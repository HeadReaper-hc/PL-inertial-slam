/*****************************************************************************
**      Stereo VO and SLAM by combining point and line segment features     **
******************************************************************************
**                                                                          **
**  Copyright(c) 2016-2018, Ruben Gomez-Ojeda, University of Malaga         **
**  Copyright(c) 2016-2018, David Zuñiga-Noël, University of Malaga         **
**  Copyright(c) 2016-2018, MAPIR group, University of Malaga               **
**                                                                          **
**  This program is free software: you can redistribute it and/or modify    **
**  it under the terms of the GNU General Public License (version 3) as     **
**  published by the Free Software Foundation.                              **
**                                                                          **
**  This program is distributed in the hope that it will be useful, but     **
**  WITHOUT ANY WARRANTY; without even the implied warranty of              **
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            **
**  GNU General Public License for more details.                            **
**                                                                          **
**  You should have received a copy of the GNU General Public License       **
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.   **
**                                                                          **
*****************************************************************************/

#include "dataset.h"

//STL
#include <algorithm>
#include <functional>
#include <limits>
#include <list>
#include <stdexcept>
#include <string>

//Boost
#include <boost/regex.hpp> //Note: using boost regex instead of C++11 regex as it isn't supported by the compiler until gcc 4.9
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>

//OpenCV
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>

//YAML
#include <yaml-cpp/yaml.h>

#include "pinholeStereoCamera.h"

//imu
#include <imudata.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

namespace StVO {

void getSortedImages(const boost::filesystem::path &img_dir, std::function<bool(const std::string &)> filter,
                     std::function<bool(const std::string &, const std::string &)> comparator, std::vector<std::string> &img_paths) {

    // get a sorted list of files in the img directories
    if (!boost::filesystem::exists(img_dir) ||
            !boost::filesystem::is_directory(img_dir))
        throw std::runtime_error("[Dataset] Invalid images subfolder");

    // get all files in the img directories
    std::list<std::string> all_imgs;
    for(auto& entry : boost::make_iterator_range(boost::filesystem::directory_iterator(img_dir), {})) {
        boost::filesystem::path filename_path = entry.path().filename();
        if (boost::filesystem::is_regular_file(entry.status()) &&
                (filename_path.extension() == ".png"  ||
                 filename_path.extension() == ".jpg"  ||
                 filename_path.extension() == ".jpeg" ||
                 filename_path.extension() == ".pnm"  ||
                 filename_path.extension() == ".tiff")) {
            all_imgs.push_back(filename_path.string());
        }
    }

    // sort
    img_paths.clear();
    img_paths.reserve(all_imgs.size());
    for (const std::string &filename : all_imgs)
        if (!filter(filename)) img_paths.push_back(filename);

    if (img_paths.empty())
        throw std::runtime_error("[Dataset] Invalid image names?");

    std::sort(img_paths.begin(), img_paths.end(), comparator);

    for (std::string &filename : img_paths)
        filename = (img_dir / filename).string();
}

Dataset::Dataset(const std::string &dataset_path, const PinholeStereoCamera &cam, int offset, int nmax, int step)
    : cam(cam) {

    boost::filesystem::path dataset_base(dataset_path);
    if (!boost::filesystem::exists(dataset_base) ||
            !boost::filesystem::is_directory(dataset_base))
        throw std::runtime_error("[Dataset] Invalid directory");

    boost::filesystem::path dataset_params = dataset_base / "dataset_params.yaml";
    if (!boost::filesystem::exists(dataset_params))
        throw std::runtime_error("[Dataset] Dataset parameters not found");
    YAML::Node dataset_config = YAML::LoadFile(dataset_params.string());

    // setup image directories
    boost::filesystem::path img_l_dir = dataset_base / dataset_config["images_subfolder_l"].as<std::string>();
    boost::filesystem::path img_r_dir = dataset_base / dataset_config["images_subfolder_r"].as<std::string>();

    boost::regex expression("^[^0-9]*([0-9]+\\.?+[0-9]*)[^0-9]*\\.[a-z]{3,4}$");
    boost::cmatch what;

    auto filename_filter = [&expression, &what](const std::string &s) {
        return !boost::regex_match(s.c_str(), what, expression);
    };

    auto sort_by_number = [&expression, &what](const std::string &a, const std::string &b) {
        double n1, n2;

        if (boost::regex_match(a.c_str(), what, expression))
            n1 = std::stod(what[1]);
        else
            throw std::runtime_error("[Dataset] Unexpected behaviour while sorting filenames");

        if (boost::regex_match(b.c_str(), what, expression))
            n2 = std::stod(what[1]);
        else
            throw std::runtime_error("[Dataset] Unexpected behaviour while sorting filenames");

        return (n1 < n2);
    };

    std::vector<std::string> img_l_paths, img_r_paths;
    getSortedImages(img_l_dir, filename_filter, sort_by_number, img_l_paths);
    getSortedImages(img_r_dir, filename_filter, sort_by_number, img_r_paths);

    if (img_l_paths.size() != img_r_paths.size())
        throw std::runtime_error("[Dataset] Left and right images");

    // decimate sequence
    offset = std::max(0, offset);
    nmax = (nmax <= 0) ? std::numeric_limits<int>::max() : nmax;
    step = std::max(1, step);
    for (int i = 0, ctr = 0; (i + offset) < img_l_paths.size() && ctr < nmax; i += step, ctr++) {
        images_l.push_back(img_l_paths[i + offset]);
        images_r.push_back(img_r_paths[i + offset]);
    }

    //read image time and imu
    std::cout<<"read image time"<<endl;
    boost::filesystem::path image_time_dir = dataset_base / dataset_config["image_time"].as<std::string>();
    std::ifstream fin(image_time_dir.string(),ios::in);

    if(!fin)
    {
        std::cerr<<"Can't open image time file, exit!"<<std::endl;
        exit(-1);
    }
    image_time_name.clear();
    vector<pair<long double,string>> img_t_name;
    string line;
    bool first_row = true;
    while(getline(fin,line)){
        if(first_row){
            first_row = false;
            continue;
        }
        istringstream sin(line);
        vector<string> fields;
        string field;
        while(getline(sin,field,',')){
            fields.push_back(field);
        }
        string s = fields[0].substr(0,10);
        string ns = fields[0].substr(10,9);
        istringstream ss(s), nss(ns);
        long double ls, lns;
        ss >> ls; nss >> lns;
        long double time = ls + lns * 1e-9;
        img_t_name.push_back(make_pair(time,fields[1]));
    }
    fin.close();
    for (int i = 0, ctr = 0; (i + offset) < img_t_name.size() && ctr < nmax; i += step, ctr++) {
        image_time_name.push_back(img_t_name[i + offset]);
    }

    //check time and image consistent
    if(image_time_name.size()!=images_l.size()){
        throw std::runtime_error("different size of image and image time");
    }

    list<pair<long double,string>>::iterator time_iter = image_time_name.begin();
    list<string>::iterator image_iter = images_l.begin();
    std::cout<<"check image time"<<std::endl;
    for(;image_iter!=images_l.end();image_iter++,time_iter++){
        string image_time_str =  string(img_l_dir.string()+(*time_iter).second);
        image_time_str = image_time_str.substr(0,(*image_iter).size());
//            std::cout<<(*image_iter)<<std::endl;
//            std::cout<<image_time_str<<std::endl;
        if(image_time_str!=(*image_iter)){
            throw std::runtime_error("wrong match between image time and image name");
        }
        img_time.push_back((*time_iter).first);
    }
    std::cout<<"check image time done...."<<std::endl;
    //test
//        list<long double>::iterator img_time_iter = img_time.begin();
//        for(int i=0;i<4000;img_time_iter++,i++){
//            cout.precision(20);
//            std::cout<<i<<"  "<<(*img_time_iter)<<std::endl;
//        }

    boost::filesystem::path imu_dir = dataset_base / dataset_config["imu_dir"].as<std::string>();
    std::ifstream imufile(imu_dir.string(),ios::in);

    if(!imufile)
    {
        std::cerr<<"Can't open imu file, exit!"<<std::endl;
        exit(-1);
    }

    first_row = true;
    imus.clear();
    while(getline(imufile,line)){
        if(first_row){
            first_row = false;
            continue;
        }
        istringstream sin(line);
        vector<string> fields;
        string field;
        while(getline(sin,field,',')){
            fields.push_back(field);
        }
        istringstream is2(fields[1]),is3(fields[2]),is4(fields[3]),is5(fields[4]),is6(fields[5]),is7(fields[6]);
        double gx,gy,gz,ax,ay,az;
        is2 >> gx; is3 >> gy; is4 >> gz; is5 >> ax; is6 >> ay; is7 >> az;
        string s = fields[0].substr(0,10);
        string ns = fields[0].substr(10,9);
        istringstream ss(s), nss(ns);
        long double ls, lns;
        ss >> ls; nss >> lns;
        long double time = ls + lns * 1e-9;
        imus.push_back(IMUData(gx,gy,gz,ax,ay,az,time));
    }
    lastImuIndex = 0;
    //test
//        vector<IMUData>::iterator imu_iter = imus.begin();
//        for(int i=0;i<4000;imu_iter++,i++){
//            IMUData temp = *imu_iter;
//            //std::cout<<i<<"  "<<temp<<std::endl;
//            std::cout<<i<<"     "<<temp._t<<std::endl;
//        }


    //get ground truth data
    boost::filesystem::path gt_dir = dataset_base / dataset_config["gt_dir"].as<std::string>();
    std::ifstream gtfile(gt_dir.string(),ios::in);
    if(!gtfile)
    {
        std::cout<<"The ground truth file not exist!"<<std::endl;
    }
    else{
        std::cout<<"The ground truth file exist, now read ground truth begin...."<<std::endl;
        gt_time.clear();
        gt_data.clear();
        string gtline;
        getline(gtfile, gtline);
        while(getline(gtfile,gtline)){
            istringstream sin(gtline);
            vector<string> fields;
            string field;
            int i=0;
            while(getline(sin,field,',')){
                if(i<8) {
                    fields.push_back(field);
                    i++;
                }
                else{
                    break;
                }
            }
            istringstream is2(fields[1]),is3(fields[2]),is4(fields[3]),is5(fields[4]),is6(fields[5]),is7(fields[6]),is8(fields[7]);
            double px,py,pz,qw,qx,qy,qz;
            is2 >> px; is3 >> py; is4 >> pz; is5 >> qw; is6 >> qx; is7 >> qy; is8 >> qz;
            string s = fields[0].substr(0,10);
            string ns = fields[0].substr(10,9);
            istringstream ss(s), nss(ns);
            long double ls, lns;
            ss >> ls; nss >> lns;
            long double time = ls + lns * 1e-9;
            gt_time.push_back(time);
            vector<double> pose;
            pose.push_back(px); pose.push_back(py); pose.push_back(pz);
            pose.push_back(qw); pose.push_back(qx); pose.push_back(qy); pose.push_back(qz);
            gt_data.push_back(pose);
        }
        have_gt = true;
        lastGtIndex = 0;
        cout<<"Read ground truth data done.... ALL size of gt is "<<gt_time.size()<<endl;
        //test
//        for(int i=0;i<10;i++){
//            cout<<"time: "<<gt_time[i]<<endl;
//            cout<<"Pose: ";
//            for(int j=0;j<7;j++){
//                cout<<gt_data[i][j]<<" ";
//            }
//            cout<<endl;
//        }
    }
}

Dataset::~Dataset() {

}

bool Dataset::nextFrame(cv::Mat &img_l, cv::Mat &img_r) {
    if (!hasNext()) return false;

    img_l = cv::imread(images_l.front(), CV_LOAD_IMAGE_UNCHANGED);
    img_r = cv::imread(images_r.front(), CV_LOAD_IMAGE_UNCHANGED);
    cam.rectifyImagesLR(img_l, img_l, img_r, img_r);
    images_l.pop_front();
    images_r.pop_front();

    return (!img_l.empty() && !img_r.empty());
}

bool Dataset::nextFrame(cv::Mat &img_l, cv::Mat &img_r, long double &t) {
    if (!hasNext()) return false;

    img_l = cv::imread(images_l.front(), CV_LOAD_IMAGE_UNCHANGED);
    img_r = cv::imread(images_r.front(), CV_LOAD_IMAGE_UNCHANGED);
    t = img_time.front();
    cam.rectifyImagesLR(img_l, img_l, img_r, img_r);
    images_l.pop_front();
    images_r.pop_front();
    img_time.pop_front();

    return (!img_l.empty() && !img_r.empty());
}

bool Dataset::hasNext() {
    return !(images_l.empty() || images_r.empty());
}

//for imu
std::vector<IMUData> Dataset::getIMUDataBetweenTwoImage(long double prevImageTime, long double currImageTime) {
    std::vector<IMUData> result;
    while(imus[lastImuIndex]._t <= prevImageTime){
        lastImuIndex++;
    }
    lastImuIndex--;
    if(lastImuIndex<0)
        lastImuIndex = 0;
    while(imus[lastImuIndex]._t< currImageTime){
        result.push_back(imus[lastImuIndex]);
        lastImuIndex++;
    }
    result.push_back(imus[lastImuIndex]);
    lastImuIndex--;
    return result;
}

//for ground truth
Matrix4d Dataset::getGroundTruthByTime(const long double kft){
    if(gt_time.size()<=0)
        return Matrix4d::Identity();

    if(kft<gt_time[0]){
        vector<double> temp = gt_data[0];
        Vector3d trans(temp[0],temp[1],temp[2]);
        Quaterniond q(temp[3],temp[4],temp[5],temp[6]);
        Matrix4d gt = Matrix4d::Identity();
        gt.block<3,3>(0,0) = q.toRotationMatrix();
        gt.block<3,1>(0,3) = trans;
        return gt;
    }

    if(kft>gt_time.back()){
        vector<double> temp = gt_data.back();
        Vector3d trans(temp[0],temp[1],temp[2]);
        Quaterniond q(temp[3],temp[4],temp[5],temp[6]);
        Matrix4d gt = Matrix4d::Identity();
        gt.block<3,3>(0,0) = q.toRotationMatrix();
        gt.block<3,1>(0,3) = trans;
        return gt;
    }

    while(gt_time[lastGtIndex]<kft){
        lastGtIndex++;
    }
    vector<double> temp = gt_data[lastGtIndex];
    Vector3d trans(temp[0],temp[1],temp[2]);
    Quaterniond q(temp[3],temp[4],temp[5],temp[6]);
    Matrix4d gt = Matrix4d::Identity();
    gt.block<3,3>(0,0) = q.toRotationMatrix();
    gt.block<3,1>(0,3) = trans;
    return gt;
}

} // namespace StVO

