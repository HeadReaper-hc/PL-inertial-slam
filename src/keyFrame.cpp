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

#include "keyFrame.h"

#include <iostream>
#include <iomanip>

#include <opencv2/imgproc.hpp>

namespace PLSLAM{

KeyFrame::KeyFrame( const StereoFrame* sf,PinholeStereoCamera* cam_)
{
    kf_idx    = -1;
    x_kf_w    = logmap_se3( T_kf_w );

    T_kf_w    = sf->Tfw ;
    xcov_kf_w = sf->Tfw_cov;
    cam = cam_;

    stereo_frame = new StereoFrame( sf->img_l, sf->img_r, kf_idx, sf->cam );
    stereo_frame->pdesc_l   = sf->pdesc_l;
    stereo_frame->pdesc_r   = sf->pdesc_r;
    stereo_frame->ldesc_l   = sf->ldesc_l;
    stereo_frame->ldesc_r   = sf->ldesc_r;
    stereo_frame->t_ = sf->t_;
    t_ = sf->t_;

    std::vector<PointFeature*> &stereo_pt = stereo_frame->stereo_pt;
    stereo_pt.resize(sf->stereo_pt.size());
    for (int idx = 0; idx < stereo_pt.size(); ++idx)
        stereo_pt[idx] = sf->stereo_pt[idx]->safeCopy();

    std::vector<LineFeature*> &stereo_ls = stereo_frame->stereo_ls;
    stereo_ls.resize(sf->stereo_ls.size());
    for (int idx = 0; idx < stereo_ls.size(); ++idx)
        stereo_ls[idx] = sf->stereo_ls[idx]->safeCopy();
}

KeyFrame::KeyFrame( const StereoFrame* sf, int kf_idx_, PinholeStereoCamera* cam_ )
{
    kf_idx    = kf_idx_;
    T_kf_w    = sf->Tfw;
    x_kf_w    = logmap_se3( T_kf_w );
    xcov_kf_w = sf->Tfw_cov;
    cam = cam_;

    stereo_frame = new StereoFrame( sf->img_l, sf->img_r, kf_idx, sf->cam );
    stereo_frame->pdesc_l   = sf->pdesc_l;
    stereo_frame->pdesc_r   = sf->pdesc_r;
    stereo_frame->ldesc_l   = sf->ldesc_l;
    stereo_frame->ldesc_r   = sf->ldesc_r;
    stereo_frame->t_ = sf->t_;
    t_ = sf->t_;
    Matrix4d tbw = T_kf_w * cam->getTbs().inverse();
    imuState.Set_Pos(tbw.block<3,1>(0,3));
    imuState.Set_Rot(tbw.block<3,3>(0,0));

    std::vector<PointFeature*> &stereo_pt = stereo_frame->stereo_pt;
    stereo_pt.resize(sf->stereo_pt.size());
    for (int idx = 0; idx < stereo_pt.size(); ++idx)
        stereo_pt[idx] = sf->stereo_pt[idx]->safeCopy();

    std::vector<LineFeature*> &stereo_ls = stereo_frame->stereo_ls;
    stereo_ls.resize(sf->stereo_ls.size());
    for (int idx = 0; idx < stereo_ls.size(); ++idx)
        stereo_ls[idx] = sf->stereo_ls[idx]->safeCopy();
}

KeyFrame::~KeyFrame() {

    delete stereo_frame;
}

Mat KeyFrame::plotKeyFrame()
{
    // create new image to modify it
    Mat img_l_aux;
    stereo_frame->img_l.copyTo( img_l_aux );
    if( img_l_aux.channels() == 1 )
        cvtColor(img_l_aux, img_l_aux, CV_GRAY2BGR, 3);
    else if (img_l_aux.channels() == 4)
        cvtColor(img_l_aux, img_l_aux, CV_BGRA2BGR, 3);
    else if (img_l_aux.channels() != 3)
        throw std::runtime_error(std::string("[KeyFrame->plotKeyFrame] unsupported image format: ") +
                                 std::to_string(img_l_aux.channels()));
    img_l_aux.convertTo(img_l_aux, CV_8UC3);

    // Variables
    unsigned int    r = 0, g, b = 0;
    Point2f         p,q;
    double          thick = 1.5;
    int             k = 0, radius  = 3;
    // plot point features
    for( vector<PointFeature*>::iterator pt_it = stereo_frame->stereo_pt.begin(); pt_it != stereo_frame->stereo_pt.end(); pt_it++)
    {
        if( (*pt_it)->idx != -1 )
        {
            g = 200;
            p = cv::Point( int((*pt_it)->pl(0)), int((*pt_it)->pl(1)) );
            circle( img_l_aux, p, radius, Scalar(b,g,r), thick);
        }
    }
    // plot line segment features
    for( vector<LineFeature*>::iterator ls_it = stereo_frame->stereo_ls.begin(); ls_it != stereo_frame->stereo_ls.end(); ls_it++)
    {
        if( (*ls_it)->idx != -1 )
        {
            g = 200;
            p = cv::Point( int((*ls_it)->spl(0)), int((*ls_it)->spl(1)) );
            q = cv::Point( int((*ls_it)->epl(0)), int((*ls_it)->epl(1)) );
            line( img_l_aux, p, q, Scalar(b,g,r), thick);
        }
    }

    return img_l_aux;
}

void KeyFrame::ComputeIMUPreIntSinceLastFrame( KeyFrame* pLastF)
{
    // Reset pre-integrator first
    IMUPreInt.reset();
    std::vector<IMUData>& imu_data = imus;
    long double prev_t = pLastF->stereo_frame->t_;
    long double curr_t = this->stereo_frame->t_;
    Vector3d bg = pLastF->GetNavState().Get_BiasGyr();
    Vector3d ba = pLastF->GetNavState().Get_BiasAcc();

    int i = 0;
    while(imu_data[i]._t < prev_t){
        i++;
    }
    double dt = imu_data[i]._t - prev_t;
    Vector3d acc = imu_data[i]._a - ba;
    Vector3d gyr = imu_data[i]._g - bg;
    IMUPreInt.update(gyr, acc, dt);
    i++;
    while(i<imu_data.size() && imu_data[i]._t<=curr_t){
        dt = imu_data[i]._t - imu_data[i-1]._t;
        acc = imu_data[i]._a - ba;
        gyr = imu_data[i]._g - bg;
        IMUPreInt.update(gyr, acc, dt);
        i++;
    }
    if(i<imu_data.size()){
        dt = curr_t - imu_data[i]._t;
        acc = imu_data[i]._a - ba;
        gyr = imu_data[i]._g - bg;
        IMUPreInt.update(gyr, acc, dt);
    }
    RefKeyframe = pLastF;
}

const NavState& KeyFrame::GetNavState(void)
{
    unique_lock<mutex> lock(mMutexNavState);
    return imuState;
}

void KeyFrame::UpdateNavState(const NavState& ns, const IMUPreintegrator& imupreint, const Vector3d& gw)
{
    unique_lock<mutex> lock(mMutexNavState);
    updateNS(ns,imupreint,gw);
}

void KeyFrame::SetNavState(const NavState& ns)
{
    unique_lock<mutex> lock(mMutexNavState);
    imuState = ns;
}

void KeyFrame::SetNavStateBiasGyr(const Vector3d &bg)
{
    unique_lock<mutex> lock(mMutexNavState);
    imuState.Set_BiasGyr(bg);
}

void KeyFrame::SetNavStateBiasAcc(const Vector3d &ba)
{
    unique_lock<mutex> lock(mMutexNavState);
    imuState.Set_BiasAcc(ba);
}

void KeyFrame::SetNavStateVel(const Vector3d &vel)
{
    unique_lock<mutex> lock(mMutexNavState);
    imuState.Set_Vel(vel);
}

void KeyFrame::SetNavStatePos(const Vector3d &pos)
{
    unique_lock<mutex> lock(mMutexNavState);
    imuState.Set_Pos(pos);
}

void KeyFrame::SetNavStateRot(const Matrix3d &rot)
{
    unique_lock<mutex> lock(mMutexNavState);
    imuState.Set_Rot(rot);
}

void KeyFrame::SetNavStateRot(const Sophus::SO3 &rot)
{
    unique_lock<mutex> lock(mMutexNavState);
    imuState.Set_Rot(rot);
}

void KeyFrame::SetNavStateDeltaBg(const Vector3d &dbg)
{
    unique_lock<mutex> lock(mMutexNavState);
    imuState.Set_DeltaBiasGyr(dbg);
}

void KeyFrame::SetNavStateDeltaBa(const Vector3d &dba)
{
    unique_lock<mutex> lock(mMutexNavState);
    imuState.Set_DeltaBiasAcc(dba);
}

void KeyFrame::updateNS(const NavState& ns, const IMUPreintegrator& imupreint, const Vector3d& gw)
{
    Matrix3d dR = imupreint.getDeltaR();
    Vector3d dP = imupreint.getDeltaP();
    Vector3d dV = imupreint.getDeltaV();
    double dt = imupreint.getDeltaTime();

    Vector3d Pwbpre = ns.Get_P();
    Matrix3d Rwbpre = ns.Get_RotMatrix();
    Vector3d Vwbpre = ns.Get_V();

    Matrix3d Rwb = Rwbpre * dR;
    Vector3d Pwb = Pwbpre + Vwbpre*dt + 0.5*gw*dt*dt + Rwbpre*dP;
    Vector3d Vwb = Vwbpre + gw*dt + Rwbpre*dV;

    // Here assume that the pre-integration is re-computed after bias updated, so the bias term is ignored
    imuState.Set_Pos(Pwb);
    imuState.Set_Vel(Vwb);
    imuState.Set_Rot(Rwb);

    // Test log
    if(ns.Get_dBias_Gyr().norm()>1e-6 || ns.Get_dBias_Acc().norm()>1e-6) std::cerr<<"delta bias in updateNS is not zero"<<ns.Get_dBias_Gyr().transpose()<<", "<<ns.Get_dBias_Acc().transpose()<<std::endl;
}

void KeyFrame::SetCameraPose(const Matrix4d &Tfw) {
    unique_lock<mutex> lock(mMutexPose);
    T_kf_w = Tfw;
}

void KeyFrame::SetCameraRot(const Matrix3d &R) {
    unique_lock<mutex> lock(mMutexPose);
    T_kf_w.block<3,3>(0,0) = R;
}

void KeyFrame::SetCameraTans(const Vector3d &t) {
    unique_lock<mutex> lock(mMutexPose);
    T_kf_w.block<3,1>(0,3);
}

Matrix4d KeyFrame::GetCameraPose(){
    unique_lock<mutex> lock(mMutexPose);
    return T_kf_w;
}

Matrix3d KeyFrame::GetCameraRot(){
    unique_lock<mutex> lock(mMutexPose);
    return T_kf_w.block<3,3>(0,0);
}

Vector3d KeyFrame::GetCameraTrans(){
    unique_lock<mutex> lock(mMutexPose);
    return T_kf_w.block<3,1>(0,3);
}

}
