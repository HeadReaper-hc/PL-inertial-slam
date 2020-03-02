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

#pragma once
#include <vector>
#include <eigen3/Eigen/Core>
#include <opencv/cv.h>

#include <DBoW2/TemplatedVocabulary.h>
#include <DBoW2/FORB.h>
#include <DBoW2/BowVector.h>
#include <DBoW2/FClass.h>
#include <DBoW2/FeatureVector.h>
#include <DBoW2/ScoringObject.h>

#include <auxiliar.h>
#include <stereoFeatures.h>
#include <stereoFrame.h>

//for imu
#include <NavState.h>
#include <imudata.h>
#include <IMUPreintegrator.h>
#include <pinholeStereoCamera.h>

using namespace cv;
using namespace std;
using namespace Eigen;
using namespace StVO;

typedef Eigen::Matrix<double,6,1> Vector6d;
typedef Eigen::Matrix<double,6,6> Matrix6d;

namespace PLSLAM{

class KeyFrame
{

public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    KeyFrame() { }
    KeyFrame( const StVO::StereoFrame* sf ,PinholeStereoCamera* cam_);
    KeyFrame( const StVO::StereoFrame* sf, int kf_idx_, PinholeStereoCamera* cam_);
    ~KeyFrame();

    Mat plotKeyFrame();

    bool     local;

    int       f_idx;
    string    img_name;

    int      kf_idx;
    Matrix4d T_kf_w;
    Vector6d x_kf_w;
    Matrix6d xcov_kf_w;

    //ground truth
    Matrix4d T_kf_w_gt;

    DBoW2::BowVector descDBoW_P, descDBoW_L;

    StVO::StereoFrame* stereo_frame;

    //for imu
    long double t_;
    NavState imuState;
    std::vector<IMUData> imus;
    IMUPreintegrator IMUPreInt;

    void ComputeIMUPreIntSinceLastFrame( KeyFrame* pLastF);
    const NavState& GetNavState(void) const;
    const IMUPreintegrator &GetIMUPreInt() { return IMUPreInt; }
    PinholeStereoCamera* cam;
    KeyFrame* RefKeyframe = nullptr;

    void UpdateNavState(const NavState& ns, const IMUPreintegrator& imupreint, const Vector3d& gw);
    void SetNavState(const NavState& ns);
    const NavState& GetNavState(void);
    void SetNavStateVel(const Vector3d &vel);
    void SetNavStatePos(const Vector3d &pos);
    void SetNavStateRot(const Matrix3d &rot);
    void SetNavStateRot(const Sophus::SO3 &rot);
    void SetNavStateBiasGyr(const Vector3d &bg);
    void SetNavStateBiasAcc(const Vector3d &ba);
    void SetNavStateDeltaBg(const Vector3d &dbg);
    void SetNavStateDeltaBa(const Vector3d &dba);
    void updateNS(const NavState& ns, const IMUPreintegrator& imupreint, const Vector3d& gw);

    void SetCameraPose(const Matrix4d& Tfw);
    void SetCameraRot(const Matrix3d& R);
    void SetCameraTans(const Vector3d& t);
    Matrix4d GetCameraPose();
    Matrix3d GetCameraRot();
    Vector3d GetCameraTrans();

protected:
    std::mutex mMutexNavState;
    std::mutex mMutexPose;
};

}
