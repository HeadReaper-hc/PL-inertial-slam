#ifndef IMUDATA_H
#define IMUDATA_H

#include <Eigen/Dense>



using namespace Eigen;

class IMUData
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // covariance of measurement
    static Matrix3d _gyrMeasCov;
    static Matrix3d _accMeasCov;
    static Matrix3d getGyrMeasCov(void) {return _gyrMeasCov;}
    static Matrix3d getAccMeasCov(void) {return _accMeasCov;}

    // covariance of bias random walk
    static Matrix3d _gyrBiasRWCov;
    static Matrix3d _accBiasRWCov;
    static Matrix3d getGyrBiasRWCov(void) {return _gyrBiasRWCov;}
    static Matrix3d getAccBiasRWCov(void) {return _accBiasRWCov;}

    static double _gyrBiasRw2;
    static double _accBiasRw2;
    static double getGyrBiasRW2(void) {return _gyrBiasRw2;}
    static double getAccBiasRW2(void) {return _accBiasRw2;}

    friend std::ostream& operator<<(std::ostream& os, const IMUData& data) {
        os << data._g.transpose() << " " << data._a.transpose() <<" "<< data._t;
        return os;
    }


    IMUData(const double& gx, const double& gy, const double& gz,
            const double& ax, const double& ay, const double& az,
            const long double& t);
    //IMUData(const IMUData& imu);

    IMUData& operator= (const IMUData& data) {
        this->_g = data._g;
        this->_a = data._a;
        this->_t = data._t;
    }

    IMUData(const IMUData& data) {
        this->_g = data._g;
        this->_a = data._a;
        this->_t = data._t;
    }

    // Raw data of imu's
    Vector3d _g;    //gyr data
    Vector3d _a;    //acc data
    long double _t;      //timestamp
};



#endif // IMUDATA_H
