#ifndef G2OTYPES_H
#define G2OTYPES_H

#include <g2o/core/base_vertex.h>
#include <g2o/core/base_unary_edge.h>
#include "so3.h"
#include "NavState.h"
#include "IMUPreintegrator.h"
#include <g2o/core/base_multi_edge.h>
#include <g2o/core/base_binary_edge.h>
#include "types_six_dof_expmap.h"
#include <g2o/solvers/eigen/linear_solver_eigen.h>
#include <g2o/core/block_solver.h>
#include <iostream>
//#include "Thirdparty/g2o/g2o/core/sparse_block_matrix.h"
typedef g2o::LinearSolverEigen<g2o::BlockSolverX::PoseMatrixType> SlamLinearSolver;
//#include "Thirdparty/g2o/g2o/core/sparse_block_matrix.h"

namespace g2o {


/**
 * @brief The VertexNavStatePVR class
 */
    class VertexNavStatePVR : public BaseVertex<9, NavState> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        VertexNavStatePVR() : BaseVertex<9, NavState>() {}

        bool read(std::istream &is) { return true; }

        bool write(std::ostream &os) const { return true; }

        virtual void setToOriginImpl() {
            _estimate = NavState();
        }

        virtual void oplusImpl(const double *update_) {
            Eigen::Map<const Vector9d> update(update_);
            _estimate.IncSmallPVR(update);
        }

    };

    class VertexNavStateBias : public BaseVertex<6, NavState> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        VertexNavStateBias() : BaseVertex<6, NavState>() {}

        bool read(std::istream &is) { return true; }

        bool write(std::ostream &os) const { return true; }

        virtual void setToOriginImpl() {
            _estimate = NavState();
        }

        virtual void oplusImpl(const double *update_) {
            Eigen::Map<const Vector6d> update(update_);
            _estimate.IncSmallBias(update);
        }

    };

    class EdgeNavStatePVR : public BaseMultiEdge<9, IMUPreintegrator> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeNavStatePVR() : BaseMultiEdge<9, IMUPreintegrator>() {
            resize(3);
        }

        bool read(std::istream &is) { return true; }

        bool write(std::ostream &os) const { return true; }

        void computeError();

        virtual void linearizeOplus();

        void SetParams(const Vector3d &gw) {
            GravityVec = gw;
        }

    protected:
        // Gravity vector in 'world' frame
        Vector3d GravityVec;
    };

    class EdgeNavStateBias : public BaseBinaryEdge<6, IMUPreintegrator, VertexNavStateBias, VertexNavStateBias> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeNavStateBias() : BaseBinaryEdge<6, IMUPreintegrator, VertexNavStateBias, VertexNavStateBias>() {}

        bool read(std::istream &is) { return true; }

        bool write(std::ostream &os) const { return true; }

        void computeError();

        virtual void linearizeOplus();

    };

    class EdgeNavStatePVRPointXYZ : public BaseBinaryEdge<2, Vector2d, VertexSBAPointXYZ, VertexNavStatePVR> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeNavStatePVRPointXYZ() : BaseBinaryEdge<2, Vector2d, VertexSBAPointXYZ, VertexNavStatePVR>() {}

        bool read(std::istream &is) { return true; }

        bool write(std::ostream &os) const { return true; }

        void computeError() {
            Vector3d Pc = computePc();
            Vector2d obs(_measurement);

            _error = obs - cam_project(Pc);
            //std::cout<<"_error: "<<_error<<std::endl;
        }

        bool isDepthPositive() {
            Vector3d Pc = computePc();
            return Pc(2) > 0.0;
        }

        Vector3d computePc() {
            const VertexSBAPointXYZ *vPoint = static_cast<const VertexSBAPointXYZ *>(_vertices[0]);
            const VertexNavStatePVR *vNavState = static_cast<const VertexNavStatePVR *>(_vertices[1]);

            const NavState &ns = vNavState->estimate();
            Matrix3d Rwb = ns.Get_RotMatrix();
            Vector3d Pwb = ns.Get_P();
            const Vector3d &Pw = vPoint->estimate();

            Matrix3d Rcb = Rbc.transpose();
            Vector3d Pc = Rcb * Rwb.transpose() * (Pw - Pwb) - Rcb * Pbc;

            return Pc;
            //Vector3d Pwc = Rwb*Pbc + Pwb;
            //Matrix3d Rcw = (Rwb*Rbc).transpose();
            //Vector3d Pcw = -Rcw*Pwc;
            //Vector3d Pc = Rcw*Pw + Pcw;
        }

        inline Vector2d project2d(const Vector3d &v) const {
            Vector2d res;
            res(0) = v(0) / v(2);
            res(1) = v(1) / v(2);
            return res;
        }

        Vector2d cam_project(const Vector3d &trans_xyz) const {
            Vector2d proj = project2d(trans_xyz);
            Vector2d res;
            res[0] = proj[0] * fx + cx;
            res[1] = proj[1] * fy + cy;
            return res;
        }

        //
        virtual void linearizeOplus();

        void SetParams(const double &fx_, const double &fy_, const double &cx_, const double &cy_,
                       const Matrix3d &Rbc_, const Vector3d &Pbc_) {
            fx = fx_;
            fy = fy_;
            cx = cx_;
            cy = cy_;
            Rbc = Rbc_;
            Pbc = Pbc_;
        }

    protected:
        // Camera intrinsics
        double fx, fy, cx, cy;
        // Camera-IMU extrinsics
        Matrix3d Rbc;
        Vector3d Pbc;
    };

    class EdgeNavStatePVRPointXYZOnlyPose : public BaseUnaryEdge<2, Vector2d, VertexNavStatePVR> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeNavStatePVRPointXYZOnlyPose() {}

        bool read(std::istream &is) { return true; }

        bool write(std::ostream &os) const { return true; }

        void computeError() {
            Vector3d Pc = computePc();
            Vector2d obs(_measurement);

            _error = obs - cam_project(Pc);
        }

        bool isDepthPositive() {
            Vector3d Pc = computePc();
            return Pc(2) > 0.0;
        }

        Vector3d computePc() {
            const VertexNavStatePVR *vNSPVR = static_cast<const VertexNavStatePVR *>(_vertices[0]);

            const NavState &ns = vNSPVR->estimate();
            Matrix3d Rwb = ns.Get_RotMatrix();
            Vector3d Pwb = ns.Get_P();
            //const Vector3d& Pw = vPoint->estimate();

            Matrix3d Rcb = Rbc.transpose();
            Vector3d Pc = Rcb * Rwb.transpose() * (Pw - Pwb) - Rcb * Pbc;

            return Pc;
            //Vector3d Pwc = Rwb*Pbc + Pwb;
            //Matrix3d Rcw = (Rwb*Rbc).transpose();
            //Vector3d Pcw = -Rcw*Pwc;
            //Vector3d Pc = Rcw*Pw + Pcw;
        }

        inline Vector2d project2d(const Vector3d &v) const {
            Vector2d res;
            res(0) = v(0) / v(2);
            res(1) = v(1) / v(2);
            return res;
        }

        Vector2d cam_project(const Vector3d &trans_xyz) const {
            Vector2d proj = project2d(trans_xyz);
            Vector2d res;
            res[0] = proj[0] * fx + cx;
            res[1] = proj[1] * fy + cy;
            return res;
        }

        //
        virtual void linearizeOplus();

        void SetParams(const double &fx_, const double &fy_, const double &cx_, const double &cy_,
                       const Matrix3d &Rbc_, const Vector3d &Pbc_, const Vector3d &Pw_) {
            fx = fx_;
            fy = fy_;
            cx = cx_;
            cy = cy_;
            Rbc = Rbc_;
            Pbc = Pbc_;
            Pw = Pw_;
        }

    protected:
        // Camera intrinsics
        double fx, fy, cx, cy;
        // Camera-IMU extrinsics
        Matrix3d Rbc;
        Vector3d Pbc;
        // Point position in world frame
        Vector3d Pw;
    };

/**
 * @brief The EdgeNavStatePrior class
 */
    class EdgeNavStatePriorPVRBias : public BaseBinaryEdge<15, NavState, VertexNavStatePVR, VertexNavStateBias> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeNavStatePriorPVRBias() : BaseBinaryEdge<15, NavState, VertexNavStatePVR, VertexNavStateBias>() {}

        bool read(std::istream &is) { return true; }

        bool write(std::ostream &os) const { return true; }

        void computeError();

        virtual void linearizeOplus();

    };

//------------------------------------------

/**
 * @brief The VertexNavState class
 * Vertex of tightly-coupled Visual-Inertial optimization
 */
    class VertexNavState : public BaseVertex<15, NavState> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        VertexNavState();

        bool read(std::istream &is);

        bool write(std::ostream &os) const;

        virtual void setToOriginImpl() {
            _estimate = NavState();
        }

        virtual void oplusImpl(const double *update_);
    };

/**
 * @brief The EdgeNavStatePrior class
 */
    class EdgeNavStatePrior : public BaseUnaryEdge<15, NavState, VertexNavState> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeNavStatePrior() {}

        bool read(std::istream &is) { return true; }

        bool write(std::ostream &os) const { return true; }

        void computeError();

        virtual void linearizeOplus();

    };

/**
 * @brief The VertexGravityW class
 */
    class VertexGravityW : public BaseVertex<2, Vector3d> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        VertexGravityW();

        bool read(std::istream &is) { return true; }

        bool write(std::ostream &os) const { return true; }

        virtual void setToOriginImpl() {
            _estimate = Vector3d(0, 0, 9.81);
        }

        virtual void oplusImpl(const double *update_);
    };

/**
 * @brief The EdgeNavStateGw class
 */
    class EdgeNavStateGw : public BaseMultiEdge<15, IMUPreintegrator> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeNavStateGw();

        bool read(std::istream &is) { return true; }

        bool write(std::ostream &os) const { return true; }

        void computeError();

        virtual void linearizeOplus();
    };

/**
 * @brief The EdgeNavState class
 * Edge between NavState_i and NavState_j, vertex[0]~i, vertex[1]~j
 * Measurement~Vector15d: 9Dof-IMUPreintegrator measurement & 6Dof-IMU bias change all Zero
 */
    class EdgeNavState : public BaseBinaryEdge<15, IMUPreintegrator, VertexNavState, VertexNavState> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeNavState();

        bool read(std::istream &is);

        bool write(std::ostream &os) const;

        void computeError();

        virtual void linearizeOplus();

        void SetParams(const Vector3d &gw) {
            GravityVec = gw;
        }

    protected:
        // Gravity vector in 'world' frame
        Vector3d GravityVec;
    };

/**
 * @brief The EdgeNavStatePointXYZ class
 * Edge between NavState and Point3D, vertex[0]~Point3D, vertex[1]~NavState
 * Measurement~Vector2d: 2Dof image feature position
 */
    class EdgeNavStatePointXYZ : public BaseBinaryEdge<2, Vector2d, VertexSBAPointXYZ, VertexNavState> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeNavStatePointXYZ();

        bool read(std::istream &is);

        bool write(std::ostream &os) const;

        void computeError() {
            Vector3d Pc = computePc();
            Vector2d obs(_measurement);

            _error = obs - cam_project(Pc);
        }

        bool isDepthPositive() {
            Vector3d Pc = computePc();
            return Pc(2) > 0.0;
        }

        Vector3d computePc() {
            const VertexSBAPointXYZ *vPoint = static_cast<const VertexSBAPointXYZ *>(_vertices[0]);
            const VertexNavState *vNavState = static_cast<const VertexNavState *>(_vertices[1]);

            const NavState &ns = vNavState->estimate();
            Matrix3d Rwb = ns.Get_RotMatrix();
            Vector3d Pwb = ns.Get_P();
            const Vector3d &Pw = vPoint->estimate();

            Matrix3d Rcb = Rbc.transpose();
            Vector3d Pc = Rcb * Rwb.transpose() * (Pw - Pwb) - Rcb * Pbc;

            return Pc;
            //Vector3d Pwc = Rwb*Pbc + Pwb;
            //Matrix3d Rcw = (Rwb*Rbc).transpose();
            //Vector3d Pcw = -Rcw*Pwc;
            //Vector3d Pc = Rcw*Pw + Pcw;
        }

        inline Vector2d project2d(const Vector3d &v) const {
            Vector2d res;
            res(0) = v(0) / v(2);
            res(1) = v(1) / v(2);
            return res;
        }

        Vector2d cam_project(const Vector3d &trans_xyz) const {
            Vector2d proj = project2d(trans_xyz);
            Vector2d res;
            res[0] = proj[0] * fx + cx;
            res[1] = proj[1] * fy + cy;
            return res;
        }

        //
        virtual void linearizeOplus();

        void SetParams(const double &fx_, const double &fy_, const double &cx_, const double &cy_,
                       const Matrix3d &Rbc_, const Vector3d &Pbc_) {
            fx = fx_;
            fy = fy_;
            cx = cx_;
            cy = cy_;
            Rbc = Rbc_;
            Pbc = Pbc_;
        }

    protected:
        // Camera intrinsics
        double fx, fy, cx, cy;
        // Camera-IMU extrinsics
        Matrix3d Rbc;
        Vector3d Pbc;
    };

    class EdgeNavStatePointXYZOnlyPose : public BaseUnaryEdge<2, Vector2d, VertexNavState> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeNavStatePointXYZOnlyPose() {}

        bool read(std::istream &is) { return true; }

        bool write(std::ostream &os) const { return true; }

        void computeError() {
            Vector3d Pc = computePc();
            Vector2d obs(_measurement);

            _error = obs - cam_project(Pc);

        }

        bool isDepthPositive() {
            Vector3d Pc = computePc();
            return Pc(2) > 0.0;
        }

        Vector3d computePc() {
            const VertexNavState *vNavState = static_cast<const VertexNavState *>(_vertices[0]);

            const NavState &ns = vNavState->estimate();
            Matrix3d Rwb = ns.Get_RotMatrix();
            Vector3d Pwb = ns.Get_P();
            //const Vector3d& Pw = vPoint->estimate();

            Matrix3d Rcb = Rbc.transpose();
            Vector3d Pc = Rcb * Rwb.transpose() * (Pw - Pwb) - Rcb * Pbc;

            return Pc;
            //Vector3d Pwc = Rwb*Pbc + Pwb;
            //Matrix3d Rcw = (Rwb*Rbc).transpose();
            //Vector3d Pcw = -Rcw*Pwc;
            //Vector3d Pc = Rcw*Pw + Pcw;
        }

        inline Vector2d project2d(const Vector3d &v) const {
            Vector2d res;
            res(0) = v(0) / v(2);
            res(1) = v(1) / v(2);
            return res;
        }

        Vector2d cam_project(const Vector3d &trans_xyz) const {
            Vector2d proj = project2d(trans_xyz);
            Vector2d res;
            res[0] = proj[0] * fx + cx;
            res[1] = proj[1] * fy + cy;
            return res;
        }

        //
        virtual void linearizeOplus();

        void SetParams(const double &fx_, const double &fy_, const double &cx_, const double &cy_,
                       const Matrix3d &Rbc_, const Vector3d &Pbc_, const Vector3d &Pw_) {
            fx = fx_;
            fy = fy_;
            cx = cx_;
            cy = cy_;
            Rbc = Rbc_;
            Pbc = Pbc_;
            Pw = Pw_;
        }

    protected:
        // Camera intrinsics
        double fx, fy, cx, cy;
        // Camera-IMU extrinsics
        Matrix3d Rbc;
        Vector3d Pbc;
        // Point position in world frame
        Vector3d Pw;
    };


/**
 * @brief The VertexGyrBias class
 * For gyroscope bias compuation in Visual-Inertial initialization
 */
    class VertexGyrBias : public BaseVertex<3, Vector3d> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        VertexGyrBias();

        bool read(std::istream &is);

        bool write(std::ostream &os) const;

        virtual void setToOriginImpl() {
            _estimate.setZero();
        }

        virtual void oplusImpl(const double *update_);
    };

/**
 * @brief The EdgeGyrBias class
 * For gyroscope bias compuation in Visual-Inertial initialization
 */
    class EdgeGyrBias : public BaseUnaryEdge<3, Vector3d, VertexGyrBias> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeGyrBias();

        bool read(std::istream &is);

        bool write(std::ostream &os) const;

        Matrix3d dRbij;
        Matrix3d J_dR_bg;
        Matrix3d Rwbi;
        Matrix3d Rwbj;

        void computeError();

        virtual void linearizeOplus();
    };

/**
 * @brief The VertexLine class
 * For Line error in pl-slam
 */
    class VertexLine : public BaseVertex<6, Vector6d> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        VertexLine();

        bool read(std::istream &is) { return true; };

        bool write(std::ostream &os) const { return true; };

        virtual void setToOriginImpl() {
            _estimate.setZero();
        }

        virtual void oplusImpl(const double *update_);
    };

/**
 * @brief The EdgeNavStateLine
 * For Line error in pl-slam
 */
    class EdgeNavStateLine : public BaseBinaryEdge<3, Vector3d, VertexLine, VertexNavStatePVR> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeNavStateLine() : BaseBinaryEdge<3, Vector3d, VertexLine, VertexNavStatePVR>() {};

        bool read(std::istream &is) { return true; };

        bool write(std::ostream &os) const { return true; };

        void computeError() {
            Vector6d Pc = computePc();
            Vector3d obs(_measurement);

            Vector2d UVs = cam_project(Pc.head(3));
            Vector2d UVe = cam_project(Pc.tail(3));

            _error(0) = obs(0) * UVs(0) + obs(1) * UVs(1) + obs(2);
            _error(1) = obs(0) * UVe(0) + obs(1) * UVe(1) + obs(2);
            _error(2) = 0;
        }

        bool isDepthPositive() {
            Vector6d Pc = computePc();
            return Pc(2) > 0.0 && Pc(5) > 0.0;
        }

        Vector6d computePc() {
            const VertexLine *vPoint = static_cast<const VertexLine *>(_vertices[0]);
            const VertexNavStatePVR *vNavState = static_cast<const VertexNavStatePVR *>(_vertices[1]);

            const NavState &ns = vNavState->estimate();
            Matrix3d Rwb = ns.Get_RotMatrix();
            Vector3d Pwb = ns.Get_P();
            const Vector6d &Pw = vPoint->estimate();

            Vector3d pw_s, pw_e;
            pw_s = Pw.head(3);
            pw_e = Pw.tail(3);

            Matrix3d Rcb = Rbc.transpose();
            Vector3d Pc_s = Rcb * Rwb.transpose() * (pw_s - Pwb) - Rcb * Pbc;
            Vector3d Pc_e = Rcb * Rwb.transpose() * (pw_e - Pwb) - Rcb * Pbc;

            Vector6d Pc;
            Pc << Pc_s, Pc_e;

            return Pc;
            //Vector3d Pwc = Rwb*Pbc + Pwb;
            //Matrix3d Rcw = (Rwb*Rbc).transpose();
            //Vector3d Pcw = -Rcw*Pwc;
            //Vector3d Pc = Rcw*Pw + Pcw;
        }

        inline Vector2d project2d(const Vector3d &v) const {
            Vector2d res;
            res(0) = v(0) / v(2);
            res(1) = v(1) / v(2);
            return res;
        }

        Vector2d cam_project(const Vector3d &trans_xyz) const {
            Vector2d proj = project2d(trans_xyz);
            Vector2d res;
            res[0] = proj[0] * fx + cx;
            res[1] = proj[1] * fy + cy;
            return res;
        }

        //
        virtual void linearizeOplus();

        void SetParams(const double &fx_, const double &fy_, const double &cx_, const double &cy_,
                       const Matrix3d &Rbc_, const Vector3d &Pbc_) {
            fx = fx_;
            fy = fy_;
            cx = cx_;
            cy = cy_;
            Rbc = Rbc_;
            Pbc = Pbc_;
        }

    protected:
        // Camera intrinsics
        double fx, fy, cx, cy;
        // Camera-IMU extrinsics
        Matrix3d Rbc;
        Vector3d Pbc;
    };

    /**
      * @brief The VertexLinePoint class
      * For Line error in pl-slam
      */
    class VertexLinePoint : public BaseVertex<3, Vector3d> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        VertexLinePoint() : BaseVertex<3, Vector3d>() {}

        bool read(std::istream &is) { return true; };

        bool write(std::ostream &os) const { return true; };

        virtual void setToOriginImpl() {
            _estimate.setZero();
        }

        virtual void oplusImpl(const double *update_);
    };

    /**
     * @brief EdgeNavStateLinePoint
     */
    class EdgeNavStateLinePoint : public BaseBinaryEdge<3, Vector3d, VertexLinePoint, VertexNavStatePVR> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeNavStateLinePoint() : BaseBinaryEdge<3, Vector3d, VertexLinePoint, VertexNavStatePVR>() {};

        bool read(std::istream &is) { return true; };

        bool write(std::ostream &os) const { return true; };

        void computeError() {
            Vector3d Pc = computePc();
            Vector3d obs(_measurement);

            Vector2d UV = cam_project(Pc);

            _error(0,0) = obs(0) * UV(0) + obs(1) * UV(1) + obs(2);
            _error(1,0) = 0;
            _error(2,0) = 0;
            // std::cout<<"error: "<<_error<<std::endl;
        }

        bool isDepthPositive() {
            Vector3d Pc = computePc();
            return Pc(2) > 0.0;
        }

        Vector3d computePc() {
            const VertexLinePoint *vPoint = static_cast<const VertexLinePoint *>(_vertices[0]);
            const VertexNavStatePVR *vNavState = static_cast<const VertexNavStatePVR *>(_vertices[1]);

            const NavState &ns = vNavState->estimate();
            Matrix3d Rwb = ns.Get_RotMatrix();
            Vector3d Pwb = ns.Get_P();
            const Vector3d &Pw = vPoint->estimate();

            Matrix3d Rcb = Rbc.transpose();
            Vector3d Pc = Rcb * Rwb.transpose() * (Pw - Pwb) - Rcb * Pbc;

            return Pc;
            //Vector3d Pwc = Rwb*Pbc + Pwb;
            //Matrix3d Rcw = (Rwb*Rbc).transpose();
            //Vector3d Pcw = -Rcw*Pwc;
            //Vector3d Pc = Rcw*Pw + Pcw;
        }

        inline Vector2d project2d(const Vector3d &v) const {
            Vector2d res;
            res(0) = v(0) / v(2);
            res(1) = v(1) / v(2);
            return res;
        }

        Vector2d cam_project(const Vector3d &trans_xyz) const {
            Vector2d proj = project2d(trans_xyz);
            Vector2d res;
            res[0] = proj[0] * fx + cx;
            res[1] = proj[1] * fy + cy;
            return res;
        }

        //
        virtual void linearizeOplus();

        void SetParams(const double &fx_, const double &fy_, const double &cx_, const double &cy_,
                       const Matrix3d &Rbc_, const Vector3d &Pbc_) {
            fx = fx_;
            fy = fy_;
            cx = cx_;
            cy = cy_;
            Rbc = Rbc_;
            Pbc = Pbc_;
        }

    protected:
        // Camera intrinsics
        double fx, fy, cx, cy;
        // Camera-IMU extrinsics
        Matrix3d Rbc;
        Vector3d Pbc;

    };
}

#endif // G2OTYPES_H
