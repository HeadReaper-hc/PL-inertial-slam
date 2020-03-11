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

#include "marginalization.h"
//#include "Thirdparty/g2o/g2o/core/sparse_block_matrix.h"
typedef g2o::LinearSolverEigen<g2o::BlockSolverX::PoseMatrixType> SlamLinearSolver;
//#include "Thirdparty/g2o/g2o/core/sparse_block_matrix.h"

namespace g2o {

    class VertexLMPointXYZ : public BaseVertex<3, Vector3d>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        VertexLMPointXYZ();
        virtual bool read(std::istream& is);
        virtual bool write(std::ostream& os) const;

        virtual void setToOriginImpl() {
            _estimate.fill(0.);
        }

        virtual void oplusImpl(const double* update)
        {
            Eigen::Map<const Vector3d> v(update);
            _estimate += v;
        }

        virtual int estimateDimension() const{
            return 3;
        }
    };


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

        virtual int estimateDimension() const{
            return 9;
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

        virtual int estimateDimension() const{
            return 6;
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

        void GetJacAddr(std::vector<double*>& addr){
            addr.clear();
            addr.push_back(_jacobianOplus[0].data());
            addr.push_back(_jacobianOplus[1].data());
            addr.push_back(_jacobianOplus[2].data());
        }

        void GetEstData(std::vector<VectorXd>& data){
            const VertexNavStatePVR* vPVRi = static_cast<const VertexNavStatePVR*>(_vertices[0]);
            const VertexNavStatePVR* vPVRj = static_cast<const VertexNavStatePVR*>(_vertices[1]);
            const VertexNavStateBias* vBiasi = static_cast<const VertexNavStateBias*>(_vertices[2]);

            const NavState& NSPVRi = vPVRi->estimate();
            Vector3d Pi = NSPVRi.Get_P();
            Vector3d Vi = NSPVRi.Get_V();
            Matrix3d Ri = NSPVRi.Get_RotMatrix();

            const NavState& NSBiasi = vBiasi->estimate();
            Vector3d dBgi = NSBiasi.Get_dBias_Gyr();
            Vector3d dBai = NSBiasi.Get_dBias_Acc();

            const NavState& NSPVRj = vPVRj->estimate();
            Vector3d Pj = NSPVRj.Get_P();
            Vector3d Vj = NSPVRj.Get_V();
            Matrix3d Rj = NSPVRj.Get_RotMatrix();

            VectorXd temp_vpi,temp_vpj,temp_bi;
            temp_vpi.resize(10);
            temp_vpj.resize(10);
            temp_bi.resize(6);

            temp_vpi.segment<3>(0) = Pi;
            temp_vpi.segment<3>(3) = Vi;
            Quaterniond quai(Ri);
            temp_vpi.segment<4>(6) = Vector4d(quai.x(),quai.y(),quai.z(),quai.w());

            temp_vpj.segment<3>(0) = Pj;
            temp_vpj.segment<3>(3) = Vj;
            Quaterniond quaj(Rj);
            temp_vpj.segment<4>(6) = Vector4d(quaj.x(),quaj.y(),quaj.z(),quaj.w());

            temp_bi.segment<3>(0) = dBgi + NSBiasi.Get_BiasGyr();
            temp_bi.segment<3>(3) = dBai + NSBiasi.Get_BiasAcc();

            data.push_back(temp_vpi);
            data.push_back(temp_vpj);
            data.push_back(temp_bi);

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

        void GetJacAddr(std::vector<double*>& addr){
            addr.clear();
            addr.push_back(_jacobianOplusXi.data());
            addr.push_back(_jacobianOplusXj.data());
        }

        void GetEstData(vector<VectorXd>& data){
            VectorXd temp_a,temp_b;
            temp_a.resize(6);
            temp_b.resize(6);

            const VertexNavStateBias* vBiasi = static_cast<const VertexNavStateBias*>(_vertices[0]);
            const VertexNavStateBias* vBiasj = static_cast<const VertexNavStateBias*>(_vertices[1]);

            const NavState& NSi = vBiasi->estimate();
            const NavState& NSj = vBiasj->estimate();

            temp_a.segment<3>(0) = NSi.Get_BiasGyr() + NSi.Get_dBias_Gyr();
            temp_a.segment<3>(3) = NSi.Get_BiasAcc() + NSi.Get_dBias_Acc();

            temp_b.segment<3>(0) = NSj.Get_BiasGyr() + NSj.Get_dBias_Gyr();
            temp_b.segment<3>(3) = NSj.Get_BiasAcc() + NSj.Get_dBias_Acc();

            data.push_back(temp_a);
            data.push_back(temp_b);


        }

    };

    class EdgeNavStatePVRPointXYZ : public BaseBinaryEdge<2, Vector2d, VertexLMPointXYZ, VertexNavStatePVR> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeNavStatePVRPointXYZ() : BaseBinaryEdge<2, Vector2d, VertexLMPointXYZ, VertexNavStatePVR>() {}

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
            const VertexLMPointXYZ *vPoint = static_cast<const VertexLMPointXYZ *>(_vertices[0]);
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

        void GetJacAddr(std::vector<double*>& addr){
            addr.clear();
            addr.push_back(_jacobianOplusXi.data());
            addr.push_back(_jacobianOplusXj.data());
        }

        void GetEstData(std::vector<VectorXd>& data){
            VectorXd temp_a,temp_b;
            temp_a.resize(3);
            temp_b.resize(10);

            const VertexLMPointXYZ* vPoint = static_cast<const VertexLMPointXYZ*>(_vertices[0]);
            const VertexNavStatePVR* vNavState = static_cast<const VertexNavStatePVR*>(_vertices[1]);

            const NavState& ns = vNavState->estimate();
            Matrix3d Rwb = ns.Get_RotMatrix();
            Vector3d Pwb = ns.Get_P();
            Vector3d Pwv = ns.Get_V();
            const Vector3d& Pw = vPoint->estimate();

            temp_a.segment<3>(0) = Pw;

            temp_b.segment<3>(0) = Pwb;
            temp_b.segment<3>(3) = Pwv;
            Quaterniond qua(Rwb);
            temp_b.segment<4>(6) = Vector4d(qua.x(),qua.y(),qua.z(),qua.w());

            data.push_back(temp_a);
            data.push_back(temp_b);
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
    class EdgeNavStatePointXYZ : public BaseBinaryEdge<2, Vector2d, VertexLMPointXYZ, VertexNavState> {
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
            const VertexLMPointXYZ *vPoint = static_cast<const VertexLMPointXYZ *>(_vertices[0]);
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

        virtual int estimateDimension() const{

            return 6;
        }
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
            const VertexLine* vPoint = static_cast<const VertexLine*>(_vertices[0]);
            const VertexNavStatePVR* vNavState = static_cast<const VertexNavStatePVR*>(_vertices[1]);

            const NavState& ns = vNavState->estimate();
            Matrix3d Rwb = ns.Get_RotMatrix();
            Vector3d Pwb = ns.Get_P();
            const Vector6d& Pw = vPoint->estimate();

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

        void GetJacAddr(std::vector<double*>& addr){
            addr.clear();
            addr.push_back(_jacobianOplusXi.data());
            addr.push_back(_jacobianOplusXj.data());
        }

        void GetEstData(std::vector<VectorXd>& data){

            VectorXd temp_a, temp_b;
            temp_a.resize(6);
            temp_b.resize(10);

            const VertexLine* vPoint = static_cast<const VertexLine*>(_vertices[0]);
            const VertexNavStatePVR* vNavState = static_cast<const VertexNavStatePVR*>(_vertices[1]);

            const NavState& ns = vNavState->estimate();
            Matrix3d Rwb = ns.Get_RotMatrix();
            Vector3d Pwb = ns.Get_P();
            Vector3d Pwv = ns.Get_V();
            const Vector6d& Pw = vPoint->estimate();

            temp_a.segment<6>(0) = Pw;

            temp_b.segment<3>(0) = Pwb;
            temp_b.segment<3>(3) = Pwv;
            Quaterniond qua(Rwb);
            temp_b.segment<4>(6) = Vector4d(qua.x(),qua.y(),qua.z(),qua.w());

            data.push_back(temp_a);
            data.push_back(temp_b);
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


    class EdgeMarginalization : public BaseMultiEdge<-1, MarginalizationInfo> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeMarginalization() : BaseMultiEdge<-1, MarginalizationInfo>() {
        }

        bool read(std::istream &is) { return true; }

        bool write(std::ostream &os) const { return true; }

        void computeError();

        virtual void linearizeOplus();

        void setDimension(int dimension_)
        {
            _dimension = dimension_;
            _information.resize(dimension_, dimension_);
            _error.resize(dimension_, 1);
        }

        void setSize(int vertices)
        {
            resize(vertices);
        }

        void GetJacAddr(std::vector<double*>& addr){
            addr.clear();
            for(int i=0;i<_vertices.size();i++){
                addr.push_back(_jacobianOplus[i].data());
            }
        }

        void GetEstData(std::vector<VectorXd>& data){
            MarginalizationInfo* marg = &(_measurement);
            if(_vertices.size()!=marg->keep_vertex_id.size()){
                std::cout<<"Wrong size between vertex and variable in MarginalizationInfo...."<<std::endl;
                exit(0);
            }
            for(int i=0;i<_vertices.size();i++){
                int N = marg->keep_vertex_size[i];
                if(N==9){
                    VectorXd temp;
                    temp.resize(10);
                    VertexNavStatePVR* pvr = dynamic_cast<VertexNavStatePVR*>(_vertices[i]);
                    const NavState& pvrState = pvr->estimate();
                    Matrix3d Rwb = pvrState.Get_RotMatrix();
                    Vector3d Pwb = pvrState.Get_P();
                    Vector3d Vwb = pvrState.Get_V();

                    temp.segment<3>(0) = Pwb;
                    temp.segment<3>(3) = Vwb;
                    Quaterniond qua(Rwb);
                    temp.segment<4>(6) = Vector4d(qua.x(),qua.y(),qua.z(),qua.w());
                    data.push_back(temp);
                }
                else if(N==6){
                    VectorXd temp;
                    temp.resize(6);
                    VertexNavStateBias* bias = dynamic_cast<VertexNavStateBias*>(_vertices[i]);
                    const NavState& biasState = bias->estimate();
                    Vector3d bias_gyr = biasState.Get_BiasGyr() + biasState.Get_dBias_Gyr();
                    Vector3d bias_acc = biasState.Get_BiasAcc() + biasState.Get_dBias_Acc();

                    temp.segment<3>(0) = bias_gyr;
                    temp.segment<3>(3) = bias_acc;
                    data.push_back(temp);
                }
                else{
                    std::cout<<"Undefined size of marginalization vertex: "<<N<<std::endl;
                    exit(0);
                }
            }
        }
    };
}

#endif // G2OTYPES_H
