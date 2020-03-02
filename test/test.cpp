//
// Created by hc on 2020/2/27.
//

#include <g2otypes.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/solvers/eigen/linear_solver_eigen.h>
#include <g2o/core/robust_kernel_impl.h>
#include "g2otypes.h"
#include <IMUPreintegrator.h>

int main() {
    NavState imuState;
    Matrix3d Rbc = Matrix3d::Identity();
    Vector3d Pbc;
    Pbc.setZero();

    g2o::SparseOptimizer optimizer;

    auto linearSolver = g2o::make_unique<SlamLinearSolver>();
//    linearSolver->setBlockOrdering(false);
    auto blockSolver = g2o::make_unique<g2o::BlockSolverX>(std::move(linearSolver));
    g2o::OptimizationAlgorithm *algorithm = new g2o::OptimizationAlgorithmLevenberg(std::move(blockSolver));

    optimizer.setAlgorithm(algorithm);

    g2o::VertexNavState * vNSPVR = new g2o::VertexNavState();
    vNSPVR->setEstimate(imuState);
    vNSPVR->setId(0);
  //  vNSPVR->setFixed(true);
    optimizer.addVertex(vNSPVR);

    g2o::VertexLinePoint* vLinePoints = new g2o::VertexLinePoint();
    vLinePoints->setEstimate(Vector3d(5,5,5));
    vLinePoints->setId(1);
    vLinePoints->setFixed(false);
    vLinePoints->setMarginalized(true);
    optimizer.addVertex(vLinePoints);

    g2o::EdgeNavStateLinePoint* es = new g2o::EdgeNavStateLinePoint();
    es->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(1)));
    es->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));

    es->setMeasurement(Vector3d(1,2,3));
    es->setInformation(Eigen::Matrix<double,3,3>::Identity());
    es->SetParams(100,100,100,100,Rbc,Pbc);

    optimizer.addEdge(es);

    optimizer.setVerbose(true);
    optimizer.initializeOptimization();
    optimizer.optimize(5);

}