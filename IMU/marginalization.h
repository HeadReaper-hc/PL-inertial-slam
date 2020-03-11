//
// Created by hc on 2020/3/5.
//

#ifndef PL_SLAM_MARGINALIZATION_H
#define PL_SLAM_MARGINALIZATION_H

#include <vector>
#include <map>

#include <keyFrame.h>
#include <mapFeatures.h>
#include <NavState.h>

#include <Eigen/Dense>

//g2o
#include <g2o/core/base_vertex.h>
#include <g2o/core/base_unary_edge.h>
#include <g2o/core/base_multi_edge.h>
#include <g2o/core/base_binary_edge.h>


using namespace g2o;
using namespace std;
using namespace Eigen;
using namespace PLSLAM;

struct ResidualBlockInfo
{
    ResidualBlockInfo(g2o::OptimizableGraph::Edge* _edge, std::vector<int> _drop_set, std::vector<VectorXd> _vertex_data,
                      std::vector<double*> _jac_addr, string _ResidualBlockType="")
        : edge(_edge), drop_set(_drop_set),  jac_addr(_jac_addr), vertex_data(_vertex_data),
          ResidualBlockType(_ResidualBlockType){
        for(int i=0;i<edge->vertices().size();i++){
            g2o::OptimizableGraph::Vertex* temp = dynamic_cast<g2o::OptimizableGraph::Vertex*>(edge->vertex(i));
            vertexs.push_back(temp);
            int N = temp->estimateDimension();
            if(N==-1){
                std::cout<<"Please check the overload function [estimateDimension()] in your vertex...."<<std::endl;
                exit(0);
            }
            vertex_local_size.push_back(N);
        }

    }

    void Evaluate();
     string ResidualBlockType;       //for debug
     g2o::OptimizableGraph::Edge* edge;
     std::vector<g2o::OptimizableGraph::Vertex*> vertexs;
     std::vector<int> vertex_local_size;
     std::vector<VectorXd> vertex_data;
     std::vector<int> drop_set;
     std::vector<double*> jac_addr;

     std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> jacobians;
     Eigen::VectorXd residuals;
   //  double** raw_jacobians;  //内存泄漏？
};

struct ThreadsStruct
{
    std::vector<ResidualBlockInfo *> sub_factors;
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    std::map<g2o::OptimizableGraph::Vertex*, int> parameter_block_size; //global size
    std::map<g2o::OptimizableGraph::Vertex*, int> parameter_block_idx; //local size
};

class MarginalizationInfo {
public:
    MarginalizationInfo(){}
    MarginalizationInfo(const MarginalizationInfo& marg_Info) = default;
    ~MarginalizationInfo(){};

    void preMarginalize();
    void marginalize();
    void marginalizeWithoutThread();
    int m,n;
    void addResidualBlockInfo(ResidualBlockInfo *residual_block_info);

    map<int, KeyFrame*> id_kf_map;
    //undefinition
    vector<int> keep_vertex_size; //local_size
    vector<int> keep_vertex_id;  //vertex idx
    vector<int> keep_vertex_idx;
    vector<VectorXd> keep_vertex_data;

    map<g2o::OptimizableGraph::Vertex*, VectorXd> vertex_addr_data;
    map<g2o::OptimizableGraph::Vertex*, int> vertex_addr_size; //all vertex and local size
    map<g2o::OptimizableGraph::Vertex*, int> vertex_addr_idx;  //all vertex and its idx in the hession matrix

    Eigen::MatrixXd linearized_jacobians;
    Eigen::VectorXd linearized_residuals;

    vector<ResidualBlockInfo*> factors;

    double eps = 1e-8;
};

#endif //PL_SLAM_MARGINALIZATION_H
