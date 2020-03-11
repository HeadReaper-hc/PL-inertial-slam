//
// Created by hc on 2020/3/5.
//

#include "marginalization.h"
#include <NavState.h>
#include <chrono>
    void* ThreadsConstructA(void* threadsstruct)
    {
        ThreadsStruct* p = ((ThreadsStruct*)threadsstruct);
        for (auto it : p->sub_factors)
        {
            for (int i = 0; i < static_cast<int>(it->vertex_local_size.size()); i++)
            {
                int idx_i = p->parameter_block_idx[it->vertexs[i]];
                int size_i = p->parameter_block_size[it->vertexs[i]];
                Eigen::MatrixXd jacobian_i = it->jacobians[i].leftCols(size_i);
                for (int j = i; j < static_cast<int>(it->vertex_local_size.size()); j++)
                {
                    int idx_j = p->parameter_block_idx[it->vertexs[j]];
                    int size_j = p->parameter_block_size[it->vertexs[j]];

                    Eigen::MatrixXd jacobian_j = it->jacobians[j].leftCols(size_j);
                    if (i == j)
                        p->A.block(idx_i, idx_j, size_i, size_j) += jacobian_i.transpose() * jacobian_j;
                    else
                    {
                        p->A.block(idx_i, idx_j, size_i, size_j) += jacobian_i.transpose() * jacobian_j;
                        p->A.block(idx_j, idx_i, size_j, size_i) = p->A.block(idx_i, idx_j, size_i, size_j).transpose();
                    }
                }
                p->b.segment(idx_i, size_i) += jacobian_i.transpose() * it->residuals;
            }
        }
        return threadsstruct;
    }

    void ResidualBlockInfo::Evaluate()
    {
       // std::cout<<"ResidualBlockType: "<<ResidualBlockType<<std::endl;
        residuals.resize(edge->dimension());
        //raw_jacobians = new double *[vertex_local_size.size()];
        jacobians.resize(vertexs.size());
        for (int i = 0; i < vertexs.size(); i++)
        {
            int vertexDimension = vertexs[i]->estimateDimension();
            if(vertexDimension==-1){
                std::cout<<"Please check the overload function [estimateDimension()] in your vertex...."<<std::endl;
                exit(0);
            }
            jacobians[i].resize(edge->dimension(), vertexDimension);
            //raw_jacobians[i] = jacobians[i].data();
            //dim += block_sizes[i] == 7 ? 6 : block_sizes[i];
        }
        edge->computeError();
        residuals = Eigen::Map<VectorXd>(edge->errorData(),edge->dimension());
        for (int i = 0; i< vertexs.size(); i++)
        {
            int D = edge->dimension();
            int vertexDimension = vertexs[i]->estimateDimension();
            if(vertexDimension==-1){
                std::cout<<"Please check the overload function [estimateDimension()] in your vertex...."<<std::endl;
                exit(0);
            }
            jacobians[i] = Eigen::Map<MatrixXd>(jac_addr[i],D,vertexDimension);
        }
        //still not consider information and huber kernel
//        if (loss_function)
//        {
//            double residual_scaling_, alpha_sq_norm_;
//
//            double sq_norm, rho[3];
//
//            sq_norm = residuals.squaredNorm();
//            loss_function->Evaluate(sq_norm, rho);
//            //printf("sq_norm: %f, rho[0]: %f, rho[1]: %f, rho[2]: %f\n", sq_norm, rho[0], rho[1], rho[2]);
//
//            double sqrt_rho1_ = sqrt(rho[1]);
//
//            if ((sq_norm == 0.0) || (rho[2] <= 0.0))
//            {
//                residual_scaling_ = sqrt_rho1_;
//                alpha_sq_norm_ = 0.0;
//            }
//            else
//            {
//                const double D = 1.0 + 2.0 * sq_norm * rho[2] / rho[1];
//                const double alpha = 1.0 - sqrt(D);
//                residual_scaling_ = sqrt_rho1_ / (1 - alpha);
//                alpha_sq_norm_ = alpha / sq_norm;
//            }
//
//            for (int i = 0; i < static_cast<int>(parameter_blocks.size()); i++)
//            {
//                jacobians[i] = sqrt_rho1_ * (jacobians[i] - alpha_sq_norm_ * residuals * (residuals.transpose() * jacobians[i]));
//            }
//
//            residuals *= residual_scaling_;
//        }
    }

    void MarginalizationInfo::addResidualBlockInfo(ResidualBlockInfo *residual_block_info)
    {
        factors.emplace_back(residual_block_info);
        std::vector<int> parameter_block_sizes = residual_block_info->vertex_local_size;
        std::vector<g2o::OptimizableGraph::Vertex*> vertexs = residual_block_info->vertexs;
        std::vector<VectorXd> parameter_data = residual_block_info->vertex_data;
        std::vector<int> drop_set = residual_block_info->drop_set;
        for (int i = 0; i < vertexs.size(); i++)
        {
            map<g2o::OptimizableGraph::Vertex*,int>::iterator iter = vertex_addr_size.begin();
            iter = vertex_addr_size.find(vertexs[i]);
            if(iter!=vertex_addr_size.end() && iter->second!=parameter_block_sizes[i]){
                std::cout<<"Wrong! The same param block with different param size....."<<std::endl;
                exit(0);
            }
            vertex_addr_size[vertexs[i]] = parameter_block_sizes[i];
            vertex_addr_data[vertexs[i]] = parameter_data[i];
        }

        for (int i = 0; i < drop_set.size(); i++)
        {
            g2o::OptimizableGraph::Vertex* vertex2drop = vertexs[drop_set[i]];
            vertex_addr_idx[vertex2drop] = 0;
        }
    }

    void MarginalizationInfo::preMarginalize()
    {
        for (auto it : factors)
        {
            it->Evaluate();

//            std::vector<int> block_sizes = it->vertex_local_size;
//            for (int i = 0; i < static_cast<int>(block_sizes.size()); i++)
//            {
//                g2o::OptimizableGraph::Vertex* temp = it->vertexs[i];
//                int size = block_sizes[i];
//                if (parameter_block_data.find(addr) == parameter_block_data.end())
//                {
//                    double *data = new double[size];
//                    memcpy(data, it->parameter_blocks[i], sizeof(double) * size);
//                    parameter_block_data[addr] = data;
//                }
//            }
        }
    }

    void MarginalizationInfo::marginalize()
    {
        chrono::steady_clock::time_point time_start=chrono::steady_clock::now();
        keep_vertex_size.clear();
        keep_vertex_id.clear();
        keep_vertex_data.clear();
        keep_vertex_idx.clear();
        int pos = 0;
        for (auto &it : vertex_addr_idx)
        {
            it.second = pos;
            pos += vertex_addr_size[it.first];

        }

        m = pos;

        for (const auto &it : vertex_addr_size)
        {
            if (vertex_addr_idx.find(it.first) == vertex_addr_idx.end())
            {
                vertex_addr_idx[it.first] = pos;
                pos += it.second;

                keep_vertex_size.push_back(vertex_addr_size[it.first]);
                keep_vertex_id.push_back( (it.first)->id() );  //only suit for PVR BIAS vertex,because landmark id add a maxKFid.
                keep_vertex_idx.push_back(vertex_addr_idx[it.first]);
                int N = vertex_addr_size[it.first];
                if(N==9){
                    N = 10; //for PVR vertex else 3
                }
                if( vertex_addr_data[it.first].size() != N ){
                    std::cout<<"Wrong size between Vertex size and Vertex data size...."<<std::endl;
                    exit(0);
                }
                keep_vertex_data.push_back( vertex_addr_data[it.first] );
            }
        }

        n = pos - m;

        //ROS_DEBUG("marginalization, pos: %d, m: %d, n: %d, size: %d", pos, m, n, (int)parameter_block_idx.size());

        Eigen::MatrixXd A(pos, pos);
        Eigen::VectorXd b(pos);
        A.setZero();
        b.setZero();
        /*
        for (auto it : factors)
        {
            for (int i = 0; i < static_cast<int>(it->parameter_blocks.size()); i++)
            {
                int idx_i = parameter_block_idx[reinterpret_cast<long>(it->parameter_blocks[i])];
                int size_i = localSize(parameter_block_size[reinterpret_cast<long>(it->parameter_blocks[i])]);
                Eigen::MatrixXd jacobian_i = it->jacobians[i].leftCols(size_i);
                for (int j = i; j < static_cast<int>(it->parameter_blocks.size()); j++)
                {
                    int idx_j = parameter_block_idx[reinterpret_cast<long>(it->parameter_blocks[j])];
                    int size_j = localSize(parameter_block_size[reinterpret_cast<long>(it->parameter_blocks[j])]);
                    Eigen::MatrixXd jacobian_j = it->jacobians[j].leftCols(size_j);
                    if (i == j)
                        A.block(idx_i, idx_j, size_i, size_j) += jacobian_i.transpose() * jacobian_j;
                    else
                    {
                        A.block(idx_i, idx_j, size_i, size_j) += jacobian_i.transpose() * jacobian_j;
                        A.block(idx_j, idx_i, size_j, size_i) = A.block(idx_i, idx_j, size_i, size_j).transpose();
                    }
                }
                b.segment(idx_i, size_i) += jacobian_i.transpose() * it->residuals;
            }
        }
        ROS_INFO("summing up costs %f ms", t_summing.toc());
        */
        //multi thread
        int NUM_THREADS = 4;

        pthread_t tids[NUM_THREADS];
        ThreadsStruct threadsstruct[NUM_THREADS];
        int i = 0;
        for (auto it : factors)
        {
            threadsstruct[i].sub_factors.push_back(it);
            i++;
            i = i % NUM_THREADS;
        }
        for (int i = 0; i < NUM_THREADS; i++)
        {
            threadsstruct[i].A = Eigen::MatrixXd::Zero(pos,pos);
            threadsstruct[i].b = Eigen::VectorXd::Zero(pos);
            threadsstruct[i].parameter_block_size = vertex_addr_size;
            threadsstruct[i].parameter_block_idx = vertex_addr_idx;
            int ret = pthread_create( &tids[i], NULL, ThreadsConstructA ,(void*)&(threadsstruct[i]));
            if (ret != 0)
            {
                cout<<"pthread_create error"<<endl;
                exit(0);
            }
        }
        for( int i = NUM_THREADS - 1; i >= 0; i--)
        {
            pthread_join( tids[i], NULL );
            A += threadsstruct[i].A;
            b += threadsstruct[i].b;
        }
        chrono::steady_clock::time_point time_mid=chrono::steady_clock::now();
        //TODO
        Eigen::MatrixXd Amm = 0.5 * (A.block(0, 0, m, m) + A.block(0, 0, m, m).transpose());
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes(Amm);
        Eigen::MatrixXd Amm_inv = saes.eigenvectors() * Eigen::VectorXd((saes.eigenvalues().array() > eps).select(saes.eigenvalues().array().inverse(), 0)).asDiagonal() * saes.eigenvectors().transpose();
        //printf("error1: %f\n", (Amm * Amm_inv - Eigen::MatrixXd::Identity(m, m)).sum());

        Eigen::VectorXd bmm = b.segment(0, m);
        Eigen::MatrixXd Amr = A.block(0, m, m, n);
        Eigen::MatrixXd Arm = A.block(m, 0, n, m);
        Eigen::MatrixXd Arr = A.block(m, m, n, n);
        Eigen::VectorXd brr = b.segment(m, n);
        A = Arr - Arm * Amm_inv * Amr;
        b = brr - Arm * Amm_inv * bmm;

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes2(A);
        Eigen::VectorXd S = Eigen::VectorXd((saes2.eigenvalues().array() > eps).select(saes2.eigenvalues().array(), 0));
        Eigen::VectorXd S_inv = Eigen::VectorXd((saes2.eigenvalues().array() > eps).select(saes2.eigenvalues().array().inverse(), 0));

        Eigen::VectorXd S_sqrt = S.cwiseSqrt();
        Eigen::VectorXd S_inv_sqrt = S_inv.cwiseSqrt();

        linearized_jacobians = S_sqrt.asDiagonal() * saes2.eigenvectors().transpose();
        linearized_residuals = S_inv_sqrt.asDiagonal() * saes2.eigenvectors().transpose() * b;
        chrono::steady_clock::time_point time_end=chrono::steady_clock::now();

        chrono::duration<double> time_mid_used=chrono::duration_cast<chrono::duration<double>>(time_mid-time_start);
        chrono::duration<double> time_end_used=chrono::duration_cast<chrono::duration<double>>(time_end-time_mid);

        cout<<"Marginalization thread construct time used:"<<time_mid_used.count()<<"s"<<endl;
        cout<<"Marginalization matrix compute time used:"<<time_end_used.count()<<"s"<<endl;
        //std::cout << A << std::endl
        //          << std::endl;
        //std::cout << linearized_jacobians << std::endl;
        //printf("error2: %f %f\n", (linearized_jacobians.transpose() * linearized_jacobians - A).sum(),
        //      (linearized_jacobians.transpose() * linearized_residuals - b).sum());
    }

    void MarginalizationInfo::marginalizeWithoutThread() {

        chrono::steady_clock::time_point time_start=chrono::steady_clock::now();
        keep_vertex_size.clear();
        keep_vertex_id.clear();
        keep_vertex_data.clear();
        keep_vertex_idx.clear();
        int pos = 0;
        for (auto &it : vertex_addr_idx)
        {
            it.second = pos;
            pos += vertex_addr_size[it.first];

        }

        m = pos;

        for (const auto &it : vertex_addr_size)
        {
            if (vertex_addr_idx.find(it.first) == vertex_addr_idx.end())
            {
                vertex_addr_idx[it.first] = pos;
                pos += it.second;

                keep_vertex_size.push_back(vertex_addr_size[it.first]);
                keep_vertex_id.push_back( (it.first)->id() );  //only suit for PVR BIAS vertex,because landmark id add a maxKFid.
                keep_vertex_idx.push_back(vertex_addr_idx[it.first]);
                int N = vertex_addr_size[it.first];
                if(N==9){
                    N = 10; //for PVR vertex else 3
                }
                if( vertex_addr_data[it.first].size() != N ){
                    std::cout<<"Wrong size between Vertex size and Vertex data size...."<<std::endl;
                    exit(0);
                }
                keep_vertex_data.push_back( vertex_addr_data[it.first] );
            }
        }

        n = pos - m;

        //ROS_DEBUG("marginalization, pos: %d, m: %d, n: %d, size: %d", pos, m, n, (int)parameter_block_idx.size());

        Eigen::MatrixXd A(pos, pos);
        Eigen::VectorXd b(pos);
        A.setZero();
        b.setZero();

        ThreadsStruct threadsstruct;
        threadsstruct.sub_factors = factors;
        threadsstruct.A = Eigen::MatrixXd::Zero(pos,pos);
        threadsstruct.b = Eigen::VectorXd::Zero(pos);
        threadsstruct.parameter_block_size = vertex_addr_size;
        threadsstruct.parameter_block_idx = vertex_addr_idx;
        ThreadsConstructA(&threadsstruct);

        A = threadsstruct.A;
        b = threadsstruct.b;
        chrono::steady_clock::time_point time_mid=chrono::steady_clock::now();
        //TODO
        Eigen::MatrixXd Amm = 0.5 * (A.block(0, 0, m, m) + A.block(0, 0, m, m).transpose());
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes(Amm);
        Eigen::MatrixXd Amm_inv = saes.eigenvectors() * Eigen::VectorXd((saes.eigenvalues().array() > eps).select(saes.eigenvalues().array().inverse(), 0)).asDiagonal() * saes.eigenvectors().transpose();
        //printf("error1: %f\n", (Amm * Amm_inv - Eigen::MatrixXd::Identity(m, m)).sum());

        Eigen::VectorXd bmm = b.segment(0, m);
        Eigen::MatrixXd Amr = A.block(0, m, m, n);
        Eigen::MatrixXd Arm = A.block(m, 0, n, m);
        Eigen::MatrixXd Arr = A.block(m, m, n, n);
        Eigen::VectorXd brr = b.segment(m, n);
        A = Arr - Arm * Amm_inv * Amr;
        b = brr - Arm * Amm_inv * bmm;

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes2(A);
        Eigen::VectorXd S = Eigen::VectorXd((saes2.eigenvalues().array() > eps).select(saes2.eigenvalues().array(), 0));
        Eigen::VectorXd S_inv = Eigen::VectorXd((saes2.eigenvalues().array() > eps).select(saes2.eigenvalues().array().inverse(), 0));

        Eigen::VectorXd S_sqrt = S.cwiseSqrt();
        Eigen::VectorXd S_inv_sqrt = S_inv.cwiseSqrt();

        linearized_jacobians = S_sqrt.asDiagonal() * saes2.eigenvectors().transpose();
        linearized_residuals = S_inv_sqrt.asDiagonal() * saes2.eigenvectors().transpose() * b;
        chrono::steady_clock::time_point time_end=chrono::steady_clock::now();

        chrono::duration<double> time_mid_used=chrono::duration_cast<chrono::duration<double>>(time_mid-time_start);
        chrono::duration<double> time_end_used=chrono::duration_cast<chrono::duration<double>>(time_end-time_mid);

        cout<<"Marginalization hession construct time used:"<<time_mid_used.count()<<"s"<<endl;
        cout<<"Marginalization hession decompose time used:"<<time_end_used.count()<<"s"<<endl;
        //std::cout << A << std::endl
        //          << std::endl;
        //std::cout << linearized_jacobians << std::endl;
        //printf("error2: %f %f\n", (linearized_jacobians.transpose() * linearized_jacobians - A).sum(),
        //      (linearized_jacobians.transpose() * linearized_residuals - b).sum());
    }

