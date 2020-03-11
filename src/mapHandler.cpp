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

#include "mapHandler.h"

#include <iostream>
#include <map>
#include <opencv2/imgproc.hpp>

#include <matching.h>
#include <timer.h>

//g2o
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/solvers/eigen/linear_solver_eigen.h>
#include <g2o/core/robust_kernel_impl.h>
#include "g2otypes.h"
#include <IMUPreintegrator.h>
#include <types_sba.h>

namespace PLSLAM
{

MapHandler::MapHandler(PinholeStereoCamera* cam_)
    : cam(cam_), threads_started(false)
{
    // load vocabulary
    if( SlamConfig::hasPoints() )
        dbow_voc_p.load( SlamConfig::dbowVocP() );
    if( SlamConfig::hasLines() )
        dbow_voc_l.load( SlamConfig::dbowVocL() );

    lc_state = LC_IDLE;

}

void MapHandler::initialize( KeyFrame *kf0 )
{
    if(datasetPtr != NULL && datasetPtr->is_have_gt()){
        SetGt(true);
    }

    if(have_gt){
        kf0->T_kf_w_gt = datasetPtr->getGroundTruthByTime(kf0->t_);
    }

    curr_kf = kf0;

    Twf = Matrix4d::Identity();
    DT = Matrix4d::Identity();

    //update first keyframe NavState
    Matrix4d tbw = curr_kf->T_kf_w * cam->getTbs().inverse();
//    cout<<"tbw: "<<tbw<<endl;
    curr_kf->imuState.Set_Pos(tbw.block<3,1>(0,3));
    curr_kf->imuState.Set_Rot(tbw.block<3,3>(0,0));

    // reset information from the map
    map_keyframes.clear();
    map_points.clear();
    map_points_kf_idx.clear();
    map_lines.clear();
    map_lines_kf_idx.clear();
    full_graph.clear();
    conf_matrix.clear();
    lc_idx_list.clear();
    lc_pose_list.clear();
    max_pt_idx = 0;
    max_ls_idx = 0;
    max_kf_idx = 0;

    // initialize graphs
    full_graph.resize(1);
    full_graph[0].resize(1);
    full_graph[0][0] = 0;
    conf_matrix.resize(1);
    conf_matrix[0].resize(1);
    conf_matrix[0][0] = 1.0;

    // reset indices
    for (PointFeature* pt : kf0->stereo_frame->stereo_pt)
        pt->idx = -1;
    for (LineFeature* ls : kf0->stereo_frame->stereo_ls)
        ls->idx = -1;

    // initialize DBoW descriptor vector and LC status
    vector<Mat> curr_desc;
    if( SlamConfig::hasPoints() )
    {
        curr_desc.reserve( kf0->stereo_frame->pdesc_l.rows );
        for ( int i = 0; i < kf0->stereo_frame->pdesc_l.rows; i++ )
            curr_desc.push_back( kf0->stereo_frame->pdesc_l.row(i) );
        dbow_voc_p.transform( curr_desc, kf0->descDBoW_P );
        curr_desc.clear();
    }
    if( SlamConfig::hasLines() )
    {
        curr_desc.reserve( kf0->stereo_frame->ldesc_l.rows );
        for ( int i = 0; i < kf0->stereo_frame->ldesc_l.rows; i++ )
            curr_desc.push_back( kf0->stereo_frame->ldesc_l.row(i) );
        dbow_voc_l.transform( curr_desc, kf0->descDBoW_L );
    }

    // insert keyframe and add to map of indexes
    vector<int> aux_vec;
    map_keyframes.push_back( kf0 );
    map_points_kf_idx.insert( std::pair<int,vector<int>>(kf0->kf_idx,aux_vec) );
    map_lines_kf_idx.insert(  std::pair<int,vector<int>>(kf0->kf_idx,aux_vec) );
    addKeyframeToSW(kf0);

    if (SlamConfig::multithreadSLAM())
        startThreads();

    time = Vector7f::Zero();
}

void MapHandler::finishSLAM()
{
    if (SlamConfig::multithreadSLAM())
        killThreads();
}

void MapHandler::addKeyFrame( KeyFrame *curr_kf )
{
    Timer timer;

    this->prev_kf = this->curr_kf;
    this->curr_kf = curr_kf;

    if(have_gt){
        curr_kf->T_kf_w_gt = datasetPtr->getGroundTruthByTime(curr_kf->t_);
    }

//    std::cout<<"prev_kf time: "<<prev_kf->stereo_frame->t_<<std::endl;
//    std::cout<<"curr_kf time: "<<curr_kf->stereo_frame->t_<<std::endl;
    std::vector<IMUData> imudata = datasetPtr->getIMUDataBetweenTwoImage(prev_kf->stereo_frame->t_,
                                                                      curr_kf->stereo_frame->t_);

    //test
//    vector<IMUData>::iterator imu_iter = imudata.begin();
//    for(int i=0;i<imudata.size();imu_iter++,i++){
//        IMUData temp = *imu_iter;
//        std::cout<<temp<<std::endl;
//    }
//    cout<<"[imus size] :"<<imudata.size()<<" "<<"[image time] : "<<curr_kf->stereo_frame->t_<<endl;
    curr_kf->imus = imudata;
    curr_kf->ComputeIMUPreIntSinceLastFrame(prev_kf);

    if( SlamConfig::multithreadSLAM() )
    {
        // expand graphs
        expandGraphs();
        // select previous keyframe
        KeyFrame* prev_kf;
        prev_kf = map_keyframes.back();
        max_kf_idx++;
        curr_kf->kf_idx = max_kf_idx;
        curr_kf->local  = true;
        // update pose of current keyframe wrt previous one (in case of LC)
        Matrix4d T_curr_w = prev_kf->T_kf_w * curr_kf->T_kf_w;
        curr_kf->x_kf_w = logmap_se3(T_curr_w);
        curr_kf->T_kf_w = expmap_se3(curr_kf->x_kf_w);
        // set current keyframe NavState
        Matrix4d tbw = curr_kf->T_kf_w * cam->getTbs().inverse();
//      cout<<"tbw: "<<tbw<<endl;
        curr_kf->imuState.Set_Pos(tbw.block<3,1>(0,3));
        curr_kf->imuState.Set_Rot(tbw.block<3,3>(0,0));
        // Estimates Twf
        Twf = expmap_se3(logmap_se3(  inverse_se3( curr_kf->T_kf_w ) ));
        // estimates pose increment
        DT = expmap_se3(logmap_se3( Twf * prev_kf->T_kf_w ));
        // reset indices
        for (PointFeature* pt : curr_kf->stereo_frame->stereo_pt)
            pt->idx = -1;
        for (LineFeature* ls : curr_kf->stereo_frame->stereo_ls)
            ls->idx = -1;
        // insert keyframe and add to map of indexes
        vector<int> aux_vec;
        map_keyframes.push_back( curr_kf );
        map_points_kf_idx.insert( std::make_pair(curr_kf->kf_idx, aux_vec) );
        map_lines_kf_idx.insert(  std::make_pair(curr_kf->kf_idx, aux_vec) );
        // call
        InterruptBA();
        addKeyFrame_multiThread(curr_kf,prev_kf);
        return;
    }


    // reset time variable
    time = Vector7f::Zero();

    // expand graphs
    timer.start();
    expandGraphs();
    time(0) = timer.stop(); //ms

    // select previous keyframe
    KeyFrame* prev_kf;
    prev_kf = map_keyframes.back();
    max_kf_idx++;
    curr_kf->kf_idx = max_kf_idx;
    curr_kf->local  = true;

    // update pose of current keyframe wrt previous one (in case of LC)
    Matrix4d T_curr_w = prev_kf->T_kf_w * curr_kf->T_kf_w;
    curr_kf->x_kf_w = logmap_se3(T_curr_w);
    curr_kf->T_kf_w = expmap_se3(curr_kf->x_kf_w);

    // Estimates Twf
    Twf = expmap_se3(logmap_se3(  inverse_se3( curr_kf->T_kf_w ) ));

    // estimates pose increment
    DT = expmap_se3(logmap_se3( Twf * prev_kf->T_kf_w ));

    // reset indices
    for (PointFeature* pt : curr_kf->stereo_frame->stereo_pt)
        pt->idx = -1;
    for (LineFeature* ls : curr_kf->stereo_frame->stereo_ls)
        ls->idx = -1;

    // look for common matches and update the full graph
    timer.start();
    lookForCommonMatches( prev_kf, curr_kf );
    time(1) = timer.stop(); //ms

    timer.start();
    if( SlamConfig::hasPoints() && SlamConfig::hasLines() )
        insertKFBowVectorPL(curr_kf);
    else if( SlamConfig::hasPoints() )
        insertKFBowVectorP(curr_kf);
    else if( SlamConfig::hasLines() )
        insertKFBowVectorL(curr_kf);

    time(2) = timer.stop(); //ms

    // insert keyframe and add to map of indexes
    vector<int> aux_vec;
    map_keyframes.push_back( curr_kf );
    map_points_kf_idx.insert( std::make_pair(curr_kf->kf_idx, aux_vec) );
    map_lines_kf_idx.insert(  std::make_pair(curr_kf->kf_idx, aux_vec) );

    // form local map
    timer.start();
    formLocalMap();
    time(3) = timer.stop(); //ms

    // perform local bundle adjustment
    timer.start();
    localBundleAdjustment();
    time(4) = timer.stop(); //ms

    // Recent map LMs culling (implement filters for line segments, which seems to be unaccurate)
    timer.start();
    removeBadMapLandmarks();
    time(5) = timer.stop(); //ms


    // LC
    timer.start();
    loopClosure();
    time(6) = timer.stop(); //ms
}

int MapHandler::matchKF2KFPoints(KeyFrame *prev_kf, KeyFrame *curr_kf) {

    int kf1_idx = prev_kf->kf_idx;
    int kf2_idx = curr_kf->kf_idx;

    StVO::StereoFrame* prev_frame = prev_kf->stereo_frame;
    StVO::StereoFrame* curr_frame = curr_kf->stereo_frame;

    matched_pt.clear();
    if (!SlamConfig::hasPoints() || curr_frame->stereo_pt.empty() || prev_frame->stereo_pt.empty())
        return 0;

    int matches = 0;
    vector<int> matches_12;

    // small window; if it fails, run standard matching
    if (SlamConfig::fastMatching()) {
        std::vector<point_2d> pj_points;
        pj_points.reserve(prev_frame->stereo_pt.size());

        for (PointFeature* pt : prev_frame->stereo_pt) {
            Vector2d point = cam->projection( DT.block(0,0,3,3) * pt->P + DT.col(3).head(3) );
            pj_points.push_back(std::make_pair(point(0) * curr_frame->inv_width, point(1) * curr_frame->inv_height));
        }

        //Fill in grid
        GridStructure grid(GRID_ROWS, GRID_COLS);
        for (int idx = 0; idx < curr_frame->stereo_pt.size(); ++idx) {
            PointFeature* point = curr_frame->stereo_pt[idx];
            grid.at(point->pl(0) * curr_frame->inv_width, point->pl(1) * curr_frame->inv_height).push_back(idx);
        }

        GridWindow w;
        int ws = SlamConfig::matchingF2FWs();
        w.width = std::make_pair(ws, ws);
        w.height = std::make_pair(ws, ws);

        matches = matchGrid(pj_points, prev_frame->pdesc_l, grid, curr_frame->pdesc_l, w, matches_12);
    }

    if (curr_frame->stereo_pt.size() > SlamConfig::minPointMatches() &&
            prev_frame->stereo_pt.size() > SlamConfig::minPointMatches() &&
            matches < SlamConfig::minPointMatches()) {
        matches = match(prev_frame->pdesc_l, curr_frame->pdesc_l, SlamConfig::minRatio12P(), matches_12);
    }

    for (int i1 = 0; i1 < matches_12.size(); ++i1) {
        const int i2 = matches_12[i1];
        if (i2 < 0) continue;

        // this shouldn't happen
        if (prev_frame->stereo_pt[i1] == nullptr)
            throw runtime_error("[MapHandler] NULL stereo point (prev)");
        if (curr_frame->stereo_pt[i2] == nullptr)
            throw runtime_error("[MapHandler] NULL stereo point (curr)");

        // new 3D landmark
        if( prev_frame->stereo_pt[i1]->idx == -1 ) {
            // assign indices
            prev_frame->stereo_pt[i1]->idx = max_pt_idx;
            curr_frame->stereo_pt[i2]->idx = max_pt_idx;
            // create new 3D landmark with the observation from previous KF
            Matrix4d Tfw = ( prev_kf->T_kf_w );
            Vector3d P3d = Tfw.block(0,0,3,3) * prev_frame->stereo_pt[i1]->P + Tfw.col(3).head(3);
            //Vector3d dir = kf0->stereo_frame->stereo_pt[lr_qdx]->P / kf0->stereo_frame->stereo_pt[lr_qdx]->P.norm();
            Vector3d dir = P3d.normalized();
            MapPoint* map_point = new MapPoint(max_pt_idx,
                                               P3d,
                                               prev_frame->pdesc_l.row(i1),
                                               kf1_idx,
                                               prev_frame->stereo_pt[i1]->pl,
                                               dir);
            // add new 3D landmark to kf_idx where it was first observed
            map_points_kf_idx.at( kf1_idx ).push_back( max_pt_idx );
            // add observation of the 3D landmark from current KF
            P3d = curr_kf->T_kf_w.block(0,0,3,3) * curr_frame->stereo_pt[i2]->P + curr_kf->T_kf_w.col(3).head(3);
            //dir = kf1->stereo_frame->stereo_pt[lr_tdx]->P / kf1->stereo_frame->stereo_pt[lr_tdx]->P.norm();
            dir = P3d / P3d.norm();
            map_point->addMapPointObservation(curr_frame->pdesc_l.row(i2),
                                              kf2_idx,
                                              curr_frame->stereo_pt[i2]->pl,
                                              dir);
            // add 3D landmark to map
            map_points.push_back(map_point);
            max_pt_idx++;
            // update full graph (new feature)
            full_graph[kf2_idx][kf1_idx]++;
            full_graph[kf1_idx][kf2_idx]++;

            // if has refine pose:
            if (SlamConfig::hasRefinement()) {
                PointFeature* pt = prev_frame->stereo_pt[i1];
                pt->pl_obs = curr_frame->stereo_pt[i2]->pl;
                pt->inlier = true;
                matched_pt.push_back(pt);
            }
        } else { // 3D landmark exists: copy idx && add observation to map landmark
            // copy idx
            int lm_idx = prev_frame->stereo_pt[i1]->idx;
            if (map_points[lm_idx] != nullptr) {
                curr_frame->stereo_pt[i2]->idx = lm_idx;
                // add observation of the 3D landmark from current KF
                Vector3d p3d = curr_kf->T_kf_w.block(0,0,3,3) * curr_frame->stereo_pt[i2]->P + curr_kf->T_kf_w.col(3).head(3);
                //Vector3d dir = kf1->stereo_frame->stereo_pt[lr_tdx]->P / kf1->stereo_frame->stereo_pt[lr_tdx]->P.norm();
                Vector3d dir = p3d.normalized();
                map_points[lm_idx]->addMapPointObservation(curr_frame->pdesc_l.row(i2),
                                                           kf2_idx,
                                                           curr_frame->stereo_pt[i2]->pl,
                                                           dir);
                // update full graph (previously observed feature)
                for (int obs : map_points[lm_idx]->kf_obs_list) {
                    if (obs != kf2_idx) {
                        full_graph[kf2_idx][obs]++;
                        full_graph[obs][kf2_idx]++;
                    }
                }

                // if has refine pose:
                if (SlamConfig::hasRefinement()) {
                    PointFeature* pt = prev_frame->stereo_pt[i1];
                    pt->pl_obs = curr_frame->stereo_pt[i2]->pl;
                    pt->inlier = true;
                    matched_pt.push_back(pt);
                }
            }
        }
    }

    return matches;
}

int MapHandler::matchKF2KFLines(KeyFrame *prev_kf, KeyFrame *curr_kf) {

    matched_ls.clear();
    if (!SlamConfig::hasLines() || curr_kf->stereo_frame->stereo_ls.empty() || prev_kf->stereo_frame->stereo_ls.empty())
        return 0;

    int kf1_idx = prev_kf->kf_idx;
    int kf2_idx = curr_kf->kf_idx;

    StVO::StereoFrame* prev_frame = prev_kf->stereo_frame;
    StVO::StereoFrame* curr_frame = curr_kf->stereo_frame;

    int matches = 0;
    vector<int> matches_12;

    // small window; if it fails, run standard matching
    if (SlamConfig::fastMatching()) {
        std::vector<line_2d> pj_lines;
        pj_lines.reserve(prev_frame->stereo_ls.size());

        for (LineFeature* ls : prev_frame->stereo_ls) {
            Vector3d sP_ = DT.block(0,0,3,3) * ls->sP + DT.col(3).head(3);
            Vector2d spl_proj = cam->projection( sP_ );

            Vector3d eP_ = DT.block(0,0,3,3) * ls->eP + DT.col(3).head(3);
            Vector2d epl_proj = cam->projection( eP_ );

            pj_lines.push_back(std::make_pair(std::make_pair(spl_proj(0), spl_proj(1)),
                                              std::make_pair(epl_proj(0), epl_proj(1))));
        }

        //Fill in grid
        list<point_2d> line_coords;
        GridStructure grid(GRID_ROWS, GRID_COLS);
        std::vector<std::pair<double, double>> directions(curr_frame->stereo_ls.size());
        for (int idx = 0; idx < curr_frame->stereo_ls.size(); ++idx) {
            LineFeature* line = curr_frame->stereo_ls[idx];

            std::pair<double, double> &v = directions[idx];
            v = std::make_pair((line->epl(0) - line->spl(0)) * curr_frame->inv_width, (line->epl(1) - line->spl(1)) * curr_frame->inv_height);
            normalize(v);

            getLineCoords(line->spl(0) * curr_frame->inv_width, line->spl(1) * curr_frame->inv_height,
                          line->epl(0) * curr_frame->inv_width, line->epl(1) * curr_frame->inv_height, line_coords);
            for (const point_2d &p : line_coords)
                grid.at(p.first, p.second).push_back(idx);
        }

        GridWindow w;
        int ws = SlamConfig::matchingF2FWs();
        w.width = std::make_pair(ws, ws);
        w.height = std::make_pair(ws, ws);

        matches = matchGrid(pj_lines, prev_frame->ldesc_l, grid, curr_frame->ldesc_l, directions, w, matches_12);
    }

    if (curr_frame->stereo_ls.size() > SlamConfig::minLineMatches() &&
            prev_frame->stereo_ls.size() > SlamConfig::minLineMatches() &&
            matches < SlamConfig::minLineMatches()) {
        matches = match(prev_frame->ldesc_l, curr_frame->ldesc_l, SlamConfig::minRatio12L(), matches_12);
//        if (matches < SlamConfig::minLineMatches()) return 0;
    }

    for (int i1 = 0; i1 < matches_12.size(); ++i1) {
        const int i2 = matches_12[i1];
        if (i2 < 0) continue;

        // this shouldn't happen
        if (prev_frame->stereo_ls[i1] == nullptr)
            throw runtime_error("[MapHandler] NULL stereo line (prev)");
        if (curr_frame->stereo_ls[i2] == nullptr)
            throw runtime_error("[MapHandler] NULL stereo line (curr)");

        // new 3D landmark
        if (prev_frame->stereo_ls[i1]->idx == -1) {
            // assign indices
            prev_frame->stereo_ls[i1]->idx = max_ls_idx;
            curr_frame->stereo_ls[i2]->idx = max_ls_idx;
            // create new 3D landmark with the observation from previous KF
            Matrix4d Tfw = ( prev_kf->T_kf_w );
            Vector3d sP3d = Tfw.block(0,0,3,3) * prev_frame->stereo_ls[i1]->sP + Tfw.col(3).head(3);
            Vector3d eP3d = Tfw.block(0,0,3,3) * prev_frame->stereo_ls[i1]->eP + Tfw.col(3).head(3);
            Vector6d L3d;
            L3d << sP3d, eP3d;
            Vector3d mP3d = 0.5*(sP3d+eP3d);
            mP3d = mP3d.normalized();
            Vector4d pts;
            pts << prev_frame->stereo_ls[i1]->spl, prev_frame->stereo_ls[i1]->epl;
            MapLine* map_line = new MapLine(max_ls_idx,
                                            L3d,
                                            prev_frame->ldesc_l.row(i1),
                                            kf1_idx,
                                            prev_frame->stereo_ls[i1]->le,
                                            mP3d,
                                            pts);
            // add new 3D landmark to kf_idx where it was first observed
            map_lines_kf_idx.at( kf1_idx ).push_back( max_ls_idx );
            // add observation of the 3D landmark from current KF
            mP3d = 0.5*( curr_frame->stereo_ls[i2]->sP + curr_frame->stereo_ls[i2]->eP );
            mP3d = curr_kf->T_kf_w.block(0,0,3,3) * mP3d + curr_kf->T_kf_w.col(3).head(3);
            mP3d = mP3d.normalized();
            pts << curr_frame->stereo_ls[i2]->spl, curr_frame->stereo_ls[i2]->epl;
            map_line->addMapLineObservation(curr_frame->ldesc_l.row(i2),
                                            kf2_idx,
                                            curr_frame->stereo_ls[i2]->le,
                                            mP3d,
                                            pts);
            // add 3D landmark to map
            map_lines.push_back(map_line);
            max_ls_idx++;
            // update full graph (new feature)
            full_graph[kf2_idx][kf1_idx]++;
            full_graph[kf1_idx][kf2_idx]++;

            // if has refine pose:
            if (SlamConfig::hasRefinement()) {
                LineFeature* ls = prev_frame->stereo_ls[i1];
                ls->sdisp_obs   = curr_frame->stereo_ls[i2]->sdisp;
                ls->edisp_obs   = curr_frame->stereo_ls[i2]->edisp;
                ls->spl_obs     = curr_frame->stereo_ls[i2]->spl;
                ls->epl_obs     = curr_frame->stereo_ls[i2]->epl;
                ls->le_obs      = curr_frame->stereo_ls[i2]->le;
                ls->inlier      = true;
                matched_ls.push_back( ls );
            }
        } else { // 3D landmark exists: copy idx && add observation to map landmark
            // copy idx
            int lm_idx = prev_frame->stereo_ls[i1]->idx;
            if (map_lines[lm_idx] != nullptr) {
                curr_frame->stereo_ls[i2]->idx = lm_idx;
                // add observation of the 3D landmark from current KF
                Vector3d mP3d = 0.5*(curr_frame->stereo_ls[i2]->sP+curr_frame->stereo_ls[i2]->eP);
                mP3d = curr_kf->T_kf_w.block(0,0,3,3) * mP3d + curr_kf->T_kf_w.col(3).head(3);
                mP3d = mP3d.normalized();
                Vector4d pts;
                pts << curr_frame->stereo_ls[i2]->spl, curr_frame->stereo_ls[i2]->epl;
                map_lines[lm_idx]->addMapLineObservation(curr_frame->ldesc_l.row(i2),
                                                         kf2_idx,
                                                         curr_frame->stereo_ls[i2]->le,
                                                         mP3d,
                                                         pts);
                // update full graph (previously observed feature)
                for (int obs : map_lines[lm_idx]->kf_obs_list) {
                    if (obs != kf2_idx) {
                        full_graph[kf2_idx][obs]++;
                        full_graph[obs][kf2_idx]++;
                    }
                }

                // if has refine pose:
                if (SlamConfig::hasRefinement()) {
                    LineFeature* ls = prev_frame->stereo_ls[i1];
                    ls->sdisp_obs   = curr_frame->stereo_ls[i2]->sdisp;
                    ls->edisp_obs   = curr_frame->stereo_ls[i2]->edisp;
                    ls->spl_obs     = curr_frame->stereo_ls[i2]->spl;
                    ls->epl_obs     = curr_frame->stereo_ls[i2]->epl;
                    ls->le_obs      = curr_frame->stereo_ls[i2]->le;
                    ls->inlier      = true;
                    matched_ls.push_back( ls );
                }
            }
        }
    }

    return matches;
}

int MapHandler::matchMap2KFPoints() {

    int kf2_idx = curr_kf->kf_idx;
    StVO::StereoFrame* curr_frame = curr_kf->stereo_frame;

    // select local map
    vector<MapPoint*> map_local_points;
    std::vector<point_2d> pj_points;
    Mat map_lpt_desc;

    if (!SlamConfig::hasPoints() || curr_frame->stereo_pt.empty())
        return 0;

    for (MapPoint* pt : map_points) {
        // if it is local and not found in the current KF
        if (pt != nullptr && pt->local && pt->kf_obs_list.back() != kf2_idx) {
            // if the LM is projected inside the current image
            Vector3d Pf = Twf.block(0,0,3,3) * pt->point3D + Twf.col(3).head(3);
            Vector2d pf = cam->projection(Pf);
            if (pf(0) > 0 && pf(0) < cam->getWidth() && pf(1) > 0 && pf(1) < cam->getHeight() && Pf(2) > 0.0) {
                // add the point and its representative descriptor
                map_local_points.push_back(pt);
                pj_points.push_back(std::make_pair(pf(0) * curr_frame->inv_width, pf(1) * curr_frame->inv_height));
                map_lpt_desc.push_back(pt->med_desc.row(0));
            }
        }
    }

    // select unmatched points
    vector<PointFeature*> unmatched_points;
    Mat unmatched_pt_desc;
    for( int idx = 0; idx < curr_kf->stereo_frame->stereo_pt.size(); ++idx) {
        PointFeature* pt = curr_kf->stereo_frame->stereo_pt[idx];
        if (pt != nullptr && pt->idx == -1) {
            unmatched_points.push_back(pt);
            unmatched_pt_desc.push_back(curr_kf->stereo_frame->pdesc_l.row(idx));
        }
    }

    if (map_local_points.empty() || unmatched_points.empty())
        return 0;

    int matches = 0;
    vector<int> matches_12;

    // track points from local map (small window; if it fails, run standard matching)
    if (SlamConfig::fastMatching()) {
        //Fill in grid
        GridStructure grid(GRID_ROWS, GRID_COLS);
        for (int idx = 0; idx < unmatched_points.size(); ++idx) {
            PointFeature* point = unmatched_points[idx];
            grid.at(point->pl(0) * curr_frame->inv_width, point->pl(1) * curr_frame->inv_height).push_back(idx);
        }

        GridWindow w;
        int ws = SlamConfig::matchingF2FWs();
        w.width = std::make_pair(ws, ws);
        w.height = std::make_pair(ws, ws);

        matches = matchGrid(pj_points, map_lpt_desc, grid, unmatched_pt_desc, w, matches_12);
    }

    if (pj_points.size() > SlamConfig::minPointMatches() &&
            map_local_points.size() > SlamConfig::minPointMatches() &&
            matches < SlamConfig::minPointMatches()) {
        matches = match(map_lpt_desc, unmatched_pt_desc, SlamConfig::minRatio12P(), matches_12);
//        if (matches < SlamConfig::minPointMatches()) return 0;
    }

    for (int i1 = 0; i1 < matches_12.size(); ++i1) {
        const int i2 = matches_12[i1];
        if (i2 < 0) continue;

        Vector3d Pf_map = Twf.block(0,0,3,3) * map_local_points[i1]->point3D + Twf.col(3).head(3);
        Vector3d Pf_kf  = unmatched_points[i2]->P;
        // check that the viewing direction condition is also satisfied
        Vector3d dir_kf  = Pf_kf.normalized();
        // check that the epipolar constraint is satisfied
        Vector2d pf_map = cam->projection( Pf_map );
        Vector2d pf_kf  = unmatched_points[i2]->pl;
        double error_epip = ( pf_map - pf_kf ).norm();
        if (error_epip < SlamConfig::maxKFEpipP()) {
            // copy idx
            int lm_idx = map_local_points[i1]->idx;
            unmatched_points[i2]->idx = lm_idx;
            // add observation of the 3D LM from current KF
            dir_kf = curr_kf->T_kf_w.block(0,0,3,3) * dir_kf + curr_kf->T_kf_w.col(3).head(3);
            map_points[lm_idx]->addMapPointObservation( unmatched_pt_desc.row(i2), kf2_idx, unmatched_points[i2]->pl, dir_kf );
            // update full graph (previously observed feature)
            for (int obs : map_points[lm_idx]->kf_obs_list) {
                if (obs != kf2_idx) {
                    full_graph[kf2_idx][obs]++;
                    full_graph[obs][kf2_idx]++;
                }
            }
        }
        else --matches;
    }

    return matches;
}

int MapHandler::matchMap2KFLines() {

    int kf2_idx = curr_kf->kf_idx;
    StVO::StereoFrame* curr_frame = curr_kf->stereo_frame;

    // select local map
    vector<MapLine*> map_local_lines;
    std::vector<line_2d> pj_lines;
    Mat map_lls_desc;

    if (!SlamConfig::hasLines() || curr_frame->stereo_ls.empty())
        return 0;

    for (MapLine* ls : map_lines) {
        if (ls != nullptr && ls->local && ls->kf_obs_list.back() != kf2_idx) {
            // if the LM is projected inside the current image
            Vector3d sPf = Twf.block(0,0,3,3) * ls->line3D.head(3) + Twf.col(3).head(3);
            Vector2d spf = cam->projection( sPf );
            Vector3d ePf = Twf.block(0,0,3,3) * ls->line3D.tail(3) + Twf.col(3).head(3);
            Vector2d epf = cam->projection( ePf );
            if (spf(0) > 0 && spf(0) < cam->getWidth() && spf(1) > 0 && spf(1) < cam->getHeight() && sPf(2) > 0.0 &&
                    epf(0) > 0 && epf(0) < cam->getWidth() && epf(1) > 0 && epf(1) < cam->getHeight() && ePf(2) > 0.0) {
                // add the point and its representative descriptor
                map_local_lines.push_back( ls );
                pj_lines.push_back(std::make_pair(std::make_pair(spf(0) * curr_frame->inv_width, spf(1) * curr_frame->inv_height),
                                                  std::make_pair(epf(0) * curr_frame->inv_width, epf(1) * curr_frame->inv_height)));
                map_lls_desc.push_back( ls->med_desc.row(0) );
            }
        }
    }

    // select unmatched line segments
    vector<LineFeature*> unmatched_lines;
    Mat unmatched_ls_desc;
    for( int idx = 0; idx < curr_frame->stereo_ls.size(); ++idx) {
        LineFeature *ls = curr_frame->stereo_ls[idx];
        if (ls != nullptr && ls->idx == -1) {
            unmatched_lines.push_back( ls );
            unmatched_ls_desc.push_back( curr_frame->ldesc_l.row(idx) );
        }
    }

    if (map_local_lines.empty() || unmatched_lines.empty())
        return 0;

    int matches = 0;
    vector<int> matches_12;

    // track lines from local map (small window; if it fails, run standard matching)
    if (SlamConfig::fastMatching()) {
        //Fill in grid
        list<point_2d> line_coords;
        GridStructure grid(GRID_ROWS, GRID_COLS);
        std::vector<std::pair<double, double>> directions(unmatched_lines.size());
        for (int idx = 0; idx < unmatched_lines.size(); ++idx) {
            LineFeature* line = unmatched_lines[idx];

            std::pair<double, double> &v = directions[idx];
            v = std::make_pair((line->epl(0) - line->spl(0)) * curr_frame->inv_width, (line->epl(1) - line->spl(1)) * curr_frame->inv_height);
            normalize(v);

            getLineCoords(line->spl(0) * curr_frame->inv_width, line->spl(1) * curr_frame->inv_height,
                          line->epl(0) * curr_frame->inv_width, line->epl(1) * curr_frame->inv_height, line_coords);
            for (const point_2d &p : line_coords)
                grid.at(p.first, p.second).push_back(idx);
        }

        GridWindow w;
        int ws = SlamConfig::matchingF2FWs();
        w.width = std::make_pair(ws, ws);
        w.height = std::make_pair(ws, ws);

        matches = matchGrid(pj_lines, map_lls_desc, grid, unmatched_ls_desc, directions, w, matches_12);
    }

    if (pj_lines.size() > SlamConfig::minLineMatches() &&
            map_local_lines.size() > SlamConfig::minLineMatches() &&
            matches < SlamConfig::minLineMatches()) {
        matches = match(map_lls_desc, unmatched_ls_desc, SlamConfig::minRatio12L(), matches_12);
//        if (matches < SlamConfig::minLineMatches()) return 0;
    }

    for (int i1 = 0; i1 < matches_12.size(); ++i1) {
        const int i2 = matches_12[i1];
        if (i2 < 0) continue;

        Vector3d sP_ = Twf.block(0,0,3,3) * map_local_lines[i1]->line3D.head(3) + Twf.col(3).head(3);
        Vector2d spl_proj = cam->projection( sP_ );
        Vector3d eP_ = Twf.block(0,0,3,3) * map_local_lines[i1]->line3D.tail(3) + Twf.col(3).head(3);
        Vector2d epl_proj = cam->projection( eP_ );
        Vector3d l_obs = unmatched_lines[i2]->le;
        // check the epipolar constraint
        Vector2d err_ls;
        err_ls(0) = l_obs(0) * spl_proj(0) + l_obs(1) * spl_proj(1) + l_obs(2);
        err_ls(1) = l_obs(0) * epl_proj(0) + l_obs(1) * epl_proj(1) + l_obs(2);
        if( err_ls(0) < SlamConfig::maxKFEpipL() && err_ls(1) < SlamConfig::maxKFEpipL() ) {
            // copy idx
            int lm_idx = map_local_lines[i1]->idx;
            unmatched_lines[i2]->idx = lm_idx;
            // add observation of the 3D landmark from current KF
            Vector3d mP3d = 0.5*(unmatched_lines[i2]->sP + unmatched_lines[i2]->eP);
            mP3d = curr_kf->T_kf_w.block(0,0,3,3) * mP3d + curr_kf->T_kf_w.col(3).head(3);
            mP3d = mP3d.normalized();
            Vector4d pts;
            pts << unmatched_lines[i2]->spl, unmatched_lines[i2]->epl;
            map_lines[lm_idx]->addMapLineObservation( unmatched_ls_desc.row(i2), kf2_idx, unmatched_lines[i2]->le, mP3d, pts );
            // update full graph (previously observed feature)
            for (int obs : map_lines[lm_idx]->kf_obs_list) {
                if (obs != kf2_idx) {
                    full_graph[kf2_idx][obs]++;
                    full_graph[obs][kf2_idx]++;
                }
            }
        }
        else --matches;
    }

    return matches;
}

void MapHandler::lookForCommonMatches( KeyFrame* kf0, KeyFrame* &kf1 )
{
    // ---------------------------------------------------
    // find matches between prev_keyframe and curr_frame
    // ---------------------------------------------------
    // points f2f tracking
    int common_pt = matchKF2KFPoints(prev_kf, curr_kf);

    // line segments f2f tracking
    int common_ls = matchKF2KFLines(prev_kf, curr_kf);

    // ---------------------------------------------------
    // refine pose between kf0 and kf1
    // ---------------------------------------------------
    if (SlamConfig::hasRefinement()) {
        StVO::StereoFrameHandler* stf = new StereoFrameHandler( cam );
        stf->matched_pt = matched_pt;
        stf->matched_ls = matched_ls;

        stf->prev_frame = kf0->stereo_frame;
        stf->curr_frame = kf1->stereo_frame;

        stf->n_inliers_pt = stf->matched_pt.size();
        stf->n_inliers_ls = stf->matched_ls.size();
        stf->n_inliers    = stf->n_inliers_pt + stf->n_inliers_ls;

        stf->optimizePose();

        double inl_ratio_pt = 100.0 * stf->n_inliers_pt / matched_pt.size();
        double inl_ratio_ls = 100.0 * stf->n_inliers_ls / matched_ls.size();

        bool condition_pt = true, condition_ls = true;
        if (SlamConfig::hasPoints())
            condition_pt = (inl_ratio_pt >= SlamConfig::kfInlierRatio());
        if (SlamConfig::hasLines())
            condition_ls = (inl_ratio_ls >= SlamConfig::kfInlierRatio());
        if (!SlamConfig::hasPoints() && !SlamConfig::hasLines()) {
            condition_pt = false;
            condition_ls = false;
        }

        bool inl_ratio_condition = (condition_pt && condition_ls);

        if (stf->n_inliers > SlamConfig::minFeatures() && inl_ratio_condition) {
            Matrix4d DT_ = stf->curr_frame->DT;
            kf1->T_kf_w = expmap_se3(logmap_se3( kf0->T_kf_w * DT_ ));
        } else
            kf1->T_kf_w = expmap_se3(logmap_se3( kf0->T_kf_w * inverse_se3(DT) ));

        // update DT & Twf
        Twf = expmap_se3(logmap_se3(  inverse_se3( curr_kf->T_kf_w ) ));
        DT = expmap_se3(logmap_se3( Twf * prev_kf->T_kf_w ));

        delete stf;
    }

    // ---------------------------------------------------
    // find point matches between local map and curr_frame
    // ---------------------------------------------------
    if (SlamConfig::hasPoints())
        common_pt += matchMap2KFPoints();

    // ---------------------------------------------------
    // find line matches between local map and curr_frame
    // ---------------------------------------------------
    if (SlamConfig::hasLines())
        common_ls += matchMap2KFLines();
}

void MapHandler::expandGraphs()
{
    int g_size = full_graph.size() + 1;
    // full graph
    full_graph.resize( g_size );
    for(unsigned int i = 0; i < g_size; i++ )
        full_graph[i].resize( g_size );
    // confusion matrix
    conf_matrix.resize( g_size );
    for(unsigned int i = 0; i < g_size; i++ )
        conf_matrix[i].resize( g_size );
}

void MapHandler::formLocalMap()
{

    // for the Single Thread version

    // reset local KFs & LMs
    for( vector<KeyFrame*>::iterator kf_it = map_keyframes.begin(); kf_it != map_keyframes.end(); kf_it++ )
    {
        if( (*kf_it) != NULL )
            (*kf_it)->local = false;
    }
    for( vector<MapPoint*>::iterator pt_it = map_points.begin(); pt_it != map_points.end(); pt_it++ )
    {
        if( (*pt_it) != NULL )
            (*pt_it)->local = false;
    }
    for( vector<MapLine*>::iterator  ls_it = map_lines.begin(); ls_it != map_lines.end(); ls_it++ )
    {
        if( (*ls_it) != NULL )
            (*ls_it)->local = false;
    }

    // set first KF and their associated LMs as local
    map_keyframes.back()->local = true;
    for( vector<PointFeature*>::iterator pt_it = map_keyframes.back()->stereo_frame->stereo_pt.begin(); pt_it != map_keyframes.back()->stereo_frame->stereo_pt.end(); pt_it++ )
    {
        if( (*pt_it) != NULL )
        {
            int lm_idx = (*pt_it)->idx;
            if( lm_idx != -1 && map_points[lm_idx] != NULL )
                map_points[lm_idx]->local = true;
        }
    }
    for( vector<LineFeature*>::iterator ls_it = map_keyframes.back()->stereo_frame->stereo_ls.begin(); ls_it != map_keyframes.back()->stereo_frame->stereo_ls.end(); ls_it++ )
    {
        if( (*ls_it) != NULL )
        {
            int lm_idx = (*ls_it)->idx;
            if( lm_idx != -1 && map_lines[lm_idx] != NULL )
                map_lines[lm_idx]->local = true;
        }
    }

    // loop over covisibility graph / full graph if we want to find more points
    int g_size = full_graph.size()-1;
    for( int i = 0; i < g_size; i++ )
    {
        if( full_graph[g_size][i] >= SlamConfig::minLMCovGraph() || abs(g_size-i) <= SlamConfig::minKFLocalMap() )
        {
            map_keyframes[i]->local = true;
            // loop over the landmarks seen by KF{i}
            for( vector<PointFeature*>::iterator pt_it = map_keyframes[i]->stereo_frame->stereo_pt.begin(); pt_it != map_keyframes[i]->stereo_frame->stereo_pt.end(); pt_it++ )
            {
                int lm_idx = (*pt_it)->idx;
                if( lm_idx != -1 && map_points[lm_idx] != NULL )
                    map_points[lm_idx]->local = true;
            }
            for( vector<LineFeature*>::iterator ls_it = map_keyframes[i]->stereo_frame->stereo_ls.begin(); ls_it != map_keyframes[i]->stereo_frame->stereo_ls.end(); ls_it++ )
            {
                int lm_idx = (*ls_it)->idx;
                if( lm_idx != -1 && map_lines[lm_idx] != NULL )
                    map_lines[lm_idx]->local = true;
            }
        }
    }

}

void MapHandler::formLocalMap( KeyFrame * kf )
{

    // reset local KFs & LMs
    for( vector<KeyFrame*>::iterator kf_it = map_keyframes.begin(); kf_it != map_keyframes.end(); kf_it++ )
    {
        if( (*kf_it) != NULL )
            (*kf_it)->local = false;
    }
    for( vector<MapPoint*>::iterator pt_it = map_points.begin(); pt_it != map_points.end(); pt_it++ )
    {
        if( (*pt_it) != NULL )
            (*pt_it)->local = false;
    }
    for( vector<MapLine*>::iterator  ls_it = map_lines.begin(); ls_it != map_lines.end(); ls_it++ )
    {
        if( (*ls_it) != NULL )
            (*ls_it)->local = false;
    }

    // set first KF and their associated LMs as local
    kf->local = true;
    for( vector<PointFeature*>::iterator pt_it = kf->stereo_frame->stereo_pt.begin(); pt_it != kf->stereo_frame->stereo_pt.end(); pt_it++ )
    {
        if( (*pt_it) != NULL )
        {
            int lm_idx = (*pt_it)->idx;
            if( lm_idx != -1 && map_points[lm_idx] != NULL )
                map_points[lm_idx]->local = true;
        }
    }
    for( vector<LineFeature*>::iterator ls_it = kf->stereo_frame->stereo_ls.begin(); ls_it != kf->stereo_frame->stereo_ls.end(); ls_it++ )
    {
        if( (*ls_it) != NULL )
        {
            int lm_idx = (*ls_it)->idx;
            if( lm_idx != -1 && map_lines[lm_idx] != NULL )
                map_lines[lm_idx]->local = true;
        }
    }

    // loop over covisibility graph / full graph if we want to find more points
    int g_size = full_graph.size()-1;
    for( int i = 0; i < g_size; i++ )
    {
        if( full_graph[g_size][i] >= SlamConfig::minLMCovGraph() || abs(g_size-i) <= SlamConfig::minKFLocalMap() )
        {
            map_keyframes[i]->local = true;
            // loop over the landmarks seen by KF{i}
            for( vector<PointFeature*>::iterator pt_it = map_keyframes[i]->stereo_frame->stereo_pt.begin(); pt_it != map_keyframes[i]->stereo_frame->stereo_pt.end(); pt_it++ )
            {
                int lm_idx = (*pt_it)->idx;
                if( lm_idx != -1 && map_points[lm_idx] != NULL )
                    map_points[lm_idx]->local = true;
            }
            for( vector<LineFeature*>::iterator ls_it = map_keyframes[i]->stereo_frame->stereo_ls.begin(); ls_it != map_keyframes[i]->stereo_frame->stereo_ls.end(); ls_it++ )
            {
                int lm_idx = (*ls_it)->idx;
                if( lm_idx != -1 && map_lines[lm_idx] != NULL )
                    map_lines[lm_idx]->local = true;
            }
        }
    }

    //test
//    MapPoint* temp = map_points.back();
//    vector<int> obs = temp->kf_obs_list;
//    for(int i=0;i<obs.size();i++){
//        cout<<"obs["<<"]: "<<obs[i]<<"  map_keyframe idx: "<<map_keyframes[obs[i]]->kf_idx<<endl;
//    }

}

// -----------------------------------------------------------------------------------------------------------------------------
// Parallelization functions
// -----------------------------------------------------------------------------------------------------------------------------

void MapHandler::addKeyFrame_multiThread(KeyFrame *curr_kf , KeyFrame *prev_kf) {

    if (!threads_started) return;

    {
        std::lock_guard<std::mutex> lk(kf_queue_mutex);
        kf_queue.push_back(std::make_pair(curr_kf,prev_kf));
    }
    new_kf.notify_one();
}

void MapHandler::handlerThread() {

    if (!threads_started) return;

    while (true) {

        std::unique_lock<std::mutex> lk(kf_queue_mutex);
        if (kf_queue.empty())
            new_kf.wait(lk, [this]{return !kf_queue.empty();});

        curr_kf_mt = kf_queue.front().first;
        prev_kf_mt = kf_queue.front().second;
        //addKeyframeToSW(curr_kf_mt);
        kf_queue.pop_front();

        lk.unlock();

        // notify threads
        {
            std::lock_guard<std::mutex> lk(lba_mutex);
            lba_thread_status = LBA_ACTIVE;
        }
        lba_start.notify_one();

        //if closing loopclosing
        {
            std::lock_guard<std::mutex> lk(lc_mutex);
            lc_thread_status = LC_ACTIVE;
        }
        lc_start.notify_one();
        //end closing loopclosing

        if( curr_kf_mt == nullptr || prev_kf_mt == nullptr ) return;

        // join localMapping and loopClosure threads
        std::unique_lock<std::mutex> lba_lk(lba_mutex);
        if( lba_thread_status != LBA_IDLE )
            lba_join.wait(lba_lk, [this]{return (lba_thread_status == LBA_IDLE);});

        //if close loopclosing
        std::unique_lock<std::mutex> lc_lk(lc_mutex);
        if ( lc_thread_status != LC_IDLE )
            lc_join.wait(lc_lk, [this]{return (lc_thread_status == LC_IDLE);});

        // loop closure
        if( lc_state == LC_READY )
        {
            loopClosureOptimizationCovGraphG2O();
            lc_state = LC_IDLE;
        }

    }


}

void MapHandler::startThreads() {

    if (threads_started) return;
    threads_started = true;

    std::thread handler(&MapHandler::handlerThread, this);
    handler.detach();

    {
        std::lock_guard<std::mutex> lk(lba_mutex);
        lba_thread_status = LBA_IDLE;
    }
    std::thread localMapping(&MapHandler::localMappingThread, this);
    localMapping.detach();

    //if closing loopclosing
    {
        std::lock_guard<std::mutex> lk(lc_mutex);
        lc_thread_status = LC_IDLE;
    }
    std::thread loopClosure(&MapHandler::loopClosureThread, this);
    loopClosure.detach();
}

void MapHandler::killThreads() {

    if (!threads_started) return;
    threads_started = false;

    {
        std::lock_guard<std::mutex> lk(kf_queue_mutex);
        kf_queue.push_back(std::make_pair(nullptr,nullptr));
    }
    new_kf.notify_one();

    print_msg("[Waiting for threads to finish...");

    std::unique_lock<std::mutex> lba_lk(lba_mutex);
    if (lba_thread_status != LBA_TERMINATED)
        lba_join.wait(lba_lk, [this]{return (lba_thread_status == LBA_TERMINATED);});

    //if closing loopclosing
    std::unique_lock<std::mutex> lc_lk(lc_mutex);
    if (lc_thread_status != LC_TERMINATED)
        lc_join.wait(lc_lk, [this]{return (lc_thread_status == LC_TERMINATED);});
}

void MapHandler::localMappingThread() {

    if (!threads_started) return;

    std::unique_lock<std::mutex> lk(lba_mutex, std::defer_lock);
    while (true) {

        lk.lock();
        if (lba_thread_status != LBA_ACTIVE)
            lba_start.wait(lk, [this]{return (lba_thread_status == LBA_ACTIVE);});
        lk.unlock();
        addKeyframeToSW(curr_kf_mt);
        if (curr_kf_mt == nullptr || prev_kf_mt == nullptr) break;

        // reset indices
        for (PointFeature* pt : curr_kf_mt->stereo_frame->stereo_pt)
            pt->idx = -1;
        for (LineFeature* ls : curr_kf_mt->stereo_frame->stereo_ls)
            ls->idx = -1;

        // look for common matches and update the full graph
        lookForCommonMatches( prev_kf_mt, curr_kf_mt );



        // perform local bundle adjustment
        if(!vioInit){
            // form local map
            formLocalMap(curr_kf_mt);
            int lba = localBundleAdjustment();
            tryVioInit();
        }
#ifndef USE_MARG
        else{
            lk_sw.lock();
            std::deque<KeyFrame*> vpKFs = slideWindow;
            lk_sw.unlock();

            //debug
//            Matrix4d T = cam->getTbs() * prev_kf_mt->T_kf_w.inverse() * curr_kf_mt->T_kf_w * cam->getTbs().inverse();
//            Quaterniond qdt(T.block<3,3>(0,0));
//            Quaterniond qdp(curr_kf_mt->GetIMUPreInt().getDeltaR());
//            std::cout<<"DT: "<<qdt.coeffs().transpose()<<std::endl;
//            std::cout<<"IMU preintegrator: "<<qdp.coeffs().transpose()<<std::endl;

           // curr_kf_mt->UpdateNavState(prev_kf_mt->GetNavState(), curr_kf_mt->GetIMUPreInt(), gw);

            //debug
//            Matrix4d T = curr_kf_mt->T_kf_w * cam->getTbs().inverse();
//            Quaterniond qdt(T.block<3,3>(0,0));
//            Vector3d tdt(T.block<3,1>(0,3));
//            Quaterniond qimu(curr_kf_mt->GetNavState().Get_RotMatrix());
//            Vector3d timu(curr_kf_mt->GetNavState().Get_P());
//            std::cout<<"Current Rot from VO:  "<<qdt.coeffs().transpose()<<"   Trans: "<<tdt.transpose()<<std::endl;
//            std::cout<<"Current Rot from IMU: "<<qimu.coeffs().transpose()<<"   Trans: "<<timu.transpose()<<std::endl;
            mbaAbort = false;
            formLocalMap(vpKFs);
            int lba = localBundleAdjustmentWithImu(vpKFs,&mbaAbort);
        }
#else
        else{
            lk_sw.lock();
            std::deque<KeyFrame*> vpKFs = slideWindow;
            lk_sw.unlock();

            mbaAbort = false;
            int lba = localBundleAdjustmentWithImuAndMarg(vpKFs,&mbaAbort);
        }

#endif
        // recent map LMs culling (implement filters for line segments, which seems to be unaccurate)
        removeBadMapLandmarks();

        while(slideWindow.size()>=max_keyframe_in_slideWindow){
            deleteKeyframeInSW();
        }

        lk.lock();
        lba_thread_status = LBA_IDLE;
        lk.unlock();

        lba_join.notify_one();
    }

    lk.lock();
    lba_thread_status = LBA_TERMINATED;
    lk.unlock();

    lba_join.notify_one();

    print_msg("[localMappingThread] terminated.");
}

void MapHandler::loopClosureThread() {

    if (!threads_started) return;

    std::unique_lock<std::mutex> lk(lc_mutex, std::defer_lock);
    while (true) {

        lk.lock();
        if (lc_thread_status != LC_ACTIVE)
            lc_start.wait(lk, [this]{return (lc_thread_status == LC_ACTIVE);});
        lk.unlock();

        if ( curr_kf_mt == nullptr && prev_kf_mt == nullptr ) break; // stop loop closure thread

        // insert BOW vector
        if( SlamConfig::hasPoints() && SlamConfig::hasLines() )
            insertKFBowVectorPL(curr_kf_mt);
        else if( SlamConfig::hasPoints() && !SlamConfig::hasLines() )
            insertKFBowVectorP(curr_kf_mt);
        else if( !SlamConfig::hasPoints() && SlamConfig::hasLines() )
            insertKFBowVectorL(curr_kf_mt);

        // look for loop closure candidates
        int lc_kf_idx = -1;
        lookForLoopCandidates(curr_kf_mt->kf_idx, lc_kf_idx);
        //bool inl_ratio_condition = false;
        if( lc_kf_idx >= 0 )
        {

            vector<Vector4i> lc_pt_idx, lc_ls_idx;
            vector<PointFeature*> lc_points;
            vector<LineFeature*>  lc_lines;
            Vector6d pose_inc;

            bool isLC = isLoopClosure( map_keyframes[lc_kf_idx], curr_kf_mt, pose_inc, lc_pt_idx, lc_ls_idx, lc_points, lc_lines );

            // if it is loop closure, add information and update status
            if( isLC )
            {
                lc_pt_idxs.push_back( lc_pt_idx );
                lc_ls_idxs.push_back( lc_ls_idx );
                lc_poses.push_back( pose_inc );
                lc_pose_list.push_back( pose_inc );
                Vector3i lc_idx;
                lc_idx(0) = map_keyframes[lc_kf_idx]->kf_idx;
                lc_idx(1) = curr_kf_mt->kf_idx;
                lc_idx(2) = 1;
                lc_idxs.push_back( lc_idx );
                lc_idx_list.push_back(lc_idx);
                if( lc_state == LC_IDLE )
                    lc_state = LC_ACTIVE;
            }
            else
            {
                if( lc_state == LC_ACTIVE )
                    lc_state = LC_READY;
            }

            for (PointFeature* pt : lc_points)
                delete pt;
            for (LineFeature* ls : lc_lines)
                delete ls;
        }
        else
        {
            if( lc_state == LC_ACTIVE )
                lc_state = LC_READY;
        }

        lk.lock();
        lc_thread_status = LC_IDLE;
        lk.unlock();

        lc_join.notify_one();

    }

    lk.lock();
    lc_thread_status = LC_TERMINATED;
    lk.unlock();

    lc_join.notify_one();

    print_msg("[LoopClosureThread] terminated.");
}

// -----------------------------------------------------------------------------------------------------------------------------
// Local Bundle Adjustment functions
// -----------------------------------------------------------------------------------------------------------------------------

int MapHandler::localBundleAdjustment()
{

    vector<double> X_aux;

    // create list of local keyframes
    vector<int> kf_list;
    for( vector<KeyFrame*>::iterator kf_it = map_keyframes.begin(); kf_it != map_keyframes.end(); kf_it++)
    {
        if( (*kf_it)!= NULL )
        {
            if( (*kf_it)->local && (*kf_it)->kf_idx != 0 )
            {
                Vector6d pose_aux = (*kf_it)->x_kf_w;
                for(int i = 0; i < 6; i++)
                    X_aux.push_back( pose_aux(i) );
                kf_list.push_back( (*kf_it)->kf_idx );
            }
        }
    }

    // create list of local point landmarks
    vector<Vector6i> pt_obs_list;
    vector<int> pt_list;
    int lm_local_idx = 0;
    for( vector<MapPoint*>::iterator pt_it = map_points.begin(); pt_it != map_points.end(); pt_it++)
    {
        if( (*pt_it)!= NULL )
        {
            if( (*pt_it)->local )
            {
                Vector3d point_aux = (*pt_it)->point3D;
                for(int i = 0; i < 3; i++)
                    X_aux.push_back( point_aux(i) );
                // gather all observations
                for( int i = 0; i < (*pt_it)->obs_list.size(); i++)
                {
                    Vector6i obs_aux;
                    obs_aux(0) = (*pt_it)->idx; // LM idx
                    obs_aux(1) = lm_local_idx;  // LM local idx
                    obs_aux(2) = i;             // LM obs idx
                    int kf_obs_list_ = (*pt_it)->kf_obs_list[i];
                    obs_aux(3) = kf_obs_list_;  // KF idx
                    obs_aux(4) = -1;            // KF local idx (-1 if not local)
                    obs_aux(5) = 1;             // 1 if the observation is an inlier
                    for( int j = 0; j < kf_list.size(); j++ )
                    {
                        if( kf_list[j] == kf_obs_list_ )
                        {
                            obs_aux(4) = j;
                            break;
                        }
                    }
                    pt_obs_list.push_back( obs_aux );
                }
                lm_local_idx++;
                // pt_list
                pt_list.push_back( (*pt_it)->idx );
            }
        }
    }

    // create list of local line segment landmarks
    vector<Vector6i> ls_obs_list;
    vector<int> ls_list;
    lm_local_idx = 0;
    for( vector<MapLine*>::iterator ls_it = map_lines.begin(); ls_it != map_lines.end(); ls_it++)
    {
        if( (*ls_it)!= NULL )
        {
            if( (*ls_it)->local )
            {
                Vector6d line_aux = (*ls_it)->line3D;
                for(int i = 0; i < 6; i++)
                    X_aux.push_back( line_aux(i) );
                // gather all observations
                for( int i = 0; i < (*ls_it)->obs_list.size(); i++)
                {
                    Vector6i obs_aux;
                    obs_aux(0) = (*ls_it)->idx; // LM idx
                    obs_aux(1) = lm_local_idx;  // LM local idx
                    obs_aux(2) = i;             // LM obs idx
                    int kf_obs_list_ = (*ls_it)->kf_obs_list[i];
                    obs_aux(3) = kf_obs_list_;  // KF idx
                    obs_aux(4) = -1;            // KF local idx (-1 if not local)
                    obs_aux(5) = 1;             // 1 if the observation is an inlier
                    for( int j = 0; j < kf_list.size(); j++ )
                    {
                        if( kf_list[j] == kf_obs_list_ )
                        {
                            obs_aux(4) = j;
                            break;
                        }
                    }
                    ls_obs_list.push_back( obs_aux );
                }
                lm_local_idx++;
                // ls_list
                ls_list.push_back( (*ls_it)->idx );
            }
        }
    }

    // Levenberg-Marquardt optimization of the local map
    if( pt_obs_list.size() + ls_obs_list.size() != 0 )
        //return levMarquardtOptimizationPoseOnlyLBA(X_aux,kf_list,pt_list,ls_list,pt_obs_list,ls_obs_list);
        return levMarquardtOptimizationLBA(X_aux,kf_list,pt_list,ls_list,pt_obs_list,ls_obs_list);
    else
        return -1;

}

int MapHandler::levMarquardtOptimizationLBA( vector<double> X_aux, vector<int> kf_list, vector<int> pt_list, vector<int> ls_list, vector<Vector6i> pt_obs_list, vector<Vector6i> ls_obs_list  )
{

    // create Levenberg-Marquardt variables
    int    Nkf = kf_list.size();
    int      N = X_aux.size();

    VectorXd X = VectorXd::Zero(N), DX = VectorXd::Zero(N);
    VectorXd g = VectorXd::Zero(N);
    MatrixXd H = MatrixXd::Zero(N,N);

    for(int i = 0; i < N; i++)
        X(i) = X_aux[i];
    SparseMatrix<double> H_(N,N);

    // create Levenberg-Marquardt parameters
    double err = 0.0, err_prev = 999999999.9;
    double lambda = SlamConfig::lambdaLbaLM(), lambda_k = SlamConfig::lambdaLbaK();
    int    max_iters = SlamConfig::maxItersLba();

    // estimate H and g to precalculate lambda
    //---------------------------------------------------------------------------------------------
    // point observations
    int Npt = 0, Npt_obs = 0;
    if( pt_obs_list.size() != 0 )
        Npt = pt_obs_list.back()(1)+1;
    for( vector<Vector6i>::iterator pt_it = pt_obs_list.begin(); pt_it != pt_obs_list.end(); pt_it++ )
    {
        int lm_idx_map = (*pt_it)(0);
        int lm_idx_loc = (*pt_it)(1);
        int lm_idx_obs = (*pt_it)(2);
        int kf_idx_map = (*pt_it)(3);
        int kf_idx_loc = (*pt_it)(4);
        if( map_points[lm_idx_map] != NULL && map_keyframes[kf_idx_map] != NULL)
        {
            // grab 3D LM (Xwj)
            Vector3d Xwj   = map_points[lm_idx_map]->point3D;
            // grab 6DoF KF (Tiw)
            Matrix4d Tiw   = map_keyframes[kf_idx_map]->T_kf_w;
            // projection error
            Tiw = inverse_se3( Tiw );
            Vector3d Xwi   = Tiw.block(0,0,3,3) * Xwj + Tiw.block(0,3,3,1);
            Vector2d p_prj = cam->projection( Xwi );
            Vector2d p_obs = map_points[lm_idx_map]->obs_list[lm_idx_obs];
            Vector2d p_err    = p_obs - p_prj;
            double p_err_norm = p_err.norm();
            // estimate useful variables
            double gx   = Xwi(0);
            double gy   = Xwi(1);
            double gz   = Xwi(2);
            double gz2  = gz*gz;
            gz2         = 1.0 / std::max(SlamConfig::homogTh(),gz2);
            double fx   = cam->getFx();
            double fy   = cam->getFy();
            double dx   = p_err(0);
            double dy   = p_err(1);
            double fxdx = fx*dx;
            double fydy = fy*dy;
            // estimate Jacobian wrt KF pose
            Vector6d Jij_Tiw = Vector6d::Zero();
            Jij_Tiw << + gz2 * fxdx * gz,
                       + gz2 * fydy * gz,
                       - gz2 * ( fxdx*gx + fydy*gy ),
                       - gz2 * ( fxdx*gx*gy + fydy*gy*gy + fydy*gz*gz ),
                       + gz2 * ( fxdx*gx*gx + fxdx*gz*gz + fydy*gx*gy ),
                       + gz2 * ( fydy*gx*gz - fxdx*gy*gz );
            Jij_Tiw = Jij_Tiw / std::max(SlamConfig::homogTh(),p_err_norm);
            // estimate Jacobian wrt LM
            Vector3d Jij_Xwj = Vector3d::Zero();
            Jij_Xwj << + gz2 * fxdx * gz,
                       + gz2 * fydy * gz,
                       - gz2 * ( fxdx*gx + fydy*gy );
            Jij_Xwj = Jij_Xwj.transpose() * Tiw.block(0,0,3,3) / std::max(SlamConfig::homogTh(),p_err_norm);
            // if employing robust cost function            
            double w = 1.0;
            w = robustWeightCauchy(p_err_norm) ;

            // update hessian, gradient, and error
            MatrixXd Haux  = MatrixXd::Zero(3,6);
            int idx = 6 * kf_idx_loc;
            int jdx = 6*Nkf + 3*lm_idx_loc;
            if( kf_idx_loc == -1 )
            {
                g.block(jdx,0,3,1) += Jij_Xwj * p_err_norm * w ;
                err += p_err_norm * p_err_norm * w ;
                H.block(jdx,jdx,3,3) += Jij_Xwj * Jij_Xwj.transpose() * w ;
            }
            else
            {
                g.block(idx,0,6,1) += Jij_Tiw * p_err_norm * w;
                g.block(jdx,0,3,1) += Jij_Xwj * p_err_norm * w;
                err += p_err_norm * p_err_norm * w;
                Haux = Jij_Xwj * Jij_Tiw.transpose() * w ;
                H.block(idx,idx,6,6) += Jij_Tiw * Jij_Tiw.transpose() * w;
                H.block(jdx,idx,3,6) += Haux;
                H.block(idx,jdx,6,3) += Haux.transpose();
                H.block(jdx,jdx,3,3) += Jij_Xwj * Jij_Xwj.transpose() * w;
            }
        }
    }
    // line segment observations
    int Nls = 0, Nls_obs = 0;
    if( ls_obs_list.size() != 0 )
        Nls = ls_obs_list.back()(1)+1;
    for( vector<Vector6i>::iterator ls_it = ls_obs_list.begin(); ls_it != ls_obs_list.end(); ls_it++ )
    {
        int lm_idx_map = (*ls_it)(0);
        int lm_idx_loc = (*ls_it)(1);
        int lm_idx_obs = (*ls_it)(2);
        int kf_idx_map = (*ls_it)(3);
        int kf_idx_loc = (*ls_it)(4);
        if( map_lines[lm_idx_map] != NULL && map_keyframes[kf_idx_map] != NULL)
        {
            // grab 3D LM (Pwj and Qwj)
            Vector3d Pwj   = map_lines[lm_idx_map]->line3D.head(3);
            Vector3d Qwj   = map_lines[lm_idx_map]->line3D.tail(3);
            // grab 6DoF KF (Tiw)
            Matrix4d Tiw   = map_keyframes[kf_idx_map]->T_kf_w;
            // projection error
            Tiw = inverse_se3( Tiw );
            Vector3d Pwi   = Tiw.block(0,0,3,3) * Pwj + Tiw.block(0,3,3,1);
            Vector3d Qwi   = Tiw.block(0,0,3,3) * Qwj + Tiw.block(0,3,3,1);
            Vector2d p_prj = cam->projection( Pwi );
            Vector2d q_prj = cam->projection( Qwi );
            Vector3d l_obs = map_lines[lm_idx_map]->obs_list[lm_idx_obs];
            Vector2d l_err;
            l_err(0) = l_obs(0) * p_prj(0) + l_obs(1) * p_prj(1) + l_obs(2);
            l_err(1) = l_obs(0) * q_prj(0) + l_obs(1) * q_prj(1) + l_obs(2);
            double l_err_norm = l_err.norm();
            // start point
            double gx   = Pwi(0);
            double gy   = Pwi(1);
            double gz   = Pwi(2);
            double gz2  = gz*gz;
            gz2         = 1.0 / std::max(SlamConfig::homogTh(),gz2);
            double fx   = cam->getFx();
            double fy   = cam->getFy();
            double lx   = l_err(0);
            double ly   = l_err(1);
            double fxlx = fx*lx;
            double fyly = fy*ly;
            // - jac. wrt. KF pose
            Vector6d Jij_Piw = Vector6d::Zero();
            Jij_Piw << + gz2 * fxlx * gz,
                       + gz2 * fyly * gz,
                       - gz2 * ( fxlx*gx + fyly*gy ),
                       - gz2 * ( fxlx*gx*gy + fyly*gy*gy + fyly*gz*gz ),
                       + gz2 * ( fxlx*gx*gx + fxlx*gz*gz + fyly*gx*gy ),
                       + gz2 * ( fyly*gx*gz - fxlx*gy*gz );
            // - jac. wrt. LM
            Vector3d Jij_Pwj = Vector3d::Zero();
            Jij_Pwj << + gz2 * fxlx * gz,
                       + gz2 * fyly * gz,
                       - gz2 * ( fxlx*gx + fyly*gy );
            Jij_Pwj = Jij_Pwj.transpose() * Tiw.block(0,0,3,3) * l_err(0) / std::max(SlamConfig::homogTh(),l_err_norm);
            // end point
            gx   = Qwi(0);
            gy   = Qwi(1);
            gz   = Qwi(2);
            gz2  = gz*gz;
            gz2         = 1.0 / std::max(SlamConfig::homogTh(),gz2);
            // - jac. wrt. KF pose
            Vector6d Jij_Qiw = Vector6d::Zero();
            Jij_Qiw << + gz2 * fxlx * gz,
                       + gz2 * fyly * gz,
                       - gz2 * ( fxlx*gx + fyly*gy ),
                       - gz2 * ( fxlx*gx*gy + fyly*gy*gy + fyly*gz*gz ),
                       + gz2 * ( fxlx*gx*gx + fxlx*gz*gz + fyly*gx*gy ),
                       + gz2 * ( fyly*gx*gz - fxlx*gy*gz );
            // - jac. wrt. LM
            Vector3d Jij_Qwj = Vector3d::Zero();
            Jij_Qwj << + gz2 * fxlx * gz,
                       + gz2 * fyly * gz,
                       - gz2 * ( fxlx*gx + fyly*gy );
            Jij_Qwj = Jij_Qwj.transpose() * Tiw.block(0,0,3,3) * l_err(1) / std::max(SlamConfig::homogTh(),l_err_norm);
            // estimate Jacobian wrt KF pose
            Vector6d Jij_Tiw = Vector6d::Zero();
            Jij_Tiw = ( Jij_Piw * l_err(0) + Jij_Qiw * l_err(1) ) / std::max(SlamConfig::homogTh(),l_err_norm);
            // estimate Jacobian wrt LM
            Vector6d Jij_Lwj = Vector6d::Zero();
            Jij_Lwj.head(3) = Jij_Pwj;
            Jij_Lwj.tail(3) = Jij_Qwj;
            // if employing robust cost function
            double w  = 1.0;
            w = robustWeightCauchy(l_err_norm) ;

            // update hessian, gradient, and error
            MatrixXd Haux  = MatrixXd::Zero(3,6);
            int idx = 6 * kf_idx_loc;
            int jdx = 6*Nkf + 3*Npt + 6*lm_idx_loc;
            if( kf_idx_loc == -1 )
            {
                g.block(jdx,0,6,1) += Jij_Lwj * l_err_norm * w;
                err += l_err_norm * l_err_norm * w;
                H.block(jdx,jdx,6,6) += Jij_Lwj * Jij_Lwj.transpose() * w;
            }
            else
            {
                g.block(idx,0,6,1) += Jij_Tiw * l_err_norm * w;
                g.block(jdx,0,6,1) += Jij_Lwj * l_err_norm * w;
                err += l_err_norm * l_err_norm * w;
                Haux = Jij_Lwj * Jij_Tiw.transpose() * w;
                H.block(idx,idx,6,6) += Jij_Tiw * Jij_Tiw.transpose() * w;
                H.block(jdx,idx,6,6) += Haux;
                H.block(idx,jdx,6,6) += Haux.transpose();
                H.block(jdx,jdx,6,6) += Jij_Lwj * Jij_Lwj.transpose() * w;
            }
        }
    }
    err /= (Npt_obs+Nls_obs);

    // initial guess of lambda
    double Hmax = 0.0;
    for( int i = 0; i < N; i++)
    {
        if( H(i,i) > Hmax || H(i,i) < -Hmax )
            Hmax = fabs( H(i,i) );
    }
    lambda *= Hmax;

    // solve the first iteration
    for(int i = 0; i < N; i++)
        H(i,i) += lambda * H(i,i) ;
    H_ = H.sparseView();
    SimplicialLDLT< SparseMatrix<double> > solver1(H_);
    DX = solver1.solve( g );

    // update KFs
    for( int i = 0; i < Nkf; i++)
    {
        Matrix4d Tprev = expmap_se3( X.block(6*i,0,6,1) );
        Matrix4d Tcurr = Tprev * inverse_se3( expmap_se3( DX.block(6*i,0,6,1) ) );
        X.block(6*i,0,6,1) = logmap_se3( Tcurr );
    }
    // update point LMs
    for( int i = 6*Nkf; i < 6*Nkf+3*Npt; i++)
        X(i) += DX(i);
    // update line LMs
    for( int i = 6*Nkf+3*Npt; i < N; i++)
        X(i) += DX(i);

    // update error
    err_prev = err;

    // LM iterations
    //---------------------------------------------------------------------------------------------
    int iters;
    for( iters = 1; iters < max_iters; iters++)
    {
        // estimate hessian and gradient (reset)
        DX = VectorXd::Zero(N);
        g  = VectorXd::Zero(N);
        H  = MatrixXd::Zero(N,N);
        err = 0.0;        
        // - point observations
        for( vector<Vector6i>::iterator pt_it = pt_obs_list.begin(); pt_it != pt_obs_list.end(); pt_it++ )
        {
            int lm_idx_map = (*pt_it)(0);
            int lm_idx_loc = (*pt_it)(1);
            int lm_idx_obs = (*pt_it)(2);
            int kf_idx_map = (*pt_it)(3);
            int kf_idx_loc = (*pt_it)(4);
            if( map_points[lm_idx_map] != NULL && map_keyframes[kf_idx_map] != NULL)
            {
                // grab 3D LM (Xwj)
                Vector3d Xwj = X.block(6*Nkf+3*lm_idx_loc,0,3,1);
                // grab 6DoF KF (Tiw)
                Matrix4d Tiw;
                if( kf_idx_loc != -1 )
                    Tiw = expmap_se3( X.block( 6*kf_idx_loc,0,6,1 ) );
                else
                    Tiw = map_keyframes[kf_idx_map]->T_kf_w;
                // projection error
                Tiw = inverse_se3( Tiw );
                Vector3d Xwi   = Tiw.block(0,0,3,3) * Xwj + Tiw.block(0,3,3,1);
                Vector2d p_prj = cam->projection( Xwi );
                Vector2d p_obs = map_points[lm_idx_map]->obs_list[lm_idx_obs];
                Vector2d p_err    = p_obs - p_prj;
                double p_err_norm = p_err.norm();
                // useful variables
                double gx   = Xwi(0);
                double gy   = Xwi(1);
                double gz   = Xwi(2);
                double gz2  = gz*gz;
                gz2         = 1.0 / std::max(SlamConfig::homogTh(),gz2);
                double fx   = cam->getFx();
                double fy   = cam->getFy();
                double dx   = p_err(0);
                double dy   = p_err(1);
                double fxdx = fx*dx;
                double fydy = fy*dy;
                // estimate Jacobian wrt KF pose
                Vector6d Jij_Tiw = Vector6d::Zero();
                Jij_Tiw << + gz2 * fxdx * gz,
                           + gz2 * fydy * gz,
                           - gz2 * ( fxdx*gx + fydy*gy ),
                           - gz2 * ( fxdx*gx*gy + fydy*gy*gy + fydy*gz*gz ),
                           + gz2 * ( fxdx*gx*gx + fxdx*gz*gz + fydy*gx*gy ),
                           + gz2 * ( fydy*gx*gz - fxdx*gy*gz );
                Jij_Tiw = Jij_Tiw / std::max(SlamConfig::homogTh(),p_err_norm);
                // estimate Jacobian wrt LM
                Vector3d Jij_Xwj = Vector3d::Zero();
                Jij_Xwj << + gz2 * fxdx * gz,
                           + gz2 * fydy * gz,
                           - gz2 * ( fxdx*gx + fydy*gy );
                Jij_Xwj = Jij_Xwj.transpose() * Tiw.block(0,0,3,3) / std::max(SlamConfig::homogTh(),p_err_norm);
                // if employing robust cost function
                double w  = 1.0;
                double s2 = map_points[lm_idx_map]->sigma_list[lm_idx_obs];
                //double w = 1.0 / ( 1.0 + p_err_norm * p_err_norm * s2 );
                w = robustWeightCauchy(p_err_norm) ;

                // update hessian, gradient, and error
                MatrixXd Haux  = MatrixXd::Zero(3,6);
                int idx = 6 * kf_idx_loc;
                int jdx = 6*Nkf + 3*lm_idx_loc;
                if( kf_idx_loc == -1 )
                {
                    g.block(jdx,0,3,1) += Jij_Xwj * p_err_norm * w;
                    err += p_err_norm * p_err_norm * w;
                    H.block(jdx,jdx,3,3) += Jij_Xwj * Jij_Xwj.transpose() * w;
                }
                else
                {
                    g.block(idx,0,6,1) += Jij_Tiw * p_err_norm * w;
                    g.block(jdx,0,3,1) += Jij_Xwj * p_err_norm * w;
                    err += p_err_norm * p_err_norm * w;
                    Haux = Jij_Xwj * Jij_Tiw.transpose() * w ;
                    H.block(idx,idx,6,6) += Jij_Tiw * Jij_Tiw.transpose() * w ;
                    H.block(jdx,idx,3,6) += Haux;
                    H.block(idx,jdx,6,3) += Haux.transpose();
                    H.block(jdx,jdx,3,3) += Jij_Xwj * Jij_Xwj.transpose() * w ;
                }
            }
        }
        // - line segment observations
        for( vector<Vector6i>::iterator ls_it = ls_obs_list.begin(); ls_it != ls_obs_list.end(); ls_it++ )
        {
            int lm_idx_map = (*ls_it)(0);
            int lm_idx_loc = (*ls_it)(1);
            int lm_idx_obs = (*ls_it)(2);
            int kf_idx_map = (*ls_it)(3);
            int kf_idx_loc = (*ls_it)(4);
            if( map_lines[lm_idx_map] != NULL && map_keyframes[kf_idx_map] != NULL)
            {
                // grab 3D LM (Pwj and Qwj)
                Vector3d Pwj = X.block(6*Nkf+3*Npt+3*lm_idx_loc,0,3,1);
                Vector3d Qwj = X.block(6*Nkf+3*Npt+3*lm_idx_loc,0,3,1);
                // grab 6DoF KF (Tiw)
                Matrix4d Tiw   = map_keyframes[kf_idx_map]->T_kf_w;
                // projection error
                Tiw = inverse_se3( Tiw );
                Vector3d Pwi   = Tiw.block(0,0,3,3) * Pwj + Tiw.block(0,3,3,1);
                Vector3d Qwi   = Tiw.block(0,0,3,3) * Qwj + Tiw.block(0,3,3,1);
                Vector2d p_prj = cam->projection( Pwi );
                Vector2d q_prj = cam->projection( Qwi );
                Vector3d l_obs = map_lines[lm_idx_map]->obs_list[lm_idx_obs];
                Vector2d l_err;
                l_err(0) = l_obs(0) * p_prj(0) + l_obs(1) * p_prj(1) + l_obs(2);
                l_err(1) = l_obs(0) * q_prj(0) + l_obs(1) * q_prj(1) + l_obs(2);
                double l_err_norm = l_err.norm();
                // start point
                double gx   = Pwi(0);
                double gy   = Pwi(1);
                double gz   = Pwi(2);
                double gz2  = gz*gz;
                gz2         = 1.0 / std::max(0.0000001,gz2);
                double fx   = cam->getFx();
                double fy   = cam->getFy();
                double lx   = l_err(0);
                double ly   = l_err(1);
                double fxlx = fx*lx;
                double fyly = fy*ly;
                // - jac. wrt. KF pose
                Vector6d Jij_Piw = Vector6d::Zero();
                Jij_Piw << + gz2 * fxlx * gz,
                           + gz2 * fyly * gz,
                           - gz2 * ( fxlx*gx + fyly*gy ),
                           - gz2 * ( fxlx*gx*gy + fyly*gy*gy + fyly*gz*gz ),
                           + gz2 * ( fxlx*gx*gx + fxlx*gz*gz + fyly*gx*gy ),
                           + gz2 * ( fyly*gx*gz - fxlx*gy*gz );
                // - jac. wrt. LM
                Vector3d Jij_Pwj = Vector3d::Zero();
                Jij_Pwj << + gz2 * fxlx * gz,
                           + gz2 * fyly * gz,
                           - gz2 * ( fxlx*gx + fyly*gy );
                Jij_Pwj = Jij_Pwj.transpose() * Tiw.block(0,0,3,3) * l_err(0) / std::max(0.0000001,l_err_norm);
                // end point
                gx   = Qwi(0);
                gy   = Qwi(1);
                gz   = Qwi(2);
                gz2  = gz*gz;
                gz2         = 1.0 / std::max(0.0000001,gz2);
                // - jac. wrt. KF pose
                Vector6d Jij_Qiw = Vector6d::Zero();
                Jij_Qiw << + gz2 * fxlx * gz,
                           + gz2 * fyly * gz,
                           - gz2 * ( fxlx*gx + fyly*gy ),
                           - gz2 * ( fxlx*gx*gy + fyly*gy*gy + fyly*gz*gz ),
                           + gz2 * ( fxlx*gx*gx + fxlx*gz*gz + fyly*gx*gy ),
                           + gz2 * ( fyly*gx*gz - fxlx*gy*gz );
                // - jac. wrt. LM
                Vector3d Jij_Qwj = Vector3d::Zero();
                Jij_Qwj << + gz2 * fxlx * gz,
                           + gz2 * fyly * gz,
                           - gz2 * ( fxlx*gx + fyly*gy );
                Jij_Qwj = Jij_Qwj.transpose() * Tiw.block(0,0,3,3) * l_err(1) / std::max(0.0000001,l_err_norm);
                // estimate Jacobian wrt KF pose
                Vector6d Jij_Tiw = Vector6d::Zero();
                Jij_Tiw = ( Jij_Piw * l_err(0) + Jij_Qiw * l_err(1) ) / std::max(0.0000001,l_err_norm);
                // estimate Jacobian wrt LM
                Vector6d Jij_Lwj = Vector6d::Zero();
                Jij_Lwj.head(3) = Jij_Pwj;
                Jij_Lwj.tail(3) = Jij_Qwj;
                // if employing robust cost function
                double w  = 1.0;
                w = robustWeightCauchy(l_err_norm) ;

                // update hessian, gradient, and error
                MatrixXd Haux  = MatrixXd::Zero(3,6);
                int idx = 6 * kf_idx_loc;
                int jdx = 6*Nkf + 3*Npt + 6*lm_idx_loc;
                if( kf_idx_loc == -1 )
                {
                    g.block(jdx,0,6,1) += Jij_Lwj * l_err_norm * w;
                    err += l_err_norm * l_err_norm * w;
                    H.block(jdx,jdx,6,6) += Jij_Lwj * Jij_Lwj.transpose() * w;
                }
                else
                {
                    g.block(idx,0,6,1) += Jij_Tiw * l_err_norm * w;
                    g.block(jdx,0,6,1) += Jij_Lwj * l_err_norm * w;
                    err += l_err_norm * l_err_norm * w;
                    Haux = Jij_Lwj * Jij_Tiw.transpose() * w;
                    H.block(idx,idx,6,6) += Jij_Tiw * Jij_Tiw.transpose() * w;
                    H.block(jdx,idx,6,6) += Haux;
                    H.block(idx,jdx,6,6) += Haux.transpose();
                    H.block(jdx,jdx,6,6) += Jij_Lwj * Jij_Lwj.transpose() * w;
                }
            }
        }
        err /= (Npt+Nls);
        // if the difference is very small stop
        if( abs(err-err_prev) < Config::minErrorChange() || err < Config::minError() )
            break;
        // add lambda to hessian
        for(int i = 0; i < N; i++)
            H(i,i) += lambda * H(i,i) ;
        // solve iteration
        H_ = H.sparseView();
        SimplicialLDLT< SparseMatrix<double> > solver1(H_);
        DX = solver1.solve( g );

        // update lambda
        if( err > err_prev ){
            lambda /= lambda_k;
        }
        else
        {
            lambda *= lambda_k;
            // update KFs
            for( int i = 0; i < Nkf; i++)
            {
                Matrix4d Tprev = expmap_se3( X.block(6*i,0,6,1) );
                Matrix4d Tcurr = Tprev * inverse_se3( expmap_se3( DX.block(6*i,0,6,1) ) );
                X.block(6*i,0,6,1) = logmap_se3( Tcurr );
            }
            // update point LMs
            for( int i = 6*Nkf; i < 6*Nkf+3*Npt; i++)
                X(i) += DX(i);
            // update line LMs
            for( int i = 6*Nkf+3*Npt; i < N; i++)
                X(i) += DX(i);
        }
        // if the parameter change is small stop
        if( DX.norm() < Config::minErrorChange() )
            break;
        // update previous values
        err_prev = err;

    }

    if( vo_status != VO_INSERTING_KF )
    {

    m_insert_kf.lock();

    // Update KFs and LMs
    //---------------------------------------------------------------------------------------------
    // update KFs
    for( int i = 0; i < Nkf; i++)
    {
        Matrix4d Test = expmap_se3( X.block( 6*i,0,6,1 ) );
        map_keyframes[ kf_list[i] ]->T_kf_w = Test;
    }
    // update point LMs
    for( int i = 0; i < Npt; i++)
    {

        Vector3d DX = X.block(6*Nkf+3*i,0,3,1) - map_points[ pt_list[i] ]->point3D;
        if( DX.norm() > 0.01 )
            map_points[ pt_list[i] ]->inlier = false;

        map_points[ pt_list[i] ]->point3D(0) = X(6*Nkf+3*i);
        map_points[ pt_list[i] ]->point3D(1) = X(6*Nkf+3*i+1);
        map_points[ pt_list[i] ]->point3D(2) = X(6*Nkf+3*i+2);

    }
    // update line segment LMs
    for( int i = 0; i < Nls; i++)
    {

        Vector6d DX = X.block(6*Nkf+3*Npt+6*i,0,6,1) - map_lines[ ls_list[i] ]->line3D;
        if( DX.norm() > 0.01 )
            map_lines[ ls_list[i] ]->inlier = false;

        map_lines[ ls_list[i] ]->line3D(0) = X(6*Nkf+3*Npt+6*i);
        map_lines[ ls_list[i] ]->line3D(1) = X(6*Nkf+3*Npt+6*i+1);
        map_lines[ ls_list[i] ]->line3D(2) = X(6*Nkf+3*Npt+6*i+2);
        map_lines[ ls_list[i] ]->line3D(3) = X(6*Nkf+3*Npt+6*i+3);
        map_lines[ ls_list[i] ]->line3D(4) = X(6*Nkf+3*Npt+6*i+4);
        map_lines[ ls_list[i] ]->line3D(5) = X(6*Nkf+3*Npt+6*i+5);

    }

    // Remove bad observations
    //---------------------------------------------------------------------------------------------
    for( vector<Vector6i>::reverse_iterator pt_it = pt_obs_list.rbegin(); pt_it != pt_obs_list.rend(); ++pt_it )
    {
        if( (*pt_it)(5) == -1 )
        {
            int lm_idx_map = (*pt_it)(0);
            int lm_idx_obs = (*pt_it)(2);
            if( map_points[lm_idx_map] != NULL )
            {
                int kf_obs = map_points[lm_idx_map]->kf_obs_list[lm_idx_obs];
                // remove observations from map_points
                if( map_points[lm_idx_map]->obs_list.size() > 1 )
                {
                    // if it is the first observation, update it from map_points_kf_idx
                    if( lm_idx_obs == 0 )
                    {
                        // delete observation from map_points_kf_idx
                        for( auto it = map_points_kf_idx.at(kf_obs).begin(); it != map_points_kf_idx.at(kf_obs).end(); it++)
                        {
                            if( (*it) == lm_idx_map )
                            {
                                int new_kf_base = map_points[(*it)]->kf_obs_list[1];
                                map_points_kf_idx.at(new_kf_base).push_back( (*it) );
                                break;
                            }
                        }
                    }
                    // remove observations from map points
                    map_points[lm_idx_map]->desc_list.erase( map_points[lm_idx_map]->desc_list.begin() + lm_idx_obs );
                    map_points[lm_idx_map]->obs_list.erase( map_points[lm_idx_map]->obs_list.begin() + lm_idx_obs );
                    map_points[lm_idx_map]->dir_list.erase( map_points[lm_idx_map]->dir_list.begin() + lm_idx_obs );
                    map_points[lm_idx_map]->kf_obs_list.erase( map_points[lm_idx_map]->kf_obs_list.begin() + lm_idx_obs );
                    // remove idx from KeyFrame stereo points
                    for(vector<PointFeature*>::iterator st_pt = map_keyframes[kf_obs]->stereo_frame->stereo_pt.begin();
                        st_pt != map_keyframes[kf_obs]->stereo_frame->stereo_pt.end(); st_pt++ )
                    {
                        if( (*st_pt)->idx == lm_idx_map )
                        {
                            (*st_pt)->idx = -1;
                            st_pt = map_keyframes[kf_obs]->stereo_frame->stereo_pt.end()-1;
                        }
                    }
                    // update main descriptor and direction
                    map_points[lm_idx_map]->updateAverageDescDir();
                    // update graphs
                    for( int i = 0; i < map_points[lm_idx_map]->kf_obs_list.size(); i++ )
                    {
                        int idx = map_points[lm_idx_map]->kf_obs_list[i];
                        if( kf_obs != idx )
                        {
                            full_graph[kf_obs][idx]--;
                            full_graph[idx][kf_obs]--;
                        }
                    }
                }
                else
                    map_points[lm_idx_map]->inlier = false;
            }
        }
    }

    for( vector<Vector6i>::reverse_iterator ls_it = ls_obs_list.rbegin(); ls_it != ls_obs_list.rend(); ++ls_it )
    {
        if( (*ls_it)(5) == -1 )
        {
            int lm_idx_map = (*ls_it)(0);
            int lm_idx_obs = (*ls_it)(2);
            if( map_lines[lm_idx_map] != NULL )
            {
                int kf_obs = map_lines[lm_idx_map]->kf_obs_list[lm_idx_obs];
                // remove observations from map_points
                if( map_lines[lm_idx_map]->obs_list.size() > 1 )
                {
                    // if it is the first observation, update it from map_points_kf_idx
                    if( lm_idx_obs == 0 )
                    {
                        // delete observation from map_points_kf_idx
                        for( auto it = map_lines_kf_idx.at(kf_obs).begin(); it != map_lines_kf_idx.at(kf_obs).end(); it++)
                        {
                            if( (*it) == lm_idx_map )
                            {
                                int new_kf_base = map_lines[(*it)]->kf_obs_list[1];
                                map_lines_kf_idx.at(new_kf_base).push_back( (*it) );
                                break;
                            }
                        }
                    }
                    // remove observations
                    map_lines[lm_idx_map]->desc_list.erase( map_lines[lm_idx_map]->desc_list.begin() + lm_idx_obs );
                    map_lines[lm_idx_map]->obs_list.erase( map_lines[lm_idx_map]->obs_list.begin() + lm_idx_obs );
                    map_lines[lm_idx_map]->pts_list.erase( map_lines[lm_idx_map]->pts_list.begin() + lm_idx_obs );
                    map_lines[lm_idx_map]->dir_list.erase( map_lines[lm_idx_map]->dir_list.begin() + lm_idx_obs );
                    map_lines[lm_idx_map]->kf_obs_list.erase( map_lines[lm_idx_map]->kf_obs_list.begin() + lm_idx_obs );
                    // remove idx from KeyFrame stereo lines
                    for(vector<LineFeature*>::iterator st_ls = map_keyframes[kf_obs]->stereo_frame->stereo_ls.begin();
                        st_ls != map_keyframes[kf_obs]->stereo_frame->stereo_ls.end(); st_ls++ )
                    {
                        if( (*st_ls)->idx == lm_idx_map )
                        {
                            (*st_ls)->idx = -1;
                            st_ls = map_keyframes[kf_obs]->stereo_frame->stereo_ls.end()-1;
                        }
                    }
                    // update main descriptor and direction
                    map_lines[lm_idx_map]->updateAverageDescDir();
                    // update graphs
                    for( int i = 0; i < map_lines[lm_idx_map]->kf_obs_list.size(); i++ )
                    {
                        int idx = map_lines[lm_idx_map]->kf_obs_list[i];
                        if( kf_obs != idx )
                        {
                            full_graph[kf_obs][idx]--;
                            full_graph[idx][kf_obs]--;
                        }
                    }
                }
                else
                    map_lines[lm_idx_map]->inlier = false;
            }

        }
    }

    m_insert_kf.unlock();

    }
    else
        return -1;

    return 0;

}

// -----------------------------------------------------------------------------------------------------------------------------
// Global Bundle Adjustment functions
// -----------------------------------------------------------------------------------------------------------------------------

void MapHandler::globalBundleAdjustment()
{

    vector<double> X_aux;

    // create list of keyframes
    vector<int> kf_list;
    for( vector<KeyFrame*>::iterator kf_it = map_keyframes.begin(); kf_it != map_keyframes.end(); kf_it++)
    {
        if( (*kf_it)!= NULL )
        {
            if( (*kf_it)->kf_idx != 0 )
            {
                Vector6d pose_aux = (*kf_it)->x_kf_w;
                for(int i = 0; i < 6; i++)
                    X_aux.push_back( pose_aux(i) );
                kf_list.push_back( (*kf_it)->kf_idx );
            }
        }
    }

    // create list of point landmarks
    vector<Vector6i> pt_obs_list;
    vector<int> pt_list;
    int lm_local_idx = 0;
    for( vector<MapPoint*>::iterator pt_it = map_points.begin(); pt_it != map_points.end(); pt_it++)
    {
        if( (*pt_it)!= NULL )
        {
            Vector3d point_aux = (*pt_it)->point3D;
            for(int i = 0; i < 3; i++)
                X_aux.push_back( point_aux(i) );
            // gather all observations
            for( int i = 0; i < (*pt_it)->obs_list.size(); i++)
            {
                Vector6i obs_aux;
                obs_aux(0) = (*pt_it)->idx; // LM idx
                obs_aux(1) = lm_local_idx;  // LM local idx
                obs_aux(2) = i;             // LM obs idx
                int kf_obs_list_ = (*pt_it)->kf_obs_list[i];
                obs_aux(3) = kf_obs_list_;  // KF idx
                obs_aux(4) = -1;            // KF local idx (-1 if not local)
                obs_aux(5) = 1;             // 1 if the observation is an inlier
                for( int j = 0; j < kf_list.size(); j++ )
                {
                    if( kf_list[j] == kf_obs_list_ )
                    {
                        obs_aux(4) = j;
                        break;
                    }
                }
                pt_obs_list.push_back( obs_aux );
            }
            lm_local_idx++;
            // pt_list
            pt_list.push_back( (*pt_it)->idx );
        }
    }

    // create list of line segment landmarks
    vector<Vector6i> ls_obs_list;
    vector<int> ls_list;
    lm_local_idx = 0;
    for( vector<MapLine*>::iterator ls_it = map_lines.begin(); ls_it != map_lines.end(); ls_it++)
    {
        if( (*ls_it)!= NULL )
        {
            Vector6d line_aux = (*ls_it)->line3D;
            for(int i = 0; i < 6; i++)
                X_aux.push_back( line_aux(i) );
            // gather all observations
            for( int i = 0; i < (*ls_it)->obs_list.size(); i++)
            {
                Vector6i obs_aux;
                obs_aux(0) = (*ls_it)->idx; // LM idx
                obs_aux(1) = lm_local_idx;  // LM local idx
                obs_aux(2) = i;             // LM obs idx
                int kf_obs_list_ = (*ls_it)->kf_obs_list[i];
                obs_aux(3) = kf_obs_list_;  // KF idx
                obs_aux(4) = -1;            // KF local idx (-1 if not local)
                obs_aux(5) = 1;             // 1 if the observation is an inlier
                for( int j = 0; j < kf_list.size(); j++ )
                {
                    if( kf_list[j] == kf_obs_list_ )
                    {
                        obs_aux(4) = j;
                        break;
                    }
                }
                ls_obs_list.push_back( obs_aux );
            }
            lm_local_idx++;
            // ls_list
            ls_list.push_back( (*ls_it)->idx );
        }
    }

    // Levenberg-Marquardt optimization
    levMarquardtOptimizationGBA(X_aux,kf_list,pt_list,ls_list,pt_obs_list,ls_obs_list);

    // -------------------------------------------------------------------------------------------------------------------

    // Recent map LMs culling (implement filters for line segments, which seems to be unaccurate)
    // Recent KFs culling
}

void MapHandler::levMarquardtOptimizationGBA( vector<double> X_aux, vector<int> kf_list, vector<int> pt_list, vector<int> ls_list, vector<Vector6i> pt_obs_list, vector<Vector6i> ls_obs_list  )
{

    // create Levenberg-Marquardt variables
    int    Nkf = kf_list.size();
    int      N = X_aux.size();
    VectorXd DX, X(N), gdense(N);
    SparseVector<double> g(N);
    SparseMatrix<double> H(N,N);
    for(int i = 0; i < N; i++)
        X.coeffRef(i) = X_aux[i];
    H.reserve( VectorXi::Constant(N,5000) );

    // create Levenberg-Marquardt parameters
    double err, err_prev = 999999999.9;
    double lambda = SlamConfig::lambdaLbaLM(), lambda_k = SlamConfig::lambdaLbaK();

    // estimate H and g to precalculate lambda
    //---------------------------------------------------------------------------------------------
    // point observations
    int Npt = 0, Npt_obs = 0;
    if( pt_obs_list.size() != 0 )
        Npt = pt_obs_list.back()(1)+1;
    for( vector<Vector6i>::iterator pt_it = pt_obs_list.begin(); pt_it != pt_obs_list.end(); pt_it++ )
    {
        int lm_idx_map = (*pt_it)(0);
        int lm_idx_loc = (*pt_it)(1);
        int lm_idx_obs = (*pt_it)(2);
        int kf_idx_map = (*pt_it)(3);
        int kf_idx_loc = (*pt_it)(4);
        if( map_points[lm_idx_map] != NULL && map_keyframes[kf_idx_map] != NULL)
        {
            // grab 3D LM (Xwj)
            Vector3d Xwj   = map_points[lm_idx_map]->point3D;
            // grab 6DoF KF (Tiw)
            Matrix4d Tiw   = map_keyframes[kf_idx_map]->T_kf_w;
            // projection error
            Tiw = inverse_se3( Tiw );
            Vector3d Xwi   = Tiw.block(0,0,3,3) * Xwj + Tiw.block(0,3,3,1);
            Vector2d p_prj = cam->projection( Xwi );
            Vector2d p_obs = map_points[lm_idx_map]->obs_list[lm_idx_obs];
            Vector2d p_err    = p_obs - p_prj;
            double p_err_norm = p_err.norm();
            double gx   = Xwi(0);
            double gy   = Xwi(1);
            double gz   = Xwi(2);
            double gz2  = gz*gz;
            gz2         = 1.0 / std::max(SlamConfig::homogTh(),gz2);
            double fx   = cam->getFx();
            double fy   = cam->getFy();
            double dx   = p_err(0);
            double dy   = p_err(1);
            double fxdx = fx*dx;
            double fydy = fy*dy;
            // estimate Jacobian wrt KF pose
            Vector6d Jij_Tiw = Vector6d::Zero();
            Jij_Tiw << + gz2 * fxdx * gz,
                       + gz2 * fydy * gz,
                       - gz2 * ( fxdx*gx + fydy*gy ),
                       - gz2 * ( fxdx*gx*gy + fydy*gy*gy + fydy*gz*gz ),
                       + gz2 * ( fxdx*gx*gx + fxdx*gz*gz + fydy*gx*gy ),
                       + gz2 * ( fydy*gx*gz - fxdx*gy*gz );
            Jij_Tiw = Jij_Tiw / std::max(SlamConfig::homogTh(),p_err_norm);
            // estimate Jacobian wrt LM
            Vector3d Jij_Xwj = Vector3d::Zero();
            Jij_Xwj << + gz2 * fxdx * gz,
                       + gz2 * fydy * gz,
                       - gz2 * ( fxdx*gx + fydy*gy );
            Jij_Xwj = Jij_Xwj.transpose() * Tiw.block(0,0,3,3) / std::max(SlamConfig::homogTh(),p_err_norm);
            // if employing robust cost function
            double w  = 1.0;
            double s2 = map_points[lm_idx_map]->sigma_list[lm_idx_obs];
            //double w = 1.0 / ( 1.0 + p_err_norm * p_err_norm * s2 );
            w = robustWeightCauchy(p_err_norm) ;

            // update hessian, gradient, and error
            int idx = 6 * kf_idx_loc;
            int jdx = 6*Nkf + 3*lm_idx_loc;
            if( kf_idx_loc == -1 )
            {
                err += p_err_norm * p_err_norm * w;
                Vector3d gi;
                Matrix3d Hjj;
                Hjj = Jij_Xwj * Jij_Xwj.transpose() * w;
                gi  = Jij_Xwj * p_err_norm * w;
                for(int i = 0; i < 3; i++)
                {
                    g.coeffRef(jdx+i) += gi(i);
                    for(int j = 0; j < 3; j++)
                        H.coeffRef(i+jdx,j+jdx) += Hjj(i,j);
                }
            }
            else
            {
                err += p_err_norm * p_err_norm * w;
                Vector3d gj;
                Vector6d gi;
                Matrix6d Hii;
                gi = Jij_Tiw * p_err_norm * w;
                gj = Jij_Xwj * p_err_norm * w;
                Hii = Jij_Tiw * Jij_Tiw.transpose() * w;
                for(int i = 0; i < 6; i++)
                {
                    g.coeffRef(i+idx) += gi(i);
                    for(int j = 0; j < 6; j++)
                        H.coeffRef(i+idx,j+idx) += Hii(i,j);
                }
                Matrix3d Hjj;
                Hjj = Jij_Xwj * Jij_Xwj.transpose() * w;
                for(int i = 0; i < 3; i++)
                {
                    g.coeffRef(i+jdx) += gj(i);
                    for(int j = 0; j < 3; j++)
                        H.coeffRef(i+jdx,j+jdx) += Hjj(i,j);
                }
                MatrixXd Hij  = MatrixXd::Zero(3,6);
                Hij = Jij_Xwj * Jij_Tiw.transpose() * w;
                for(int i = 0; i < 3; i++)
                {
                    for(int j = 0; j < 6; j++)
                    {
                        H.coeffRef(i+jdx,j+idx) += Hij(i,j);
                        H.coeffRef(j+idx,i+jdx) += Hij(i,j);
                    }
                }
            }
        }
    }
    // line segment observations
    int Nls = 0, Nls_obs = 0;
    if( ls_obs_list.size() != 0 )
        Nls = ls_obs_list.back()(1)+1;
    for( vector<Vector6i>::iterator ls_it = ls_obs_list.begin(); ls_it != ls_obs_list.end(); ls_it++ )
    {
        int lm_idx_map = (*ls_it)(0);
        int lm_idx_loc = (*ls_it)(1);
        int lm_idx_obs = (*ls_it)(2);
        int kf_idx_map = (*ls_it)(3);
        int kf_idx_loc = (*ls_it)(4);
        if( map_lines[lm_idx_map] != NULL && map_keyframes[kf_idx_map] != NULL)
        {
            // grab 3D LM (Pwj and Qwj)
            Vector3d Pwj   = map_lines[lm_idx_map]->line3D.head(3);
            Vector3d Qwj   = map_lines[lm_idx_map]->line3D.tail(3);
            // grab 6DoF KF (Tiw)
            Matrix4d Tiw   = map_keyframes[kf_idx_map]->T_kf_w;
            // projection error
            Tiw = inverse_se3( Tiw );
            Vector3d Pwi   = Tiw.block(0,0,3,3) * Pwj + Tiw.block(0,3,3,1);
            Vector3d Qwi   = Tiw.block(0,0,3,3) * Qwj + Tiw.block(0,3,3,1);
            Vector2d p_prj = cam->projection( Pwi );
            Vector2d q_prj = cam->projection( Qwi );
            Vector3d l_obs = map_lines[lm_idx_map]->obs_list[lm_idx_obs];
            Vector2d l_err;
            l_err(0) = l_obs(0) * p_prj(0) + l_obs(1) * p_prj(1) + l_obs(2);
            l_err(1) = l_obs(0) * q_prj(0) + l_obs(1) * q_prj(1) + l_obs(2);
            double l_err_norm = l_err.norm();
            // start point
            double gx   = Pwi(0);
            double gy   = Pwi(1);
            double gz   = Pwi(2);
            double gz2  = gz*gz;
            gz2         = 1.0 / std::max(SlamConfig::homogTh(),gz2);
            double fx   = cam->getFx();
            double fy   = cam->getFy();
            double lx   = l_err(0);
            double ly   = l_err(1);
            double fxlx = fx*lx;
            double fyly = fy*ly;
            // - jac. wrt. KF pose
            Vector6d Jij_Piw = Vector6d::Zero();
            Jij_Piw << + gz2 * fxlx * gz,
                       + gz2 * fyly * gz,
                       - gz2 * ( fxlx*gx + fyly*gy ),
                       - gz2 * ( fxlx*gx*gy + fyly*gy*gy + fyly*gz*gz ),
                       + gz2 * ( fxlx*gx*gx + fxlx*gz*gz + fyly*gx*gy ),
                       + gz2 * ( fyly*gx*gz - fxlx*gy*gz );
            // - jac. wrt. LM
            Vector3d Jij_Pwj = Vector3d::Zero();
            Jij_Pwj << + gz2 * fxlx * gz,
                       + gz2 * fyly * gz,
                       - gz2 * ( fxlx*gx + fyly*gy );
            Jij_Pwj = Jij_Pwj.transpose() * Tiw.block(0,0,3,3) * l_err(0) / std::max(SlamConfig::homogTh(),l_err_norm);
            // end point
            gx   = Qwi(0);
            gy   = Qwi(1);
            gz   = Qwi(2);
            gz2  = gz*gz;
            gz2  = 1.0 / std::max(SlamConfig::homogTh(),gz2);
            // - jac. wrt. KF pose
            Vector6d Jij_Qiw = Vector6d::Zero();
            Jij_Qiw << + gz2 * fxlx * gz,
                       + gz2 * fyly * gz,
                       - gz2 * ( fxlx*gx + fyly*gy ),
                       - gz2 * ( fxlx*gx*gy + fyly*gy*gy + fyly*gz*gz ),
                       + gz2 * ( fxlx*gx*gx + fxlx*gz*gz + fyly*gx*gy ),
                       + gz2 * ( fyly*gx*gz - fxlx*gy*gz );
            // - jac. wrt. LM
            Vector3d Jij_Qwj = Vector3d::Zero();
            Jij_Qwj << + gz2 * fxlx * gz,
                       + gz2 * fyly * gz,
                       - gz2 * ( fxlx*gx + fyly*gy );
            Jij_Qwj = Jij_Qwj.transpose() * Tiw.block(0,0,3,3) * l_err(1) / std::max(SlamConfig::homogTh(),l_err_norm);
            // estimate Jacobian wrt KF pose
            Vector6d Jij_Tiw = Vector6d::Zero();
            Jij_Tiw = ( Jij_Piw * l_err(0) + Jij_Qiw * l_err(1) ) / std::max(SlamConfig::homogTh(),l_err_norm);
            // estimate Jacobian wrt LM
            Vector6d Jij_Lwj = Vector6d::Zero();
            Jij_Lwj.head(3) = Jij_Pwj;
            Jij_Lwj.tail(3) = Jij_Qwj;
            // if employing robust cost function
            double w  = 1.0;
            w = robustWeightCauchy(l_err_norm) ;
            // update hessian, gradient, and error
            VectorXd gi = VectorXd::Zero(N), gj = VectorXd::Zero(N);
            int idx = 6 * kf_idx_loc;
            int jdx = 6*Nkf + 3*Npt + 6*lm_idx_loc;
            if( kf_idx_loc == -1 )
            {
                gj = Jij_Lwj * l_err_norm * w;
                err += l_err_norm * l_err_norm * w;
                Matrix6d Hjj;
                Hjj = Jij_Lwj * Jij_Lwj.transpose() * w;
                for(int i = 0; i < 6; i++)
                {
                    g.coeffRef(jdx+i) += gj(i);
                    for(int j = 0; j < 6; j++)
                        H.coeffRef(i+jdx,j+jdx) += Hjj(i,j);
                }
            }
            else
            {
                gi = Jij_Tiw * l_err_norm * w;
                gj = Jij_Lwj * l_err_norm * w;
                err += l_err_norm * l_err_norm * w;
                Matrix6d Hii, Hjj;
                MatrixXd Hij   = MatrixXd::Zero(3,6);
                Hii = Jij_Tiw * Jij_Tiw.transpose() * w;
                Hjj = Jij_Lwj * Jij_Lwj.transpose() * w;
                Hij = Jij_Lwj * Jij_Tiw.transpose() * w;
                for(int i = 0; i < 6; i++)
                {
                    g.coeffRef(i+idx) += gi(i);
                    g.coeffRef(i+jdx) += gj(i);
                    for(int j = 0; j < 6; j++)
                    {
                        H.coeffRef(idx+i,idx+j) += Hii(i,j);
                        H.coeffRef(jdx+i,jdx+j) += Hjj(i,j);
                        H.coeffRef(idx+i,jdx+j) += Hij(i,j);
                        H.coeffRef(jdx+j,idx+i) += Hij(i,j);
                    }
                }
            }
        }
    }
    err /= (Npt_obs+Nls_obs);

    // initial guess of lambda
    int Hmax = 0.0;
    for( int i = 0; i < N; i++)
    {
        if( H.coeffRef(i,i) > Hmax || H.coeffRef(i,i) < -Hmax )
            Hmax = fabs( H.coeffRef(i,i) );
    }
    lambda *= Hmax;
    // solve the first iteration
    H.makeCompressed();
    for(int i = 0; i < N; i++)
    {
        H.coeffRef(i,i) += lambda * H.coeffRef(i,i);
        gdense(i) = g.coeffRef(i);
    }
    SimplicialLDLT< SparseMatrix<double> > solver1(H);
    DX = solver1.solve( gdense );
    // update KFs
    for( int i = 0; i < Nkf; i++)
    {
        Matrix4d Tprev = expmap_se3( X.block(6*i,0,6,1) );
        Matrix4d Tcurr = Tprev * inverse_se3( expmap_se3( DX.block(6*i,0,6,1) ) );       
        X.block(6*i,0,6,1) = logmap_se3( Tcurr );
    }
    // update LMs
    for( int i = 6*Nkf; i < N; i++)
        X(i) += DX(i);
    // update error
    err_prev = err;

    // LM iterations
    //---------------------------------------------------------------------------------------------
    int iters;
    for( iters = 1; iters < SlamConfig::maxItersLba(); iters++ )
    {
        // estimate hessian and gradient (reset)
        DX.setZero();
        g.setZero();
        H.setZero();
        H.reserve( VectorXi::Constant(N,5000) );
        err = 0.0;
        // - point observations
        for( vector<Vector6i>::iterator pt_it = pt_obs_list.begin(); pt_it != pt_obs_list.end(); pt_it++ )
        {
            int lm_idx_map = (*pt_it)(0);
            int lm_idx_loc = (*pt_it)(1);
            int lm_idx_obs = (*pt_it)(2);
            int kf_idx_map = (*pt_it)(3);
            int kf_idx_loc = (*pt_it)(4);
            if( map_points[lm_idx_map] != NULL && map_keyframes[kf_idx_map] != NULL)
            {
                // grab 3D LM (Xwj)
                Vector3d Xwj = X.block(6*Nkf+3*lm_idx_loc,0,3,1);
                // grab 6DoF KF (Tiw)
                Matrix4d Tiw;
                if( kf_idx_loc != -1 )
                    Tiw = expmap_se3( X.block( 6*kf_idx_loc,0,6,1 ) );
                else
                    Tiw = map_keyframes[kf_idx_map]->T_kf_w;
                // projection error
                Tiw = inverse_se3( Tiw );
                Vector3d Xwi   = Tiw.block(0,0,3,3) * Xwj + Tiw.block(0,3,3,1);
                Vector2d p_prj = cam->projection( Xwi );
                Vector2d p_obs = map_points[lm_idx_map]->obs_list[lm_idx_obs];
                Vector2d p_err    = p_obs - p_prj;
                double p_err_norm = p_err.norm();
                // useful variables
                double gx   = Xwi(0);
                double gy   = Xwi(1);
                double gz   = Xwi(2);
                double gz2  = gz*gz;
                gz2         = 1.0 / std::max(SlamConfig::homogTh(),gz2);
                double fx   = cam->getFx();
                double fy   = cam->getFy();
                double dx   = p_err(0);
                double dy   = p_err(1);
                double fxdx = fx*dx;
                double fydy = fy*dy;
                // estimate Jacobian wrt KF pose
                Vector6d Jij_Tiw = Vector6d::Zero();
                Jij_Tiw << + gz2 * fxdx * gz,
                           + gz2 * fydy * gz,
                           - gz2 * ( fxdx*gx + fydy*gy ),
                           - gz2 * ( fxdx*gx*gy + fydy*gy*gy + fydy*gz*gz ),
                           + gz2 * ( fxdx*gx*gx + fxdx*gz*gz + fydy*gx*gy ),
                           + gz2 * ( fydy*gx*gz - fxdx*gy*gz );
                Jij_Tiw = Jij_Tiw / std::max(SlamConfig::homogTh(),p_err_norm);
                // estimate Jacobian wrt LM
                Vector3d Jij_Xwj = Vector3d::Zero();
                Jij_Xwj << + gz2 * fxdx * gz,
                           + gz2 * fydy * gz,
                           - gz2 * ( fxdx*gx + fydy*gy );
                Jij_Xwj = Jij_Xwj.transpose() * Tiw.block(0,0,3,3) / std::max(SlamConfig::homogTh(),p_err_norm);
                // if employing robust cost function
                double  w = 1.0;
                w = robustWeightCauchy(p_err_norm) ;
                // update hessian, gradient, and error
                int idx = 6 * kf_idx_loc;
                int jdx = 6*Nkf + 3*lm_idx_loc;
                if( kf_idx_loc == -1 )
                {
                    err += p_err_norm * p_err_norm * w;
                    Vector3d gi;
                    Matrix3d Hjj;
                    Hjj = Jij_Xwj * Jij_Xwj.transpose() * w;
                    gi  = Jij_Xwj * p_err_norm * w;
                    for(int i = 0; i < 3; i++)
                    {
                        g.coeffRef(jdx+i) += gi(i);
                        for(int j = 0; j < 3; j++)
                            H.coeffRef(i+jdx,j+jdx) += Hjj(i,j);
                    }
                }
                else
                {
                    err += p_err_norm * p_err_norm * w;
                    Vector3d gj;
                    Vector6d gi;
                    Matrix6d Hii;
                    gi = Jij_Tiw * p_err_norm * w;
                    gj = Jij_Xwj * p_err_norm * w;
                    Hii = Jij_Tiw * Jij_Tiw.transpose() * w;
                    for(int i = 0; i < 6; i++)
                    {
                        g.coeffRef(i+idx) += gi(i);
                        for(int j = 0; j < 6; j++)
                            H.coeffRef(i+idx,j+idx) += Hii(i,j);
                    }
                    Matrix3d Hjj;
                    Hjj = Jij_Xwj * Jij_Xwj.transpose() * w;
                    for(int i = 0; i < 3; i++)
                    {
                        g.coeffRef(i+jdx) += gj(i);
                        for(int j = 0; j < 3; j++)
                            H.coeffRef(i+jdx,j+jdx) += Hjj(i,j);
                    }
                    MatrixXd Hij  = MatrixXd::Zero(3,6);
                    Hij = Jij_Xwj * Jij_Tiw.transpose() * w;
                    for(int i = 0; i < 3; i++)
                    {
                        for(int j = 0; j < 6; j++)
                        {
                            H.coeffRef(i+jdx,j+idx) += Hij(i,j);
                            H.coeffRef(j+idx,i+jdx) += Hij(i,j);
                        }
                    }
                }
            }
        }
        // - line segment observations
        for( vector<Vector6i>::iterator ls_it = ls_obs_list.begin(); ls_it != ls_obs_list.end(); ls_it++ )
        {
            int lm_idx_map = (*ls_it)(0);
            int lm_idx_loc = (*ls_it)(1);
            int lm_idx_obs = (*ls_it)(2);
            int kf_idx_map = (*ls_it)(3);
            int kf_idx_loc = (*ls_it)(4);
            if( map_lines[lm_idx_map] != NULL && map_keyframes[kf_idx_map] != NULL)
            {
                // grab 3D LM (Pwj and Qwj)
                //Vector3d Pwj   = map_lines[lm_idx_map]->line3D.head(3);
                //Vector3d Qwj   = map_lines[lm_idx_map]->line3D.tail(3);
                Vector3d Pwj = X.block(6*Nkf+3*Npt+3*lm_idx_loc,0,3,1);
                Vector3d Qwj = X.block(6*Nkf+3*Npt+3*lm_idx_loc,0,3,1);
                // grab 6DoF KF (Tiw)
                Matrix4d Tiw   = map_keyframes[kf_idx_map]->T_kf_w;
                // projection error
                Tiw = inverse_se3( Tiw );
                Vector3d Pwi   = Tiw.block(0,0,3,3) * Pwj + Tiw.block(0,3,3,1);
                Vector3d Qwi   = Tiw.block(0,0,3,3) * Qwj + Tiw.block(0,3,3,1);
                Vector2d p_prj = cam->projection( Pwi );
                Vector2d q_prj = cam->projection( Qwi );
                Vector3d l_obs = map_lines[lm_idx_map]->obs_list[lm_idx_obs];
                Vector2d l_err;
                l_err(0) = l_obs(0) * p_prj(0) + l_obs(1) * p_prj(1) + l_obs(2);
                l_err(1) = l_obs(0) * q_prj(0) + l_obs(1) * q_prj(1) + l_obs(2);
                double l_err_norm = l_err.norm();
                // start point
                double gx   = Pwi(0);
                double gy   = Pwi(1);
                double gz   = Pwi(2);
                double gz2  = gz*gz;
                gz2         = 1.0 / std::max(SlamConfig::homogTh(),gz2);
                double fx   = cam->getFx();
                double fy   = cam->getFy();
                double lx   = l_err(0);
                double ly   = l_err(1);
                double fxlx = fx*lx;
                double fyly = fy*ly;
                // - jac. wrt. KF pose
                Vector6d Jij_Piw = Vector6d::Zero();
                Jij_Piw << + gz2 * fxlx * gz,
                           + gz2 * fyly * gz,
                           - gz2 * ( fxlx*gx + fyly*gy ),
                           - gz2 * ( fxlx*gx*gy + fyly*gy*gy + fyly*gz*gz ),
                           + gz2 * ( fxlx*gx*gx + fxlx*gz*gz + fyly*gx*gy ),
                           + gz2 * ( fyly*gx*gz - fxlx*gy*gz );
                // - jac. wrt. LM
                Vector3d Jij_Pwj = Vector3d::Zero();
                Jij_Pwj << + gz2 * fxlx * gz,
                           + gz2 * fyly * gz,
                           - gz2 * ( fxlx*gx + fyly*gy );
                Jij_Pwj = Jij_Pwj.transpose() * Tiw.block(0,0,3,3) * l_err(0) / std::max(SlamConfig::homogTh(),l_err_norm);
                // end point
                gx   = Qwi(0);
                gy   = Qwi(1);
                gz   = Qwi(2);
                gz2  = gz*gz;
                gz2  = 1.0 / std::max(SlamConfig::homogTh(),gz2);
                // - jac. wrt. KF pose
                Vector6d Jij_Qiw = Vector6d::Zero();
                Jij_Qiw << + gz2 * fxlx * gz,
                           + gz2 * fyly * gz,
                           - gz2 * ( fxlx*gx + fyly*gy ),
                           - gz2 * ( fxlx*gx*gy + fyly*gy*gy + fyly*gz*gz ),
                           + gz2 * ( fxlx*gx*gx + fxlx*gz*gz + fyly*gx*gy ),
                           + gz2 * ( fyly*gx*gz - fxlx*gy*gz );
                // - jac. wrt. LM
                Vector3d Jij_Qwj = Vector3d::Zero();
                Jij_Qwj << + gz2 * fxlx * gz,
                           + gz2 * fyly * gz,
                           - gz2 * ( fxlx*gx + fyly*gy );
                Jij_Qwj = Jij_Qwj.transpose() * Tiw.block(0,0,3,3) * l_err(1) / std::max(SlamConfig::homogTh(),l_err_norm);
                // estimate Jacobian wrt KF pose
                Vector6d Jij_Tiw = Vector6d::Zero();
                Jij_Tiw = ( Jij_Piw * l_err(0) + Jij_Qiw * l_err(1) ) / std::max(SlamConfig::homogTh(),l_err_norm);
                // estimate Jacobian wrt LM
                Vector6d Jij_Lwj = Vector6d::Zero();
                Jij_Lwj.head(3) = Jij_Pwj;
                Jij_Lwj.tail(3) = Jij_Qwj;
                // if employing robust cost function
                double  w = 1.0;
                double s2 = map_lines[lm_idx_map]->sigma_list[lm_idx_obs];
                //double w = 1.0 / ( 1.0 + l_err_norm * l_err_norm * s2 );
                w = robustWeightCauchy(l_err_norm) ;
                // update hessian, gradient, and error
                VectorXd gi = VectorXd::Zero(N), gj = VectorXd::Zero(N);
                int idx = 6 * kf_idx_loc;
                int jdx = 6*Nkf + 3*Npt + 6*lm_idx_loc;
                if( kf_idx_loc == -1 )
                {
                    gj = Jij_Lwj * l_err_norm * w;
                    err += l_err_norm * l_err_norm * w;
                    Matrix6d Hjj;
                    Hjj = Jij_Lwj * Jij_Lwj.transpose() * w;
                    for(int i = 0; i < 6; i++)
                    {
                        g.coeffRef(jdx+i) += gj(i);
                        for(int j = 0; j < 6; j++)
                            H.coeffRef(i+jdx,j+jdx) += Hjj(i,j);
                    }
                }
                else
                {
                    gi = Jij_Tiw * l_err_norm * w;
                    gj = Jij_Lwj * l_err_norm * w;
                    err += l_err_norm * l_err_norm * w;
                    Matrix6d Hii, Hjj;
                    MatrixXd Hij   = MatrixXd::Zero(3,6);
                    Hii = Jij_Tiw * Jij_Tiw.transpose() * w;
                    Hjj = Jij_Lwj * Jij_Lwj.transpose() * w;
                    Hij = Jij_Lwj * Jij_Tiw.transpose() * w;
                    for(int i = 0; i < 6; i++)
                    {
                        g.coeffRef(i+idx) += gi(i);
                        g.coeffRef(i+jdx) += gj(i);
                        for(int j = 0; j < 6; j++)
                        {
                            H.coeffRef(idx+i,idx+j) += Hii(i,j);
                            H.coeffRef(jdx+i,jdx+j) += Hjj(i,j);
                            H.coeffRef(idx+i,jdx+j) += Hij(i,j);
                            H.coeffRef(jdx+j,idx+i) += Hij(j,i);
                        }
                    }
                }
            }
        }
        err /= (Npt_obs+Nls_obs);
        // if the difference is very small stop
        if( abs(err-err_prev) < numeric_limits<double>::epsilon() || err < numeric_limits<double>::epsilon() )
            break;
        // add lambda to diagonal
        for(int i = 0; i < N; i++)
            H.coeffRef(i,i) += lambda * H.coeffRef(i,i) ;
        // solve iteration
        H.makeCompressed();
        for( int i = 0; i < N; i++)
            gdense(i) = g.coeffRef(i);
        SimplicialLDLT< SparseMatrix<double> > solver1(H);
        DX = solver1.solve( gdense );
        // update lambda
        if( err > err_prev ){
            lambda /= lambda_k;
        }
        else
        {
            lambda *= lambda_k;
            // update KFs
            for( int i = 0; i < Nkf; i++)
            {
                Matrix4d Tprev = expmap_se3( X.block(6*i,0,6,1) );
                Matrix4d Tcurr = Tprev * inverse_se3( expmap_se3( DX.block(6*i,0,6,1) ) );
                X.block(6*i,0,6,1) = logmap_se3( Tcurr );
            }
            // update LMs
            for( int i = 6*Nkf; i < N; i++)
                X(i) += DX(i);
        }
        // if the parameter change is small stop
        if( DX.norm() < numeric_limits<double>::epsilon() )
            break;
        // update previous values
        err_prev = err;
    }

    // Update KFs and LMs
    //---------------------------------------------------------------------------------------------
    // update KFs
    for( int i = 0; i < Nkf; i++)
    {
        Matrix4d Test = expmap_se3( X.block( 6*i,0,6,1 ) );
        map_keyframes[ kf_list[i] ]->T_kf_w = Test;
    }
    // update point LMs
    for( int i = 0; i < Npt; i++)
    {
        map_points[ pt_list[i] ]->point3D(0) = X(6*Nkf+3*i);
        map_points[ pt_list[i] ]->point3D(1) = X(6*Nkf+3*i+1);
        map_points[ pt_list[i] ]->point3D(2) = X(6*Nkf+3*i+2);
    }
    // update line segment LMs
    for( int i = 0; i < Nls; i++)
    {
        map_lines[ ls_list[i] ]->line3D(0) = X(6*Nkf+3*Npt+6*i);
        map_lines[ ls_list[i] ]->line3D(1) = X(6*Nkf+3*Npt+6*i+1);
        map_lines[ ls_list[i] ]->line3D(2) = X(6*Nkf+3*Npt+6*i+2);
        map_lines[ ls_list[i] ]->line3D(3) = X(6*Nkf+3*Npt+6*i+3);
        map_lines[ ls_list[i] ]->line3D(4) = X(6*Nkf+3*Npt+6*i+4);
        map_lines[ ls_list[i] ]->line3D(5) = X(6*Nkf+3*Npt+6*i+5);
    }

}

// -----------------------------------------------------------------------------------------------------------------------------
// Culling functions
// -----------------------------------------------------------------------------------------------------------------------------

void MapHandler::removeBadMapLandmarks()
{

    // point features
    for( vector<MapPoint*>::iterator pt_it = map_points.begin(); pt_it != map_points.end(); pt_it++)
    {
        if( (*pt_it)!=NULL )
        {
            if( (*pt_it)->local == false && max_kf_idx - (*pt_it)->kf_obs_list[0] > 10 )
            {
                if( (*pt_it)->inlier == false || (*pt_it)->obs_list.size() < SlamConfig::minLMObs() )
                {
                    int kf_obs = (*pt_it)->kf_obs_list[0];
                    int lm_idx = (*pt_it)->idx;
                    // remove idx from KeyFrame stereo points
                    for(vector<PointFeature*>::iterator st_pt = map_keyframes[kf_obs]->stereo_frame->stereo_pt.begin();
                        st_pt != map_keyframes[kf_obs]->stereo_frame->stereo_pt.end(); st_pt++ )
                    {
                        if( (*st_pt)->idx == lm_idx )
                        {
                            (*st_pt)->idx = -1;
                            break;
                        }
                    }
                    // remove from map_points_kf_idx
                    int iter = 0;
                    for( auto it = map_points_kf_idx.at(kf_obs).begin(); it != map_points_kf_idx.at(kf_obs).end(); it++, iter++)
                    {
                        if( (*it) == lm_idx )
                        {
                            map_points_kf_idx.at(kf_obs).erase( map_points_kf_idx.at(kf_obs).begin() + iter );
                            break;
                        }
                    }
                    // remove LM
                    delete (*pt_it);
                    (*pt_it) = nullptr;
                }
            }
        }
    }

    // line features
    for( vector<MapLine*>::iterator ls_it = map_lines.begin(); ls_it != map_lines.end(); ls_it++)
    {
        if((*ls_it) != NULL)
        {
            if( (*ls_it)->local == false && max_kf_idx-(*ls_it)->kf_obs_list[0] > 10 )
            {
                if( (*ls_it)->inlier == false || (*ls_it)->obs_list.size() < SlamConfig::minLMObs() )
                {
                    int kf_obs = (*ls_it)->kf_obs_list[0];
                    int lm_idx = (*ls_it)->idx;
                    // remove idx from KeyFrame stereo points
                    for(vector<LineFeature*>::iterator st_ls = map_keyframes[kf_obs]->stereo_frame->stereo_ls.begin();
                        st_ls != map_keyframes[kf_obs]->stereo_frame->stereo_ls.end(); st_ls++ )
                    {
                        if( (*st_ls)->idx == lm_idx )
                        {
                            (*st_ls)->idx = -1;
                            break;
                        }
                    }
                    // remove from map_points_kf_idx
                    int iter = 0;
                    for( auto it = map_lines_kf_idx.at(kf_obs).begin(); it != map_lines_kf_idx.at(kf_obs).end(); it++, iter++)
                    {
                        if( (*it) == lm_idx )
                        {
                            map_lines_kf_idx.at(kf_obs).erase( map_lines_kf_idx.at(kf_obs).begin() + iter );
                            break;
                        }
                    }

                    // remove LM
                    delete (*ls_it);
                    (*ls_it) = nullptr;
                }
            }
        }
    }
}

void MapHandler::removeRedundantKFs()
{

    // select which KFs are going to remove
    vector<int> kf_idxs;
    for( vector<KeyFrame*>::iterator kf_it = map_keyframes.begin(); kf_it != map_keyframes.end(); kf_it++)
    {
        if((*kf_it) != NULL)
        {
            int kf_idx = (*kf_it)->kf_idx;
            if( !(*kf_it)->local && kf_idx > 1 && kf_idx < max_kf_idx )
            {
                // estimate number of landmarks observed by this KF
                int n_feats = 0;
                for( vector<PointFeature*>::iterator pt_it = (*kf_it)->stereo_frame->stereo_pt.begin();
                     pt_it != (*kf_it)->stereo_frame->stereo_pt.end(); pt_it++ )
                {
                    if( (*pt_it)->idx != -1 )
                        n_feats++;
                }
                for( vector<LineFeature*>::iterator ls_it = (*kf_it)->stereo_frame->stereo_ls.begin();
                     ls_it != (*kf_it)->stereo_frame->stereo_ls.end(); ls_it++)
                {
                    if( (*ls_it)->idx != -1 )
                        n_feats++;
                }
                int max_n_feats = int( SlamConfig::maxCommonFtsKF() * double(n_feats) );
                // check if the KF is redundant
                int n_graph = full_graph.size();
                for(int i = 0; i < n_graph-1, i != kf_idx; i++)
                {
                    if( map_keyframes[i] != NULL )
                    {
                        if( full_graph[kf_idx][i] > n_feats )
                            kf_idxs.push_back(kf_idx);
                    }
                }
            }
        }
    }

    // eliminate KFs, LMs observed only by this KFs, and all observations from this KF
    for( int i = 0; i < kf_idxs.size(); i++)
    {
        int kf_idx = kf_idxs[i];
        if( map_keyframes[kf_idx] != NULL )
        {
            // delete observation from map_points_kf_idx
            if( !(map_points_kf_idx.find(kf_idx) == map_points_kf_idx.end()) )
            {
                for( auto it = map_points_kf_idx.at(kf_idx).begin(); it != map_points_kf_idx.at(kf_idx).end(); it++)
                {
                    if( map_points[(*it)] != NULL &&  map_points[(*it)]->kf_obs_list.size() <= 1 )
                    {
                        bool found = false;
                        for( int k = 1; k < map_points[(*it)]->kf_obs_list.size(); k++)
                        {
                            map_points_kf_idx.find( map_points[(*it)]->kf_obs_list[k] );
                            int new_kf_base = map_points[(*it)]->kf_obs_list[k];    // if the second dont exist...
                            if( !(map_points_kf_idx.find(new_kf_base)==map_points_kf_idx.end()) )
                            {
                                map_points_kf_idx.at(new_kf_base).push_back( (*it) );
                                found = true;
                            }
                        }
                        if( !found )
                            map_points[(*it)]->inlier = false;
                    }
                }
                map_points_kf_idx.erase(kf_idx);
            }

            // delete observation from map_lines_kf_idx
            for( auto it = map_lines_kf_idx.at(kf_idx).begin(); it != map_lines_kf_idx.at(kf_idx).end(); it++)
            {
                if( map_points[(*it)]!= NULL )
                {
                    int new_kf_base = map_lines[(*it)]->kf_obs_list[1];
                    map_lines_kf_idx.at(new_kf_base).push_back( (*it) );
                }
            }
            map_lines_kf_idx.erase(kf_idx);

            // iterate over point features
            for( vector<PointFeature*>::iterator pt_it = map_keyframes[kf_idx]->stereo_frame->stereo_pt.begin();
                 pt_it != map_keyframes[kf_idx]->stereo_frame->stereo_pt.end(); pt_it++)
            {
                int pt_idx = (*pt_it)->idx;
                if( pt_idx != -1 )
                {
                    if( map_points[pt_idx] != NULL )
                    {
                        for( int j = 0; j < map_points[pt_idx]->obs_list.size(); j++)
                        {
                            if( map_points[pt_idx]->kf_obs_list[j] == kf_idx )
                            {
                                // delete observation
                                map_points[pt_idx]->desc_list.erase( map_points[pt_idx]->desc_list.begin() + j );
                                map_points[pt_idx]->obs_list.erase( map_points[pt_idx]->obs_list.begin() + j );
                                map_points[pt_idx]->dir_list.erase( map_points[pt_idx]->dir_list.begin() + j );
                                map_points[pt_idx]->kf_obs_list.erase( map_points[pt_idx]->kf_obs_list.begin() + j );
                                // update main descriptor and direction
                                map_points[pt_idx]->updateAverageDescDir();
                            }
                        }
                    }
                }
            }

            // iterate over line segment features
            for( vector<LineFeature*>::iterator ls_it = map_keyframes[kf_idx]->stereo_frame->stereo_ls.begin();
                 ls_it != map_keyframes[kf_idx]->stereo_frame->stereo_ls.end(); ls_it++)
            {
                int ls_idx = (*ls_it)->idx;
                if( ls_idx != -1 )
                {
                    if( map_lines[ls_idx] != NULL )
                    {
                        for( int j = 0; j < map_lines[ls_idx]->obs_list.size(); j++)
                        {
                            if( map_lines[ls_idx]->kf_obs_list[j] == kf_idx )
                            {
                                // delete observation
                                map_lines[ls_idx]->desc_list.erase( map_lines[ls_idx]->desc_list.begin() + j );
                                map_lines[ls_idx]->obs_list.erase( map_lines[ls_idx]->obs_list.begin() + j );
                                map_lines[ls_idx]->dir_list.erase( map_lines[ls_idx]->dir_list.begin() + j );
                                map_lines[ls_idx]->kf_obs_list.erase( map_lines[ls_idx]->kf_obs_list.begin() + j );
                                map_lines[ls_idx]->pts_list.erase( map_lines[ls_idx]->pts_list.begin() + j );
                                // update main descriptor and direction
                                map_lines[ls_idx]->updateAverageDescDir();
                            }
                        }
                    }
                }
            }

            // update full graph
            for( int k = 0; k < full_graph.size()-1; k++ )
            {
                full_graph[kf_idx][k] = 0;
                full_graph[k][kf_idx] = 0;
            }

            // erase KF
            delete map_keyframes[kf_idx];
            map_keyframes[kf_idx] = nullptr;
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------
// Loop Closure functions
// -----------------------------------------------------------------------------------------------------------------------------

void MapHandler::loopClosure()
{

    Timer timer;

    // look for loop closure candidates
    int kf_prev_idx, kf_curr_idx;
    kf_curr_idx = max_kf_idx;
    timer.start();
    bool is_lc_candidate = lookForLoopCandidates(kf_curr_idx,kf_prev_idx);
    time(4) = timer.stop(); // ms

    // compute relative transformation if it is LC candidate
    if( is_lc_candidate )
    {
        vector<Vector4i> lc_pt_idx, lc_ls_idx;
        vector<PointFeature*> lc_points;
        vector<LineFeature*>  lc_lines;
        Vector6d pose_inc;
        timer.start();
        bool isLC = isLoopClosure( map_keyframes[kf_prev_idx], map_keyframes[kf_curr_idx], pose_inc, lc_pt_idx, lc_ls_idx, lc_points, lc_lines );
        time(5) = timer.stop(); //ms
        // if it is loop closure, add information and update status
        if( isLC )
        {
            lc_pt_idxs.push_back( lc_pt_idx );
            lc_ls_idxs.push_back( lc_ls_idx );
            lc_poses.push_back( pose_inc );
            lc_pose_list.push_back( pose_inc );
            Vector3i lc_idx;
            lc_idx(0) = kf_prev_idx;
            lc_idx(1) = kf_curr_idx;
            lc_idx(2) = 1;
            lc_idxs.push_back( lc_idx );
            lc_idx_list.push_back(lc_idx);
            if( lc_state == LC_IDLE )
                lc_state = LC_ACTIVE;
        }
        else
        {
            if( lc_state == LC_ACTIVE )
                lc_state = LC_READY;
        }

        for (PointFeature* pt : lc_points)
            delete pt;
        for (LineFeature* ls : lc_lines)
            delete ls;
    }
    else
    {
        if( lc_state == LC_ACTIVE )
            lc_state = LC_READY;
    }

    // LC computation
    if( lc_state == LC_READY )  // add condition indicating that the LC has finished
    {
        timer.start();
        loopClosureOptimizationCovGraphG2O();
        time(6) = timer.stop(); //ms
        lc_state = LC_IDLE;
    }
}

void MapHandler::insertKFBowVectorP( KeyFrame* kf  )
{

    // transform Mat to vector<Mat>
    vector<Mat> curr_desc;
    curr_desc.reserve( kf->stereo_frame->pdesc_l.rows );
    for ( int i = 0; i < kf->stereo_frame->pdesc_l.rows; i++ )
        curr_desc.push_back( kf->stereo_frame->pdesc_l.row(i) );
    // transform to DBoW2::BowVector
    dbow_voc_p.transform( curr_desc, kf->descDBoW_P );

    // fill the confusion matrix for the new KF
    int idx = kf->kf_idx;
    for( int i = 0; i < idx; i++)
    {
        if( map_keyframes[i] != NULL )
        {
            double score = dbow_voc_p.score( kf->descDBoW_P, map_keyframes[i]->descDBoW_P );
            conf_matrix[idx][i] = score;
            conf_matrix[i][idx] = score;
        }
    }
    conf_matrix[idx][idx] = dbow_voc_p.score( kf->descDBoW_P, kf->descDBoW_P );
}

void MapHandler::insertKFBowVectorL( KeyFrame* kf  )
{

    // transform Mat to vector<Mat>
    vector<Mat> curr_desc;
    curr_desc.reserve( kf->stereo_frame->ldesc_l.rows );
    for ( int i = 0; i < kf->stereo_frame->ldesc_l.rows; i++ )
        curr_desc.push_back( kf->stereo_frame->ldesc_l.row(i) );
    // transform to DBoW2::BowVector
    dbow_voc_l.transform( curr_desc, kf->descDBoW_L );

    // fill the confusion matrix for the new KF
    int idx = kf->kf_idx;
    for( int i = 0; i < idx; i++)
    {
        if( map_keyframes[i] != NULL )
        {
            double score = dbow_voc_l.score( kf->descDBoW_L, map_keyframes[i]->descDBoW_L );
            conf_matrix[idx][i] = score;
            conf_matrix[i][idx] = score;
        }
    }
    conf_matrix[idx][idx] = dbow_voc_l.score( kf->descDBoW_L, kf->descDBoW_L );
}

void MapHandler::insertKFBowVectorPL( KeyFrame* kf  )
{

    // Point Features
    // --------------------------------------------------------------
    // transform Mat to vector<Mat>
    vector<Mat> curr_desc;
    curr_desc.reserve( kf->stereo_frame->pdesc_l.rows );
    for ( int i = 0; i < kf->stereo_frame->pdesc_l.rows; i++ )
        curr_desc.push_back( kf->stereo_frame->pdesc_l.row(i) );
    // transform to DBoW2::BowVector
    dbow_voc_p.transform( curr_desc, kf->descDBoW_P );

    // estimate dispersion of point features
    vector<double> pt_x, pt_y;
    for (PointFeature* pt : kf->stereo_frame->stereo_pt)
    {
        pt_x.push_back( pt->pl(0) );
        pt_y.push_back( pt->pl(1) );
    }
    double std_pt = vector_stdv( pt_x ) + vector_stdv( pt_y );
    int    n_pt   = pt_x.size();

    // Line Segment Features
    // --------------------------------------------------------------
    // transform Mat to vector<Mat>
    curr_desc.clear();
    curr_desc.reserve( kf->stereo_frame->ldesc_l.rows );
    for ( int i = 0; i < kf->stereo_frame->ldesc_l.rows; i++ )
        curr_desc.push_back( kf->stereo_frame->ldesc_l.row(i) );
    // transform to DBoW2::BowVector
    dbow_voc_l.transform( curr_desc, kf->descDBoW_L );

    // estimate dispersion of point features
    vector<double> ls_x, ls_y;
    for (LineFeature* ls : kf->stereo_frame->stereo_ls)
    {
        Vector2d mp;
        mp << (ls->spl + ls->epl)*0.5;
        ls_x.push_back( mp(0) );
        ls_y.push_back( mp(1) );
    }
    double std_ls = vector_stdv( ls_x ) + vector_stdv( ls_y );
    double std_pl = std_ls + std_pt;
    int    n_ls  = ls_x.size();
    int    n_pl  = n_pt + n_ls;

    // Combined Approach
    // --------------------------------------------------------------
    // fill the confusion matrix for the new KF
    double score, score_p, score_l;
    int idx = kf->kf_idx;
    for( int i = 0; i < idx; i++)
    {
        if( map_keyframes[i] != NULL )
        {
            score_p = dbow_voc_p.score( kf->descDBoW_P, map_keyframes[i]->descDBoW_P );
            score_l = dbow_voc_l.score( kf->descDBoW_L, map_keyframes[i]->descDBoW_L );
            score = 0.0;
            score += ( score_p * n_pt   + score_l * n_ls   ) / n_pl;    // strategy#1
            score += ( score_p * std_pt + score_l * std_ls ) / std_pl;  // strategy#2
            conf_matrix[idx][i] = score;
            conf_matrix[i][idx] = score;
        }
    }
    score_p = dbow_voc_p.score( kf->descDBoW_P, kf->descDBoW_P );
    score_l = dbow_voc_l.score( kf->descDBoW_L, kf->descDBoW_L);
    score = 0.0;
    score += ( score_p * n_pt   + score_l * n_ls   ) / n_pl;    // strategy#1
    score += ( score_p * std_pt + score_l * std_ls ) / std_pl;  // strategy#2
    conf_matrix[idx][idx] = score;
}

bool MapHandler::lookForLoopCandidates( int kf_curr_idx, int &kf_prev_idx )
{
    bool is_lc_candidate = false;
    kf_prev_idx = -1;

    // find the best matches
    vector<Vector2d> max_confmat;
    for(int i = 0; i < kf_curr_idx - SlamConfig::lcKFDist(); i++)
    {
        if( map_keyframes[i] != NULL )
        {
            Vector2d aux;
            aux(0) = i;
            aux(1) = conf_matrix[i][kf_curr_idx];
            max_confmat.push_back( aux );
        }
    }

    // if there are enough matches..
    if( max_confmat.size() > SlamConfig::lcKFMaxDist() )
    {
        sort( max_confmat.begin(), max_confmat.end(), sort_confmat_by_score() );

        // find the minimum score in the covisibility graph
        double lc_min_score = 1.0;
        for( int i = 0; i < kf_curr_idx; i++ )
        {
            if( full_graph[i][kf_curr_idx] >= SlamConfig::minLMCovGraph() || kf_curr_idx - i <= SlamConfig::minKFLocalMap()+3 )
            {
                double score_i = conf_matrix[i][kf_curr_idx];
                if( score_i < lc_min_score && score_i > 0.001 )
                    lc_min_score = score_i;
            }
        }

        // the best match must has an score above lc_dbow_score_max
        int idx_max = int( max_confmat[0](0) );
        int Nkf_closest = 0;
        if( max_confmat[0](1) >= lc_min_score )
        {
            // there must be at least lc_nkf_closest KFs conected to the LC candidate with a score above lc_dbow_score_min
            for( int i = 1; i < max_confmat.size(); i++ )
            {
                int idx = int( max_confmat[i](0) );
                // frame closest && connected by the cov_graph && score > lc_dbow_score_min
                if( abs(idx-idx_max) <= SlamConfig::lcKFMaxDist() &&
                    max_confmat[i](1) >= lc_min_score * 0.8 )
                    Nkf_closest++;
            }

            // update in case of being loop closure candidate
            if( Nkf_closest >= SlamConfig::lcNKFClosest() )
            {
                is_lc_candidate = true;
                kf_prev_idx     = idx_max;
            }
        }
    }

    return is_lc_candidate;
}

bool MapHandler::isLoopClosure( const KeyFrame* kf0, const KeyFrame* kf1, Vector6d &pose_inc,
                                vector<Vector4i> &lc_pt_idx, vector<Vector4i> &lc_ls_idx,
                                vector<PointFeature*> &lc_points, vector<LineFeature*>  &lc_lines)
{

    // grab frame number
    int kf0_idx = kf0->kf_idx;
    int kf1_idx = kf1->kf_idx;

    // number of stereo matches
    int n_pt_0 = kf0->stereo_frame->stereo_pt.size();
    int n_pt_1 = kf1->stereo_frame->stereo_pt.size();
    int n_ls_0 = kf0->stereo_frame->stereo_ls.size();
    int n_ls_1 = kf1->stereo_frame->stereo_ls.size();

    // initial pose increment
    Matrix4d DT = Matrix4d::Identity();

    // structures to optimize camera drift
    lc_pt_idx.clear();
    lc_ls_idx.clear();
    lc_points.clear();
    lc_lines.clear();

    // find matches between both KFs
    // ---------------------------------------------------
    // points f2f tracking
    int common_pt = 0;
    if( SlamConfig::hasPoints() && !kf1->stereo_frame->stereo_pt.empty() && !kf0->stereo_frame->stereo_pt.empty() )
    {
        vector<int> matches_12;
        common_pt = match(kf0->stereo_frame->pdesc_l, kf1->stereo_frame->pdesc_l, SlamConfig::minRatio12P(), matches_12);

        for (int i1 = 0; i1 < matches_12.size(); ++i1) {
            const int i2 = matches_12[i1];
            if (i2 < 0) continue;

            // save data for optimization
            Vector3d P       = kf0->stereo_frame->stereo_pt[i1]->P;
            Vector2d pl_obs  = kf1->stereo_frame->stereo_pt[i2]->pl;
            PointFeature* pt = new PointFeature( P, pl_obs );
            lc_points.push_back(pt);
            // save indices for fusing LMs
            Vector4i idx;
            idx(0) = kf0->stereo_frame->stereo_pt[i1]->idx;
            idx(1) = i1;
            idx(2) = kf1->stereo_frame->stereo_pt[i2]->idx;
            idx(3) = i2;
            lc_pt_idx.push_back(idx);
        }
    }

    // line segments f2f tracking
    int common_ls = 0;
    if( SlamConfig::hasLines() && !kf1->stereo_frame->stereo_ls.empty() && !kf0->stereo_frame->stereo_ls.empty() )
    {
        vector<int> matches_12;
        common_ls = match(kf0->stereo_frame->ldesc_l, kf1->stereo_frame->ldesc_l, SlamConfig::minRatio12L(), matches_12);

        for (int i1 = 0; i1 < matches_12.size(); ++i1) {
            const int i2 = matches_12[i1];
            if (i2 < 0) continue;

            // save data for optimization
            Vector3d sP     = kf0->stereo_frame->stereo_ls[i1]->sP;
            Vector3d eP     = kf0->stereo_frame->stereo_ls[i1]->eP;
            Vector3d le_obs = kf1->stereo_frame->stereo_ls[i2]->le;
            LineFeature* ls = new LineFeature( sP, eP, le_obs, kf1->stereo_frame->stereo_ls[i2]->spl, kf1->stereo_frame->stereo_ls[i2]->epl );
            lc_lines.push_back(ls);
            // save indices for fusing LMs
            Vector4i idx;
            idx(0) = kf0->stereo_frame->stereo_ls[i1]->idx;
            idx(1) = i1;
            idx(2) = kf1->stereo_frame->stereo_ls[i2]->idx;
            idx(3) = i2;
            lc_ls_idx.push_back(idx);
        }
    }

    // estimate relative pose between both KFs
    // ---------------------------------------------------
    double inl_ratio_pt = max( 100.0 * common_pt / n_pt_0, 100.0 * common_pt / n_pt_1 );
    double inl_ratio_ls = max( 100.0 * common_ls / n_ls_0, 100.0 * common_ls / n_ls_1 );
    bool inl_ratio_condition = false;

    if( SlamConfig::hasPoints() && SlamConfig::hasLines() )
    {
        if( inl_ratio_pt > SlamConfig::lcInlierRatio() && inl_ratio_ls > SlamConfig::lcInlierRatio() )
            inl_ratio_condition = true;
    }
    else if( SlamConfig::hasPoints() )
    {
        if( inl_ratio_pt > SlamConfig::lcInlierRatio() )
            inl_ratio_condition = true;
    }
    else if( SlamConfig::hasLines() )
    {
        if( inl_ratio_ls > SlamConfig::lcInlierRatio() )
            inl_ratio_condition = true;
    }


    //inl_ratio_condition = true;
    if( inl_ratio_condition )
        return computeRelativePoseRobustGN(lc_points,lc_lines,lc_pt_idx,lc_ls_idx,pose_inc);
    else
        return false;

}

bool MapHandler::computeRelativePoseGN( vector<PointFeature*> &lc_points, vector<LineFeature*> &lc_lines,
                                        vector<Vector4i>      &lc_pt_idx, vector<Vector4i>     &lc_ls_idx,
                                        Vector6d &pose_inc) const
{

    // create GN variables
    Vector6d x_inc = Vector6d::Zero(), x_prev = Vector6d::Zero();
    Matrix4d T_inc = Matrix4d::Identity(), T_prev = Matrix4d::Zero();
    Matrix6d H_l, H_p, H;
    Vector6d g_l, g_p, g;
    H = Matrix6d::Zero();
    g = Vector6d::Zero();
    H_l = H; H_p = H;
    g_l = g; g_p = g;
    double   e_l = 0.0, e_p = 0.0, e = 0.0, S_l, S_p;

    // create GN parameters
    double err, err_prev = 999999999.9;
    int    max_iters_first = SlamConfig::maxIters(), max_iters = SlamConfig::maxItersRef();

    // GN iterations
    //---------------------------------------------------------------------------------------------
    for( int iters = 0; iters < max_iters_first; iters++)
    {
        H = Matrix6d::Zero(); H_l = H; H_p = H;
        g = Vector6d::Zero(); g_l = g; g_p = g;
        e = 0.0; e_p = 0.0; e_l = 0.0;
        // point observations
        int N_p = 0;
        for( vector<PointFeature*>::iterator pt_it = lc_points.begin(); pt_it != lc_points.end(); pt_it++ )
        {
            if( (*pt_it)->inlier )
            {
                Vector3d P_ = T_inc.block(0,0,3,3) * (*pt_it)->P + T_inc.col(3).head(3);
                Vector2d pl_proj = cam->projection( P_ );
                // projection error
                Vector2d err_i    = pl_proj - (*pt_it)->pl_obs;
                double err_i_norm = err_i.norm();
                double gx   = P_(0);
                double gy   = P_(1);
                double gz   = P_(2);
                double gz2  = gz*gz;
                double fgz2 = cam->getFx() / std::max(SlamConfig::homogTh(),gz2);
                double dx   = err_i(0);
                double dy   = err_i(1);
                // jacobian
                Vector6d J_aux;
                J_aux << + fgz2 * dx * gz,
                         + fgz2 * dy * gz,
                         - fgz2 * ( gx*dx + gy*dy ),
                         - fgz2 * ( gx*gy*dx + gy*gy*dy + gz*gz*dy ),
                         + fgz2 * ( gx*gx*dx + gz*gz*dx + gx*gy*dy ),
                         + fgz2 * ( gx*gz*dy - gy*gz*dx );
                J_aux = J_aux / std::max(SlamConfig::homogTh(),err_i_norm);
                // if employing robust cost function
                double  w = 1.0;
                double s2 = (*pt_it)->sigma2;
                //double w = 1.0 / ( 1.0 + err_i_norm * err_i_norm * s2 );
                w = robustWeightCauchy(err_i_norm) ;

                // update hessian, gradient, and error
                H_p += J_aux * J_aux.transpose() * w;
                g_p += J_aux * err_i_norm * w;
                e_p += err_i_norm * err_i_norm * w;
                N_p++;
            }
        }
        // line segment features
        int N_l = 0;
        for( vector<LineFeature*>::iterator ls_it = lc_lines.begin(); ls_it != lc_lines.end(); ls_it++ )
        {
            if( (*ls_it)->inlier )
            {
                Vector3d sP_ = T_inc.block(0,0,3,3) * (*ls_it)->sP + T_inc.col(3).head(3);
                Vector2d spl_proj = cam->projection( sP_ );
                Vector3d eP_ = T_inc.block(0,0,3,3) * (*ls_it)->eP + T_inc.col(3).head(3);
                Vector2d epl_proj = cam->projection( eP_ );
                Vector3d l_obs = (*ls_it)->le_obs;
                // projection error
                Vector2d err_i;
                err_i(0) = l_obs(0) * spl_proj(0) + l_obs(1) * spl_proj(1) + l_obs(2);
                err_i(1) = l_obs(0) * epl_proj(0) + l_obs(1) * epl_proj(1) + l_obs(2);
                double err_i_norm = err_i.norm();
                // start point
                double gx   = sP_(0);
                double gy   = sP_(1);
                double gz   = sP_(2);
                double gz2  = gz*gz;
                double fgz2 = cam->getFx() / std::max(SlamConfig::homogTh(),gz2);
                double ds   = err_i(0);
                double de   = err_i(1);
                double lx   = l_obs(0);
                double ly   = l_obs(1);
                Vector6d Js_aux;
                Js_aux << + fgz2 * lx * gz,
                          + fgz2 * ly * gz,
                          - fgz2 * ( gx*lx + gy*ly ),
                          - fgz2 * ( gx*gy*lx + gy*gy*ly + gz*gz*ly ),
                          + fgz2 * ( gx*gx*lx + gz*gz*lx + gx*gy*ly ),
                          + fgz2 * ( gx*gz*ly - gy*gz*lx );
                // end point
                gx   = eP_(0);
                gy   = eP_(1);
                gz   = eP_(2);
                gz2  = gz*gz;
                fgz2 = cam->getFx() / std::max(SlamConfig::homogTh(),gz2);
                Vector6d Je_aux, J_aux;
                Je_aux << + fgz2 * lx * gz,
                          + fgz2 * ly * gz,
                          - fgz2 * ( gx*lx + gy*ly ),
                          - fgz2 * ( gx*gy*lx + gy*gy*ly + gz*gz*ly ),
                          + fgz2 * ( gx*gx*lx + gz*gz*lx + gx*gy*ly ),
                          + fgz2 * ( gx*gz*ly - gy*gz*lx );
                // jacobian
                J_aux = ( Js_aux * ds + Je_aux * de ) / std::max(SlamConfig::homogTh(),err_i_norm);
                // if employing robust cost function
                double  w = 1.0;
                w = robustWeightCauchy(err_i_norm) ;
                // update hessian, gradient, and error
                H_l += J_aux * J_aux.transpose() * w;
                g_l += J_aux * err_i_norm * w;
                e_l += err_i_norm * err_i_norm * w;
                N_l++;
            }
        }
        // sum H, g and err from both points and lines
        H = H_p + H_l;
        g = g_p + g_l;
        e = e_p + e_l;
        // normalize error
        e /= (N_l+N_p);
        // if the difference is very small stop
        if( abs(e-err_prev) < numeric_limits<double>::epsilon() || e < numeric_limits<double>::epsilon() )
            break;
        // solve
        ColPivHouseholderQR<MatrixXd> solver(H);
        x_inc = solver.solve( g );
        T_inc = T_inc * inverse_se3( expmap_se3(x_inc) );
        // if the parameter change is small stop
        if( x_inc.norm() < numeric_limits<double>::epsilon() )
            break;
        // update error
        err_prev = e;
    }
    x_inc = logmap_se3(T_inc);

    // Remove outliers
    //---------------------------------------------------------------------------------------------
    // point observations
    for( vector<PointFeature*>::iterator pt_it = lc_points.begin(); pt_it != lc_points.end(); pt_it++ )
    {
        if( (*pt_it)->inlier )
        {
            Vector3d P_ = T_inc.block(0,0,3,3) * (*pt_it)->P + T_inc.col(3).head(3);
            Vector2d pl_proj = cam->projection( P_ );
            // projection error
            Vector2d err_i    = pl_proj - (*pt_it)->pl_obs;
            double s2 = (*pt_it)->sigma2;
            if( err_i.norm() > sqrt(7.815) )
                (*pt_it)->inlier = false;
        }
    }
    // line segments observations
    for( vector<LineFeature*>::iterator ls_it = lc_lines.begin(); ls_it != lc_lines.end(); ls_it++ )
    {
        if((*ls_it)->inlier)
        {
            Vector3d sP_ = T_inc.block(0,0,3,3) * (*ls_it)->sP + T_inc.col(3).head(3);
            Vector2d spl_proj = cam->projection( sP_ );
            Vector3d eP_ = T_inc.block(0,0,3,3) * (*ls_it)->eP + T_inc.col(3).head(3);
            Vector2d epl_proj = cam->projection( eP_ );
            Vector3d l_obs = (*ls_it)->le_obs;
            // projection error
            Vector2d err_i;
            err_i(0) = l_obs(0) * spl_proj(0) + l_obs(1) * spl_proj(1) + l_obs(2);
            err_i(1) = l_obs(0) * epl_proj(0) + l_obs(1) * epl_proj(1) + l_obs(2);
            double s2 = sqrt((*ls_it)->sigma2);
            if( err_i.norm() > sqrt(7.815) )
                (*ls_it)->inlier = false;
        }
    }

    // Check whether it is Loop Closure or not
    //---------------------------------------------------------------------------------------------
    // Residue value
    bool lc_res = ( e < SlamConfig::lcRes() );

    // Uncertainty value
    Matrix6d DT_cov;
    DT_cov = H.inverse();
    SelfAdjointEigenSolver<Matrix6d> eigensolver(DT_cov);
    Vector6d DT_cov_eig = eigensolver.eigenvalues();
    bool lc_unc = ( DT_cov_eig(5) < SlamConfig::lcUnc() );

    // Ratio of outliers
    int N = lc_points.size() + lc_lines.size(), N_inl = 0;
    for( vector<PointFeature*>::iterator pt_it = lc_points.begin(); pt_it != lc_points.end(); pt_it++ )
    {
        if( (*pt_it)->inlier )
            N_inl++;
    }
    for( vector<LineFeature*>::iterator  ls_it = lc_lines.begin();  ls_it != lc_lines.end();  ls_it++ )
    {
        if( (*ls_it)->inlier )
            N_inl++;
    }
    double ratio_inliers = double(N_inl) / double(N);
    bool lc_inl = ( ratio_inliers > SlamConfig::lcInl() );
    //lc_inl = true;

    // Geometry condition
    double t = x_inc.head(3).norm();
    double r = x_inc.tail(3).norm() * 180.f / CV_PI;
    bool lc_trs = ( t < SlamConfig::lcTrs() );
    bool lc_rot = ( r < SlamConfig::lcRot() );

    // Decision
    if( lc_res && lc_unc && lc_inl && lc_trs && lc_rot )
    {
        // erase outliers from lc_pt_idx & lc_ls_idx
        int iter = 0;
        vector<Vector4i> lc_pt_idx_, lc_ls_idx_;
        vector<PointFeature*> lc_points_;
        vector<LineFeature*>  lc_lines_;
        for( vector<PointFeature*>::iterator pt_it = lc_points.begin(); pt_it != lc_points.end(); pt_it++, iter++ )
        {
            if( (*pt_it)!= NULL )
            {
                if( (*pt_it)->inlier )
                {
                    lc_points_.push_back( (*pt_it) );
                    lc_pt_idx_.push_back( lc_pt_idx[iter] );
                }
            }
        }
        iter = 0;
        for( vector<LineFeature*>::iterator ls_it = lc_lines.begin(); ls_it != lc_lines.end(); ls_it++, iter++ )
        {
            if( (*ls_it)!= NULL )
            {
                if( (*ls_it)->inlier )
                {
                    lc_lines_.push_back( (*ls_it) );
                    lc_ls_idx_.push_back( lc_ls_idx[iter] );
                }
            }
        }
        lc_pt_idx.clear();
        lc_pt_idx = lc_pt_idx_;
        lc_ls_idx.clear();
        lc_ls_idx = lc_ls_idx_;
        lc_points.clear();
        lc_points = lc_points_;
        lc_lines.clear();
        lc_points = lc_points_;
        // assign pose increment
        pose_inc = logmap_se3( inverse_se3( T_inc ) );
        return true;
    }
    else
        return false;

}

bool MapHandler::computeRelativePoseRobustGN( vector<PointFeature*> &lc_points, vector<LineFeature*> &lc_lines,
                                        vector<Vector4i>      &lc_pt_idx, vector<Vector4i>     &lc_ls_idx,
                                        Vector6d &pose_inc) const
{

    // create GN variables
    Vector6d x_inc = Vector6d::Zero(), x_prev = Vector6d::Zero();
    Matrix4d T_inc = Matrix4d::Identity(), T_prev = Matrix4d::Zero();
    Matrix6d H_l, H_p, H;
    Vector6d g_l, g_p, g;
    H = Matrix6d::Zero();
    g = Vector6d::Zero();
    H_l = H; H_p = H;
    g_l = g; g_p = g;
    double   e_l = 0.0, e_p = 0.0, e = 0.0, S_l, S_p;

    // create GN parameters
    double err, err_prev = 999999999.9;
    int    max_iters_first = SlamConfig::maxIters(), max_iters = SlamConfig::maxItersRef();

    // GN iterations
    //---------------------------------------------------------------------------------------------
    for( int iters = 0; iters < max_iters_first; iters++)
    {
        H = Matrix6d::Zero(); H_l = H; H_p = H;
        g = Vector6d::Zero(); g_l = g; g_p = g;
        e = 0.0; e_p = 0.0; e_l = 0.0;
        // point observations
        int N_p = 0;
        for( vector<PointFeature*>::iterator pt_it = lc_points.begin(); pt_it != lc_points.end(); pt_it++ )
        {
            if( (*pt_it)->inlier )
            {
                Vector3d P_ = T_inc.block(0,0,3,3) * (*pt_it)->P + T_inc.col(3).head(3);
                Vector2d pl_proj = cam->projection( P_ );
                // projection error
                Vector2d err_i    = pl_proj - (*pt_it)->pl_obs;
                double err_i_norm = err_i.norm();
                double gx   = P_(0);
                double gy   = P_(1);
                double gz   = P_(2);
                double gz2  = gz*gz;
                double fgz2 = cam->getFx() / std::max(SlamConfig::homogTh(),gz2);
                double dx   = err_i(0);
                double dy   = err_i(1);
                // jacobian
                Vector6d J_aux;
                J_aux << + fgz2 * dx * gz,
                         + fgz2 * dy * gz,
                         - fgz2 * ( gx*dx + gy*dy ),
                         - fgz2 * ( gx*gy*dx + gy*gy*dy + gz*gz*dy ),
                         + fgz2 * ( gx*gx*dx + gz*gz*dx + gx*gy*dy ),
                         + fgz2 * ( gx*gz*dy - gy*gz*dx );
                J_aux = J_aux / std::max(SlamConfig::homogTh(),err_i_norm);
                // if employing robust cost function
                double  w = 1.0;
                double s2 = (*pt_it)->sigma2;
                w = robustWeightCauchy(err_i_norm) ;
                // update hessian, gradient, and error
                H_p += J_aux * J_aux.transpose() * w;
                g_p += J_aux * err_i_norm * w;
                e_p += err_i_norm * err_i_norm * w;
                N_p++;
            }
        }
        // line segment features
        int N_l = 0;
        for( vector<LineFeature*>::iterator ls_it = lc_lines.begin(); ls_it != lc_lines.end(); ls_it++ )
        {
            if( (*ls_it)->inlier )
            {
                Vector3d sP_ = T_inc.block(0,0,3,3) * (*ls_it)->sP + T_inc.col(3).head(3);
                Vector2d spl_proj = cam->projection( sP_ );
                Vector3d eP_ = T_inc.block(0,0,3,3) * (*ls_it)->eP + T_inc.col(3).head(3);
                Vector2d epl_proj = cam->projection( eP_ );
                Vector3d l_obs = (*ls_it)->le_obs;
                // projection error
                Vector2d err_i;
                err_i(0) = l_obs(0) * spl_proj(0) + l_obs(1) * spl_proj(1) + l_obs(2);
                err_i(1) = l_obs(0) * epl_proj(0) + l_obs(1) * epl_proj(1) + l_obs(2);
                double err_i_norm = err_i.norm();
                // start point
                double gx   = sP_(0);
                double gy   = sP_(1);
                double gz   = sP_(2);
                double gz2  = gz*gz;
                double fgz2 = cam->getFx() / std::max(SlamConfig::homogTh(),gz2);
                double ds   = err_i(0);
                double de   = err_i(1);
                double lx   = l_obs(0);
                double ly   = l_obs(1);
                Vector6d Js_aux;
                Js_aux << + fgz2 * lx * gz,
                          + fgz2 * ly * gz,
                          - fgz2 * ( gx*lx + gy*ly ),
                          - fgz2 * ( gx*gy*lx + gy*gy*ly + gz*gz*ly ),
                          + fgz2 * ( gx*gx*lx + gz*gz*lx + gx*gy*ly ),
                          + fgz2 * ( gx*gz*ly - gy*gz*lx );
                // end point
                gx   = eP_(0);
                gy   = eP_(1);
                gz   = eP_(2);
                gz2  = gz*gz;
                fgz2 = cam->getFx() / std::max(SlamConfig::homogTh(),gz2);
                Vector6d Je_aux, J_aux;
                Je_aux << + fgz2 * lx * gz,
                          + fgz2 * ly * gz,
                          - fgz2 * ( gx*lx + gy*ly ),
                          - fgz2 * ( gx*gy*lx + gy*gy*ly + gz*gz*ly ),
                          + fgz2 * ( gx*gx*lx + gz*gz*lx + gx*gy*ly ),
                          + fgz2 * ( gx*gz*ly - gy*gz*lx );
                // jacobian
                J_aux = ( Js_aux * ds + Je_aux * de ) / std::max(SlamConfig::homogTh(),err_i_norm);
                // if employing robust cost function
                double  w = 1.0;
                double s2 = (*ls_it)->sigma2;
                w = robustWeightCauchy(err_i_norm) ;
                // update hessian, gradient, and error
                H_l += J_aux * J_aux.transpose() * w;
                g_l += J_aux * err_i_norm * w;
                e_l += err_i_norm * err_i_norm * w;
                N_l++;
            }
        }
        // sum H, g and err from both points and lines
        H = H_p + H_l;
        g = g_p + g_l;
        e = e_p + e_l;
        // normalize error
        e /= (N_l+N_p);

        // if the difference is very small stop
        if( abs(e-err_prev) < numeric_limits<double>::epsilon() || e < numeric_limits<double>::epsilon() )
            break;

        // solve
        ColPivHouseholderQR<MatrixXd> solver(H);
        x_inc = solver.solve( g );
        T_inc = T_inc * inverse_se3( expmap_se3(x_inc) );

        // if the parameter change is small stop
        if( x_inc.norm() < numeric_limits<double>::epsilon() )
            break;

        // update error
        err_prev = e;

    }

    // Remove outliers
    //---------------------------------------------------------------------------------------------
    // point observations
    for( vector<PointFeature*>::iterator pt_it = lc_points.begin(); pt_it != lc_points.end(); pt_it++ )
    {
        if( (*pt_it)->inlier )
        {
            Vector3d P_ = T_inc.block(0,0,3,3) * (*pt_it)->P + T_inc.col(3).head(3);
            Vector2d pl_proj = cam->projection( P_ );
            // projection error
            Vector2d err_i    = pl_proj - (*pt_it)->pl_obs;
            double s2 = (*pt_it)->sigma2;
            if( err_i.norm() > sqrt(7.815) )
                (*pt_it)->inlier = false;
        }
    }
    // line segments observations
    for( vector<LineFeature*>::iterator ls_it = lc_lines.begin(); ls_it != lc_lines.end(); ls_it++ )
    {
        if((*ls_it)->inlier)
        {
            Vector3d sP_ = T_inc.block(0,0,3,3) * (*ls_it)->sP + T_inc.col(3).head(3);
            Vector2d spl_proj = cam->projection( sP_ );
            Vector3d eP_ = T_inc.block(0,0,3,3) * (*ls_it)->eP + T_inc.col(3).head(3);
            Vector2d epl_proj = cam->projection( eP_ );
            Vector3d l_obs = (*ls_it)->le_obs;
            // projection error
            Vector2d err_i;
            err_i(0) = l_obs(0) * spl_proj(0) + l_obs(1) * spl_proj(1) + l_obs(2);
            err_i(1) = l_obs(0) * epl_proj(0) + l_obs(1) * epl_proj(1) + l_obs(2);
            double s2 = sqrt((*ls_it)->sigma2);
            if( err_i.norm() > sqrt(7.815) )
                (*ls_it)->inlier = false;
        }
    }

    // GN refinement
    //---------------------------------------------------------------------------------------------
    for( int iters = 0; iters < max_iters; iters++)
    {
        H = Matrix6d::Zero(); H_l = H; H_p = H;
        g = Vector6d::Zero(); g_l = g; g_p = g;
        e = 0.0; e_p = 0.0; e_l = 0.0;
        // point observations
        int N_p = 0;
        for( vector<PointFeature*>::iterator pt_it = lc_points.begin(); pt_it != lc_points.end(); pt_it++ )
        {
            if( (*pt_it)->inlier )
            {
                Vector3d P_ = T_inc.block(0,0,3,3) * (*pt_it)->P + T_inc.col(3).head(3);
                Vector2d pl_proj = cam->projection( P_ );
                // projection error
                Vector2d err_i    = pl_proj - (*pt_it)->pl_obs;
                double err_i_norm = err_i.norm();
                double gx   = P_(0);
                double gy   = P_(1);
                double gz   = P_(2);
                double gz2  = gz*gz;
                double fgz2 = cam->getFx() / std::max(SlamConfig::homogTh(),gz2);
                double dx   = err_i(0);
                double dy   = err_i(1);
                // jacobian
                Vector6d J_aux;
                J_aux << + fgz2 * dx * gz,
                         + fgz2 * dy * gz,
                         - fgz2 * ( gx*dx + gy*dy ),
                         - fgz2 * ( gx*gy*dx + gy*gy*dy + gz*gz*dy ),
                         + fgz2 * ( gx*gx*dx + gz*gz*dx + gx*gy*dy ),
                         + fgz2 * ( gx*gz*dy - gy*gz*dx );
                J_aux = J_aux / std::max(SlamConfig::homogTh(),err_i_norm);
                // if employing robust cost function
                double  w = 1.0;
                w = robustWeightCauchy(err_i_norm) ;
                // update hessian, gradient, and error
                H_p += J_aux * J_aux.transpose() * w;
                g_p += J_aux * err_i_norm * w;
                e_p += err_i_norm * err_i_norm * w;
                N_p++;
            }
        }
        // line segment features
        int N_l = 0;
        for( vector<LineFeature*>::iterator ls_it = lc_lines.begin(); ls_it != lc_lines.end(); ls_it++ )
        {
            if( (*ls_it)->inlier )
            {
                Vector3d sP_ = T_inc.block(0,0,3,3) * (*ls_it)->sP + T_inc.col(3).head(3);
                Vector2d spl_proj = cam->projection( sP_ );
                Vector3d eP_ = T_inc.block(0,0,3,3) * (*ls_it)->eP + T_inc.col(3).head(3);
                Vector2d epl_proj = cam->projection( eP_ );
                Vector3d l_obs = (*ls_it)->le_obs;
                // projection error
                Vector2d err_i;
                err_i(0) = l_obs(0) * spl_proj(0) + l_obs(1) * spl_proj(1) + l_obs(2);
                err_i(1) = l_obs(0) * epl_proj(0) + l_obs(1) * epl_proj(1) + l_obs(2);
                double err_i_norm = err_i.norm();
                // start point
                double gx   = sP_(0);
                double gy   = sP_(1);
                double gz   = sP_(2);
                double gz2  = gz*gz;
                double fgz2 = cam->getFx() / std::max(SlamConfig::homogTh(),gz2);
                double ds   = err_i(0);
                double de   = err_i(1);
                double lx   = l_obs(0);
                double ly   = l_obs(1);
                Vector6d Js_aux;
                Js_aux << + fgz2 * lx * gz,
                          + fgz2 * ly * gz,
                          - fgz2 * ( gx*lx + gy*ly ),
                          - fgz2 * ( gx*gy*lx + gy*gy*ly + gz*gz*ly ),
                          + fgz2 * ( gx*gx*lx + gz*gz*lx + gx*gy*ly ),
                          + fgz2 * ( gx*gz*ly - gy*gz*lx );
                // end point
                gx   = eP_(0);
                gy   = eP_(1);
                gz   = eP_(2);
                gz2  = gz*gz;
                fgz2 = cam->getFx() / std::max(SlamConfig::homogTh(),gz2);
                Vector6d Je_aux, J_aux;
                Je_aux << + fgz2 * lx * gz,
                          + fgz2 * ly * gz,
                          - fgz2 * ( gx*lx + gy*ly ),
                          - fgz2 * ( gx*gy*lx + gy*gy*ly + gz*gz*ly ),
                          + fgz2 * ( gx*gx*lx + gz*gz*lx + gx*gy*ly ),
                          + fgz2 * ( gx*gz*ly - gy*gz*lx );
                // jacobian
                J_aux = ( Js_aux * ds + Je_aux * de ) / std::max(SlamConfig::homogTh(),err_i_norm);
                // if employing robust cost function
                double  w = 1.0;
                double s2 = (*ls_it)->sigma2;
                w = robustWeightCauchy(err_i_norm) ;
                // update hessian, gradient, and error
                H_l += J_aux * J_aux.transpose() * w;
                g_l += J_aux * err_i_norm * w;
                e_l += err_i_norm * err_i_norm * w;
                N_l++;
            }
        }
        // sum H, g and err from both points and lines
        H = H_p + H_l;
        g = g_p + g_l;
        e = e_p + e_l;
        // normalize error
        e /= (N_l+N_p);
        // if the difference is very small stop
        if( abs(e-err_prev) < numeric_limits<double>::epsilon() || e < numeric_limits<double>::epsilon() )
            break;
        // solve
        ColPivHouseholderQR<MatrixXd> solver(H);
        x_inc = solver.solve( g );
        T_inc = T_inc * inverse_se3( expmap_se3(x_inc) );
        // if the parameter change is small stop
        if( x_inc.norm() < numeric_limits<double>::epsilon() )
            break;
        // update error
        err_prev = e;
    }

    x_inc = logmap_se3(T_inc);

    // Check whether it is Loop Closure or not
    //---------------------------------------------------------------------------------------------
    // Residue value
    bool lc_res = ( e < SlamConfig::lcRes() );

    // Uncertainty value
    Matrix6d DT_cov;
    DT_cov = H.inverse();
    SelfAdjointEigenSolver<Matrix6d> eigensolver(DT_cov);
    Vector6d DT_cov_eig = eigensolver.eigenvalues();
    bool lc_unc = ( DT_cov_eig(5) < SlamConfig::lcUnc() );

    // Ratio of outliers
    int N = lc_points.size() + lc_lines.size(), N_inl = 0;
    for( vector<PointFeature*>::iterator pt_it = lc_points.begin(); pt_it != lc_points.end(); pt_it++ )
    {
        if( (*pt_it)->inlier )
            N_inl++;
    }
    for( vector<LineFeature*>::iterator  ls_it = lc_lines.begin();  ls_it != lc_lines.end();  ls_it++ )
    {
        if( (*ls_it)->inlier )
            N_inl++;
    }
    double ratio_inliers = double(N_inl) / double(N);
    bool lc_inl = ( ratio_inliers > SlamConfig::lcInl() );

    lc_inl = true;

    // Geometry condition
    double t = x_inc.head(3).norm();
    double r = x_inc.tail(3).norm() * 180.f / CV_PI;
    bool lc_trs = ( t < SlamConfig::lcTrs() );
    bool lc_rot = ( r < SlamConfig::lcRot() );

    // Decision
    if( lc_res && lc_unc && lc_inl && lc_trs && lc_rot )
    {
        // erase outliers from      lc_pt_idx & lc_ls_idx
        int iter = 0;
        vector<Vector4i> lc_pt_idx_, lc_ls_idx_;
        vector<PointFeature*> lc_points_;
        vector<LineFeature*>  lc_lines_;
        for( vector<PointFeature*>::iterator pt_it = lc_points.begin(); pt_it != lc_points.end(); pt_it++, iter++ )
        {
            if( (*pt_it)!= NULL )
            {
                if( (*pt_it)->inlier )
                {
                    lc_points_.push_back( (*pt_it) );
                    lc_pt_idx_.push_back( lc_pt_idx[iter] );
                }
            }
        }
        iter = 0;
        for( vector<LineFeature*>::iterator ls_it = lc_lines.begin(); ls_it != lc_lines.end(); ls_it++, iter++ )
        {
            if( (*ls_it)!= NULL )
            {
                if( (*ls_it)->inlier )
                {
                    lc_lines_.push_back( (*ls_it) );
                    lc_ls_idx_.push_back( lc_ls_idx[iter] );
                }
            }
        }
        lc_pt_idx.clear();
        lc_pt_idx = lc_pt_idx_;
        lc_ls_idx.clear();
        lc_ls_idx = lc_ls_idx_;
        lc_points.clear();
        lc_points = lc_points_;
        lc_lines.clear();
        lc_points = lc_points_;
        // assign pose increment
        pose_inc = logmap_se3( inverse_se3( expmap_se3(x_inc) ) );
        return true;
    }
    else
        return false;

}

bool MapHandler::loopClosureOptimizationEssGraphG2O()
{
    typedef g2o::BlockSolver<g2o::BlockSolverTraits<6,3>> BlockSolverType;
    typedef g2o::LinearSolverCholmod<BlockSolverType::PoseMatrixType> LinearSolverType;
    auto solver = new g2o::OptimizationAlgorithmLevenberg(g2o::make_unique<BlockSolverType>(
            g2o::make_unique<LinearSolverType>()));
    g2o::SparseOptimizer optimizer;

    // define G2O variables
//    g2o::SparseOptimizer optimizer;
//    optimizer.setVerbose(true);
//    g2o::BlockSolver_6_3::LinearSolverType* linearSolver;
//    linearSolver = new g2o::LinearSolverCholmod<g2o::BlockSolver_6_3::PoseMatrixType>();
//    g2o::BlockSolver_6_3* solver_ptr = new g2o::BlockSolver_6_3(linearSolver);
//    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    solver->setUserLambdaInit(1e-10);
    optimizer.setAlgorithm(solver);

    // select min and max KF indices
    int kf_prev_idx =  2 * max_kf_idx;
    int kf_curr_idx = -1;
    for( auto it = lc_idxs.begin(); it != lc_idxs.end(); it++)
    {
        if( (*it)(0) < kf_prev_idx )
            kf_prev_idx = (*it)(0);
        if( (*it)(1) > kf_curr_idx )
            kf_curr_idx = (*it)(1);
    }

    // grab the KFs included in the optimization
    vector<int> kf_list;
    for( int i = kf_prev_idx; i <= kf_curr_idx; i++)
    {
        if( map_keyframes[i] != NULL )
        {
            // check if it is a LC vertex
            bool is_lc_i = false;
            bool is_lc_j = false;
            int id = 0;
            for( auto it = lc_idxs.begin(); it != lc_idxs.end(); it++, id++ )
            {
                if( (*it)(0) == i )
                {
                    is_lc_i = true;
                    break;
                }
                if( (*it)(1) == i )
                {
                    is_lc_j = true;
                    break;
                }
            }
            kf_list.push_back(i);
            // create SE3 vertex
            g2o::VertexSE3* v_se3 = new g2o::VertexSE3();
            v_se3->setId(i);
            v_se3->setMarginalized( false );
            if( is_lc_j )
            {
                // update pose of LC vertex
                v_se3->setFixed(true);
                v_se3->setEstimate( g2o::SE3Quat::exp( reverse_se3(logmap_se3( (expmap_se3(lc_pose_list[id])) * map_keyframes[lc_idxs[id](0)]->T_kf_w )) ) );
            }
            else
            {
                v_se3->setEstimate( g2o::SE3Quat::exp( reverse_se3(map_keyframes[i]->x_kf_w) ) );
                if( is_lc_i || i == 0 )
                    v_se3->setFixed(true);
            }
            optimizer.addVertex( v_se3 );
        }
    }

    // introduce edges
    for( int i = kf_prev_idx; i <= kf_curr_idx; i++ )
    {
        for( int j = i+1; j <= kf_curr_idx; j++ )
        {
            if( map_keyframes[i] != NULL && map_keyframes[j] != NULL &&
                ( full_graph[i][j] >= SlamConfig::minLMEssGraph() || abs(i-j) == 1  ) )
            {
                // kf2kf constraint
                Matrix4d T_ji_constraint = inverse_se3( map_keyframes[i]->T_kf_w ) * map_keyframes[j]->T_kf_w;
                // add edge
                g2o::EdgeSE3* e_se3 = new g2o::EdgeSE3();
                e_se3->setVertex( 0, optimizer.vertex(i) );
                e_se3->setVertex( 1, optimizer.vertex(j) );
                Vector6d x;
                x = reverse_se3(logmap_se3(T_ji_constraint) );
                e_se3->setMeasurement( g2o::SE3Quat::exp(x) );
                e_se3->setInformation( Matrix6d::Identity() );
                optimizer.addEdge( e_se3 );
            }
        }
    }

    // introduce loop closure edges
    int id = 0;
    for( auto it = lc_idx_list.begin(); it != lc_idx_list.end(); it++, id++ )
    {
        // add edge
        g2o::EdgeSE3* e_se3 = new g2o::EdgeSE3();
        e_se3->setVertex( 0, optimizer.vertex((*it)(0)) );
        e_se3->setVertex( 1, optimizer.vertex((*it)(1)) );
        Vector6d x;
        x = reverse_se3( lc_pose_list[id] );
        e_se3->setMeasurement( g2o::SE3Quat::exp(x) );
        e_se3->information() = Matrix6d::Identity();
        optimizer.addEdge( e_se3 );
    }

    // optimize graph
    optimizer.initializeOptimization();
    optimizer.computeInitialGuess();
    optimizer.computeActiveErrors();
    optimizer.optimize(SlamConfig::maxItersPGO());

    // recover pose and update map
    Matrix4d Tkfw_corr;
    for( auto kf_it = kf_list.begin(); kf_it != kf_list.end(); kf_it++)
    {
        g2o::VertexSE3* v_se3 = static_cast<g2o::VertexSE3*>(optimizer.vertex( (*kf_it) ));
        g2o::SE3Quat Tiw_corr =  v_se3->estimateAsSE3Quat();
        Vector6d x;
        Matrix4d Tkfw, Tkfw_prev;
        x = reverse_se3(Tiw_corr.log());
        Tkfw = expmap_se3( x );
        Tkfw_prev = map_keyframes[ (*kf_it) ]->T_kf_w;
        map_keyframes[ (*kf_it) ]->T_kf_w = Tkfw;
        map_keyframes[ (*kf_it) ]->x_kf_w = logmap_se3(Tkfw);
        // update map
        Tkfw_corr = Tkfw * inverse_se3( Tkfw_prev );
        for( auto it = map_points_kf_idx.at((*kf_it)).begin(); it != map_points_kf_idx.at((*kf_it)).end(); it++ )
        {
           if( map_points[(*it)] != NULL )
           {
               // update 3D position
               Vector3d point3D = map_points[(*it)]->point3D;
               map_points[(*it)]->point3D = Tkfw_corr.block(0,0,3,3) * point3D + Tkfw_corr.block(0,3,3,1);
               // update direction of observation
               Vector3d obs_dir = map_points[(*it)]->med_obs_dir;
               map_points[(*it)]->med_obs_dir = Tkfw_corr.block(0,0,3,3) * obs_dir + Tkfw_corr.block(0,3,3,1);
               // update direction of each observation
               for( auto dir_it = map_points[(*it)]->dir_list.begin(); dir_it != map_points[(*it)]->dir_list.end(); dir_it++)
               {
                   Vector3d dir_list_ = (*dir_it);
                   (*dir_it) = Tkfw_corr.block(0,0,3,3) * dir_list_ + Tkfw_corr.block(0,3,3,1);
               }
           }
        }
        for( auto it = map_lines_kf_idx.at((*kf_it)).begin(); it != map_lines_kf_idx.at((*kf_it)).end(); it++ )
        {
           if( map_lines[(*it)] != NULL )
           {
               // update 3D position
               Vector3d sP3D = map_lines[(*it)]->line3D.head(3);
               Vector3d eP3D = map_lines[(*it)]->line3D.tail(3);
               map_lines[(*it)]->line3D.head(3) = Tkfw_corr.block(0,0,3,3) * sP3D + Tkfw_corr.block(0,3,3,1);
               map_lines[(*it)]->line3D.tail(3) = Tkfw_corr.block(0,0,3,3) * eP3D + Tkfw_corr.block(0,3,3,1);
               // update direction of observation
               Vector3d obs_dir = map_lines[(*it)]->med_obs_dir;
               map_lines[(*it)]->med_obs_dir = Tkfw_corr.block(0,0,3,3) * obs_dir + Tkfw_corr.block(0,3,3,1);
               // update direction of each observation
               for( auto dir_it = map_lines[(*it)]->dir_list.begin(); dir_it != map_lines[(*it)]->dir_list.end(); dir_it++)
               {
                   Vector3d dir_list_ = (*dir_it);
                   (*dir_it) = Tkfw_corr.block(0,0,3,3) * dir_list_ + Tkfw_corr.block(0,3,3,1);
               }
           }
        }
    }

    // update pose and map of the rest of frames
    for( int i = kf_curr_idx + 1; i < map_keyframes.size(); i++ )
    {
        // update pose
        map_keyframes[i]->T_kf_w = Tkfw_corr * map_keyframes[i]->T_kf_w;
        map_keyframes[i]->x_kf_w = logmap_se3(map_keyframes[i]->T_kf_w);
        // update landmarks
        for( auto it = map_points_kf_idx.at(i).begin(); it != map_points_kf_idx.at(i).end(); it++ )
        {
           if( map_points[(*it)] != NULL )
           {
               // update 3D position
               Vector3d point3D = map_points[(*it)]->point3D;
               map_points[(*it)]->point3D = Tkfw_corr.block(0,0,3,3) * point3D + Tkfw_corr.block(0,3,3,1);
               // update direction of observation
               Vector3d obs_dir = map_points[(*it)]->med_obs_dir;
               map_points[(*it)]->med_obs_dir = Tkfw_corr.block(0,0,3,3) * obs_dir + Tkfw_corr.block(0,3,3,1);
               // update direction of each observation
               for( auto dir_it = map_points[(*it)]->dir_list.begin(); dir_it != map_points[(*it)]->dir_list.end(); dir_it++)
               {
                   Vector3d dir_list_ = (*dir_it);
                   (*dir_it) = Tkfw_corr.block(0,0,3,3) * dir_list_ + Tkfw_corr.block(0,3,3,1);
               }
           }
        }
        for( auto it = map_lines_kf_idx.at(i).begin(); it != map_lines_kf_idx.at(i).end(); it++ )
        {
           if( map_lines[(*it)] != NULL )
           {
               // update 3D position
               Vector3d sP3D = map_lines[(*it)]->line3D.head(3);
               Vector3d eP3D = map_lines[(*it)]->line3D.tail(3);
               map_lines[(*it)]->line3D.head(3) = Tkfw_corr.block(0,0,3,3) * sP3D + Tkfw_corr.block(0,3,3,1);
               map_lines[(*it)]->line3D.tail(3) = Tkfw_corr.block(0,0,3,3) * eP3D + Tkfw_corr.block(0,3,3,1);
               // update direction of observation
               Vector3d obs_dir = map_lines[(*it)]->med_obs_dir;
               map_lines[(*it)]->med_obs_dir = Tkfw_corr.block(0,0,3,3) * obs_dir + Tkfw_corr.block(0,3,3,1);
               // update direction of each observation
               for( auto dir_it = map_lines[(*it)]->dir_list.begin(); dir_it != map_lines[(*it)]->dir_list.end(); dir_it++)
               {
                   Vector3d dir_list_ = (*dir_it);
                   (*dir_it) = Tkfw_corr.block(0,0,3,3) * dir_list_ + Tkfw_corr.block(0,3,3,1);
               }
           }
        }
    }

    // mark as optimized the lc_idx_list edges
    for( auto it = lc_idx_list.begin(); it != lc_idx_list.end(); it++)
        (*it)(2) = 0;

    // fuse local map from both sides of the loop and update graphs
    loopClosureFuseLandmarks();

    lc_state = LC_IDLE;

    return true;
}

bool MapHandler::loopClosureOptimizationCovGraphG2O()
{
    typedef g2o::BlockSolver<g2o::BlockSolverTraits<6,3>> BlockSolverType;
    typedef g2o::LinearSolverCholmod<BlockSolverType::PoseMatrixType> LinearSolverType;
    auto solver = new g2o::OptimizationAlgorithmLevenberg(g2o::make_unique<BlockSolverType>(
            g2o::make_unique<LinearSolverType>()));
    g2o::SparseOptimizer optimizer;
    // define G2O variables
//    g2o::SparseOptimizer optimizer;
    optimizer.setVerbose(false);
//    g2o::BlockSolver_6_3::LinearSolverType* linearSolver;
//    linearSolver = new g2o::LinearSolverCholmod<g2o::BlockSolver_6_3::PoseMatrixType>();
//    g2o::BlockSolver_6_3* solver_ptr = new g2o::BlockSolver_6_3(linearSolver);
//    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    solver->setUserLambdaInit(1e-10);
    optimizer.setAlgorithm(solver);

    // select min and max KF indices
    int kf_prev_idx =  2 * max_kf_idx;
    int kf_curr_idx = -1;
    for( auto it = lc_idx_list.begin(); it != lc_idx_list.end(); it++)
    {
        if( (*it)(0) < kf_prev_idx )
            kf_prev_idx = (*it)(0);
        if( (*it)(1) > kf_curr_idx )
            kf_curr_idx = (*it)(1);
    }
    kf_prev_idx = 0;

    // grab the KFs included in the optimization
    vector<int> kf_list;
    for( int i = kf_prev_idx; i <= kf_curr_idx; i++)
    {
        if( map_keyframes[i] != NULL )
        {
            // check if it is a LC vertex
            bool is_lc_i = false;
            bool is_lc_j = false;
            int id = 0;
            for( auto it = lc_idx_list.begin(); it != lc_idx_list.end(); it++, id++ )
            {
                if( (*it)(0) == i )
                {
                    is_lc_i = true;
                    break;
                }
                if( (*it)(1) == i )
                {
                    is_lc_j = true;
                    break;
                }
            }
            kf_list.push_back(i);
            // create SE3 vertex
            g2o::VertexSE3* v_se3 = new g2o::VertexSE3();
            v_se3->setId(i);
            v_se3->setMarginalized( false );
            if( is_lc_j )
            {
                // update pose of LC vertex
                v_se3->setFixed(false);
                v_se3->setEstimate( g2o::SE3Quat::exp( reverse_se3(logmap_se3( (expmap_se3(lc_pose_list[id])) * map_keyframes[lc_idx_list[id](0)]->T_kf_w )) ) );
            }
            else
            {
                v_se3->setEstimate( g2o::SE3Quat::exp( reverse_se3(map_keyframes[i]->x_kf_w) ) );
                if( i == 0 )
                    v_se3->setFixed(true);
            }
            optimizer.addVertex( v_se3 );
        }
    }

    // introduce edges
    for( int i = kf_prev_idx; i <= kf_curr_idx; i++ )
    {
        for( int j = i+1; j <= kf_curr_idx; j++ )
        {
            if( map_keyframes[i] != NULL && map_keyframes[j] != NULL &&
                ( full_graph[i][j] >= SlamConfig::minLMEssGraph() || full_graph[i][j] >= SlamConfig::minLMCovGraph() || abs(i-j) == 1  ) )
            {
                // kf2kf constraint
                Matrix4d T_ji_constraint = inverse_se3( map_keyframes[i]->T_kf_w ) * map_keyframes[j]->T_kf_w;
                // add edge
                g2o::EdgeSE3* e_se3 = new g2o::EdgeSE3();
                e_se3->setVertex( 0, optimizer.vertex(i) );
                e_se3->setVertex( 1, optimizer.vertex(j) );
                Vector6d x;
                x = reverse_se3(logmap_se3(T_ji_constraint) );
                e_se3->setMeasurement( g2o::SE3Quat::exp(x) );
                e_se3->setInformation( Matrix6d::Identity() );
                optimizer.addEdge( e_se3 );
            }
        }
    }

    // introduce loop closure edges
    int id = 0;
    for( auto it = lc_idx_list.begin(); it != lc_idx_list.end(); it++, id++ )
    {
        // add edge
        g2o::EdgeSE3* e_se3 = new g2o::EdgeSE3();
        e_se3->setVertex( 0, optimizer.vertex((*it)(0)) );
        e_se3->setVertex( 1, optimizer.vertex((*it)(1)) );
        Vector6d x;
        x = reverse_se3( lc_pose_list[id] );
        e_se3->setMeasurement( g2o::SE3Quat::exp(x) );
        e_se3->information() = Matrix6d::Identity();
        optimizer.addEdge( e_se3 );
    }

    // optimize graph
    optimizer.initializeOptimization();
    optimizer.computeInitialGuess();
    optimizer.computeActiveErrors();
    optimizer.optimize(SlamConfig::maxItersPGO());

    // recover pose and update map
    Matrix4d Tkfw_corr;
    for( auto kf_it = kf_list.begin(); kf_it != kf_list.end(); kf_it++)
    {
        g2o::VertexSE3* v_se3 = static_cast<g2o::VertexSE3*>(optimizer.vertex( (*kf_it) ));
        g2o::SE3Quat Tiw_corr =  v_se3->estimateAsSE3Quat();
        Vector6d x;
        Matrix4d Tkfw, Tkfw_prev;
        x = reverse_se3(Tiw_corr.log());
        Tkfw = expmap_se3( x );
        Tkfw_prev = map_keyframes[ (*kf_it) ]->T_kf_w;
        map_keyframes[ (*kf_it) ]->T_kf_w = Tkfw;
        map_keyframes[ (*kf_it) ]->x_kf_w = logmap_se3(Tkfw);
        // update map
        Tkfw_corr = Tkfw * inverse_se3( Tkfw_prev );
        for( auto it = map_points_kf_idx.at((*kf_it)).begin(); it != map_points_kf_idx.at((*kf_it)).end(); it++ )
        {
           if( map_points[(*it)] != NULL )
           {
               // update 3D position
               Vector3d point3D = map_points[(*it)]->point3D;
               map_points[(*it)]->point3D = Tkfw_corr.block(0,0,3,3) * point3D + Tkfw_corr.block(0,3,3,1);
               // update direction of observation
               Vector3d obs_dir = map_points[(*it)]->med_obs_dir;
               map_points[(*it)]->med_obs_dir = Tkfw_corr.block(0,0,3,3) * obs_dir + Tkfw_corr.block(0,3,3,1);
               // update direction of each observation
               for( auto dir_it = map_points[(*it)]->dir_list.begin(); dir_it != map_points[(*it)]->dir_list.end(); dir_it++)
               {
                   Vector3d dir_list_ = (*dir_it);
                   (*dir_it) = Tkfw_corr.block(0,0,3,3) * dir_list_ + Tkfw_corr.block(0,3,3,1);
               }
           }
        }
        for( auto it = map_lines_kf_idx.at((*kf_it)).begin(); it != map_lines_kf_idx.at((*kf_it)).end(); it++ )
        {
           if( map_lines[(*it)] != NULL )
           {
               // update 3D position
               Vector3d sP3D = map_lines[(*it)]->line3D.head(3);
               Vector3d eP3D = map_lines[(*it)]->line3D.tail(3);
               map_lines[(*it)]->line3D.head(3) = Tkfw_corr.block(0,0,3,3) * sP3D + Tkfw_corr.block(0,3,3,1);
               map_lines[(*it)]->line3D.tail(3) = Tkfw_corr.block(0,0,3,3) * eP3D + Tkfw_corr.block(0,3,3,1);
               // update direction of observation
               Vector3d obs_dir = map_lines[(*it)]->med_obs_dir;
               map_lines[(*it)]->med_obs_dir = Tkfw_corr.block(0,0,3,3) * obs_dir + Tkfw_corr.block(0,3,3,1);
               // update direction of each observation
               for( auto dir_it = map_lines[(*it)]->dir_list.begin(); dir_it != map_lines[(*it)]->dir_list.end(); dir_it++)
               {
                   Vector3d dir_list_ = (*dir_it);
                   (*dir_it) = Tkfw_corr.block(0,0,3,3) * dir_list_ + Tkfw_corr.block(0,3,3,1);
               }
           }
        }
    }

    // update pose and map of the rest of frames
    for( int i = kf_curr_idx + 1; i < map_keyframes.size(); i++ )
    {
        // update pose
        map_keyframes[i]->T_kf_w = Tkfw_corr * map_keyframes[i]->T_kf_w;
        map_keyframes[i]->x_kf_w = logmap_se3(map_keyframes[i]->T_kf_w);
        // update landmarks
        for( auto it = map_points_kf_idx.at(i).begin(); it != map_points_kf_idx.at(i).end(); it++ )
        {
           if( map_points[(*it)] != NULL )
           {
               // update 3D position
               Vector3d point3D = map_points[(*it)]->point3D;
               map_points[(*it)]->point3D = Tkfw_corr.block(0,0,3,3) * point3D + Tkfw_corr.block(0,3,3,1);
               // update direction of observation
               Vector3d obs_dir = map_points[(*it)]->med_obs_dir;
               map_points[(*it)]->med_obs_dir = Tkfw_corr.block(0,0,3,3) * obs_dir + Tkfw_corr.block(0,3,3,1);
               // update direction of each observation
               for( auto dir_it = map_points[(*it)]->dir_list.begin(); dir_it != map_points[(*it)]->dir_list.end(); dir_it++)
               {
                   Vector3d dir_list_ = (*dir_it);
                   (*dir_it) = Tkfw_corr.block(0,0,3,3) * dir_list_ + Tkfw_corr.block(0,3,3,1);
               }
           }
        }
        for( auto it = map_lines_kf_idx.at(i).begin(); it != map_lines_kf_idx.at(i).end(); it++ )
        {
           if( map_lines[(*it)] != NULL )
           {
               // update 3D position
               Vector3d sP3D = map_lines[(*it)]->line3D.head(3);
               Vector3d eP3D = map_lines[(*it)]->line3D.tail(3);
               map_lines[(*it)]->line3D.head(3) = Tkfw_corr.block(0,0,3,3) * sP3D + Tkfw_corr.block(0,3,3,1);
               map_lines[(*it)]->line3D.tail(3) = Tkfw_corr.block(0,0,3,3) * eP3D + Tkfw_corr.block(0,3,3,1);
               // update direction of observation
               Vector3d obs_dir = map_lines[(*it)]->med_obs_dir;
               map_lines[(*it)]->med_obs_dir = Tkfw_corr.block(0,0,3,3) * obs_dir + Tkfw_corr.block(0,3,3,1);
               // update direction of each observation
               for( auto dir_it = map_lines[(*it)]->dir_list.begin(); dir_it != map_lines[(*it)]->dir_list.end(); dir_it++)
               {
                   Vector3d dir_list_ = (*dir_it);
                   (*dir_it) = Tkfw_corr.block(0,0,3,3) * dir_list_ + Tkfw_corr.block(0,3,3,1);
               }
           }
        }
    }

    // mark as optimized the lc_idx_list edges
    for( auto it = lc_idx_list.begin(); it != lc_idx_list.end(); it++)
        (*it)(2) = 0;

    // fuse local map from both sides of the loop and update graphs
    loopClosureFuseLandmarks();

    lc_state = LC_IDLE;

    return true;
}

void MapHandler::loopClosureFuseLandmarks()
{

    // point matches
    int lc_idx = 0;
    for( auto idx_it = lc_pt_idxs.begin(); idx_it != lc_pt_idxs.end(); idx_it++, lc_idx++ )
    {
        if( lc_idx_list[lc_idx](2) == 1 )   // if not already optimized
        {
            int kf_prev_idx = lc_idx_list[lc_idx](0);
            int kf_curr_idx = lc_idx_list[lc_idx](1);
            for( auto lm_it = (*idx_it).begin(); lm_it != (*idx_it).end(); lm_it++ )
            {
                // grab indices
                int lm_idx0  = (*lm_it)(0);
                int lm_ldx0  = (*lm_it)(1); // lr_qdx
                int lm_idx1  = (*lm_it)(2);
                int lm_ldx1  = (*lm_it)(3); // lr_tdx
                // if the LM exists just once, add observation
                if( lm_idx0 == -1 && lm_idx1 != -1 )
                {
                    if( map_keyframes[kf_prev_idx]->stereo_frame->stereo_pt[lm_ldx0] != NULL && map_points[lm_idx1] != NULL )
                    {
                        map_keyframes[kf_prev_idx]->stereo_frame->stereo_pt[lm_ldx0]->idx = lm_idx1;
                        Vector3d dir  = map_keyframes[kf_prev_idx]->stereo_frame->stereo_pt[lm_ldx0]->P / map_keyframes[kf_prev_idx]->stereo_frame->stereo_pt[lm_ldx0]->P.norm();
                        map_points[lm_idx1]->addMapPointObservation( map_keyframes[kf_prev_idx]->stereo_frame->pdesc_l.row(lm_ldx0), map_keyframes[kf_prev_idx]->kf_idx, map_keyframes[kf_prev_idx]->stereo_frame->stereo_pt[lm_ldx0]->pl, dir );
                        // increase full graph for each KF that has already observed this LM
                        for( auto kf_it = map_points[lm_idx1]->kf_obs_list.begin(); kf_it != map_points[lm_idx1]->kf_obs_list.end(); kf_it++)
                        {
                            full_graph[(*kf_it)][kf_curr_idx]++;
                            full_graph[kf_curr_idx][(*kf_it)]++;
                        }
                    }
                }
                if( lm_idx0 != -1 && lm_idx1 == -1 )
                {
                    if( map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1] != NULL && map_points[lm_idx0] != NULL )
                    {
                        map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1]->idx = lm_idx0;
                        Vector3d dir  = map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1]->P / map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1]->P.norm();
                        map_points[lm_idx0]->addMapPointObservation( map_keyframes[kf_curr_idx]->stereo_frame->pdesc_l.row(lm_ldx1),
                                                                     map_keyframes[kf_curr_idx]->kf_idx,
                                                                     map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1]->pl, dir );
                        // increase full graph for each KF that has already observed this LM
                        for( auto kf_it = map_points[lm_idx0]->kf_obs_list.begin(); kf_it != map_points[lm_idx0]->kf_obs_list.end(); kf_it++)
                        {
                            full_graph[(*kf_it)][kf_prev_idx]++;
                            full_graph[kf_prev_idx][(*kf_it)]++;
                        }
                    }
                }
                // if not, create LM and add both observations
                if( lm_idx0 == -1 && lm_idx1 == -1 )
                {
                    if( map_keyframes[kf_prev_idx]->stereo_frame->stereo_pt[lm_ldx0] != NULL && map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1] != NULL )
                    {
                        // assign indices
                        map_keyframes[kf_prev_idx]->stereo_frame->stereo_pt[lm_ldx0]->idx = max_pt_idx;
                        map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1]->idx = max_pt_idx;
                        // create new 3D landmark with the observation from previous KF
                        Matrix4d Tfw = ( map_keyframes[kf_prev_idx]->T_kf_w );
                        Vector3d P3d = Tfw.block(0,0,3,3) * map_keyframes[kf_prev_idx]->stereo_frame->stereo_pt[lm_ldx0]->P + Tfw.col(3).head(3);
                        Vector3d dir = P3d / P3d.norm();
                        MapPoint* map_point = new MapPoint(max_pt_idx,P3d,map_keyframes[kf_prev_idx]->stereo_frame->pdesc_l.row(lm_ldx0),map_keyframes[kf_prev_idx]->kf_idx,map_keyframes[kf_prev_idx]->stereo_frame->stereo_pt[lm_ldx0]->pl,dir);
                        // add new 3D landmark to kf_idx where it was first observed
                        map_points_kf_idx.at( kf_prev_idx ).push_back( max_pt_idx );
                        // add observation of the 3D landmark from current KF
                        P3d = map_keyframes[kf_curr_idx]->T_kf_w.block(0,0,3,3) *  map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1]->P + map_keyframes[kf_curr_idx]->T_kf_w.col(3).head(3);
                        dir = P3d / P3d.norm();
                        map_point->addMapPointObservation(map_keyframes[kf_curr_idx]->stereo_frame->pdesc_l.row(lm_ldx1),map_keyframes[kf_curr_idx]->kf_idx,map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1]->pl,dir);
                        // add 3D landmark to map
                        map_points.push_back(map_point);
                        // update full graph (new feature)
                        max_pt_idx++;
                        full_graph[kf_prev_idx][kf_curr_idx]++;
                        full_graph[kf_curr_idx][kf_prev_idx]++;
                    }
                }
                // if the LM observed is different in each KF, then fuse them and erase the old one
                if( lm_idx0 != -1 && lm_idx1 != -1 )
                {
                    if( map_points[lm_idx0] != NULL && map_points[lm_idx1] != NULL
                        && map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1] != NULL )
                    {
                        int Nobs_lm_prev = map_points[lm_idx0]->kf_obs_list.size();
                        // fuse LMs while updating the full graph
                        int iter = 0;
                        for( auto it = map_points[lm_idx1]->desc_list.begin(); it != map_points[lm_idx1]->desc_list.end(); it++, iter++)
                        {
                            // concatenate desc, obs, dir, and kf_obs lists
                            map_points[lm_idx0]->desc_list.push_back(   (*it) );
                            map_points[lm_idx0]->obs_list.push_back(    map_points[lm_idx1]->obs_list[iter]    );
                            map_points[lm_idx0]->dir_list.push_back(    map_points[lm_idx1]->dir_list[iter]    );
                            map_points[lm_idx0]->kf_obs_list.push_back( map_points[lm_idx1]->kf_obs_list[iter] );
                            // update full graph
                            int jdx = map_points[lm_idx1]->kf_obs_list[iter];
                            for( int i = 0; i < Nobs_lm_prev; i++ )
                            {
                                int idx = map_points[lm_idx0]->kf_obs_list[i];
                                full_graph[idx][jdx]++;
                                full_graph[jdx][idx]++;
                            }
                            // update average descriptor and direction of observation
                            map_points[lm_idx0]->updateAverageDescDir();
                            // change idx in stereo_pt
                            map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1]->idx = lm_idx0;
                        }
                        // remove from map_points_kf_idx
                        iter = 0;
                        int kf_lm_obs = map_points[lm_idx1]->kf_obs_list[0];
                        for( auto it = map_points_kf_idx.at(kf_lm_obs).begin(); it != map_points_kf_idx.at(kf_lm_obs).end(); it++, iter++)
                        {
                            if( (*it) == lm_idx1 )
                            {
                                map_points_kf_idx.at(kf_lm_obs).erase( map_points_kf_idx.at(kf_lm_obs).begin() + iter );
                                break;
                            }
                        }
                        // erase old landmark
                        delete map_points[lm_idx1];
                        map_points[lm_idx1] = nullptr;
                    }
                }
            }
        }
    }

    // line segment matches
    lc_idx = 0;
    for( auto idx_it = lc_ls_idxs.begin(); idx_it != lc_ls_idxs.end(); idx_it++, lc_idx++ )
    {
        if( lc_idx_list[lc_idx](2) == 1 )   // if not already optimized
        {
            int kf_prev_idx = lc_idx_list[lc_idx](0);
            int kf_curr_idx = lc_idx_list[lc_idx](1);
            for( auto lm_it = (*idx_it).begin(); lm_it != (*idx_it).end(); lm_it++ )
            {
                // grab indices
                int lm_idx0  = (*lm_it)(0);
                int lm_ldx0  = (*lm_it)(1); // lr_qdx
                int lm_idx1  = (*lm_it)(2);
                int lm_ldx1  = (*lm_it)(3); // lr_tdx
                // if the LM exists just once, add observation
                if( lm_idx0 == -1 && lm_idx1 != -1 )
                {
                    if( map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0] != NULL && map_lines[lm_idx1] != NULL )
                    {
                        map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->idx = lm_idx1;
                        Vector3d dir  = ( map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->sP + map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->eP )
                                      / ( map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->sP + map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->eP ).norm();
                        Vector4d pts;
                        pts.head(2) = map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->spl_obs;
                        pts.tail(2) = map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->epl_obs;
                        map_lines[lm_idx1]->addMapLineObservation( map_keyframes[kf_prev_idx]->stereo_frame->ldesc_l.row(lm_ldx0),
                                                                   map_keyframes[kf_prev_idx]->kf_idx,
                                                                   map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->le,
                                                                   dir, pts );
                        // increase full graph for each KF that has already observed this LM
                        for( auto kf_it = map_lines[lm_idx1]->kf_obs_list.begin(); kf_it != map_lines[lm_idx1]->kf_obs_list.end(); kf_it++)
                        {
                            full_graph[(*kf_it)][kf_curr_idx]++;
                            full_graph[kf_curr_idx][(*kf_it)]++;
                        }
                    }
                }
                if( lm_idx0 != -1 && lm_idx1 == -1 )
                {
                    if( map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1] != NULL && map_lines[lm_idx0] != NULL )
                    {
                        map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->idx = lm_idx0;
                        Vector3d dir  = (map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->sP + map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->eP)
                                / (map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->sP+map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->eP).norm();
                        Vector4d pts;
                        pts.head(2) = map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->spl_obs;
                        pts.tail(2) = map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->epl_obs;
                        map_lines[lm_idx0]->addMapLineObservation( map_keyframes[kf_curr_idx]->stereo_frame->ldesc_l.row(lm_ldx1),
                                                                   map_keyframes[kf_curr_idx]->kf_idx,
                                                                   map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->le, dir, pts );
                        // increase full graph for each KF that has already observed this LM
                        for( auto kf_it = map_lines[lm_idx0]->kf_obs_list.begin(); kf_it != map_lines[lm_idx0]->kf_obs_list.end(); kf_it++)
                        {
                            full_graph[(*kf_it)][kf_prev_idx]++;
                            full_graph[kf_prev_idx][(*kf_it)]++;
                        }
                    }
                }
                // if not, create LM and add both observations
                if( lm_idx0 == -1 && lm_idx1 == -1 )
                {
                    if( map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0] != NULL && map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1] != NULL )
                    {
                        // assign indices
                        map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->idx = max_ls_idx;
                        map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->idx = max_ls_idx;
                        // create new 3D landmark with the observation from previous KF
                        Matrix4d Tfw  = ( map_keyframes[kf_prev_idx]->T_kf_w );
                        Vector3d sP3d = Tfw.block(0,0,3,3) * map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->sP + Tfw.col(3).head(3);
                        Vector3d eP3d = Tfw.block(0,0,3,3) * map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->eP + Tfw.col(3).head(3);
                        Vector3d mP3d = 0.5 * ( sP3d + eP3d );
                        Vector3d dir = mP3d / mP3d.norm();
                        Vector6d L3D; L3D.head(3) = sP3d; L3D.tail(3) = eP3d;
                        Vector4d pts;
                        pts.head(2) = map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->spl;
                        pts.tail(2) = map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->epl;
                        MapLine* map_line = new MapLine(max_ls_idx,L3D,map_keyframes[kf_prev_idx]->stereo_frame->ldesc_l.row(lm_ldx0),map_keyframes[kf_prev_idx]->kf_idx,
                                                        map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->le,dir,pts);
                        // add new 3D landmark to kf_idx where it was first observed
                        map_lines_kf_idx.at( kf_prev_idx ).push_back( max_ls_idx );
                        // add observation of the 3D landmark from current KF
                        sP3d = map_keyframes[kf_curr_idx]->T_kf_w.block(0,0,3,3) *  map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->sP + map_keyframes[kf_curr_idx]->T_kf_w.col(3).head(3);
                        eP3d = map_keyframes[kf_curr_idx]->T_kf_w.block(0,0,3,3) *  map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->eP + map_keyframes[kf_curr_idx]->T_kf_w.col(3).head(3);
                        mP3d = 0.5 * ( sP3d + eP3d );
                        dir = mP3d / mP3d.norm();
                        pts.head(2) = map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->spl;
                        pts.tail(2) = map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx1]->epl;
                        map_line->addMapLineObservation(map_keyframes[kf_curr_idx]->stereo_frame->ldesc_l.row(lm_ldx1),map_keyframes[kf_curr_idx]->kf_idx,
                                                        map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->le,dir,pts);
                        // add 3D landmark to map
                        map_lines.push_back(map_line);
                        // update full graph (new feature)
                        max_ls_idx++;
                        full_graph[kf_prev_idx][kf_curr_idx]++;
                        full_graph[kf_curr_idx][kf_prev_idx]++;
                    }
                }
                // if the LM observed is different in each KF, then fuse them and erase the old one
                if( lm_idx0 != -1 && lm_idx1 != -1 )
                {
                    if( map_lines[lm_idx0] != NULL && map_lines[lm_idx1] != NULL
                        && map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1] != NULL )
                    {
                        int Nobs_lm_prev = map_lines[lm_idx0]->kf_obs_list.size();
                        // fuse LMs while updating the full graph
                        int iter = 0;
                        for( auto it = map_lines[lm_idx1]->desc_list.begin(); it != map_lines[lm_idx1]->desc_list.end(); it++, iter++)
                        {
                            // concatenate desc, obs, dir, pts, and kf_obs lists
                            map_lines[lm_idx0]->desc_list.push_back(   (*it) );
                            map_lines[lm_idx0]->obs_list.push_back(    map_lines[lm_idx1]->obs_list[iter]    );
                            map_lines[lm_idx0]->dir_list.push_back(    map_lines[lm_idx1]->dir_list[iter]    );
                            map_lines[lm_idx0]->pts_list.push_back(    map_lines[lm_idx1]->pts_list[iter]    );
                            map_lines[lm_idx0]->kf_obs_list.push_back( map_lines[lm_idx1]->kf_obs_list[iter] );
                            // update full graph
                            int jdx = map_lines[lm_idx1]->kf_obs_list[iter];
                            for( int i = 0; i < Nobs_lm_prev; i++ )
                            {
                                int idx = map_lines[lm_idx0]->kf_obs_list[i];
                                full_graph[idx][jdx]++;
                                full_graph[jdx][idx]++;
                            }
                            // update average descriptor and direction of observation
                            map_lines[lm_idx0]->updateAverageDescDir();
                            // change idx in stereo_ls
                            map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->idx = lm_idx0;
                        }
                        // remove from map_points_kf_idx
                        iter = 0;
                        int kf_lm_obs = map_lines[lm_idx1]->kf_obs_list[0];
                        for( auto it = map_lines_kf_idx.at(kf_lm_obs).begin(); it != map_lines_kf_idx.at(kf_lm_obs).end(); it++, iter++)
                        {
                            if( (*it) == lm_idx1 )
                            {
                                map_lines_kf_idx.at(kf_lm_obs).erase( map_lines_kf_idx.at(kf_lm_obs).begin() + iter );
                                break;
                            }
                        }
                        // erase old landmark
                        delete map_lines[lm_idx1];
                        map_lines[lm_idx1] = nullptr;
                    }
                }
            }
        }
    }

}

void MapHandler::print_msg(const std::string &msg) {

    {
        std::lock_guard<std::mutex> lock(cout_mutex);
        cout << msg << endl;
    }
}

void MapHandler::addKeyframeToSW(PLSLAM::KeyFrame *kf) {
    lk_sw.lock();
    slideWindow.push_back(kf);
    lk_sw.unlock();
}

void MapHandler::deleteKeyframeInSW() {
    lk_sw.lock();
    slideWindow.pop_front();
    lk_sw.unlock();
}

bool MapHandler::tryVioInit() {

    lk_sw.lock();
    std::deque<KeyFrame*> vpKFs = slideWindow;
    lk_sw.unlock();

    int N = vpKFs.size();
    if (N < 10)    // 初始化需要若干个关键帧
        return false;

    cout<<"[begin IMU bias_g init...."<<endl;
    Vector3d bgest = IMUInitEstBg(vpKFs);

    cout<<"[bias_g init done]: "<<bgest.transpose()<<endl;

    for(auto pkf : vpKFs){
        pkf->SetNavStateBiasGyr(bgest);
    }

    KeyFrame* kf0 = vpKFs.front();
    for(auto pkf : vpKFs){
        if(pkf==kf0)
            continue;
        pkf->ComputeIMUPreIntSinceLastFrame(pkf->RefKeyframe);
    }

    MatrixXd C(3 * (N-2), 3);
    C.setZero();

    VectorXd D(3 * (N-2));
    D.setZero();

    Matrix3d I3 = Matrix3d::Identity();
    Matrix4d Tcb = cam->getTbs().inverse();
    Vector3d pcb = Tcb.block<3,1>(0,3);
    Matrix3d Rcb = Tcb.block<3,3>(0,0);
    for(int i=0; i<N-2;i++){
        KeyFrame* pKF1 = vpKFs[i];
        KeyFrame* pKF2 = vpKFs[i+1];
        KeyFrame* pKF3 = vpKFs[i+2];
        //Delta time between frames
        double dt12 = pKF2->GetIMUPreInt().getDeltaTime();
        double dt23 = pKF3->GetIMUPreInt().getDeltaTime();
        //Pre-integrated measurements
        Vector3d dp12 = pKF2->GetIMUPreInt().getDeltaP();
        Vector3d dv12 = pKF2->GetIMUPreInt().getDeltaV();
        Vector3d dp23 = pKF3->GetIMUPreInt().getDeltaP();

        //Rotation of camera in world frame
        Matrix3d Rc1 = pKF1->T_kf_w.block<3,3>(0,0);
        Matrix3d Rc2 = pKF2->T_kf_w.block<3,3>(0,0);
        Matrix3d Rc3 = pKF3->T_kf_w.block<3,3>(0,0);
        //Position of camera center
        Vector3d pc1 = pKF1->T_kf_w.block<3,1>(0,3);
        Vector3d pc2 = pKF2->T_kf_w.block<3,1>(0,3);
        Vector3d pc3 = pKF3->T_kf_w.block<3,1>(0,3);

        //Stack to C/D matrix
        Vector3d lambda = (pc2-pc1)*dt23 + (pc2-pc3)*dt12;
        Matrix3d beta = 0.5*I3*(dt12*dt12*dt23 + dt12*dt23*dt23);
        Vector3d gamma = (Rc3-Rc2)*pcb*dt12 + (Rc1-Rc2)*pcb*dt23 + Rc1*Rcb*dp12*dt23 - Rc2*Rcb*dp23*dt12 - Rc1*Rcb*dv12*dt12*dt23;

        C.block<3,3>(3*i,0) = beta;
        D.segment<3>(3*i) = gamma - lambda;
    }

    JacobiSVD<MatrixXd> svd2(C, ComputeThinU | ComputeThinV);
    VectorXd y = svd2.solve(D);
    Vector3d gpre = y.head(3);

    cout<<"[Gravity init is done...]"<<"   G: "<<gpre.transpose()<<endl;

    //normalize g
    Vector3d g0 = gpre / gpre.norm() * 9.810;

    //estimate the bias from g
    MatrixXd A(3*(N-2),3);
    A.setZero();
    VectorXd B(3*(N-2));
    B.setZero();

    for(int i = 0;i < N-2;i++) {
        KeyFrame* pKF1 = vpKFs[i];
        KeyFrame* pKF2 = vpKFs[i+1];
        KeyFrame* pKF3 = vpKFs[i+2];

        //Delta time between frames
        double dt12 = pKF2->GetIMUPreInt().getDeltaTime();
        double dt23 = pKF3->GetIMUPreInt().getDeltaTime();
        //Pre-integrated measurements
        Vector3d dp12 = pKF2->GetIMUPreInt().getDeltaP();
        Vector3d dv12 = pKF2->GetIMUPreInt().getDeltaV();
        Vector3d dp23 = pKF3->GetIMUPreInt().getDeltaP();

        Matrix3d Jpba12 = pKF2->GetIMUPreInt().getJPBiasa();
        Matrix3d Jvba12 = pKF2->GetIMUPreInt().getJVBiasa();
        Matrix3d Jpba23 = pKF3->GetIMUPreInt().getJPBiasa();

        //Rotation of camera in world frame
        Matrix3d Rc1 = pKF1->T_kf_w.block<3,3>(0,0);
        Matrix3d Rc2 = pKF2->T_kf_w.block<3,3>(0,0);
        Matrix3d Rc3 = pKF3->T_kf_w.block<3,3>(0,0);
        //Position of camera center
        Vector3d pc1 = pKF1->T_kf_w.block<3,1>(0,3);
        Vector3d pc2 = pKF2->T_kf_w.block<3,1>(0,3);
        Vector3d pc3 = pKF3->T_kf_w.block<3,1>(0,3);

        Vector3d lambda = pc2*dt23 - pc3*dt12 -pc1*dt23 + pc2*dt12 + 0.5*g0*(dt12*dt12*dt23+dt23*dt23*dt12) -
                          Rc1*Rcb*dp12*dt23 + Rc2*Rcb*dp23*dt12 - (Rc1-Rc2)*pcb*dt23 + (Rc2-Rc3)*pcb*dt12 +
                          Rc1*Rcb*dv12*dt12*dt23;
        Matrix3d fai = Rc1*Rcb*Jpba12*dt23 - Rc2*Rcb*Jpba23*dt12 - Rc1*Rcb*Jvba12*dt12*dt23;
        A.block<3,3>(i*3,0) = fai;
        B.segment<3>(3*i) = lambda;
    }

    JacobiSVD<MatrixXd> svd(A, ComputeThinV | ComputeThinU);
    VectorXd y2 = svd.solve(B);
    Vector3d baest = y2.head(3);
    cout<<"[bias acc init done...."<<"   ba: "<<baest.transpose()<<endl;

    for(int i=0;i<N;i++){
        vpKFs[i]->SetNavStateBiasAcc(baest);
    }
    for(int i=1;i<N;i++){
        vpKFs[i]->ComputeIMUPreIntSinceLastFrame(vpKFs[i]->RefKeyframe);
    }

    //estimate the velocity
    for(int i=0;i<N;i++){
        auto pKF = vpKFs[i];
        if(i != N-1){
            auto pKF2 = vpKFs[i+1];
            const Vector3d pb2 = pKF2->imuState.Get_P();
            const Vector3d pb1 = pKF->imuState.Get_P();
            const Matrix3d R1 = pKF->imuState.Get_R().matrix();
            const double dt12 = pKF2->GetIMUPreInt().getDeltaTime();
            const Vector3d dp12 = pKF2->GetIMUPreInt().getDeltaP();

            Vector3d v1 = (pb2-pb1-0.5*g0*dt12*dt12-R1*dp12) / dt12;
            pKF->imuState.Set_Vel(v1);
            cout<<"Velocity ["<<i<<"]: "<<v1.transpose()<<endl;
        }
        else{
            auto pKF0 = vpKFs[i-1];
            const Matrix3d R0 = pKF0->imuState.Get_R().matrix();
            const Vector3d v0 = pKF0->imuState.Get_V();
            const double dt01 = pKF->GetIMUPreInt().getDeltaTime();
            const Vector3d dv01 = pKF->GetIMUPreInt().getDeltaV();

            Vector3d v1 = v0 + g0*dt01 + R0*dv01;
            pKF->SetNavStateVel(v1);
            cout<<"Velocity ["<<i<<"]: "<<v1.transpose()<<endl;
        }
    }

    gw = g0;

    vioInit = true;

    return true;
}

Vector3d MapHandler::IMUInitEstBg(const std::deque<KeyFrame*> &vpKFs) {

    // Setup optimizer
    g2o::SparseOptimizer optimizer;

    auto linearSolver = g2o::make_unique<SlamLinearSolver>();
    linearSolver->setBlockOrdering(false);
    auto blockSolver = g2o::make_unique<g2o::BlockSolverX>(std::move(linearSolver));
    g2o::OptimizationAlgorithm *algorithm = new g2o::OptimizationAlgorithmLevenberg(std::move(blockSolver));

    optimizer.setAlgorithm(algorithm);
    // Add vertex of gyro bias, to optimizer graph
    g2o::VertexGyrBias *vBiasg = new g2o::VertexGyrBias();
    vBiasg->setEstimate(Eigen::Vector3d::Zero());
    vBiasg->setId(0);
    optimizer.addVertex(vBiasg);

    // Add unary edges for gyro bias vertex
    KeyFrame* pPrevKF0 = vpKFs.front();
    for (auto pKF : vpKFs) {
        // Ignore the first KF
        if (pKF == vpKFs.front())
            continue;

        KeyFrame* pPrevKF = pKF->RefKeyframe;

        const IMUPreintegrator &imupreint = pKF->GetIMUPreInt();
        g2o::EdgeGyrBias *eBiasg = new g2o::EdgeGyrBias();
        eBiasg->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(0)));

        // measurement is not used in EdgeGyrBias
        eBiasg->dRbij = imupreint.getDeltaR();
        eBiasg->J_dR_bg = imupreint.getJRBiasg();
        eBiasg->Rwbi = pPrevKF->GetNavState().Get_RotMatrix();
        eBiasg->Rwbj = pKF->GetNavState().Get_RotMatrix();
        eBiasg->setInformation(imupreint.getCovPVPhi().bottomRightCorner(3, 3).inverse());
        optimizer.addEdge(eBiasg);
    }

    // It's actualy a linear estimator, so 1 iteration is enough.
    optimizer.initializeOptimization();
    optimizer.optimize(1);

    // update bias G
    g2o::VertexGyrBias *vBgEst = static_cast<g2o::VertexGyrBias *>(optimizer.vertex(0));

    return vBgEst->estimate();
}

void MapHandler::formLocalMap( deque<KeyFrame *> localKeyFrames) {
    // reset local KFs & LMs
    for( vector<KeyFrame*>::iterator kf_it = map_keyframes.begin(); kf_it != map_keyframes.end(); kf_it++ )
    {
        if( (*kf_it) != NULL )
            (*kf_it)->local = false;
    }
    for( vector<MapPoint*>::iterator pt_it = map_points.begin(); pt_it != map_points.end(); pt_it++ )
    {
        if( (*pt_it) != NULL )
            (*pt_it)->local = false;
    }
    for( vector<MapLine*>::iterator  ls_it = map_lines.begin(); ls_it != map_lines.end(); ls_it++ )
    {
        if( (*ls_it) != NULL )
            (*ls_it)->local = false;
    }

    //Set slidewindow keyframe's [local\ to true
    for(deque<KeyFrame*>::iterator kf_sw_it = localKeyFrames.begin(); kf_sw_it!=localKeyFrames.end(); kf_sw_it++ ){
        if( (*kf_sw_it) != NULL ) {
            (*kf_sw_it)->local = true;
            KeyFrame* kf = (*kf_sw_it);

            for( vector<PointFeature*>::iterator pt_it = kf->stereo_frame->stereo_pt.begin(); pt_it != kf->stereo_frame->stereo_pt.end(); pt_it++ )
            {
                if( (*pt_it) != NULL )
                {
                    int lm_idx = (*pt_it)->idx;
                    if( lm_idx != -1 && map_points[lm_idx] != NULL )
                        map_points[lm_idx]->local = true;
                }
            }
            for( vector<LineFeature*>::iterator ls_it = kf->stereo_frame->stereo_ls.begin(); ls_it != kf->stereo_frame->stereo_ls.end(); ls_it++ )
            {
                if( (*ls_it) != NULL )
                {
                    int lm_idx = (*ls_it)->idx;
                    if( lm_idx != -1 && map_lines[lm_idx] != NULL )
                        map_lines[lm_idx]->local = true;
                }
            }

        }
    }

}

int MapHandler::localBundleAdjustmentWithImu(deque<KeyFrame *> vpKFs, bool* mbaAbort){

    cout<<"Perform Local Bundle Adjustment With Imu...."<<endl;
    //camera param
    double fx = cam->getFx();
    double fy = cam->getFy();
    double cx = cam->getCx();
    double cy = cam->getCy();

    Matrix4d Tcb = cam->getTbs();  //from camera to imu
    Matrix4d Tbc = cam->getTbs().inverse();
    Matrix3d Rbc = Tcb.block<3,3>(0,0);
    Vector3d tbc = Tcb.block<3,1>(0,3);

    vector<MapPoint*> local_pt;
    vector<MapLine*> local_ls;
    map<int, KeyFrame*> idx_fix_kfs;
    map<int, KeyFrame*> all_idx_kfs;
    local_pt.clear();
    local_ls.clear();
    idx_fix_kfs.clear();

    vector<MapPoint*> all_map_points = GetMapPoints();
    vector<MapLine*> all_map_lines = GetMapLines();
    vector<KeyFrame*> all_map_keyframes = GetMapKeyFrames();

    for(vector<MapPoint*>::iterator pt_it = all_map_points.begin(); pt_it != all_map_points.end(); pt_it++){
        if( (*pt_it) != NULL && (*pt_it)->local == true ){
            local_pt.push_back( (*pt_it) );
        }
    }
    for(vector<MapLine*>::iterator ls_it = all_map_lines.begin(); ls_it != all_map_lines.end(); ls_it++){
        if( (*ls_it) != NULL && (*ls_it)->local == true ){
            local_ls.push_back( (*ls_it) );
        }
    }

    for(vector<MapPoint*>::iterator lpt_it = local_pt.begin(); lpt_it != local_pt.end(); lpt_it++){
        vector<int> obs = (*lpt_it)->kf_obs_list;
        for(int i=0;i<obs.size();i++){
            if(all_map_keyframes[obs[i]]->kf_idx != obs[i]){
                cout<<"obs: "<<obs[i]<<"  "<<"keyframe idx: "<<all_map_keyframes[obs[i]]->kf_idx<<endl;
                cerr<<"[Wrong index in the map_keyframes and MapPoint obs.....]"<<endl;
                exit(0);
            }
            else{
                if(all_map_keyframes[obs[i]]->local != true){
                    idx_fix_kfs.insert(make_pair(obs[i],all_map_keyframes[obs[i]]));
                    all_map_keyframes[obs[i]]->local = true;
                }
            }
        }
    }

    for(vector<MapLine*>::iterator lls_it = local_ls.begin(); lls_it != local_ls.end(); lls_it++){
        vector<int> obs = (*lls_it)->kf_obs_list;
        for(int i=0;i<obs.size();i++){
            if(all_map_keyframes[obs[i]]->kf_idx != obs[i]){
                cerr<<"[Wrong index in the map_keyframes and MapLine obs.....]"<<endl;
                exit(0);
            }
            else{
                if(all_map_keyframes[obs[i]]->local != true){
                    idx_fix_kfs.insert(make_pair(obs[i],all_map_keyframes[obs[i]]));
                    all_map_keyframes[obs[i]]->local = true;
                }
            }
        }
    }

    //Ref KeyFrame
    KeyFrame* firstKF = vpKFs.front();
    KeyFrame* RefKF = NULL;
    bool firstFix = true;
    if(firstKF->RefKeyframe != NULL){
        RefKF = firstKF->RefKeyframe;
        firstFix = false;
        map<int, KeyFrame*>::iterator iter = idx_fix_kfs.find(RefKF->kf_idx);
        if(iter==idx_fix_kfs.end())
            idx_fix_kfs.insert(make_pair(RefKF->kf_idx, RefKF));
    }

    all_idx_kfs = idx_fix_kfs;
    for(deque<KeyFrame*>::const_iterator lit=vpKFs.begin(), lend=vpKFs.end(); lit!=lend; lit++){
        all_idx_kfs.insert(make_pair((*lit)->kf_idx,(*lit)));
    }


    g2o::SparseOptimizer optimizer;

    auto linearSolver = g2o::make_unique<SlamLinearSolver>();
//    linearSolver->setBlockOrdering(false);
    auto blockSolver = g2o::make_unique<g2o::BlockSolverX>(std::move(linearSolver));
    g2o::OptimizationAlgorithm *algorithm = new g2o::OptimizationAlgorithmLevenberg(std::move(blockSolver));

    optimizer.setAlgorithm(algorithm);

    if(mbaAbort)
        optimizer.setForceStopFlag(mbaAbort);

    int maxKFid = 0;

    //slidewindow keyframe vertex
    for(deque<KeyFrame*>::const_iterator lit=vpKFs.begin(), lend=vpKFs.end(); lit!=lend; lit++)
    {
        KeyFrame* pKFi = *lit;
        int idKF = pKFi->kf_idx*2;
        // Vertex of PVR
        {
            g2o::VertexNavStatePVR * vNSPVR = new g2o::VertexNavStatePVR();
            vNSPVR->setEstimate(pKFi->GetNavState());
            vNSPVR->setId(idKF);
            vNSPVR->setFixed(false);
            if(pKFi==firstKF && firstFix){
                vNSPVR->setFixed(true);
            }
            optimizer.addVertex(vNSPVR);
        }
        // Vertex of Bias
        {
            g2o::VertexNavStateBias * vNSBias = new g2o::VertexNavStateBias();
            vNSBias->setEstimate(pKFi->GetNavState());
            vNSBias->setId(idKF+1);
            vNSBias->setFixed(false);
            if(pKFi==firstKF && firstFix){
                vNSBias->setFixed(true);
            }
            optimizer.addVertex(vNSBias);
        }
        if(idKF+1>maxKFid)
            maxKFid = idKF+1;
    }

    //fix keyframe vertex
    for( map<int,KeyFrame*>::iterator lit=idx_fix_kfs.begin(), lend=idx_fix_kfs.end(); lit!=lend; lit++)
    {
        KeyFrame* pKFi = (*lit).second;
        int idKF = pKFi->kf_idx*2;
        // For common fixed KeyFrames, only add PVR vertex
        {
            g2o::VertexNavStatePVR * vNSPVR = new g2o::VertexNavStatePVR();
            vNSPVR->setEstimate(pKFi->GetNavState());
            vNSPVR->setId(idKF);
            vNSPVR->setFixed(true);
            optimizer.addVertex(vNSPVR);
        }
        // For Local-Window-Previous KeyFrame, add Bias vertex
        if( RefKF != NULL && pKFi == RefKF)
        {
            g2o::VertexNavStateBias * vNSBias = new g2o::VertexNavStateBias();
            vNSBias->setEstimate(pKFi->GetNavState());
            vNSBias->setId(idKF+1);
            vNSBias->setFixed(true);
            optimizer.addVertex(vNSBias);
        }
        if(idKF+1>maxKFid)
            maxKFid = idKF+1;
    }

    vector<g2o::EdgeNavStatePVR*> vpEdgesNavStatePVR;
    vector<g2o::EdgeNavStateBias*> vpEdgesNavStateBias;
    const float thHuberNavStatePVR = sqrt(21.666);
    const float thHuberNavStateBias = sqrt(16.812);
    Matrix<double,6,6> InvCovBgaRW = Matrix<double,6,6>::Identity();
    InvCovBgaRW.topLeftCorner(3,3) = Matrix3d::Identity()/IMUData::getGyrBiasRW2();       // Gyroscope bias random walk, covariance INVERSE
    InvCovBgaRW.bottomRightCorner(3,3) = Matrix3d::Identity()/IMUData::getAccBiasRW2();   // Accelerometer bias random walk, covariance INVERSE

    for(deque<KeyFrame*>::const_iterator lit=vpKFs.begin(), lend=vpKFs.end(); lit!=lend; lit++)
    {
        KeyFrame* pKF1 = *lit;                // Current KF, store the IMU pre-integration between previous-current
        KeyFrame* pKF0 = pKF1->RefKeyframe;   // Previous KF

        if(pKF1 == firstKF && pKF0 == NULL)
            continue;

        // PVR edge
        {
            g2o::EdgeNavStatePVR * epvr = new g2o::EdgeNavStatePVR();
            epvr->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(2*pKF0->kf_idx)));
            epvr->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(2*pKF1->kf_idx)));
            epvr->setVertex(2, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(2*pKF0->kf_idx+1)));
            epvr->setMeasurement(pKF1->GetIMUPreInt());

            Matrix9d InvCovPVR = pKF1->GetIMUPreInt().getCovPVPhi().inverse();
            epvr->setInformation(InvCovPVR);
            epvr->SetParams(gw);

            g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
            epvr->setRobustKernel(rk);
            rk->setDelta(thHuberNavStatePVR);

            optimizer.addEdge(epvr);
            vpEdgesNavStatePVR.push_back(epvr);
        }
        // Bias edge
        {
            g2o::EdgeNavStateBias * ebias = new g2o::EdgeNavStateBias();
            ebias->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(2*pKF0->kf_idx+1)));
            ebias->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(2*pKF1->kf_idx+1)));
            ebias->setMeasurement(pKF1->GetIMUPreInt());

            ebias->setInformation(InvCovBgaRW/pKF1->GetIMUPreInt().getDeltaTime());

            g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
            ebias->setRobustKernel(rk);
            rk->setDelta(thHuberNavStateBias);

            optimizer.addEdge(ebias);
            vpEdgesNavStateBias.push_back(ebias);
        }
    }

    //MapPoint vertex
    vector<g2o::EdgeNavStatePVRPointXYZ*> vpEdgesMono;
//    vpEdgesMono.reserve(nExpectedSize);
    vector<KeyFrame*> vpEdgeKFMono;
//    vpEdgeKFMono.reserve(nExpectedSize);
    vector<MapPoint*> vpMapPointEdgeMono;
//    vpMapPointEdgeMono.reserve(nExpectedSize);
    vector<int> vpLmObsIdx;
    int maxPointId = maxKFid;
    const float thHuberMono = sqrt(5.991);
    for(vector<MapPoint*>::iterator lit=local_pt.begin(), lend=local_pt.end(); lit!=lend; lit++)
    {
        MapPoint* pMP = *lit;
        g2o::VertexLMPointXYZ* vPoint = new g2o::VertexLMPointXYZ();
        vPoint->setEstimate(pMP->point3D);
        int id = pMP->idx+maxKFid+1;
        vPoint->setId(id);
        vPoint->setFixed(false);
        vPoint->setMarginalized(true);
        optimizer.addVertex(vPoint);

        const vector<int> observations = pMP->kf_obs_list;

        // Set edges between KeyFrame and MapPoint
        for(int i=0;i<observations.size();i++)
        {
            int kf_id = observations[i];
            KeyFrame* pKF = NULL;
            map<int,KeyFrame*>::iterator kf_it = all_idx_kfs.find(kf_id);
            if(kf_it==all_idx_kfs.end()){
                cerr<<"[can't find KeyFrame in all_idx_xfs.....]"<<endl;
                exit(0);
            }
            pKF = kf_it->second;

            const Vector2d obs = pMP->obs_list[i];

            g2o::EdgeNavStatePVRPointXYZ* e = new g2o::EdgeNavStatePVRPointXYZ();

            e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
            e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(2*kf_id)));
            e->setMeasurement(obs);
            const float &invSigma2 = 1.0 / pMP->sigma_list[i];
            e->setInformation(Eigen::Matrix2d::Identity()*invSigma2);

            g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
            e->setRobustKernel(rk);
            rk->setDelta(thHuberMono);

            e->SetParams(fx,fy,cx,cy,Rbc,tbc);

            optimizer.addEdge(e);
            vpEdgesMono.push_back(e);
            vpEdgeKFMono.push_back(pKF);
            vpMapPointEdgeMono.push_back(pMP);
            vpLmObsIdx.push_back(i);
        }

        if(maxPointId <id+1);
        maxPointId = id+1;
    }

    //MapLine Vertex
    vector<g2o::EdgeNavStateLine*> vlEdgesMono;
    vector<KeyFrame*> vlEdgeKFMono;
    vector<MapLine*> vlMapLineEdgeMono;
    vector<int> vlLmObsIdx;

    const float thHuberLine = sqrt(5.991);
    for(vector<MapLine*>::iterator lit=local_ls.begin(), lend=local_ls.end(); lit!=lend; lit++)
    {
        MapLine* lML = *lit;
        g2o::VertexLine* vLine = new g2o::VertexLine();
        vLine->setEstimate(lML->line3D);
        int id = lML->idx+maxPointId+1;
        vLine->setId(id);
        vLine->setMarginalized(true);
        optimizer.addVertex(vLine);

        const vector<int> observations = lML->kf_obs_list;

        // Set edges between KeyFrame and MapPoint
        for(int i=0;i<observations.size();i++)
        {
            int kf_id = observations[i];
            KeyFrame* pKF = NULL;
            map<int,KeyFrame*>::iterator kf_it = all_idx_kfs.find(kf_id);
            if(kf_it==all_idx_kfs.end()){
                cerr<<"[can't find KeyFrame in all_idx_xfs.....]"<<endl;
                exit(0);
            }
            pKF = kf_it->second;

            const Vector3d obs = lML->obs_list[i];

            g2o::EdgeNavStateLine* e = new g2o::EdgeNavStateLine();

            e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
            e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(2*kf_id)));
            e->setMeasurement(obs);
            const float &invSigma2 = 1.0 / lML->sigma_list[i];
            e->setInformation(Eigen::Matrix3d::Identity()*invSigma2);

            g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
            e->setRobustKernel(rk);
            rk->setDelta(thHuberLine);

            e->SetParams(fx,fy,cx,cy,Rbc,tbc);

            optimizer.addEdge(e);
            vlEdgesMono.push_back(e);
            vlEdgeKFMono.push_back(pKF);
            vlMapLineEdgeMono.push_back(lML);
            vlLmObsIdx.push_back(i);
        }
    }


    //MapLinePoint Vertex
//    vector<g2o::EdgeNavStateLinePoint*> vlpEdgesMono;
//    vector<KeyFrame*> vlpEdgeKFMono;
//    vector<MapLine*> vlpMapLineEdgeMono;
//
//    const float thHuberLine = sqrt(5.991);
//    for(vector<MapLine*>::iterator lit=local_ls.begin(), lend=local_ls.end(); lit!=lend; lit++)
//    {
//        MapLine* lML = *lit;
//        g2o::VertexLinePoint* vLinePoints = new g2o::VertexLinePoint();
//        g2o::VertexLinePoint* vLinePointe = new g2o::VertexLinePoint();
//        vLinePoints->setEstimate((lML->line3D).head(3));
//        vLinePointe->setEstimate((lML->line3D).tail(3));
//        int ids = lML->idx*2+maxPointId+1;
//        int ide = lML->idx*2+1+maxPointId+1;
////        std::cout<<"ids: "<<ids<<"  ide: "<<ide<<std::endl;
//        vLinePoints->setId(ids);
//        vLinePointe->setId(ide);
//        vLinePoints->setFixed(false);
//        vLinePointe->setFixed(false);
//        vLinePoints->setMarginalized(true);
//        vLinePointe->setMarginalized(true);
//        optimizer.addVertex(vLinePoints);
//        optimizer.addVertex(vLinePointe);
//
//        const vector<int> observations = lML->kf_obs_list;
//
//       // Set edges between KeyFrame and MapLinePoint
//        for(int i=0;i<observations.size();i++)
//        {
//            int kf_id = observations[i];
//            KeyFrame* pKF = NULL;
//            map<int,KeyFrame*>::iterator kf_it = all_idx_kfs.find(kf_id);
//            if(kf_it==all_idx_kfs.end()){
//                cerr<<"[can't find KeyFrame in all_idx_xfs.....]"<<endl;
//                exit(0);
//            }
//            pKF = kf_it->second;
//
//            const Vector3d obs = lML->obs_list[i];
//
//            const float &invSigma2 = 1.0 / lML->sigma_list[i];
//            g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
//            rk->setDelta(thHuberLine);
//
//            g2o::EdgeNavStateLinePoint* es = new g2o::EdgeNavStateLinePoint();
//            g2o::EdgeNavStateLinePoint* ee = new g2o::EdgeNavStateLinePoint();
//
//            es->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(ids)));
//            es->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(2*kf_id)));
//            es->setMeasurement(obs);
//            es->setInformation(Eigen::Matrix<double,3,3>::Identity()*invSigma2);
//            es->setRobustKernel(rk);
//            es->SetParams(fx,fy,cx,cy,Rbc,tbc);
//
//            ee->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(ide)));
//            ee->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(2*kf_id)));
//            ee->setMeasurement(obs);
//            ee->setInformation(Eigen::Matrix<double,3,3>::Identity()*invSigma2);
//            ee->setRobustKernel(rk);
//            ee->SetParams(fx,fy,cx,cy,Rbc,tbc);
//
//            optimizer.addEdge(es);
//            optimizer.addEdge(ee);
//            vlpEdgesMono.push_back(es);
//            vlpEdgesMono.push_back(ee);
//            vlpEdgeKFMono.push_back(pKF);
//            vlpMapLineEdgeMono.push_back(lML);
//        }
//    }

    if(mbaAbort)
        if(*mbaAbort)
            return 1;

    //  cout<<"Begin optimize...."<<endl;
    //   optimizer.setVerbose(true);
    optimizer.initializeOptimization();
    optimizer.optimize(5);
    //  cout<<"End optimize"<<endl;

    bool bDoMore= true;

    if(mbaAbort)
        if(*mbaAbort)
            bDoMore = false;
    if(bDoMore) {
        //MapPoint error
        for (size_t i = 0, iend = vpEdgesMono.size(); i < iend; i++) {
            g2o::EdgeNavStatePVRPointXYZ *e = vpEdgesMono[i];
            //MapPoint* pMP = vpMapPointEdgeMono[i];

            if (e->chi2() > 5.991 || !e->isDepthPositive()) {
                e->setLevel(1);
            }

            e->setRobustKernel(0);
        }

        //MapLine error
        for (size_t i = 0; i < vlEdgesMono.size(); i++) {
            g2o::EdgeNavStateLine *e = vlEdgesMono[i];
            MapLine *lML = vlMapLineEdgeMono[i];

            if (e->chi2() > 5.991 || !e->isDepthPositive()) {
                e->setLevel(1);
            }
            e->setRobustKernel(0);
        }

        //MapLinePoint error
//    for(size_t i=0; i<vlpEdgesMono.size();i++) {
//        g2o::EdgeNavStateLinePoint *e = vlpEdgesMono[i];
//
//        if (e->chi2() > 5.991 || !e->isDepthPositive()) {
//            e->setLevel(1);
//        }
//        e->setRobustKernel(0);
//    }

        //   optimizer.setVerbose(true);
        optimizer.initializeOptimization(0);
        optimizer.optimize(10);
    }

    //remove bad point observations
    //std::cout<<"Begin romove bad point obs"<<std::endl;
    int bad_point_obs = 0;
    int actually_bad_point_obs = 0;
    for(int i = vpEdgesMono.size()-1; i>=0; i--){
        g2o::EdgeNavStatePVRPointXYZ* e = vpEdgesMono[i];
        if(e->level()==1){
            e->computeError();
        }
        if(e->chi2() > 5.991 || !e->isDepthPositive()){
            bad_point_obs++;
            KeyFrame* kf = vpEdgeKFMono[i];
            MapPoint* pMP = vpMapPointEdgeMono[i];
            if( pMP->obs_list.size() > 1 )
            {
                actually_bad_point_obs++;
                int kf_obs = kf->kf_idx;
                int lm_idx_map = pMP->idx;
                int lm_idx_obs = vpLmObsIdx[i];
                // if it is the first observation, update it from map_points_kf_idx
                if( lm_idx_obs == 0 )
                {
                    // delete observation from map_points_kf_idx
                    for( auto it = map_points_kf_idx.at(kf_obs).begin(); it != map_points_kf_idx.at(kf_obs).end(); it++)
                    {
                        if( (*it) == lm_idx_map )
                        {
                            int new_kf_base = pMP->kf_obs_list[1];
                            map_points_kf_idx.at(new_kf_base).push_back( (*it) );
                            break;
                        }
                    }
                }
                // remove observations from map points
                pMP->desc_list.erase( pMP->desc_list.begin() + lm_idx_obs );
                pMP->obs_list.erase( pMP->obs_list.begin() + lm_idx_obs );
                pMP->dir_list.erase( pMP->dir_list.begin() + lm_idx_obs );
                pMP->kf_obs_list.erase( pMP->kf_obs_list.begin() + lm_idx_obs );
                // remove idx from KeyFrame stereo points
                for(vector<PointFeature*>::iterator st_pt = kf->stereo_frame->stereo_pt.begin();
                    st_pt != kf->stereo_frame->stereo_pt.end(); st_pt++ )
                {
                    if( (*st_pt)->idx == lm_idx_map )
                    {
                        (*st_pt)->idx = -1;
                        st_pt = kf->stereo_frame->stereo_pt.end()-1;
                    }
                }
                // update main descriptor and direction
                pMP->updateAverageDescDir();
                // update graphs
                for( int i = 0; i < pMP->kf_obs_list.size(); i++ )
                {
                    int idx = pMP->kf_obs_list[i];
                    if( kf_obs != idx )
                    {
                        full_graph[kf_obs][idx]--;
                        full_graph[idx][kf_obs]--;
                    }
                }
            }
            else{
                pMP->inlier = false;
            }
        }
    }
    std::cout<<"All Point Edge is "<<vpEdgesMono.size()<<", Bad Point Edge is "<<bad_point_obs<<", Actually delete "
             <<actually_bad_point_obs<<" Edges."<<std::endl;

    //remove bad Line observations
    int bad_line_obs = 0;
    int actually_bad_line_obs = 0;
    for(int i =vlEdgesMono.size()-1; i>=0; i--){
        g2o::EdgeNavStateLine* e = vlEdgesMono[i];
        if(e->level()==1){
            e->computeError();
        }

        if(e->chi2() > 5.991 || !e->isDepthPositive()){
            bad_line_obs++;

            KeyFrame* kf = vlEdgeKFMono[i];
            MapLine* lML = vlMapLineEdgeMono[i];

            if( lML->obs_list.size() > 1 )
            {
                actually_bad_line_obs++;
                int kf_obs = kf->kf_idx;
                int lm_idx_map = lML->idx;
                int lm_idx_obs = vlLmObsIdx[i];
                // if it is the first observation, update it from map_points_kf_idx
                if( lm_idx_obs == 0 )
                {
                    // delete observation from map_points_kf_idx
                    for( auto it = map_points_kf_idx.at(kf_obs).begin(); it != map_points_kf_idx.at(kf_obs).end(); it++)
                    {
                        if( (*it) == lm_idx_map )
                        {
                            int new_kf_base = lML->kf_obs_list[1];
                            map_points_kf_idx.at(new_kf_base).push_back( (*it) );
                            break;
                        }
                    }
                }

                // remove observations from map points
                lML->desc_list.erase( lML->desc_list.begin() + lm_idx_obs );
                lML->obs_list.erase( lML->obs_list.begin() + lm_idx_obs );
                lML->dir_list.erase( lML->dir_list.begin() + lm_idx_obs );
                lML->kf_obs_list.erase( lML->kf_obs_list.begin() + lm_idx_obs );
                lML->pts_list.erase( lML->pts_list.begin() + lm_idx_obs );

                // remove idx from KeyFrame stereo points
                for(vector<LineFeature*>::iterator st_pt = kf->stereo_frame->stereo_ls.begin();
                    st_pt != kf->stereo_frame->stereo_ls.end(); st_pt++ )
                {
                    if( (*st_pt)->idx == lm_idx_map )
                    {
                        (*st_pt)->idx = -1;
                        st_pt = kf->stereo_frame->stereo_ls.end()-1;
                    }
                }
                // update main descriptor and direction
                lML->updateAverageDescDir();
                // update graphs
                for( int i = 0; i < lML->kf_obs_list.size(); i++ )
                {
                    int idx = lML->kf_obs_list[i];
                    if( kf_obs != idx )
                    {
                        full_graph[kf_obs][idx]--;
                        full_graph[idx][kf_obs]--;
                    }
                }
            }
            else{
                lML->inlier = false;
            }

        }
    }
    std::cout<<"All Point Edge is "<<vpEdgesMono.size()<<", Bad Point Edge is "<<bad_line_obs<<", Actually delete "
             <<actually_bad_line_obs<<" Edges."<<std::endl;


    // recover Keyframes
    for(deque<KeyFrame*>::const_iterator lit=vpKFs.begin(), lend=vpKFs.end(); lit!=lend; lit++)
    {
        KeyFrame* pKFi = *lit;
        g2o::VertexNavStatePVR* vNSPVR = static_cast<g2o::VertexNavStatePVR*>(optimizer.vertex(2*pKFi->kf_idx));
        g2o::VertexNavStateBias* vNSBias = static_cast<g2o::VertexNavStateBias*>(optimizer.vertex(2*pKFi->kf_idx+1));
        // In optimized navstate, bias not changed, delta_bias not zero, should be added to bias
        const NavState& optPVRns = vNSPVR->estimate();
        const NavState& optBiasns = vNSBias->estimate();
        NavState primaryns = pKFi->GetNavState();
        // Update NavState
        pKFi->SetNavStatePos(optPVRns.Get_P());
        pKFi->SetNavStateVel(optPVRns.Get_V());
        pKFi->SetNavStateRot(optPVRns.Get_R());

        pKFi->SetNavStateDeltaBg(optBiasns.Get_dBias_Gyr());
        pKFi->SetNavStateDeltaBa(optBiasns.Get_dBias_Acc());

        //update Pose
        Matrix3d Rot = pKFi->GetNavState().Get_RotMatrix() * Tcb.block<3,3>(0,0);
        Vector3d trans = pKFi->GetNavState().Get_P() + pKFi->GetNavState().Get_RotMatrix()*Tcb.block<3,1>(0,3);
        pKFi->SetCameraRot(Rot);
        pKFi->SetCameraTans(trans);
    }

    //recover MapPoints
    for(vector<MapPoint*>::const_iterator lit=local_pt.begin(), lend=local_pt.end(); lit!=lend; lit++)
    {
        MapPoint* pMP = *lit;
        g2o::VertexLMPointXYZ* vPoint = static_cast<g2o::VertexLMPointXYZ*>(optimizer.vertex(pMP->idx+maxKFid+1));
        pMP->point3D = vPoint->estimate();
    }

    //recover MapLines
    for(vector<MapLine*>::const_iterator lit=local_ls.begin(), lend=local_ls.end(); lit!=lend; lit++){
        MapLine* lML = *lit;
        g2o::VertexLine* vLine = static_cast<g2o::VertexLine*>(optimizer.vertex(lML->idx+maxPointId+1));
        lML->line3D = vLine->estimate();
    }

    //recover MapLinePoints
//    for(vector<MapLine*>::const_iterator lit=local_ls.begin(), lend=local_ls.end(); lit!=lend; lit++){
//        MapLine* lML = *lit;
//        g2o::VertexLinePoint* vLinePoints = static_cast<g2o::VertexLinePoint*>(optimizer.vertex(lML->idx*2+maxPointId+1));
//        g2o::VertexLinePoint* vLinePointe = static_cast<g2o::VertexLinePoint*>(optimizer.vertex(lML->idx*2+1+maxPointId+1));
//        Vector6d est;
//        est.head(3) = vLinePoints->estimate();
//        est.tail(3) = vLinePointe->estimate();
//        lML->line3D = est;
//    }

    cout<<"Finish Local Bundle Adjustment With IMU!"<<endl;
    return 0;
}

int MapHandler::localBundleAdjustmentWithImuAndMarg(deque<KeyFrame *> vpKFs, bool *mbaAbort) {

    cout<<"Perform Local Bundle Adjustment With Imu And Marg...."<<endl;
    //camera param
    double fx = cam->getFx();
    double fy = cam->getFy();
    double cx = cam->getCx();
    double cy = cam->getCy();

    Matrix4d Tcb = cam->getTbs();  //from camera to imu
    Matrix4d Tbc = cam->getTbs().inverse();
    Matrix3d Rbc = Tcb.block<3,3>(0,0);
    Vector3d tbc = Tcb.block<3,1>(0,3);

    vector<MapPoint*> local_pt;
    vector<MapLine*> local_ls;
    map<int, KeyFrame*> idx_fix_kfs;
    map<int, KeyFrame*> all_idx_kfs;
    local_pt.clear();
    local_ls.clear();
    idx_fix_kfs.clear();

    vector<MapPoint*> all_map_points = GetMapPoints();
    vector<MapLine*> all_map_lines = GetMapLines();
    vector<KeyFrame*> all_map_keyframes = GetMapKeyFrames();

    //form local map
    int first_kf_idx = vpKFs.front()->kf_idx;
    for(int i=0;i<all_map_points.size();i++){
        MapPoint* pt = all_map_points[i];
        if(pt!=NULL && pt->kf_obs_list.size()>0 && pt->kf_obs_list[0]>=first_kf_idx){
            local_pt.push_back(pt);
        }
    }

    for(int i=0;i<all_map_lines.size();i++){
        MapLine* ls = all_map_lines[i];
        if(ls!=NULL && ls->kf_obs_list.size()>0 && ls->kf_obs_list[0]>=first_kf_idx){
            local_ls.push_back(ls);
        }
    }

    for(int i=0;i<vpKFs.size();i++){
        all_idx_kfs.insert(make_pair(vpKFs[i]->kf_idx,vpKFs[i]));
    }

    g2o::SparseOptimizer optimizer;

    auto linearSolver = g2o::make_unique<SlamLinearSolver>();
//    linearSolver->setBlockOrdering(false);
    auto blockSolver = g2o::make_unique<g2o::BlockSolverX>(std::move(linearSolver));
    g2o::OptimizationAlgorithm *algorithm = new g2o::OptimizationAlgorithmLevenberg(std::move(blockSolver));

    optimizer.setAlgorithm(algorithm);

    if(mbaAbort)
        optimizer.setForceStopFlag(mbaAbort);

    int maxKFid = 0;
    bool firstFix = true;
    //slidewindow keyframe vertex
    for(deque<KeyFrame*>::const_iterator lit=vpKFs.begin(), lend=vpKFs.end(); lit!=lend; lit++)
    {
        KeyFrame* pKFi = *lit;
        int idKF = pKFi->kf_idx*2;
        // Vertex of PVR
        {
            g2o::VertexNavStatePVR * vNSPVR = new g2o::VertexNavStatePVR();
            vNSPVR->setEstimate(pKFi->GetNavState());
            vNSPVR->setId(idKF);
            vNSPVR->setFixed(false);
            if(pKFi==vpKFs[0] && firstFix){
                vNSPVR->setFixed(true);
            }
            optimizer.addVertex(vNSPVR);
        }
        // Vertex of Bias
        {
            g2o::VertexNavStateBias * vNSBias = new g2o::VertexNavStateBias();
            vNSBias->setEstimate(pKFi->GetNavState());
            vNSBias->setId(idKF+1);
            vNSBias->setFixed(false);
            if(pKFi==vpKFs[0] && firstFix){
                vNSBias->setFixed(true);
            }
            optimizer.addVertex(vNSBias);
        }
        if(idKF+1>maxKFid)
            maxKFid = idKF+1;
    }

    vector<g2o::EdgeNavStatePVR*> vpEdgesNavStatePVR;
    vector<g2o::EdgeNavStateBias*> vpEdgesNavStateBias;
    vpEdgesNavStatePVR.clear();
    vpEdgesNavStateBias.clear();
    const float thHuberNavStatePVR = sqrt(21.666);
    const float thHuberNavStateBias = sqrt(16.812);
    Matrix<double,6,6> InvCovBgaRW = Matrix<double,6,6>::Identity();
    InvCovBgaRW.topLeftCorner(3,3) = Matrix3d::Identity()/IMUData::getGyrBiasRW2();       // Gyroscope bias random walk, covariance INVERSE
    InvCovBgaRW.bottomRightCorner(3,3) = Matrix3d::Identity()/IMUData::getAccBiasRW2();   // Accelerometer bias random walk, covariance INVERSE

    for(deque<KeyFrame*>::const_iterator lit=vpKFs.begin(), lend=vpKFs.end(); lit!=lend; lit++)
    {
        KeyFrame* pKF1 = *lit;                // Current KF, store the IMU pre-integration between previous-current
        KeyFrame* pKF0 = pKF1->RefKeyframe;   // Previous KF

        if(pKF1 == vpKFs[0])
            continue;

        // PVR edge
        {
            g2o::EdgeNavStatePVR * epvr = new g2o::EdgeNavStatePVR();
            epvr->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(2*pKF0->kf_idx)));
            epvr->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(2*pKF1->kf_idx)));
            epvr->setVertex(2, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(2*pKF0->kf_idx+1)));
            epvr->setMeasurement(pKF1->GetIMUPreInt());

            Matrix9d InvCovPVR = pKF1->GetIMUPreInt().getCovPVPhi().inverse();
            epvr->setInformation(InvCovPVR);
            epvr->SetParams(gw);

            g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
            epvr->setRobustKernel(rk);
            rk->setDelta(thHuberNavStatePVR);

            optimizer.addEdge(epvr);
            vpEdgesNavStatePVR.push_back(epvr);
        }
        // Bias edge
        {
            g2o::EdgeNavStateBias * ebias = new g2o::EdgeNavStateBias();
            ebias->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(2*pKF0->kf_idx+1)));
            ebias->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(2*pKF1->kf_idx+1)));
            ebias->setMeasurement(pKF1->GetIMUPreInt());

            ebias->setInformation(InvCovBgaRW/pKF1->GetIMUPreInt().getDeltaTime());

            g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
            ebias->setRobustKernel(rk);
            rk->setDelta(thHuberNavStateBias);

            optimizer.addEdge(ebias);
            vpEdgesNavStateBias.push_back(ebias);
        }
    }

    //MapPoint vertex
    vector<g2o::EdgeNavStatePVRPointXYZ*> vpEdgesMono;
//    vpEdgesMono.reserve(nExpectedSize);
    vector<KeyFrame*> vpEdgeKFMono;
//    vpEdgeKFMono.reserve(nExpectedSize);
    vector<MapPoint*> vpMapPointEdgeMono;
//    vpMapPointEdgeMono.reserve(nExpectedSize);
    vector<int> vpLmObsIdx;
    int maxPointId = maxKFid;
    const float thHuberMono = sqrt(5.991);
    for(vector<MapPoint*>::iterator lit=local_pt.begin(), lend=local_pt.end(); lit!=lend; lit++)
    {
        MapPoint* pMP = *lit;
        g2o::VertexLMPointXYZ* vPoint = new g2o::VertexLMPointXYZ();
        vPoint->setEstimate(pMP->point3D);
        int id = pMP->idx+maxKFid+1;
        vPoint->setId(id);
        vPoint->setFixed(false);
        vPoint->setMarginalized(true);
        optimizer.addVertex(vPoint);

        const vector<int> observations = pMP->kf_obs_list;

        // Set edges between KeyFrame and MapPoint
        for(int i=0;i<observations.size();i++)
        {
            int kf_id = observations[i];
            KeyFrame* pKF = NULL;
            map<int,KeyFrame*>::iterator kf_it = all_idx_kfs.find(kf_id);
            if(kf_it==all_idx_kfs.end()){
                cout<<"[can't find KeyFrame in all_idx_xfs.....]"<<endl;
//                exit(0);
                continue;
            }
            pKF = kf_it->second;

            const Vector2d obs = pMP->obs_list[i];

            g2o::EdgeNavStatePVRPointXYZ* e = new g2o::EdgeNavStatePVRPointXYZ();

            e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
            e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(2*kf_id)));
            e->setMeasurement(obs);
            const float &invSigma2 = 1.0 / pMP->sigma_list[i];
            e->setInformation(Eigen::Matrix2d::Identity()*invSigma2);

            g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
            e->setRobustKernel(rk);
            rk->setDelta(thHuberMono);

            e->SetParams(fx,fy,cx,cy,Rbc,tbc);

            optimizer.addEdge(e);
            vpEdgesMono.push_back(e);
            vpEdgeKFMono.push_back(pKF);
            vpMapPointEdgeMono.push_back(pMP);
            vpLmObsIdx.push_back(i);
        }

        if(maxPointId <id+1);
        maxPointId = id+1;
    }

    //MapLine Vertex
    vector<g2o::EdgeNavStateLine*> vlEdgesMono;
    vector<KeyFrame*> vlEdgeKFMono;
    vector<MapLine*> vlMapLineEdgeMono;
    vector<int> vlLmObsIdx;

    const float thHuberLine = sqrt(5.991);
    for(vector<MapLine*>::iterator lit=local_ls.begin(), lend=local_ls.end(); lit!=lend; lit++)
    {
        MapLine* lML = *lit;
        g2o::VertexLine* vLine = new g2o::VertexLine();
        vLine->setEstimate(lML->line3D);
        int id = lML->idx+maxPointId+1;
        vLine->setId(id);
        vLine->setMarginalized(true);
        optimizer.addVertex(vLine);

        const vector<int> observations = lML->kf_obs_list;

        // Set edges between KeyFrame and MapPoint
        for(int i=0;i<observations.size();i++)
        {
            int kf_id = observations[i];
            KeyFrame* pKF = NULL;
            map<int,KeyFrame*>::iterator kf_it = all_idx_kfs.find(kf_id);
            if(kf_it==all_idx_kfs.end()){
                cout<<"[can't find KeyFrame in all_idx_xfs.....]"<<endl;
//                exit(0);
                continue;
            }
            pKF = kf_it->second;

            const Vector3d obs = lML->obs_list[i];

            g2o::EdgeNavStateLine* e = new g2o::EdgeNavStateLine();

            e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
            e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(2*kf_id)));
            e->setMeasurement(obs);
            const float &invSigma2 = 1.0 / lML->sigma_list[i];
            e->setInformation(Eigen::Matrix3d::Identity()*invSigma2);

            g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
            e->setRobustKernel(rk);
            rk->setDelta(thHuberLine);

            e->SetParams(fx,fy,cx,cy,Rbc,tbc);

            optimizer.addEdge(e);
            vlEdgesMono.push_back(e);
            vlEdgeKFMono.push_back(pKF);
            vlMapLineEdgeMono.push_back(lML);
            vlLmObsIdx.push_back(i);
        }
    }
    vector<g2o::EdgeMarginalization*> vpEdgesMarg;
    vpEdgesMarg.clear();
    if( marg_info ){
        std::cout<<"Begin add Marg edge..."<<std::endl;
        vector<int> keep_vertex_id = marg_info->keep_vertex_id;
        int dimension = marg_info->n;
        if(keep_vertex_id.size()>0){
            int N = keep_vertex_id.size();

            //debug
            std::cout<<"Keep vertex id: ";
            for(int i=0;i<keep_vertex_id.size();i++){
                std::cout<<keep_vertex_id[i]<<" ";
            }
            std::cout<<std::endl;
            std::cout<<"Dimension: "<<dimension<<std::endl;

            g2o::EdgeMarginalization* edgeMarg = new g2o::EdgeMarginalization();
            edgeMarg->setSize(N);
            edgeMarg->setDimension(dimension);
            for(int i=0;i<keep_vertex_id.size();i++){
                int VertexId = keep_vertex_id[i];
                edgeMarg->setVertex(i, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(VertexId)));
            }
            edgeMarg->setMeasurement(*marg_info);
            edgeMarg->setInformation(MatrixXd::Identity(dimension,dimension));
            optimizer.addEdge(edgeMarg);
            vpEdgesMarg.push_back(edgeMarg);
        }
    }

   // cout<<"Begin optimize...."<<endl;
   // optimizer.setVerbose(true);
    optimizer.initializeOptimization();
    optimizer.optimize(5);
   // cout<<"End optimize"<<endl;

    bool bDoMore= true;

    if(mbaAbort)
        if(*mbaAbort)
            bDoMore = false;
    if(bDoMore) {
        //MapPoint error
        for (size_t i = 0, iend = vpEdgesMono.size(); i < iend; i++) {
            g2o::EdgeNavStatePVRPointXYZ *e = vpEdgesMono[i];
            if (e->chi2() > 5.991 || !e->isDepthPositive()) {
                e->setLevel(1);
            }

            e->setRobustKernel(0);
        }

        //MapLine error
        for (size_t i = 0; i < vlEdgesMono.size(); i++) {
            g2o::EdgeNavStateLine *e = vlEdgesMono[i];
            if (e->chi2() > 5.991 || !e->isDepthPositive()) {
                e->setLevel(1);
            }
            e->setRobustKernel(0);
        }

        //   optimizer.setVerbose(true);
        optimizer.initializeOptimization(0);
        optimizer.optimize(10);
    }

    //do marginalization, marginalize the oldest frame
    int NUM = 50; //the number of Point and Line to marginalize, reduce the time of marginalize
    MarginalizationInfo* new_marg_info = new MarginalizationInfo();
    if(vpKFs.size()>=max_keyframe_in_slideWindow){
        std::cout<<"Begin Build Marg...."<<std::endl;
        {
            g2o::OptimizableGraph::Edge *edgePVR = dynamic_cast<g2o::OptimizableGraph::Edge *>(vpEdgesNavStatePVR[0]);
            vector<int> drop_set;
            vector<double*> addr;
            vector<VectorXd> estimate_data;
            drop_set.clear();
            addr.clear();
            estimate_data.clear();
            drop_set.push_back(0);
            vpEdgesNavStatePVR[0]->GetJacAddr(addr);
            vpEdgesNavStatePVR[0]->GetEstData(estimate_data);

            ResidualBlockInfo* ImuRes = new ResidualBlockInfo(edgePVR,drop_set,estimate_data,addr,"PVR");
            new_marg_info->addResidualBlockInfo(ImuRes);
        }

        {
            g2o::OptimizableGraph::Edge *edgeBias = dynamic_cast<g2o::OptimizableGraph::Edge*>(vpEdgesNavStateBias[0]);
            vector<int> drop_set;
            vector<double*> addr;
            vector<VectorXd> estimate_data;
            drop_set.clear();
            addr.clear();
            estimate_data.clear();
            drop_set.push_back(0);
            vpEdgesNavStateBias[0]->GetJacAddr(addr);
            vpEdgesNavStateBias[0]->GetEstData(estimate_data);

            ResidualBlockInfo* BiasRes = new ResidualBlockInfo(edgeBias,drop_set,estimate_data,addr,"BIAS");
            new_marg_info->addResidualBlockInfo(BiasRes);
        }

        {
            int num = 0;
            for(int i=0;i<vpEdgesMono.size();i++){
                MapPoint* vp = vpMapPointEdgeMono[i];
                if(vp->kf_obs_list[0] == first_kf_idx){
                    g2o::OptimizableGraph::Edge *edgePoint = dynamic_cast<g2o::OptimizableGraph::Edge*>(vpEdgesMono[i]);
                    vector<int> drop_set;
                    vector<double*> addr;
                    vector<VectorXd> estimate_data;
                    drop_set.clear();
                    addr.clear();
                    estimate_data.clear();
                    drop_set.push_back(0);
                    if(edgePoint->vertex(1)->id()==2*first_kf_idx){
                        drop_set.push_back(1);
                    }
                    vpEdgesMono[i]->GetJacAddr(addr);
                    vpEdgesMono[i]->GetEstData(estimate_data);

                    ResidualBlockInfo* PointRes = new ResidualBlockInfo(edgePoint,drop_set,estimate_data,addr,"POINT");
                    new_marg_info->addResidualBlockInfo(PointRes);
                    num++;
                    if(num>NUM){
                        break;
                    }
                }
            }
        }

        {
            int num = 0;
            for(int i=0;i<vlEdgesMono.size();i++){
                MapLine* ls = vlMapLineEdgeMono[i];
                if(ls->kf_obs_list[0] == first_kf_idx){
                    g2o::OptimizableGraph::Edge *edgeLine = dynamic_cast<g2o::OptimizableGraph::Edge*>(vlEdgesMono[i]);
                    vector<int> drop_set;
                    vector<double*> addr;
                    vector<VectorXd> estimate_data;
                    drop_set.clear();
                    addr.clear();
                    estimate_data.clear();
                    drop_set.push_back(0);
                    if(edgeLine->vertex(1)->id()==2*first_kf_idx){
                        drop_set.push_back(1);
                    }
                    vlEdgesMono[i]->GetJacAddr(addr);
                    vlEdgesMono[i]->GetEstData(estimate_data);

                    ResidualBlockInfo* LineRes = new ResidualBlockInfo(edgeLine,drop_set,estimate_data,addr,"LINE");
                    new_marg_info->addResidualBlockInfo(LineRes);
                    num++;
                    if(num>NUM){
                        break;
                    }
                }
            }
        }

        {
            if(vpEdgesMarg.size()!=0){
                g2o::OptimizableGraph::Edge* edgeMarg = dynamic_cast<g2o::OptimizableGraph::Edge*>(vpEdgesMarg[0]);
                vector<int> drop_set;
                vector<double*> addr;
                vector<VectorXd> estimate_data;
                drop_set.clear();
                addr.clear();
                estimate_data.clear();
                for(int i=0;i<vpEdgesMarg[0]->vertices().size();i++){
                    g2o::OptimizableGraph::Vertex* vertex = dynamic_cast<g2o::OptimizableGraph::Vertex*>(edgeMarg->vertex(i));
                    if(vertex->id() == 2*first_kf_idx || vertex->id()==2*first_kf_idx+1){
                        drop_set.push_back(i);
                    }
                }
                vpEdgesMarg[0]->GetJacAddr(addr);
                vpEdgesMarg[0]->GetEstData(estimate_data);

                ResidualBlockInfo* margRes = new ResidualBlockInfo(edgeMarg,drop_set,estimate_data,addr,"MARG");
                new_marg_info->addResidualBlockInfo(margRes);
            }
        }
        cout<<"begin preMarginlize...."<<endl;
        new_marg_info->preMarginalize();
        cout<<"begin marginalize...."<<endl;
        new_marg_info->marginalizeWithoutThread();

        if(marg_info!= nullptr)
            delete marg_info;

        marg_info = new_marg_info;
        std::cout<<"Marginalization, Keep vertex size: "<<marg_info->keep_vertex_size.size()<<std::endl;
    }

    // recover Keyframes
    for(deque<KeyFrame*>::const_iterator lit=vpKFs.begin(), lend=vpKFs.end(); lit!=lend; lit++)
    {
        KeyFrame* pKFi = *lit;
        g2o::VertexNavStatePVR* vNSPVR = static_cast<g2o::VertexNavStatePVR*>(optimizer.vertex(2*pKFi->kf_idx));
        g2o::VertexNavStateBias* vNSBias = static_cast<g2o::VertexNavStateBias*>(optimizer.vertex(2*pKFi->kf_idx+1));
        // In optimized navstate, bias not changed, delta_bias not zero, should be added to bias
        const NavState& optPVRns = vNSPVR->estimate();
        const NavState& optBiasns = vNSBias->estimate();
        NavState primaryns = pKFi->GetNavState();
        // Update NavState
        pKFi->SetNavStatePos(optPVRns.Get_P());
        pKFi->SetNavStateVel(optPVRns.Get_V());
        pKFi->SetNavStateRot(optPVRns.Get_R());

        pKFi->SetNavStateDeltaBg(optBiasns.Get_dBias_Gyr());
        pKFi->SetNavStateDeltaBa(optBiasns.Get_dBias_Acc());

        //update Pose
        Matrix3d Rot = pKFi->GetNavState().Get_RotMatrix() * Tcb.block<3,3>(0,0);
        Vector3d trans = pKFi->GetNavState().Get_P() + pKFi->GetNavState().Get_RotMatrix()*Tcb.block<3,1>(0,3);
        pKFi->SetCameraRot(Rot);
        pKFi->SetCameraTans(trans);
    }

    //recover MapPoints
    for(vector<MapPoint*>::const_iterator lit=local_pt.begin(), lend=local_pt.end(); lit!=lend; lit++)
    {
        MapPoint* pMP = *lit;
        g2o::VertexLMPointXYZ* vPoint = static_cast<g2o::VertexLMPointXYZ*>(optimizer.vertex(pMP->idx+maxKFid+1));
        pMP->point3D = vPoint->estimate();
    }

    //recover MapLines
    for(vector<MapLine*>::const_iterator lit=local_ls.begin(), lend=local_ls.end(); lit!=lend; lit++){
        MapLine* lML = *lit;
        g2o::VertexLine* vLine = static_cast<g2o::VertexLine*>(optimizer.vertex(lML->idx+maxPointId+1));
        lML->line3D = vLine->estimate();
    }

    //recover MapLinePoints
//    for(vector<MapLine*>::const_iterator lit=local_ls.begin(), lend=local_ls.end(); lit!=lend; lit++){
//        MapLine* lML = *lit;
//        g2o::VertexLinePoint* vLinePoints = static_cast<g2o::VertexLinePoint*>(optimizer.vertex(lML->idx*2+maxPointId+1));
//        g2o::VertexLinePoint* vLinePointe = static_cast<g2o::VertexLinePoint*>(optimizer.vertex(lML->idx*2+1+maxPointId+1));
//        Vector6d est;
//        est.head(3) = vLinePoints->estimate();
//        est.tail(3) = vLinePointe->estimate();
//        lML->line3D = est;
//    }

    cout<<"Finish Local Bundle Adjustment With IMU And Marg!"<<endl;
    return 0;
}

void MapHandler::SaveKeyFrameTrajectoryTUM(const string &filename)
{
    cout << endl << "Saving keyframe trajectory to " << filename << " ..." << endl;

    vector<KeyFrame*> vpKFs = GetMapKeyFrames();
    sort(vpKFs.begin(),vpKFs.end(),MapHandler::lId);

    ofstream f;
    f.open(filename.c_str());
    f << fixed;

    for(size_t i=0; i<vpKFs.size(); i++)
    {
        KeyFrame* pKF = vpKFs[i];

        Matrix3d R = pKF->GetCameraRot();
        Eigen::Quaterniond q(R);
        std::vector<float> v(4);
        v[0] = q.x();
        v[1] = q.y();
        v[2] = q.z();
        v[3] = q.w();
        Vector3f t = pKF->GetCameraTrans().cast<float>();
        long double tt = pKF->t_;
        f << setprecision(6) << tt<< setprecision(7) << " " << t(0) << " " << t(1) << " " << t(2)
         << " " << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << endl;

    }

    f.close();
    cout << endl << "trajectory saved!" << endl;
}



} //end namespace


