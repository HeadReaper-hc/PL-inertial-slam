project( PL-SLAM )

cmake_minimum_required(VERSION 2.7)

#gdb
#SET(CMAKE_BUILD_TYPE "Debug")
#SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g2 -ggdb")
#SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")
set(CMAKE_BUILD_TYPE "Release")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_modules")

find_package(OpenCV 3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS regex thread system filesystem)
find_package(G2O     REQUIRED)
#find_package(Sophus REQUIRED)
#find_package(Eigen3   REQUIRED)
INCLUDE_DIRECTORIES("/usr/local/include/eigen3")
if(COMMAND cmake_policy)
  	cmake_policy(SET CMP0003 NEW)	
endif(COMMAND cmake_policy)
link_directories(${OpenCV_LIBS_DIR})
include_directories(${OpenCV2_INCLUDE_DIRS})
#include_directories(${Sophus_INCLUDE_DIRS})

# Odometry library (Future TODO: integrate with this library)
SET( StVO_LIBRARY 		"" CACHE FILEPATH "Visual Odometry Library" )
SET( StVO_INCLUDE_DIR 	"" CACHE PATH "Visual Odometry Include folder" )	

set(DEFAULT_HAS_MRPT ON)
set(HAS_MRPT ${DEFAULT_HAS_MRPT} CACHE BOOL "Build the PointGrey Bumblebee2 SVO application that employs the MRPT library")

# choose whether use marginalization
set(DEFAULT_USE_MARG ON)
set(USE_MARG ${DEFAULT_USE_MARG} CACHE BOOL "Use marginalization in the localBundleAdjustment")
if(USE_MARG)
add_definitions(-DUSE_MARG)
endif(USE_MARG)

SET(BUILD_SHARED_LIBS ON)
SET(CMAKE_MODULE_PATH $ENV{CMAKE_MODULE_PATH})
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -mtune=native")

add_definitions(-DBOOST_NO_CXX11_SCOPED_ENUMS)

# MRPT library (optional, only with representation purposes)
if(HAS_MRPT)
FIND_PACKAGE(MRPT REQUIRED base opengl gui hwdrivers)
set(MRPT_DONT_USE_DBG_LIBS 1) #use release libraries for linking even if "Debug" CMake build
add_definitions(-DHAS_MRPT)
endif(HAS_MRPT)


# YAML library
FIND_PACKAGE(yaml-cpp REQUIRED CONFIG PATHS ${YAML_PATHS})

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

SET(G2O_LIBS g2o_cli g2o_ext_freeglut_minimal g2o_simulator g2o_solver_slam2d_linear g2o_types_icp g2o_types_slam2d g2o_core g2o_interface g2o_solver_csparse g2o_solver_structure_only g2o_types_sba g2o_types_slam3d g2o_csparse_extension g2o_opengl_helper g2o_solver_dense g2o_stuff g2o_types_sclam2d g2o_parser g2o_solver_pcg g2o_types_data g2o_types_sim3 cxsparse g2o_ext_csparse cholmod )

# Include dirs
include_directories(
  include
  stvo-pl/include
  IMU
#  ${Eigen3_INCLUDE_DIRS}
  ${OpenCV_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${YAML_CPP_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/../stvo-pl/include/
  ${PROJECT_SOURCE_DIR}/3rdparty/DBoW2/include/
  ${PROJECT_SOURCE_DIR}/3rdparty/line_descriptor/include/
  ${G2O_INCLUDE_DIR}
  /usr/include/suitesparse # for cholmod
)

# Set link libraries
list(APPEND LINK_LIBS 
  ${OpenCV_LIBS}
  ${Boost_LIBRARIES}
  ${YAML_CPP_LIBRARIES}
  ${G2O_LIBS}
#  ${SOPHUS_LIBRARIES}
#  ${PROJECT_SOURCE_DIR}/../stvo-pl/lib/libstvo.so
  ${PROJECT_SOURCE_DIR}/3rdparty/DBoW2/lib/libDBoW2.so
  ${PROJECT_SOURCE_DIR}/3rdparty/line_descriptor/lib/liblinedesc.so
)

# Set source files 
if(HAS_MRPT)
list(APPEND SOURCEFILES
  src/mapHandler.cpp
  src/mapFeatures.cpp
  src/keyFrame.cpp
  src/slamConfig.cpp
  src/slamScene.cpp
  stvo-pl/src/auxiliar.cpp
  stvo-pl/src/config.cpp
  stvo-pl/src/dataset.cpp
  stvo-pl/src/gridStructure.cpp
  stvo-pl/src/lineIterator.cpp
  stvo-pl/src/matching.cpp
  stvo-pl/src/pinholeStereoCamera.cpp
  stvo-pl/src/sceneRepresentation.cpp
  stvo-pl/src/stereoFeatures.cpp
  stvo-pl/src/stereoFrame.cpp
  stvo-pl/src/stereoFrameHandler.cpp
  stvo-pl/src/timer.cpp
  IMU/so3.cpp
  IMU/g2otypes.cpp
  IMU/imudata.cpp
  IMU/IMUPreintegrator.cpp
  IMU/NavState.cpp
  IMU/types_six_dof_expmap.cpp
  IMU/types_sba.cpp
  IMU/se3_ops.hpp
  IMU/marginalization.cpp
)
else()
list(APPEND SOURCEFILES
  src/mapHandler.cpp
  src/mapFeatures.cpp
  src/keyFrame.cpp
  src/slamConfig.cpp
  stvo-pl/src/auxiliar.cpp
  stvo-pl/src/config.cpp
  stvo-pl/src/dataset.cpp
  stvo-pl/src/gridStructure.cpp
  stvo-pl/src/lineIterator.cpp
  stvo-pl/src/matching.cpp
  stvo-pl/src/pinholeStereoCamera.cpp
 # stvo-pl/src/sceneRepresentation.cpp
  stvo-pl/src/stereoFeatures.cpp
  stvo-pl/src/stereoFrame.cpp
  stvo-pl/src/stereoFrameHandler.cpp
  stvo-pl/src/timer.cpp
  IMU/so3.cpp
  IMU/g2otypes.cpp
  IMU/imudata.cpp
  IMU/IMUPreintegrator.cpp
  IMU/NavState.cpp
  IMU/types_six_dof_expmap.cpp
  IMU/types_sba.cpp
  IMU/se3_ops.hpp
  IMU/marginalization.cpp
)
endif()

# List all files (headers) contained by StVO-PL library
file(GLOB_RECURSE all_include_files RELATIVE "${CMAKE_SOURCE_DIR}" *.h *.hpp)

# Visualize the files of this directory in IDE creating an custom empty target
add_custom_target( plslam_includes DEPENDS ${all_include_files} SOURCES ${all_include_files} )

# Create StVO-PL library
add_library(plslam SHARED ${SOURCEFILES})

if(HAS_MRPT)
target_link_libraries(plslam ${LINK_LIBS} ${MRPT_LIBS}  )
else()
target_link_libraries(plslam ${LINK_LIBS})
endif()

# Applications 
if(HAS_MRPT)
add_executable       ( plslam_dataset app/plslam_dataset.cpp )
target_link_libraries( plslam_dataset plslam )

#add_executable( test test/test.cpp)
#target_link_libraries(test plslam)
else(HAS_MRPT)
add_executable(plslam_dataset app/plslam_dataset.cpp)
	target_link_libraries(plslam_dataset plslam)
endif(HAS_MRPT)









