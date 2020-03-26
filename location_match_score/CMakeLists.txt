cmake_minimum_required(VERSION 3.1)
project(location_match_score)

if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 11)
endif()

find_package(catkin REQUIRED
  COMPONENTS
    cmake_modules
    message_filters
    nav_msgs
    roscpp
    sensor_msgs
    std_srvs
    tf
    tf2
    tf2_msgs
    tf2_ros
)

find_package(Boost REQUIRED )
find_package(Eigen3 REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(OpenCV 3 REQUIRED) 

include_directories(include)
include_directories(${catkin_INCLUDE_DIRS} 
                    ${EIGEN_INCLUDE_DIRS} 
                    ${Boost_INCLUDE_DIRS}
                    ${OpenCV_INCLUDE_DIRS})

add_definitions(${EIGEN_DEFINITIONS})

catkin_package(
  CATKIN_DEPENDS
    roscpp
    tf
    tf2
    tf2_msgs
    tf2_ros
    sensor_msgs
    message_filters
    nav_msgs 
    std_srvs
  INCLUDE_DIRS include
)

add_library(scan_match src/scan_match.cpp)
target_link_libraries(scan_match ${OpenCV_LIBS} yaml-cpp)

install(TARGETS
    scan_match
    DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)


add_executable(location_match_score src/location_match_score.cpp)

target_link_libraries(location_match_score 
                      scan_match
                      ${catkin_LIBRARIES}
                      ${Boost_LIBRARIES}
                      )


install(TARGETS location_match_score
    RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)


