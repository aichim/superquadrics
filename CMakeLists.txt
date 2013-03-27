cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(superquadrics_visualize)

find_package(PCL 1.7 REQUIRED io common visualization)

include_directories(${PCL_INCLUDE_DIRS} include/)
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable (superquadrics_sample MACOSX_BUNDLE src/main_sample.cpp)
target_link_libraries (superquadrics_sample ${PCL_LIBRARIES})

