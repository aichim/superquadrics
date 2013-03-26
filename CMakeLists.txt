cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(superquadrics_visualize)

find_package(PCL 1.7 REQUIRED io common visualization)

include_directories(${PCL_INCLUDE_DIRS} include/)
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable (superquadrics_visualize MACOSX_BUNDLE src/main_visualize.cpp)
target_link_libraries (superquadrics_visualize ${PCL_LIBRARIES})

