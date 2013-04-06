cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(superquadrics_visualize)

find_package(PCL 1.7 REQUIRED io common visualization)

include_directories(${PCL_INCLUDE_DIRS} include/)
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})


set (framework_incs include/fit_superquadric_lm.h include/superquadric_formulas.h)
set (framework_impls include/impl/fit_superquadric_lm.hpp)
set (framework_srcs src/fit_superquadric_lm.cpp)

ADD_LIBRARY (superquadrics_framework ${framework_incs} ${framework_impls} ${framework_srcs})
TARGET_LINK_LIBRARIES (superquadrics_framework ${libs_to_link})


add_executable (superquadrics_sample MACOSX_BUNDLE src/main_sample.cpp)
target_link_libraries (superquadrics_sample ${PCL_LIBRARIES})

add_executable (superquadrics_sample_novtk MACOSX_BUNDLE src/main_sample_novtk.cpp)
target_link_libraries (superquadrics_sample_novtk ${PCL_LIBRARIES})

add_executable (fit_superquadric_to_samples MACOSX_BUNDLE src/main_fit_superquadric_to_samples.cpp)
target_link_libraries (fit_superquadric_to_samples ${PCL_LIBRARIES} superquadrics_framework)


add_executable (split_and_merge MACOSX_BUNDLE src/main_split_and_merge.cpp)
target_link_libraries (split_and_merge ${PCL_LIBRARIES} superquadrics_framework)
