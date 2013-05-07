cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(superquadrics_visualize)

find_package(PCL 1.7 REQUIRED io common visualization)
include_directories(${PCL_INCLUDE_DIRS} include/)
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

find_package (Ceres REQUIRED)
include_directories (${CERES_INCLUDES})


set (framework_incs include/fit_superquadric_lm.h 
					include/fit_superquadric_ceres.h
                    include/superquadric_formulas.h
                    include/sample_superquadric.h)
set (framework_impls include/impl/fit_superquadric_lm.hpp
					 include/impl/fit_superquadric_ceres.hpp
                     include/impl/sample_superquadric.hpp)
set (framework_srcs src/fit_superquadric_lm.cpp
                    src/fit_superquadric_ceres.cpp
                    src/sample_superquadric.cpp)

ADD_LIBRARY (superquadrics_framework ${framework_incs} ${framework_impls} ${framework_srcs})
TARGET_LINK_LIBRARIES (superquadrics_framework)

set (libs_to_link ${PCL_LIBRARIES} ${CERES_LIBRARIES} superquadrics_framework)



add_executable (superquadrics_sample MACOSX_BUNDLE src/main_sample.cpp)
target_link_libraries (superquadrics_sample ${PCL_LIBRARIES})

add_executable (superquadrics_sample_novtk MACOSX_BUNDLE src/main_sample_novtk.cpp)
target_link_libraries (superquadrics_sample_novtk ${PCL_LIBRARIES})

add_executable (fit_superquadric_to_samples MACOSX_BUNDLE src/main_fit_superquadric_to_samples.cpp)
target_link_libraries (fit_superquadric_to_samples ${PCL_LIBRARIES} superquadrics_framework)


add_executable (split_and_merge MACOSX_BUNDLE src/main_split_and_merge.cpp)
target_link_libraries (split_and_merge ${PCL_LIBRARIES} superquadrics_framework)

add_executable (fit_superquadric_ceres src/main_fit_superquadric_ceres.cpp)
target_link_libraries (fit_superquadric_ceres ${libs_to_link})

add_executable (ceres_test MACOSX_BUNDLE src/main_ceres_test.cpp)
target_link_libraries (ceres_test ${libs_to_link})
