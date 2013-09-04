cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(superquadrics_visualize)

find_package(PCL 1.7 REQUIRED io common segmentation sample_consensus surface visualization)
include_directories(${PCL_INCLUDE_DIRS} include/)
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

find_package (Ceres REQUIRED)
include_directories (${CERES_INCLUDES})


find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set (framework_incs include/fit_superquadric_lm.h 
					include/fit_superquadric_ceres.h
                    include/superquadric_formulas.h
                    include/sample_superquadric.h
                    include/sample_superquadric_uniform.h
                    include/segmentation_utils.h
                    include/tabletop_segmentation.h)
set (framework_impls include/impl/fit_superquadric_lm.hpp
					 include/impl/fit_superquadric_ceres.hpp
					 include/impl/superquadric_formulas.hpp
                     include/impl/sample_superquadric.hpp
                     include/impl/sample_superquadric_uniform.hpp
                     include/impl/segmentation_utils.hpp
                     include/impl/tabletop_segmentation.hpp)
set (framework_srcs src/fit_superquadric_lm.cpp
                    src/fit_superquadric_ceres.cpp
                    src/sample_superquadric.cpp
                    src/sample_superquadric_uniform.cpp
                    src/segmentation_utils.cpp
                    src/tabletop_segmentation.cpp)

ADD_LIBRARY (superquadrics_framework ${framework_incs} ${framework_impls} ${framework_srcs})
TARGET_LINK_LIBRARIES (superquadrics_framework)

set (libs_to_link ${PCL_LIBRARIES} ${CERES_LIBRARIES} superquadrics_framework)



add_executable (sample_superquadric_vtk MACOSX_BUNDLE src/main_sample_superquadric_vtk.cpp)
target_link_libraries (sample_superquadric_vtk ${PCL_LIBRARIES})

add_executable (sample_superquadric src/main_sample_superquadric.cpp)
target_link_libraries (sample_superquadric ${libs_to_link})

add_executable (fit_superquadric_to_samples MACOSX_BUNDLE src/main_fit_superquadric_to_samples.cpp)
target_link_libraries (fit_superquadric_to_samples ${PCL_LIBRARIES} superquadrics_framework)


add_executable (split_and_merge MACOSX_BUNDLE src/main_split_and_merge.cpp)
target_link_libraries (split_and_merge ${libs_to_link})

add_executable (fit_superquadric_ceres src/main_fit_superquadric_ceres.cpp)
target_link_libraries (fit_superquadric_ceres ${libs_to_link})

add_executable (ceres_test MACOSX_BUNDLE src/main_ceres_test.cpp)
target_link_libraries (ceres_test ${libs_to_link})


add_executable (tabletop_segmentation MACOSX_BUNDLE src/main_tabletop_segmentation.cpp)
target_link_libraries (tabletop_segmentation ${libs_to_link})

add_executable (segmentation_n_fitting_pipeline MACOSX_BUNDLE src/main_segmentation_n_fitting_pipeline.cpp)
target_link_libraries (segmentation_n_fitting_pipeline ${libs_to_link})
