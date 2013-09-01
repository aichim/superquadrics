#pragma once

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/pcl_base.h>
#include <pcl/PolygonMesh.h>

namespace sq
{

template<typename PointT> bool
findLargestPlane (typename pcl::PointCloud<PointT>::ConstPtr cloud,
                  pcl::IndicesPtr &indices,
                  std::vector<int> &inlier_indices,
                  Eigen::Vector4f &plane_params);




template<typename PointT> void
computeInliersContour (typename pcl::PointCloud<PointT>::ConstPtr cloud,
                       pcl::IndicesPtr &inlier_indices,
                       const Eigen::Vector4f &plane_params,
                       pcl::PointCloud<pcl::PointXYZ> &contour);


void
triangulizeContour (pcl::PointCloud<pcl::PointXYZ>::ConstPtr contour,
                    pcl::PolygonMesh &mesh);
}

#include "impl/segmentation_utils.hpp"
