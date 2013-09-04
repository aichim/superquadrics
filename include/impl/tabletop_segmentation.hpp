#pragma once

#include "tabletop_segmentation.h"


template<typename PointT>
sq::TabletopSegmentation<PointT>::TabletopSegmentation ()
{

}



template<typename PointT> void
sq::TabletopSegmentation<PointT>::process ()
{

  pcl::IndicesPtr indices (new std::vector<int> ());
  for (size_t i = 0; i < cloud_input_->size (); ++i)
    if (pcl::isFinite ((*cloud_input_)[i]))
      indices->push_back (i);

  PCL_INFO ("[TabletopSegmentation] Cloud size %zu, indices size %zu\n", cloud_input_->size (), indices->size ());

  /// Segment the largest plane
  pcl::IndicesPtr plane_indices (new std::vector<int> ());
  Eigen::Vector4f plane_params;
  sq::findLargestPlane<PointT> (cloud_input_, indices, *plane_indices, plane_params);

  PCL_INFO ("[TabletopSegmentation] Found plane of size %zu\n", plane_indices->size ());

  /// Extract the plane cloud
  std::vector<bool> is_planar_point (cloud_input_->size (), false);
  for (size_t i = 0; i < plane_indices->size (); ++i)
    is_planar_point[(*plane_indices)[i]] = true;

  result_plane_.reset (new Cloud ());
  result_non_plane_.reset (new Cloud ());
  for (size_t i = 0; i < cloud_input_->size (); ++i)
    if (is_planar_point[i])
      result_plane_->push_back ((*cloud_input_)[i]);
    else
      result_non_plane_->push_back ((*cloud_input_)[i]);


  CloudPtr plane_contour (new Cloud ());
  sq::computeInliersContour<PointT> (cloud_input_, plane_indices, plane_params, *plane_contour);
  PCL_INFO ("[TabletopSegmentation] Plane contour has %zu points.\n", plane_contour->size ());

  plane_contour->push_back (plane_contour->front ());
  result_plane_mesh_.reset (new pcl::PolygonMesh ());
  sq::triangulizeContour (plane_contour, *result_plane_mesh_);


  /// TODO throw away the points that are not above/below the plane

  /// Cluster the non-planar points
  pcl::EuclideanClusterExtraction<PointT> cluster_extraction;
  cluster_extraction.setClusterTolerance (0.025);
  cluster_extraction.setMinClusterSize (1000);
  cluster_extraction.setMaxClusterSize (std::numeric_limits<int>::max ());
  cluster_extraction.setInputCloud (result_non_plane_);
  std::vector<pcl::PointIndices> clusters;
  cluster_extraction.extract (clusters);
  PCL_INFO ("[TabletopSegmentation] Found %zu clusters.\n", clusters.size ());

  /// Save clusters separately
  result_clusters_.clear ();
  for (size_t c_i = 0; c_i < clusters.size (); ++c_i)
  {
    CloudPtr cloud_cluster (new Cloud ());
    cloud_cluster->reserve (clusters[c_i].indices.size ());
    for (size_t i = 0; i < clusters[c_i].indices.size (); ++i)
      cloud_cluster->push_back ((*result_non_plane_)[clusters[c_i].indices[i]]);
    cloud_cluster->height = 1;
    cloud_cluster->width = cloud_cluster->size ();
    result_clusters_.push_back (cloud_cluster);
  }
}
