#pragma once


#include "segmentation_utils.h"

namespace sq
{

template<typename PointT>
class TabletopSegmentation
{
public:
  typedef pcl::PointCloud<PointT> Cloud;
  typedef typename Cloud::Ptr CloudPtr;
  typedef typename Cloud::ConstPtr CloudConstPtr;

  TabletopSegmentation ();

  void
  setInputCloud (CloudConstPtr cloud_input)
  { cloud_input_ = cloud_input; }


  void
  process ();

  void
  getObjectClusters (std::vector<CloudPtr> &clusters)
  { clusters = result_clusters_; }

  void
  getPlane (CloudPtr &plane)
  { plane = result_plane_; }

  void
  getPlaneMesh (pcl::PolygonMesh::Ptr &plane_mesh)
  { plane_mesh = result_plane_mesh_; }

  void
  getNonPlane (CloudPtr &non_plane)
  { non_plane = result_non_plane_; }

  void
  getOthers (CloudPtr &others)
  { others = result_others_; }


protected:
  CloudConstPtr cloud_input_;

  std::vector<CloudPtr> result_clusters_;
  CloudPtr result_plane_;
  CloudPtr result_non_plane_;
  CloudPtr result_others_;
  pcl::PolygonMesh::Ptr result_plane_mesh_;
};
}


#include "impl/segmentation_utils.hpp"
