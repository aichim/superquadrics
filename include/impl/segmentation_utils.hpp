#pragma once

#include "segmentation_utils.h"

#include <pcl/sample_consensus/sac_model_normal_plane.h>
#include <pcl/sample_consensus/ransac.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/segmentation/extract_clusters.h>
#include <pcl/features/normal_3d.h>

#include <pcl/surface/concave_hull.h>
#include <pcl/geometry/polygon_operations.h>
#include <pcl/filters/project_inliers.h>

#include <pcl/io/vtk_lib_io.h>



template<typename PointT> bool
sq::findLargestPlane (typename pcl::PointCloud<PointT>::ConstPtr cloud,
                      pcl::IndicesPtr &indices,
                      std::vector<int> &inlier_indices,
                      Eigen::Vector4f &plane_params)
{
  /// HACK: it seems that passing indices to the sac_model does not influence the results
  typename pcl::PointCloud<PointT>::Ptr cloud_new (new pcl::PointCloud<PointT> ());
  pcl::ExtractIndices<PointT> ei;
  ei.setKeepOrganized (true);
  ei.setInputCloud (cloud);
  ei.setIndices (indices);
  ei.filter (*cloud_new);

  //    typename pcl::SampleConsensusModelPlane<PointT>::Ptr model (new pcl::SampleConsensusModelPlane<PointT> (cloud, indices));
  typename pcl::SampleConsensusModelPlane<PointT>::Ptr model (new pcl::SampleConsensusModelPlane<PointT> (cloud_new));
  pcl::RandomSampleConsensus<PointT> sac (model, 0.025);
  sac.setMaxIterations (1500);

  bool res = sac.computeModel ();

  if (!res)
  {
    PCL_ERROR ("No planar model found. Relax thresholds or select object again and continue.\n");
    return (false);
  }

  res = sac.refineModel (3., 500);

  if (!res)
  {
    PCL_ERROR ("No planar model found. Relax thresholds or select object again and continue.\n");
    return (false);
  }

  sac.getInliers (inlier_indices);
  Eigen::VectorXf aux;
  sac.getModelCoefficients (aux);
  plane_params = aux.head<4> ();



  /// Cluster the inlier indices to take only the contiguous plane
  std::vector<pcl::PointIndices> clusters;
  typename pcl::search::KdTree<PointT>::Ptr tree (new pcl::search::KdTree<PointT> ());
  tree->setInputCloud (cloud, pcl::IndicesPtr (new std::vector<int> (inlier_indices)));
  pcl::extractEuclideanClusters<PointT> (*cloud, inlier_indices, tree, 0.03, clusters);

  /// Keep only the largest cluster
  size_t max_cluster_size = 0, max_cluster_i = 0;
  for (size_t i = 0; i < clusters.size (); ++i)
    if (clusters[i].indices.size () > max_cluster_size)
    {
      max_cluster_size = clusters[i].indices.size ();
      max_cluster_i = i;
    }
  inlier_indices = clusters[max_cluster_i].indices;

  /// Recompute the plane equation for this cluster only
  float curvature_aux;
  pcl::computePointNormal<PointT> (*cloud, inlier_indices, plane_params, curvature_aux);

  /// Orient the normal towards the viewpoint
  if (plane_params.head<3> ().dot (Eigen::Vector3f (0.f, 0.f, 1.f)) < 0.f)
    plane_params *= -1.0;


  return (res);
}



template<typename PointT> void
sq::computeInliersContour (typename pcl::PointCloud<PointT>::ConstPtr cloud,
                           pcl::IndicesPtr &inlier_indices,
                           const Eigen::Vector4f &plane_params,
                           pcl::PointCloud<pcl::PointXYZ> &contour)
{
  /// First project the inliers to the plane => kill one eigenvalue of the covariance matrix
  pcl::ProjectInliers<PointT> project_inliers;
  project_inliers.setInputCloud (cloud);
  project_inliers.setIndices (inlier_indices);
  pcl::ModelCoefficientsPtr model_coeffs (new pcl::ModelCoefficients ());
  model_coeffs->values.resize (4);
  for (size_t i = 0; i < 4; ++i)
    model_coeffs->values[i] = plane_params[i];
  project_inliers.setModelCoefficients (model_coeffs);
  typename pcl::PointCloud<PointT>::Ptr points_projected (new pcl::PointCloud<PointT> ());
  project_inliers.filter (*points_projected);

  pcl::ConcaveHull<PointT> concave_hull_alg;
  concave_hull_alg.setInputCloud (points_projected);
  concave_hull_alg.setDimension (2);
  concave_hull_alg.setAlpha (0.2);
  pcl::PolygonMesh mesh_hull;
  concave_hull_alg.reconstruct (mesh_hull);

  pcl::io::savePolygonFileVTK ("plane_cv.vtk", mesh_hull);

  pcl::PointCloud<pcl::PointXYZ>::Ptr hull_points (new pcl::PointCloud<pcl::PointXYZ> ());
  pcl::fromPCLPointCloud2 (mesh_hull.cloud, *hull_points);

  /// Assume all of the points are already in the plane
  /// Compute the uv-coordinates for the 3d points in the plane
  Eigen::Vector3d normal = plane_params.head<3> ().cast<double> ().normalized ();
  Eigen::Vector3d u_vec = normal.unitOrthogonal ().normalized ();
  Eigen::Vector3d v_vec = normal.cross (u_vec).normalized ();
  Eigen::Matrix3d proj_matrix;
  for (size_t i = 0; i < 3; ++i)
  {
    proj_matrix (0, i) = u_vec[i];
    proj_matrix (1, i) = v_vec[i];
    proj_matrix (2, i) = normal[i];
  }


  pcl::PointCloud<pcl::PointXYZ>::Ptr points_uv (new pcl::PointCloud<pcl::PointXYZ> ());
  points_uv->resize (hull_points->size ());
  for (size_t i = 0; i < points_uv->size (); ++i)
  {
    Eigen::Vector3d hp ((*hull_points)[i].x, (*hull_points)[i].y, (*hull_points)[i].z);
    Eigen::Vector3d aux = proj_matrix * (hp);

    (*points_uv)[i].x = aux[0];
    (*points_uv)[i].y = aux[1];
    (*points_uv)[i].z = aux[2];//0.f; //aux[2];//0.f;  // aux[2] + plane_params[3];
  }

  /// Simplify the contour
//  contour.points = points_uv->points;
  pcl::approximatePolygon2D<pcl::PointXYZ> (points_uv->points, contour.points, 0.05, false);


  /// Un-project points back to 3D
  proj_matrix = Eigen::Matrix3d (proj_matrix.inverse ());
  for (size_t i = 0; i < contour.size (); ++i)
  {
    Eigen::Vector3d aux = proj_matrix * Eigen::Vector3d (contour[i].x, contour[i].y, -plane_params[3]); //contour[i].z);//-plane_params[3]);
    contour[i].x = aux[0];
    contour[i].y = aux[1];
    contour[i].z = aux[2];
  }
}
