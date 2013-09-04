#include <pcl/console/parse.h>
#include <pcl/console/print.h>

#include <pcl/io/pcd_io.h>
#include <pcl/io/vtk_lib_io.h>
#include <pcl/segmentation/extract_clusters.h>

#include "segmentation_utils.h"

using namespace pcl;

int
main (int argc,
      char **argv)
{
  /// Read in command line parameters
  std::string cloud_path = "";
  console::parse_argument (argc, argv, "-cloud", cloud_path);

  if (cloud_path == "")
  {
    PCL_ERROR ("Syntax: %s -cloud <path_to_rgbd_cloud.pcd>\n", argv[0]);
    return (-1);
  }


  /// Read in point cloud
  PointCloud<PointXYZ>::Ptr cloud_input (new PointCloud<PointXYZ> ());
  io::loadPCDFile (cloud_path, *cloud_input);

  IndicesPtr indices (new std::vector<int> ());
  for (size_t i = 0; i < cloud_input->size (); ++i)
    if (isFinite ((*cloud_input)[i]))
      indices->push_back (i);

  PCL_INFO ("Cloud size %zu, indices size %zu\n", cloud_input->size (), indices->size ());

  /// Segment the largest plane
  IndicesPtr plane_indices (new std::vector<int> ());
  Eigen::Vector4f plane_params;
  sq::findLargestPlane<PointXYZ> (cloud_input, indices, *plane_indices, plane_params);

  PCL_INFO ("Found plane of size %zu\n", plane_indices->size ());

  /// Extract the plane cloud
  std::vector<bool> is_planar_point (cloud_input->size (), false);
  for (size_t i = 0; i < plane_indices->size (); ++i)
    is_planar_point[(*plane_indices)[i]] = true;

  PointCloud<PointXYZ>::Ptr cloud_plane (new PointCloud<PointXYZ> ());
  PointCloud<PointXYZ>::Ptr cloud_non_plane (new PointCloud<PointXYZ> ());
  for (size_t i = 0; i < cloud_input->size (); ++i)
    if (is_planar_point[i])
      cloud_plane->push_back ((*cloud_input)[i]);
    else
      cloud_non_plane->push_back ((*cloud_input)[i]);


  io::savePCDFile ("plane.pcd", *cloud_plane, true);
  io::savePCDFile ("non-plane.pcd", *cloud_non_plane, true);

  PointCloud<PointXYZ>::Ptr plane_contour (new PointCloud<PointXYZ> ());
  sq::computeInliersContour<PointXYZ> (cloud_input, plane_indices, plane_params, *plane_contour);
  PCL_INFO ("Plane contour has %zu points.\n", plane_contour->size ());
  io::savePCDFile ("plane_contour.pcd", *plane_contour, true);

  plane_contour->push_back (plane_contour->front ());
  PolygonMesh plane_mesh;
  sq::triangulizeContour (plane_contour, plane_mesh);
  io::savePolygonFileVTK ("plane_mesh.vtk", plane_mesh);

  /// TODO throw away the points that are not above/below the plane

  /// Cluster the non-planar points
  EuclideanClusterExtraction<PointXYZ> cluster_extraction;
  cluster_extraction.setClusterTolerance (0.025);
  cluster_extraction.setMinClusterSize (1000);
  cluster_extraction.setMaxClusterSize (std::numeric_limits<int>::max ());
  cluster_extraction.setInputCloud (cloud_non_plane);
  std::vector<PointIndices> clusters;
  cluster_extraction.extract (clusters);

  /// Save clusters separately
  PCL_INFO ("Found %zu clusters.\n", clusters.size ());
  for (size_t c_i = 0; c_i < clusters.size (); ++c_i)
  {
    PointCloud<PointXYZ>::Ptr cloud_cluster (new PointCloud<PointXYZ> ());
    cloud_cluster->reserve (clusters[c_i].indices.size ());
    for (size_t i = 0; i < clusters[c_i].indices.size (); ++i)
      cloud_cluster->push_back ((*cloud_non_plane)[clusters[c_i].indices[i]]);
    cloud_cluster->height = 1;
    cloud_cluster->width = cloud_cluster->size ();

    char str[512];
    sprintf (str, "cluster_%03zu.pcd", c_i);
    io::savePCDFile (str, *cloud_cluster, true);
  }


  return (0);
}
