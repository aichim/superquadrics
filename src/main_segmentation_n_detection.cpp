#include <pcl/console/parse.h>
#include <pcl/console/print.h>
#include <pcl/visualization/pcl_visualizer.h>


#include "tabletop_segmentation.h"
#include "impl/tabletop_segmentation.hpp"

#include "fit_superquadric_ceres.h"
#include "sample_superquadric_uniform.h"
#include "superquadric_detection.h"


using namespace pcl;

int
main (int argc,
      char **argv)
{
  google::InitGoogleLogging (argv[0]);

  /// Read in command line parameters
  std::string cloud_path = "";
  console::parse_argument (argc, argv, "-cloud", cloud_path);

  if (cloud_path == "")
  {
    PCL_ERROR ("Syntax: %s -cloud <path_to_rgbd_cloud.pcd>\n", argv[0]);
    return (-1);
  }

  sq::SuperquadricParameters<double> params_to_detect;
  console::parse_argument (argc, argv, "-e1", params_to_detect.e1);
  console::parse_argument (argc, argv, "-e2", params_to_detect.e2);
  console::parse_argument (argc, argv, "-a", params_to_detect.a);
  console::parse_argument (argc, argv, "-b", params_to_detect.b);
  console::parse_argument (argc, argv, "-c", params_to_detect.c);



  /// Read in point cloud
  PointCloud<PointXYZ>::Ptr cloud_input (new PointCloud<PointXYZ> ());
  io::loadPCDFile (cloud_path, *cloud_input);

  sq::TabletopSegmentation<PointXYZ> tabletop_segmentation;
  tabletop_segmentation.setInputCloud (cloud_input);
  tabletop_segmentation.process ();

  std::vector<PointCloud<PointXYZ>::Ptr> cloud_clusters;
  tabletop_segmentation.getObjectClusters (cloud_clusters);

  visualization::PCLVisualizer visualizer;
  visualizer.setBackgroundColor (255., 255., 255.);

  for (size_t c_i = 0; c_i < cloud_clusters.size (); ++c_i)
  {
    char str[512];
    sprintf (str, "cluster_%03zu", c_i);

    visualization::PointCloudColorHandlerRandom<PointXYZ> color_random (cloud_clusters[c_i]);
    visualizer.addPointCloud<PointXYZ> (cloud_clusters[c_i], color_random, str);
  }

  PolygonMesh::Ptr plane_mesh;
  tabletop_segmentation.getPlaneMesh (plane_mesh);
  visualizer.addPolygonMesh (*plane_mesh, "plane_mesh");
  PCL_INFO ("Tabletop clustering done, press 'q' to detect superquadrics\n");
  visualizer.spin ();

  for (size_t c_i = 0; c_i < cloud_clusters.size (); ++c_i)
  {
    sq::SuperquadricDetection<PointXYZ, double> detection;
    detection.setInputCloud (cloud_clusters[c_i]);
    detection.setSuperquadricParams (params_to_detect);
    std::vector<sq::SuperquadricDetection<PointXYZ, double>::SuperquadricDetectionHypothesis> hypotheses;
    detection.process (hypotheses);

    /// TODO need to sort the hypotheses
    /// Could just hack it and take the best fitting error instead
/*
    double min_fit_error = std::numeric_limits<double>::max ();
    size_t min_fit_error_i = 0;
    for (size_t h_i = 0; h_i < hypotheses.size (); ++h_i)
    {
      if (min_fit_error > hypotheses[h_i].fit_error)
      {
        min_fit_error = hypotheses[h_i].fit_error;
        min_fit_error_i = h_i;
      }
    }*/

    int max_num_surface_points = std::numeric_limits<int>::min ();
    size_t min_fit_error_i = 0;
    for (size_t h_i = 0; h_i < hypotheses.size (); ++h_i)
    {
      if (max_num_surface_points < hypotheses[h_i].num_surface_points)
      {
        max_num_surface_points = hypotheses[h_i].num_surface_points;
        min_fit_error_i = h_i;
      }
    }


    sq::SuperquadricSampling<PointXYZ, double> sampling;
    sq::SuperquadricParameters<double> params = params_to_detect;
    params.transform = hypotheses[min_fit_error_i].transform;
    sampling.setParameters (params);

    PolygonMesh mesh;
    sampling.generateMesh (mesh);
    char str[512];
    sprintf (str, "hypothesis%03zu", c_i);
    visualizer.addPolygonMesh (mesh, str);

    /*
    for (size_t h_i = 0; h_i < hypotheses.size (); ++h_i)
    {
      sq::SuperquadricSampling<PointXYZ, double> sampling;
      sq::SuperquadricParameters<double> params = params_to_detect;
      params.transform = hypotheses[h_i].transform;
      sampling.setParameters (params);

      PolygonMesh mesh;
      sampling.generateMesh (mesh);
      char str[512];
      sprintf (str, "hypothesis%03zu%03zu", c_i, h_i);
      visualizer.addPolygonMesh (mesh, str);
    }*/

    visualizer.spin ();
  }



  return (0);
}
