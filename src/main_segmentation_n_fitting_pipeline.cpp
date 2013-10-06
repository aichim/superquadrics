#include <pcl/console/parse.h>
#include <pcl/console/print.h>
#include <pcl/visualization/pcl_visualizer.h>


#include "tabletop_segmentation.h"
#include "impl/tabletop_segmentation.hpp"

#include "fit_superquadric_ceres.h"
#include "sample_superquadric_uniform.h"

#include <pcl/console/time.h>

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


  /// Read in point cloud
  PointCloud<PointXYZ>::Ptr cloud_input (new PointCloud<PointXYZ> ());
  io::loadPCDFile (cloud_path, *cloud_input);

  pcl::console::TicToc timer;
  timer.tic ();
  sq::TabletopSegmentation<PointXYZ> tabletop_segmentation;
  tabletop_segmentation.setInputCloud (cloud_input);
  tabletop_segmentation.process ();
  double time = timer.toc ();

  PCL_ERROR ("Segmentation time: %f\n", time);


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
  PCL_INFO ("Tabletop clustering done, press 'q' to fit superquadrics\n");
  visualizer.spin ();

  pcl::console::TicToc timer_fitting;
  timer_fitting.tic ();
  for (size_t c_i = 0; c_i < cloud_clusters.size (); ++c_i)
  {
    timer.tic ();

    double min_fit = std::numeric_limits<double>::max ();
    sq::SuperquadricParameters<double> min_params;

    for (size_t i = 0; i < 3; ++i)
    {
      sq::SuperquadricFittingCeres<PointXYZ> fitting;
      fitting.setInputCloud (cloud_clusters[c_i]);
      fitting.setPreAlign (true, i);

      sq::SuperquadricParameters<double> params;
      double fit = fitting.fit (params);
      printf ("pre_align axis %d, fit %f\n", i, fit);

      if (fit < min_fit)
      {
        min_fit = fit;
        min_params = params;
      }
    }

    time = timer.toc ();
    PCL_ERROR ("3 fittings: %f\n", time);

    if (min_fit > std::numeric_limits<double>::max () - 100.)
      continue;

//    sq::SuperquadricSamplingUniform<PointXYZ, double> sampling;
//    sampling.setSpatialSampling (0.001);
    sq::SuperquadricSampling<PointXYZ, double> sampling;
    sampling.setParameters (min_params);

    PolygonMesh mesh;
    sampling.generateMesh (mesh);
    char str[512];
    sprintf (str, "cluster_fitted_mesh_%03zu", c_i);
    visualizer.addPolygonMesh (mesh, str);

    PCL_INFO ("--- Displaying superquadric with error: %f and params: -e1 %f -e2 %f -a %f -b %f -c %f\n",
              min_fit, min_params.e1, min_params.e2, min_params.a, min_params.b, min_params.c);

//    PointCloud<PointXYZ>::Ptr cloud_fitted (new PointCloud<PointXYZ> ());
//    sampling.generatePointCloud (*cloud_fitted);

//    char str[512];
//    sprintf (str, "cluster_fitted_cloud_%03zu", c_i);
//    visualization::PointCloudColorHandlerCustom<PointXYZ> color_cloud (cloud_fitted, 255., 0., 0.);
//    visualizer.addPointCloud<PointXYZ> (cloud_fitted, color_cloud, str);
//    io::savePCDFile ("caca.pcd", *cloud_fitted, true);
  }

  time = timer_fitting.toc ();
  PCL_ERROR ("Total fitting time: %f\n", time);

  PCL_INFO ("Superquadric fitting done. Press 'q' to quit.\n");
  visualizer.spin ();



  return (0);
}
