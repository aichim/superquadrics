#include <pcl/console/parse.h>
#include <pcl/console/print.h>
#include <pcl/visualization/pcl_visualizer.h>


#include "tabletop_segmentation.h"
#include "impl/tabletop_segmentation.hpp"

#include "fit_superquadric_ceres.h"
#include "sample_superquadric_uniform.h"
#include "superquadric_detection.h"

#include <pcl/io/vtk_lib_io.h>

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

  std::vector<double> transf_vec;
  console::parse_x_arguments (argc, argv, "-transform", transf_vec);

  if (transf_vec.size () != 0 && transf_vec.size () != 16)
  {
    PCL_ERROR ("Transformation vector/matrix must have 16 entries\n");
    return (-1);
  }

  Eigen::Matrix4d transformation (Eigen::Matrix4d::Identity ());
  if (transf_vec.size () == 16)
  {
    transformation.setZero ();
    transformation (0, 0) = transf_vec[0];
    transformation (0, 1) = transf_vec[1];
    transformation (0, 2) = transf_vec[2];
    transformation (0, 3) = transf_vec[3];

    transformation (1, 0) = transf_vec[4];
    transformation (1, 1) = transf_vec[5];
    transformation (1, 2) = transf_vec[6];
    transformation (1, 3) = transf_vec[7];

    transformation (2, 0) = transf_vec[8];
    transformation (2, 1) = transf_vec[9];
    transformation (2, 2) = transf_vec[10];
    transformation (2, 3) = transf_vec[11];

    transformation (3, 0) = transf_vec[12];
    transformation (3, 1) = transf_vec[13];
    transformation (3, 2) = transf_vec[14];
    transformation (3, 3) = transf_vec[15];

    std::cout << "using transformation:\n" << transformation << std::endl;
  }



  /// Read in point cloud
  PointCloud<PointXYZ>::Ptr cloud_input (new PointCloud<PointXYZ> ());
  io::loadPCDFile (cloud_path, *cloud_input);

  sq::SuperquadricRigidRegistration<PointXYZ, double> sq_reg;
  sq_reg.setInputCloud (cloud_input);
  sq_reg.setParameters (params_to_detect.e1, params_to_detect.e2,
                        params_to_detect.a, params_to_detect.b, params_to_detect.c);


  double min_fit = std::numeric_limits<double>::max ();
  Eigen::Matrix4d min_transf;

  /// TODO 0 - 3
  for (int i = 0; i < 3; ++i)
  {
    sq_reg.setPreAlign (true, i);
    Eigen::Matrix4d transf;
    double fit = sq_reg.fit (transf);
    printf ("pre_align axis %d, fit %f\n", i, fit);

    if (fit < min_fit)
    {
      min_fit = fit;
      min_transf = transf;
    }


    PointCloud<PointXYZ>::Ptr cloud_tr (new PointCloud<PointXYZ> ());
    transformPointCloud (*cloud_input, *cloud_tr, transf);
    char str[512];
    sprintf (str, "cloud_trial_%d.pcd", i);
    io::savePCDFileBinaryCompressed (str, *cloud_tr);
  }


  std::cout << "transform result:\n" << min_transf << std::endl;


  PointCloud<PointXYZ>::Ptr cloud_tr (new PointCloud<PointXYZ> ());
  transformPointCloud (*cloud_input, *cloud_tr, min_transf);
  io::savePCDFileBinaryCompressed ("cloud_tr.pcd", *cloud_tr);

  sq::SuperquadricSampling<PointXYZ, double> sampling;
//  params_to_detect.transform = Eigen::Matrix4d::Identity ();
  params_to_detect.transform = min_transf;
  sampling.setParameters (params_to_detect);
  PointCloud<PointXYZ>::Ptr cloud_obj (new PointCloud<PointXYZ> ());
  sampling.generatePointCloud (*cloud_obj);
  io::savePCDFileBinaryCompressed ("cloud_obj.pcd", *cloud_obj);




  //// using the detection pipeline


  sq::SuperquadricDetection<PointXYZ, double> detection;
  detection.setInputCloud (cloud_input);
  detection.setSuperquadricParams (params_to_detect);
  std::vector<sq::SuperquadricDetection<PointXYZ, double>::SuperquadricDetectionHypothesis> hypotheses;
  detection.process (hypotheses);

  double min_fit_error = std::numeric_limits<double>::max ();
  size_t min_fit_error_i = 0;
  for (size_t h_i = 0; h_i < hypotheses.size (); ++h_i)
  {
    if (min_fit_error > hypotheses[h_i].fit_error)
    {
      min_fit_error = hypotheses[h_i].fit_error;
      min_fit_error_i = h_i;
    }
  }

  sq::SuperquadricSampling<PointXYZ, double> sampling_mesh;
  sq::SuperquadricParameters<double> params = params_to_detect;
  params.transform = hypotheses[min_fit_error_i].transform;
  sampling_mesh.setParameters (params);
  PolygonMesh mesh;
  sampling_mesh.generateMesh (mesh);
  io::savePolygonFileVTK ("cloud_detected.vtk", mesh);


  return (0);
}
