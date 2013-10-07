#include <pcl/io/pcd_io.h>
#include <pcl/io/obj_io.h>
#include <pcl/io/vtk_io.h>

#include <pcl/console/print.h>
#include <pcl/console/parse.h>

#include <boost/algorithm/string.hpp>

#include "sample_superquadric.h"
#include "sample_superquadric_uniform.h"
#include "superquadric_formulas.h"

using namespace pcl;


int
main (int argc,
      char **argv)
{
  PCL_INFO ("Syntax: %s [-e1 x.x] [-e2 y.y] [-a z.z] [-b t.t] [-c u.u] [-output xxx.pcd]\n", argv[0]);

  /// Parse command line arguments
  double epsilon_1 = 0.;
  console::parse_argument (argc, argv, "-e1", epsilon_1);

  double epsilon_2 = 0.;
  console::parse_argument (argc, argv, "-e2", epsilon_2);

  double a = 1.;
  console::parse_argument (argc, argv, "-a", a);

  double b = 1.;
  console::parse_argument (argc, argv, "-b", b);

  double c = 1.;
  console::parse_argument (argc, argv, "-c", c);

  std::string output_file = "sample.pcd";
  console::parse_argument (argc, argv, "-output", output_file);

  std::vector<double> transf_vec;
  console::parse_x_arguments (argc, argv, "-transform", transf_vec);

  if (transf_vec.size () != 0 && transf_vec.size () != 16)
  {
    PCL_ERROR ("Transformation vector/matrix must have 16 entries\n");
    return (-1);
  }

  PCL_INFO ("### Superquadric parameters:\n   epsilon_1: %f\n   epsilon_2: %f\n   a: %f\n   b: %f\n   c: %f\n",
            epsilon_1, epsilon_2, a, b, c);

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


  /// Prepare the sampling instance
  sq::SuperquadricSampling<PointXYZ, double> sampling;
  sq::SuperquadricParameters<double> params;
  params.e1 = epsilon_1;
  params.e2 = epsilon_2;
  params.a = a;
  params.b = b;
  params.c = c;
  params.transform = transformation;

  sampling.setParameters (params);


  /// Sample the superquadric to a pcd file
  std::string extension = output_file.substr (output_file.size () - 3);
  boost::algorithm::to_lower (extension);
  if (extension == "pcd")
  {
    PointCloud<PointXYZ> cloud;
    sampling.generatePointCloud (cloud);
    io::savePCDFile (output_file, cloud, true);
  }
  else if (extension == "obj" ||
           extension == "vtk")
  {
    PolygonMesh mesh;
    sampling.generateMesh (mesh);

    if (extension == "obj")
      io::saveOBJFile (output_file, mesh);
    else if (extension == "vtk")
      io::saveVTKFile (output_file, mesh);
  }
  else
  {
    PCL_ERROR ("Unknown file format, use .pcd or .obj\n");
    return (-1);
  }




  /// Generate superquadric with the new smart sampling scheme
  sq::SuperquadricSamplingUniform<PointXYZ, double> sampling_uniform;
  sampling_uniform.setParameters (params);
  sampling_uniform.setSpatialSampling (0.07);
  PointCloud<PointXYZ> cloud_uniform;
  sampling_uniform.generatePointCloud (cloud_uniform);
  io::savePCDFile ("uniform.pcd", cloud_uniform, true);



  return EXIT_SUCCESS;
}
