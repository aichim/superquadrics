#include <pcl/io/pcd_io.h>

#include <pcl/console/print.h>
#include <pcl/console/parse.h>

using namespace pcl;


int
main (int argc,
      char **argv)
{
  PCL_INFO ("Syntax: %s [-e1 x.x] [-e2 y.y] [-a z.z] [-b t.t] [-c u.u] [-visualize 0/1] [-output xxx.pcd]\n", argv[0]);

  /// Parse command line arguments
  bool visualize_enabled = true;
  console::parse_argument (argc, argv, "-visualize", visualize_enabled);

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

  std::string output_file = "";
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


  /// Sample the superquadric to a pcd file
  if (output_file != "")
  {
    PointCloud<PointXYZ> cloud;
    for (double eta = -M_PI / 2.0; eta < M_PI / 2.; eta += 0.01)
    {
      for (double mu = -M_PI; mu < M_PI; mu += 0.01)
      {
        PointXYZ point;
        Eigen::Vector4d p;
        p[0] = a * (cos(eta) / fabs(cos(eta))) * pow (fabs(cos (eta)), epsilon_1) * (cos(mu) / fabs(cos(mu))) * pow (fabs (cos (mu)), epsilon_2);
        p[1] = b * (cos(eta) / fabs(cos(eta))) * pow (fabs(cos (eta)), epsilon_1) * (sin(mu) / fabs(sin(mu))) * pow (fabs (sin (mu)), epsilon_2);
        p[2] = c * (sin(eta) / fabs(sin(eta))) * pow (fabs(sin (eta)), epsilon_1);
        p[3] = 1.;

        Eigen::Vector4d p_tr = transformation.inverse () * p;
        point.x = p_tr[0];
        point.y = p_tr[1];
        point.z = p_tr[2];

        if (isFinite (point))
          cloud.push_back (point);
      }
    }

    cloud.height = 1.;
    cloud.width = cloud.size ();
    io::savePCDFile (output_file, cloud, true);
  }



  return EXIT_SUCCESS;
}
