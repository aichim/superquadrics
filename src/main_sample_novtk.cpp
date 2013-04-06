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

  if (transf_vec.size () != 0 && transf_vec.size () != 6)
  {
    PCL_ERROR ("Transformation vector must have 6 entries\n");
    return (-1);
  }

  PCL_INFO ("### Superquadric parameters:\n   epsilon_1: %f\n   epsilon_2: %f\n   a: %f\n   b: %f\n   c: %f\n",
            epsilon_1, epsilon_2, a, b, c);

  Eigen::Matrix4d transformation (Eigen::Matrix4d::Identity ());
  if (transf_vec.size () == 6)
  {
    transformation.setZero ();
    transformation (0, 3) = transf_vec[0];
    transformation (1, 3) = transf_vec[1];
    transformation (2, 3) = transf_vec[2];
    transformation (3, 3) = 1.;

    double angle_x = transf_vec[3],
        angle_y = transf_vec[4],
        angle_z = transf_vec[5];
    double aux_a = cos (angle_x),
        aux_b = sin (angle_x),
        aux_c = cos (angle_y),
        aux_d = sin (angle_y),
        aux_e = cos (angle_z),
        aux_f = sin(angle_z),
        aux_ad = aux_a * aux_d,
        aux_bd = aux_b * aux_d;

    transformation (0, 0) = aux_c * aux_e;
    transformation (0, 1) = -aux_c * aux_f;
    transformation (0, 2) = -aux_d;
    transformation (1, 0) = -aux_bd * aux_e + aux_a * aux_f;
    transformation (1, 1) = aux_bd * aux_f + aux_a * aux_e;
    transformation (1, 2) = -aux_b * aux_c;
    transformation (2, 0) = aux_ad * aux_e + aux_b * aux_f;
    transformation (2, 1) = -aux_ad * aux_f + aux_b * aux_e;
    transformation (2, 2) = aux_a * aux_c;

    std::cout << "using transformation:\n" << transformation << std::endl;
  }


  /// Sample the superquadric to a pcd file
  if (output_file != "")
  {
    PointCloud<PointXYZ> cloud;
    for (double eta = -M_PI / 2.0; eta < M_PI / 2.; eta += 0.05)
    {
      for (double mu = -M_PI; mu < M_PI; mu += 0.05)
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


//        double term_1 = pow (fabs(p.x/a), 2./epsilon_2);
//        double term_2 = pow (fabs(p.y/b), 2./epsilon_2);
//        double term_3 = pow (fabs(p.z/c), 2./epsilon_1);
//        double superellipsoid_f = pow (fabs(term_1 + term_2), epsilon_2/epsilon_1) + term_3;

//        PCL_INFO ("f %f\n", superellipsoid_f);

//        double e1 = epsilon_1,
//            e2 = epsilon_2,
//            tx = transf_vec[0],
//            ty = transf_vec[1],
//            tz = transf_vec[2],
//            ax = transf_vec[3],
//            ay = transf_vec[4],
//            az = transf_vec[5],
//            x = p[0],
//            y = p[1],
//            z = p[2];

        // Formula from Maple
//        double val = pow(pow((cos(ay) * cos(az) * x - cos(ay) * sin(az) * y + sin(ay) * z + tx) / a, (double) (2 / e2)) + pow(((-sin(ax) * sin(ay) * cos(az) + cos(ax) * sin(az)) * x + (sin(ax) * sin(ay) * sin(az) + cos(ax) * cos(az)) * y - sin(ax) * cos(ay) * z + ty) / b, (double) (2 / e2)), (double) (e2 / e1)) + pow(((cos(ax) * sin(ay) * cos(az) + sin(ax) * sin(az)) * x + (-cos(ax) * sin(ay) * sin(az) + sin(ax) * cos(az)) * y + cos(ax) * cos(ay) * z + tz) / c, (double) (2 / e1));
//        PCL_INFO ("val = %f\n", val);
      }
    }

    cloud.height = 1.;
    cloud.width = cloud.size ();
    io::savePCDFile (output_file, cloud);
  }



  return EXIT_SUCCESS;
}
