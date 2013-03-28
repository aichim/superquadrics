#include <pcl/io/pcd_io.h>

#include <pcl/console/print.h>
#include <pcl/console/parse.h>

using namespace pcl;


double signum (const double &x)
{
  return ((x > 0) - (x < 0));
}


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

  PCL_INFO ("### Superquadric parameters:\n   epsilon_1: %f\n   epsilon_2: %f\n   a: %f\n   b: %f\n   c: %f\n",
            epsilon_1, epsilon_2, a, b, c);


  /// Sample the superquadric to a pcd file
  if (output_file != "")
  {
    PointCloud<PointXYZ> cloud;
    for (double eta = -M_PI / 2.0; eta < M_PI / 2.; eta += 0.1)
    {
      for (double mu = -M_PI; mu < M_PI; mu += 0.1)
      {
        PointXYZ p;
        p.x = a * (cos(eta) / fabs(cos(eta))) * pow (fabs(cos (eta)), epsilon_1) * (cos(mu) / fabs(cos(mu))) * pow (fabs (cos (mu)), epsilon_2);
        p.y = b * (cos(eta) / fabs(cos(eta))) * pow (fabs(cos (eta)), epsilon_1) * (sin(mu) / fabs(sin(mu))) * pow (fabs (sin (mu)), epsilon_2);
        p.z = c * (sin(eta) / fabs(sin(eta))) * pow (fabs(sin (eta)), epsilon_1);

        if (isFinite (p))
          cloud.push_back (p);


        double term_1 = pow (fabs(p.x/a), 2./epsilon_2);
        double term_2 = pow (fabs(p.y/b), 2./epsilon_2);
        double term_3 = pow (fabs(p.z/c), 2./epsilon_1);
        double superellipsoid_f = pow (fabs(term_1 + term_2), epsilon_2/epsilon_1) + term_3;

        PCL_INFO ("f %f\n", superellipsoid_f);

      }
    }

    cloud.height = 1.;
    cloud.width = cloud.size ();
    io::savePCDFile (output_file, cloud);
  }



  return EXIT_SUCCESS;
}
