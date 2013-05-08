#include <pcl/io/pcd_io.h>

#include "fit_superquadric_lm.h"

using namespace pcl;

int
main (int argc,
      char **argv)
{
  PointCloud<PointXYZ>::Ptr cloud_in (new PointCloud<PointXYZ> ());
  io::loadPCDFile (argv[1], *cloud_in);

//  Eigen::Vector4d centroid;
//  pcl::compute3DCentroid (*cloud_in, centroid);
//  PointCloud<PointXYZ>::Ptr cloud_centered (new PointCloud<PointXYZ> ());
//  demeanPointCloud (*cloud_in, centroid, *cloud_centered);

//  io::savePCDFile ("cloud_input.pcd", *cloud_centered, true);

  sq::SuperquadricFittingLM<PointXYZ, double> sq_fit;
  sq_fit.setInputCloud (cloud_in);

  double min_fit = std::numeric_limits<double>::max ();
  sq::SuperquadricParameters<double> min_params;
  for (int i = 0; i < 3; ++i)
  {
    sq::SuperquadricParameters<double> params;
    sq_fit.setPreAlign (true, i);
    double fit = sq_fit.fit (params);
    printf ("pre_align axis %d, fit %f\n", i, fit);

    if (fit < min_fit)
    {
      min_fit = fit;
      min_params = params;
    }
  }


//  params[8] = params[8] - floor (params[8] / M_PI) * M_PI;
//  params[9] = params[9] - floor (params[9] / M_PI) * M_PI;
//  params[10] = params[10] - floor (params[10] / M_PI) * M_PI;


//  PCL_INFO ("Command for sampler:\n-e1 %f -e2 %f -a %f -b %f -c %f -transform %f,%f,%f,%f,%f,%f\n",
//            params[0], params[1], params[2], params[3], params[4], params[5], params[6], params[7], params[8], params[9], params[10]);


  PCL_INFO ("Command for sampler:\n-e1 %f -e2 %f -a %f -b %f -c %f -transform %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
            min_params.e1, min_params.e2, min_params.a, min_params.b, min_params.c,
            min_params.transform (0, 0), min_params.transform (0, 1), min_params.transform (0, 2), min_params.transform (0, 3),
            min_params.transform (1, 0), min_params.transform (1, 1), min_params.transform (1, 2), min_params.transform (1, 3),
            min_params.transform (2, 0), min_params.transform (2, 1), min_params.transform (2, 2), min_params.transform (2, 3),
            min_params.transform (3, 0), min_params.transform (3, 1), min_params.transform (3, 2), min_params.transform (3, 3));




  return (0);
}
