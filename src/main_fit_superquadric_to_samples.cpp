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

  SuperquadricFittingLM<PointXYZ, double> sq_fit;
  sq_fit.setInputCloud (cloud_in);
  Eigen::Matrix<double, Eigen::Dynamic, 1> params (11);
  Eigen::Matrix4d transformation;
  sq_fit.fit (params, transformation);


//  params[8] = params[8] - floor (params[8] / M_PI) * M_PI;
//  params[9] = params[9] - floor (params[9] / M_PI) * M_PI;
//  params[10] = params[10] - floor (params[10] / M_PI) * M_PI;
  std::cout << "optimized parameters: " << params << std::endl;


//  PCL_INFO ("Command for sampler:\n-e1 %f -e2 %f -a %f -b %f -c %f -transform %f,%f,%f,%f,%f,%f\n",
//            params[0], params[1], params[2], params[3], params[4], params[5], params[6], params[7], params[8], params[9], params[10]);


  PCL_INFO ("Command for sampler:\n-e1 %f -e2 %f -a %f -b %f -c %f -transform %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
            params[0], params[1], params[2], params[3], params[4],
            transformation (0, 0), transformation (0, 1), transformation (0, 2), transformation (0, 3),
            transformation (1, 0), transformation (1, 1), transformation (1, 2), transformation (1, 3),
            transformation (2, 0), transformation (2, 1), transformation (2, 2), transformation (2, 3),
            transformation (3, 0), transformation (3, 1), transformation (3, 2), transformation (3, 3));




  return (0);
}
