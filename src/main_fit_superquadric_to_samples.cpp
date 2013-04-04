#include <pcl/io/pcd_io.h>

#include "fit_superquadric_lm.h"

using namespace pcl;

int
main (int argc,
      char **argv)
{
  PointCloud<PointXYZ>::Ptr cloud (new PointCloud<PointXYZ> ());
  io::loadPCDFile (argv[1], *cloud);

//  for (size_t i = 0; i < cloud->size (); ++i)
//    (*cloud)[i].getVector3fMap () /= 200.;

  SuperquadricFittingLM<PointXYZ, double> sq_fit;
  sq_fit.setInputCloud (cloud);
  Eigen::Matrix<double, Eigen::Dynamic, 1> params (11);
  Eigen::Matrix4d transformation;
  sq_fit.fit (params, transformation);


  params[8] = params[8] - floor (params[8] / M_PI) * M_PI;
  params[9] = params[9] - floor (params[9] / M_PI) * M_PI;
  params[10] = params[10] - floor (params[10] / M_PI) * M_PI;
  std::cout << "optimized parameters: " << params << std::endl;



  return (0);
}
