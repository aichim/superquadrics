#include <pcl/io/pcd_io.h>

#include "fit_superquadric_lm.h"

using namespace pcl;

int
main (int argc,
      char **argv)
{
  PointCloud<PointXYZ>::Ptr cloud (new PointCloud<PointXYZ> ());
  io::loadPCDFile (argv[1], *cloud);

  /*cloud->clear ();
  for (size_t i = 0; i < 50; ++i)
  {
    PointXYZ p;
    p.x = i;
    p.y = i*2;
    p.z = i *3;
    cloud->push_back (p);
  }*/

  SuperquadricFittingLM<PointXYZ, double> sq_fit;
  sq_fit.setInputCloud (cloud);
  Eigen::Matrix<double, Eigen::Dynamic, 1> params (6);
  Eigen::Matrix4d transformation;
  sq_fit.fit (params, transformation);



  std::cout << "optimized parameters: " << params << std::endl;

  return (0);
}
