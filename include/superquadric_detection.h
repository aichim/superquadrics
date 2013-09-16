#pragma once

#include <boost/random.hpp>

#include <fit_superquadric_ceres.h>


namespace sq
{

template<typename PointT, typename Scalar>
class SuperquadricDetection
{
public:
  typedef pcl::PointCloud<PointT> Cloud;
  typedef typename Cloud::Ptr CloudPtr;
  typedef typename Cloud::ConstPtr CloudConstPtr;

  struct SuperquadricDetectionHypothesis
  {
    double fit_error;
    int num_interior_points;
    int num_surface_points;
    Eigen::Matrix<Scalar, 4, 4> transform;
  };


  SuperquadricDetection ();

  void
  setInputCloud (CloudConstPtr cloud_input)
  { cloud_input_ = cloud_input; }

  void
  setSuperquadricParams (SuperquadricParameters<Scalar> &params)
  { superquadric_params_ = params; }

  void
  process (typename std::vector<SuperquadricDetectionHypothesis> &hypotheses);

protected:
  CloudConstPtr cloud_input_;
  boost::random::mt19937 rng_;


  int ransac_starts_;
  int ransac_hypotheses_;
  int ransac_model_points_;
  double gamma_;
  SuperquadricParameters<Scalar> superquadric_params_;
};
}


#include "impl/superquadric_detection.hpp"
