#pragma once

#include <pcl/point_cloud.h>
#include <pcl/common/common.h>



namespace sq
{
template <typename T>
struct SuperquadricParameters;

template <typename PointT, typename MatScalar = double>
class SuperquadricFittingCeres
{
  typedef pcl::PointCloud<PointT> Cloud;
  typedef typename Cloud::Ptr CloudPtr;
  typedef typename Cloud::ConstPtr CloudConstPtr;


public:
  SuperquadricFittingCeres ();

  void
  setPreAlign (bool pre_align,
               int pre_align_axis = 2)
  {
    pre_align_ = pre_align;
    pre_align_axis_ = pre_align_axis;
  }

  void
  setInputCloud (const CloudConstPtr &cloud)
  { input_ = cloud; }

  void
  setInitParameters (SuperquadricParameters<MatScalar> &init_params)
  { init_parameters_ = init_params; }

  void
  preAlign (Eigen::Matrix<MatScalar, 4, 4> &transformation_prealign,
            Eigen::Matrix<MatScalar, 3, 1> &variances);



  double
  fit (SuperquadricParameters<MatScalar> &parameters);


  struct SuperquadricCostFunctor
  {
    SuperquadricCostFunctor (const PointT &point)
    { point_ = point; }

    template <typename T> bool
    operator () (const T* const x, T* residual) const;

    PointT point_;
  };


protected:
  CloudConstPtr input_;
  CloudPtr input_prealigned_;

  SuperquadricParameters<MatScalar> init_parameters_;

  bool pre_align_;
  int pre_align_axis_;
};
}

#include "impl/fit_superquadric_ceres.hpp"
