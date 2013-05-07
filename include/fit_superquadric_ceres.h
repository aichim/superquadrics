#pragma once

#include <pcl/point_cloud.h>
#include <pcl/common/common.h>

#include <ceres/ceres.h>


template<typename Scalar>
struct SuperquadricParams
{
  SuperquadricParams ()
    : e1 (0.), e2 (0.), a (1.), b (1.), c (1.)
    , transform (Eigen::Matrix<Scalar, 4, 4>::Identity ())
  {}

  Scalar e1, e2, a, b, c;
  Eigen::Matrix<Scalar, 4, 4> transform;
};


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

//  void
//  setIndices (const pcl::IndicesConstPtr &indices)
//  { indices_ = indices; }

  void
  preAlign (Eigen::Matrix<MatScalar, 4, 4> &transformation_prealign,
            Eigen::Matrix<MatScalar, 3, 1> &variances);



  double
  fit (SuperquadricParams<MatScalar> &parameters);


  struct SuperquadricCostFunctor
  {
    template <typename T> bool
    operator () (const T* const x, T* residual) const
    {
      residual[0] = T (10.0) - x[0];
      return (true);
    }
  };


protected:
  CloudConstPtr input_;
  CloudPtr input_prealigned_;

  bool pre_align_;
  int pre_align_axis_;
};


#include "impl/fit_superquadric_ceres.hpp"
