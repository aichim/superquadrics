#pragma once

#include <pcl/point_cloud.h>
#include <pcl/common/common.h>



namespace sq
{
template <typename T>
struct SuperquadricParameters;

template <typename PointT, typename MatScalar = double>
class SuperquadricRigidRegistration
{
  typedef pcl::PointCloud<PointT> Cloud;
  typedef typename Cloud::Ptr CloudPtr;
  typedef typename Cloud::ConstPtr CloudConstPtr;


public:
  SuperquadricRigidRegistration ();

  void
  setInputCloud (const CloudConstPtr &cloud)
  { input_ = cloud; }

  void
  setParameters (MatScalar e1, MatScalar e2, MatScalar a, MatScalar b, MatScalar c)
  {
    e1_ = e1; e2_ = e2;
    a_ = a; b_ = b; c_ = c;
  }

  void
  setInitTransform (Eigen::Matrix<MatScalar, 4, 4> &tr)
  { init_transform_ = tr; }

  void
  setPreAlign (bool pre_align,
               int pre_align_axis = 2)
  {
    pre_align_ = pre_align;
    pre_align_axis_ = pre_align_axis;
  }

  void
  preAlign (Eigen::Matrix<MatScalar, 4, 4> &transformation_prealign);


  double
  fit (Eigen::Matrix<MatScalar, 4, 4> &transform);


  struct SuperquadricCostFunctor
  {
    SuperquadricCostFunctor (const PointT &point,
                             MatScalar e1, MatScalar e2, MatScalar a, MatScalar b, MatScalar c)
    {
      point_ = point;
      e1_ = e1; e2_ = e2;
      a_ = a; b_ = b; c_ = c;
    }

    template <typename T> bool
    operator () (const T* const x, T* residual) const;

    PointT point_;
    MatScalar e1_, e2_, a_, b_, c_;
  };


protected:
  CloudConstPtr input_;
  CloudPtr input_prealigned_;

  MatScalar e1_, e2_, a_, b_, c_;
  Eigen::Matrix<MatScalar, 4, 4> init_transform_;

  bool pre_align_;
  int pre_align_axis_;
};
}

#include "impl/superquadric_rigid_registration.hpp"
