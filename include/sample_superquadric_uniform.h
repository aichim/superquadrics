#pragma once


#include "sample_superquadric.h"


namespace sq
{
  template<typename T>
  struct SuperquadricParameters;

  /**
   * Smart superquadric sampling, should result in almost uniform sampling in 3D space
   */
  template <typename PointT, typename Scalar>
  class SuperquadricSamplingUniform
  {
  public:
    SuperquadricSamplingUniform ();

    inline void
    setParameters (const SuperquadricParameters<Scalar> &params)
    { params_ = params; }

    inline void
    setSpatialSampling (const double D)
    { D_ = D; }

    void
    generatePointCloud (pcl::PointCloud<PointT> &cloud);

    void
    generateMesh (pcl::PolygonMesh &mesh);

    double
    updateTheta (const double a1,
                 const double a2,
                 const double eps,
                 const double theta);


  protected:
    SuperquadricParameters<Scalar> params_;
    double D_;
  };
}


#include "impl/sample_superquadric_uniform.hpp"
