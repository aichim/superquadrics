#pragma once


namespace sq
{

/** \brief Structure containing the parameters used to define a superquadric */
template<typename Scalar>
struct SuperquadricParameters
{
  SuperquadricParameters ()
    : e1 (0.), e2 (0.), a (1.), b (1.), c (1.)
    , transform (Eigen::Matrix<Scalar, 4, 4>::Identity ())
  {}

  Scalar e1, e2, a, b, c;
  Eigen::Matrix<Scalar, 4, 4> transform;
};


template <typename Scalar> inline Scalar
superquadric_function (const Scalar &x,
                       const Scalar &y,
                       const Scalar &z,
                       const Scalar &e1,
                       const Scalar &e2,
                       const Scalar &a,
                       const Scalar &b,
                       const Scalar &c);

template <typename Scalar> inline Scalar
superquadric_function_scale_weighting (const Scalar &x,
                                       const Scalar &y,
                                       const Scalar &z,
                                       const Scalar &e1,
                                       const Scalar &e2,
                                       const Scalar &a,
                                       const Scalar &b,
                                       const Scalar &c);


/** Computed with Maple, used for the Eigen LM, Ceres uses Automatic Differentiation*/
template<typename Scalar> inline void
superquadric_derivative (const Scalar &x, const Scalar &y, const Scalar &z,
                         const Scalar &e1, const Scalar &e2,
                         const Scalar &a, const Scalar &b, const Scalar &c,
                         const Scalar &tx, const Scalar &ty, const Scalar &tz,
                         const Scalar &ax, const Scalar &ay, const Scalar &az,
                         Scalar &dS_de1, Scalar &dS_de2,
                         Scalar &dS_da, Scalar &dS_db, Scalar &dS_dc,
                         Scalar &dS_dtx, Scalar &dS_dty, Scalar &dS_dtz,
                         Scalar &dS_dax, Scalar &dS_day, Scalar &dS_daz);

/** \brief Compute the error. */
template<typename PointT, typename Scalar> Scalar
computeSuperQuadricError (typename pcl::PointCloud<PointT>::ConstPtr cloud,
                          const Scalar &e1, const Scalar &e2,
                          const Scalar &a, const Scalar &b, const Scalar &c,
                          const Eigen::Matrix<Scalar, 4, 4> &transform);


/// Stirling's approximation
template<typename Scalar> Scalar
betaFunction (Scalar x, Scalar y);

template<typename Scalar> Scalar
computeSuperQuadricVolume (const Scalar &e1, const Scalar &e2,
                           const Scalar &a, const Scalar &b, const Scalar &c);
}

#include "impl/superquadric_formulas.hpp"
