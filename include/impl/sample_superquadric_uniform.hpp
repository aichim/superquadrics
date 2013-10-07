#pragma once

#include <pcl/point_cloud.h>
#include <pcl/PolygonMesh.h>

#include "sample_superquadric_uniform.h"
#include "superquadric_formulas.h"

template <typename PointT, typename Scalar>
sq::SuperquadricSamplingUniform<PointT, Scalar>::SuperquadricSamplingUniform ()
  : D_ (1)
{
}


template <typename PointT, typename Scalar> double
sq::SuperquadricSamplingUniform<PointT, Scalar>::updateTheta (const double a1,
                                                              const double a2,
                                                              const double eps,
                                                              const double theta)
{
  const double thresh = 1e-2;
  if (theta > thresh)
//      theta < M_PI/2. - 1e-2)
  {
    PCL_INFO ("first: %f, inside sqrt %f\n",
              D_ / eps * cos (theta) * sin (theta),
              1./ (a1*a1 * pow (cos(theta), 2*eps) * pow(sin(theta), 4.) + a2*a2 * pow (sin(theta), 2*eps) * pow(cos(theta), 4.)));

    return (D_ / eps * cos (theta) * sin (theta) * sqrt (1./ (a1*a1 * pow (cos(theta), 2*eps) * pow(sin(theta), 4.) +
                                                              a2*a2 * pow (sin(theta), 2*eps) * pow(cos(theta), 4.))));
  }
  else if (theta <= thresh && theta > 0.)
  {
    PCL_INFO ("special option: theta %f, theta^eps %f,   %f\n",
              theta, pow(theta, eps),
              D_ / a2  - pow (theta, eps));
    return (pow (fabs (D_ / a2  - pow (theta, eps)), 1. / eps) - theta);
  }
  else
  {
    return (pow (fabs (D_ / a2), 1. / eps));
  }
//  else //if (theta > M_PI/2. - 1e-2)
//  {
//    return (pow (fabs (D_ / a1  + pow (M_PI/2. - theta, eps)), 1. / eps) - (M_PI/2. - theta));
//  }
}



template <typename PointT, typename Scalar> void
sq::SuperquadricSamplingUniform<PointT, Scalar>::generatePointCloud (pcl::PointCloud<PointT> &cloud)
{
  const double start_value = 0.; //std::numeric_limits<double>::epsilon ();
  double eta = start_value;
  double omega = start_value;
  bool eta_increasing = true;
  bool omega_increasing = true;

  while (eta > -std::numeric_limits<double>::epsilon ())
  {
    PCL_INFO ("Eta %f, omega %f at sample %zu\n", eta, omega, cloud.size ());
    /// Generate sample
    double x1, x2, y1, y2;
    if (omega > M_PI/2. - 6e-1)
    {
      omega = std::max (M_PI/2. - omega, 0.);
      omega_increasing = false;
    }

    if (omega_increasing)
    {
      x1 = pow (cos (omega), params_.e1);
      x2 = params_.c * pow (sin (omega), params_.e1);
    }
    else
    {
      /// Swap them
      x1 = pow (sin (omega), params_.e1);
      x2 = params_.c * pow (cos (omega), params_.e1);
    }

    /// Now update the parameters
    double update_omega = std::max (updateTheta (params_.a, params_.b, params_.e1, omega), 1e-4);
    PCL_INFO ("Update omega: %f\n", update_omega);
    if (omega_increasing)
      omega += update_omega;
    else
      omega -= update_omega;


    if (eta_increasing)
    {
      y1 = params_.a * pow (cos (eta), params_.e2);
      y2 = params_.b * pow (sin (eta), params_.e2);
    }
    else
    {
      y1 = params_.a * pow (sin (eta), params_.e2);
      y2 = params_.b * pow (cos (eta), params_.e2);
    }

    if (omega < 0
        && !omega_increasing)
    {
      omega = start_value;
      omega_increasing = true;

      if (eta > M_PI/2. - 6e-1)
      {
        eta = std::max (M_PI/2. - eta, 0.);
        eta_increasing = false;
      }


      double update_eta = std::max (updateTheta (1., params_.c, params_.e2, eta), 1e-4);
      PCL_INFO ("update eta: %f\n", update_eta);
      if (eta_increasing)
        eta += update_eta;
      else
      {
        eta -= update_eta;
        PCL_INFO ("Eta decreasing new eta %f, update_eta %f\n", eta, update_eta);
      }
    }




    for (int k = 0; k < 8; ++k)
    {
      Eigen::Matrix<Scalar, 4, 1> p;
      p[0] = x1 * y1 * pow (-1., k / 4);
      p[1] = x1 * y2 * pow (-1., (k % 4) / 2);
      p[2] = x2 * pow (-1, k % 2);
      p[3] = 1.;

      PointT point;
      Eigen::Matrix<Scalar, 4, 1> p_tr = params_.transform.inverse () * p;
      point.x = p_tr[0];
      point.y = p_tr[1];
      point.z = p_tr[2];
      if (isFinite (point))
        cloud.push_back (point);
    }


//    PCL_INFO ("Update eta: %f, update omega %f\n", update_eta, update_omega);
  }

  cloud.height = 1.;
  cloud.width = cloud.size ();
}



template <typename PointT, typename Scalar> void
sq::SuperquadricSamplingUniform<PointT, Scalar>::generateMesh (pcl::PolygonMesh &mesh)
{
  /*  /// First generate the point cloud
  pcl::PointCloud<PointT> cloud;
  generatePointCloud (cloud);


  /// Then look into the connectivity
  for (size_t v = 1; v < eta_samples_; ++v)
  {
    for (size_t u = 1; u < mu_samples_; ++u)
    {
      pcl::Vertices polygon;
      polygon.vertices.push_back (v * eta_samples_ + u);
      polygon.vertices.push_back (v * eta_samples_ + u - 1);
      polygon.vertices.push_back ((v - 1) * eta_samples_ + u);
      mesh.polygons.push_back (polygon);

      polygon.vertices.clear ();
      polygon.vertices.push_back ((v - 1) * eta_samples_ + u);
      polygon.vertices.push_back (v * eta_samples_ + u - 1);
      polygon.vertices.push_back ((v - 1) * eta_samples_ + u - 1);
      mesh.polygons.push_back (polygon);
    }

    /// And connect the last column with the first one
    pcl::Vertices polygon;
    polygon.vertices.push_back (v * eta_samples_ + 0);
    polygon.vertices.push_back (v * eta_samples_ + eta_samples_ - 1);
    polygon.vertices.push_back ((v - 1) * eta_samples_ + 0);
    mesh.polygons.push_back (polygon);

    polygon.vertices.clear ();
    polygon.vertices.push_back ((v - 1) * eta_samples_ + 0);
    polygon.vertices.push_back (v * eta_samples_ + eta_samples_ - 1);
    polygon.vertices.push_back ((v - 1) * eta_samples_ + eta_samples_ - 1);
    mesh.polygons.push_back (polygon);
  }
*/
  /// TODO
  /*
  /// Take care of the top: add a new point as the centroid of the first row
  /// and do a triangle strip there
  Eigen::Vector3f centroid (0.f, 0.f, 0.f);
  for (size_t i = 0; i < eta_samples_; ++i)
    centroid += cloud[i].getVector3fMap ();
  centroid /= mu_samples_;
  PointT point_centroid;
  point_centroid.getVector3fMap () = centroid;
  cloud.push_back (point_centroid);

  for (size_t i = 1; i < eta_samples_; ++i)
  {
    pcl::Vertices polygon;
    polygon.vertices.push_back (cloud.size () - 1);
    polygon.vertices.push_back (i - 1);
    polygon.vertices.push_back (i);

    mesh.polygons.push_back (polygon);
  }
  pcl::Vertices polygon;
  polygon.vertices.push_back (cloud.size () - 1);
  polygon.vertices.push_back (mu_samples_);
  polygon.vertices.push_back (0);
  mesh.polygons.push_back (polygon);
*/
  //  toROSMsg (cloud, mesh.cloud);
}
