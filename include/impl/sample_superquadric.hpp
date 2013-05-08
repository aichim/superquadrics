#pragma once

#include <pcl/point_cloud.h>
#include <pcl/PolygonMesh.h>

#include "sample_superquadric.h"
#include "superquadric_formulas.h"

template <typename PointT, typename Scalar>
sq::SuperquadricSampling<PointT, Scalar>::SuperquadricSampling ()
{
  transform_ = Eigen::Matrix<Scalar, 4, 4>::Identity ();
  eta_samples_ = 100;
  mu_samples_ = 100;
}


template <typename PointT, typename Scalar> void
sq::SuperquadricSampling<PointT, Scalar>::generatePointCloud (pcl::PointCloud<PointT> &cloud)
{
  Scalar eta_sample_rate = static_cast<Scalar> (M_PI) / static_cast<Scalar> (eta_samples_);
  Scalar mu_sample_rate = static_cast<Scalar> (M_PI) / static_cast<Scalar> (mu_samples_);

  for (double eta = -M_PI / 2.0; eta < M_PI / 2.; eta += eta_sample_rate)
  {
    for (double mu = -M_PI; mu < M_PI; mu += mu_sample_rate)
    {
      Eigen::Matrix<Scalar, 4, 1> p;
      p[0] = params_.a * (cos(eta) / fabs(cos(eta))) * pow (fabs(cos (eta)), params_.e1) * (cos(mu) / fabs(cos(mu))) * pow (fabs (cos (mu)), params_.e2);
      p[1] = params_.b * (cos(eta) / fabs(cos(eta))) * pow (fabs(cos (eta)), params_.e1) * (sin(mu) / fabs(sin(mu))) * pow (fabs (sin (mu)), params_.e2);
      p[2] = params_.c * (sin(eta) / fabs(sin(eta))) * pow (fabs(sin (eta)), params_.e1);
      p[3] = 1.;

      PointT point;
      Eigen::Matrix<Scalar, 4, 1> p_tr = transform_.inverse () * p;
      point.x = p_tr[0];
      point.y = p_tr[1];
      point.z = p_tr[2];


      if (isFinite (point))
        cloud.push_back (point);
    }
  }

  cloud.height = 1.;
  cloud.width = cloud.size ();
}



template <typename PointT, typename Scalar> void
sq::SuperquadricSampling<PointT, Scalar>::generateMesh (pcl::PolygonMesh &mesh)
{
  /// First generate the point cloud
  pcl::PointCloud<PointT> cloud;
  generatePointCloud (cloud);

  toROSMsg (cloud, mesh.cloud);

  /// Then look into the connectivity
}
