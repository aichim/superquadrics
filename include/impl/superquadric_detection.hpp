#pragma once

#include "superquadric_detection.h"
#include "superquadric_rigid_registration.h"
#include <pcl/search/kdtree.h>


template<typename PointT, typename Scalar>
sq::SuperquadricDetection<PointT, Scalar>::SuperquadricDetection ()
  : ransac_starts_ (10)
  , ransac_hypotheses_ (20) //(12)
  , ransac_model_points_ (30)
  , gamma_ (0.2)
{
  rng_.seed (std::time (NULL));
}



template<typename PointT, typename Scalar> void
sq::SuperquadricDetection<PointT, Scalar>::process (typename std::vector<SuperquadricDetectionHypothesis> &hypotheses)
{
  /// Set up a search tree
  pcl::search::KdTree<PointT> tree;
  tree.setInputCloud (cloud_input_);

  double min_radius = std::min (superquadric_params_.a,
                                std::min (superquadric_params_.b, superquadric_params_.c));
  double max_radius = 1.2 * std::max (superquadric_params_.a,
                                      std::max (superquadric_params_.b, superquadric_params_.c));

  boost::random::uniform_int_distribution<> random_index_input (0, cloud_input_->size () - 1);


  /// Run ransac_hypotheses_ instances of RANSAC, i.e., ransac_hypotheses_ random seed points
  for (size_t hypothesis_i = 0; hypothesis_i < ransac_hypotheses_; ++hypothesis_i)
  {
    /// Choose a seed point at random
    int seed_index = random_index_input (rng_);

    /// Find the min_radius neighborhood around this point
    std::vector<int> nn_indices;
    std::vector<float> nn_dists;
    tree.radiusSearch (seed_index, min_radius, nn_indices, nn_dists);

    boost::random::uniform_int_distribution<> random_index_model (0, nn_indices.size () - 1);


    double min_fit_error = std::numeric_limits<double>::max ();
    Eigen::Matrix<Scalar, 4, 4> min_transf = Eigen::Matrix<Scalar, 4, 4>::Identity ();

    /// For each seed point, choose ransac_model_points_ around the seed point ransac_starts_ times
    for (size_t ransac_run = 0; ransac_run < ransac_starts_; ++ransac_run)
    {
      CloudPtr cloud (new Cloud ());
      cloud->reserve (ransac_model_points_);
      /// Choose ransac_model_points_ at random from the neighborhood
      for (size_t i = 0; i < ransac_model_points_; ++i)
      {
        int index = nn_indices[random_index_model (rng_)];
        cloud->push_back ((*cloud_input_)[index]);
      }
      cloud->width = cloud->size ();
      cloud->height = 1;

      /// Fit a superquadric
      SuperquadricRigidRegistration<PointT, Scalar> sq_reg;
      sq_reg.setInputCloud (cloud);
      sq_reg.setParameters (superquadric_params_.e1, superquadric_params_.e2, superquadric_params_.a, superquadric_params_.b, superquadric_params_.c);

      for (int i = 0; i < 3; ++i)
      {
        sq_reg.setPreAlign (true, i);
        Eigen::Matrix4d transf;
        double fit = sq_reg.fit (transf);
        printf ("pre_align axis %d, fit %f\n", i, fit);

        if (fit < min_fit_error)
        {
          min_fit_error = fit;
          min_transf = transf;
        }
      }
    }


    /// Now fit the superquadric using the full resolution neighborhood
    CloudPtr cloud (new Cloud ());
    cloud->reserve (nn_indices.size ());
    for (size_t i = 0; i < nn_indices.size (); ++i)
      cloud->push_back ((*cloud_input_)[nn_indices[i]]);

    SuperquadricRigidRegistration<PointT, Scalar> sq_reg;
    sq_reg.setInputCloud (cloud);
    sq_reg.setPreAlign (false);
    sq_reg.setParameters (superquadric_params_.e1, superquadric_params_.e2, superquadric_params_.a, superquadric_params_.b, superquadric_params_.c);
    sq_reg.setInitTransform (min_transf);
    Eigen::Matrix<Scalar, 4, 4> transf;
    double fitting_error = sq_reg.fit (transf);
    PCL_INFO ("Fitting error on full res neighborhood: %f\n", fitting_error);

    /// TODO should set this to a reasonable value or return false from fitting class
//    if (fitting_error > 10000)
//      continue;


    /// Compute the number of interior and surface points
    std::vector<int> nn_large_indices;
    std::vector<float> nn_large_dists;
    PointT p;
    Eigen::Matrix4d tr = transf.template cast<double> ().inverse ();
    p.x = tr (0, 3); p.y = tr (1, 3); p.z = tr (2, 3);
    tree.radiusSearch (p, max_radius, nn_large_indices, nn_large_dists);
//    tree.radiusSearch (p, min_radius, nn_large_indices, nn_large_dists);

    PCL_INFO ("Checking out large neighborhood of radius %f\n", max_radius);
    PCL_INFO ("Seed center %f %f %f, sq center %f %f %f\n",
              (*cloud_input_)[seed_index].x, (*cloud_input_)[seed_index].y, (*cloud_input_)[seed_index].z,
              p.x, p.y, p.z);

    int num_interior_points = 0;
    int num_surface_points = 0;
    for (size_t i = 0; i < nn_large_indices.size (); ++i)
    {
      Eigen::Vector4d point ((*cloud_input_)[nn_large_indices[i]].x,
                             (*cloud_input_)[nn_large_indices[i]].y,
                             (*cloud_input_)[nn_large_indices[i]].z,
                             1.);
      point = transf * point;

      double value = superquadric_function<Scalar> (point[0], point[1], point[2],
                                                    superquadric_params_.e1, superquadric_params_.e2, superquadric_params_.a, superquadric_params_.b, superquadric_params_.c);

      printf ("value = %f\n", value);

      if (value < - gamma_)
          num_interior_points ++;
      if (value > - gamma_ &&
          value < + gamma_)
        num_surface_points ++;
    }

    SuperquadricDetectionHypothesis hypo;
    hypo.transform = transf;
    hypo.fit_error = fitting_error;
    hypo.num_interior_points = num_interior_points;
    hypo.num_surface_points = num_surface_points;
    hypotheses.push_back (hypo);

    PCL_INFO ("Generated hypothesis with:\n   fitting_error: %f\n   #interior: %d\n   #surface: %d\n",
              hypo.fit_error, hypo.num_interior_points, hypo.num_surface_points);
  }
}
