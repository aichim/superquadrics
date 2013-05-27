#include <pcl/io/pcd_io.h>
#include <pcl/common/pca.h>

#include "fit_superquadric_ceres.h"
#include "superquadric_formulas.h"


using namespace pcl;

const double cluster_sq_error = 6000.;
const size_t min_cluster_points_ = 20;


struct SuperquadricCluster
{
  PointCloud<PointXYZ>::Ptr cloud;
  sq::SuperquadricParameters<double> params;
};


double
fitBestSuperquadric (PointCloud<PointXYZ>::ConstPtr cloud,
                     sq::SuperquadricParameters<double> &min_params)
{
  sq::SuperquadricFittingCeres<PointXYZ, double> sq_fit;
  sq_fit.setInputCloud (cloud);

  double min_fit = std::numeric_limits<double>::max ();
  for (int i = 0; i < 3; ++i)
  {
    sq::SuperquadricParameters<double> params;
    sq_fit.setPreAlign (true, i);
    double fit = sq_fit.fit (params);
    printf ("pre_align axis %d, fit %f\n", i, fit);

    if (fit < min_fit)
    {
      min_fit = fit;
      min_params = params;
    }
  }

  return (min_fit);
}

void
splitCluster (PointCloud<PointXYZ>::ConstPtr cloud,
              PointCloud<PointXYZ> &cluster_one,
              PointCloud<PointXYZ> &cluster_two,
              Eigen::Vector4d &splitting_plane_params)
{
  /// Compute PCA for the input cloud
  PCA<PointXYZ> pca;
  pca.setInputCloud (cloud);
  Eigen::Vector3f eigenvalues = pca.getEigenValues ();
  Eigen::Matrix3f eigenvectors = pca.getEigenVectors ();

  PCL_INFO ("eigenvalues: %f %f %f\n", eigenvalues[0]/eigenvalues[0],
      eigenvalues[1] / eigenvalues[0], eigenvalues[2]/eigenvalues[0]);
  Eigen::Vector3d main_pca_axis = eigenvectors.col (0).cast<double> ();

  /// Compute the centroid
  Eigen::Vector4d centroid;
  pcl::compute3DCentroid (*cloud, centroid);

  /// The plane equation
  double d = (-1) * (centroid[0] * main_pca_axis[0] +
                     centroid[1] * main_pca_axis[1] +
                     centroid[2] * main_pca_axis[2]);


  /// Now divide the input cloud into two clusters based on the splitting plane
  for (size_t i = 0; i < cloud->size (); ++i)
  {
    const Eigen::Vector3d point = (*cloud)[i].getVector3fMap ().cast<double> ();
    if (point.dot (main_pca_axis) + d < 0.f)
      cluster_one.push_back ((*cloud)[i]);
    else
      cluster_two.push_back ((*cloud)[i]);
  }


  /// Set the splitting plane parameters
  splitting_plane_params[0] = main_pca_axis[0];
  splitting_plane_params[1] = main_pca_axis[1];
  splitting_plane_params[2] = main_pca_axis[2];
  splitting_plane_params[3] = d;
}


void
splitClusters (PointCloud<PointXYZ>::ConstPtr cloud_in,
               std::vector<SuperquadricCluster> &clusters_result,
               std::vector<Eigen::Vector4d> &splitting_planes)
{
  std::deque<PointCloud<PointXYZ>::Ptr> cluster_queue;
  cluster_queue.push_front (PointCloud<PointXYZ>::Ptr (new PointCloud<PointXYZ> (*cloud_in)));

  while (cluster_queue.size () != 0)
  {
    /// Get the last cluster
    PointCloud<PointXYZ>::Ptr cloud_current = cluster_queue.back ();
    cluster_queue.pop_back ();

    /// Do not try to fit a superquadric if there are not enough points
    if (cloud_current->size () < min_cluster_points_)
    {
      SuperquadricCluster cluster_done;
      cluster_done.cloud = cloud_current;
      clusters_result.push_back (cluster_done);
      continue;
    }

    PCL_INFO ("Fitting a superquadric to a cloud with %ld points.\n", cloud_current->size ());
    sq::SuperquadricParameters<double> params;
    double fit_error = fitBestSuperquadric (cloud_current, params);

    PCL_INFO ("Superquadric fit error = %f\n", fit_error);

    if (fit_error > cluster_sq_error)
    {
      PCL_INFO ("Splitting ...\n");
      PointCloud<PointXYZ>::Ptr cluster_one (new PointCloud<PointXYZ> ());
      PointCloud<PointXYZ>::Ptr cluster_two (new PointCloud<PointXYZ> ());
      Eigen::Vector4d splitting_plane_params;
      splitCluster (cloud_current,
                    *cluster_one, *cluster_two,
                    splitting_plane_params);

      cluster_queue.push_front (cluster_one);
      cluster_queue.push_front (cluster_two);

      splitting_planes.push_back (splitting_plane_params);
    }
    else
    {
      SuperquadricCluster cluster_done;
      cluster_done.cloud = cloud_current;
      cluster_done.params = params;
      clusters_result.push_back (cluster_done);
    }
  }

  PCL_INFO ("Number of clusters after splitting: %ld\n   with sizes: ", clusters_result.size ());
  for (size_t i = 0; i < clusters_result.size (); ++i)
    PCL_INFO ("%ld ", clusters_result[i].cloud->size ());
  PCL_INFO ("\n");
}



void
mergeClusters (std::vector<SuperquadricCluster> &clusters_in,
               std::vector<SuperquadricCluster> &clusters_merged,
               std::vector<Eigen::Vector4d> &splitting_planes)
{
  /// TODO would need to carry them over from the splitting step
  /// Compute the centroids in order to find the neighbors
  std::vector<Eigen::Vector3d> centroids;
  for (size_t cluster_i = 0; cluster_i < clusters_in.size (); ++cluster_i)
  {
    Eigen::Vector4d centroid;
    compute3DCentroid (*clusters_in[cluster_i].cloud, centroid);
    centroids.push_back (Eigen::Vector3d (centroid[0], centroid[1], centroid[2]));
  }


  int merge_count = 1;
  while (merge_count != 0)
  {
    merge_count = 0;

    std::vector<std::pair<size_t, size_t> > possible_merges;
    std::vector<double> possible_errors;
    /// For each cluster pair
    for (size_t c_i = 0; c_i < clusters_in.size (); ++c_i)
    {
      for (size_t c_j = c_i + 1; c_j < clusters_in.size (); ++c_j)
      {
        int planes_intersected = 0;
        for (size_t p_i = 0; p_i < splitting_planes.size (); ++p_i)
        {
          bool side_i = (centroids[c_i][0] * splitting_planes[p_i][0] +
                        centroids[c_i][1] * splitting_planes[p_i][1] +
                        centroids[c_i][2] * splitting_planes[p_i][2] +
                        splitting_planes[p_i][3] > 0.);
          bool side_j = (centroids[c_i][0] * splitting_planes[p_i][0] +
                        centroids[c_i][1] * splitting_planes[p_i][1] +
                        centroids[c_i][2] * splitting_planes[p_i][2] +
                        splitting_planes[p_i][3] > 0.);

          if (side_i ^ side_j)
            planes_intersected ++;
        }

        if (planes_intersected <= 2)
        {
          /// Compute the union of the clusters
          PointCloud<PointXYZ>::Ptr clusters_union (new PointCloud<PointXYZ> ());
          *clusters_union = *(clusters_in[c_i].cloud) + *(clusters_in[c_j].cloud);

          /// Center the cloud
          Eigen::Vector4d centroid;
          pcl::compute3DCentroid (*clusters_union, centroid);
          PointCloud<PointXYZ>::Ptr cloud_centered (new PointCloud<PointXYZ> ());
          demeanPointCloud (*clusters_union, centroid, *cloud_centered);

          /// Fit a superquadric in the union
          sq::SuperquadricFittingCeres<PointXYZ, double> sq_fit;
          sq_fit.setInputCloud (cloud_centered);
          sq::SuperquadricParameters<double> params;
          double fit_error = sq_fit.fit (params);
          PCL_INFO ("Superquadric fit error = %f\n", fit_error);

          if (fit_error < cluster_sq_error)
          {
            possible_merges.push_back (std::make_pair (c_i, c_j));
            possible_errors.push_back (fit_error);
          }
        }
      }
    }

    PCL_INFO ("Possible merges: ");
    for (size_t i = 0; i < possible_merges.size (); ++i)
      PCL_INFO ("(%ld (%ld points), %ld (%ld points) - %f)\n",
                possible_merges[i].first, clusters_in[possible_merges[i].first].cloud->size (),
                possible_merges[i].second, clusters_in[possible_merges[i].second].cloud->size (),
                possible_errors[i]);
  }


}


int
main (int argc,
      char **argv)
{
  google::InitGoogleLogging (argv[0]);

  PointCloud<PointXYZ>::Ptr cloud (new PointCloud<PointXYZ> ());
  io::loadPCDFile (argv[1], *cloud);

  PCL_INFO ("####### SPLITTING ######\n");

  std::vector<SuperquadricCluster> clusters_split;
  std::vector<Eigen::Vector4d> splitting_planes;
  splitClusters (cloud, clusters_split, splitting_planes);

  PCL_INFO ("Clusters after splitting:\n");
  for (size_t i = 0; i < clusters_split.size (); ++i)
  {
    PCL_INFO ("Command for sampler:\n-e1 %f -e2 %f -a %f -b %f -c %f -transform %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
              clusters_split[i].params.e1, clusters_split[i].params.e2, clusters_split[i].params.a, clusters_split[i].params.b, clusters_split[i].params.c,
              clusters_split[i].params.transform (0, 0), clusters_split[i].params.transform (0, 1), clusters_split[i].params.transform (0, 2), clusters_split[i].params.transform (0, 3),
              clusters_split[i].params.transform (1, 0), clusters_split[i].params.transform (1, 1), clusters_split[i].params.transform (1, 2), clusters_split[i].params.transform (1, 3),
              clusters_split[i].params.transform (2, 0), clusters_split[i].params.transform (2, 1), clusters_split[i].params.transform (2, 2), clusters_split[i].params.transform (2, 3),
              clusters_split[i].params.transform (3, 0), clusters_split[i].params.transform (3, 1), clusters_split[i].params.transform (3, 2), clusters_split[i].params.transform (3, 3));
    PCL_INFO ("\n");
  }

  PCL_INFO ("####### MERGING ######\n");

  std::vector<SuperquadricCluster> clusters_merge;
  mergeClusters (clusters_split, clusters_merge, splitting_planes);


  return (0);
}
