#include <pcl/io/pcd_io.h>
#include <pcl/common/pca.h>
#include <pcl/io/vtk_lib_io.h>

#include "fit_superquadric_ceres.h"
#include "superquadric_formulas.h"
#include "sample_superquadric.h"

#include <pcl/console/parse.h>


using namespace pcl;

const double cluster_sq_error = 500;
const size_t min_cluster_points_ = 300;


struct SuperquadricCluster
{
  typedef boost::shared_ptr<SuperquadricCluster> Ptr;
  SuperquadricCluster ()
  {
    cloud.reset (new PointCloud<PointXYZ> ());
  }

  SuperquadricCluster (const SuperquadricCluster& copy)
  {
    if (copy.cloud)
    {
      cloud.reset (new PointCloud<PointXYZ> ());
      *cloud = *copy.cloud;
    }
    params = copy.params;
    fit_error = copy.fit_error;
  }

  PointCloud<PointXYZ>::Ptr cloud;
  sq::SuperquadricParameters<double> params;
  double fit_error;
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
               std::vector<SuperquadricCluster::Ptr> &clusters_result,
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
      SuperquadricCluster::Ptr cluster_done (new SuperquadricCluster ());
      cluster_done->cloud = cloud_current;
      cluster_done->fit_error = std::numeric_limits<double>::max ();
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
      SuperquadricCluster::Ptr cluster_done (new SuperquadricCluster ());
      cluster_done->cloud = cloud_current;
      cluster_done->params = params;
      cluster_done->fit_error = fit_error;
      clusters_result.push_back (cluster_done);
    }
  }

  PCL_INFO ("Number of clusters after splitting: %ld\n   with sizes: ", clusters_result.size ());
  for (size_t i = 0; i < clusters_result.size (); ++i)
    PCL_INFO ("%ld ", clusters_result[i]->cloud->size ());
  PCL_INFO ("\n");
}




void
computeMergeCandidates (const std::vector<SuperquadricCluster::Ptr> &clusters,
                        const std::vector<Eigen::Vector4d> &splitting_planes,
                        std::set<std::pair<size_t, size_t> > &pairs)
{
  /// Compute the centroids for all of the clusters
  std::vector<Eigen::Vector4d> centroids (clusters.size ());
  for (size_t i = 0; i < clusters.size (); ++i)
    pcl::compute3DCentroid (*clusters[i]->cloud, centroids[i]);


  /// Compute by checking if the centroids are on opposite sides of each plane
  for (size_t c_i = 0; c_i < clusters.size (); ++c_i)
  {
    for (size_t c_j = c_i + 1; c_j < clusters.size (); ++c_j)
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

      if (planes_intersected <= 3)
      {
        if (clusters[c_i]->cloud->size () > clusters[c_j]->cloud->size ())
          pairs.insert (std::make_pair (c_i, c_j));
        else
          pairs.insert (std::make_pair (c_j, c_i));
      }
    }
  }


  /* // Compute by expliciting calculating the intersection


  /// Take all pairs of clusters
  for (size_t src_i = 0; src_i < clusters.size (); ++src_i)
    for (size_t tgt_i = src_i + 1; tgt_i < clusters.size (); ++tgt_i)
    {
      /// Count how many planes the line segment connecting the two clusters intersect
      size_t intersected_planes = 0;
      for (size_t plane_i = 0; plane_i < splitting_planes.size (); ++plane_i)
      {
        /// Compute the intersection of the segment connecting the two clusters and the splitting plane
        Eigen::Vector3d A (centroids[src_i].head<3> ()),
            B (centroids[tgt_i].head<3> ());
        Eigen::Vector3d n (splitting_planes[plane_i].head<3> ());
        double d = splitting_planes[plane_i][3];

        double t = - (A.dot (n) + d) / (B - A).dot (n);

        if (t > 0. && t < 1.)
          intersected_planes ++;
      }

      /// If we intersected less than 3 planes, then we consider the clusters as being neighbors
      if (intersected_planes <= 3)
        pairs.insert (std::make_pair (src_i, tgt_i));
    }
    */
}




void
mergeClusters (const std::vector<SuperquadricCluster::Ptr> &clusters_in,
               std::vector<SuperquadricCluster::Ptr> &clusters_merged,
               std::vector<Eigen::Vector4d> &splitting_planes)
{
  size_t merge_count = 0;

  clusters_merged = clusters_in;

  /// Repeat until we cannot do any more merges
  do
  {
    merge_count = 0;

    std::set<std::pair<size_t, size_t> > candidate_pairs;
    computeMergeCandidates (clusters_merged, splitting_planes, candidate_pairs);

    PCL_INFO ("Candidate pairs: %zu\n", candidate_pairs.size ());

    std::vector<std::set<std::pair<size_t, SuperquadricCluster::Ptr> > > possible_merges (clusters_merged.size ());

    for (std::set<std::pair<size_t, size_t> >::const_iterator p_it = candidate_pairs.begin (); p_it != candidate_pairs.end (); ++p_it)
    {
      size_t cluster_src_index = p_it->first,
          cluster_tgt_index = p_it->second;
      /// Compute the union of the clusters
      PointCloud<PointXYZ>::Ptr clusters_union (new PointCloud<PointXYZ> ());
      *clusters_union = *(clusters_merged[cluster_src_index]->cloud) + *(clusters_merged[cluster_tgt_index]->cloud);

      /// Center the cloud
      Eigen::Vector4d centroid;
      pcl::compute3DCentroid (*clusters_union, centroid);

      /// Fit a superquadric in the union
      sq::SuperquadricParameters<double> params;
      double fit_error = fitBestSuperquadric (clusters_union, params);
      PCL_INFO ("Superquadric fit error = %f\n", fit_error);

      if (fit_error < cluster_sq_error)
      {
        SuperquadricCluster::Ptr cluster (new SuperquadricCluster ());
        cluster->cloud = clusters_union;
        cluster->params = params;
        cluster->fit_error = fit_error;
        possible_merges[cluster_src_index].insert (std::make_pair (cluster_tgt_index, cluster));
      }
    }

    PCL_INFO ("Possible merges: ");
    for (size_t i = 0; i < possible_merges.size (); ++i)
    {
      PCL_INFO ("Merges for cluster %zu: ", i);
      for (std::set<std::pair<size_t, SuperquadricCluster::Ptr> >::const_iterator s_it = possible_merges[i].begin ();
           s_it != possible_merges[i].end (); ++s_it)
        PCL_INFO ("(%zu, %f)   ", s_it->first, s_it->second->fit_error);
      PCL_INFO ("\n");
    }


    std::vector<SuperquadricCluster::Ptr> clusters_new;
    std::vector<bool> taken (clusters_merged.size (), false);
    for (size_t i = 0; i < possible_merges.size (); ++i)
    {
      if (taken[i])
        continue;

      double min_error = std::numeric_limits<double>::max ();
      SuperquadricCluster::Ptr min_cluster;
      int min_index = -1;
      for (std::set<std::pair<size_t, SuperquadricCluster::Ptr> >::const_iterator s_it = possible_merges[i].begin ();
           s_it != possible_merges[i].end (); ++s_it)
        if (min_error > s_it->second->fit_error &&
            s_it->second->fit_error < cluster_sq_error &&
            !taken[s_it->first])
        {
          min_index = s_it->first;
          min_error = s_it->second->fit_error;
          min_cluster = s_it->second;
        }

      /// Combine the two clusters
      if (min_index != -1)
      {
        clusters_new.push_back (min_cluster);
        taken[min_index] = true;
        taken[i] = true;
        merge_count ++;
      }
    }

    for (size_t i = 0; i < clusters_merged.size (); ++i)
      if (!taken[i])
        clusters_new.push_back (clusters_merged[i]);


    PCL_INFO ("Clusters merged: %zu - old size %zu, new size %zu\n",
              merge_count, clusters_merged.size (), clusters_new.size ());

    clusters_merged = clusters_new;
  } while (merge_count != 0);
}




int
main (int argc,
      char **argv)
{
  google::InitGoogleLogging (argv[0]);

  std::string cloud_path = "";
  console::parse_argument (argc, argv, "-cloud", cloud_path);

  std::string output_dir = "";
  console::parse_argument (argc, argv, "-output", output_dir);

  if (cloud_path == "" ||
      output_dir == "")
  {
    PCL_ERROR ("Syntax error. Correct: %s -cloud <path_to_input_cloud> -output <output folder to write the results in>\n", argv[0]);
    return (-1);
  }

  PointCloud<PointXYZ>::Ptr cloud (new PointCloud<PointXYZ> ());
  io::loadPCDFile (cloud_path, *cloud);


  PCL_INFO ("####### SPLITTING ######\n");

  std::vector<SuperquadricCluster::Ptr> clusters_split;
  std::vector<Eigen::Vector4d> splitting_planes;
  splitClusters (cloud, clusters_split, splitting_planes);

  PCL_INFO ("Clusters after splitting:\n");
  for (size_t i = 0; i < clusters_split.size (); ++i)
  {
    PCL_INFO ("Command for sampler:\n-e1 %f -e2 %f -a %f -b %f -c %f -transform %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
              clusters_split[i]->params.e1, clusters_split[i]->params.e2, clusters_split[i]->params.a, clusters_split[i]->params.b, clusters_split[i]->params.c,
              clusters_split[i]->params.transform (0, 0), clusters_split[i]->params.transform (0, 1), clusters_split[i]->params.transform (0, 2), clusters_split[i]->params.transform (0, 3),
              clusters_split[i]->params.transform (1, 0), clusters_split[i]->params.transform (1, 1), clusters_split[i]->params.transform (1, 2), clusters_split[i]->params.transform (1, 3),
              clusters_split[i]->params.transform (2, 0), clusters_split[i]->params.transform (2, 1), clusters_split[i]->params.transform (2, 2), clusters_split[i]->params.transform (2, 3),
              clusters_split[i]->params.transform (3, 0), clusters_split[i]->params.transform (3, 1), clusters_split[i]->params.transform (3, 2), clusters_split[i]->params.transform (3, 3));
    PCL_INFO ("\n");
  }

  PCL_INFO ("####### MERGING ######\n");


  /// HACK: filter out the clusters that are below 500 points
  std::vector<SuperquadricCluster::Ptr> clusters_split_filtered;
  for (size_t i = 0; i < clusters_split.size (); ++i)
    if (clusters_split[i]->cloud->size () > 500)
      clusters_split_filtered.push_back (clusters_split[i]);

  std::vector<SuperquadricCluster::Ptr> clusters_merge;
  mergeClusters (clusters_split_filtered, clusters_merge, splitting_planes);


  PCL_INFO ("Clusters after merging:\n");
  for (size_t i = 0; i < clusters_merge.size (); ++i)
  {
    PCL_INFO ("Command for sampler:\n-e1 %f -e2 %f -a %f -b %f -c %f -transform %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
              clusters_merge[i]->params.e1, clusters_merge[i]->params.e2, clusters_merge[i]->params.a, clusters_merge[i]->params.b, clusters_merge[i]->params.c,
              clusters_merge[i]->params.transform (0, 0), clusters_merge[i]->params.transform (0, 1), clusters_merge[i]->params.transform (0, 2), clusters_merge[i]->params.transform (0, 3),
              clusters_merge[i]->params.transform (1, 0), clusters_merge[i]->params.transform (1, 1), clusters_merge[i]->params.transform (1, 2), clusters_merge[i]->params.transform (1, 3),
              clusters_merge[i]->params.transform (2, 0), clusters_merge[i]->params.transform (2, 1), clusters_merge[i]->params.transform (2, 2), clusters_merge[i]->params.transform (2, 3),
              clusters_merge[i]->params.transform (3, 0), clusters_merge[i]->params.transform (3, 1), clusters_merge[i]->params.transform (3, 2), clusters_merge[i]->params.transform (3, 3));
    PCL_INFO ("\n");
  }


  /// Generate the cluster point clouds for easy result visualization
  boost::filesystem::create_directory (output_dir.c_str ());
  for (size_t i = 0; i < clusters_merge.size (); ++i)
  {
    sq::SuperquadricSampling<PointXYZ, double> sampling;
    sampling.setParameters (clusters_merge[i]->params);

    PointCloud<PointXYZ> cluster_cloud;
    sampling.generatePointCloud (cluster_cloud);

    char str[512];
    sprintf (str, "%s/cluster_%zu.pcd", output_dir.c_str (), i);
    io::savePCDFileBinaryCompressed (str, cluster_cloud);

    PolygonMesh cluster_mesh;
    sampling.generateMesh (cluster_mesh);
    sprintf (str, "%s/cluster_%zu.vtk", output_dir.c_str (), i);
    io::savePolygonFileVTK (str, cluster_mesh);
  }



  return (0);
}
