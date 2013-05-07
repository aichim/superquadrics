#include <ceres/ceres.h>
#include <glog/logging.h>


struct CostFunctor
{
  template <typename T> bool
  operator () (const T* const x, T* residual) const
  {
    residual[0] = T (10.0) - x[0];
    return (true);
  }
};



int
main (int argc,
      char **argv)
{
  google::InitGoogleLogging (argv[0]);

  double initial_x = 5.0;
  double x = initial_x;

  ceres::Problem problem;

  ceres::CostFunction *cost_function = new ceres::AutoDiffCostFunction<CostFunctor, 1, 1> (new CostFunctor);
  problem.AddResidualBlock (cost_function, NULL, &x);

  ceres::Solver::Options options;
  options.minimizer_type = ceres::TRUST_REGION;
  options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
  options.minimizer_progress_to_stdout = true;

  ceres::Solver::Summary summary;
  ceres::Solve (options, &problem, &summary);

  std::cout << summary.BriefReport () << "\n";
  std::cout << "x : " << initial_x << " -> " << x << "\n";

  return (0);
}
