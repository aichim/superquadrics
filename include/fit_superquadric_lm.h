#pragma once

#include <pcl/point_cloud.h>
#include <pcl/common/common.h>


namespace sq
{

template <typename t>
struct SuperquadricParameters;

template <typename PointT, typename MatScalar = double>
class SuperquadricFittingLM
{
  typedef pcl::PointCloud<PointT> Cloud;
  typedef typename Cloud::Ptr CloudPtr;
  typedef typename Cloud::ConstPtr CloudConstPtr;
  typedef Eigen::Matrix<MatScalar, Eigen::Dynamic, 1> VectorX;

public:
  SuperquadricFittingLM ();

  SuperquadricFittingLM (const SuperquadricFittingLM &src)
  { }

  SuperquadricFittingLM& operator = (const SuperquadricFittingLM &src)
  { }

  virtual
  ~SuperquadricFittingLM () {}

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

  void
  setIndices (const pcl::IndicesConstPtr &indices)
  { indices_ = indices; }

  void
  preAlign (Eigen::Matrix<MatScalar, 4, 4> &transformation_prealign,
            Eigen::Matrix<MatScalar, 3, 1> &variances);

  double
  fit (SuperquadricParameters<MatScalar> &parameters);


protected:
  CloudConstPtr input_;
  CloudPtr input_prealigned_;
  pcl::IndicesConstPtr indices_;
  bool pre_align_;
  int pre_align_axis_;

  /** \brief The vector of residual weights. Used internally in the LM loop. */
  std::vector<double> weights_;


  /** Base functor all the models that need non linear optimization must
    * define their own one and implement operator() (const Eigen::VectorXd& x, Eigen::VectorXd& fvec)
    * or operator() (const Eigen::VectorXf& x, Eigen::VectorXf& fvec) dependening on the choosen _Scalar
    */
  template<typename _Scalar, int NX=Eigen::Dynamic, int NY=Eigen::Dynamic>
  struct Functor
  {
    typedef _Scalar Scalar;
    enum
    {
      InputsAtCompileTime = NX,
      ValuesAtCompileTime = NY
    };
    typedef Eigen::Matrix<_Scalar,InputsAtCompileTime,1> InputType;
    typedef Eigen::Matrix<_Scalar,ValuesAtCompileTime,1> ValueType;
    typedef Eigen::Matrix<_Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

    /** \brief Empty Construtor. */
    Functor () : m_data_points_ (ValuesAtCompileTime) {}

    /** \brief Constructor
            * \param[in] m_data_points number of data points to evaluate.
            */
    Functor (int m_data_points) : m_data_points_ (m_data_points) {}

    /** \brief Destructor. */
    virtual ~Functor () {}

    /** \brief Get the number of values. */
    int
    values () const { return (m_data_points_); }

  protected:
    int m_data_points_;
  };

  struct OptimizationFunctor : public Functor<MatScalar>
  {
    using Functor<MatScalar>::values;

    /** Functor constructor
            * \param[in] m_data_points the number of data points to evaluate
            * \param[in,out] estimator pointer to the estimator object
            */
    OptimizationFunctor (int m_data_points,
                         const SuperquadricFittingLM *estimator)
      :  Functor<MatScalar> (m_data_points), estimator_ (estimator)
    {}

    /** Copy constructor
            * \param[in] the optimization functor to copy into this
            */
    inline OptimizationFunctor (const OptimizationFunctor &src) :
      Functor<MatScalar> (src.m_data_points_), estimator_ ()
    {
      *this = src;
    }

    /** Copy operator
      * \param[in] the optimization functor to copy into this
      */
    inline OptimizationFunctor&
    operator = (const OptimizationFunctor &src)
    {
      Functor<MatScalar>::operator=(src);
      estimator_ = src.estimator_;
      return (*this);
    }

    /** \brief Destructor. */
    virtual ~OptimizationFunctor () {}

    /** Fill fvec from x. For the current state vector x fill the f values
            * \param[in] x state vector
            * \param[out] fvec f values vector
            */
    int
    operator () (const VectorX &x,
                 VectorX &fvec) const;

    int
    df (const VectorX &x,
        Eigen::Matrix<MatScalar, Eigen::Dynamic, Eigen::Dynamic> &fjac) const;

    const SuperquadricFittingLM<PointT, MatScalar> *estimator_;
  };

  struct OptimizationFunctorWithIndices : public Functor<MatScalar>
  {
    using Functor<MatScalar>::values;

    /** Functor constructor
            * \param[in] m_data_points the number of data points to evaluate
            * \param[in,out] estimator pointer to the estimator object
            */
    OptimizationFunctorWithIndices (int m_data_points,
                                    const SuperquadricFittingLM *estimator)
      : Functor<MatScalar> (m_data_points), estimator_ (estimator)
    {}

    /** Copy constructor
            * \param[in] the optimization functor to copy into this
            */
    inline OptimizationFunctorWithIndices (const OptimizationFunctorWithIndices &src)
      : Functor<MatScalar> (src.m_data_points_), estimator_ ()
    {
      *this = src;
    }

    /** Copy operator
      * \param[in] the optimization functor to copy into this
      */
    inline OptimizationFunctorWithIndices&
    operator = (const OptimizationFunctorWithIndices &src)
    {
      Functor<MatScalar>::operator=(src);
      estimator_ = src.estimator_;
      return (*this);
    }

    /** \brief Destructor. */
    virtual ~OptimizationFunctorWithIndices () {}

    /** Fill fvec from x. For the current state vector x fill the f values
      * \param[in] x state vector
      * \param[out] fvec f values vector
      */
    int
    operator () (const VectorX &x,
                 VectorX &fvec) const;

    const SuperquadricFittingLM<PointT, MatScalar> *estimator_;
  };

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}

#include "impl/fit_superquadric_lm.hpp"


