#include <iostream>
#include <eigen3/Eigen/Dense>
using namespace Eigen;
template <class ArgType> class Circulant;
namespace Eigen {
  namespace internal {
    template <class ArgType>
    struct traits<Circulant<ArgType> >
    {
      typedef Eigen::Dense StorageKind;
      typedef Eigen::MatrixXpr XprKind;
      typedef typename ArgType::StorageIndex StorageIndex;
      typedef typename ArgType::Scalar Scalar;
      enum { 
        Flags = Eigen::ColMajor,
        RowsAtCompileTime = ArgType::RowsAtCompileTime,
        ColsAtCompileTime = ArgType::RowsAtCompileTime,
        MaxRowsAtCompileTime = ArgType::MaxRowsAtCompileTime,
        MaxColsAtCompileTime = ArgType::MaxRowsAtCompileTime
      };
    };
  }
}


template <class ArgType>
class Circulant : public Eigen::MatrixBase<Circulant<ArgType> >
{
public:
  Circulant(const ArgType& arg, const ArgType& arg2, const ArgType& arg3)
    : m_arg(arg), s_arg(arg2), o_arg(arg3)
  { 
    EIGEN_STATIC_ASSERT(ArgType::ColsAtCompileTime == 1,
                        YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_A_MATRIX);
  }
  typedef typename Eigen::internal::ref_selector<Circulant>::type Nested; 
  typedef Eigen::Index Index;
  Index Rows() const { return m_arg.rows(); }
  Index Cols() const { return m_arg.rows(); }
  Index bRows() const { return s_arg.rows(); }
  Index bCols() const { return s_arg.rows(); }
  Index cRows() const { return o_arg.rows(); }
  Index cCols() const { return o_arg.rows(); }
  typedef typename Eigen::internal::ref_selector<ArgType>::type ArgTypeNested;
  ArgTypeNested m_arg;
  ArgTypeNested s_arg;
  ArgTypeNested o_arg;
};

namespace Eigen {
  namespace internal {
    template<typename ArgType>
    struct evaluator<Circulant<ArgType> >
      : evaluator_base<Circulant<ArgType> >
    {
      typedef Circulant<ArgType> XprType;
      typedef typename nested_eval<ArgType, XprType::ColsAtCompileTime>::type ArgTypeNested;
      typedef typename remove_all<ArgTypeNested>::type ArgTypeNestedCleaned;
      typedef typename XprType::CoeffReturnType CoeffReturnType;
      enum { 
        CoeffReadCost = evaluator<ArgTypeNestedCleaned>::CoeffReadCost,
        Flags = Eigen::ColMajor 
      };
      
      evaluator(const XprType& xpr)
        : m_argImpl(xpr.m_arg), m_rows(xpr.rows())
      { }
      CoeffReturnType coeff(Index row, Index col) const
      {
        Index index = row - col;
        
	if (index < 0) index += m_rows;
	
        return m_argImpl.coeff(index);
      }
      evaluator<ArgTypeNestedCleaned> m_argImpl;
      evaluator<ArgTypeNestedCleaned> m_argImp2;
      const Index m_rows;
    };
  }
}


template <class ArgType>
Circulant<ArgType> makeCirculant(const Eigen::MatrixBase<ArgType>& arg, const Eigen::MatrixBase<ArgType>& arg2, const Eigen::MatrixBase<ArgType>& arg3)
{
  return Circulant<ArgType>(arg.derived(), arg2.derived(), arg3.derived());
}


int main()
{
  Eigen::VectorXd vec(4);
  vec << 1, 2, 4, 8;
  Eigen::VectorXd vec2(4);
  Eigen::MatrixXd mat;
  mat = makeCirculant(vec, vec, vec2);
  std::cout << mat << std::endl;
}
