#include <eigen3/Eigen/Dense>
#include <iostream>

using Eigen::MatrixXd;
using namespace std;
using namespace Eigen;

struct logical_xor {
  unsigned int operator() (unsigned int a, unsigned int b) const
  {
    return a != b;
  }
};


template<class ArgType, class RowIndexType, class ColIndexType>
class indexing_functor {
  const ArgType &m_arg;
  const RowIndexType &m_rowIndices;
  const ColIndexType &m_colIndices;
public:
  typedef Matrix<typename ArgType::Scalar,
                 RowIndexType::SizeAtCompileTime,
                 ColIndexType::SizeAtCompileTime,
                 ArgType::Flags&RowMajorBit?RowMajor:ColMajor,
                 RowIndexType::MaxSizeAtCompileTime,
                 ColIndexType::MaxSizeAtCompileTime> MatrixType;
  indexing_functor(const ArgType& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
    : m_arg(arg), m_rowIndices(row_indices), m_colIndices(col_indices)
  {}
  const typename ArgType::Scalar& operator() (Index row, Index col) const {
    return m_arg(m_rowIndices[row], m_colIndices[col]);
  }
};


template <class ArgType, class RowIndexType, class ColIndexType>
CwiseNullaryOp<indexing_functor<ArgType,RowIndexType,ColIndexType>, typename indexing_functor<ArgType,RowIndexType,ColIndexType>::MatrixType>
indexing(const Eigen::MatrixBase<ArgType>& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
{
  typedef indexing_functor<ArgType,RowIndexType,ColIndexType> Func;
  typedef typename Func::MatrixType MatrixType;
  return MatrixType::NullaryExpr(row_indices.size(), col_indices.size(), Func(arg.derived(), row_indices, col_indices));
}

int main(){
/*int array[8];
int b[8];
for(int i = 0; i < 8; ++i) array[i] = 1;
for(int i=0; i < 8; ++i) b[i] = 1;
cout << "Column-major:\n" << Map<Matrix<int,2,4> >(array) << endl;
cout << "Row-major:\n" << Map<Matrix<int,2,4,RowMajor> >(array) << endl;
cout << "Row-major using stride:\n" <<
  *///Map<Matrix<int,2,4>, Unaligned, Stride<1,4> >(array) << endl;
  Eigen::MatrixXi A = Eigen::MatrixXi::Random(4,4);
  Array3i ri(1,1,1);
  ArrayXi ci(6); ci << 1,1,1, 1,1, 1;
  Eigen::MatrixXi B = indexing(A, ri, ci);
  std::cout << "A =" << std::endl;
  std::cout << A << std::endl << std::endl;
  std::cout << "A([" << ri.transpose() << "], [" << ci.transpose() << "]) =" << std::endl;
  std::cout << B << std::endl;


  B =  indexing(A, ri+1, ci);
  std::cout << "A(ri+1,ci) =" << std::endl;
  std::cout << B << std::endl << std::endl;
#if __cplusplus >= 201103L
  B =  indexing(A, ArrayXi::LinSpaced(13,0,12).unaryExpr([](int x){return x%4;}), ArrayXi::LinSpaced(4,0,3));
  std::cout << "A(ArrayXi::LinSpaced(13,0,12).unaryExpr([](int x){return x%4;}), ArrayXi::LinSpaced(4,0,3)) =" << std::endl;
  std::cout << B << std::endl << std::endl;
#endif

}
