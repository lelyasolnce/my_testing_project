#ifndef __LINBOX_bm_addsub_INL
#define __LINBOX_bm_addsub_INL

namespace LinBox
{
	template<class Field, class Matrix1, class Matrix2, class Matrix3>
	class BlasMatrixDomainAdd
	{
	public:
		Matrix1& operator() (Field& F, Matrix1& C, Matrix2& A, Matrix3& B)
		{
			typename Matrix2::Element el2; typename Matrix3::Element el3;
			for (size_t i = 0; i < C.rowdim(); i++)
			{
				for (size_t j = 0; j < C.coldim(); j++)
				{
					el2 = A.getEntry(i,j); el3 = B.getEntry(i,j); el2 = el2+el3;
					C.setEntry(i, j, el2);
				}
			}
			return C;
		}
	};

	template<class Field, class Matrix1, class Matrix2, class Matrix3>
	class BlasMatrixDomainSub
	{
	public:
		Matrix1& operator() (Field& F, Matrix1& C, Matrix2& A, Matrix3& B)
		{
			typename Matrix2::Element el2; typename Matrix3::Element el3;
			for (size_t i = 0; i < C.rowdim(); i++)
			{
				for (size_t j = 0; j < C.coldim(); j++)
				{
					el2 = A.getEntry(i,j); el3 = B.getEntry(i,j); el2 = el2-el3;
					C.setEntry(i, j, el2);
				}
			}
			return C;
		}
	};

	template<class Field, class Matrix1, class Matrix2>
	class BlasMatrixDomainAddin
	{
	public:
		Matrix1& operator() (Field& F, Matrix1& C, Matrix2& A)
		{
			typename Matrix2::Element el2; typename Matrix1::Element el3;
			for (size_t i = 0; i < C.rowdim(); i++)
			{
				for (size_t j = 0; j < C.coldim(); j++)
				{
					el2 = A.getEntry(i,j); el3 = C.getEntry(i,j); el2 = el2+el3;
					C.setEntry(i, j, el2);
				}
			}
			return C;
		}
	};

	template<class Field, class Matrix1, class Matrix2>
	class BlasMatrixDomainSubin
	{
	public:
		Matrix1& operator() (Field& F, Matrix1& C, Matrix2& A)
		{
			typename Matrix2::Element el2; typename Matrix1::Element el3;
			for (size_t i = 0; i < C.rowdim(); i++)
			{
				for (size_t j = 0; j < C.coldim(); j++)
				{
					el2 = A.getEntry(i,j); el3 = C.getEntry(i,j); el2 = el3 - el2;
					C.setEntry(i, j, el2);
				}
			}
			return C;
		}
	};
} // LinBox

#endif
