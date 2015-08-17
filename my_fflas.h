#ifndef __my_FFLAS_H
#define __my_FFLAS_H
namespace FFLAS
{
	enum FFLAS_TRANSPOSE {FflasNoTrans, FfslasTrans};

	template<class Field>
	typename Field::Element_ptr
	fgemm(Field& F, FFLAS_TRANSPOSE ta, FFLAS_TRANSPOSE tb, size_t m,
		size_t n, size_t q, typename Field::Element alpha,
		typename Field::Element_ptr A, size_t lda,
		typename Field::Element_ptr B, size_t ldb,
		typename Field::Element beta,
		typename Field::Element_ptr C, size_t ldc)
	{
		typename Field::Element el1;
		typename Field::Element el2;
		typename Field::Element el3;
		typename Field::Element p = F.characteristic();
		for (size_t i = 0; i < m; i++)
		{
			for (size_t j = 0; j < n; j++)
			{
				el1 = 0;
				for (size_t k = 0; k < q; k++)
				{
					el2 = A[i*lda + k]; el3 = B[k*ldb + j];
					el1 += (el2 * el3);
				}					
				el1 = el1 % p;
				C[i*ldc + j] = el1;
			}
		}
	}
}

#endif
