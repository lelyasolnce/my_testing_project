#ifndef __my_FFPACK_H
#define __my_FFPACK_H

namespace FFPACK
{
	template<class Field>
	typename Field::Element_ptr
	Invert(Field& F, size_t M,
		typename Field::Element_ptr A, size_t lda,
		typename Field::Element_ptr C, size_t ldc,
		int& nullity)
	{

		typedef typename Field::Element Element;
		size_t i = 0;
		size_t j = 0;
		size_t s = 0;
		size_t t = 0;
		typename Field::Element_ptr D = new Element [M * M];
		for (size_t i = 0; i < M; i++)
			for (size_t j = 0; j < M; j++)
				D[i*M + j] = A[i*lda + j];
		for (size_t i = 0; i < M; i++)
			for (size_t j = 0; j < M; j++)
				C[i*ldc +j] = (i == j) ? F.one : F.zero;
		while ((i < M) && (j < M))
		{
			s = i;
			t = j;
			bool nonzero = false;
			while ((!nonzero) && (t < M))
			{
				if (D[s*M+t] != F.zero)
				{
					nonzero = true;
				}
				else
				{
					if (s == M - 1)
					{
						s = i;
						t++;
					}
					else
					{
						s++;
					}
				}
			}
			if (!nonzero) return C;
			else
			{
				typename Field::Element el;
				for (size_t u = 0; u < M; u++) 
				{
					el = D[i*M + u];
					D[i*M + u] = D[s*M + u];
					D[s*M + u] = el;
					el = C[i*ldc + u];
					C[i*ldc + u] = C[s*ldc + u];
					C[s*ldc + u] = el;
				}
				el = D[i*M + t];
				el = F.invert(el);
				for (size_t u = 0; u < M; u++) 
				{
					D[i*M + u] = F.mul(D[i*M + u], el);
					C[i*ldc + u] = F.mul(C[i*ldc + u], el);
				}
				for (size_t v = 0; v < M; v++) 
				{
					if (v != i)
					{
Element ell = D[v*M + t];
						for (size_t u = 0; u < M; u++) 
						{
							D[v*M +u] = F.minus(D[v*M + u], F.mul(ell, D[i*M + u]));
							C[v*ldc +u] = F.minus(C[v*ldc + u], F.mul(ell, C[i*ldc + u]));
						}
					}
				}
				i++;
				j++;
			}
		}
		return C;
	}
}

#endif
