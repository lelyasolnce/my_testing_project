#ifndef __GIVARO_myGF_H
#define __GIVARO_myGF_H

#include <vector>

template<class XXX> struct Signed_Trait
{
	typedef XXX signed_type;
	typedef XXX unsigned_type;
};
template<> struct Signed_Trait<int64_t>
{
	typedef int64_t signed_type;
	typedef uint64_t unsigned_type;
};
template<> struct Signed_Trait<double>
{
	typedef double signed_type;
	typedef double unsigned_type;
};

namespace Givaro
{
	template<class XXX>
	class GFqDom
	{
	public:
		typedef XXX Element;
		typedef typename Signed_Trait<XXX>::unsigned_type Residu_t;
	private:
		int p;
		int n;
		std::vector<XXX> polynomial;
		Residu_t _irred;
	public:
		GFqDom (XXX p0, XXX n0)
		{
			p = p0;
			n = n0;
			_irred = 1; for (int k = 0; k < n; k++) _irred *= p;
		}
		int& exponent() {return n;}
		int& characteristic() {return p;}
		Residu_t irreducible() { return _irred;}
		GFqDom<XXX>& operator=(const GFqDom<XXX>& F)
		{
			this->p = F.p;
			this->n = F.n;
			this->polynomial = F.polynomial;
			return *this;
		}
	};
}

#endif
