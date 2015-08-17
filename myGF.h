#ifndef __GIVARO_myGF_H
#define __GIVARO_myGF_H

#include <vector>

namespace Givaro
{
	template<class XXX>
	class GFqDom
	{
	public:
		typedef XXX Element;
	private:
		int p;
		int n;
		std::vector<XXX> polynomial;
	public:
		GFqDom (XXX p0, XXX n0)
		{
			p = p0;
			n = n0;
			//polynomial;
		}
		int& exponent() {return n;}
		int& characteristic() {return p;}
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
