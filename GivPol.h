//!typedef Givaro::Poly1Dom<IntField, Dense>::Rep polynomial;
#ifndef GivPol_H
#define GivPol_H

#include <vector>

namespace Givaro
{
	template <typename T>
	class givvector : public std::vector<T>
	{};

	template <class Domain>
	class Poly1Dom //Dense
	{
	private:
		typedef typename Domain::Element Element;
		Domain *domain;
	public:
		Poly1Dom(Domain& d) {domain = &d;}
		typedef givvector<Element> Rep;
		Rep& mod(Rep &q, Rep &a, Rep &b)
		{
			Rep theta;
			theta = a;
			size_t theta_size = theta.size(); size_t b_size = b.size();
			while (theta_size >= b_size)
			{
				Rep delta;
				size_t k;
				for (k = 0; k < theta_size - b_size; k++)
				{
					delta.push_back(0);
				}
				for (; k < theta_size; k++)
				{
					delta.push_back(domain->mul(theta[theta_size - 1], b[k - theta_size + b_size]));
				}
				subin(theta, delta);
				theta_size--;
			}
			for (size_t k = 0; k < b_size - 1; k++)
			{
				q.push_back(theta[k]);
			}
			return q;
		}

		Rep& modin(Rep &a, Rep &b) {Rep q; a = mod(q, a, b); return a;}

		Rep& subin(Rep &q, Rep &b)
		{
			for (size_t i = 0; i < b.size(); i++) q[i] = domain->minus(q[i], b[i]);
			return q;
		}
	};
}

#endif
