#ifndef __GIVARO_modular_H
#define __GIVARO_modular_H

namespace Givaro
{
	template <class XXX>
	class Modular<XXX>
	{
	public:
		typedef XXX Element;
		typedef Modular<XXX> Self_t;
		typedef uint64_t Residu_t;
		// ----- Constantes
		const Element zero;
		const Element one;
		const Element mOne;
		const Element p;
		// ----- Constructors
		Modular();
		Modular(Residu_t p);
		Modular(const Self_t& F);
	};

	template <>
	class Modular<double>
	{
	public:
		typedef XXX Element;
		typedef Modular<XXX> Self_t;
		typedef uint64_t Residu_t;
		// ----- Constantes
		const Element zero;
		const Element one;
		const Element mOne;
		const Element _p;
		// ----- Constructors
		Modular() : zero(0.0), one(1.0), mOne(-1.0), _p(0.0) {}
		Modular(const XXX& p) : zero(0.0), one(1.0), mOne((Element)p - 1.0), _p((double)p) {}
		Modular(const Self_t& F) : zero(F.zero), one(F.one), mOne(F.mOne), _p(F._p) {}
	};

	template<>
	class Modular<int64_t> 
	{
	public:
		typedef XXX Element;
		typedef Modular<XXX> Self_t;
		typedef uint64_t Residu_t;
		// ----- Constantes
		const Element zero;
		const Element one;
		const Element mOne;
		const Element p;
		// ----- Constructors
		Modular() : zero(0), one(1), mOne(-1), _p(0) {}
		Modular(Residu_t p) : zero(0), one(1), mOne((Element)p-1), _p(p) {}
		Modular(const Self_t& F) : zero(F.zero), one(F.one), mOne(F.mOne), _p(F._p) {}
	};
}

#endif
