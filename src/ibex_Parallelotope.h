/* ============================================================================
 * I B E X 
 * ============================================================================
 * Copyright   : Ecole Polytechnique (FRANCE)
 * License     : LGPL, see the file COPYING.LESSER.
 * Author(s)   : Benjamin Martin
 * Created     : 04.02.2019
 * Last Update : 
 * ---------------------------------------------------------------------------- */

#ifndef __IBEX_PARALLELOTOPE_H__
#define __IBEX_PARALLELOTOPE_H__

#include "ibex_IntervalMatrix.h"


namespace ibex
{



/**
 * \brief Class representing a parallelotope: a tilted box.
 * 
 * The attributes of the class are:
 * - a transformation matrix C
 * - a box w
 * - a vector xtilde representing the center
 * 
 * The concretization of the parallelotope is then given by C w + xtilde
 * 
 * There is an extra attribute Cinv: an IntervalMatrix containing the inverse
 * of C. This is used for set interesection of inclusion tests.
 * 
 * TODO: Think about how we can optimize the storage of Parallelotope.
 * There are two dense matrices which can be cumbersome for large dimension.
 * Cinv is only used for testing: once we know that no more tests will happen
 * on the parallelotope, it can be supressed. C and Cinv could be shared between
 * parallelotopes having the same transformation. 
 * 
 **/
class Parallelotope { 
	
	public:

		/**
		 * 	\brief Computes a transformation matrix given an approximate jacobian
		 * 	and a tangent vector. 
		 **/
		static Matrix transformation(const Matrix& J, const Vector& tan);
		
		/**
		 *  \brief Computes the box w of size n with only 0 excepts w[n-1] = Interval(0,h)
		 **/
		static IntervalVector auxiliary_box(int n, double h);

		/**
		 * 	\brief Computes an IntervalMatrix enclosing the inverse of C
		 **/
		static IntervalMatrix computeInvMatrix(const Matrix& C);

		/**
		 * 	\brief Parallelotope constructor.
		 * 
		 * 	The parallelotope is constructed given:
		 * 	- a center xtilde,
		 * 	- a jacobian matrix jac (at xtilde),
		 * 	- a step length h,
		 * 	- a vector tan approximately tangent to jac.
		 **/		
		Parallelotope(const Vector& xtilde, const Matrix& jac, double h, const Vector& tan)
		:	C(transformation(jac, tan)),
			w(auxiliary_box(xtilde.size(),h)),
			xtilde(xtilde),
			Cinv(computeInvMatrix(C))
		{}
		
		/**
		 * 	\brief Partial copy constructor.
		 * 
		 * 	The transition matrix, its inverse and the center of p are copied.
		 * 	The characteristic box is set to wn.
		 **/
		Parallelotope(const Parallelotope& p, const IntervalVector& wn)
		:	C(p.C),
			w(wn),
			xtilde(p.xtilde),
			Cinv(p.Cinv)
		{}
		
		/**
		 * 	\brief Returns the interval hull of the parallelotope
		 **/
		IntervalVector hull() const;
		
		/**
		 * 	\brief Empty intersection test w.r.t another parallelotope p.
		 **/
		bool not_intersects(const Parallelotope& p) const;
		
		/**
		 * 	\brief Empty intersection test w.r.t a box x.
		 **/
		bool not_intersects(const IntervalVector& x) const;
		
		/**
		 * 	\brief Tests whether the parallelotope contains the box x.
		 **/
		bool is_superset(const IntervalVector& x) const;
		
		/**
		 * 	\brief The size of the parallelotope: the dimension of the space
		 * 	in which it is embedded
		 **/
		size_t size() const;
		
		
		// The characteristic matrix
		Matrix C;
		
		// The characteristic box
		IntervalVector w;
		
		// The center
		Vector xtilde;
		
		// The enclosure of the inverse of C. Used by the tests.
		IntervalMatrix Cinv;
		
}; // Parallelotope


inline IntervalVector Parallelotope::hull() const
{
	return C*(w) + xtilde;
}

inline bool Parallelotope::not_intersects(const Parallelotope& p) const
{
	// TODO: the test can be more accurate by developing the matrix product
	// involved in p.hull
	return ! (Cinv*(p.hull() - xtilde)).intersects(w);
}

inline bool Parallelotope::not_intersects(const IntervalVector& x) const
{		
	return ! w.intersects(Cinv*(x-xtilde));
}

inline bool Parallelotope::is_superset(const IntervalVector& x) const
{
	return (Cinv*(x - xtilde)).is_subset(w);	
}

inline size_t Parallelotope::size() const
{
	return xtilde.size();
}

} // ibex

#endif /* __IBEX_PARALLELOTOPE_H__ */
