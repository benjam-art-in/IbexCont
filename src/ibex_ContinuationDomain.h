/* ============================================================================
 * I B E X 
 * ============================================================================
 * Copyright   : Ecole Polytechnique (FRANCE)
 * License     : LGPL, see the file COPYING.LESSER.
 * Author(s)   : Benjamin Martin
 * Created     : 04.02.2019
 * Last Update : 
 * ---------------------------------------------------------------------------- */

#ifndef __IBEX_CONTINUATION_DOMAIN_H__
#define __IBEX_CONTINUATION_DOMAIN_H__

#include "ibex_Parallelotope.h"
#include "ibex_Function.h"

namespace ibex {


/**
 *	\brief Wrapper class for a domain paving a manifold.
 *
 *	Domain specifically used for continuation of 1-dimensional manifolds.
 *  Once a domain is certified, the manifold traverse the domain with
 *  respect to a given direction.
 **/
class ContinuationDomain
{
	public:
		/**
		 * \brief Certify the existence and unicity of a piece
		 * of manifold of solution to f = 0.
		 **/
		virtual bool certify(Function& f) = 0;
		
		/**
		 * 	\brief Contract on the side of the domain where the manifold "enters"
		 *  the domain. Returns the corresponding domain.
		 **/
		virtual ContinuationDomain* contractIn(Function& f) = 0;
		
		/**
		 * 	\brief Contract on the side of the domain where the manifold "exits"
		 *  the domain. Returns the corresponding domain
		 **/
		virtual ContinuationDomain* contractOut(Function& f) = 0;
		
		/**
		 * 	\brief Empty intersection test between the domain and a box.
		 **/
		virtual bool not_intersects(const IntervalVector& x) const = 0;
		
		/**
		 *	\brief Empty intersection test between the domain and another domain.
		 *
		 *  The hull of the other domain is taken.
		 **/
		bool not_intersects(const ContinuationDomain& d) const;
		
		/**
		 * 	\brief Tests whether the domain contains a box.
		 **/
		virtual bool is_superset(const IntervalVector& b) const = 0;
		
		/**
		 * 	\brief Tests whether the domain contains the domain d.
		 * 
		 *	The hull of d is taken.
		 **/
		bool is_superset(const ContinuationDomain& d) const;

		/**
		 * 	\brief Builds an interval hull of the domain.
		 **/
		virtual IntervalVector hull() const = 0;
		
		virtual void print() const = 0;
		
}; // ContinuationDomain

inline bool ContinuationDomain::not_intersects(const ContinuationDomain& d) const
{
	return not_intersects(d.hull());
}

inline bool ContinuationDomain::is_superset(const ContinuationDomain& d) const
{
	return is_superset(d.hull());
}
		
/**
 * 	\brief Class implementing the Box domain for continuation. 
 * 
 * Used for implementing the BoxCont continuation algorithm.
 * 
 * TODO: not necessary to store a varset explicitly ? 
 **/
class ContinuationDomainBox : public ContinuationDomain
{
	public:
	
		/**
		 * \brief Basic constructor
		 * 
		 * Construct a box continuation domain with box b, a varset v containing 
		 * a single paramter being the parametric dimension of the box domain
		 * and a sign s indicating if the parameter dimension is crossed in the 
		 * positive direction.
		 **/
		ContinuationDomainBox(const IntervalVector& b, const VarSet& v, bool s)
		: 	x(b),
			vs(v),
			sign(s)
		{}
		
		/**
		 * \brief Certify the existence and unicity of a piece
		 * of manifold of solution to f = 0.
		 **/
		virtual bool certify(Function& f);
		
		/**
		 * 	\brief Contract on the side of the domain where the manifold "enters"
		 *  the domain. Returns the corresponding domain.
		 **/
		virtual ContinuationDomain* contractIn(Function& f);
		
		/**
		 * 	\brief Contract on the side of the domain where the manifold "exits"
		 *  the domain. Returns the corresponding domain
		 **/
		virtual ContinuationDomain* contractOut(Function& f);

		/**
		 * 	\brief Empty intersection test between the domain and a box.
		 **/
		virtual bool not_intersects(const IntervalVector& b) const;
		
		/**
		 * 	\brief Tests whether the domain contains a box.
		 **/
		virtual bool is_superset(const IntervalVector& b) const;
		
		/**
		 * 	\brief Builds an interval hull of the domain.
		 **/
		virtual IntervalVector hull() const;
		
		virtual void print() const;
	
	
		/** \brief the box **/
		IntervalVector x;
		
		/** \brief the VarSet identifying the "crossed" parameter **/
		VarSet vs;
		
		/** \brief the direction of the manifold within the box **/
		bool sign;
		
	private:
		
		/**
		 * 	\brief Used by contractIn and contractOut
		 **/
		ContinuationDomain* contract(Function& f, bool out = true);
	
}; // ContinuationDomainBox

inline bool ContinuationDomainBox::not_intersects(const IntervalVector& b) const
{
	return !x.intersects(b);
}

inline bool ContinuationDomainBox::is_superset(const IntervalVector& b) const
{
	return x.is_superset(b);
}

inline IntervalVector ContinuationDomainBox::hull() const
{
	return x;
}


/**
 * 	\brief Class implementing the Parallelotope continuation domain. 
 * 
 * 	Used for implementing the ParCont continuation algorithm. 
 *  Contrary to the box domain, we consider here that parallelotopes can
 *  be build so as to implicitly store the parameter "crossed" by the manifold
 *  and its direction.
 **/ 
class ContinuationDomainParallelotope : public ContinuationDomain
{
	public:
	
		/**
		 *	\brief Basic constructor from a Parallelotope.
		 **/
		ContinuationDomainParallelotope(const Parallelotope& p)
		: 	x(p)			
		{}
		
		/**
		 * \brief Certify the existence and unicity of a piece
		 * of manifold of solution to f = 0.
		 **/
		virtual bool certify(Function& f);
		
		/**
		 * 	\brief Contract on the side of the domain where the manifold "enters"
		 *  the domain. Returns the corresponding domain.
		 **/
		virtual ContinuationDomain* contractIn(Function& f);
		
		/**
		 * 	\brief Contract on the side of the domain where the manifold "exits"
		 *  the domain. Returns the corresponding domain
		 **/
		virtual ContinuationDomain* contractOut(Function& f);
		
		/**
		 * 	\brief Empty intersection test between the domain and a box.
		 **/
		virtual bool not_intersects(const IntervalVector& b) const;
		
		/**
		 * 	\brief Tests whether the domain contains a box.
		 **/
		virtual bool is_superset(const IntervalVector& b) const;
		
		/**
		 * 	\brief Builds an interval hull of the domain.
		 **/
		virtual IntervalVector hull() const;
		
		virtual void print() const;
		
		// the Parallelotope
		Parallelotope x;
		
	private:

		/**
		 * 	\brief Used by contractIn and contractOut
		 **/	
		ContinuationDomain* contract(Function& f, bool out = true);
}; // ContinuationDomainParallelotope

inline bool ContinuationDomainParallelotope::not_intersects(const IntervalVector& b) const
{
	return x.not_intersects(b);
}

inline bool ContinuationDomainParallelotope::is_superset(const IntervalVector& b) const
{
	return x.is_superset(b);
}

inline IntervalVector ContinuationDomainParallelotope::hull() const
{
	return x.hull();
}



/**
 * 	\brief Factory for constructing new instances of ContinuationDomain.
 **/
class FactoryContinuationDomain
{
	public:
		
		/**
		 * 	\brief Builds a ContinuationDomainBox for continuation.
		 * 
		 *  The box is constructed from:
		 *  - a vector dir (~ tangent to the manifold), 
		 *  - a point xmid close to the manifold (where to start by default),
		 *  - a step length h
		 *  - an input domain in (that the box must contain if it is not null)
		 *  - a previously computed domain (default nullptr) for using the
		 * 	  domain initialisation heuristic.
		 * 
		 *  The parameter and direction are computed from the vector dir, tangent to the manifold.
		 **/
		static ContinuationDomain* constructDomainBox(const Vector& dir, const Vector& xmid, double h, const ContinuationDomain* in, const ContinuationDomain* previous = nullptr);
		
		/**
		 * 	\brief Builds a ContinuationDomainParallelotope for continuation.
		 * 
		 *  The parallelotope is constructed from:
		 * 	- a vector dir (~ tangent to the manifold),
		 *  - the jacobian of the system (at xmid)
		 *  - a point xmid close to the manifold
		 * 	- an input domain in (that the parallelotope must contain if it is not null)
		 * 	- a previously computed domain (default nullptr) for using the
		 * 	  domain initialisation heuristic.
		 **/
		static ContinuationDomain* constructDomainParallelotope(const Vector& dir, const IntervalMatrix& jac, const Vector& xmid, double h, const ContinuationDomain* in, ContinuationDomain* previous = nullptr);
		
}; // FactoryContinuationDomain


} // namespace ibex
#endif
