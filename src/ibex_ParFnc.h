/* ============================================================================
 * I B E X 
 * ============================================================================
 * Copyright   : Ecole Polytechnique (FRANCE)
 * License     : LGPL, see the file COPYING.LESSER.
 * Author(s)   : Benjamin Martin
 * Created     : 04.02.2019
 * Last Update : 
 * ---------------------------------------------------------------------------- */
 
#ifndef __IBEX_PARFNC_H__
#define __IBEX_PARFNC_H__

#include <ibex_Fnc.h>
#include "ibex_Parallelotope.h"

namespace ibex {

/**
 * 	\brief Fnc class implementing the evaluation of a function according to a Parallelotope.
 * 
 *  Considers that the input box to evaluate lies in the auxiliary space of some parallelotope.
 **/
class ParFnc : public Fnc {
	public:

		/**
		 * 	\brief Basic constructor
		 **/
		ParFnc(Fnc& f, const Parallelotope& p) 
		: 	Fnc(f.nb_var(),f.image_dim()), 
			f(f), 
			p(p) 
		{ }

		virtual Interval eval(const IntervalVector& w) const 
		{	
			return f.eval(p.xtilde + p.C*w);
		}

		virtual IntervalVector eval_vector(const IntervalVector& w, const BitSet& components) const 
		{	
			return f.eval_vector(p.xtilde + p.C*w, components);
		}

		virtual void gradient(const IntervalVector& w, IntervalVector& g) const
		{
			f.gradient(p.xtilde + p.C*w, g);
			g = g*p.C;	
		}
		
		virtual void jacobian(const IntervalVector& w, IntervalMatrix& J, const BitSet& components, int v=-1) const
		{				
			f.jacobian(p.xtilde + p.C*w, J, components, v);
			J *= p.C;
		}


		// the function
		Fnc& f;
		
		// the parallelotope
		const Parallelotope& p;	
}; // ParFnc

} // ibex

#endif /* __IBEX_PARFNC_H__ */
