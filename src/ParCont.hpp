/* ============================================================================
 * I B E X 
 * ============================================================================
 * Copyright   : Ecole Polytechnique (FRANCE)
 * License     : LGPL, see the file COPYING.LESSER.
 * Author(s)   : Benjamin Martin
 * Created     : 04.02.2019
 * Last Update : 
 * ---------------------------------------------------------------------------- */

#ifndef __IBEX_CONTINUATION_SOLVER_H__
#define __IBEX_CONTINUATION_SOLVER_H__

#include <vector>
#include <ibex.h>

#include "ContinuationDomain.hpp"

namespace ibex {


/**
 *	\brief Implementation of the ParCont continuation algorithm.
 *
 *  The implementation follows:
 *	Benjamin Martin, Alexandre Goldsztejn, Laurent Granvilliers and Christophe Jermann,
 *	Certified Parallelotope Continuation for One-Manifolds,
 *	SIAM Journal on Numerical Analysis, Volume 51(6), Pages 3373-3401, 2013
 *
 *	The same class can be used for BoxCont. The two methods operates
 * 	the same algorithm only with small changes.
 **/
class ContinuationSolver
{
	public:


		/**
		 * 	\brief Sublcass for managing the construction of tangent vectors.
		 * 
		 *  \note Should it be an external class ?
		 **/
		class ContinuationDirectionBuilder
		{
			public:
				
				/**
				 * 	\brief Basic constructor taking the dimension as input.
				 **/
				ContinuationDirectionBuilder(int n)
				:	sign(true),
					lastJacobian(n-1,n),
					lastTangent(n)
				{}
				
				/**
				 * 	\brief Builds a vector tangent to eq = 0 at p.
				 **/
				bool buildNewDirection(Function& eq, const IntervalVector& p);
				
				/**
				 * 	\brief Change the direction of continuation.
				 **/
				void changeSign()
				{
					sign = !sign;
				}
				
				/**
				 * 	\brief returns the previously evaluated jacobian.
				 **/
				const IntervalMatrix& getLastJacobian() const
				{
					return lastJacobian;
				}
				
				/**
				 * 	\brief returns the previously evaluated tangent.
				 **/
				const Vector& getLastTangent() const
				{
					return lastTangent;
				}
				
				// direction of continuation
				bool sign;
				
				// the last evaluated jacobian
				IntervalMatrix lastJacobian;
				
				// the last constructed tangent vector
				Vector lastTangent;
		}; // ContinuationDirectionBuilder

		// default parameters
		static double default_hmin;
		static double default_alpha;
		static double default_beta;
		
		/**
		 * 	\brief Constructor for the continuation solver on the system
		 * 	f = 0 and universe box u.
		 * 
		 * 	boxcont (optional) indicates to use BoxCont.
		 **/
		ContinuationSolver(	Function& f, 
							const IntervalVector& u, 
							bool boxcont=false,
							double hmin = default_hmin,
							double alpha = default_alpha,
							double beta = default_beta);
	
		/**
		 * 	\brief Execution of the continuation algorithm
		 **/
		void solve(const Vector& init, double hstart = 1.0);
		
		/**
		 * 	\brief Resets the outputs of the solver.
		 **/
		void reset();
		
		
		// The n-1 x n system of equations
		Function& equations;
		
		// The boundary box
		IntervalVector universe;
		
		// flag for using boxcont
		bool boxcont;
		
		// minimal step length
		double hmin;
		
		// step length decresing factor
		double alpha;
		
		// step length incresing factor
		double beta;
		
		// flag for using the domain initialisation heuristic
		bool flag_heuristic_init;
		
		// For the output: the domains covering the manifold and their 
		// exit sides.
		std::vector<ContinuationDomain*> manifold_out;
		std::vector<ContinuationDomain*> checkpoints_out;
		
		

}; // ContinuationSolver

} // namespace ibex

#endif /* __IBEX_CONTINUATION_SOLVER_H__ */
