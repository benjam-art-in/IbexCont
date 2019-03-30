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
#include "ibex_Function.h"
#include "ibex_ContinuationDomain.h"
#include "ibex_CovContinuation.h"
#include "ibex_Timer.h"

namespace ibex {


/**
 *	\brief Implementation of the ParCont continuation algorithm.
 *
 *  The implementation follows:
 *	> Benjamin Martin, Alexandre Goldsztejn, Laurent Granvilliers and Christophe Jermann,
 *	> Certified Parallelotope Continuation for One-Manifolds,
 *	> SIAM Journal on Numerical Analysis, Volume 51(6), Pages 3373-3401, 2013 
 *
 *	The same class can be used for BoxCont. The two methods operates
 * 	the same algorithm only with small changes.
 **/
class ContinuationSolver
{
	public:

		/**
		 *	\brief Status of the solver.
		 *	- TANGENT_FAILURE: failure when computing a tangent direction.
		 *	- BACKTRACK: the recently constructed domains goes backward.
		 *	(may indicate an error in tangent computation)
		 *	- LOW_STEP: stoped due to the step size reaching the minimal value.
		 *  - MAX_STEP_NUMBER: stoped due to exceeding the maximal number of steps
		 *  - TIME_OUT: stoped due to time out
		 *	- EXITS_UNIVERSE: the solutions have been successfully followed and crossed
		 * 	the boundary of the 'universe' box.
		 * 	- LOOP: the manifold is periodic within the universe.
		 **/
		enum class Status { TANGENT_FAILURE, BACKTRACK, LOW_STEP, MAX_STEP_NUMBER, TIME_OUT, EXITS_UNIVERSE, LOOP }; 

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
				void changeSign();
				
				/**
				 * 	\brief returns the previously evaluated jacobian.
				 **/
				const IntervalMatrix& getLastJacobian() const;
				
				/**
				 * 	\brief returns the previously evaluated tangent.
				 **/
				const Vector& getLastTangent() const;
				
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
		 * \brief Destructor.
		 **/
		~ContinuationSolver();
		
		/**
		 * 	\brief Execution of the continuation algorithm
		 **/
		void solve(const Vector& init, double hstart = 1.0);
		
		
		
		/**
		 * 	\brief Resets the outputs of the solver.
		 **/
		void flush();
		
		
		/**
		 * \brief Display in the prompt the current state of the solver
		 **/
		void report();
		
		/**
		 * \brief get the number of (succesfully) produced domains
		 **/
		size_t get_nb_domains() const;
		
		void save_cov(const char* filename) const;
		
		void load_cov(const char* filename);
		
	protected:
		
		void init_solving(size_t n);
		
		void stop_solving(Status stat);
		
		void check_time();
		
		void check_domain_limit();
		
		/** \brief The n-1 x n system of equations. **/
		Function& equations;
		
		/** \brief The boundary box. **/
		IntervalVector universe;
		
		/** \brief flag for using boxcont. **/
		bool boxcont;
		
		/** \brief minimal step length. **/
		double hmin;
		
		/** \brief step length decreasing factor. **/
		double alpha;
		
		/** \brief step length incresing factor. **/
		double beta;
		
		/** \brief flag for using the domain initialisation heuristic. **/
		bool flag_heuristic_init;
		
		/** \brief elapsed time. **/
		double time;
		
		/** \brief timer. **/
		Timer timer;
		
		/** \brief number of iterations. **/
		unsigned int nb_iterations;
		
		/** \brief number of ContinuationDomain objects successfully produced. **/
		unsigned int nb_domains;
		
		/** \brief number of connected components. **/
		unsigned int nb_components;
		
		/** \brief status of the solver. **/
		Status solving_status;
		
	public:
		
		/** \brief time limitation (-1 by default, meaning no limit) **/
		double time_limit;
		
		/** \brief limit on the number of ContinuationDomain objects (0 by default, meaning no limit) **/
		unsigned int domain_limit;
		
		// For the output: the domains covering the manifold and their 
		// exit sides.
		std::vector<ContinuationDomain*> manifold;
		std::vector<ContinuationDomain*> checkpoints;
		
		CovContinuation* cov;
		
}; // ContinuationSolver

inline void ContinuationSolver::ContinuationDirectionBuilder::changeSign()
{
	sign = !sign;
}

inline const IntervalMatrix& ContinuationSolver::ContinuationDirectionBuilder::getLastJacobian() const
{
	return lastJacobian;
}

inline const Vector& ContinuationSolver::ContinuationDirectionBuilder::getLastTangent() const
{
	return lastTangent;
}


inline size_t ContinuationSolver::get_nb_domains() const
{
	return nb_domains;
}

inline void ContinuationSolver::save_cov(const char* filename) const
{
	cov->save(filename);
}

inline void ContinuationSolver::load_cov(const char* filename)
{
	if(cov) delete cov;
	cov = new CovContinuation(filename);
}

} // namespace ibex

#endif /* __IBEX_CONTINUATION_SOLVER_H__ */
