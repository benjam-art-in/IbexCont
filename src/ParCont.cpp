/* ============================================================================
 * I B E X 
 * ============================================================================
 * Copyright   : Ecole Polytechnique (FRANCE)
 * License     : LGPL, see the file COPYING.LESSER.
 * Author(s)   : Benjamin Martin
 * Created     : 04.02.2019
 * Last Update : 
 * ---------------------------------------------------------------------------- */

#include <stdlib.h>
#include <vector>

#include "ParCont.hpp"
#include "ContinuationDomain.hpp"

namespace ibex {
	
	
namespace 
{
	// Helper function for constructing the next domain.
	ContinuationDomain* buildNewContinuationDomain(bool boxcont, const ContinuationSolver::ContinuationDirectionBuilder& cbuild,  double h, const Vector& xmid, ContinuationDomain* previousPoint, ContinuationDomain* previousDom = nullptr)
	{
		if (boxcont)
		{
			return FactoryContinuationDomain::constructDomainBox(cbuild.getLastTangent(), xmid, h, previousPoint, previousDom);
		}
		else
		{
			return FactoryContinuationDomain::constructDomainParallelotope(cbuild.getLastTangent(), cbuild.getLastJacobian(), xmid, h, previousPoint, previousDom);
		}
	}
	
	
}

/** ContinuationSolver::ContinuationDirectionBuilder **/

bool ContinuationSolver::ContinuationDirectionBuilder::buildNewDirection(Function& eq, const IntervalVector& p)
{
	// assume that eq.image_dim = eq.nb_var -1
	size_t n = p.size();
	
	// compute the Jacobian of the system.
	// there are extra rows in the matrix for storing the kernel
	lastJacobian = eq.jacobian(p);
	
	// compute the kernel
	Matrix K(kernel(lastJacobian.mid()));
	lastTangent = K.row(0);
	
	IntervalMatrix C(n,n);
	for(int i = 0; i < lastJacobian.nb_rows(); ++i)
	{
		C.set_row(i,lastJacobian.row(i));
	}
	
	C.set_row(n-1, K.row(0));
	
	
	// check (approximately) the correct direction of the tangent	
	Interval d = det(C);
	
	if (d.contains(0.0))
	{
		// Too close to singularity. The tangent vector is not reliable
		return false;
	}
	
	if ((!sign && d.lb() > 0) || (sign && d.ub() < 0))
	{
		// determinant is positive while its sign is set negative
		// or vice-versa. Take the opposite direction
		lastTangent = -lastTangent; 
	} 
	
	return true;
}


/** ContinuationSolver **/

double ContinuationSolver::default_hmin = 1e-9;
double ContinuationSolver::default_alpha = 0.5;
double ContinuationSolver::default_beta = 1.1;

ContinuationSolver::ContinuationSolver(Function& f, const IntervalVector& u, bool boxcont, double hmin, double alpha, double beta)
:	equations(f),
	universe(u),
	boxcont(boxcont),
	hmin(hmin),
	alpha(alpha),
	beta(beta),
	flag_heuristic_init(false)
{}

void ContinuationSolver::solve(const Vector& init, double hstart)
{
	//~ std::cout <<"ParCont solver" << std::endl;
	//~ std::cout << "alpha " << alpha << std::endl;
	//~ std::cout << "beta " << beta << std::endl;
	//~ std::cout << "hmin " << hmin << std::endl;
	size_t n = equations.nb_var();

	ContinuationDirectionBuilder contbuild(n);
		
	size_t k = 1;
	bool stop = false;
	double h = hstart; // Initial step length
	
	std::vector<ContinuationDomain*> checkpoints(0);
	checkpoints.push_back(nullptr);
	std::vector<ContinuationDomain*> manifold(0);
	
	// Note: we count iterations starting from k = 1
	// at the begining of each iteration, manifold is always of size k-1
	// checkpoints is always of size k
	
	
	Vector xtildek(init);
	
	try{
		contbuild.buildNewDirection(equations, xtildek);
	}
	catch(SingularMatrixException& e)
	{
		std::cout << "SingularMatrixException at first tangent construction." << std::endl;
		return;
	}
	
	bool lastFail = true;
	while (! stop)
	{
		//std::cout << "iteration " << k << std::endl;
		
		// Build a new domain
		ContinuationDomain* xcurrent = buildNewContinuationDomain(boxcont, contbuild, h, xtildek, checkpoints[k-1], (flag_heuristic_init && !lastFail) ? manifold[k-2] : nullptr);
		//xcurrent->print();
		// Prove xcurrent contains the manifold
		bool success = false;
		try{
			success = xcurrent->certify(equations);
		}
		catch(SingularMatrixException& e)
		{
			std::cout << "SingularMatrixException when certifying." << std::endl;
			return;
		}
		
		//~ std::cout << "certified ?" << success << std::endl;
		// If so, then contract the "output" side of the parallelotope
		ContinuationDomain* newcheck = nullptr;
		IntervalVector newcheck_hull(n);
		
		bool no_backtrack = true;
		if (success)
		{
			//std::cout << "validate the success" << std::endl;
			// contract Out
			try{
				newcheck = xcurrent->contractOut(equations);
			}
			catch(SingularMatrixException& e)
			{
				std::cout << "SingularMatrixException when contracting out" << std::endl;
				return;
			}
			
			newcheck_hull = newcheck->hull();
			
			// Validate the sucess: 1) No backtrack 2) loop from the begining 3) Crossing the domain or remaining inside
			no_backtrack = (k==1) || (manifold[k-2]->not_intersects(*newcheck) && xcurrent->not_intersects(*checkpoints[k-2]));
			//std::cout << "no_backtrack " << no_backtrack << std::endl;
			bool no_loop = (k==1) || (xcurrent->not_intersects(*checkpoints[0])) || xcurrent->is_superset(*checkpoints[0]);
			//std::cout << "no loop " << no_loop << std::endl;
			//if (k>1)
			//	std::cout << "xk not touch y0 " << (xcurrent->not_intersects(*checkpoints[0])) << " or xk contains y0 " << xcurrent->is_superset(*checkpoints[0]) << std::endl;
			bool no_domain = !(universe.intersects(newcheck_hull)) || universe.is_superset(xcurrent->hull());
			//std::cout << "no domain " << no_domain << std::endl;
			//std::cout << "universe not contains yk " << !(universe.intersects(newcheck_hull)) << " or xk within universe " <<  universe.is_superset(xcurrent->hull()) << std::endl;
			
			// update the success status accordingly
			success = no_backtrack  && no_loop && no_domain;	
		}
		//~ else
		//~ {	
			//~ int t;
			//~ std::cout << "Fail ! h = " << h << ". Abort ?" << std::endl;
			
			//~ std::cin >> t;
			//~ if(t == 1)
			//~ {
				//~ manifold_out = manifold;
				//~ checkpoints_out = checkpoints;
				//~ return;
			//~ }
		//~ }

		
		
		// update the step length and check stopping condition
		if (success)
		{
			//std::cout << "success validated" << std::endl;
			if (k == 1)
			{
				// contract In and considers it as the first checkpoint
				try{
					checkpoints[0] = xcurrent->contractIn(equations);
				}
				catch(SingularMatrixException& e)
				{
					std::cout << "SingularMatrixException when contracting in" << std::endl;
					return;
				}
			}
			manifold.push_back(xcurrent);
			checkpoints.push_back(newcheck);
			h *= beta;
			
			// after a success: , check loop or leaving the domain
			
			bool is_loop = (k>=2) && xcurrent->is_superset(*checkpoints[0]);
			bool is_leaving_universe = ! universe.intersects(newcheck_hull);
			
			// bool is_low_effective_step = distance(checkpoints[k-1]->hull(), newcheck->hull()) <= hmin; // Was necessary for proving termination. In practice ?
			// TODO: if is_leaving_universe, one can restart by reversing the checkpoints and manifold vector and stop at the next is_leaving_universe
			stop = is_loop || is_leaving_universe;
			
			if(!stop)
			{
				lastFail = false;
				k = k+1;
				xtildek = newcheck_hull.mid();
				try{
				contbuild.buildNewDirection(equations,xtildek);
				}
				catch(SingularMatrixException& e)
				{
						std::cout << "singular matrix encountered when computing next tangent: stop" << std::endl;
						break;
				}
			}
		}
		else
		{
			lastFail = true;
			//std::cout << "Failure" << std::endl;
			h *= alpha;
			delete xcurrent;
			if (newcheck)
				delete newcheck;
				
			// after a failure: no backtrack or small step
			// Note that the ContinuationDirectionBuilder is used to avoid cases of backtrack
			bool is_low_step = h <= hmin;
			stop = !(no_backtrack) || is_low_step;
			
		}	
	}
	
	manifold_out = manifold;
	checkpoints_out = checkpoints;
}	

void ContinuationSolver::reset()
{
	for(auto m : manifold_out)
	{
		delete m;
	}
	while(manifold_out.size()>0)
		manifold_out.pop_back();
		
	for(auto m : checkpoints_out)
	{
		delete m;
	}
	while(checkpoints_out.size()>0)
		checkpoints_out.pop_back();
}

} // namespace ibex
