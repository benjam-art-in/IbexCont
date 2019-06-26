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

#include "ibex_ContinuationSolver.h"
#include "ibex_Kernel.h"
#include "ibex_Linear.h"

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
	
	void addSuccesfulDomain(CovContinuation* cov, ContinuationDomain* domain, bool boxcont)
	{
		// TODO: better method than dynamic cast ?
		if(boxcont)
		{
			ContinuationDomainBox* dombox = dynamic_cast<ContinuationDomainBox*>(domain);
			if(dombox) cov->add_solution(dombox->x, dombox->vs);
		}
		else
		{
			ContinuationDomainParallelotope* dompar = dynamic_cast<ContinuationDomainParallelotope*>(domain);
			if(dompar) cov->add_solution_parallelotope(dompar->x);
		}
	}
	
	class MaxNumberDomainException : Exception { };
} // anonymous namespace

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
	flag_heuristic_init(false),
	time(0.0),
	timer(),
	nb_iterations(0),
	nb_domains(0),
	nb_components(0),
	time_limit(-1.0), 
	domain_limit(-1),
	solving_status(Status::LOOP), // TODO: have an extra status for initialisation ?
	cov(nullptr)
{}

ContinuationSolver::~ContinuationSolver()
{
	if(cov) delete cov;
	flush();
}

void ContinuationSolver::init_solving(size_t n)
{
	flush();
	// TODO: check whether the initial point belongs to an existing component ?
	if(!cov) cov = new CovContinuation(n, n-1,0);
	
	if(cov) cov->begin_component();
	
	time = 0;
	nb_iterations = 0;
	nb_domains = 0;
	nb_components = 1;
}

void ContinuationSolver::stop_solving(Status stat)
{
	timer.stop();
	time = timer.get_time();
	solving_status = stat;
	
	cov->end_component();
	cov->set_solver_status((unsigned int) stat);
	cov->set_time(cov->get_time() + time);
	cov->set_iterations(cov->get_iterations() + nb_iterations);
}

void ContinuationSolver::check_time()
{
	timer.check(time_limit); // Throw TimeOutException
}

void ContinuationSolver::check_domain_limit()
{
	if (get_nb_domains() > domain_limit)
		throw MaxNumberDomainException();
}


//~ void ContinuationSolver::solvebis(double hstart)
//~ {
	//~ size_t n = equations.nb_var();
	
	//~ init_solving(n);
	
	
	//~ try
	//~ {
		//~ while(has_next_component())
		//~ {
			//~ // get the next component with a free connector
			//~ size_t c = get_next_component();
			//~ try
			//~ {
				//~ continuation(c, hstart);
			//~ }
			//~ catch(SingularMatrixException& e)
			//~ {
				//~ // maybe we can avoid to stop solving and continue with another component
				//~ stop_solving(Status::TANGENT_FAILURE);
				//~ return;
			//~ } // TODO: add an exception in case of backtrack
		//~ }
	//~ }
	//~ catch(TimeOutException& e)
	//~ {
		//~ stop_solving(Status::TIME_OUT);
	//~ }
	//~ catch(MaxNumberDomainException& e)
	//~ {
		//~ stop_solving(Status::MAX_STEP_NUMBER);
	//~ }
//~ }

//~ void ContinuationSolver::continuation(size_t c, double hstart)
//~ {
	//~ size_t n = equations.nb_var();
	
	//~ ContinuationDirectionBuilder contbuild(n);
	
	//~ ContinuationDomain* connector = get_right_connector(c);
	//~ Vector xtildek(connector->hull().mid());
	
	//~ // TODO: the sense of continuation has to be limited
	//~ contbuild.setSign(get_component_sense(c));
	//~ contbuild.buildNewDirection(equations, xtildek); 
	
	//~ bool stop = false;
	//~ double h = hstart;
	
	//~ while(! stop)
	//~ {
		//~ // check stopping conditions (if they are set).
		//~ if(time_limit > 0.0) check_time();
		//~ if(domain_limit > 0) check_domain_limit();
		
		//~ // Build a new domain.
		//~ ContinuationDomain* xcurrent = buildNewContinuationDomain(boxcont, contbuild, h, xtildek, connector, (flag_heuristic_init && !lastFail) ? manifold[k-2] : nullptr);

		//~ bool success = certify_and_check(c, xcurrent);
		
		//~ if(success)
		//~ {
			//~ // successful step
			
			//~ nb_domains++;
			
			//~ // check if continuation is finished
			//~ if(!continuable_component(c))
			//~ {
				//~ return;
			//~ }
			
			//~ h *= beta;
			//~ k = k+1;
			
			//~ connector = get_right_connector(c);
			//~ xtildek = connector->hull().mid();
			
			//~ contbuild.buildNewDirection(equations,xtildek);
		//~ }
		//~ else
		//~ {
			//~ lastFail = true;

			//~ h *= alpha;
			//~ delete xcurrent;
	
			//~ if(h <= hmin)
			//~ {
				//~ // TODO: check singularity / close the right connector
				
				//~ return;
			//~ }
		//~ }
		
		//~ nb_iterations++;
	//~ }
	
//~ }

//~ bool ContinuationSolver::certify_and_check(size_t& c, ContinuationDomain* xc)
//~ {
	//~ assert(c < components.size());
	
	//~ if(xc->certify((Fnc&)equations))
	//~ {
		//~ // validate the success
		
		//~ // check loop
		//~ bool is_loop = false;
		//~ if(! components[c].is_empty())
		//~ {
			//~ ContinuationDomain* lc = components[c].get_left_connector();
			//~ if(! xc->not_intersects(*lc))
			//~ {
				//~ // there is an intersection, check if everything is alright
				//~ is_loop = xc->is_superset(*lc);
				//~ if(!is_loop) // TEST
				//~ {
					//~ // xc is not enclosing lc, force a union and certify again
					//~ xc->join(*lc);
					//~ is_loop = xc->certify((Fnc&)equations);
					//~ if(!is_loop) return false;
				//~ }
			//~ }			
		//~ }
		
		//~ // check join
		
		
		//~ ContinuationDomain* newcheck = xc->contractOut((Fnc&)equations);
		//~ IntervalVector newcheck_hull(newcheck->hull());
		
		
		//~ // check backtrack
		//~ bool is_backtrack = false;
		//~ if(! components[c].is_empty())
		//~ {
			//~ // TODO
		//~ }
		//~ if(is_backtrack) // TODO: thow exception
		
		//~ // check domain
		//~ bool is_domain = ! universe.intersects(newcheck_hull); 
		//~ if(!is_domain && ! universe.is_superset(xcurrent->hull())) return false;
		
		//~ // successful
		//~ components[c].push_right(xc, newcheck);
		//~ if(is_loop) components[c].close_loop();
		//~ else if(is_domain) components[c].close_right();
	//~ }
	//~ else
	//~ {
		//~ return false;
	//~ }
//~ }

void ContinuationSolver::solve(const Vector& init, double hstart)
{
	
	
	//~ std::cout <<"ParCont solver" << std::endl;
	//~ std::cout << "alpha " << alpha << std::endl;
	//~ std::cout << "beta " << beta << std::endl;
	//~ std::cout << "hmin " << hmin << std::endl;
	size_t n = equations.nb_var();
	
	init_solving(n);

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
		flush();
		stop_solving(Status::TANGENT_FAILURE);
		return;
	}
	
	// boolean that sets to false if the last continuation step was a success
	bool lastFail = true;
	
	// if the continuation leaves the universe, restart from the beggining and go the opposite direction
	bool other_direction = false;
	while (! stop)
	{
		// check time
		if(time_limit > 0)
		{
			try
			{
				check_time();
			}
			catch(TimeOutException& e)
			{
				flush();
				stop_solving(Status::TIME_OUT);
				return;
			}
		}
		
		// check domain
		if(domain_limit > 0)
		{
			try
			{
				check_domain_limit();
			}
			catch(MaxNumberDomainException& e)
			{
				flush();
				stop_solving(Status::MAX_STEP_NUMBER);
				return;
			}
		}
			
		// Build a new domain
		ContinuationDomain* xcurrent = buildNewContinuationDomain(boxcont, contbuild, h, xtildek, checkpoints[k-1], (flag_heuristic_init && !lastFail) ? manifold[k-2] : nullptr);

		// Prove xcurrent contains the manifold
		bool success = false;
		try{
			success = xcurrent->certify((Fnc&)equations);
		}
		catch(SingularMatrixException& e)
		{
			std::cout << "SingularMatrixException when certifying." << std::endl;
			flush();
			stop_solving(Status::TANGENT_FAILURE);
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
				newcheck = xcurrent->contractOut((Fnc&)equations);
			}
			catch(SingularMatrixException& e)
			{
				std::cout << "SingularMatrixException when contracting out" << std::endl;
				flush();
				stop_solving(Status::TANGENT_FAILURE);
				return;
			}
			
			newcheck_hull = newcheck->hull();
			
			// Validate the sucess: 1) No backtrack 2) loop from the begining 3) Crossing the domain or remaining inside
			no_backtrack = (k==1) || (manifold[k-2]->not_intersects(*newcheck) && xcurrent->not_intersects(*checkpoints[k-2]));

			//~ std::cout << "loop test " << std::endl;
			bool no_loop = (k==1) || (xcurrent->not_intersects(*checkpoints[0])) || xcurrent->is_superset(*checkpoints[0]);
			if(! no_loop)
			{
					// TODO: try to make the union and certify again !
			}

			//if (k>1)
			bool no_domain = !(universe.intersects(newcheck_hull)) || universe.is_superset(xcurrent->hull());
			
			
			//~ std::cout << no_backtrack << " " << no_loop << " " << no_domain << std::endl;
			//~ std::cout << xcurrent->hull() << std::endl;
			//~ if(k>1) std::cout << checkpoints[0]->hull() << std::endl;
			// update the success status accordingly
			success = no_backtrack  && no_loop && no_domain;	
		}

		
		
		// update the step length and check stopping condition
		if (success)
		{
			//~ std::cout << "sucessful iteration" << std::endl;

			if (k == 1)
			{
				// contract In and considers it as the first checkpoint
				try{
					checkpoints[0] = xcurrent->contractIn((Fnc&)equations);
				}
				catch(SingularMatrixException& e)
				{
					std::cout << "SingularMatrixException when contracting in" << std::endl;
					flush();
					stop_solving(Status::TANGENT_FAILURE);
					return;
				}
			}
			manifold.push_back(xcurrent);
			checkpoints.push_back(newcheck);
			h *= beta;
						
			nb_domains++;
			addSuccesfulDomain(cov, xcurrent, boxcont);
			
			
			// after a success: , check loop or leaving the domain
			
			bool is_loop = (k>=2) && xcurrent->is_superset(*checkpoints[0]);
			bool is_leaving_universe = ! universe.intersects(newcheck_hull);
			
			// bool is_low_effective_step = distance(checkpoints[k-1]->hull(), newcheck->hull()) <= hmin; // Was necessary for proving termination. In practice ?

			stop = is_loop || (is_leaving_universe && other_direction);
			
			
			if(!is_loop && is_leaving_universe && !other_direction)
			{
				other_direction = true;
				// reverse the continuation domains, direction of continuation and restart
				contbuild.switchSign();
				std::reverse(manifold.begin(), manifold.end());
				std::reverse(checkpoints.begin(), checkpoints.end());

				newcheck_hull = checkpoints.back()->hull();
			}
			
			if(!stop)
			{
				// prepare the next iteration
				lastFail = false;
				k = k+1;
				xtildek = newcheck_hull.mid();
				try{
					contbuild.buildNewDirection(equations,xtildek);
				}
				catch(SingularMatrixException& e)
				{
					std::cout << "singular matrix encountered when computing next tangent: stop" << std::endl;
					flush();
					stop_solving(Status::TANGENT_FAILURE);
					return;
				}
			}
			else
			{
				flush();
				// stopping
				if(is_loop)
				{
					stop_solving(Status::LOOP);
				}
				else
				{
					stop_solving(Status::EXITS_UNIVERSE);
				}
				return;
			}
		}
		else
		{
						//~ std::cout << "fail iteration" << std::endl;
						
			lastFail = true;

			h *= alpha;
			delete xcurrent;
			if (newcheck)
				delete newcheck;
				
			// after a failure: no backtrack or small step
			// Note that the ContinuationDirectionBuilder is used to avoid cases of backtrack
			bool is_low_step = h <= hmin;
			stop = !(no_backtrack) || is_low_step;
			
			if(stop)
			{
				// stopping
				flush();
				if(is_low_step)
				{
					stop_solving(Status::LOW_STEP);
				}
				else
				{
					stop_solving(Status::BACKTRACK);
				}
				return;
			}
		}
		
		nb_iterations++;
	}
}	

void ContinuationSolver::flush()
{
	for(auto m : manifold)
	{
		delete m;
	}
	while(manifold.size()>0)
		manifold.pop_back();
		
	for(auto m : checkpoints)
	{
		delete m;
	}
	while(checkpoints.size()>0)
		checkpoints.pop_back();
		
	//~ if(cov) delete cov;
}


// To match ibexsolve report .
namespace {
const char* green() {
#ifndef _WIN32
	return "\033[32m";
#else
	return "";
#endif
}

const char* red(){
#ifndef _WIN32
	return "\033[31m";
#else
	return "";
#endif
}

const char* white() {
#ifndef _WIN32
	return "\033[0m";
#else
	return "";
#endif
}

} // anonymous namespace


void ContinuationSolver::report() 
{

	switch ((Status) cov->solver_status()) 
	{
	case Status::LOOP: 
		std::cout << green() << " solving successfully found a periodic manifold !" << std::endl;
		break;
	case Status::EXITS_UNIVERSE: 
		std::cout << green() << " solving successfully found a manifold traversing the init_box !" << std::endl;
		break;
	case Status::MAX_STEP_NUMBER: 
		std::cout << red() << " number of successulf iterations " << domain_limit << " reached " << std::endl;
		break;
	case Status::TIME_OUT: 
		std::cout << red() << " time limit " << time_limit << "s. reached " << std::endl;
		break;
	case Status::LOW_STEP:
		std::cout << red() << " minimal step length " << hmin << " reached " << std::endl;
		break;
	case Status::BACKTRACK: 
		std::cout << red() << " fail ! Backtracking observed " << std::endl;
		break;
	case Status::TANGENT_FAILURE:
		std::cout << red() << " fail ! Error when computing tangent to jacobian " << std::endl; 
	}

	std::cout << white() << std::endl;

	std::cout << " number of domains   :\t";
	if (cov->size()==0) std::cout << "--"; else std::cout << cov->size();
	std::cout << std::endl;
	std::cout << " number of components:\t";
	if (cov->nb_components()==0) std::cout << "--"; else std::cout << cov->nb_components();
	std::cout << std::endl;
	std::cout << " cpu time used       :\t" << time << "s";
	if (cov->get_time()!=time)
		std::cout << " [total=" << cov->get_time() << "]";
	std::cout << std::endl;
	std::cout << " number of iterations:\t" << nb_iterations;
	if (cov->get_iterations()!=nb_iterations)
		std::cout << " [total=" << cov->get_iterations() << "]";
	std::cout << std::endl << std::endl;
}


} // namespace ibex
