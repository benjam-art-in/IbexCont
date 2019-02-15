/* ============================================================================
 * I B E X 
 * ============================================================================
 * Copyright   : Ecole Polytechnique (FRANCE)
 * License     : LGPL, see the file COPYING.LESSER.
 * Author(s)   : Benjamin Martin
 * Created     : 04.02.2019
 * Last Update : 
 * ---------------------------------------------------------------------------- */

#include "ContinuationDomain.hpp"
#include "ParFnc.hpp"

namespace ibex {
	
	/** ContinuationDomainBox **/
	
	bool ContinuationDomainBox::certify(Function& f)
	{
		IntervalVector xe(x);
		IntervalVector xu(x);
		
		bool cert = inflating_newton(f, vs, x, xe, xu);
		
		if(cert)
		{
			// if it is certified, we modify the box
			x = xu;
		} 
		
		return cert;
	}
	
	ContinuationDomain* ContinuationDomainBox::contractIn(Function& f)
	{
		return contract(f, false);
	}

	ContinuationDomain* ContinuationDomainBox::contractOut(Function& f)
	{
		return contract(f, true);
	}
	
	ContinuationDomain* ContinuationDomainBox::contract(Function& f, bool out)
	{
		IntervalVector p(vs.param_box(x));
		
		if (sign)
		{
			p[0] = out ? p[0].ub() : p[0].lb();
		}
		else
		{
			p[0] = out ? p[0].lb() : p[0].ub();
		}
		
		IntervalVector xn(x);
		vs.set_param_box(xn, p); 
		
		ContinuationDomainBox* db = new ContinuationDomainBox(xn, vs, sign);
		
		bool cert = newton(f, vs, db->x);
		
		// TODO: check the result of cert ? technically, we only call contract
		// on a certified domain.
		return db;
	}

	void ContinuationDomainBox::print() const
	{
		std::cout << "BoxDomain : " << std::endl;
		std::cout << "Box : " << x << std::endl;
		std::cout << "VarSet : " << vs << std::endl;
		std::cout << "Sign : " << sign << std::endl;
	}

	/** ContinuationDomainParallelotope **/
		
	bool ContinuationDomainParallelotope::certify(Function& f)
	{
		// Creation of the VarSet for parametric Newton
		// NOTE: we could avoid this by having a newton automatically
		// assuming the last components being the parameters
		BitSet bs = BitSet::all(x.size());
		bs.remove(x.size()-1);
		VarSet vs(x.size(), bs);
		
		// The evaluation object for a parallelotope
		ParFnc pf(f, x);
		
		// The output boxes
		IntervalVector xe(x.w);
		IntervalVector xu(x.w);
		
		// Attempt to certify the parallelotope
		bool cert = inflating_newton(pf, vs, x.w, xe, xu);//, 15, 1, 1.1, 0.0);
		
		if (cert)
		{
			// If it is certified, we change the characteristic box of the parallelotope.
			x.w = xu;
		}
		return cert;
	}
	
	ContinuationDomain* ContinuationDomainParallelotope::contractIn(Function& f)
	{
		return contract(f, false);
	}

	ContinuationDomain* ContinuationDomainParallelotope::contractOut(Function& f)
	{
		return contract(f, true);
	}	
	
	ContinuationDomain* ContinuationDomainParallelotope::contract(Function& f, bool out)
	{
		size_t n = x.size();
		// Creation of the VarSet for parametric Newton
		// NOTE: we could avoid this by having a newton auomatically
		// assuming the last components being the parameters
		BitSet bs = BitSet::all(n);
		bs.remove(n-1);
		VarSet vs(n, bs);
		
		// The evaluation object for a parallelotope
		ParFnc pf(f, x);
		
		// The new parallelotope
		IntervalVector wn(x.w);
		wn[n-1] = out ? wn[n-1].ub() : wn[n-1].lb();
		ContinuationDomainParallelotope* dp = new ContinuationDomainParallelotope(Parallelotope(x, wn));
		
		// Attempt to contract the parallelotope
		bool cert = newton(pf, vs, dp->x.w);//, 1e-20, 1e-16);
		
		return dp;
	}
	
	void ContinuationDomainParallelotope::print() const
	{
		std::cout << "ParallelotopeDomain : " << std::endl;
		std::cout << "Transition : " << x.C << std::endl;
		std::cout << "Box : " << x.w << std::endl;
		std::cout << "xtilde : " << x.xtilde << std::endl;
		std::cout << "hull : " << x.hull() << std::endl;
	}
	
	/** FactoryContinuationDomain **/
	
	ContinuationDomain* FactoryContinuationDomain::constructDomainBox(const Vector& dir, const Vector& xmid, double h, const ContinuationDomain* in, const ContinuationDomain* previous)
	{
		size_t n = dir.size();
		IntervalVector x(xmid);
		
		// get the largest dimension contribution
		int maxi = 0;
		double maxd = std::abs(dir[0]);
		for(int i = 1; i < n; ++i)
		{
			if(maxd < std::abs(dir[i]))
			{
				maxi = i;
				maxd = std::abs(dir[i]);
			}
		}
		
		// build the varset
		BitSet bs = BitSet::all(n);
		bs.remove(maxi);
		VarSet vs(n, bs);
		 
		// The box to consider
		bool sign = dir[maxi] >= 0;
		
		x[maxi] = (sign) ? Interval(x[maxi].lb(),x[maxi].lb()+h) : Interval(x[maxi].ub()-h, x[maxi].ub());
		
		// Enforce to contain the last output
		if(in)
			x |= in->hull();
			
		if(previous)
		{
			// heuristic
			IntervalVector h(previous->hull());
			for(int i = 0; i < n; ++i)
			{
				if(i!=maxi)
					x[i] |= h[i];
			}
			
		}
		
		return new ContinuationDomainBox(x, vs, sign);
	}
		
	ContinuationDomain* FactoryContinuationDomain::constructDomainParallelotope(const Vector& dir, const IntervalMatrix& jac, const Vector& xmid, double h, const ContinuationDomain* in, ContinuationDomain* previous)
	{
		Vector xtilde(xmid);
		
		// TEST: normalisation of the rows of the jacobian
		Matrix jmid(jac.mid());
		for(int i = 0; i < jac.nb_rows(); ++i)
		{
			Vector& r = jmid.row(i);
			double norminv = 1.0 / norm(r);
			r *= norminv;
		}
		
		
		Parallelotope p(xtilde, jmid, h, dir);
		
		// Enforce to contain the last output
		if (in)
			p.w |= p.Cinv*(in->hull()-xtilde);
			
		if(previous)
		{	
			// heuristic
			// TODO: check the cast is correct
			ContinuationDomainParallelotope* parprevious = (ContinuationDomainParallelotope*) previous;
			for(int i = 0; i < p.size()-1; ++i)
			{
				p.w[i] |= parprevious->x.w[i];
			}
			
		}
		ContinuationDomain* dom = new ContinuationDomainParallelotope(p);
		
		return dom;
	}


} // namespace ibex
