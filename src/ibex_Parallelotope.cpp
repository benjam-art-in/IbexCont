/* ============================================================================
 * I B E X 
 * ============================================================================
 * Copyright   : Ecole Polytechnique (FRANCE)
 * License     : LGPL, see the file COPYING.LESSER.
 * Author(s)   : Benjamin Martin
 * Created     : 04.02.2019
 * Last Update : 
 * ---------------------------------------------------------------------------- */

#include "ibex_Parallelotope.h"
#include "ibex_Linear.h"

namespace ibex {
	
	
Matrix Parallelotope::transformation(const Matrix& J, const Vector& tan)
{
	Matrix C(J.nb_cols(), J.nb_cols());
	for(int i = 0; i < J.nb_rows(); ++i)
	{
		C.set_row(i, J.row(i));
	}
	C.set_row(J.nb_cols()-1, tan);
	return real_inverse(C); 
}


IntervalMatrix Parallelotope::computeInvMatrix(const Matrix& C)
{
	IntervalMatrix Cinv(C);
	
	//TODO: check that we do not need a copy of C
	neumaier_inverse(Cinv, Cinv);
	return Cinv;
}

IntervalVector Parallelotope::auxiliary_box(int n, double h)
{
	IntervalVector w(n,Interval::ZERO);
	w[n-1]=Interval(0,h);
	return w;	
}



} // namespace ibex
