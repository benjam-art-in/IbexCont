/* ============================================================================
 * I B E X 
 * ============================================================================
 * Copyright   : Ecole Polytechnique (FRANCE)
 * License     : LGPL, see the file COPYING.LESSER.
 * Author(s)   : Benjamin Martin
 * Created     : 04.02.2019
 * Last Update : 
 * ---------------------------------------------------------------------------- */

//#include <ibex.h>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <ctime>

#include "ibex_ContinuationSolver.h"
#include "ibex_ContinuationDomain.h"
#include "ibex_Parallelotope.h"
#include "ibex_ParFnc.h"


//#define __USE_VIBES__

#ifdef __USE_VIBES__
#include "vibes.h"
#endif

using namespace std;
using namespace ibex;

/** Prints and vibes draw **/

#ifdef __USE_VIBES__
void draw_paral(const Parallelotope& p, string s = "r[]")
{ 
	Vector p1(2), p2(2), p3(2), p4(2);
	
	p1[0] = p.w[0].lb(); p1[1] = p.w[1].lb();
	p2[0] = p.w[0].lb(); p2[1] = p.w[1].ub();
	p3[0] = p.w[0].ub(); p3[1] = p.w[1].ub();
	p4[0] = p.w[0].ub(); p4[1] = p.w[1].lb();
	
	p1= (p.C*p1 + p.xtilde);
	p2= (p.C*p2 + p.xtilde);
	p3 =  (p.C*p3 + p.xtilde); 
	p4 = (p.C*p4 + p.xtilde);
	
	std::vector<double> x(4),y(4);
	x[0] = p1[0]; x[1] = p2[0]; x[2] = p3[0]; x[3] = p4[0];
	y[0] = p1[1]; y[1] = p2[1]; y[2] = p3[1]; y[3] = p4[1];
	vibes::drawPolygon(x,y,s);
}

void draw_parals(const vector<ContinuationDomain*> v, const string& c = "r[]")
{
	for(auto p : v)
	{
		ContinuationDomainParallelotope* dp = (ContinuationDomainParallelotope*) p;
		draw_paral(dp->x, c);
	}
}

void draw_boxes(const vector<ContinuationDomain*> v, const string& c = "r[]")
{
	for(auto p : v)
	{
		ContinuationDomainBox* db = (ContinuationDomainBox*) p;
		vibes::drawBox(db->x, c);
	}
}
#endif

void print_paral(const Parallelotope& p)
{
	Vector p1(2), p2(2), p3(2), p4(2);
	
	p1[0] = p.w[0].lb(); p1[1] = p.w[1].lb();
	p2[0] = p.w[0].lb(); p2[1] = p.w[1].ub();
	p3[0] = p.w[0].ub(); p3[1] = p.w[1].ub();
	p4[0] = p.w[0].ub(); p4[1] = p.w[1].lb();
	
	cout << (p.C*p1 + p.xtilde) << ", " << (p.C*p2 + p.xtilde) << ", " << (p.C*p3 + p.xtilde) << ", " << (p.C*p4 + p.xtilde) << endl;
}

void print_paral_mathematica(const Parallelotope& p)
{
	Vector p1(2), p2(2), p3(2), p4(2);
	
	p1[0] = p.w[0].lb(); p1[1] = p.w[1].lb();
	p2[0] = p.w[0].lb(); p2[1] = p.w[1].ub();
	p3[0] = p.w[0].ub(); p3[1] = p.w[1].ub();
	p4[0] = p.w[0].ub(); p4[1] = p.w[1].lb();
	
	p1= (p.C*p1 + p.xtilde);
	p2= (p.C*p2 + p.xtilde);
	p3 =  (p.C*p3 + p.xtilde); 
	p4 = (p.C*p4 + p.xtilde);
	
	cout << "{";
	cout << "{" << p1[0] << ","<< p1[1] << "},";
	cout << "{" << p2[0] << ","<< p2[1] << "},";
	cout << "{" << p3[0] << ","<< p3[1] << "},";
	cout << "{" << p4[0] << ","<< p4[1] << "}";
	cout << "}" << endl;
}

void print_parals_mathematica(const vector<ContinuationDomainParallelotope*> manifold)
{
	cout << "{";
	for(auto m : manifold)
	{
		print_paral_mathematica(m->x);
		cout << ",";
	}
	cout << "}";
}

void print_box_mathematica(const IntervalVector& b)
{
	// to use with Rectangle
	cout << "{";
	cout << "{" << b[0].lb() <<", "<< b[1].lb()<< "},{"<< b[0].ub() <<", "<< b[1].ub()<<"}";
	cout << "}" << endl; 
}

void print_boxes_mathematica(const vector<ContinuationDomainParallelotope*> manifold)
{
	cout << "{";
	for(auto m : manifold)
	{
		print_box_mathematica(m->hull());
		cout << ",";
	}
	cout << "}";
}


void print_box3D_mathematica(const IntervalVector& b)
{
	// to use with Cuboid
	cout << "{";
	cout << "{" << b[0].lb() <<", "<< b[1].lb()<< "," << b[2].lb() << "},{"<< b[0].ub() <<", "<< b[1].ub()<< "," << b[2].ub() <<"}";
	cout << "}" << endl; 
}

void print_boxes3D_mathematica(const vector<ContinuationDomain*> manifold)
{	
	cout << "{";
	for(auto m : manifold)
	{
		print_box3D_mathematica(m->hull());
		cout << ",";
	}
	cout << "}";
}

void print_paral3D_mathematica(const Parallelotope& p)
{
	// Matrix, box, center form
	cout << "{";
	
	// Matrix
	cout <<"{";
	for(int i = 0; i < p.size(); ++i)
	{
		cout << "{" << p.C[i][0];
		for(int j = 1; j < p.size(); ++j)
		{
			cout << "," << p.C[i][j];
		}
		cout << "}";
		if (i != p.size()-1)
			cout << ",";
	}
	cout << "},";
	
	// box (cuboid form)
	IntervalVector b(p.w);
	cout << "{";
	cout << "{" << b[0].lb() <<", "<< b[1].lb()<< "," << b[2].lb() << "},{"<< b[0].ub() <<", "<< b[1].ub()<< "," << b[2].ub() <<"}";
	cout << "},"; 
	
	// center
	cout << "{" << p.xtilde[0] << "," << p.xtilde[1] << "," << p.xtilde[2] << "}";
	cout << "}" << endl;
}

void print_parals3D_mathematica(const vector<ContinuationDomain*> manifold)
{	
	cout << "{";
	for(auto m : manifold)
	{
		ContinuationDomainParallelotope* dp = (ContinuationDomainParallelotope*) m;
		print_paral3D_mathematica(dp->x);
		cout << ",";
	}
	cout << "}";
}


/** Test the continuation solver **/

void test_cont_solver()
{
	Function f("x","y","x^2 + y^2 -1");
	
	Vector xtilde(2);
	xtilde[0]=0;
	xtilde[1]=1;
	
	IntervalVector universe(2, Interval(-2,2));
	
	ContinuationSolver solver(f, universe,false);
	
	solver.solve(xtilde);
	
	#ifdef __USE_VIBES__
	vibes::beginDrawing();
	vibes::newFigure("ParCont");
	vibes::drawCircle(0.0,0.0,1.0,"k");
	draw_parals(solver.manifold);
	vibes::endDrawing();
	#endif
}

void test_cont_flower()
{
	
	double eps = 0.999;
	Variable x("x"),y("y");
	Function f(x,y,
		pow(x,8) -(1-eps)*pow(x,6) + 4*pow(x,6)*pow(y,2) - (3+15*eps)*pow(x,4)*pow(y,2) + 6*pow(x,4)*pow(y,4) 
	- (3-15*eps)*pow(x,2)*pow(y,4) + pow(y,8) -(1+eps)*pow(y,6) + 4*pow(y,6)*pow(x,2)
	);
	
	Vector xtilde(2);
	xtilde[0]=0;
	xtilde[1]=sqrt(eps+1);
	
	xtilde[0] = sqrt(1-eps);
	xtilde[1] = 0;
	
	IntervalVector universe(2, Interval(-2,2));
	
	bool box = false;
	for(int i = 0; i < 2; ++i)
	{
		if(box)
			cout << "BoxCont" << endl;
		else
			cout << "ParCont" << endl;
			
		ContinuationSolver solver(f, universe,box);
		
		
		clock_t start = clock();
		solver.solve(xtilde);
		cout << "Solving during " << ((double)(clock()-start))/CLOCKS_PER_SEC << " s" << endl;
		cout << "Produced " << solver.manifold.size() << " elements" << endl;
		
		#ifdef __USE_VIBES__
		vibes::beginDrawing();
		if(box)
		{
			//vibes::newFigure("BoxCont Flower");
			draw_boxes(solver.manifold, "b[]");
		}
		else
		{
			vibes::newFigure("ParCont Flower");
			draw_parals(solver.manifold, "r[]");
		}
		vibes::endDrawing();
		#endif
		
		box = !box;
	}
}

void test_cont_heart()
{
	double eps = 0.000005;
	Variable x("x"),y("y");
	Function f(x,y,
		pow(sqr(x)+sqr(y)-1,3)-sqr(x)*pow(y,3)-eps
		//(sqr(x)+sqr(y)-1)*(sqr(x-eps)+sqr(y)-1)*(sqr(x+eps)+sqr(y)-1)-sqr(x)*pow(y,3)
	);
	
	Vector xtilde(2);
	xtilde[0]=0;
	xtilde[1]=1.1;
	
	IntervalVector universe(2, Interval(-2,2));
	
	bool box = false;
	for(int i = 0; i < 2; ++i)
	{
		if(box)
			cout << "BoxCont" << endl;
		else
			cout << "ParCont" << endl;
			
		ContinuationSolver solver(f, universe,box);
		
		
		clock_t start = clock();
		solver.solve(xtilde);
		cout << "Solving during " << ((double)(clock()-start))/CLOCKS_PER_SEC << " s" << endl;
		cout << "Produced " << solver.manifold.size() << " elements" << endl;
		
		#ifdef __USE_VIBES__
		vibes::beginDrawing();
		if(box)
		{
			//vibes::newFigure("BoxCont Heart");
			draw_boxes(solver.manifold, "b[]");
		}
		else
		{
			vibes::newFigure("ParCont Heart");
			draw_parals(solver.manifold, "r[]");
		}
		vibes::endDrawing();
		#endif
		
		box = !box;
	}
}

void test_cont_condition()
{
	Variable x(3,"x");
	double eps = 1.0/50.0;
	Function f(x, Return(sqr(x[0] + eps) + sqr(x[1]) + sqr(x[2]) - (1+eps*eps), sqr(x[0] - eps) + sqr(x[1]) + sqr(x[2]) - (1+eps*eps)));
	
	Vector xtilde(3);
	xtilde[0] = 0; xtilde[1] = 1; xtilde[2] = 0; 
	
	IntervalVector universe(3, Interval(-2,2));
	
	bool box = false;
	for(int i = 0; i < 2; ++i)
	{
		if(box)
			cout << "BoxCont" << endl;
		else
			cout << "ParCont" << endl;
			
		ContinuationSolver solver(f, universe,box);
		
		
		clock_t start = clock();
		solver.solve(xtilde);
		cout << "Solving during " << ((double)(clock()-start))/CLOCKS_PER_SEC << " s" << endl;
		cout << "Produced " << solver.manifold.size() << " elements" << endl;
		
		
		box = !box;
	}
}


void test_cont_dimension()
{
	size_t n = 32;
	double eps = 1;//e-3;
	
	
	Matrix Q = Matrix::eye(n);
	Q[n-1][n-1] = eps;
	Matrix A(n-2,n, 0.0);
	for(int i = 0; i < n-2; ++i)
		A[i][i] = 1;
	
	Matrix M = Matrix::rand(n);
	gram_schmidt(M);
	
	Q = M.transpose() * Q * M;
	A = A*M;
	
	
	Variable x(n, "x");
	Function f(x, Return(x * Q * x -1, A*x));
	cout << M << endl;
	cout << f << endl;
	
	Vector xtilde(n,0.0);
	xtilde[n-2] = 1;
	xtilde = M.transpose()*xtilde;
	cout << xtilde << endl;
	
	IntervalVector universe(n, Interval(-100000, 100000));
	bool box = false;
	for(int i = 0; i < 2; ++i)
	{
		if(box)
			cout << "BoxCont" << endl;
		else
			cout << "ParCont" << endl;
			
		ContinuationSolver solver(f, universe,box);
		
		
		clock_t start = clock();
		solver.solve(xtilde);
		cout << "Solving during " << ((double)(clock()-start))/CLOCKS_PER_SEC << " s" << endl;
		cout << "Produced " << solver.manifold.size() << " elements" << endl;
		
		
		box = !box;
	}
}

void test_cont_param_synthesis_RRRRRR()
{
	Variable x(2,"x");
	Variable u(2,"u");
	
	double p1 = 2.0, p2 = 2.0;
	double l1 = 1.5, l2 = 1.5;
	double a1 = -1.0, a2 = 0.0;
	double b1 = 1.0, b2 = 0.0;
	
	Function f(x,u, Return(sqr(x[0] - (cos(u[0])*p1 + a1)) + sqr(x[1] - (sin(u[0])*p1 + a2)) -l1*l1, sqr(x[0] - (cos(u[1])*p2 + b1)) + sqr(x[1] - (sin(u[1])*p2 + b2)) -l2*l2));
	
	Variable t("t");
	Function xt(t, Return(sin(Interval::TWO_PI*t), cos(3*(Interval::TWO_PI*t)) + 1.5));
	
	
	cout << f << endl;
	cout << xt << endl;
	
	Function F(u,t, f(xt(t),u));
	cout << F << endl;
	
	
	Vector xtilde(3);
	xtilde[0] = 0.6085;
	xtilde[1] = 1.3699;
	xtilde[2] = 0;
	
	IntervalVector universe(3, Interval::TWO_PI * Interval(-1.1,1.1));
	universe[2] = Interval(-1,1);
	
	bool box = false;
	ContinuationSolver solver(F, universe, box);
	clock_t start = clock();
	solver.solve(xtilde);
	cout << "Solving during " << ((double)(clock()-start))/CLOCKS_PER_SEC << " s" << endl;
	cout << "Produced " << solver.manifold.size() << " elements" << endl;

	//*
	if(box)
		print_boxes3D_mathematica(solver.manifold);
	else
		print_parals3D_mathematica(solver.manifold);
	//*/
	
}

void test_normalize_paral()
{
	double eps = 0.999;
	Variable x("x"),y("y");
	Function f(x,y,
		pow(x,8) -(1-eps)*pow(x,6) + 4*pow(x,6)*pow(y,2) - (3+15*eps)*pow(x,4)*pow(y,2) + 6*pow(x,4)*pow(y,4) 
	- (3-15*eps)*pow(x,2)*pow(y,4) + pow(y,8) -(1+eps)*pow(y,6) + 4*pow(y,6)*pow(x,2)
	);
	
	Vector xtilde(2);
	xtilde[0]=0;
	xtilde[1]=sqrt(eps+1);
	
	xtilde[0] = sqrt(1-eps);
	xtilde[1] = 0;
	
	IntervalMatrix J(f.jacobian(xtilde));
	Vector g(J.row(0).mid());
	g = (1.0/ norm(g)) * g;
	Matrix C(2,2);
	C.set_row(0, g);
	
	Matrix K(kernel(J.mid()));
	C.set_row(1, K.row(0));
	
	cout << "J" << endl;
	cout << J << endl;
	cout << "C" << endl;
	cout << C << endl;
	cout << "C-" << endl;
	Matrix InvC(real_inverse(C));
	cout << InvC << endl;
	cout << J*InvC << endl;
	
	IntervalVector w(2, Interval(0.0,0.0));
	w[1] = Interval(0.0, 1e-5);
	
	J = f.jacobian(InvC*w + xtilde)*InvC;
	cout << J << endl;
	
	
}

//~ void test_minibex()
//~ {
	//~ System sys("examples/flower.txt");
	//~ Vector xinit(2);
	//~ xinit[0] = 0; xinit[1] = 1.41;
	
	//~ System sys("examples/ellipsoidal_3_1.0.txt")
	//~ Vector xinit(3);
	
	
	
	//~ // Get the equations
	//~ System eq(sys, System::EQ_ONLY);
//~ }

int main() {

	// test_cont_flower();
	test_cont_heart();
	//	test_cont_condition();
	//	test_cont_dimension();
	//test_cont_param_synthesis_RRRRRR();
	//~ test_normalize_paral();
}
