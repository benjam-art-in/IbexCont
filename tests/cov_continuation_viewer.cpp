#include <stdio.h>
#include <ibex_CovContinuation.h>

#include "vibes.h"

using namespace ibex;
using namespace std;



#define COL_OUTER "black[#ff413c]"
#define COL_UNKNOWN "black[]"
#define COL_INNER "black[#006dff]"
#define COL_BOUNDARY "black[#f5dc00]"
#define COL_SOLUTION "black[#c869ff]"


void vibes_draw_paral(const Parallelotope& p, string s)
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

void vibes_cov_continuation_viewer(const char* fname)
{
	CovContinuation cov(fname);
	
	const string str_fname(fname);
	
	vibes::beginDrawing();
	vibes::newFigure(str_fname);
	
	// TOOD: assumes a CovContinuation has been correctly loaded.
	for(unsigned int i = 0; i < cov.size(); ++i)
	{
		if(cov.is_box(i))
		{
			vibes::drawBox(cov[i], COL_SOLUTION);
		}
	}
	
	for(unsigned int i = 0; i < cov.nb_parallelotopes(); ++i)
	{
		vibes_draw_paral(cov.get_parallelotope(i), COL_SOLUTION);
	}
	
	vibes::axisAuto(str_fname);
	vibes::endDrawing();
}




int main(int argc, char** argv) {
	if(argc < 2)
	{
		cout << "Usage: covviewer [covfile]. Opens vibes viewer before usage." << endl;
		return 0;
	}
	
	cout << "Open the vibes viewer before usage." << endl;
	vibes_cov_continuation_viewer(argv[1]);
	return 0;
}
