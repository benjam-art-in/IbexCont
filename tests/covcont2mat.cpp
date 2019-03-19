#include <stdio.h>
#include <string.h>
#include <sstream>

#include <ibex_CovContinuation.h>

using namespace ibex;
using namespace std;



void print_double_mat(double d)
{
	stringstream ss;
	ss << d;
	
	string strd(ss.str());
	
	size_t pos = strd.find_first_of("e"); // find e character if existing
	if(pos == string::npos)
	{
		// no match, print regularly
		cout << strd;
	}
	else
	{
		strd.replace(pos,1," 10^");
		cout << strd;
	}
		
}

void print_matrix_mathematica(const Matrix& m)
{
	cout <<"{";
	for(int i = 0; i < m.nb_rows(); ++i)
	{
		cout << "{";
		print_double_mat(m[i][0]);
		for(int j = 1; j < m.nb_cols(); ++j)
		{
			cout << ",";
			print_double_mat(m[i][j]);
		}
		cout << "}";
		if (i != m.nb_rows()-1)
			cout << ",";
	}
	cout << "}";
}

void print_vector_mathematica(const Vector& v)
{
	cout << "{";
	print_double_mat(v[0]);
	for(int i = 1; i < v.size(); ++i)
	{
		cout << ", ";
		print_double_mat(v[i]);
	}
	cout << "}";
}

void print_box_mathematica(const IntervalVector& b)
{
	// to use with Rectangle
	cout << "{";
	cout << "{";
	print_double_mat(b[0].lb());
	cout << ", ";
	print_double_mat(b[1].lb());
	cout << "},{";
	print_double_mat(b[0].ub()); 
	cout <<", ";
	print_double_mat(b[1].ub());
	cout <<"}}" << endl; 
}

void print_box3D_mathematica(const IntervalVector& b)
{
	// to use with Cuboid
	cout << "{";
	cout << "{";
	print_double_mat(b[0].lb());
	cout << ", ";
	print_double_mat(b[1].lb());
	cout << ", ";
	print_double_mat(b[2].lb());
	cout << "},{";
	print_double_mat(b[0].ub()); 
	cout <<", ";
	print_double_mat(b[1].ub());
	cout <<", ";
	print_double_mat(b[2].ub());
	cout <<"}}" << endl; 
}

void print_paral_mathematica(const Parallelotope& p)
{
	// Matrix, box, center form
	cout << "{";
	
	// Matrix
	print_matrix_mathematica(p.C);
	cout << ", ";
	
	// box (Rectangle form)
	print_box_mathematica(p.w);
	cout << ", ";
	
	// center
	print_vector_mathematica(p.xtilde);
	cout << "}" << endl;
}


void print_paral3D_mathematica(const Parallelotope& p)
{
	// Matrix, box, center form
	cout << "{";
	
	// Matrix
	print_matrix_mathematica(p.C);
	cout << ", ";
	
	// box (cuboid form)
	print_box3D_mathematica(p.w);
	cout << ", ";
	
	// center
	print_vector_mathematica(p.xtilde);
	cout << "}" << endl;
}


void cov2mat(const char* filename)
{
	CovContinuation cov(filename);
	
	bool mod3D;
	if(cov.n == 2) mod3D = false;
	else if (cov.n == 3) mod3D = true;
	else
	{
		cerr << "covcont2mat error. Dimension of elements in the cov file larger than 3." << endl;
	}
	
	// TODO: check we have correctly loaded a CovContinuation file
	
	bool first_box = true;
	for(unsigned int i = 0; i < cov.size(); ++i)
	{
		if(cov.is_box(i))
		{
			if(first_box)
			{
				cout << "lb = {";
				first_box = false;
			}
			else
			{
				cout << ", ";
			}
			if(mod3D) print_box3D_mathematica(cov[i]);
			else print_box_mathematica(cov[i]);
		}
	}
	if(!first_box) cout << "};" << endl; 
	
	bool first_paral = true;
	for(unsigned int i = 0; i < cov.nb_parallelotopes(); ++i)
	{
		if(first_paral)
		{
			cout << "lp = {";
			first_paral = false;
		}
		else
		{
			cout << ", ";
		}
		if(mod3D) print_paral3D_mathematica(cov.get_parallelotope(i));
		else print_paral_mathematica(cov.get_parallelotope(i));
	}
	if(!first_paral) cout << "};" << endl;
}

int main(int argc, char** argv) {
	if(argc < 2)
	{
		cout << "Usage: covcont2mat [covfile]." << endl;
		return 0;
	}
	
	cov2mat(argv[1]);
	return 0;
}
