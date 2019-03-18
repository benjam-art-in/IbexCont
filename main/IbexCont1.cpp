/* ============================================================================
 *                                I B E X 
 * 
 *                               IbexCont1
 * 
 * 
 * ============================================================================
 * Copyright   : Ecole Polytechnique (FRANCE)
 * License     : LGPL, see the file COPYING.LESSER.
 * Author(s)   : Benjamin Martin
 * Created     : 18.02.2019
 * Last Update : 
 * ---------------------------------------------------------------------------- */

#include <ibex.h>
#include <ctime> // Temporary

#include "args.hxx"
#include "ibex_ContinuationSolver.h"
#include <sstream>

//~ #ifndef _IBEX_WITH_CONT1_
//~ #error "You need to install the IbexCont1 plugin."
//~ #endif

using namespace std;
using namespace ibex;

int main(int argc, char** argv) {
	
	
	args::ArgumentParser parser("********** IbexCont1 **********", "Continuation on a n x n-1 system of equations in Minibex file.");
	args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
	
	// to be deprecated...
	args::ValueFlag<std::string> output_file(parser, "filename", "COV output file. The file will contain the "
			"description of the manifold with boxes in the COV (binary) format. See --format", {'o',"output"});
	args::Flag format(parser, "format", "Give a description of the COV format used by IbexSolve", {"format"});
	
	args::ValueFlag<double> hmin(parser, "float", "Minimal step length of continuation.", {"hmin"});
	args::ValueFlag<double> hstart(parser, "float", "Starting step length of continuation.", {"hstart"});
	args::ValueFlag<double> alpha(parser, "float", "Decreasing step length factor.", {"alpha"});
	args::ValueFlag<double> beta(parser, "float", "Increasing step length factor.", {"beta"});
	args::Flag boxcont(parser, "boxcont", "Use the box continuation algorithm (parallelotope by default)",{"boxcont"});
	

	// TODO Adds verbose, trace, quiet, output, ...
	args::Flag quiet(parser, "quiet", "Print no report on the standard output.",{'q',"quiet"});
	
	args::Positional<std::string> filename(parser, "filename", "The name of the MINIBEX file.");
	args::PositionalList<double> initpoint(parser, "initpoint", "The starting point (blank separated float coordinates). Missing coordinates are treated as 0.");
	
	
	try
	{
		parser.ParseCLI(argc, argv);
	}
	catch (args::Help&)
	{
		std::cout << parser;
		return 0;
	}
	catch (args::ParseError& e)
	{
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}
	catch (args::ValidationError& e)
	{
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}


	if (filename.Get()=="") {
		ibex_error("no input file (try ibexcont1 --help)");
		exit(1);
	}
	
	if(format)
	{
		cout << CovContinuation::format() << endl;
		exit(0);
	}
	
	try
	{
		if (!quiet) 
		{
			cout << endl << "***************************** setup *****************************" << endl;
			cout << "  file loaded:\t\t" << filename.Get() << endl;
		}
		
		// output_file management taken from ibexsolve
		string output_manifold_file;
		bool overwitten=false;       // is it overwritten?
		string manifold_copy;
		
		if (output_file) 
		{
			output_manifold_file = output_file.Get();
		} 
		else 
		{
			// got from stackoverflow.com:
			string::size_type const p(filename.Get().find_last_of('.'));
			// filename without extension
			string filename_no_ext=filename.Get().substr(0, p);
			stringstream ss;
			ss << filename_no_ext << ".cov";
			output_manifold_file=ss.str();

			ifstream file;
			file.open(output_manifold_file.c_str(), ios::in); // to check if it exists

			if (file.is_open()) 
			{
				overwitten = true;
				stringstream ss;
				ss << output_manifold_file << "~";
				manifold_copy=ss.str();
				// got from stackoverflow.com:
				ofstream dest(manifold_copy, ios::binary);

			    istreambuf_iterator<char> begin_source(file);
			    istreambuf_iterator<char> end_source;
			    ostreambuf_iterator<char> begin_dest(dest);
			    copy(begin_source, end_source, begin_dest);
			}
			file.close();
		}
		
		if (!quiet)
		{
			cout << "  output file:\t\t" << output_manifold_file << "\n";
		}
		
		// check parameters
		double min_h = hmin ? hmin.Get() : ContinuationSolver::default_hmin;
		if(min_h <= 0.0)
		{
			// just a warning if hmin is set to zero (or negative).
			std::cout << "\nWarning: hmin has been set to 0 or a negative value.\n";
		}
		
		if(alpha && (alpha.Get() < 0.0 || alpha.Get() >= 1.0))
		{
			// alpha must lie within [0, 1)
			std::cerr << "\nError: alpha must be within [0,1).\n";
			exit(0);
		}
		double alpha_param = (alpha) ? alpha.Get() : ContinuationSolver::default_alpha;
		
		if(beta && (beta.Get() < 1.0))
		{
			// beta must be > 1
			std::cerr << "\nError: beta must be greater than 1.\n";
			exit(0);
		}
		double beta_param = (beta) ? beta.Get() : ContinuationSolver::default_beta;
		
		double start_h = hstart ? hstart.Get() : 1.0;
		if(start_h < min_h)
		{
			std::cerr << "\nError: starting step length smaller than minimal one.\n";
			exit(0);
		}
		
		if (!quiet)
		{
			cout << "  box   :\t\t" << ((boxcont) ? " yes" : " no") << endl;
			cout << "  hstart:\t\t" << start_h << endl;
			cout << "  hmin  :\t\t" << min_h << endl;
			cout << "  alpha :\t\t" << alpha_param << endl;
			cout << "  beta  :\t\t" << beta_param << endl;
		}
				
		// Load the system
		System sys(filename.Get().c_str());
		
		// Get the equations
		System eq(sys, System::EQ_ONLY);
		if(eq.nb_var != eq.nb_ctr +1)
		{
			std::cerr << "\nError: " << "System of equations with " << eq.nb_var << " variables and " << eq.nb_ctr << " equations (it should be n x n-1).\n";
			exit(0); 
		}
		
		// Get the initial point
		Vector init(eq.nb_var, 0.0);
		auto p = initpoint.Get();
		if(p.size() < eq.nb_var)
		{
			std::cout << "\nWarning: " << "Not all coordinates of the initial point have been provided. Missing coordinates are set to 0.\n";
		}
		else if(p.size() > eq.nb_var)
		{
			std::cerr << "\nError: " << "The provided initial point has more coordinates than variables.\n";
			exit(0);
		}
		
		int i =0;
		for(auto pi : p)
		{
			init[i++] = pi;
		}
		//cout << "initial " << init << endl;
		
		if (!quiet)
		{
			cout << "  equations:" <<  eq.f_ctrs << endl;
			cout << "*****************************************************************" << endl << endl;
		}
		
		
		// Create the continuation solver
		ContinuationSolver solver(	eq.f_ctrs,
									sys.box,
									boxcont,
									min_h,
									alpha_param,
									beta_param
									);
									
		// Solve it
		if(!quiet) std::cout << "Solving...";
		clock_t start_time = clock();
		solver.solve(init, start_h);
		if(!quiet) std::cout << " done !\n";
		
		// TODO: output
		if(!quiet)
		{
			solver.report();
		}
		
		solver.save_cov(output_manifold_file.c_str());

		if (!quiet) 
		{
			cout << " results written in " << output_manifold_file << "\n";
			if (overwitten)
				cout << " (old file saved in " << manifold_copy << ")\n";
		}
	}
	catch(ibex::UnknownFileException& e) {
		cerr << "Error: cannot read file '" << filename.Get() << "'" << endl;
	}
	catch(ibex::SyntaxError& e) {
		cout << e << endl;
	}
	catch(ibex::DimException& e) {
		cout << e << endl;
	}
	
}
