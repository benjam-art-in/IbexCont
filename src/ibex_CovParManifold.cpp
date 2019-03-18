/* ============================================================================
 * I B E X 
 * ============================================================================
 * Copyright   : Ecole Polytechnique (FRANCE)
 * License     : LGPL, see the file COPYING.LESSER.
 * Author(s)   : Benjamin Martin
 * Created     : 14.03.2019
 * Last Update : 
 * ---------------------------------------------------------------------------- */

#include <sstream>
#include <algorithm>
#include <limits>

#include "ibex_CovParManifold.h"

namespace ibex
{
	const unsigned int CovParManifold::FORMAT_VERSION = 1;

	const unsigned int CovParManifold::subformat_level = 5;

	const unsigned int CovParManifold::subformat_number = 0;
	
	CovParManifold::CovParManifold(size_t n, size_t m, size_t nb_ineq, ParMode par_mode)
	:	CovManifold(n,m, nb_ineq, CovManifold::EQU_ONLY),
		data(new Data()),
		own_data(true)
	{
		data->_par_manifold_par_mode = par_mode;
	}

	CovParManifold::CovParManifold(const char* filename)
	:	CovParManifold(0,0,0)
	{
		std::stack<unsigned int> format_id;
		std::stack<unsigned int> format_version;
		std::ifstream* f = CovParManifold::read(filename, *this, format_id, format_version);
		f->close();
		delete f;
	}

	CovParManifold::CovParManifold(const Cov& cov, bool copy)
	:	CovManifold(cov, copy)
	{
		const CovParManifold* covpar = dynamic_cast<const CovParManifold*> (&cov);
		
		if(covpar)
		{
			if(copy)
			{
				data = new Data();
				data->_par_manifold_par_mode = covpar->data->_par_manifold_par_mode;
				data->_par_manifold_parallelotopes = covpar->data->_par_manifold_parallelotopes;
				data->_par_manifold_aux_boxes = covpar->data->_par_manifold_aux_boxes;
				data->_par_manifold_centers = covpar->data->_par_manifold_centers;
				data->_par_manifold_transition_matrices = covpar->data->_par_manifold_transition_matrices;
				
				own_data = true;
			}
			else
			{
				data = covpar->data;
				own_data = false;
			}
		}
		else
		{
			data = new Data();
			own_data = true;
		}
	}
	
	CovParManifold::~CovParManifold()
	{
		if(own_data) delete data;
	}
		
	void CovParManifold::add_solution_parallelotope(const Parallelotope& p)
	{
		IntervalVector hp(p.hull());
		
		CovIBUList::add_boundary(p.hull());
		
		// set this box to be unknown for the CovManifold level
		CovManifold::data->_manifold_unknown.push_back(size()-1);
		CovManifold::data->_manifold_status.push_back(UNKNOWN);
		
		data->_par_manifold_parallelotopes.push_back(size()-1);
		data->_par_manifold_aux_boxes.push_back(p.w);
		data->_par_manifold_transition_matrices.push_back(p.C);
		data->_par_manifold_centers.push_back(p.xtilde);
	}
	
	void CovParManifold::add_solution(const IntervalVector& box, const VarSet& varset)
	{
		CovManifold::add_solution(box, box, varset);
	}

	void CovParManifold::save(const char* filename) const
	{
		std::stack<unsigned int> format_id;
		std::stack<unsigned int> format_version;
		std::ofstream* of=CovParManifold::write(filename, *this, format_id, format_version);
		of->close();
		delete of;	
	}
	
	Vector CovParManifold::read_vector(std::ifstream& f, size_t n)
	{
		Vector v(n);
		
		for( unsigned int i = 0; i < n; ++i)
		{
			v[i] = read_double(f);
		}
		
		return v;
	}
	
	void CovParManifold::write_vector(std::ofstream& f, const Vector& v)
	{
		for( unsigned int i = 0; i < v.size(); ++i)
		{
			write_double(f, v[i]);
		}
	}
	
	Matrix CovParManifold::read_sq_matrix(std::ifstream& f, size_t n)
	{
		Matrix m(n,n);
		for(unsigned int i = 0; i < n; ++i)
		{
			for( unsigned int j = 0; j < n; ++j)
			{
				m[i][j] = read_double(f);
			}
		}
		
		return m;
	}
	
	void CovParManifold::write_sq_matrix(std::ofstream& f, const Matrix& m)
	{
		for(unsigned int i = 0; i < m.nb_rows(); ++i)
		{
			for( unsigned int j = 0; j < m.nb_cols(); ++j)
			{
				write_double(f,m[i][j]);
			}
		}
	}
	
	std::ifstream* CovParManifold::read(const char* filename, CovParManifold& cov, std::stack<unsigned int>& format_id, std::stack<unsigned int>& format_version)
	{
		std::ifstream* f = CovManifold::read(filename, cov, format_id, format_version);
		
		size_t nb_parallelotopes;
		
		if(format_id.empty() || format_id.top()!=subformat_number || format_version.top()!=FORMAT_VERSION)
		{
			nb_parallelotopes = 0;
		}
		else
		{
			format_id.pop();
			format_version.pop();
			
			unsigned int par_mode = read_pos_int(*f);

			switch(par_mode)
			{
				case 0  : (ParMode&) cov.data->_par_manifold_par_mode = FULL_MATRIX; break;
				case 1  : (ParMode&) cov.data->_par_manifold_par_mode = NO_MATRIX; break;
				default : ibex_error("[CovParManifold]: unknown storage mode.");
			}
			
			nb_parallelotopes = read_pos_int(*f);
			
			// read indexes
			for(unsigned int i = 0; i < nb_parallelotopes; ++i)
			{
				uint32_t j = read_pos_int(*f);
				
				// check ordering (do we need to check this in ibexcont ?)
				if(!cov.data->_par_manifold_parallelotopes.empty())
				{
					if(j < cov.data->_par_manifold_parallelotopes.back())
						ibex_error("[CovParManifold]: indices of parallelotopes are not in increasing order.");
					if (j==cov.data->_par_manifold_parallelotopes.back())
						ibex_error("[CovParManifold]: duplicated index of parallelotopes.");
				}
				cov.data->_par_manifold_parallelotopes.push_back(j);
			}
			
			// read auxiliary boxes
			for(unsigned int i = 0; i < cov.nb_parallelotopes(); ++i)
			{
				cov.data->_par_manifold_aux_boxes.push_back(read_box(*f, cov.n));
			}
						
			// read centers
			for(unsigned int i = 0; i < cov.nb_parallelotopes(); ++i)
			{
				cov.data->_par_manifold_centers.push_back(read_vector(*f, cov.n));
			}
						
			// read matrices if necessary
			if(cov.data->_par_manifold_par_mode == FULL_MATRIX)
				for(unsigned int i = 0; i < cov.nb_parallelotopes(); ++i)
				{
					cov.data->_par_manifold_transition_matrices.push_back(read_sq_matrix(*f, cov.n));
				}
		}
		
		return f;
	}


	std::ofstream* CovParManifold::write(const char* filename, const CovParManifold& cov, std::stack<unsigned int>& format_id, std::stack<unsigned int>& format_version)
	{
		format_id.push(subformat_number);
		format_version.push(FORMAT_VERSION);

		std::ofstream* f = CovManifold::write(filename, cov, format_id, format_version);
		
		switch(cov.par_mode())
		{
			case FULL_MATRIX : write_pos_int(*f, (uint32_t) 0); break;
			case NO_MATRIX   : write_pos_int(*f, (uint32_t) 1); break;
			default          : assert(false);
		}
		
		write_pos_int(*f, cov.nb_parallelotopes());
		for(auto i : cov.data->_par_manifold_parallelotopes)
		{
			assert(i < std::numeric_limits<uint32_t>::max());
			write_pos_int(*f, (uint32_t) i);
		}

		for(auto b : cov.data->_par_manifold_aux_boxes)
		{
			write_box(*f, b);
		}
		
		for(auto v : cov.data->_par_manifold_centers)
		{
			write_vector(*f, v);
		}
		
		if(cov.data->_par_manifold_par_mode == FULL_MATRIX)
			for(auto m : cov.data->_par_manifold_transition_matrices)
			{
				write_sq_matrix(*f, m);
			}
			
		return f;
	} 

	void CovParManifold::format(std::stringstream& ss, const string& title, std::stack<unsigned int>& format_id, std::stack<unsigned int>& format_version) {
		format_id.push(subformat_number);
		format_version.push(FORMAT_VERSION);

		CovManifold::format(ss, title, format_id, format_version);

		ss
		<< space << " - 1 integer:     the type of storage of parallelotopes\n"
		<< space << "                  - 0 = Transition matrices are stored alongside\n"
		<< space << "                        (FULL_MATRIX mode)\n"
		<< space << "                  - 1 = Transition matrices are not stored.\n"
		<< space << "                        Saves some space but the matrices have\n"
		<< space << "                        to be rebuild carefully\n"
		<< space << "                        (NO_MATRIX mode)\n"
		<< "|   CovParManifold  |" 
		<< space << " - 1 integer:     the number Np of parallelotopes\n"
		<< space << " - Np integers:   the index of parallelotopes\n"
		<< space << " - Np boxes :     the auxiliary boxes of the parallelotopes\n"
		<< space << " - Np vectors:    the centers o\n"
		<< space << " +----[if mode = FULL_MATRIX]----\n"
		<< space << " - Np matrices:   the transition matrices (if FULL_MATRIX mode).\n"
		<< separator;
	}

	string CovParManifold::format() {
		std::stringstream ss;
		std::stack<unsigned int> format_id;
		std::stack<unsigned int> format_version;
		format(ss, "CovParManifold", format_id, format_version);
		return ss.str();
	}
	
	
	std::ostream& operator<<(std::ostream& os, const CovParManifold& cov)
	{
		os << " mode " << cov.par_mode() << std::endl;
		for(int i = 0; i < cov.nb_parallelotopes(); ++i)
		{
			os << " n° " << (cov.id_parallelotope(i)+1) << " is a parallelotope." << std::endl;
		}
		
		return os;
	}
	
} // namespace ibex
