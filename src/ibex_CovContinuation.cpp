/* ============================================================================
 * I B E X 
 * ============================================================================
 * Copyright   : Ecole Polytechnique (FRANCE)
 * License     : LGPL, see the file COPYING.LESSER.
 * Author(s)   : Benjamin Martin
 * Created     : 15.03.2019
 * Last Update : 
 * ---------------------------------------------------------------------------- */

#include <sstream>
#include <algorithm>

#include "ibex_CovContinuation.h"


namespace ibex
{
	const unsigned int CovContinuation::FORMAT_VERSION = 1;
		
	const unsigned int CovContinuation::subformat_level = 6;
	 
	const unsigned int CovContinuation::subformat_number = 0;
	
	CovContinuation::CovContinuation(size_t n, size_t m, size_t nb_ineq)
	:	CovParManifold(n, m, nb_ineq, CovParManifold::FULL_MATRIX),
		data(new Data()),
		own_data(true)
	{}
			
	CovContinuation::CovContinuation(const char* filename)
	:	CovContinuation(0,0,0)
	{
		std::stack<unsigned int> format_id;
		std::stack<unsigned int> format_version;
		std::ifstream* f = CovContinuation::read(filename, *this, format_id, format_version);
		f->close();
		delete f;
	}
	
	CovContinuation::CovContinuation(const Cov& cov, bool copy)
	:	CovParManifold(cov, copy)
	{
		// TODO
	}
	
	CovContinuation::~CovContinuation()
	{
		if(own_data) delete data;
	}
	
	void CovContinuation::add_solution_parallelotope(const Parallelotope& p)
	{
		CovParManifold::add_solution_parallelotope(p);
		
		data->_continuation_solver_domain_type.push_back(PARALLELOTOPE);
	}

	void CovContinuation::add_solution(const IntervalVector& box, const VarSet& varset)
	{
		CovParManifold::add_solution(box, varset);
		
		data->_continuation_solver_domain_type.push_back(BOX);
	}
			
	void CovContinuation::begin_component()
	{
		start_component =size();
	}
			
	void CovContinuation::end_component()
	{
		// no element have been added since the last call to begin_component.
		if(start_component == size())
			return;
			
		data->_continuation_solver_components.push_back(std::make_pair(start_component, size()));
	}
	
	void CovContinuation::save(const char* filename) const
	{
		std::stack<unsigned int> format_id;
		std::stack<unsigned int> format_version;
		std::ofstream* of=CovContinuation::write(filename, *this, format_id, format_version);
		of->close();
		delete of;	
	}
	
	std::ifstream* CovContinuation::read(const char* filename, CovContinuation& cov, std::stack<unsigned int>& format_id, std::stack<unsigned int>& format_version)
	{
		std::ifstream* f = CovParManifold::read(filename, cov, format_id, format_version);
		
		if(format_id.empty() || format_id.top()!=subformat_number || format_version.top()!=FORMAT_VERSION)
		{
			cov.data->_continuation_solver_status = 0; // ?
			cov.data->_continuation_solver_time = -1.0;
			cov.data->_continuation_solver_iterations = 0;
			
			// TODO: Not safe but it should be initialise according to whether we loaded a CovManifold
			// or a CovParManifold
			for(int i = 0; i < cov.size(); ++i)
				cov.data->_continuation_solver_domain_type.push_back(BOX);
		}
		else
		{
			format_id.pop();
			format_version.pop();
			
			// TODO: check ?
			unsigned int status = read_pos_int(*f);
			cov.data->_continuation_solver_status = status;
			
			cov.data->_continuation_solver_time = read_double(*f);
			
			cov.data->_continuation_solver_iterations = read_pos_int(*f);
			
			// components
			unsigned int nb_compo = read_pos_int(*f);
			for(unsigned int i = 0; i < nb_compo; ++i)
			{
				size_t start = read_pos_int(*f);
				size_t end = read_pos_int(*f);
				cov.data->_continuation_solver_components.push_back(std::make_pair(start, end));
			}
			
			for(unsigned int i = 0; i < cov.size(); ++i)
			{
				unsigned int type = read_pos_int(*f);
				switch(type)
				{
					case 0 : cov.data->_continuation_solver_domain_type.push_back(DomainType::PARALLELOTOPE); break;
					case 1 : cov.data->_continuation_solver_domain_type.push_back(DomainType::BOX); break;
					default : ibex_error("[CovContinuation]: Unkown domain type when reading.");
				}
			}
		}
		
		return f;
	}
				
	std::ofstream* CovContinuation::write(const char* filename, const CovContinuation& cov, std::stack<unsigned int>& format_id, std::stack<unsigned int>& format_version)
	{
		format_id.push(subformat_number);
		format_version.push(FORMAT_VERSION);

		std::ofstream* f = CovParManifold::write(filename, cov, format_id, format_version);
		
		write_pos_int(*f, cov.data->_continuation_solver_status);
		
		write_double(*f, cov.data->_continuation_solver_time);
		
		write_pos_int(*f, cov.data->_continuation_solver_iterations);
		
		write_pos_int(*f, cov.nb_components());
		
		for(unsigned int i = 0; i < cov.nb_components(); ++i)
		{
			auto compo = cov.data->_continuation_solver_components[i];
			write_pos_int(*f, compo.first);
			write_pos_int(*f, compo.second);
		}
		
		for(unsigned int i = 0; i < cov.size(); ++i)
		{
			switch(cov.data->_continuation_solver_domain_type[i])
			{
				case PARALLELOTOPE : write_pos_int(*f, 0); break;
				case BOX : write_pos_int(*f, 1); break;
				default : ibex_error("[CovContinuation]: Unkown domain type when writing.");
			}

		}
		
		return f;
	}

	void CovContinuation::format(std::stringstream& ss, const std::string& title, std::stack<unsigned int>& format_id, std::stack<unsigned int>& format_version)
	{
		format_id.push(subformat_number);
		format_version.push(FORMAT_VERSION);

		CovManifold::format(ss, title, format_id, format_version);

		ss << space << " - 1 integer:     the continuation solver status\n"
		<< space << "                  - 0 = Transition matrices are stored alongside\n"
		<< space << "                  - 1 = Transition matrices are not stored.\n"
		<< space << " - 1 double:      time (in seconds)\n"
		<< "|  CovContinuation  |" 
		<< space << " - 1 integer:     the number of iterations (including failures)\n"
		<< space << " - 1 integer:     the number Nc of connected components\n"
		<< space << " - Nc int pairs:  the components seen as pairs of starting and ending index\n"
		<< space << " - N integer:     type of the stored element\n"
		<< separator;
	}

	std::string CovContinuation::format()
	{
		std::stringstream ss;
		std::stack<unsigned int> format_id;
		std::stack<unsigned int> format_version;
		format(ss, "CovContinuation", format_id, format_version);
		return ss.str();
	}
	
	std::ostream& operator<<(std::ostream& os, const CovContinuation& cov)
	{
		
		os << " status " << cov.solver_status() << std::endl;
		os << " timer: " << cov.get_time() << " s" << std::endl;
		os << " # iterations: " << cov.get_iterations() << std::endl;
		os << " # components: " << cov.nb_components() << std::endl;
		for(int i = 0; i < cov.nb_components(); ++i)
		{
			auto c = cov.get_component(i);
			os << "\t component " << (i+1) << " from " << (c.first+1) << " to " << c.second << std::endl;
		}

		for(int i = 0; i < cov.size(); ++i)
		{
			
			if(cov.is_box(i))
				os << " element n° " << (i+1) << " is a box " << std::endl;
		}
		
		return os;
	}
} // namespace ibex
