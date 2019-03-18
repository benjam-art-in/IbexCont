/* ============================================================================
 * I B E X 
 * ============================================================================
 * Copyright   : Ecole Polytechnique (FRANCE)
 * License     : LGPL, see the file COPYING.LESSER.
 * Author(s)   : Benjamin Martin
 * Created     : 15.03.2019
 * Last Update : 
 * ---------------------------------------------------------------------------- */

#ifndef __COV_CONTINUATION_H__
#define __COV_CONTINUATION_H__

#include "ibex_CovParManifold.h"

namespace ibex
{
	
	class CovContinuation : public CovParManifold
	{
		public:
		
			/**
			 * \brief The type of domains added
			 * 
			 * TODO: when BOX, adds indication on sense of continuation ?
			 **/
			typedef enum {PARALLELOTOPE, BOX} DomainType;
		
			CovContinuation(size_t n, size_t m, size_t nb_ineq = 0);
			
			CovContinuation(const char* filename);
			
			CovContinuation(const Cov& cov, bool copy = false);
			
			~CovContinuation();
			
			/**
			 * \brief Adds a parallelotope to the CovParManifold.
			 **/
			virtual void add_solution_parallelotope(const Parallelotope& p);

			virtual void add_solution(const IntervalVector& box, const VarSet& varset);
			
			bool is_parallelotope(int j) const;
			
			bool is_box(int j) const;
			
			void set_time(double t);
			
			double get_time() const;
			
			void set_iterations(unsigned int it);
			
			unsigned int get_iterations() const;
			
			unsigned int nb_components() const;
			
			std::pair<size_t, size_t> get_component(int j) const;
			
			void set_solver_status(unsigned int s);
			
			unsigned int solver_status() const;

			void begin_component();
			
			void end_component();
			
			/**
			 * \brief Save this as a COV file.
			 */
			void save(const char* filename) const;

			/**
			 * \brief Display the format of a CovParManifold file.
			 */
			static std::string format();

			/**
			 * \brief CovParManifold file format version.
			 */
			static const unsigned int FORMAT_VERSION;
		
		protected:
		
			/**
			 * \brief Load a manifold covering from a COV file.
			 */
			static std::ifstream* read(const char* filename, CovContinuation& cov, std::stack<unsigned int>& format_id, std::stack<unsigned int>& format_version);

			/**
			 * \brief Write a manifold covering into a COV file.
			 */
			static std::ofstream* write(const char* filename, const CovContinuation& cov, std::stack<unsigned int>& format_id, std::stack<unsigned int>& format_version);

			static void format(std::stringstream& ss, const std::string& title, std::stack<unsigned int>& format_id, std::stack<unsigned int>& format_version);

			
			/**
			 * \brief Subformat level.
			 */
			static const unsigned int subformat_level;

			/**
			 * \brief Subformat identifying number.
			 */
			static const unsigned int subformat_number;
			
			size_t start_component;

			struct Data {
				unsigned int _continuation_solver_status;
				double		 _continuation_solver_time;
				unsigned int _continuation_solver_iterations;
				std::vector<std::pair<size_t,size_t>> _continuation_solver_components;
				std::vector<DomainType> _continuation_solver_domain_type;
			} *data;

			bool own_data;
	}; // CovContinuation
	
	std::ostream& operator<<(std::ostream& os, const CovContinuation& cov);
	
	inline bool CovContinuation::is_parallelotope(int j) const
	{
		return data->_continuation_solver_domain_type[j] == PARALLELOTOPE;
	}
	
	inline bool CovContinuation::is_box(int j) const
	{
		return data->_continuation_solver_domain_type[j] == BOX;
	}
	
	inline void CovContinuation::set_time(double t)
	{
		data->_continuation_solver_time = t;
	}
	
	inline double CovContinuation::get_time() const
	{
		return data->_continuation_solver_time;
	}
	
	inline void CovContinuation::set_iterations(unsigned int it)
	{
		data->_continuation_solver_iterations = it;
	}
	
	inline unsigned int CovContinuation::get_iterations() const
	{
		return data->_continuation_solver_iterations;
	}
	
	inline unsigned int CovContinuation::nb_components() const
	{
		return data->_continuation_solver_components.size();
	}
	
	inline std::pair<size_t, size_t> CovContinuation::get_component(int j) const
	{
		return data->_continuation_solver_components[j];
	}

	inline void CovContinuation::set_solver_status(unsigned int s)
	{
		data->_continuation_solver_status = s;
	}

	inline unsigned int CovContinuation::solver_status() const
	{
		return data->_continuation_solver_status;
	}
} // namespace ibex

#endif // __COV_CONTINUATION_H__
