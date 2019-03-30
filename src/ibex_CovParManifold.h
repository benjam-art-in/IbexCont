/* ============================================================================
 * I B E X 
 * ============================================================================
 * Copyright   : Ecole Polytechnique (FRANCE)
 * License     : LGPL, see the file COPYING.LESSER.
 * Author(s)   : Benjamin Martin
 * Created     : 14.03.2019
 * Last Update : 
 * ---------------------------------------------------------------------------- */

#ifndef __IBEX_COV_PAR_MANIFOLD_H__
#define __IBEX_COV_PAR_MANIFOLD_H__

#include "ibex_Parallelotope.h"
#include "ibex_CovManifold.h"

namespace ibex
{
	
	/**
	 * \brief Cov structure for the covering of a manifold by parallelotopes.
	 * 
	 *  TODO: we can store inside the data the parallelotopes directly. Only write the 
	 * (matrix, aux, xtilde) format to the file. This requires rebuilding the parallelotopes
	 * when loading a file.
	 **/
	class CovParManifold : public CovManifold
	{
		public:
		
			/**
			 * \brief Mode for storing the parallelotope
			 **/
			typedef enum {FULL_MATRIX, NO_MATRIX} ParMode;
			
			
			/**
			 * \brief Create an empty manifold covering.
			 */
			CovParManifold(size_t n, size_t m, size_t nb_ineq=0, ParMode par_mode=FULL_MATRIX);

			/**
			 * \brief Load a manifold covering from a COV file.
			 */
			CovParManifold(const char* filename);

			/**
			 * \brief Conversion from a COV.
			 */
			CovParManifold(const Cov& cov, bool copy=false);

			/**
			 * \brief Delete this
			 */
			~CovParManifold();
			
			/**
			 * \brief Adds a parallelotope to the CovParManifold.
			 **/
			virtual void add_solution_parallelotope(const Parallelotope& p);

			/**
			 * \brief Adds a solution box to the CovParManifold.
			 * \note adds nothing to this level, but pass to the manifold level.
			 **/
			virtual void add_solution(const IntervalVector& box, const VarSet& varset);

			/**
			 * \brief Save this as a COV file.
			 */
			void save(const char* filename) const;
		
			/**
			 * \brief Type of parallelotopes.
			 *
			 * By default: FULL_MATRIX.
			 */
			ParMode par_mode() const;

			/**
			 * \brief Number of parallelotopes.
			 */
			size_t nb_parallelotopes() const;
			
			size_t id_parallelotope(int j) const;
			
			/**
			 * \brief Get the j'th parallelotope.
			 * 
			 * \note requires to reconstruct the Cinv IntervalMatrix of the parallelotope.
			 * This operation has O(n^3) complexity.
			 */
			Parallelotope get_parallelotope(int j) const;

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
			 * \brief Load a parallelotopic manifold covering from a COV file.
			 */
			static std::ifstream* read(const char* filename, CovParManifold& cov, std::stack<unsigned int>& format_id, std::stack<unsigned int>& format_version);

			/**
			 * \brief Write a parallelotopic manifold covering into a COV file.
			 */
			static std::ofstream* write(const char* filename, const CovParManifold& cov, std::stack<unsigned int>& format_id, std::stack<unsigned int>& format_version);

			static void format(std::stringstream& ss, const std::string& title, std::stack<unsigned int>& format_id, std::stack<unsigned int>& format_version);

			static Vector read_vector(std::ifstream& f, size_t n);
			
			static void write_vector(std::ofstream& f, const Vector& v);
			
			static Matrix read_sq_matrix(std::ifstream& f, size_t n);
			
			/**
			 * \brief write a square matrix to the file (row wise)
			 */
			static void write_sq_matrix(std::ofstream& f, const Matrix& m);
			
			
			/**
			 * \brief Subformat level.
			 */
			static const unsigned int subformat_level;

			/**
			 * \brief Subformat identifying number.
			 */
			static const unsigned int subformat_number;

			struct Data {
				ParMode                 	 _par_manifold_par_mode;  // storage mode
				std::vector<size_t>          _par_manifold_parallelotopes;  // indices of parallelotopes
				std::vector<IntervalVector>  _par_manifold_aux_boxes;   // all the auxiliary boxes
				std::vector<Vector>			 _par_manifold_centers;	// all the centers
				std::vector<Matrix>			 _par_manifold_transition_matrices;	// all the transition matrices
			} *data;

			bool own_data;
		
	}; // CovParManifold
	
	std::ostream& operator<<(std::ostream& os, const CovParManifold& cov);
	
	
	inline size_t CovParManifold::nb_parallelotopes() const
	{
		return data->_par_manifold_parallelotopes.size();
	}
	
	inline size_t CovParManifold::id_parallelotope(int j) const
	{
		return data->_par_manifold_parallelotopes[j];
	}
	
	inline Parallelotope CovParManifold::get_parallelotope(int j) const
	{
		assert(j >= 0 && j <= nb_parallelotopes());
		return Parallelotope( data->_par_manifold_transition_matrices[j],
							  data->_par_manifold_aux_boxes[j],
							  data->_par_manifold_centers[j]);
					
	}
	
	inline CovParManifold::ParMode CovParManifold::par_mode() const
	{
		return data->_par_manifold_par_mode;
	}

} // namespace ibex

#endif // __IBEX_COV_PAR_MANIFOLD_H__
