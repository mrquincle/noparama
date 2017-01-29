
/*!
 * A cluster is represented by a cluster_t structure and contains parameters in the form of sufficient statistics.
 *
 * There are no dependencies between clusters assumed (hence no references to other clusters or hierarchy).
 */

#include <np_data.h>

class cluster_t {
	private: 
		// The parameters
		SufficientStatistics _statistics;

		// The data is a set of individual data items
		//t_data _data;

	public:
		cluster_t(SufficientStatistics & statistics);
		
		~cluster_t();

		// Get sufficient statistics
		SufficientStatistics & getSufficientStatistics();
		
		// Set sufficient statistics
		void setSufficientStatistics(SufficientStatistics & statistics);

		// Add a data item
		//size_t add(t_datum & item);
		
		// Remove a data item
		//void erase(size_t position);

		// Number of data items
		//int count();

		// Get data
		//t_data & data();
};
