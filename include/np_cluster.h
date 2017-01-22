
/**
 * A Cluster exists out of two things:
 *
 *  - sufficient statistics / parameters
 *  - data items
 *
 * There are no dependencies between clusters assumed (hence no references to other clusters or hierarchy).
 */

typedef MatrixXd t_datum;

typedef std::vector<t_datum> t_data;

class Cluster {
	private: 
		// The parameters
		SufficientStatistics _statistics;

		// The data is a set of individual data items
		t_data _data;

	public:
		Cluster(SufficientStatistics & statistics);
		
		~Cluster();

		// Get sufficient statistics
		SufficientStatistics & getSufficientStatistics();

		// Add a data item
		size_t add(t_datum & item);
		
		// Remove a data item
		void erase(size_t position);

		// Number of data items
		int count();
};
