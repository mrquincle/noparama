
#include <Eigen/dense>

enum error_t { error_none, error_already_assigned, error_assignment_remaining, error_assignment_absent };
	
/*!
 * The binary matrix is currently defined as a dense matrix. This needs actually some profiling to know if a sparse 
 * matrix is faster. Namely, the content of the matrix is sparse with one non-false value per row. However, for now
 * a dense matrix is used, because it is assumed that how rows and columns are accessed favors a dense matrix.
 */
typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> binary_matrix_t;

typedef std::vector<cluster_t> clusters_t;

//! A dataset
typedef std::set<data_t> dataset_t;

//! A dataset per cluster
typedef std::vector< dataset_t > clusters_dataset_t;

/*!
 * The membertrix data structure is a binary matrix optimimized for storing membership information.
 * The membership is asymmetric. A data item can be assigned to only one cluster. In contrary, a cluster can have
 * multiple data points as members.
 *
 * The structure is stored with clusters as columns, and data items as rows. Reasons:
 *   + Eigen by default uses column-major storage [1]
 *   + The getData() function retrieves multiple data items
 *   + The getCluster() function can return on first item found
 *
 * References:
 *   [1] https://eigen.tuxfamily.org/dox/group__TopicStorageOrders.html
 */
class membertrix {
	private:
		// membership matrix
		binary_matrix_t _membership_matrix;

		// map from cluster pointers to column entries
		clusters_t _column_labels; 

		// map from data pointers to row entries
		std::vector<data_t> _row_labels; 

		// store data items per cluster
		clusters_dataset_t _clusters_dataset;
	public:

		membertrix();

		id_cluster_t addCluster(cluster_t & cluster);

		id_data_t addData(data_t & data);

		error_t assign(id_cluster_t & cluster, id_data_t & data);
		
		error_t retract(id_cluster_t & id_cluster, id_data_t & id_data);

		error_t retract(id_data_t & id_data);

		id_cluster_t getCluster(id_data_t & id_data);

		clusters_t & getClusters();

		dataset_t & getData(id_cluster_t & id_cluster);
};

