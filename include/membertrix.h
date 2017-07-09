#pragma once

#include <Eigen/Dense>
#include "np_data.h"
#include "np_cluster.h"

#include <map>
#include <unordered_map>

//! Hashmap for the cluster indices and cluster objects
typedef std::unordered_map<cluster_id_t, cluster_t*> clusters_t;

//! A dataset per cluster
typedef std::unordered_map<cluster_id_t, dataset_t*> clusters_dataset_t;

enum np_error_t { error_none, error_already_assigned, error_assignment_remaining, error_assignment_absent };

static std::map< np_error_t, const char * > np_error_str = {
	{error_none,                       "none"},
	{error_already_assigned,           "already assigned"},
	{error_assignment_remaining,       "assignment remaining"},
	{error_assignment_absent,          "assignment absent"}
};

/*!
 * The binary matrix is currently defined as a dense matrix. This needs actually some profiling to know if a sparse 
 * matrix is faster. Namely, the content of the matrix is sparse with one non-false value per row. However, for now
 * a dense matrix is used, because it is assumed that how rows and columns are accessed favors a dense matrix.
 */
typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> binary_matrix_t;

/*!
 * The membertrix data structure is a binary matrix optimimized for storing membership information.
 * The membership is asymmetric. A data item can be assigned to only one cluster. In contrary, a cluster can have
 * multiple data points as members.
 *
 * The data points and clusters are stored in separate vectors.
 *
 * The structure is stored with clusters as columns, and data items as rows. Reasons:
 *   + Eigen by default uses column-major storage [1]
 *   + The getData() function retrieves multiple data items (columns)
 *   + The getClusterId() function can return on first item found (rows)
 *
 * Usage:
 *   + Call cluster_id = addCluster(cluster),   store cluster id somewhere safe
 *   + Call data_id = addData(data),            store data id somewhere safe
 *   + Call assign(cluster_id, data_id),        assign data to clusters
 *
 * References:
 *   [1] https://eigen.tuxfamily.org/dox/group__TopicStorageOrders.html
 */
class membertrix {
	private:
		// membership matrix
		binary_matrix_t _membership_matrix;

		// map from cluster pointers to column entries
		clusters_t _cluster_objects; 

		// map from data pointers to row entries
		dataset_t _data_objects; 

		// store data items per cluster
		clusters_dataset_t _clusters_dataset;

		// verbosity
		char _verbosity;

	protected:
		
		bool exists(cluster_id_t cluster_id);
		
	public:
		/*!
		 * The default constructor.
		 */
		membertrix();

		/*!
		 * The copy constructor. This is not a true copy constructor. If you copy a membership matrix we assume you
		 * want to optimize the internal structures. This invalidates all cluster_id's. If an exact clone is required
		 * you will need to add a clone() member function.
		 *
		 * This constructor calls addData and addCluster to have all internal datastructures consistent and reduce
		 * the matrix to the minimum size. The alternative would be all kind of book-keeping swapping columns in the
		 * matrix, moving datasets from own cluster to the next, etc.
		 */
		membertrix(const membertrix &other);

		/*!
		 * The destructor.
		 */
		~membertrix();

		/*!
		 * Add a cluster to the membership matrix. The cluster is not physically stored, only a reference is kept. If
		 * the memory is deallocated, errors can be expected.
		 *
		 * The returned index should be kept as a reference for use in the functions assign() and retract().
		 *
		 * @param[in] cluster          A cluster object
		 * @return                     An index to the given cluster object
		 */
		cluster_id_t addCluster(cluster_t *cluster);

		/*!
		 * Get cluster given cluster id.
		 *
		 * @param[in] cluster_id       An index to a cluster object
		 * @return                     The cluster object itself
		 */
		cluster_t * getCluster(cluster_id_t cluster_id);

		/*!
		 * Add a data point to the membership matrix. The data are not physically stored, only a reference is kept. If
		 * the memory is deallocated, errors can be expected.
		 *
		 * The returned index should be kept as a reference for use in the functions assign() and retract().
		 *
		 * @param[in] data             A data object
		 * @return                     An index to the given data object
		 */
		data_id_t addData(data_t & data);

		/*!
		 * Return data point with given index.
		 *
		 * @param[in] data_id          An index to a particular data point
		 * @return                     A data point that has been set previously through addData
		 */
		data_t * getDatum(data_id_t data_id);

		/*!
		 * Assign a previously added data item (through addData) to a previously added cluster (through addCluster).
		 *
		 * @param[in] cluster_id       An index to a cluster object
		 * @param[in] data_id          An index to a data point
		 * @return                     Error if data item does not exist (for example)
		 */
		np_error_t assign(cluster_id_t cluster_id, data_id_t data_id);
		
		/*!
		 * Retract a previously assigned data-cluster pair (through assign). If the cluster does not have any data
		 * points left, also the object will be deallocated depending on auto_remove setting.
		 *
		 * @param[in] cluster_id       An index to a cluster object
		 * @param[in] data_id          An index to a data point
		 * @param[in] auto_remove      Automatically deallocate cluster object if there is no data assigned anymore
		 */
		np_error_t retract(cluster_id_t cluster_id, data_id_t data_id, bool auto_remove = true);

		/*!
		 * Retract a previously assigned data-cluster pair (through assign) where the search for this particular 
		 * cluster is left to getClusterId(data_id). If the cluster does not have any data points left, also the object 
		 * will be deallocated (depending on auto_remove flag). This has the same effect as:
		 * 
		 *   retract(getClusterId(data_id), data_id);
		 *
		 * @param[in] data_id          An index to a data point
		 * @param[in] auto_remove      Automatically deallocate cluster object if there is no data assigned anymore
		 */
		np_error_t retract(data_id_t data_id, bool auto_remove = true);

		/*!
		 * Remove a cluster (should be empty). This function should be called if auto_remove is set to false in the
		 * retract() functions. If auto_remove is set to true (default) there is no need to call remove().
		 *
		 * @param[in] cluster_id       Index to cluster to be removed
		 * @return                     Error if cluster is still non-empty for example
		 */
		np_error_t remove(cluster_id_t cluster_id);

		/*!
		 * If the data item is assigned to any cluster this function will return true. In all other cases it returns
		 * false.
		 *
		 * @param[in] data_id          An index to a data point
		 * @return                     Boolean representing any assignment
		 */
		bool assigned(data_id_t data_id) const;

		/*!
		 * Get all assignments to given cluster.
		 *
		 * @param[in] cluster_id       An index to a cluster object
		 * @param[out] data_ids        Set of data ids
		 */
		void getAssignments(cluster_id_t cluster_id, data_ids_t & data_ids) const;

		/*!
		 * Get the cluster id given a particular data id.
		 * @param[in] data_id          An index to a data point
		 * @return                     An index to a cluster
		 */
		cluster_id_t getClusterId(data_id_t data_id) const;

		/*!
		 * Get all clusters to iterate over them. The cluster set is const to protect the user from accidentally 
		 * removing clusters in a for-loop in a way that destroys the user iterator.
		 * 
		 * This function can be used to adjust the parameters of the cluster_t objects. The function only reads
		 * cluster information and does not change the membertrix instance, hence it is const.
		 *
		 * @return                     Set of clusters
		 */
		const clusters_t & getClusters() const;

		/*!
		 * Get number of clusters. Note that retract and assign adjust the number of clusters!
		 *
		 * @return                     Total number of clusters
		 */
		size_t getClusterCount() const;

		/*!
		 * Aggressive restructuring of all data structures. This will relabel all cluster_id's to consecutive numbers.
		 * The assignments are still valid but with different cluster_id's. 
		 *
		 * This is called automatically on assignments!!
		 */
		void relabel();

		/*!
		 * Print cluster to stream.
		 */
		void print(cluster_id_t cluster_id, std::ostream &os) const;
		
		/*!
		 * Print entire membership matrix to stream.
		 */
		void print(std::ostream& os) const;

		/*!
		 * Allow a membership to be printed to a standard stream using the << operator.
		 */
		friend std::ostream& operator<<(std::ostream& os, const membertrix& m);  

		/*!
		 * Return all data points.
		 *
		 * @return                     The entire dataset (do not need to be assigned)
		 */
		dataset_t* getData();
		
		/*!
		 * Return all data points that are assigned to a particular cluster.
		 *
		 * @param[in] cluster_id       An index to a particular cluster
		 * @param[out] dataset         A dataset (vector) of data points that have been assigned through assign()
		 */
		dataset_t* getData(const cluster_id_t cluster_id) const;

		void getData(const cluster_id_t cluster_id, dataset_t & dataset) const;

		/*!
		 * Return a particular subset of data points. Can belong to a particular cluster or not. As long as they have
		 * been added through addData(). They do not have to be assigned to a cluster yet.
		 *
		 * @param[in] data_ids         A set of data point ids
		 * @return                     The data points themselves
		 */
		void getData(const data_ids_t data_ids, dataset_t & dataset) const;

		/*!
		 * Return count of data points within the given cluster.
		 * @param[in] cluster_id       An index to a particular cluster
		 * @return                     Number of data points (should be the same as getData(cluster_id).size()).
		 */
		size_t count(cluster_id_t cluster_id) const;

		/*!
		 * Return total number of data points. This should be the same as calling count(cluster_id_t) for each
		 * cluster returned by getClusters().
		 * @return                     The total number of data points
		 */
		size_t count() const;
	
		/*!
		 * Indicate if a cluster is empty or non-empty (one or more data items assigned to it).
		 * @return                     True or false depending on zero or nonzero data items in the cluster
		 */
		bool empty(cluster_id_t cluster_id);

		/*!
		 * Clean up internal data structures after data items and clusters have been assigned to remove all clusters
		 * that didn't get assigned.
		 */
		int cleanup();

		/*!
		 * The assignment operator is implemented by not passing by reference, but having the argument as a copy.
		 * Subsequently, only a swap operation needs to be called.
		 *
		 * @param[in] membertrix       Another membertrix object
		 * @return                     A copy of this membertrix object, optimized.
		 */
		membertrix &operator=(membertrix other);

		/*!
		 * Swap all member fields of the two objects. This is a very lightweight implementation only swapping the
		 * five member fields on the level of references, nothing is copied.
		 */
		friend void swap(membertrix& first, membertrix& second) {
			using std::swap;

			swap(first._membership_matrix, second._membership_matrix);
			swap(first._cluster_objects, second._cluster_objects);
			swap(first._data_objects, second._data_objects);
			swap(first._clusters_dataset, second._clusters_dataset);
			swap(first._verbosity, second._verbosity);
		}
};
