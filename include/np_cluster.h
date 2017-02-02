#pragma once

#include <vector>

#include <np_data.h>
#include <np_suffies.h>

// Forward references
class cluster_t;

//! Index to a cluster
typedef int cluster_id_t;

typedef std::vector<cluster_t> clusters_t;

//! A dataset per cluster
typedef std::vector< dataset_t > clusters_dataset_t;

/*!
 * A cluster is represented by a cluster_t structure and contains parameters in the form of sufficient statistics.
 *
 * There are no dependencies between clusters assumed (hence no references to other clusters or hierarchy).
 */
class cluster_t {
	private: 
		//! The parameters are sufficient statistics of the related probability density function
		Suffies _statistics;

		// The data is a set of individual data items
		//t_data _data;

	public:
		cluster_t(Suffies & statistics);
		
		~cluster_t();

		//! Get sufficient statistics
		Suffies & getSuffies();
		
		//! Set sufficient statistics
		void setSuffies(Suffies & statistics);

		// Add a data item
		//size_t add(t_datum & item);
		
		// Remove a data item
		//void erase(size_t position);

		// Number of data items
		//int count();

		// Get data
		//t_data & data();
};
