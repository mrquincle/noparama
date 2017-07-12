#pragma once

#include <vector>

#include "np_suffies.h"

#include <iostream>

// Forward reference to cluster_t
class cluster_t;

//! Index to a cluster
typedef int cluster_id_t;

/*!
 * A set of cluster ids, not the clusters themselves.
 */
typedef std::vector<cluster_id_t> cluster_ids_t;

/*!
 * A cluster is represented by a cluster_t structure and contains parameters in the form of sufficient statistics.
 *
 * There are no dependencies between clusters assumed (hence no references to other clusters or hierarchy).
 *
 * The data is not stored in the cluster itself, but in the membertrix object.
 */
class cluster_t {
	private: 
		//! The parameters are sufficient statistics of the related probability density function
		//! They are stored by reference, or else the construction will be a copy by value operation.
		Suffies & _suffies;
	public:
		cluster_t(Suffies & suffies): _suffies(suffies) {
		}

		cluster_t(const cluster_t &other): _suffies(other._suffies) {
			std::cout << "Copy cluster " << std::endl;
		}

		~cluster_t() {
		};

		//! Get sufficient statistics
		Suffies & getSuffies() const {
			return _suffies;
		}	
		
		//! Set sufficient statistics
		void setSuffies(Suffies & suffies) {
			_suffies = suffies;
		}

		void updateSuffies(Suffies suffies) {
			// overwrite 
			_suffies = suffies;
		}
		
		friend void swap(cluster_t& first, cluster_t& second) {
			using std::swap;

			swap(first._suffies, second._suffies);
		}
		
		cluster_t &operator=(cluster_t other) {
			swap(*this, other);

			return *this;
		}
};
