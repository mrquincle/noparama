#pragma once

#include <vector>
#include <random>

#include <membertrix.h>

#include <np_data.h>

#include <np_init_clusters.h>
#include <np_update_clusters.h>
#include <np_update_cluster_population.h>

/*!
 * The MCMC class can be equiped with InitCluster, UpdateCluster, and UpdateClusterPopulation objects that each will
 * be used in a Monte-Carlo Monte Chain algorithm. This algorithm is run for T steps and returns a membership matrix
 * in the form of a membertrix object.
 */
class MCMC {
	private:
		//! Todo: set generator
		std::default_random_engine & _generator;
	
		//! Membertrix object itself
		membertrix _membertrix;
		
		//! Reference to InitCluster object
		InitClusters & _init_clusters;

		//! Reference to UpdateCluster object
		UpdateClusters & _update_clusters;

		//! Reference to UpdateClusterPopulation object
		UpdateClusterPopulation & _update_cluster_population;

		//! Verbosity level for the class
		char _verbosity;

		//! Fixed scan or random scan
		bool _random_scan;

		//! Number of simultaneous data points under consideration
		int _subset_count;
	public:
		/*!
		 * Constructor for MCMC.
		 * 
		 * @param[in] generator                  A random number generator
		 * @param[in] init_clusters              An object that initializes clusters
		 * @param[in] update_clusters            An object that updates the cluster parameters
		 * @param[in] update_cluster_population  An object that generates and removes clusters
		 */
		MCMC(
				std::default_random_engine & generator,
				InitClusters & init_clusters, 
				UpdateClusters & update_clusters, 
				UpdateClusterPopulation & update_cluster_population,
				int subset_count
			);

		/*!
		 * Run the MCMC method (can be one or multiple chains) for T steps.
		 *
		 * @param[in] dataset                    The dataset to run the MCMC on
		 * @param[in] T                          The number of steps 
		 */
		void run(dataset_t & dataset, int T);

		/*!
		 * The results of the MCMC algorithm is a membership matrix.
		 *
		 * @return                               The membership matrix (const)
		 */
		const membertrix & getMembershipMatrix() const;
};
