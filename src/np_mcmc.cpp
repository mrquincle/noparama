#include <assert.h>

#include <Eigen/Dense>

#include <np_mcmc.h>
#include <np_update_clusters.h>
#include <np_update_cluster_population.h>

#include <iostream>

#include <pretty_print.hpp>

#include <chrono>

using namespace std;

MCMC::MCMC(
		UpdateClusters & update_clusters, 
		UpdateClusterPopulation & update_cluster_population,
		std::default_random_engine & generator,
		distribution_t & nonparametrics
	):
	_update_clusters(update_clusters), 
	_update_cluster_population(update_cluster_population),
	_generator(generator),
	_nonparametrics(nonparametrics)
{
	_verbosity = 4;
}

void MCMC::run(dataset_t & dataset) {
	int T = 20;
	int K = 30;

	int N = dataset.size();

	int number_mh_steps = 20;

	fout << "Add data to membership matrix" << endl;
	for (int i = 0; i < N; ++i) {
		data_t & datum = *dataset[i];
		data_id_t j = _membertrix.addData(datum);
//		fout << "Data " << datum << " got index " << j <<endl;
		assert (i == j);
	}

	// initialize clusters
	fout << "Add clusters to membership matrix" << endl;
	for (int k = 0; k < K; ++k) {
		// sample sufficient statistics from nonparametrics
		//Suffies & suffies = _nonparametrics.sample(_generator);
		Suffies *suffies = _nonparametrics(_generator);

		cluster_t *cluster = new cluster_t(*suffies);
		int kTest = _membertrix.addCluster(cluster);

		assert (kTest == k);
	}

	std::vector<double> weights(K);
	for (int k = 0; k < K; ++k) {
		weights[k] = 1/(double)K;
	}
	std::discrete_distribution<int> distribution(weights.begin(), weights.end());

	// randomly assign data to clusters
	for (int i = 0; i < N; ++i) {
		int k = distribution(_generator);
	//	fout << "Assign data item " << i << " to cluster " << k <<endl;
		np_error_t err = _membertrix.assign(k, i);
		if (err) fout << "Error: " << np_error_str[err] << endl;
	}	

	// clean up clusters that have been created, but didn't get data assigned
	int removed = _membertrix.cleanup();
	fout << "Removed " << removed << " empty clusters" << endl;
	removed = _membertrix.cleanup();
	fout << "Removed " << removed << " empty clusters" << endl;

	// update clusters
	for (int t = 0; t < T; ++t) {
			
		foutvar(5) << "Metropolis Hastings update " << t << endl;
		
		// iterator over all observations (observations get shifted around, but should not be counted twice)
		for (int i = 0; i < N; ++i) {
			fout << "Retract observation " << i << " from membertrix" << endl;
			
			// remove observation under consideration from cluster
			_membertrix.retract(i);
			
			// update cluster assignments, delete and create clusters
			_update_cluster_population.update(_membertrix, i);
		}
			
		fout << "Update cluster " << endl;

		// update cluster parameters
		_update_clusters.update(_membertrix, number_mh_steps);
	}

	clusters_t &clusters = _membertrix.getClusters();
	for (auto cluster_pair: clusters) {
		auto const &key = cluster_pair.first;
		foutvar(5) << "Cluster ";
		_membertrix.print(key, std::cout);
		std::cout << endl;
	}
}
