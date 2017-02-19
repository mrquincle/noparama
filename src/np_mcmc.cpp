#include <assert.h>
#include <sstream>
#include <fstream>

#include <Eigen/Dense>

#include <np_mcmc.h>
#include <np_init_clusters.h>
#include <np_update_clusters.h>
#include <np_update_cluster_population.h>

#include <iostream>

#include <pretty_print.hpp>

#include <chrono>

using namespace std;

MCMC::MCMC(
		default_random_engine & generator,
		InitClusters & init_clusters, 
		UpdateClusters & update_clusters, 
		UpdateClusterPopulation & update_cluster_population
	):
	_generator(generator),
	_init_clusters(init_clusters), 
	_update_clusters(update_clusters), 
	_update_cluster_population(update_cluster_population)
{
	_verbosity = 4;
}

void MCMC::run(dataset_t & dataset, int T) {
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

	_init_clusters.init(_membertrix, K);

	vector<double> weights(K);
	for (int k = 0; k < K; ++k) {
		weights[k] = 1/(double)K;
	}
	discrete_distribution<int> distribution(weights.begin(), weights.end());

	// randomly assign data to clusters
	for (int i = 0; i < N; ++i) {
		int k = distribution(_generator);
	//	fout << "Assign data item " << i << " to cluster " << k <<endl;
		np_error_t err = _membertrix.assign(k, i);
		if (err) fout << "Error: " << np_error_str[err] << endl;
	}	
	foutvar(7) << "We have in total " << _membertrix.count() << " assigned data items" << endl;

	// clean up clusters that have been created, but didn't get data assigned
	int removed = _membertrix.cleanup();
	fout << "Removed " << removed << " empty clusters" << endl;
	removed = _membertrix.cleanup();
	fout << "Removed " << removed << " empty clusters" << endl;
	foutvar(7) << "We have in total " << _membertrix.count() << " assigned data items" << endl;

	// update clusters
	for (int t = 0; t < T; ++t) {
		
		if (t % 100 == 0) { foutvar(5) << "Metropolis Hastings update " << t << endl; }

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

}

const membertrix & MCMC::getMembershipMatrix() const {
	return _membertrix;
}
