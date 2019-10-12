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

#include <dim1algebra.hpp>

using namespace std;
using namespace algebra;

MCMC::MCMC(
		default_random_engine & generator,
		InitClusters & init_clusters, 
		UpdateClusters & update_clusters, 
		UpdateClusterPopulation & update_cluster_population,
		int subset_count,
		distribution_t & likelihood
	):
	_generator(generator),
	_init_clusters(init_clusters), 
	_update_clusters(update_clusters), 
	_update_cluster_population(update_cluster_population),
	_likelihood(likelihood)
{
	_verbosity = Information;

	// the number of sets we split / merge into, e.g. 3 for the triadic sampler
	_subset_count = subset_count;

	_max_likelihood = -std::numeric_limits<double>::infinity();
}

void MCMC::run(dataset_t & dataset, int T) {
	int K = 20;

	int N = dataset.size();

	int number_mh_steps = 40; // was 20
	number_mh_steps = 20;
//	number_mh_steps = 2;

	foutvar(Debug) << "Add data to membership matrix" << endl;
	for (int i = 0; i < N; ++i) {
		data_t & datum = *dataset[i];
		data_id_t j = _membertrix.addData(datum);
		foutvar(Debug) << "Data " << datum << " got index " << j <<endl;
		assert (i == j);
	}

	// initialize parameters of a first set of clusters 
	_init_clusters.init(_membertrix, K);

	// weights are all the same (does actually not need to be normalized: assigning 1 to each would be fine as well)
	vector<double> weights(K);
	for (int k = 0; k < K; ++k) {
		weights[k] = 1/(double)K;
	}
	discrete_distribution<int> distribution(weights.begin(), weights.end());

	// randomly assign data to clusters
	for (int i = 0; i < N; ++i) {
		// a cluster is selected uniformly (between 0 and K-1)
		int k = distribution(_generator);
#ifdef DEBUG
	    fout << "Assign data item " << i << " to cluster " << k <<endl;
#endif
		// assign to matrix 
		np_error_t err = _membertrix.assign(k, i);
		if (err) fout << "Error: " << np_error_str[err] << endl;
	}	
	foutvar(7) << "We have in total " << _membertrix.count() << " assigned data items" << endl;

	// clean up clusters that have been created, but didn't get data assigned
	int removed = _membertrix.cleanup();
	fout << "Removed " << removed << " empty clusters" << endl;
	removed = _membertrix.cleanup();
	fout << "Removed " << removed << " empty clusters" << endl;
	foutvar(7) << "After cleanup we have " << _membertrix.count() << " assigned data items (should be the same)" << endl;

	// just only pick a few, set M=1 for debugging, but M<N can also be used to perform fewer "large" steps
	int M = N;
//	M = 100;
//	M = 20;
	if (M > N) M = N;

	int cycle_print = 10;
	int cycle_max_likelihood = 5;

	foutvar(Notice) << "There will be " << T << " MH updates" << endl;
	foutvar(Notice) << "Within each update there will be " << M << " cluster updates" << endl;
	foutvar(Notice) << "Total number of steps is " << T*M << endl;

	// update clusters
	for (int t = 0; t < T; ++t) {
		
		if (t % cycle_print == 0) { 
			foutvar(Notice) << "Metropolis Hastings update " << t << endl; 
			_membertrix.relabel();
		}

		/* We iterate over all observations. We are using a fixed scan where we randomize all items once and then loop
		 * over them. For e.g. a triadic sampler we need three of such vectors.
		 */
		foutvar(Debug) << "Create multiple vectors with indices and randomly order each of them" << endl;
		std::vector< std::vector<int> > indices(_subset_count);
		for (int i = 0; i < _subset_count; ++i) {
			indices[i].resize(N);
			fill_successively(indices[i].begin(), indices[i].end());
			random_order(indices[i].begin(), indices[i].end());
		}

//#define TEST_ADDITIONAL_DUPLICATES
#ifdef TEST_ADDITIONAL_DUPLICATES
		int Q = 3; 
		if (Q > M) Q = M;
		for (int i = 0; i < Q; ++i) {
			double copy = indices[0][i];
			std::vector<int> subset(_subset_count);
			for (int j = 1; j < _subset_count; ++j) {
				indices[j][i] = copy;
			}
		}
#endif

		// We assume M to be very large. Hence sampling 3 items randomly and check on uniqueness by creating a set is
		// a fast way that makes sure there is a unique set sampled. On average we update the cluster assignments
		// much more often than that we update the cluster parameters. 
		// Note that the indices are randomly ordered. This means that M can be made smaller than N without the risk
		// that particular data points are never seen.
		foutvar(Debug) << "Create subset of size " << _subset_count << endl;
		for (int i = 0; i < M; ++i) {
			// create subset vector of size _subset_count
			std::vector<int> subset(_subset_count);
			std::set<int> uniq;
			for (int j = 0; j < _subset_count; ++j) {
				foutvar(Debug) << "Add " << indices[j][i] << " to subset" << endl;
				subset[j] = indices[j][i];
				uniq.insert(subset[j]);
			}
			if ((int)uniq.size() != _subset_count) {
				foutvar(Debug) << "Oops... Not all " << _subset_count << "indices are unique. Try Again." << endl;
				continue;
			}
	
			// update cluster assignments, delete and create clusters
			foutvar(Debug) << "Update cluster assignment i=" << i << std::endl;
			_update_cluster_population.update(_membertrix, subset);

		}
		//_update_cluster_population.printStatistics();

		fout << "Update cluster t=" << t << endl;

		// update cluster parameters
		_update_clusters.update(_membertrix, number_mh_steps);

		if (t % cycle_max_likelihood == 0) {
			considerMaxLikelihood();
		}
	}

}

const membertrix & MCMC::getMembershipMatrix() const {
	return _membertrix;
}

const membertrix & MCMC::getMaxLikelihoodMatrix() const {
	return _max_likelihood_membertrix;
}

void MCMC::considerMaxLikelihood() {
//	fout << "Check max likelihood " << endl;
	const clusters_t &clusters = _membertrix.getClusters();
	double current_likelihood = .0;
	for (auto cluster_pair: clusters) {
		auto const key = cluster_pair.first;
		auto const cluster = cluster_pair.second;
		dataset_t *dataset = _membertrix.getData(key);
		_likelihood.init(cluster->getSuffies());
		current_likelihood += _likelihood.logprobability(*dataset) ;
	}
	if (current_likelihood > _max_likelihood) {
		_max_likelihood = current_likelihood;
		_max_likelihood_membertrix = *_membertrix.clone();
	}
	foutvar(Notice) << "Loglikelihood now: " << current_likelihood << endl;
}

