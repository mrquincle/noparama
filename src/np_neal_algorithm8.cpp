#include <np_neal_algorithm8.h>

#include <iostream>
#include <iomanip>
#include <pretty_print.hpp>
#include <dim1algebra.hpp>

using namespace std;

NealAlgorithm8::NealAlgorithm8(
			random_engine_t & generator,
			distribution_t & likelihood,
			dirichlet_process & nonparametrics
		): 
			_generator(generator),
			_likelihood(likelihood),
			_nonparametrics(nonparametrics)
{
	_alpha = _nonparametrics.getSuffies().alpha;

	_verbosity = 5;

	fout << "likelihood: " << _likelihood.getSuffies() << endl;
	
	// the number of auxiliary variables, M=1 should be the same as algorithm 2
	_M = 3;
}

/*
 * This function assigns the given data item and updates the number of clusters if it is assigned to a new cluster.
 * It does not reassign other data items. Hence, it is a single update step. This naturally has disadvantages for
 * problems where we might converge faster if we would be able to re-assign a set of data points at once. The latter
 * is the idea behind split and merge samplers.
 *
 * The partition under consideration emerges from calling this update() function regularly. The data points only
 * enter the equation through the likelihood function where the number of data points per cluster determines their
 * relative weight. This can be seen as some kind of voting mechanism. However, there is no "guidance" of the chain
 * to partitions of data that might make more sense. The MCMC is not data-driven. Sampling of the cluster parameters
 * is a not data-driven either. We just sample from the base distribution "_nonparametrics.sample_base" without 
 * regard for that part of parameter space that makes sense from the perspective of the data.
 */
void NealAlgorithm8::update(
			membertrix & cluster_matrix,
			const data_ids_t & data_ids
		) {
	
	assert (data_ids.size() == 1);
	data_id_t data_id = data_ids[0];

	static int step = 0;
	fout << "Update step " << ++step << " in algorithm 8 for data item " << data_id << endl;
	
	// remove observation under consideration from cluster
	for (auto index: data_ids) {
		fout << "Retract observation " << index << " from cluster matrix" << endl;
		cluster_matrix.retract(index);
	}

	// existing clusters
	// if data item removed one of the clusters, this is still in there with a particular weight
	auto clusters = cluster_matrix.getClusters();

	// Create temporary data structure for K clusters
	size_t K = clusters.size();
	std::vector<int> cluster_ids(K);

	fout << "Calculate likelihood for K=" << K << " existing clusters" << endl;

	fout << "Calculate parameters for M new clusters" << endl;
	// Calculate parameters for M new clusters
	clusters_t new_clusters(_M);
	for (int m = 0; m < _M; ++m) {
		Suffies *suffies = _nonparametrics.sample_base(_generator);
		cluster_t *temp = new cluster_t(*suffies);
		new_clusters[m] = temp;
		fout << "Suffies generated: " << new_clusters[m]->getSuffies() << endl;
	}	
	
	data_t & observation = *cluster_matrix.getDatum(data_id);
	fout << "Current observation " << data_id << " under consideration: " << observation << endl;

	// Get (weighted) likelihoods for each existing and new cluster given the above observation
	std::vector<double> weighted_likelihood(K+_M);

	int k = 0;
	for (auto cluster_pair: clusters) {
		auto const &key = cluster_pair.first;
		auto const &cluster = cluster_pair.second;
	
		cluster_ids[k] = key;

		fout << "Obtained suffies from cluster " << key << ": " << cluster->getSuffies() << endl;

		// here _likelihood is sudden a niw distribution...
		// TODO: is this indeed the case? or an old debug statement?
		_likelihood.init(cluster->getSuffies());
		fout << "calculate probability " << endl;
		double ll = _likelihood.probability(observation);
		fout << "which is: " << ll << endl;
		weighted_likelihood[k] = _likelihood.probability(observation) * cluster_matrix.count(key);
		k++;
	}

	/*
	 * The Metropolis-Hastings step with generating M proposals at the same time does not look like the typical 
	 * alpha < u rule, however, it is. With a single proposal we can first sample for accept when it falls under
	 * a uniform random variable. And then select one of the existing proposals with a probability proportional to
	 * their (weighted) likelihoods. To remind you, the weight is determined by the number of observations tied to
	 * that cluster for existing clusters and alpha for proposed clusters. With multiple proposals we can just directly
	 * sample from the complete vector (existing and proposed clusters) with weighted likelihoods.
	 */
	for (int m = 0; m < _M; ++m) {
		_likelihood.init(new_clusters[m]->getSuffies());
		weighted_likelihood[K+m] = _likelihood.probability(observation) * _alpha / (double)_M;
		fout << "New cluster " << m << 
			"[~#" << _alpha/_M << "]: " << \
			_likelihood.probability(observation) << \
			new_clusters[m]->getSuffies() << endl;
	}

	// sample uniformly from the vector "weighted_likelihood" according to the weights in the vector
	// hence [0.5 0.25 0.25] will be sampled in the ratio 2:1:1, and return index 0 1 2 accordingly
	size_t index = algebra::random_weighted_pick(weighted_likelihood.begin(), weighted_likelihood.end(), _generator);

	/*
	 * Pick either a new or old cluster using calculated values. With a new cluster generate new sufficient 
	 * statistics and update the membership matrix.
	 */
	if (index >= K) {
		// pick new cluster
		fout << "New cluster: " << index-K << endl;
		// lazy, just allocate a new one with proper sufficient statistics and deallocate all temporary ones at the end
		cluster_t *new_cluster = new cluster_t(new_clusters[index-K]->getSuffies());
		fout << "Add to membership matrix" << endl;
		cluster_id_t cluster_index = cluster_matrix.addCluster(new_cluster);
		fout << "Assign data to cluster with id " << cluster_index << endl;
		cluster_matrix.assign(cluster_index, data_id);
		_statistics.step[0].cluster_events_accept++;
	} else {
		// pick existing cluster with given cluster_id
		fout << "Existing cluster" << endl;
		
		cluster_id_t cluster_id = cluster_ids[index];

		fout << "Assign data to cluster with id " << cluster_id << endl;
		assert (cluster_matrix.count(cluster_id) != 0);
		
		cluster_matrix.assign(cluster_id, data_id);
		_statistics.step[0].cluster_events_reject++;
	}
	
	// deallocate new clusters that have not been used
	fout << "Deallocate temporary clusters" << endl;
	for (int m = 0; m < _M; ++m) {
		delete new_clusters[m];
	}

	// check
	assert(cluster_matrix.assigned(data_id));
}

void NealAlgorithm8::printStatistics() {
	int verbosity = _verbosity;
	_verbosity = 0;
	fout << "Statistics:" << endl;
	fout << " # of new cluster events accepted: " << _statistics.step[0].cluster_events_accept << endl;
	fout << " # of new cluster events rejected: " << _statistics.step[0].cluster_events_reject << endl;
	_verbosity = verbosity;
} 
