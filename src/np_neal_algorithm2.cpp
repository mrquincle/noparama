#include <np_neal_algorithm2.h>

#include <iostream>
#include <iomanip>
#include <pretty_print.hpp>
#include <dim1algebra.hpp>

using namespace std;

NealAlgorithm2::NealAlgorithm2(
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
}

void NealAlgorithm2::update(
			membertrix & cluster_matrix,
			const data_ids_t & data_ids
		) {

	static int step = 0;
	fout << "Update step " << ++step << " in algorithm 2" << endl;

	assert (data_ids.size() == 1);
	data_id_t data_id = data_ids[0];
	fout << "For data item " << data_id << endl;
	
	// remove observation under consideration from cluster
	for (auto index: data_ids) {
		fout << "Retract observation " << index << " from cluster matrix" << endl;
		cluster_matrix.retract(index);
	}
	
	// existing clusters
	// if data item removed one of the clusters, this is still in there with a particular weight
	auto clusters = cluster_matrix.getClusters();

	// Create temporary data structure for K clusters
	fout << "Create temp data structures to store per cluster likelihood and indices" << endl;
	size_t K = clusters.size();
	std::vector<double> weighted_likelihood(K);
	std::vector<int> cluster_ids(K);

	data_t & observation = *cluster_matrix.getDatum(data_id);
	fout << "Current observation " << data_id << " under consideration: " << observation << endl;

	// Check all statistics given observation
	int k = 0;
	for (auto cluster_pair: clusters) {
		auto const &key = cluster_pair.first;
		auto const &cluster = cluster_pair.second;
	
		cluster_ids[k] = key;
		fout << "Obtained suffies from cluster " << key << ": " << cluster->getSuffies() << endl;
		
		_likelihood.init(cluster->getSuffies());
		weighted_likelihood[k] = _likelihood.probability(observation) * cluster_matrix.count(key);
		k++;
	}

	// Calculate parameters for new cluster
	Suffies *suffies = _nonparametrics.sample_base(_generator);
	cluster_t *new_cluster = new cluster_t(*suffies);
	
	// posterior predictive distribution is the distribution of unobserved values given observed ones
	// in this case, this probability is independent of the observed values and can be drawn from the Dirichlet
	_likelihood.init(new_cluster->getSuffies());
	double weighted_posterior_predictive = _likelihood.probability(observation) * _alpha;

	double sum_likelihoods = accumulate(weighted_likelihood.begin(), weighted_likelihood.end(), 0.0);
	
	double normalization_constant = sum_likelihoods + weighted_posterior_predictive;

	double prob_new_sample = weighted_posterior_predictive / normalization_constant;

	double reject = _distribution(_generator);

	if (reject < prob_new_sample) {
		// note, acceptance of new cluster does not depend on its likelihood
		fout << "New cluster " << endl;
		cluster_id_t cluster_index = cluster_matrix.addCluster(new_cluster);
		fout << "Add to membership matrix" << endl;
		cluster_matrix.assign(cluster_index, data_id);
		_statistics.step[0].cluster_events_accept++;
	} else {
		// pick existing cluster with weights 
		fout << "Existing cluster" << endl;
		size_t index = algebra::random_weighted_pick(weighted_likelihood.begin(), weighted_likelihood.end(), _generator);
		cluster_id_t cluster_id = cluster_ids[index];
		
		fout << "Assign data to cluster with id " << cluster_id << endl;
		assert (cluster_matrix.count(cluster_id) != 0);

		cluster_matrix.assign(cluster_id, data_id);
		//cluster_matrix.reassign(index, data_id);
		_statistics.step[0].cluster_events_reject++;
	}
}

void NealAlgorithm2::printStatistics() {
	int verbosity = _verbosity;
	_verbosity = 0;
	fout << "Statistics:" << endl;
	fout << " # of new cluster events accepted: " << _statistics.step[0].cluster_events_accept << endl;
	fout << " # of new cluster events rejected: " << _statistics.step[0].cluster_events_reject << endl;
	_verbosity = verbosity;
} 
