#include <np_jain_neal_algorithm.h>

#include <iostream>
#include <iomanip>
#include <pretty_print.hpp>
#include <dim1algebra.hpp>
#include <cmath>

using namespace std;

// note, you'll need g++-7 and C++17 to obtain special_math function beta.

#define USE_LOGARITHM

/*
 * This algorithm is a split-merge sampler. 
 */
JainNealAlgorithm::JainNealAlgorithm(
			random_engine_t & generator,
			distribution_t & likelihood,
			dirichlet_process & nonparametrics
		): 
			_generator(generator),
			_likelihood(likelihood),
			_nonparametrics(nonparametrics)
{
	_alpha = _nonparametrics.getSuffies().alpha;

	_verbosity = Debug;
	_verbosity = Warning;
	
	_statistics.step[0].type = "merge";
	_statistics.step[1].type = "split";

	fout << "likelihood: " << _likelihood.getSuffies() << endl;

	_split_method = sams_prior;
	//_split_method = simple_random_split;
}

JainNealAlgorithm::~JainNealAlgorithm() {
}

double JainNealAlgorithm::ratioStateProb(bool split, int nc0, int nc1, int nc) {
#ifdef USE_LOGARITHM
	assert (nc == nc0 + nc1);
	double result = std::log(_alpha) + lgamma(nc0) + lgamma(nc1) - lgamma(nc);
	if (split)
		return result;
	else
		return -result;
#else
	double result = _alpha * beta(nc0, nc1);
	if (split)
		return result;
	else
		return 1.0 / result;
#endif
}

double JainNealAlgorithm::ratioProposal(bool split, int N, int C) {
#ifdef USE_LOGARITHM
	if (_split_method == sams_prior) {
		return 0.0;
	}
	double result1 = (C - N - 1) * std::log(C - 1) - (C - N) * std::log(C);
	assert (C == 2);
	double result = std::log(pow(2, N - 2));

	assert(result1 == result);

	if (split)
		return result;
	else
		return -result;
#else
	if (_split_method == sams_prior) {
		return 1.0;
	}
	// when assuming C=2:
	// double result = pow(2, N - 2);
	assert (C == 2);
	double result = pow(2, N - 2);
	if (split) 
		return result;
	else
		return 1.0 / result;
#endif
}

void JainNealAlgorithm::propose_split(data_ids_t & move, data_ids_t & remain, data_id_t data_i, data_id_t data_j, 
		cluster_id_t current_cluster_id, cluster_t & new_cluster, split_method_t split_method) {

	// assign data_i and data_j to separate clusters
	move.push_back(data_i);
	remain.push_back(data_j);
	
	// current cluster
	auto current_cluster = _cluster_matrix->getCluster(current_cluster_id);

	// get assignments (just ids, not the data itself)
	data_ids_t all_data;
	_cluster_matrix->getAssignments(current_cluster_id, all_data);
	fout << "The number of data points assigned to the cluster now is " << all_data.size() << endl;
	
	// shuffle the data items
	algebra::random_order(all_data.begin(), all_data.end());

	// for random_mixing method prepare cut_off parameter
	int cut_off;
	if (split_method == random_mixing) {
		double u = _distribution(_generator);
		cut_off = all_data.size() * u;
	}

	// assign all data points
	for (int i = 0; i < (int)all_data.size(); ++i) {
		int data_id = all_data[i];
		if (data_id == data_i) continue;
		if (data_id == data_j) continue;

		switch (split_method) {
			case random_mixing:
				{
					if (i < cut_off) {
						move.push_back(data_id);
					} else {
						remain.push_back(data_id);
					}
					break;
				}
			case simple_random_split:  
				{
					// sample uniform random variable to randomly split over clusters
					double toss = _distribution(_generator);
					if (toss <= 0.5) {
						// move these points
						move.push_back(data_id);
					} else {
						// keep this points at same cluster
						remain.push_back(data_id);
					}
					break;
				} 
			case sams_prior: 
				{
					// Observation: a lot of times there will be no single extra data point attached to the new cluster
					// because it is so obviously bad w.r.t. likelihood
					// shouldn't we at least sample more options as with the m auxiliary proposals?

					// SAMS
					std::vector<double> ldata(2);
					data_t * datum = _cluster_matrix->getDatum(data_id);
					
					_likelihood.init(current_cluster->getSuffies());
					//cout << "Get probability for " << *datum << endl;
#ifdef USE_LOGARITHM
					ldata[0] = _likelihood.logprobability(*datum) + remain.size();
#else
					ldata[0] = _likelihood.probability(*datum) * remain.size();
#endif

					_likelihood.init(new_cluster.getSuffies());
#ifdef USE_LOGARITHM
					ldata[1] = _likelihood.logprobability(*datum) + move.size();
#else
					ldata[1] = _likelihood.probability(*datum) * move.size();
#endif
					int index = algebra::random_weighted_pick(ldata.begin(), ldata.end(), _generator);
					if (index == 0) {
						remain.push_back(data_id);
					} else {
						move.push_back(data_id);
					}
					break;
				}
			case sams_random_walk: 
				{
					// TODO
					assert(false);
					break;
				}
		}
	}
}

void JainNealAlgorithm::checkLikelihoods(double lsrc, double ldest, step_t statistics_step, double &lratio, bool &accept, bool &overwrite) {
	overwrite = false;
	if (ldest == 0 && lsrc == 0) {
		fout << "Both likelihood of existing and new cluster is zero. Reject it." << endl;
		overwrite = true;
		accept = false;
		statistics_step.likelihood_both_zero++;
	} else if (lsrc == 0) {
		fout << "Likelihood of source cluster is zero. Accept new cluster." << endl;
		overwrite = true;
		accept = true;
		statistics_step.source_likelihood_zero++;
	} else if (ldest == 0) {
		fout << "Likelihood of target cluster is " << ldest << ". Reject it in all cases." << endl;
		fout << "Just for the record, likelihood of source is: " << lsrc << endl;
		overwrite = true;
		accept = false;
		statistics_step.target_likelihood_zero++;
	} else {
		lratio = ldest / lsrc; 
		fout << "Likelihood ratio is something useful " << lratio << endl;
		statistics_step.likelihood_both_nonzero++;
	}
}

/**
 * Here we sample the parameters from the new cluster that we need for the split. The split is not yet accepted!
 */
bool JainNealAlgorithm::split(
		data_ids_t data_ids,
		cluster_id_t current_cluster_id
		) {

	data_ids_t move, remain;
	dataset_t data_to_move;
	int nc0, nc1, nc;
	double p21, q12;
	double lsrc, ldest, lratio;
	bool accept = true, overwrite = false;

	_statistics.step[1].attempts++;
	
	// current cluster
	auto current_cluster = _cluster_matrix->getCluster(current_cluster_id);
	
	// new cluster
	Suffies *suffies = _nonparametrics.sample_base(_generator);
	cluster_t *new_cluster = new cluster_t(*suffies);

	propose_split(move, remain, data_ids[0], data_ids[1], current_cluster_id, *new_cluster, _split_method);

	nc0 = move.size();
	nc1 = remain.size();
	nc = nc0 + nc1; 
	fout << "Cluster is size " << nc << " and after split would become size " << nc0 << " and " << nc1 << endl;

	p21 = ratioStateProb(true, nc0, nc1, nc);
	assert ( ratioStateProb(true, nc1, nc0, nc) == p21 );
	q12 = ratioProposal(true, nc);

#ifdef USE_LOGARITHM
	fout << "The q(.|.)p(.) ratio for the MH split is " << q12 << " + " << p21 << " = " << q12 + p21 << endl;
#else
	fout << "The q(.|.)p(.) ratio for the MH split is " << q12 << " * " << p21 << " = " << q12 * p21 << endl;
#endif
	// only consider the data points that move to the new one, but at the old cluster
	_cluster_matrix->getData(move, data_to_move);

	// likelihood if data potentially to be moved would also remain at current cluster
	_likelihood.init(current_cluster->getSuffies());
#ifdef USE_LOGARITHM
	lsrc = _likelihood.logprobability(data_to_move);
#else
	lsrc = _likelihood.probability(data_to_move);
#endif

	// likelihood if data would move to new cluster
	_likelihood.init(new_cluster->getSuffies());
#ifdef USE_LOGARITHM
	ldest = _likelihood.logprobability(data_to_move);
#else
	ldest = _likelihood.probability(data_to_move);
#endif

#ifdef USE_LOGARITHM
	lratio = ldest - lsrc;
#else
	// handle weird corner cases
	checkLikelihoods(lsrc, ldest, _statistics.step[1], lratio, accept, overwrite);
#endif
	// normal MH step
	if (!overwrite) {
#ifdef USE_LOGARITHM
		double a_split = std::exp(q12 + p21 + lratio);
#else
		double a_split = q12 * p21 * lratio;
#endif
		fout << "Acceptance for split: " << a_split << endl;

		// sample uniform random variable
		double u = _distribution(_generator);
		fout << "Sampled u: " << u << endl;

		if (a_split < u) {
			_statistics.step[1].likelihood_both_nonzero_reject++;
			accept = false;
		} else {
			_statistics.step[1].likelihood_both_nonzero_accept++;
			accept = true;
		}
	}

	if (accept) {
		fout << "Accept split!" << endl;
		// actually perform the move: only points that move have to be moved...
		cluster_id_t new_cluster_id = _cluster_matrix->addCluster(new_cluster);
		fout << "Add cluster with id " << new_cluster_id << endl;				
		for (auto data: move) {
			fout << "Add data " << data << " to the new cluster " << new_cluster_id << endl;
			_cluster_matrix->retract(data);
			_cluster_matrix->assign(new_cluster_id, data);
		}
	} else {
		fout << "Reject split!" << endl;
		delete new_cluster;
	}
	return accept;
}

bool JainNealAlgorithm::merge(
		data_ids_t data_ids,
		std::vector<cluster_id_t> & current_clusters
		) {

	data_ids_t data_ids0, data_ids1; 
	dataset_t *data_to_move;
	int nc0, nc1, nc;
	double p12, q21;
	double lsrc, ldest, lratio;
	bool accept = true, overwrite = false;

	_statistics.step[0].attempts++;

	// get assignments (just ids, not the data itself)
	_cluster_matrix->getAssignments(current_clusters[0], data_ids0);
	_cluster_matrix->getAssignments(current_clusters[1], data_ids1);

	nc0 = data_ids0.size();
	nc1 = data_ids1.size();
	nc = nc0 + nc1;
	fout << "Clusters are size " << nc0 << " and " << nc1 << " and after merge become " << nc << endl;

	p12 = ratioStateProb(false, nc0, nc1, nc);
	assert ( ratioStateProb(false, nc1, nc0, nc) == p12 );
	q21 = ratioProposal(false, nc);
#ifdef USE_LOGARITHM
	fout << "The q(.|.)p(.) ratio for the MH merge is " << q21 << " + " << p12 << " = " << q21 + p12 << endl;
#else
	fout << "The q(.|.)p(.) ratio for the MH merge is " << q21 << " * " << p12 << " = " << q21 * p12 << endl;
#endif	
	// data to move sits in cluster 0
	data_to_move = _cluster_matrix->getData(current_clusters[0]);

	// likelihood of data in cluster 0 before merge
	auto cluster0 = _cluster_matrix->getCluster(current_clusters[0]);
	_likelihood.init(cluster0->getSuffies());
#ifdef USE_LOGARITHM
	lsrc = _likelihood.logprobability(*data_to_move);
#else
	lsrc = _likelihood.probability(*data_to_move);
#endif
//	fout << "Likelihood of data in cluster 0 before merge: l2_0 = " << l2_0 << endl;

	// likelihood of data if moved to cluster 1 
	auto cluster1 = _cluster_matrix->getCluster(current_clusters[1]);
	_likelihood.init(cluster1->getSuffies());
#ifdef USE_LOGARITHM
	ldest = _likelihood.logprobability(*data_to_move);
#else
	ldest = _likelihood.probability(*data_to_move);
#endif
//	fout << "Likelihood of same data in case of moving them to cluster 1 (after merge): l1_0 = " << l1_0 << endl;

#ifdef USE_LOGARITHM
	lratio = ldest - lsrc;
#else
	checkLikelihoods(lsrc, ldest, _statistics.step[0], lratio, accept, overwrite);
#endif

	if (!overwrite) {
#ifdef USE_LOGARITHM
		double a_merge = std::exp(q21 + p12 + lratio); 
#else
		double a_merge = q21 * p12 * lratio; 
#endif
		fout << "Acceptance for merge: " << a_merge << endl;

		// sample uniform random variable
		double u = _distribution(_generator);
		fout << "Sampled u: " << u << endl;

		if (a_merge < u) {
			_statistics.step[0].likelihood_both_nonzero_reject++;
			accept = false;
		} else {
			_statistics.step[0].likelihood_both_nonzero_accept++;
			accept = true;
		}
	}

	if (accept) {
		fout << "Accept merge!" << endl;
		// actually perform the move: all items from i go to j
		for (auto data: data_ids0) {
			_cluster_matrix->retract(data);
			_cluster_matrix->assign(current_clusters[1], data);
		}

	} else {
		fout << "Reject merge!" << endl;
	}
	return accept;
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
void JainNealAlgorithm::update(
			membertrix & cluster_matrix,
			const data_ids_t & data_ids
		) {
	// there should be exactly two data points sampled
	assert (data_ids.size() == 2);

	// the data ids should be unique
	assert (data_ids[0] != data_ids[1]);

	_cluster_matrix = &cluster_matrix;

	static int step = 0;
	fout << "Update step " << ++step << " in Jain & Neal's split-merge algorithm" << endl;

	// store cluster indices of current observations
	std::vector<cluster_id_t> current_clusters; 
	current_clusters.clear();
	for (auto data_id: data_ids) {
		cluster_id_t cluster_id = _cluster_matrix->getClusterId(data_id);
		fout << "We found data item " << data_id << " at cluster " << cluster_id << endl;
		current_clusters.push_back( cluster_id );
	}

	int uniq_cluster_count = algebra::count_unique(current_clusters.begin(), current_clusters.end());
	fout << "The data comes from " << uniq_cluster_count << " unique cluster(s)" << endl;
	
	// a single cluster
	switch (uniq_cluster_count) {
		case 0: default:
			assert(false);
			break;
		case 1: 
			{ // split step
				fout << "We are gonna split this single cluster" << endl;
				assert (current_clusters[0] == current_clusters[1]);

				size_t Ka = _cluster_matrix->getClusterCount();

				int current_cluster = current_clusters[0];
				bool accept = split(data_ids, current_cluster);

				if (accept) {
					size_t Kb = _cluster_matrix->getClusterCount();
					assert (Kb == Ka + 1);

					_statistics.step[1].cluster_events_accept++;
				} else {
					_statistics.step[1].cluster_events_reject++;
				}
				break;
			}
		case 2: 
			{ // merge step
				fout << "Merge these clusters" << endl;
				assert (current_clusters[0] != current_clusters[1]);

				size_t Ka = _cluster_matrix->getClusterCount();
				
				bool accept = merge(data_ids, current_clusters);

				if (accept) {
					size_t Kb = _cluster_matrix->getClusterCount();
					assert (Kb == Ka - 1);
					
					_statistics.step[0].cluster_events_accept++;
				} else {
					_statistics.step[0].cluster_events_reject++;
				}
				break;
			}
	} 

	// check
	for (auto data_id: data_ids) {
		fout << "Check data id " << data_id << endl;
		assert(_cluster_matrix->assigned(data_id));
	}
}

void JainNealAlgorithm::printStatistics() {
	int verbosity = _verbosity;
	_verbosity = Debug;

	fout << endl;
	fout << "Statistics:" << endl;
	for (int i = 0; i < 2; ++i) {
		step_t step = _statistics.step[i];
		fout << " # of " << step.type << " attempts: " << step.attempts << endl;
		fout << "   o of accepted " << step.type << " cluster events: " << step.cluster_events_accept << endl;
		fout << "     - of accept 0->? overrides: " << step.source_likelihood_zero << endl;
		fout << "     - of events accepted after consideration: " << step.likelihood_both_nonzero_accept << endl;

		fout << "   o of rejected " << step.type << " cluster events: " << step.cluster_events_reject << endl;
		fout << "   # of reject ?->0 overrides: " << step.target_likelihood_zero << endl;

		fout << "   # of metropolis-hastings decisions: " << step.likelihood_both_nonzero << endl;
		fout << "     - accepted: " << step.likelihood_both_nonzero_accept << endl;
		fout << "     - rejected: " << step.likelihood_both_nonzero_reject << endl;
	}

	// set back to zero
	_statistics = {};
	_statistics.step[0].type = "merge";
	_statistics.step[1].type = "split";

	_verbosity = verbosity;
} 
