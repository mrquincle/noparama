#include <np_jain_neal_algorithm.h>

#include <iostream>
#include <iomanip>
#include <pretty_print.hpp>
#include <dim1algebra.hpp>
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <cmath>

using namespace std;

enum { simple_random_split, sams_prior, sams_random_walk };

#define JAIN_NEAL_Q_CALCULATION

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
	
	_statistics = (statistics_t){0};

	fout << "likelihood: " << _likelihood.getSuffies() << endl;
		
}

bool JainNealAlgorithm::split(
		membertrix & cluster_matrix,
		data_ids_t data_ids,
		std::vector<cluster_id_t> & prev_clusters
		) {
	bool accept;
	//_verbosity = Warning;

	_statistics.split_attempts++;

	// get assignments (just ids, not the data itself)
	data_ids_t data_ids1;
	data_ids1.clear();
	cluster_matrix.getAssignments(prev_clusters[0], data_ids1);
	fout << "The number of data points assigned to the cluster now is " << data_ids1.size() << endl;
	fout << "The data items assigned to the cluster are: " << data_ids1 << endl;

	data_ids_t data_ids2_0 = { data_ids[0] };
	data_ids_t data_ids2_1 = { data_ids[1] };

//	int split_method = sams_prior;
	int split_method = simple_random_split;

	// shuffle the data items
	algebra::random_order(data_ids1.begin(), data_ids1.end());

	// old and new cluster
	auto cluster0 = cluster_matrix.getCluster(prev_clusters[0]);

	int M = 5; 
	int m = 0;
	std::vector<cluster_t*> new_clusters;
	for (int i = 0; i < M; ++i) {
		Suffies *suffies = _nonparametrics.sample_base(_generator);
		cluster_t *new_cluster = new cluster_t(*suffies);
		new_clusters.push_back(new_cluster);
	}

	switch (split_method) {
		case simple_random_split:  
			{
				fout << "We will use a simple random selection to divide the data over the existing and new cluster." << endl;
				break;

			} 
		case sams_prior: case sams_random_walk:
			{
				fout << "We will use SAMS to divide the data over the existing and new cluster intelligently." << endl;
				fout << "The existing cluster has as parameters " << cluster0->getSuffies() << endl;
				fout << "The new cluster m=" << m << " has as parameters " << new_clusters[m]->getSuffies() << endl;
				break;
			}
	}
	
	// set m using only first data point that is not i nor j
	bool set_m = false;
	for (auto data_id: data_ids1) {
		if (data_id == data_ids[0]) continue;
		if (data_id == data_ids[1]) continue;

		switch (split_method) {
			case sams_prior: 
				{
					std::vector<double> likelihoods(M);
					for (int i = 0; i < M; ++i) {
						_likelihood.init(new_clusters[i]->getSuffies());
						data_t * datum = cluster_matrix.getDatum(data_id);
						likelihoods[i] = _likelihood.probability(*datum);
					}
					m = algebra::random_weighted_pick(likelihoods.begin(), likelihoods.end(), _generator);
					//cout << "Likelihoods: " << likelihoods << endl;
					//cout << "Picked: " << m << endl;
					set_m = true;
					break;
				}
		}
		if (set_m) break;
	}

	for (auto data_id: data_ids1) {
		if (data_id == data_ids[0]) continue;
		if (data_id == data_ids[1]) continue;

		switch (split_method) {
			case simple_random_split:  
				{
					// sample uniform random variable to randomly split over clusters
					double toss = _distribution(_generator);
					if (toss <= 0.5) {
						// move these points
						data_ids2_0.push_back(data_id);
					} else {
						data_ids2_1.push_back(data_id);
					}
					break;
				} 
			case sams_prior: 
				{
					// Observation: a lot of times there will be no single extra data point attached to the new cluster
					// because it is so obviously bad w.r.t. likelihood
					// shouldn't we at least sample more options as with the m auxiliary proposals?

					// SAMS
					_likelihood.init(new_clusters[m]->getSuffies());
					data_t * datum = cluster_matrix.getDatum(data_id);
					//fout << "Consider datum: " << data_id << ": " << *datum << endl;
					double l2_0 = _likelihood.probability(*datum) * data_ids2_0.size();

					//fout << "It has a likelihood at the new cluster of " << l2_0 << endl;
					_likelihood.init(cluster0->getSuffies());

					double l2_1 = _likelihood.probability(*datum) * data_ids2_1.size();
					//fout << "It has a likelihood at the existing cluster of " << l2_1 << endl;

					double frac;
					if (l2_0 + l2_1 == 0) {
						frac = 1/2;
					} else {
						frac = l2_0 / (l2_0 + l2_1);
					}
					//cout << "frac: " << frac << endl;
					double u = _distribution(_generator);

					if (frac > u) {
						//fout << "Add " << data_id << " to new cluster" << endl;
						data_ids2_0.push_back(data_id);
					} else {
						//fout << "Keep " << data_id << " at existing cluster" << endl;
						data_ids2_1.push_back(data_id);
					}

					break;
				}
			case sams_random_walk: 
				{
					assert(false);
					break;
				}
		}
	}

	int nc0 = data_ids2_0.size();
	assert (nc0 > 0);
	int nc1 = data_ids2_1.size();
	assert (nc1 > 0);
	int nc = nc0 + nc1; // TODO: check again size
	fout << "Cluster is size " << nc << " and after split would become size " << nc0 << " and " << nc1 << endl;

	// note, you'll need g++-7 and C++17 to obtain special_math function beta.
	double p21 = _alpha * beta(nc0, nc1);
	// q12 is quite large for a split move
	double q12;
#ifdef JAIN_NEAL_Q_CALCULATION
	q12 = pow(2, nc - 2);
#else
	// q12 is much smaller, but this takes into account n_i and n_j. I's not equally likely to sample (.. .. | .. .. ..) 
	// as to sample (.. | .. .. .. ..) isn't it?
	int f1 = 2;
	if (nc0 == nc1) f1 = 1;
	q12 = 1 / (f1 * (nc0*nc1)/( 6.0 * (nc0+nc1-1)*(nc0+nc1)*(nc0+nc1+1)  ) );
#endif

	fout << "The ratio to be considered for a Metropolis-Hastings split " << nc0 << ", " << nc1 << " without the likelihood ratio is " << q12 << " * " << p21 << " = " << q12 * p21 << endl;

	// only consider the data points that move to the new one, but at the old cluster
	dataset_t data0;
	data0.clear();
	cluster_matrix.getData(data_ids2_0, data0);
	_likelihood.init(cluster0->getSuffies());
	double l1 = _likelihood.probability(data0);

	// at new cluster consider same dataset there
	_likelihood.init(new_clusters[m]->getSuffies());
	double l2_0 = _likelihood.probability(data0);

	double l21; 
	bool overwrite = false;

	if (l2_0 == 0 && l1 == 0) {
		fout << "Both likelihood of existing and new cluster is zero. Reject it." << endl;
		overwrite = true;
		accept = false;
		_statistics.split_likelihood_both_zero++;
	} else if (l1 == 0) {
		fout << "Likelihood of source cluster is zero. Accept new cluster." << endl;
		overwrite = true;
		accept = true;
		_statistics.split_source_likelihood_zero++;
	} else if (l2_0 == 0) {
		fout << "Likelihood of target cluster is " << l2_0 << ". Reject it in all cases." << endl;
		fout << "Just for the record, likelihood of source is: " << l1 << endl;
		overwrite = true;
		accept = false;
		_statistics.split_target_likelihood_zero++;
	} else {
		l21 = l2_0 / l1; 
		fout << "Likelihood ratio is something useful " << l21 << endl;
		_statistics.split_likelihood_both_nonzero++;
	}

	/*
	// if the point we randomly picked is the only one under consideration, overwrite everything and reject it
	if (data_ids2_0.size() <= 1) {
		overwrite = true;
		accept = false;
	} */

	if (!overwrite) {
		double a_split = q12 * p21 * l21; // or actually min(1, ...) but we compare with u ~ U(0,1) anyway
		fout << "Acceptance for split: " << a_split << endl;

		// sample uniform random variable
		double u = _distribution(_generator);
		fout << "Sampled u: " << u << endl;

		if (a_split < u) {
			_statistics.split_likelihood_both_nonzero_reject++;
			//cout << "Reject: clusters sizes " << data_ids2_0.size() << " and " << data_ids2_1.size() << endl;
			accept = false;
		} else {
			//cout << "Accept: clusters sizes " << data_ids2_0.size() << " and " << data_ids2_1.size() << endl;
			_statistics.split_likelihood_both_nonzero_accept++;
			accept = true;
		}
	}

	if (accept) {
		fout << "--------------------------------------------------------------------------------------------" << endl;
		fout << "Accept split!" << endl;
		// actually perform the move: only points that move have to be moved...
		cluster_id_t new_cluster_id = cluster_matrix.addCluster(new_clusters[m]);
		fout << "Add cluster with id " << new_cluster_id << endl;				
		for (auto data: data_ids2_0) {
			fout << "Add data " << data << " to the new cluster " << new_cluster_id << endl;
			cluster_matrix.retract(data);
			cluster_matrix.assign(new_cluster_id, data);
		}
		// remove clusters that have not been picked
		for (int i = 0; i < M; ++i) {
			if (i != m) delete new_clusters[i];
		}
		// points in data_ids2_1 need not be moved
	} else {
		fout << "Reject split!" << endl;
		// remove all clusters
		for (int i = 0; i < M; ++i) {
			delete new_clusters[i];
		}
	}
	// TODO: anything that needs to be re-established on rejection?
	//_verbosity = Debug;
	return accept;
}

bool JainNealAlgorithm::merge(
		membertrix & cluster_matrix,
		data_ids_t data_ids,
		std::vector<cluster_id_t> & prev_clusters
		) {
	bool accept;

	_statistics.merge_attempts++;

	// get assignments (just ids, not the data itself)
	data_ids_t data_ids0; 
	data_ids0.clear();
	cluster_matrix.getAssignments(prev_clusters[0], data_ids0);
	data_ids_t data_ids1; 
	data_ids1.clear();
	cluster_matrix.getAssignments(prev_clusters[1], data_ids1);

	int nc0 = data_ids0.size();
	int nc1 = data_ids1.size();
	int nc = nc0 + nc1;
	fout << "Clusters are size " << nc0 << " and " << nc1 << " and after merge become " << nc << endl;

	// check if type is correct (no rounding off to ints by accident)
	double p12 = 1.0/(_alpha * beta(nc0, nc1));
	double q21;
#ifdef JAIN_NEAL_Q_CALCULATION
	q21 = pow(0.5, nc - 2);
#else
	int f1 = 2;
	if (nc0 == nc1) f1 = 1;
	q21 = f1 * (nc0*nc1)/( 6.0 * (nc0+nc1-1)*(nc0+nc1)*(nc0+nc1+1)  );
#endif
	fout << "The ratio to be considered for Metropolis-Hastings without the likelihood ratio is " << q21 << " * " << p12 << " = " << q21 * p12 << endl;

	// calculate likelihoods
	// if we would have stored the likelihoods somewhere of the original clusters we might have reused them
	// now we have to calculate them from scratch

	// likelihood of data in cluster 0 before merge
	auto cluster0 = cluster_matrix.getCluster(prev_clusters[0]);
	_likelihood.init(cluster0->getSuffies());
	dataset_t * data0 = cluster_matrix.getData(prev_clusters[0]);
	double l2_0 = _likelihood.probability(*data0);
	fout << "Likelihood of data in cluster 0 before merge " << l2_0 << endl;

	// likelihood of data if moved to cluster 1 
	auto cluster1 = cluster_matrix.getCluster(prev_clusters[1]);
	_likelihood.init(cluster1->getSuffies());
	double l1_0 = _likelihood.probability(*data0);
	fout << "Likelihood of same data in case of moving them to cluster 1 (after merge) " << l1_0 << endl;

	// likelihood of data in cluster 1 
	dataset_t * data1 = cluster_matrix.getData(prev_clusters[1]);
	double l2_1 = _likelihood.probability(*data1);
	fout << "Likelihood of data in cluster 1 before merge " << l2_1 << endl;

	fout << "Likelihood of cluster 1 and cluster 2 before merge " << l2_0 * l2_1 << endl;
	// likelihood of all data in cluster 1
	dataset_t data01;
	for (auto d: *data0) {
		data01.push_back(d);
	}
	for (auto d: *data1) {
		data01.push_back(d);
	}
	double l1 = _likelihood.probability(data01);
	fout << "Likelihood of all data after merge " << l1 << endl;

	double l12;
	bool overwrite = false;
	if (l1_0 == 0 && l2_0 == 0) {
		// likelihood is zero before and after merge... we should just walk...
		overwrite = true;
		accept = false;
		_statistics.merge_likelihood_both_zero++;
	} else if (l2_0 == 0) {
		// just accept, source is 0, nothing to loose
		overwrite = true;
		accept = true;
		_statistics.merge_source_likelihood_zero++;
	} else if (l1_0 == 0) {
		// just reject
		overwrite = true;
		accept = false;
		_statistics.merge_target_likelihood_zero++;
	} else {
		l12 = l1_0 / l2_0;
		fout << "Likelihood ratio becomes: " << l12 << endl;
		_statistics.merge_likelihood_both_nonzero++;
	}

	if (!overwrite) {
		double a_merge = q21 * p12 * l12; // or actually min(1, ...) but we compare with u ~ U(0,1) anyway
		fout << "Acceptance for merge: " << a_merge << endl;

		// sample uniform random variable
		double u = _distribution(_generator);
		fout << "Sampled u: " << u << endl;

		if (a_merge < u) {
			_statistics.merge_likelihood_both_nonzero_reject++;
			accept = false;
		} else {
			_statistics.merge_likelihood_both_nonzero_accept++;
			accept = true;
		}
	}

	if (accept) {
		fout << "--------------------------------------------------------------------------------------------" << endl;
		fout << "Accept merge!" << endl;
		// actually perform the move: all items from i go to j
		for (auto data: data_ids0) {
			cluster_matrix.retract(data);
			cluster_matrix.assign(prev_clusters[1], data);
		}

	} else {
		fout << "Reject merge!" << endl;
	}
	// TODO: anything that needs to be re-established on rejection?
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
			data_ids_t data_ids
		) {
	// there should be exactly two data points sampled
	assert (data_ids.size() == 2);

	// the data ids should be unique
	assert (data_ids[0] != data_ids[1]);

	static int step = 0;
	fout << "Update step " << ++step << " in Jain & Neal's split-merge algorithm" << endl;

	// store cluster indices of current observations
	std::vector<cluster_id_t> prev_clusters; 
	prev_clusters.clear();
	for (auto data_id: data_ids) {
		cluster_id_t cluster_id = cluster_matrix.getClusterId(data_id);
		fout << "We found data item " << data_id << " at cluster " << cluster_id << endl;
		prev_clusters.push_back( cluster_id );
	}

	int uniq_cluster_count = algebra::count_unique(prev_clusters.begin(), prev_clusters.end());
	fout << "The data comes from " << uniq_cluster_count << " unique cluster(s)" << endl;

	// a single cluster
	switch (uniq_cluster_count) {
		case 0: default:
			assert(false);
			break;
		case 1: 
			{ // split step
				fout << "We are gonna split this single cluster" << endl;
				assert (prev_clusters[0] == prev_clusters[1]);

				size_t Ka = cluster_matrix.getClusterCount();

				bool accept = split(cluster_matrix, data_ids, prev_clusters);

				if (accept) {
					size_t Kb = cluster_matrix.getClusterCount();
					assert (Kb == Ka + 1);

					_statistics.split_cluster_events_accept++;
				} else {
					_statistics.split_cluster_events_reject++;
				}
				break;
			}
		case 2: 
			{ // merge step
				fout << "Merge these clusters" << endl;
				assert (prev_clusters[0] != prev_clusters[1]);

				size_t Ka = cluster_matrix.getClusterCount();
				
				bool accept = merge(cluster_matrix, data_ids, prev_clusters);

				if (accept) {
					size_t Kb = cluster_matrix.getClusterCount();
					assert (Kb == Ka - 1);
					
					_statistics.merge_cluster_events_accept++;
				} else {
					_statistics.merge_cluster_events_reject++;
				}
				break;
			}
	} 

	// check
	for (auto data_id: data_ids) {
		fout << "Check data id " << data_id << endl;
		assert(cluster_matrix.assigned(data_id));
	}
}

void JainNealAlgorithm::printStatistics() {
	int verbosity = _verbosity;
	_verbosity = Debug;

	fout << endl;
	fout << "Statistics:" << endl;
	fout << " # of split attempts: " << _statistics.split_attempts << endl;
	fout << "   o of accepted split cluster events: " << _statistics.split_cluster_events_accept << endl;
	fout << "     - of accept 0->? overrides: " << _statistics.split_source_likelihood_zero << endl;
//	fout << "     - of likelihood ratio 0/0 (accepted anyhow): " << _statistics.split_likelihood_both_zero << endl;
	fout << "     - of events accepted after consideration: " << _statistics.split_likelihood_both_nonzero_accept << endl;
	
	fout << "   o of rejected split cluster events: " << _statistics.split_cluster_events_reject << endl;
	fout << "   # of reject 0->0 overrides: " << _statistics.split_likelihood_both_zero << endl;
	fout << "   # of reject ?->0 overrides: " << _statistics.split_target_likelihood_zero << endl;

	fout << "   # of metropolis-hastings decisions: " << _statistics.split_likelihood_both_nonzero << endl;
	fout << "     - accepted: " << _statistics.split_likelihood_both_nonzero_accept << endl;
	fout << "     - rejected: " << _statistics.split_likelihood_both_nonzero_reject << endl;

	fout << " # of merge attempts: " << _statistics.merge_attempts << endl;
	fout << "   o of accepted merge cluster events: " << _statistics.merge_cluster_events_accept << endl;
	fout << "     - of accept 0->? overrides: " << _statistics.merge_source_likelihood_zero << endl;
//	fout << "     - of likelihood ratio 0/0 (accepted anyhow): " << _statistics.merge_likelihood_both_zero << endl;
	fout << "     - of events accepted after consideration: " << _statistics.merge_likelihood_both_nonzero_accept << endl;
	
	fout << "   o of rejected merge cluster events: " << _statistics.merge_cluster_events_reject << endl;
	fout << "   # of reject 0->0 overrides: " << _statistics.merge_likelihood_both_zero << endl;
	fout << "   # of reject ?->0 overrides: " << _statistics.merge_target_likelihood_zero << endl;

	fout << "   # of metropolis-hastings decisions: " << _statistics.merge_likelihood_both_nonzero << endl;
	fout << "     - accepted: " << _statistics.merge_likelihood_both_nonzero_accept << endl;
	fout << "     - rejected: " << _statistics.merge_likelihood_both_nonzero_reject << endl;


	// set back to zero
	_statistics = (statistics_t){0};

	_verbosity = verbosity;
} 
