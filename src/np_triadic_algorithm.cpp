#include <np_triadic_algorithm.h>

#include <iostream>
#include <iomanip>
#include <pretty_print.hpp>
#include <dim1algebra.hpp>
#include <cmath>

using namespace std;

// note, you'll need g++-7 and C++17 to obtain special_math function beta.

TriadicAlgorithm::TriadicAlgorithm(
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
	
	_statistics = (statistics_t){{0}};

	fout << "likelihood: " << _likelihood.getSuffies() << endl;

	_beta = 0.5;

	_split_method = sams_prior;
	_split_method = simple_random_split;
}

double TriadicAlgorithm::ratioStateProb(bool split, const std::vector<int> & more, const std::vector<int> & less) {
	double logfraction = 0.0;
	for (int i = 0; i < (int)more.size(); ++i) {
		logfraction += lgamma(more[i]);
	}
	for (int i = 0; i < (int)less.size(); ++i) {
		logfraction -= lgamma(less[i]);
	}
	double result = _alpha * std::exp(logfraction);
	if (split)
		return result;
	else
		return 1.0 / result;
}

/**
 * @param[in] N                                  number of data points assigned to clustesr
 * @param[in] C                                  maximum number of clusters (in merge/split between 2 to 3, it is 3)
 */
double TriadicAlgorithm::ratioProposal(bool split, int N, int C) {
	if (_split_method == sams_prior) {
		return 1.0;
	}
	double f = C / (C - 1.0);
	//double result = pow(3, nc - 3) * pow(2, 2 - nc);
	// we can generalize to  (c/(c-1))^(n-c+1)/c
	double result = pow(f, N - C + 1) / C;
	if (split) 
		return result;
	else
		return 1.0 / result;
}

// always merge first Q
void TriadicAlgorithm::propose_merge(std::vector<data_ids_t> &pdata, const data_ids_t &data_picks,
		cluster_ids_t &cluster_ids, split_method_t split_method) {
	
	//int D = data_picks.size();
	int C = cluster_ids.size();
	int Q = C - 1;
	std::vector<cluster_t*> cluster(Q);

	// assign Q of the D uniformly sampled data points to separate clusters
	for (int i = 0; i < Q; ++i) {
		pdata[i].push_back(data_picks[i]);
	}
	
	// get assignments (just ids, not the data itself)
	data_ids_t data_ids;
	for (int i = 0; i < C; ++i) {
		cluster[i] = _cluster_matrix->getCluster(cluster_ids[i]);
		_cluster_matrix->getAssignments(cluster_ids[i], data_ids);
	}
	fout << "The number of data points assigned to the cluster now is " << data_ids.size() << endl;
	
	// shuffle the data_id items
	algebra::random_order(data_ids.begin(), data_ids.end());
	
	std::vector<double> ldata(Q, 1.0); 

	// assign all data points
	for (int i = 0; i < (int)data_ids.size(); ++i) {
		int data_id = data_ids[i];
		for (int i = 0; i < Q; ++i) {
			if (data_picks[i] == data_id) 
				continue;
		}

		switch (split_method) {
			case random_mixing:
				{
					assert(false);
					break;
				}
			case simple_random_split:  
				{
					int index = algebra::random_weighted_pick(ldata.begin(), ldata.end(), _generator);
					pdata[index].push_back(data_id);
					break;
				} 
			case sams_prior: 
				{
					data_t * datum = _cluster_matrix->getDatum(data_id);
					for (int i = 0; i < Q; i++) {
						_likelihood.init(cluster[i]->getSuffies());
						ldata[i] = _likelihood.probability(*datum) * pdata[i].size();
					}
					int index = algebra::random_weighted_pick(ldata.begin(), ldata.end(), _generator);
					pdata[index].push_back(data_id);
					break;
				}
			case sams_random_walk: 
				{
					assert(false);
					break;
				}
		}
	}
}

void TriadicAlgorithm::propose_split(std::vector<data_ids_t> &pdata, const data_ids_t &data_picks, 
		cluster_ids_t &cluster_ids, cluster_t *new_cluster, split_method_t split_method) {

	int D = data_picks.size();
	int C = cluster_ids.size();
	int Q = C + 1;
	std::vector<cluster_t*> cluster(Q);
	data_ids_t data_ids;

	// assign uniformly sampled data points to separate clusters
	for (int i = 0; i < D; ++i) {
		pdata[i].push_back(data_picks[i]);
	}
	
	// current cluster
	for (int i = 0; i < C; ++i) {
		cluster[i] = _cluster_matrix->getCluster(cluster_ids[i]);
	}
	cluster[C] = new_cluster;

	for (int i = 0; i < C; ++i) {
		_cluster_matrix->getAssignments(cluster_ids[i], data_ids);
	}
	fout << "The number of data points assigned to the cluster now is " << data_ids.size() << endl;
	
	// shuffle the data items
	algebra::random_order(data_ids.begin(), data_ids.end());

	std::vector<double> ldata(D, 1.0);

	// assign all data points
	for (int i = 0; i < (int)data_ids.size(); ++i) {
		int data_id = data_ids[i];
		for (int i = 0; i < D; ++i) {
			if (data_picks[i] == data_id) 
				continue;
		}

		switch (split_method) {
			case random_mixing:
				{
					assert(false);
					break;
				}
			case simple_random_split:  
				{
					int index = algebra::random_weighted_pick(ldata.begin(), ldata.end(), _generator);
					pdata[index].push_back(data_id);
					break;
				} 
			case sams_prior: 
				{
					data_t * datum = _cluster_matrix->getDatum(data_id);
					for (int i = 0; i < Q; i++) {
						_likelihood.init(cluster[i]->getSuffies());
						ldata[i] = _likelihood.probability(*datum) * pdata[i].size();
					}
					int index = algebra::random_weighted_pick(ldata.begin(), ldata.end(), _generator);
					pdata[index].push_back(data_id);
					break;
				}
			case sams_random_walk: 
				{
					assert(false);
					break;
				}
		}
	}
}

void TriadicAlgorithm::checkLikelihoods(double lsrc, double ldest, step_t statistics_step, bool &accept, bool &overwrite) {
	if (ldest == 0 && lsrc == 0) {
		fout << "Both likelihood of existing and new cluster is zero. Reject it." << endl;
		overwrite = true;
		accept = accept && false;
		statistics_step.likelihood_both_zero++;
	} else if (lsrc == 0) {
		fout << "Likelihood of source cluster is zero. Accept new cluster." << endl;
		overwrite = true;
		accept = accept && true;
		statistics_step.source_likelihood_zero++;
	} else if (ldest == 0) {
		fout << "Likelihood of target cluster is " << ldest << ". Reject it in all cases." << endl;
		fout << "Just for the record, likelihood of source is: " << lsrc << endl;
		overwrite = true;
		accept = accept && false;
		statistics_step.target_likelihood_zero++;
	} else {
		statistics_step.likelihood_both_nonzero++;
	}
}

bool TriadicAlgorithm::split(
		const data_ids_t & data_picks,
		cluster_ids_t & cluster_ids
		) {
	
	int C = cluster_ids.size();
	int Q = C + 1;
	assert (Q == (int)data_picks.size());

	// dimension of dataX is C, pdataX is Q, ldataX is set both to Q, cluster as well
	// data can be obtained as pointers, pdata has to be constructed by copying elements
	std::vector<data_ids_t> data_ids(C), pdata_ids(Q);
	std::vector<dataset_t*> data(C);
	std::vector<dataset_t> pdata(Q);
	std::vector<cluster_t*> cluster(Q);
	std::vector<int> ndata(C), npdata(Q);
	std::vector<double> ldata(Q), lpdata(Q);
	int N = 0;
	double rP, rQ, rL = 1.0;
	bool accept = true, overwrite = false;

	_statistics.split.attempts++;
	
	for (int i = 0; i < C; i++) {
		cluster[i] = _cluster_matrix->getCluster(cluster_ids[i]);
		_cluster_matrix->getAssignments(cluster_ids[i], data_ids[i]);
	}
	
	// new cluster
	Suffies *suffies = _nonparametrics.sample_base(_generator);
	cluster[C] = new cluster_t(*suffies);
	
	for (int i = 0; i < C; ++i) {
		ndata[i] = data_ids[i].size();
		N += ndata[i];
	}

	// split also concerns moves from 1 to 2 and 2 to 1, not only from 1 to 3 and 2 to 3.
	propose_split(pdata_ids, data_picks, cluster_ids, cluster[C], _split_method);

	for (int i = 0; i < Q; ++i) {
		npdata[i] = pdata_ids[i].size();
	}

	fout << "Clusters sizes: from " << ndata << " to " << npdata << endl;

	rP = ratioStateProb(true, npdata, ndata);
	rQ = ratioProposal(true, N, C + 1);

	fout << "The q(.|.)p(.) ratio for the MH split is " << rQ << " * " << rP << " = " << rQ * rP << endl;
	
	// before split
	for (int i = 0; i < C; ++i) {
		cluster[i] = _cluster_matrix->getCluster(cluster_ids[i]);
		_likelihood.init(cluster[i]->getSuffies());
		data[i] = _cluster_matrix->getData(cluster_ids[i]);
		ldata[i] = _likelihood.probability(*data[i]);
	}
	ldata[C] = 1.0;

	// after split
	for (int i = 0; i < C; ++i) {
		cluster[i] = _cluster_matrix->getCluster(cluster_ids[i]);
		_likelihood.init(cluster[i]->getSuffies());
		_cluster_matrix->getData(pdata_ids[i], pdata[i]);
		lpdata[i] = _likelihood.probability(pdata[i]);
	}
	_likelihood.init(cluster[C]->getSuffies());
	_cluster_matrix->getData(pdata_ids[C], pdata[C]);
	lpdata[C] = _likelihood.probability(pdata[C]);

	for (int i = 0; i < C; ++i) {
		checkLikelihoods(ldata[i], lpdata[i], _statistics.split, accept, overwrite);
	}

	// normal MH step
	if (!overwrite) {
		for (int i = 0; i < C; ++i) {
			rL *= lpdata[i] / ldata[i];
		}
		double a_split = rQ * rP * rL;
		fout << "Acceptance for split: " << a_split << endl;

		double u = _distribution(_generator);

		if (a_split < u) {
			_statistics.split.likelihood_both_nonzero_reject++;
			accept = false;
		} else {
			_statistics.split.likelihood_both_nonzero_accept++;
			accept = true;
		}
	}

	if (accept) {
		fout << "Accept split!" << endl;
		cluster_id_t new_cluster_id = _cluster_matrix->addCluster(cluster[C]);
		cluster_ids.push_back(new_cluster_id);
		// assign or reassign (even to same cluster data points)
		for (int i = 0; i < Q; ++i) {
			for (auto data_id: pdata_ids[i]) {
				fout << "Retract and reassign data " << data_id << " to cluster " << cluster_ids[i] << endl;
				_cluster_matrix->retract(data_id);
				_cluster_matrix->assign(cluster_ids[i], data_id);
			}
		}
	} else {
		fout << "Reject split!" << endl;
		delete cluster[C];
	}
	return accept;
}

bool TriadicAlgorithm::merge(
		const data_ids_t & data_picks,
		cluster_ids_t & cluster_ids
		) {

	int C = cluster_ids.size();
	int D = data_picks.size();

	assert(C==D);

	std::vector<data_ids_t> data_ids(C), pdata_ids(C-1);
	std::vector<dataset_t*> data(C);
	std::vector<dataset_t> pdata(C-1);
	std::vector<cluster_t*> cluster(C);
	std::vector<int> ndata(C), npdata(C-1);
	std::vector<double> ldata(C), lpdata(C);
	int N = 0;
	double rP, rQ, rL = 1.0;
	bool accept = true, overwrite = false;

	_statistics.merge.attempts++;

	for (int i = 0; i < C; ++i) {
		cluster[i] = _cluster_matrix->getCluster(cluster_ids[i]);
		_cluster_matrix->getAssignments(cluster_ids[i], data_ids[i]);
	}
	
	for (int i = 0; i < C; ++i) {
		ndata[i] = data_ids[i].size();
		N += ndata[i];
	}

	propose_merge(pdata_ids, data_picks, cluster_ids, _split_method);
	
	for (int i = 0; i < (C - 1); ++i) {
		npdata[i] = pdata_ids[i].size();
	}
	
	fout << "Clusters sizes: from " << ndata << " to " << npdata << endl;

	rP = ratioStateProb(false, npdata, ndata);
	rQ = ratioProposal(false, N, C);

	fout << "The q(.|.)p(.) ratio for the MH merge is " << rQ << " * " << rP << " = " << rQ * rP << endl;

	// before merge
	for (int i = 0; i < C; ++i) {
		data[i] = _cluster_matrix->getData(cluster_ids[i]);
		cluster[i] = _cluster_matrix->getCluster(cluster_ids[i]);
		_likelihood.init(cluster[i]->getSuffies());
		ldata[i] = _likelihood.probability(*data[i]);
	}
	
	// after merge
	for (int i = 0; i < C - 1; ++i) {
		cluster[i] = _cluster_matrix->getCluster(cluster_ids[i]);
		_likelihood.init(cluster[i]->getSuffies());
		_cluster_matrix->getData(pdata_ids[i], pdata[i]);
		lpdata[i] = _likelihood.probability(pdata[i]);
	}
	lpdata[C-1] = 1.0;

	for (int i = 0; i < C; ++i) {
		checkLikelihoods(ldata[i], lpdata[i], _statistics.merge, accept, overwrite);
	}

	if (!overwrite) {
		for (int i = 0; i < C; ++i) {
			rL *= lpdata[i] / ldata[i];
		}

		double a_merge = rQ * rP * rL;
		fout << "Acceptance for merge: " << a_merge << endl;

		double u = _distribution(_generator);
		fout << "Sampled u: " << u << endl;

		if (a_merge < u) {
			_statistics.merge.likelihood_both_nonzero_reject++;
			accept = false;
		} else {
			_statistics.merge.likelihood_both_nonzero_accept++;
			accept = true;
		}
	}

	if (accept) {
		fout << "Accept merge!" << endl;
		// actually perform the move: all items from i go to j
		for (int i = 0; i < C; ++i) {
			for (auto data_id: pdata_ids[i]) {
				fout << "Retract and reassign data " << data_id << " to cluster " << cluster_ids[i] << endl;
				_cluster_matrix->retract(data_id);
				_cluster_matrix->assign(cluster_ids[i], data_id);
			}
		}
		// TODO: actually delete now empty cluster
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
void TriadicAlgorithm::update(
			membertrix & cluster_matrix,
			const data_ids_t & cdata_picks
		) {
	// copy data points (because we will need to be removing points for the dyadic sampler steps)
	data_ids_t data_picks = cdata_picks;

	// there should be exactly three data points sampled
	assert (data_picks.size() == 3);

	// the data ids should be unique
	assert (data_picks[0] != data_picks[1]);
	assert (data_picks[0] != data_picks[2]);

	_cluster_matrix = &cluster_matrix;

	static int step = 0;
	fout << "Update step " << ++step << " in Triadic split-merge algorithm" << endl;

	// store cluster indices of current observations
	cluster_ids_t cluster_ids; 
	cluster_ids.clear();
	for (auto data_id: data_picks) {
		cluster_id_t cluster_id = _cluster_matrix->getClusterId(data_id);
		fout << "We found data item " << data_id << " at cluster " << cluster_id << endl;
		cluster_ids.push_back( cluster_id );
	}

	int uniq_cluster_count = algebra::count_unique(cluster_ids.begin(), cluster_ids.end());
	fout << "The data comes from " << uniq_cluster_count << " unique cluster(s)" << endl;

	bool jain_neal_split = false;
	double u = _distribution(_generator);
	if (uniq_cluster_count == 1) {
		fout << "We are gonna split this one cluster" << endl;

		// remove one of the data points and corresponding cluster
		data_picks.pop_back();
		cluster_ids.pop_back();

		assert(data_picks.size() == 2);
		assert(cluster_ids.size() == 2);

		size_t Ka = _cluster_matrix->getClusterCount();

		bool accept = split(data_picks, cluster_ids);

		if (accept) {
			size_t Kb = _cluster_matrix->getClusterCount();
			assert (Kb == Ka + 1);

			_statistics.split.cluster_events_accept++;
		} else {
			_statistics.split.cluster_events_reject++;
		}
		jain_neal_split = true;
	} else if (u < _beta) {
		fout << "Merge these two clusters" << endl;

		// remove one of the data points that points to the same cluster as one of the other data points
		int index = algebra::duplicate_pick(cluster_ids.begin(), cluster_ids.end(), _generator);
		data_picks.erase(data_picks.begin() + index);
		cluster_ids.erase(cluster_ids.begin() + index);

		size_t Ka = _cluster_matrix->getClusterCount();

		bool accept = merge(data_picks, cluster_ids);

		if (accept) {
			size_t Kb = _cluster_matrix->getClusterCount();
			assert (Kb == Ka - 1);

			_statistics.merge.cluster_events_accept++;
		} else {
			_statistics.merge.cluster_events_reject++;
		}

		jain_neal_split = true;
	}

	if (jain_neal_split) {
		// a single cluster
		switch (uniq_cluster_count) {
			case 0: default:
				assert(false);
				break;
			case 1:
				assert(false);
				break;
			case 2: 
				{
					fout << "We are gonna split these two clusters" << endl;

					size_t Ka = _cluster_matrix->getClusterCount();

					bool accept = split(data_picks, cluster_ids);

					if (accept) {
						size_t Kb = _cluster_matrix->getClusterCount();
						assert (Kb == Ka + 1);

						_statistics.split.cluster_events_accept++;
					} else {
						_statistics.split.cluster_events_reject++;
					}
					break;
				}
			case 3: 
				{ // merge step
					fout << "Merge these clusters" << endl;

					size_t Ka = _cluster_matrix->getClusterCount();

					bool accept = merge(data_picks, cluster_ids);

					if (accept) {
						size_t Kb = _cluster_matrix->getClusterCount();
						assert (Kb == Ka - 1);

						_statistics.merge.cluster_events_accept++;
					} else {
						_statistics.merge.cluster_events_reject++;
					}
					break;
				}
		} 
	}
	// check
	for (auto data_id: data_picks) {
		fout << "Check data id " << data_id << endl;
		assert(_cluster_matrix->assigned(data_id));
	}
}

void TriadicAlgorithm::printStatistics() {
	int verbosity = _verbosity;
	_verbosity = Debug;

	fout << endl;
	fout << "Statistics:" << endl;
	fout << " # of split attempts: " << _statistics.split.attempts << endl;
	fout << "   o of accepted split cluster events: " << _statistics.split.cluster_events_accept << endl;
	fout << "     - of accept 0->? overrides: " << _statistics.split.source_likelihood_zero << endl;
//	fout << "     - of likelihood ratio 0/0 (accepted anyhow): " << _statistics.split.likelihood_both_zero << endl;
	fout << "     - of events accepted after consideration: " << _statistics.split.likelihood_both_nonzero_accept << endl;
	
	fout << "   o of rejected split cluster events: " << _statistics.split.cluster_events_reject << endl;
	fout << "   # of reject 0->0 overrides: " << _statistics.split.likelihood_both_zero << endl;
	fout << "   # of reject ?->0 overrides: " << _statistics.split.target_likelihood_zero << endl;

	fout << "   # of metropolis-hastings decisions: " << _statistics.split.likelihood_both_nonzero << endl;
	fout << "     - accepted: " << _statistics.split.likelihood_both_nonzero_accept << endl;
	fout << "     - rejected: " << _statistics.split.likelihood_both_nonzero_reject << endl;

	fout << " # of merge attempts: " << _statistics.merge.attempts << endl;
	fout << "   o of accepted merge cluster events: " << _statistics.merge.cluster_events_accept << endl;
	fout << "     - of accept 0->? overrides: " << _statistics.merge.source_likelihood_zero << endl;
//	fout << "     - of likelihood ratio 0/0 (accepted anyhow): " << _statistics.merge.likelihood_both_zero << endl;
	fout << "     - of events accepted after consideration: " << _statistics.merge.likelihood_both_nonzero_accept << endl;
	
	fout << "   o of rejected merge cluster events: " << _statistics.merge.cluster_events_reject << endl;
	fout << "   # of reject 0->0 overrides: " << _statistics.merge.likelihood_both_zero << endl;
	fout << "   # of reject ?->0 overrides: " << _statistics.merge.target_likelihood_zero << endl;

	fout << "   # of metropolis-hastings decisions: " << _statistics.merge.likelihood_both_nonzero << endl;
	fout << "     - accepted: " << _statistics.merge.likelihood_both_nonzero_accept << endl;
	fout << "     - rejected: " << _statistics.merge.likelihood_both_nonzero_reject << endl;


	// set back to zero
	_statistics = (statistics_t){0};

	_verbosity = verbosity;
} 
