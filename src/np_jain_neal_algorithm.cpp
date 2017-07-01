#include <np_jain_neal_algorithm.h>

#include <iostream>
#include <iomanip>
#include <pretty_print.hpp>
#include <dim1algebra.hpp>
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <cmath>

using namespace std;

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

//	_verbosity = Debug;
	_verbosity = Warning;

	fout << "likelihood: " << _likelihood.getSuffies() << endl;
		
	_statistics.new_clusters_events = 0;
}

int factorial(int n) {
	int ret = 1;
	for(int i = 1; i <= n; ++i)
		ret *= i;
	return ret;
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
	fout << "Update step " << ++step << " in Jain & Neal's split-merge algorithm for several data items at once" << endl;

	// store cluster indices and retract only current observations
	std::vector<cluster_id_t> prev_clusters; //data_ids.size());
	prev_clusters.clear();
	for (auto data_id: data_ids) {
		//fout << "Retract data item " << data_id << endl;
		//assert (cluster_matrix.assigned(data_id) );
		cluster_id_t cluster_id = cluster_matrix.getClusterId(data_id);
		fout << "Found data " << data_id << " at cluster " << cluster_id << endl;
		prev_clusters.push_back( cluster_id );
		// Todo: check c_j, should it be retracted?
		//assert (cluster_matrix.assigned(data_id) );
		//fout << "Retract data " << data_id << " from cluster " << cluster_id << endl;
		//cluster_matrix.retract(cluster_id, data_id);
		//fout << "Retracted" << endl;
	}
	int uniq_cluster_count = algebra::count_unique(prev_clusters.begin(), prev_clusters.end());
	fout << "There are " << uniq_cluster_count << " unique clusters" << endl;
	// a single cluster
	switch (uniq_cluster_count) {
		case 0: default:
			assert(false);
		break;
		case 1: { // split step
			fout << "Split this cluster" << endl;
			assert (prev_clusters[0] == prev_clusters[1]);

			size_t Ka = cluster_matrix.getClusterCount();

			// get assignments (just ids, not the data itself)
			data_ids_t *data_ids1 = cluster_matrix.getAssignments(prev_clusters[0]);

			data_ids_t data_ids2_0 = { data_ids[0] };
			data_ids_t data_ids2_1 = { data_ids[1] };

			for (auto data_id: *data_ids1) {
				if (data_id == data_ids[0]) continue;
				if (data_id == data_ids[1]) continue;
				
				// sample uniform random variable to randomly split over clusters
				double toss = _distribution(_generator);
				if (toss <= 0.5) {
					// move these points
					data_ids2_0.push_back(data_id);
				} else {
					// keep these points
					data_ids2_1.push_back(data_id);
				}
			}

			int nc0 = data_ids2_0.size();
			assert (nc0 > 0);
			int nc1 = data_ids2_1.size();
			assert (nc1 > 0);
			int nc = nc0 + nc1; // TODO: check again size
			fout << "Cluster is size " << nc << " and after split becomes " << nc0 << " and " << nc1 << endl;

			// note, you'll need g++-7 and C++17 to obtain special_math function beta.
			double p21 = _alpha * beta(nc0-1, nc1-1);
			double q12 = pow(2, nc - 2);
			
			// only consider the data points that move to the new one, but at the old cluster
			auto cluster0 = cluster_matrix.getCluster(prev_clusters[0]);
			_likelihood.init(cluster0->getSuffies());
			dataset_t *data0 = cluster_matrix.getData(data_ids2_0);
			double l2_0 = _likelihood.probability(*data0);

			// create new cluster, and consider same dataset there
			Suffies *suffies = _nonparametrics.sample_base(_generator);
			cluster_t *new_cluster = new cluster_t(*suffies);
			_likelihood.init(new_cluster->getSuffies());
			double l1 = _likelihood.probability(*data0);

			double l21; 
	
			if (l2_0 == 0 && l1 == 0) {
				// both ways zero, so just random walk
				l21 = 1;
				//std::cout << "|";
			} else if (l1 == 0) {
				// likelihood is zero where we are going
				l21 = 0;
				//std::cout << "|";
			} else {
				// if we do not have a stream of . then there is nothing to accept!
				l21 = l2_0 / l1; 
				//std::cout << ".";
			}

			double a_split = q12 * p21 * l21; // or actually min(1, ...) but we compare with u ~ U(0,1) anyway
			
			// sample uniform random variable
			double u = _distribution(_generator);

			bool accept;
			if (a_split < u) {
				fout << "Reject!" << endl;
				accept = false;
			} else {
				fout << "Accept!" << endl;
				accept = true;
			}

			if (accept) {
				// actually perform the move: only points that move have to be moved...
				cluster_id_t new_cluster_id = cluster_matrix.addCluster(new_cluster);
				fout << "Add cluster with id " << new_cluster_id << endl;				
				for (auto data: data_ids2_0) {
					fout << "Add data " << data << " to the new cluster " << new_cluster_id << endl;
					cluster_matrix.retract(data);
					cluster_matrix.assign(new_cluster_id, data);
				}
				// points in data_ids2_1 need not be moved

				size_t Kb = cluster_matrix.getClusterCount();
				assert (Kb == Ka + 1);
			}
			// TODO: anything that needs to be re-established on rejection?
		}
		break;
		case 2: { // merge step
			fout << "Merge these clusters" << endl;
			assert (prev_clusters[0] != prev_clusters[1]);

			size_t Ka = cluster_matrix.getClusterCount();

			// get assignments (just ids, not the data itself)
			data_ids_t *data_ids0 = cluster_matrix.getAssignments(prev_clusters[0]);
			data_ids_t *data_ids1 = cluster_matrix.getAssignments(prev_clusters[1]);

			int nc0 = data_ids0->size();
			int nc1 = data_ids1->size();
			int nc = nc0 + nc1;
			fout << "Clusters are size " << nc0 << " and " << nc1 << " and after merge become " << nc << endl;

			// check if type is correct (no rounding off to ints by accident)
			double p12 = 1.0/(_alpha * beta(nc0-1, nc1-1));
			fout << "Fraction P(1)/P(2) becomes: " << p12 << endl;
			double q21 = pow(0.5, nc - 2);
			fout << "Fraction q(1|2)/q(2|1) becomes: " << q21 << endl;

			// calculate likelihoods
			// if we would have stored the likelihoods somewhere of the original clusters we might have reused them
			// now we have to calculate them from scratch
			
			// likelihood of cluster 0 before merge
			auto cluster0 = cluster_matrix.getCluster(prev_clusters[0]);
			_likelihood.init(cluster0->getSuffies());
			dataset_t * data0 = cluster_matrix.getData(prev_clusters[0]);
			double l2_0 = _likelihood.probability(*data0);
			fout << "Likelihood of data in cluster 0 before merge " << l2_0 << endl;

#define CLUSTER_SUFFICIENT_STATISTICS_UNCHANGED_AND_LIKELIHOOD_DATA_INDEPENDENT
#ifdef CLUSTER_SUFFICIENT_STATISTICS_UNCHANGED_AND_LIKELIHOOD_DATA_INDEPENDENT
			// 
			auto cluster1 = cluster_matrix.getCluster(prev_clusters[1]);
			_likelihood.init(cluster1->getSuffies());
			double l1_0 = _likelihood.probability(*data0);
			
			fout << "Likelihood of same data in cluster 1 after merge " << l1_0 << endl;
			double l12;
			if (l1_0 == 0 && l2_0 == 0) {
				// likelihood is zero before and after merge... we should just walk...
				l12 = 1;
				//std::cout << "]";
			} else if (l2_0 == 0) {
				// just reject
				l12 = 0;
				//std::cout << "-";
			} else {
//				assert (l2_0 != 0);
				l12 = l1_0 / l2_0;
				//std::cout << "o";
			}
#else
			// likelihood of cluster 1 before merge
			auto cluster1 = cluster_matrix.getCluster(prev_clusters[1]);
			_likelihood.init(cluster1->getSuffies());
			dataset_t * data1 = cluster_matrix.getData(prev_clusters[1]);
			double l2_1 = _likelihood.probability(data1);

			// likelihood of cluster 1 [c_j] after merge (cluster 0 [c_i] is removed)
			// likelihood is still initialized to sufficient statistics of second cluster
			double l1_0 = _likelihood.probability(data0);
			double l1_1 = l2_1;
			double l1 = l1_0 * l1_1;
			double l12 = l1/(l2_0*l2_1);
#endif
			double a_merge = q21 * p12 * l12; // or actually min(1, ...) but we compare with u ~ U(0,1) anyway

			// sample uniform random variable
			double u = _distribution(_generator);

			bool accept;
			if (a_merge < u) {
				fout << "Reject!" << endl;
				accept = false;
			} else {
				fout << "Accept!" << endl;
				accept = true;
			}

			if (accept) {
				// actually perform the move: all items from i go to j
				for (auto data: *data_ids0) {
					cluster_matrix.retract(data);
					cluster_matrix.assign(prev_clusters[1], data);
				}

				fout << "Check cluster count again" << endl;
				size_t Kb = cluster_matrix.getClusterCount();
				assert (Kb == Ka - 1);
			}
			// TODO: anything that needs to be re-established on rejection?
		}
		break;
	} 

	// check
	for (auto data_id: data_ids) {
		fout << "Check data id " << data_id << endl;
		assert(cluster_matrix.assigned(data_id));
	}
}

void JainNealAlgorithm::printStatistics() {
	int verbosity = _verbosity;
	_verbosity = 0;
	fout << "Statistics:" << endl;
	fout << " # of new cluster events: " << _statistics.new_clusters_events << endl;
	_verbosity = verbosity;
} 
