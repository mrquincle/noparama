#include <np_jain_neal_algorithm.h>

#include <iostream>
#include <iomanip>
#include <pretty_print.hpp>

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

	static int step = 0;
	fout << "Update step " << ++step << " in Jain & Neal's split-merge algorithm for several data items at once" << endl;

	// store cluster indices and retract only current observations
	std::vector<cluster_id_t> prev_clusters(data_ids.size());
	for (auto index: data_ids) {
		cluster_id_t cluster_id = cluster_matrix.getClusterId(index);
		prev_clusters.push_back( cluster_id );
		// Todo: check c_j, should it be retracted?
		cluster_matrix.retract(cluster_id, index);
	}
	std::set<cluster_id_t> prev_uniq_clusters;
	for (auto cl: prev_clusters) {
		prev_uniq_clusters.insert(cl);
	}
	int uniq_cluster_count = prev_uniq_clusters.size();

	// a single cluster
	switch (uniq_cluster_count) {
		case 0: default:
			assert(false);
		break;
		case 1: { // split step
			assert (prev_clusters[0] == prev_clusters[1]);

			size_t Ka = cluster_matrix.getClusterCount();

			// get assignments (just ids, not the data itself)
			data_ids_t *data_ids1 = cluster_matrix.getAssignments(prev_clusters[0]);

			data_ids_t data_ids2_0 = { data_ids[0] };
			data_ids_t data_ids2_1 = { data_ids[1] };

			for (auto data_id: *data_ids1) {
				
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
			int nc1 = data_ids2_1.size();
			int nc = nc0 + nc1; // TODO: check again size

			// TODO: check if factorials do not become too big
			// TODO: check if we can have a lookup for the Beta function
			double p21 = _alpha * (factorial(nc0 - 1) * factorial(nc1 - 1)) / factorial(nc - 1);
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

			double l21 = l2_0 / l1; 

			double a_split = q12 * p21 * l21; // or actually min(1, ...) but we compare with u ~ U(0,1) anyway

			// sample uniform random variable
			double u = _distribution(_generator);

			bool accept;
			if (a_split < u) {
				accept = false;
			} else {
				accept = true;
			}

			if (accept) {
				// actually perform the move: only points that move have to be moved...
				cluster_id_t new_cluster_id = cluster_matrix.addCluster(new_cluster);
				for (auto data: data_ids2_0) {
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
			assert (prev_clusters[0] != prev_clusters[1]);

			size_t Ka = cluster_matrix.getClusterCount();

			// get assignments (just ids, not the data itself)
			data_ids_t *data_ids0 = cluster_matrix.getAssignments(prev_clusters[0]);
			data_ids_t *data_ids1 = cluster_matrix.getAssignments(prev_clusters[1]);

			int nc0 = data_ids0->size();
			int nc1 = data_ids1->size();
			int nc = nc0 + nc1;

			// check if type is correct (no rounding off to ints by accident)
			double p12 = 1.0/_alpha * factorial(nc - 1) / (  factorial(nc0 - 1) * factorial(nc1 - 1) );
			double q21 = pow(0.5, nc - 2);

			// calculate likelihoods
			// if we would have stored the likelihoods somewhere of the original clusters we might have reused them
			// now we have to calculate them from scratch
			
			// likelihood of cluster 0 before merge
			auto cluster0 = cluster_matrix.getCluster(prev_clusters[0]);
			_likelihood.init(cluster0->getSuffies());
			dataset_t * data0 = cluster_matrix.getData(prev_clusters[0]);
			double l2_0 = _likelihood.probability(*data0);

#define CLUSTER_SUFFICIENT_STATISTICS_UNCHANGED_AND_LIKELIHOOD_DATA_INDEPENDENT
#ifdef CLUSTER_SUFFICIENT_STATISTICS_UNCHANGED_AND_LIKELIHOOD_DATA_INDEPENDENT
			// 
			auto cluster1 = cluster_matrix.getCluster(prev_clusters[1]);
			_likelihood.init(cluster1->getSuffies());
			double l1_0 = _likelihood.probability(*data0);
			double l12 = l1_0 / l2_0;
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
				accept = false;
			} else {
				accept = true;
			}

			if (accept) {
				// actually perform the move: all items from i go to j
				for (auto data: *data_ids0) {
					cluster_matrix.retract(data);
					cluster_matrix.assign(prev_clusters[1], data);
				}

				size_t Kb = cluster_matrix.getClusterCount();
				assert (Kb == Ka - 1);
			}
			// TODO: anything that needs to be re-established on rejection?
		}
		break;
	} 

	// check
	for (auto data_id: data_ids) {
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
