#include <assert.h>

#include <np_mcmc.h>
#include <np_update_clusters.h>
#include <np_update_cluster_population.h>

#include <Eigen/Dense>

using Eigen;

MCMC::MCMC(UpdateClusters & update_clusters, UpdateClusterPopulation & update_cluster_population):
	_update_clusters(update_cluster), _update_cluster_population(update_cluster_population)
{

}

void MCMC::run(dataset_t & dataset) {

	_dataset = dataset;

	int T = 10;
	int K = 30;

	int N = _dataset.size();
	
	for (int i = 0; i < N; ++i) {
		data_id_t j = _membertrix.addData(_dataset[i]);
		assert (i == j);
	}

	// initialize clusters
	for (int k = 0; k < K; ++k) {
		// sample sufficient statistics from nonparametrics
		//Suffies & suffies = sample_pdf(nonparametrics.base_distribution);

		Suffies & suffies = _nonparametrics.base_distribution(_generator);

		cluster_t *cluster = new cluster_t(suffies);
		_membertrix.addCluster(*cluster);
	}
	
	// update clusters
	for (int t = 0; t < T; ++t) {
		
		// iterator over all observations (observations get shifted around, but should not be counted twice)
		for (int i = 0; i < N; ++i) {
			// remove observation under consideration from cluster
			_membertrix.retract(i);

			// update cluster assignments, delete and create clusters
			update_cluster_population.update(membertrix.getClusters(), observations[i], nonparametrics, sample_pdf);

			// 

		}

		// update cluster parameters
		update_clusters.update(clusters, nonparametrics, prior, number_mh_steps);

	}


}
