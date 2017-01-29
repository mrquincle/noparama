#include <np_mcmc.h>
#include <assert.h>

MCMC::MCMC() {
}

void MCMC::run(ordered_data_t & ordered_data) {

	_dataset = ordered_data;

	int T = 10;
	
	for (int i = 0; i < N; ++i) {
		data_id_t j = _membertrix.addData(_dataset[i]);
		assert (i == j);
	}

	// initialize clusters
	for (int k = 0; t < K; ++k) {
		// sample sufficient statistics from nonparametrics
		sufficient_statistics <- nonparametrics;

		cluster_t *cluster = new cluster_t(sufficient_statistics);
		_membertrix.addCluster(cluster);
	}
	
	// update clusters
	for (int t = 0; t < T; ++t) {
		
		// iterator over all observations (observations get shifted around, but should not be counted twice)
		for (int i = 0; i < N; ++i) {
			// remove observation under consideration from cluster
			_membertrix.retract(i);

			// update cluster assignments, delete and create clusters
			UpdateClusterPopulation->update(membertrix.getClusters(), observations[i], nonparametrics, sample_pdf);

			// 

		}

		// update cluster parameters
		UpdateClusters->update(clusters, nonparametrics, prior, number_mh_steps);

	}


}
