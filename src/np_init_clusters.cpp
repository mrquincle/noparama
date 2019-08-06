#include <np_init_clusters.h>

#include <pretty_print.hpp>

#include <statistics/dirichlet.h>

using namespace std;

InitClusters::InitClusters(
			random_engine_t & generator,
			distribution_t & nonparametrics
		): 
			_generator(generator),
			_nonparametrics(nonparametrics),
			_distribution(0.0, 1.0)
{
	_verbosity = 4;
}

/**
 * This implementation samples the sufficient statistics for each cluster (mean, variance) from a Dirichlet
 * distribution.
 */
void InitClusters::init(
			membertrix & cluster_matrix,
			int K
		) 
{
	// initialize clusters
	fout << "Add clusters to membership matrix" << endl;
	for (int k = 0; k < K; ++k) {
		// sample sufficient statistics from nonparametrics
		Suffies *suffies = ((dirichlet_process&)_nonparametrics).sample_base(_generator);
		
		cluster_t *cluster = new cluster_t(*suffies);
		int kTest = cluster_matrix.addCluster(cluster);

		assert (kTest == k);
	}

}
