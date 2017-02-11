#include <iostream>
#include <fstream>
#include <random>

#include <np_suffies.h>
#include <membertrix>
#include <pretty_print.hpp>

using namespace std;
using namespace Eigen;

/**
 * Test the multivariate normal distribution. The test uses a predefined mean and covariance matrix and checks with
 * particular data items if the calculated probability is the same as the one calculated beforehand using octave.
 */
int main() {
	std::cout << "Start test membertrix" << std::endl;
	std::random_device r;
	default_random_engine generator(r());

	membertrix trix;

	int N = 5;
	int K = 4;

	std::cout << "Dataset" << std::endl;
	dataset_t dataset(N);
	for (int i = 0; i < N; ++i) {
		dataset[i] = new data_t({1, 1});
		
		// add to membertrix
		trix.addData(*dataset[i]);
	}

	std::cout << "Clusters" << std::endl;
	clusters_t clusters(K);
	std::vector<cluster_id_t> cluster_ids(K);
	for (int k = 0; k < K; ++k) {
		Vector2d mean(1,1);
		Matrix2d covar;
		covar << 2,0,1,2;
		Suffies_MultivariateNormal *suffies = new Suffies_MultivariateNormal(mean.size());
		suffies->mu = mean;
		suffies->sigma = covar;
		clusters[k] = new cluster_t(*suffies);
	
		// add to membertrix
		cluster_ids[k] = trix.addCluster(clusters[k]);

	}

	std::cout << "Assignments" << std::endl;
	std::vector<std::vector<data_id_t> > assignments(K);
	std::uniform_int_distribution<> dist(0, K-1);
	for (int i = 0; i < N; ++i) {
		int k = dist(generator);
		std::cout << "For data " << i << " sample cluster " << k << std::endl;
		assignments[k].push_back(i);
		trix.assign(k, i);
	}
	std::cout << "Print assignments" << std::endl;
	for (int k = 0; k < K; ++k) {
		std::cout << assignments[k] << std::endl;
	}

	std::cout << trix << std::endl;
	
	trix.cleanup();

	clusters_t check_clusters = trix.getClusters();

	std::cout << "Number of clusters: " << check_clusters.size() << std::endl;
//	assert ((int)check_clusters.size() == K);

	int removeClusters = 0;

	int l = dist(generator);
	std::vector<data_id_t> &members = assignments[l];
	std::cout << "Retract data items of cluster " << l << std::endl;
	for (auto member: members) {
		removeClusters = 1;
		std::cout << "Remove data item " << member << " from cluster " << l << std::endl;
		trix.retract(l, member);
	}
	
	std::cout << trix << std::endl;

	clusters_t check_clusters1 = trix.getClusters();

	std::cout << "Number of clusters: " << check_clusters1.size() << std::endl;
	assert ((int)check_clusters1.size() == ((int)check_clusters.size()-removeClusters));

}
