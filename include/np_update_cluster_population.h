typedef MatrixXd t_cluster;
typedef std:set<t_cluster> t_cluster_population;

/**
 * This class UpdateClusterPopulation deletes, adds, and adjusts clusters. This in contrast with UpdateCluster which 
 * only adjusts the parameters assigned to a cluster and will leave the number of clusters invariant. 
 */
class UpdateClusterPopulation {
	public:
		UpdateClusterPopulation();

		void update(t_cluster_population & clusters);
};
