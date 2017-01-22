
class Likelihood {
	private:
		SufficientStatistics _sufficient_statistics;
	public:
		// Constructor
		Likelihood(SufficientStatistics & sufficient_statistics);

		// Get likelihood given data
		double get(t_data & observation);

		// Get likelihood given set of data
		std::vector<double> get(t_cluster & cluster);
};
