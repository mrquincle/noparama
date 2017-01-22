#include <Eigen/Dense>
#include <random>

/*!
 * The base class t_distribution allows using this class, there were needs to be sampled from a probability 
 * distribution. The sufficient statistics are distribution dependent and are stored within this derived class.
 */
class normalinvwishart_distribution: t_distribution {
	private:
		// assumed to be in t_distribution: SufficientStatistics _statistics;
		NIWStatistics _normalinvwishart_statistics;

		// prepare statistics to be used
		NormalStatistics _normal_statistics; 
		IWStatistics _invwishart_statistics;
		
		// result
		NormalStatistics _result_statistics; 
	public:
		normalinvwishart_distribution(NIWStatistics & niw_statistics): _normalinvwishart_statistics(niw_statistics)
		{
			_invwishart_statistics.nu = niw_statistics.nu;
			_invwishart_statistics.lambda = niw_statistics.lambda;
			
			_normal_statistics.mu = niw_statistics.mu;

		}

		template<typename _UniformRandomNumberGenerator >
		NormalStatistics & operator()(_UniformRandomNumberGenerator & generator) const
		{
			// first sample covariance matrix from IW
			invwishart_distribution invwishart_dist(_invwishart_statistics);
			_result_statistics.sigma = invwishart_dist(generator);

			// next sample mu from normal distribution
			_normal_statistics.sigma = _result_statistics.sigma / _normalinvwishart_statistics.kappa;
			multivariate_normal_distribution multivariate_normal_dist(_normal_statistics);
			_result_statistics.mu = multivariate_normal_dis(generator);
			
			return _result_statistics;
		}
};
