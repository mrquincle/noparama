#pragma once

#include <Eigen/Dense>
#include <random>

#include <statistics/invwishart.h>
#include <statistics/multivariatenormal.h>

/*!
 * The base class distribution_t allows using this class, there were needs to be sampled from a probability 
 * distribution. The sufficient statistics are distribution dependent and are stored within this derived class.
 */
class normal_inverse_wishart_distribution: public distribution_t {
	private:
		Suffies_NormalInvWishart _suffies_normalinvwishart;

		// prepare statistics to be used
		Suffies_MultivariateNormal _suffies_mvn; 
		Suffies_InvWishart _suffies_iw;
		
		// result
		Suffies_MultivariateNormal _suffies_result; 
	public:
		normal_inverse_wishart_distribution(Suffies_NormalInvWishart & suffies_niw): 
			_suffies_normalinvwishart(suffies_niw)
		{
			_suffies_iw.nu = suffies_niw.nu;
			_suffies_iw.Lambda = suffies_niw.Lambda;
		
			// _suffies_mvn.sigma will be set at operator()
			_suffies_mvn.mu = suffies_niw.mu;
		}

		template<typename _UniformRandomNumberGenerator >
		Suffies_MultivariateNormal & operator()(_UniformRandomNumberGenerator & generator) const
		{
			// first sample covariance matrix from IW
			inverse_wishart_distribution invwishart_dist(_suffies_iw);
			_suffies_result.sigma = invwishart_dist(generator);

			// next sample mu from normal distribution
			_suffies_mvn.sigma = _suffies_result.sigma / _suffies_normalinvwishart.kappa;
			multivariate_normal_distribution multivariate_normal_dist(_suffies_mvn);
			_suffies_result.mu = multivariate_normal_dis(generator);
			
			return _suffies_result;
		}
};
