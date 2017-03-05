#pragma once

#include <Eigen/Dense>
#include <random>

#include <statistics/invwishart.h>
#include <statistics/multivariatenormal.h>

#include <pretty_print.hpp>

/*!
 * The base class distribution_t allows using this class, there were needs to be sampled from a probability 
 * distribution. The sufficient statistics are distribution dependent and are stored within this derived class.
 */
class normal_inverse_wishart_distribution: public distribution_t {
	private:
		Suffies_NormalInvWishart & _suffies_normalinvwishart;

		// prepare statistics to be used
		Suffies_MultivariateNormal _suffies_mvn; 
		Suffies_InvWishart _suffies_iw;
		int _D;

		// result
		Suffies_MultivariateNormal * _suffies_result; 
	public:
		normal_inverse_wishart_distribution(Suffies_NormalInvWishart & suffies_niw): 
			_suffies_normalinvwishart(suffies_niw), 
			_suffies_mvn(suffies_niw.D), 
			_suffies_iw(suffies_niw.D),
			_suffies_result(NULL)
		{
			_distribution_type = NormalInvWishart;

			_D = _suffies_mvn.mu.size();

			_suffies_iw.nu = suffies_niw.nu;
			_suffies_iw.Lambda = suffies_niw.Lambda;
		
			// _suffies_mvn.sigma will be set at operator()
			_suffies_mvn.mu = suffies_niw.mu;
		}

		Suffies_MultivariateNormal* operator()(random_engine_t & generator) 
		{
			_suffies_result = new Suffies_MultivariateNormal(_D);

			// first sample covariance matrix from IW
			inverse_wishart_distribution invwishart_dist(_suffies_iw);
			Suffies_ZeroCentered_MultivariateNormal *suffies_mvn0 = invwishart_dist(generator);
			std::cout << "Sigma comes from: " << *suffies_mvn0 << std::endl;
			_suffies_result->sigma = suffies_mvn0->sigma;
//			delete suffies_mvn0;

			// next sample mu from normal distribution
			_suffies_mvn.sigma = _suffies_result->sigma / _suffies_normalinvwishart.kappa;
			std::cout << "Mu comes from: " << _suffies_mvn << std::endl;
			multivariate_normal_distribution multivariate_normal_dist(_suffies_mvn);
			Suffies_Unity_MultivariateNormal *suffies_mvn1 = multivariate_normal_dist(generator);

			_suffies_result->mu = suffies_mvn1->mu;
			
			return _suffies_result;
		}
};
