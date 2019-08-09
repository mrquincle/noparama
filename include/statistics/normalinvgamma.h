#pragma once

#include <Eigen/Dense>
#include <random>

#include <statistics/distribution.h>
#include <statistics/gamma.h>
#include <statistics/multivariatenormal.h>

#include <pretty_print.hpp>

/*!
 * The base class distribution_t allows using this class, there were needs to be sampled from a probability 
 * distribution. The sufficient statistics are distribution dependent and are stored within this derived class.
 *
 * This is a multivariate form of the NIG distribution. We expect _suffies_nig to contain a Lambda matrix (not just
 * a gamma scalar).
 */
class normal_inverse_gamma_distribution: public distribution_t {
	private:
		Suffies_NormalInvGamma & _suffies_nig;

		// prepare statistics to be used
		Suffies_MultivariateNormal _suffies_mvn; 
		Suffies_Gamma _suffies_g;
		int _D;

		// result
		Suffies_ScalarNoise_MultivariateNormal * _suffies_result; 
	public:
		normal_inverse_gamma_distribution(Suffies_NormalInvGamma & suffies_nig): 
			_suffies_nig(suffies_nig), 
			_suffies_mvn(suffies_nig.D), 
			_suffies_result(NULL)
		{
			_distribution_type = NormalInvGamma;

			_D = _suffies_mvn.mu.size();

			_suffies_g.alpha = suffies_nig.alpha;
			_suffies_g.beta = suffies_nig.beta;
		
			_suffies_mvn.mu = suffies_nig.mu;
		}

		Suffies_ScalarNoise_MultivariateNormal* operator()(random_engine_t & generator) 
		{
			_suffies_result = new Suffies_ScalarNoise_MultivariateNormal(_D);

			// first sample variance from IG
			gamma_distribution gamma_dist(_suffies_g);
			Suffies_Double *sd = gamma_dist(generator);
			assert (sd->val != 0);
			double sigma_pow2 = 1.0/sd->val;
			_suffies_result->sigma = sqrt(sigma_pow2);

			// next sample mu from normal distribution
			_suffies_mvn.sigma = sigma_pow2 * _suffies_nig.Lambda.inverse();
			multivariate_normal_distribution multivariate_normal_dist(_suffies_mvn);
			Suffies_Unity_MultivariateNormal *suffies_mvn1 = multivariate_normal_dist(generator);

			_suffies_result->mu = suffies_mvn1->mu;
			
			return _suffies_result;
		}
};
