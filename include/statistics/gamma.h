#pragma once

#include <Eigen/Dense>

#include <random>

#include <statistics/distribution.h>

/*!
 * To sample from a gamma distribution, we need to initialize it with mean and variance
 * and then get samples out using a uniform random generator as input.
 */
class gamma_distribution: public distribution_t {
	private:
		double _alpha;
		double _beta;

		Suffies_Double *_suffies_result;
	public:
		/**
		 * Constructor for a gamma distribution.
		 *
		 * @param[in] suffies		mean and variance in one structure
		 */
		gamma_distribution(Suffies_Gamma & suffies)
			: _alpha(suffies.alpha), _beta(suffies.beta)
		{
		}

		/**
		 * Sample a random value according to the gamma distribution by calling the () operator.
		 *
		 * @param[in] generator		a uniform random number generator
		 * @return					a double, representing a single sample from the distribution
		 */
		Suffies_Double* operator()(random_engine_t & generator) 
		{
			_suffies_result = new Suffies_Double();

			static std::gamma_distribution<> dist(_alpha, _beta);
			
			_suffies_result->val = dist(generator);

			return _suffies_result;
		}
		
		/**
		 * Calculate the probability of a value given mean and variance.
		 *
		 * Mathematically, with mu, the mean and s the variance, this is:
		 *
		 *   p(x|mu,s)
		 *
		 * @param[in] generator		a uniform random number generator
		 * @return					probability (value between 0 and 1)
		 */
		double operator()(double & value) const
		{
			assert(false);
			return -1;
		}
};
