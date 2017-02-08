#pragma once

#include <Eigen/Dense>

#include <random>

#include <statistics/distribution.h>

/*!
 * To sample from a normal distribution, we need to initialize it with mean and variance
 * and then get samples out using a uniform random generator as input.
 */
class normal_distribution: public distribution_t {
	private:
		double _mean;
		double _variance;

		Suffies_Unity_Normal *_suffies_result;
	public:
		/**
		 * Constructor for a normal distribution. The mean is assumed to be a zero.
		 *
		 * @param[in] variance		double
		 */
		normal_distribution(double variance)
			: normal_distribution(0, variance) 
		{
		}

		/**
		 * Constructor for a normal distribution.
		 *
		 * @param[in] mean			double
		 * @param[in] variance		double
		 */
		normal_distribution(double mean, double variance)
			: distribution_t(), _mean(mean), _variance(variance)
		{
		}

		/**
		 * Constructor for a normal distribution.
		 *
		 * @param[in] suffies		mean and variance in one structure
		 */
		normal_distribution(Suffies_Normal & suffies)
			: normal_distribution(suffies.mu, suffies.sigma)
		{
		}

		/**
		 * Sample a random value according to the normal distribution by calling the () operator.
		 *
		 * @param[in] generator		a uniform random number generator
		 * @return					a double, representing a single sample from the distribution
		 */
		Suffies_Unity_Normal* operator()(random_engine_t & generator) 
		{
			_suffies_result = new Suffies_Unity_Normal();

			static std::normal_distribution<> dist(_mean, _variance);
			
			_suffies_result->mu = dist(generator);

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
			return -1;
			/*
			size_t D = value.size();

			const auto diff = value - _mean;
			auto inverse = _variance.inverse();
			auto exponent = -0.5 * diff.transpose() * inverse * diff;
			auto det = _variance.determinant();
			auto constant = std::pow( 2*M_PI * det, -0.5*D ) ; // check how to get det out of llt
			return constant * std::exp(exponent);		 */
		}
};
