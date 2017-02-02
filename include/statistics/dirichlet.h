#pragma once

#include <Eigen/Dense>

#include <random>

#include <statistics/distribution.h>

/*!
 * To sample from a Dirichlet distribution, we need to initialize it with a dispersion factor and a base distribution
 * and then get samples out using a uniform random generator as input.
 */
class dirichlet_distribution: public distribution_t {
	private:
		double _alpha;
		distribution_t _base_distribution;
	public:
		dirichlet_distribution(double alpha, distribution_t & base_distribution):
			_alpha(alpha), _base_distribution(base_distribution) 
		{
		}

		//! TODO: Here we run into trouble... The above has a base distribution as parameter, from which can be
		dirichlet_distribution(Suffies_Dirichlet suffies) //:
//			dirichlet_distribution(suffies.alpha, suffies.base_suffies) 
		{
		}
};

