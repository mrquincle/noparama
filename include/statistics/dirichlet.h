#pragma once

#include <Eigen/Dense>

#include <random>

#include <statistics/distribution.h>

#include <np_suffies.h>

/*!
 * To sample from a Dirichlet distribution, we need to initialize it with a dispersion factor and a base distribution
 * and then get samples out using a uniform random generator as input.
 */
class dirichlet_distribution: public distribution_t {
	private:
		/*! 
		 * The sufficient statistics for the Dirichlet, should actually store a reference to a distribution. If we 
		 * just store the sufficient statistics of the base distribution, it is not enough to sample from them because
		 * the sample() function is defined for the distribution_t object, not the suffies. 
		 *
		 * Hence, currently it doesn't store anything about the base distribution. This is kept into this distribution
		 * object.
		 */
		Suffies_Dirichlet _suffies_dirichlet;

		distribution_t _base_distribution;
	public:
		dirichlet_distribution(Suffies_Dirichlet suffies_dirichlet, distribution_t & base_distribution):
			_suffies_dirichlet(suffies_dirichlet), _base_distribution(base_distribution)
		{

		}

		void init(Suffies & suffies) {
			_suffies_dirichlet = (Suffies_Dirichlet&)suffies;
		}

		Suffies_Dirichlet & getSuffies() {
			return _suffies_dirichlet;
		}
};

