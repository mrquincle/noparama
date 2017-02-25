#pragma once

#include <Eigen/Dense>

#include <random>
#include <iostream>

#include <statistics/distribution.h>

#include <np_suffies.h>
#include <pretty_print.hpp>



/*!
 * To sample from a Dirichlet process, we need to initialize it with a dispersion factor and a base distribution
 * and then get samples out using a uniform random generator as input.
 *
 * The sufficient statistics that are used as parameters are stored by reference. The base distribution is also
 * stored by reference.
 */
class dirichlet_process: public distribution_t {
	private:
		/*! 
		 * The sufficient statistics for the Dirichlet, should actually store a reference to a distribution. If we 
		 * just store the sufficient statistics of the base distribution, it is not enough to sample from them because
		 * the sample() function is defined for the distribution_t object, not the suffies. 
		 *
		 * Hence, currently it doesn't store anything about the base distribution. This is kept into this distribution
		 * object.
		 */
		Suffies_Dirichlet & _suffies_dirichlet;

		distribution_t & _base_distribution;

		int _verbosity;
	public:
		dirichlet_process(Suffies_Dirichlet & suffies_dirichlet, distribution_t & base_distribution):
			_suffies_dirichlet(suffies_dirichlet), _base_distribution(base_distribution)
		{
			_distribution_type = Dirichlet;
			_verbosity = 7;
		}

		void init(Suffies & suffies) {
			_suffies_dirichlet = (Suffies_Dirichlet&)suffies;
		}

		Suffies_Dirichlet & getSuffies() {
			return _suffies_dirichlet;
		}

		/**
		 * Return sample from base distribution. It is not known what type the sufficient statistics will be. 
		 */
		Suffies* operator()(random_engine_t & generator) {
			std::cout << "It is not possible to just sample from a Dirichlet Process like this." << std::endl;
			std::cout << "We would need to store an infinite number of samples to be able to draw independently from them." << std::endl;

			assert(false);
		};

		/**
		 * Draw X_0 from sample(generator). Draw X_i with i > 0 from: 
		 *  sample(generator) with probability alpha / (alpha + N).
		 *  previous sample with 1 / (alpha + N)
		 */
		Suffies_MultivariateNormal* sample(
				random_engine_t & generator, 
				const std::vector<Suffies_MultivariateNormal*> & history) {
			int N = history.size();
			if (!N) return sample_base(generator);

			std::uniform_real_distribution<double> dist_uniform(0.0,1.0);
			double threshold = _suffies_dirichlet.alpha / (_suffies_dirichlet.alpha + N);
			double pick = dist_uniform(generator);
			if (pick < threshold) {
				return sample_base(generator);
			} else {
				std::uniform_int_distribution<> dist_categorical(0,N-1);
				int k = dist_categorical(generator);

				return history[k];
			}
		}

		Suffies_MultivariateNormal* sample_base(random_engine_t & generator) {
			return (Suffies_MultivariateNormal*)_base_distribution(generator);
		}
};

