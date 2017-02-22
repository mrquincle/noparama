#pragma once

#include <Eigen/Dense>

#include <random>
#include <iostream>

#include <statistics/distribution.h>

#include <np_suffies.h>
#include <pretty_print.hpp>

/*!
 * To sample from a Dirichlet distribution, we need to initialize it with a dispersion factor and a base distribution
 * and then get samples out using a uniform random generator as input.
 *
 * The sufficient statistics that are used as parameters are stored by reference. The base distribution is also
 * stored by reference.
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
		Suffies_Dirichlet & _suffies_dirichlet;

		distribution_t & _base_distribution;

	public:
		dirichlet_distribution(Suffies_Dirichlet & suffies_dirichlet, distribution_t & base_distribution):
			_suffies_dirichlet(suffies_dirichlet), _base_distribution(base_distribution)
		{
			_distribution_type = Dirichlet;
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
		//	fout << "Are you sure you need to sample from the Dirichlet?" << std::endl;
			return _base_distribution(generator);
		};
};

