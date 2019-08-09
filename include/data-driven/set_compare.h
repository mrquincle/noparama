#pragma once

#include <Eigen/Dense>

#include <random>

#include <statistics/distribution.h>

#include <iostream>

#include <pretty_print.hpp>
#include <typeinfo>

/*!
 * The set_compare class is not a true distribution. It can just return likelihood values.
 */
class set_compare: public distribution_t {
	private:
		
		Eigen::VectorXd _mean;
		
		int _D;
		
		dataset_t _dataset_reference;

		float * _dataset_reference_raw;
		
		float * _dataset_reference_mean;
		
		float * _dataset_compare_raw;
		
		float * _dataset_compare_mean;
		
		float * _dataset_match;

		Suffies_ScalarNoise_MultivariateNormal * _suffies_mvn;

	protected:


		void calc_mean(float * data, int size, int dim, float * result) const;

		void get_shape(const dataset_t & dataset, int & size, int & dim) const;

		void get_raw(dataset_t & source, float* target) const;

	public:

		/**
		 * Constructor for a set comparison.
		 *
		 * The suffies argument is by reference to prevent unnecessary copy actions from one suffies to the next. 
		 * However, be very careful! Do not pass accidentally a temporary variable towards the distribution that goes
		 * out of scope after calling the constructor.
		 *
		 * @param[in] suffies          mean 
		 */
		set_compare(Suffies_ScalarNoise_MultivariateNormal & suffies, dataset_t & dataset_reference);
		
		void init(Suffies & suffies);

		Suffies & getSuffies();

		//! Do not use
		double probability(data_t & data) const;
		
		//! Do not use
		double logprobability(data_t & data) const;
	
		/**
		 * Calculate the probability of a dataset given the original dataset and the mean.
		 *
		 * @param[in] dataset          multiple data items
		 * @return                     probability (value between 0 and 1)
		 */
		double probability(dataset_t & dataset) const;
	
		/**
		 * Log probability, see function above.
		 */
		double logprobability(dataset_t & dataset) const;
};
