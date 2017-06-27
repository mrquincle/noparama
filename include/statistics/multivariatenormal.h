#pragma once

#include <Eigen/Dense>

#include <random>

#include <statistics/distribution.h>

#include <iostream>

#include <pretty_print.hpp>
#include <typeinfo>

/*!
 * To sample from a multivariate normal distribution, we need to initialize it with mean and covariance matrix
 * and then get samples out using a uniform random generator as input.
 *
 * There is no particular exception handling for the cases in which the user uses mean and covariance vectors and 
 * matrices with unmatching dimensions. Eigen will throw errors in that case.
 * 
 * Usage with mean [1 1] and covariance matrix [1 0; 0.2 1]:
 *
 *   int main() {
 *   	Vector2d mean(1, 1);
 *   	Matrix2d covar; covar << 1, 0, 0.2, 1;
 *   			
 *   	multivariate_normal_distribution distribution(mean, covar);
 *   	default_random_engine generator;
 *   
 *   	const int nRolls = 100;
 *   	for (int i = 0; i < nRolls; ++i) {
 *   		std::fout << distribution(generator).transpose() << std::endl;
 *   	}
 *   }
 */
class multivariate_normal_distribution: public distribution_t {
	private:
		//! Store result in the form 
		Eigen::VectorXd _mean;
		Eigen::MatrixXd _covar;
		Eigen::MatrixXd _transform;
		int _D;

		Suffies_MultivariateNormal * _suffies_mvn;
			
		Suffies_Unity_MultivariateNormal * _suffies_result;
	
		//! Regression flag
		bool _regression;
	public:

		/**
		 * Constructor for a multivariate normal distribution.
		 *
		 * The suffies argument is by reference to prevent unnecessary copy actions from one suffies to the next. 
		 * However, be very careful! Do not pass accidentally a temporary variable towards the distribution that goes
		 * out of scope after calling the constructor.
		 *
		 * @param[in] suffies          mean and covariance matrix in one structure
		 * @param[in] regression       handle data as regression data [y = X*d]
		 *
		 * TODO: Find out a way to communicate that the caller has to take care of the lifetime of suffies 
		 */
		multivariate_normal_distribution(Suffies_MultivariateNormal & suffies, bool regression = false);
		
		inline void setRegression(bool regression);

		void init(Suffies & suffies);

		Suffies & getSuffies();

		/**
		 * Assume _suffies_mvn is set. Prepare the other fields.
		 */
		void prepare();

		/**
		 * Sample a random value according to the multivariate normal distribution by calling the () operator.
		 *
		 * @param[in] generator        a uniform random number generator
		 * @return                     a suffies with a vector 1xD with doubles, representing a single sample from the
		 *                             distribution
		 */
		Suffies_Unity_MultivariateNormal* operator()(random_engine_t &generator);

		/**
		 * Calculate the probability of a value given mean and covariance.
		 *
		 * Mathematically, with mu, the mean and S the covariance matrix, this is:
		 *
		 *   p(x|mu,S)
		 *
		 * @param[in] data             individual data item
		 * @return                     probability (value between 0 and 1)
		 */
		double probability(data_t & data) const;
	
		/**
		 * Calculate the probability of multiple values given mean and covariance. 
		 *
		 * This is equivalent to the product of individual data items for the multivariate normal distribution.
		 * @param[in] dataset          multiple data items
		 * @return                     probability (value between 0 and 1)
		 */
		double probability(dataset_t & dataset) const;

		/**
		 * Calculate the likelihood of mean and covariance given a series of observations.
		 *
		 * Mathematically, with mu, the mean and S the covariance matrix, this is:
		 *
		 *   L(mu,S|x0,x1,...,xN)
		 *
		 * @param[in] values           a vector of supposedly normally distribution multivariate variables
		 * @return                     likelihood (unnormalized probability).
		 */

/*		double operator()(std::vector< Eigen::VectorXd&> & values) const
		{
			return 0;
		}
*/

};
