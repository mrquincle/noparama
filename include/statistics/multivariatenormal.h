#pragma once

#include <Eigen/Dense>

#include <random>

#include <statistics/distribution.h>

#include <iostream>

#include <pretty_print.hpp>

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
		Eigen::VectorXd _mean;
		Eigen::MatrixXd _covar;
		Eigen::MatrixXd _transform;
		int _D;

		Suffies_MultivariateNormal & _suffies_mvn;
			
		Suffies_Unity_MultivariateNormal * _suffies_result;
	public:

		/**
		 * Constructor for a multivariate normal distribution. The mean is assumed to be a zero-vector with 1xD values.
		 *
		 * @param[in] covar            matrix DxD with doubles (make sure it is a covariance matrix!)
		 */
		//multivariate_normal_distribution(Eigen::MatrixXd const& covar)
	//		: multivariate_normal_distribution(Eigen::VectorXd::Zero(covar.rows()), covar) 
	//	{
	//	}

		/**
		 * Constructor for a multivariate normal distribution.
		 *
		 * @param[in] mean             vector 1xD with doubles
		 * @param[in] covar            matrix DxD with doubles (make sure it is a covariance matrix!)
		 */
	//	multivariate_normal_distribution(Eigen::VectorXd const& mean, Eigen::MatrixXd const& covar): 
	//		distribution_t(), 
	//		_mean(mean), 
	//		_covar(covar), 
	//		_suffies_mvn1(mean.size())
		/**
		 * Constructor for a multivariate normal distribution.
		 *
		 * @param[in] suffies          mean and covariance matrix in one structure
		 */
		multivariate_normal_distribution(Suffies_MultivariateNormal & suffies)
			//: multivariate_normal_distribution(suffies.mu, suffies.sigma)
			: _suffies_mvn(suffies), _suffies_result(NULL)
		{
			_distribution_type = MultivariateNormal;
			prepare();
		}

		void init(Suffies & suffies) {
			_suffies_mvn = (Suffies_MultivariateNormal&)suffies;
			prepare();
		}

		/**
		 * Assume _suffies_mvn is set. Prepare the other fields.
		 */
		void prepare() {
			_mean = _suffies_mvn.mu;
			_covar = _suffies_mvn.sigma;
			_D = _mean.size();
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(_covar);
			_transform = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
		}

		/**
		 * Sample a random value according to the multivariate normal distribution by calling the () operator.
		 *
		 * @param[in] generator        a uniform random number generator
		 * @return                     a suffies with a vector 1xD with doubles, representing a single sample from the
		 *                             distribution
		 */
		Suffies_Unity_MultivariateNormal* operator()(random_engine_t &generator) 
		{
			_suffies_result = new Suffies_Unity_MultivariateNormal(_D);

			static std::normal_distribution<> dist;

			_suffies_result->mu = _mean + _transform * 
				Eigen::VectorXd{ _mean.size() }.unaryExpr([&](auto x) { return dist(generator); });
			return _suffies_result;
		}

		/**
		 * Calculate the probability of a value given mean and covariance.
		 *
		 * Mathematically, with mu, the mean and S the covariance matrix, this is:
		 *
		 *   p(x|mu,S)
		 *
		 * @param[in] generator        a uniform random number generator
		 * @return                     probability (value between 0 and 1)
		 */
		double probability(data_t & data) const
		{
			double* ptr = &data[0];
			Eigen::Map< Eigen::VectorXd > value(ptr, data.size());
			
			size_t D = value.size();

			const auto diff = value - _mean;
//			fout << diff.transpose() << std::endl;
			auto inverse = _covar.inverse();
			auto exponent = -0.5 * diff.transpose() * inverse * diff;
//			std::cout << "Exponent: " << exponent << std::endl;
			double det = _covar.determinant();
//			std::cout << "Determinant: " << det << std::endl;
			//auto constant = std::pow( 2*M_PI * det, -0.5*D ) ; 
			auto constant = std::sqrt( std::pow( 2*M_PI, D ) * det); // check how to get det out of llt
//			std::cout << "Constant: " << constant << std::endl;

			return std::exp(exponent) / constant;
		}
		
		double probability(dataset_t & dataset) const
		{
			assert (dataset.size() > 0);
			double result = 1.0; 
			for (int i = 0; i < (int)dataset.size(); ++i) {
				result *= probability(*dataset[i]);
			}
			return result;
		}

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
