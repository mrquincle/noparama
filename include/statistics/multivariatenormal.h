#pragma once

#include <Eigen/Dense>

#include <random>

#include <statistics/distribution.h>

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
 *   		std::cout << distribution(generator).transpose() << std::endl;
 *   	}
 *   }
 */
class multivariate_normal_distribution: public distribution_t {
	private:
		Eigen::VectorXd _mean;
		Eigen::MatrixXd _covar;
		Eigen::MatrixXd _transform;
	public:
		/**
		 * Constructor for a multivariate normal distribution. The mean is assumed to be a zero-vector with 1xD values.
		 *
		 * @param[in] covar		matrix DxD with doubles (make sure it is a covariance matrix!)
		 */
		multivariate_normal_distribution(Eigen::MatrixXd const& covar)
			: multivariate_normal_distribution(Eigen::VectorXd::Zero(covar.rows()), covar) 
		{
		}

		/**
		 * Constructor for a multivariate normal distribution.
		 *
		 * @param[in] mean		vector 1xD with doubles
		 * @param[in] covar		matrix DxD with doubles (make sure it is a covariance matrix!)
		 */
		multivariate_normal_distribution(Eigen::VectorXd const& mean, Eigen::MatrixXd const& covar)
			: distribution_t(), _mean(mean), _covar(covar)
		{
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(_covar);
			_transform = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
		}

		/**
		 * Constructor for a multivariate normal distribution.
		 *
		 * @param[in] suffies		mean and covariance matrix in one structure
		 */
		multivariate_normal_distribution(Suffies_MultivariateNormal & suffies)
			: multivariate_normal_distribution(suffies.mu, suffies.sigma)
		{
		}

		/**
		 * Sample a random value according to the multivariate normal distribution by calling the () operator.
		 *
		 * @param[in] generator		a uniform random number generator
		 * @return					a vector 1xD with doubles, representing a single sample from the distribution
		 */
		template<typename _UniformRandomNumberGenerator >
		Eigen::VectorXd operator()(_UniformRandomNumberGenerator & generator) const
		{
			static std::normal_distribution<> dist;

			return _mean + _transform * 
				Eigen::VectorXd{ _mean.size() }.unaryExpr([&](auto x) { return dist(generator); });
		}
		
		template<typename _UniformRandomNumberGenerator >
		Suffies sample(_UniformRandomNumberGenerator &generator) const
		{
			Suffies_Unity_MultivariateNormal suffies;
			static std::normal_distribution<> dist;

			suffies.mu = _mean + _transform * 
				Eigen::VectorXd{ _mean.size() }.unaryExpr([&](auto x) { return dist(generator); });
			return suffies;
		}

		/**
		 * Calculate the probability of a value given mean and covariance.
		 *
		 * Mathematically, with mu, the mean and S the covariance matrix, this is:
		 *
		 *   p(x|mu,S)
		 *
		 * @param[in] generator		a uniform random number generator
		 * @return					probability (value between 0 and 1)
		 */
		double probability(Eigen::VectorXd & value) const
		{
			size_t D = value.size();

			const auto diff = value - _mean;
			auto inverse = _covar.inverse();
			auto exponent = -0.5 * diff.transpose() * inverse * diff;
			auto det = _covar.determinant();
			auto constant = std::pow( 2*M_PI * det, -0.5*D ) ; // check how to get det out of llt
			return constant * std::exp(exponent);		
		}

		/**
		 * Calculate the likelihood of mean and covariance given a series of observations.
		 *
		 * Mathematically, with mu, the mean and S the covariance matrix, this is:
		 *
		 *   L(mu,S|x0,x1,...,xN)
		 *
		 * @param[in] values		a vector of supposedly normally distribution multivariate variables
		 * @return					likelihood (unnormalized probability).
		 */
/*		double operator()(std::vector< Eigen::VectorXd&> & values) const
		{
			return 0;
		}
*/

};
