#include <statistics/scalarnoise_multivariatenormal.h>

scalarnoise_multivariate_normal_distribution::scalarnoise_multivariate_normal_distribution(
		Suffies_ScalarNoise_MultivariateNormal & suffies, 
		representation_mode_t mode): _suffies_result(NULL), _mode(mode) 
{
	_suffies_mvn = &suffies;
	_distribution_type = ScalarNoise_MultivariateNormal;
	prepare();
}

void scalarnoise_multivariate_normal_distribution::setMode(representation_mode_t mode) {
	_mode = mode;
}

void scalarnoise_multivariate_normal_distribution::init(Suffies & suffies) {
	if (typeid(suffies)!=typeid(*_suffies_mvn)) {
		//std::cout << "init mvn" << std::endl;
		//std::cout << "current suffies: " << *_suffies_mvn << std::endl;
		//std::cout << "incoming suffies: " << suffies << std::endl;
		assert (typeid(suffies)==typeid(*_suffies_mvn));
	}
	_suffies_mvn = &dynamic_cast<Suffies_ScalarNoise_MultivariateNormal&>(suffies);
	prepare();
}

Suffies & scalarnoise_multivariate_normal_distribution::getSuffies() {
	return *_suffies_mvn;
}

void scalarnoise_multivariate_normal_distribution::prepare() {
	_mean = _suffies_mvn->mu;
	_var = _suffies_mvn->sigma;
	_D = _mean.size();
	switch (_mode) {
		case angular_mode: 
			{
				// first array element is a distance: >=0
				_mean[0] = abs(_mean[0]);
				// second array element is an angle 0 <= alpha < 2pi
				_mean[1] = fmod(abs(_mean[1]), 2*M_PI);
				break;
			}
		default:
			{
			// nothing
			}
	}
}

/**
 * We return the results of a sampling action by a Suffies struct again. It would also be possible to just return
 * data_t here. The returned struct only has the mu field set to nonzero values. The covariance matrix is the
 * identity matrix times a variance constant. This means we can sample each dimension independently. 
 */
Suffies_Unity_MultivariateNormal* scalarnoise_multivariate_normal_distribution::operator()(random_engine_t &generator) 
{
	_suffies_result = new Suffies_Unity_MultivariateNormal(_D);

	static std::normal_distribution<> dist;

	std::cout << "Before: " << _mean.transpose() << std::endl;
	_suffies_result->mu = _mean;
	for (int i = 0; i < _D; ++i) {
		double shift = _var * dist(generator);
		_suffies_result->mu[i] += shift;
	}
	return _suffies_result;
}

/**
 * If _mode == regression, the data is organized with [x_0 x_1 ... x_n y] with x_i the independent variables and y the
 * dependent variable.
 * If _mode == angular, the data is organized as [theta_0 theta_1 ... theta_n D] with theta_i the angles w.r.t. 
 * subsequent dimensions and D the distance to the origin in an R^{n+1} dimensional space.
 */
double scalarnoise_multivariate_normal_distribution::probability(data_t & data) const
{
	double* ptr = &data[0];
	size_t D;
	switch (_mode) {
		case regression_mode: 
			{
				// assume that the last item is y
				D = data.size() - 1;
				double yd = data[D]; 
				Eigen::VectorXd y(1); y << yd;
				Eigen::Map< Eigen::VectorXd > value(ptr, D);
				auto y_proj = value.transpose() * _mean;
				assert(y_proj.size() == 1);
				Eigen::VectorXd y_diff = y - y_proj;
				double diff = y_diff(0);

				double inverse = 1/(_var*_var);
				double exponent = -0.5 * (diff * diff) * inverse; 
				auto constant = std::sqrt( 2*M_PI*_var*_var ); 
				return std::exp(exponent) / constant;
			}
		case clustering_mode:
			{
				D = data.size();
				Eigen::Map< Eigen::VectorXd > value(ptr, D);
				const auto diff = value - _mean;

				Eigen::MatrixXd covar = Eigen::MatrixXd::Identity(D,D).array() * _var;

				auto inverse = covar.inverse();
				assert (inverse.cols() == diff.rows());
				auto exponent = -0.5 * diff.transpose() * inverse * diff;
				double det = covar.determinant();
				auto constant = std::sqrt( std::pow( 2*M_PI, D ) * det); 
				return std::exp(exponent) / constant;
			}
		case angular_mode: 
			{
				D = data.size();
				assert (D == 2);
				double d = _mean[0];
				double theta = _mean[1];
				// Clockwise rotation matrix
				Eigen::Matrix2d m;
				m << cos(theta), sin(theta),
					-sin(theta), cos(theta);
			
				Eigen::Map< Eigen::VectorXd > value(ptr, D);
				const auto q = m * value;
				double diff = std::abs( d - q(1));
				
				double inverse = 1/(_var*_var);
				double exponent = -0.5 * (diff * diff) * inverse; 
				auto constant = std::sqrt( 2*M_PI*_var*_var ); 
				return std::exp(exponent) / constant;
			}
		default:
			{
				std::cerr << "Unknown mode" << std::endl;
				return -1;
			}
	}
}

// Numbers are always compared in a ratio where the number of data points is exactly the same in the nominator and
// the denominator, so multiplying with the number of data points is unnecessary. This is especially easy to see 
// with loglikelihoods, where we just add some constants above and below the fraction bar.

// Do or do not factor in size.
//#define FACTOR_IN_SIZE

double scalarnoise_multivariate_normal_distribution::probability(dataset_t & dataset) const
{
	assert (dataset.size() > 0);
	double result = 1.0; 
	for (int i = 0; i < (int)dataset.size(); ++i) {
		result *= probability(*dataset[i]);
	}
#ifdef FACTOR_IN_SIZE
	return result * dataset.size();
#else
	return result; 
#endif
}

double scalarnoise_multivariate_normal_distribution::logprobability(data_t & data) const
{
	double* ptr = &data[0];
	size_t D;
	switch (_mode) {
		case regression_mode: 
			{
				// in 2D, data.size() should be 3, y = data[2], ptr = data[0-1], with data[0] = 1.
				D = data.size() - 1;
				double yd = data[D]; 
				Eigen::VectorXd y(1); y << yd;
				Eigen::Map< Eigen::VectorXd > value(ptr, D);
				auto y_proj = value.transpose() * _mean;
				assert(y_proj.size() == 1);
				Eigen::VectorXd y_diff = y - y_proj;
				double diff = y_diff(0);

				double inverse = 1/(_var*_var);
				double exponent = -0.5 * (diff * diff) * inverse; 
				auto constant = 2*M_PI*_var*_var; 
				return exponent - std::log(constant)/2;
				//auto constant = std::sqrt( 2*M_PI*_var*_var ); 
				//return exponent - std::log(constant);
			} 
		case clustering_mode:
			{
				D = data.size();
				Eigen::Map< Eigen::VectorXd > value(ptr, D);
				const auto diff = value - _mean;

				Eigen::MatrixXd covar = Eigen::MatrixXd::Identity(D,D).array() * _var;

				auto inverse = covar.inverse();
				assert (inverse.cols() == diff.rows());
				auto exponent = -0.5 * diff.transpose() * inverse * diff;
				double det = covar.determinant();
				auto constant = std::sqrt( std::pow( 2*M_PI, D ) * det); 
				return exponent - std::log(constant);
			}
		case angular_mode:
			{
				// angular mode, polar mode, normal form, etc. is not easy to extend to higher dimensional spaces
				// already in 4D is rotation not easy to define: each rotation has at least two invariant axis-planes
				D = data.size();
				assert (D == 2);
				double d = _mean[0];
				double theta = _mean[1];
				//std::cout << "theta: " << theta << std::endl; 
				//std::cout << "d: " << d << std::endl; 
				// Clockwise rotation matrix
				Eigen::Matrix2d m;
				m << cos(theta), sin(theta),
					-sin(theta), cos(theta);
				//std::cout << "m: " << m << std::endl; 
			
				Eigen::Map< Eigen::VectorXd > value(ptr, D);
				//std::cout << "value: " << value << std::endl; 
				const auto q = m * value;
				//std::cout << "q: " << q << std::endl; 
				double diff = std::abs( d - q(1));
				//std::cout << "diff: " << diff << std::endl; 
				
				double inverse = 1/(_var*_var);
				double exponent = -0.5 * (diff * diff) * inverse; 
				auto constant = 2*M_PI*_var*_var; 
				return exponent - std::log(constant)/2;
			}
		default:
			{
				std::cerr << "Unknown mode" << std::endl;
				return -1;
			}
	}
}

double scalarnoise_multivariate_normal_distribution::logprobability(dataset_t & dataset) const
{
	assert (dataset.size() > 0);
	double result = 0.0; 
	for (int i = 0; i < (int)dataset.size(); ++i) {
		result += logprobability(*dataset[i]);
	}
#ifdef FACTOR_IN_SIZE
	return result + dataset.size();
#else
	return result;
#endif
}
