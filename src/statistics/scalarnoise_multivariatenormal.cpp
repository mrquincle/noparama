#include <statistics/scalarnoise_multivariatenormal.h>

scalarnoise_multivariate_normal_distribution::scalarnoise_multivariate_normal_distribution(
		Suffies_ScalarNoise_MultivariateNormal & suffies, 
		bool regression): _suffies_result(NULL), _regression(regression) 
{
	_suffies_mvn = &suffies;
	_distribution_type = ScalarNoise_MultivariateNormal;
	prepare();
}

void scalarnoise_multivariate_normal_distribution::setRegression(bool regression) {
	_regression = regression;
}

void scalarnoise_multivariate_normal_distribution::init(Suffies & suffies) {
	if (typeid(suffies)!=typeid(*_suffies_mvn)) {
		std::cout << "init mvn" << std::endl;
		std::cout << "current suffies: " << *_suffies_mvn << std::endl;
		std::cout << "incoming suffies: " << suffies << std::endl;
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
}

Suffies_Unity_MultivariateNormal* scalarnoise_multivariate_normal_distribution::operator()(random_engine_t &generator) 
{
	_suffies_result = new Suffies_Unity_MultivariateNormal(_D);

	static std::normal_distribution<> dist;

	double scale = _var * dist(generator);
	_suffies_result->mu = _mean.array() * scale;
	return _suffies_result;
}

double scalarnoise_multivariate_normal_distribution::probability(data_t & data) const
{
	double* ptr = &data[0];
	size_t D;
	if (_regression) {
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
	} else {
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
}

double scalarnoise_multivariate_normal_distribution::probability(dataset_t & dataset) const
{
	assert (dataset.size() > 0);
	double result = 1.0; 
	for (int i = 0; i < (int)dataset.size(); ++i) {
		result *= probability(*dataset[i]);
	}
	return result;
}

double scalarnoise_multivariate_normal_distribution::logprobability(data_t & data) const
{
	double* ptr = &data[0];
	size_t D;
	if (_regression) {
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
		return exponent - std::log(constant);
	} else {
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
}

double scalarnoise_multivariate_normal_distribution::logprobability(dataset_t & dataset) const
{
	assert (dataset.size() > 0);
	double result = 0.0; 
	for (int i = 0; i < (int)dataset.size(); ++i) {
		result += logprobability(*dataset[i]);
	}
	return result;
}
