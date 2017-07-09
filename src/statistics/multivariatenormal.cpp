#include <statistics/multivariatenormal.h>

multivariate_normal_distribution::multivariate_normal_distribution(
		Suffies_MultivariateNormal & suffies, 
		bool regression): _suffies_result(NULL), _regression(regression) 
{
	_suffies_mvn = &suffies;
	_distribution_type = MultivariateNormal;
	prepare();
}

void multivariate_normal_distribution::setRegression(bool regression) {
	_regression = regression;
}

void multivariate_normal_distribution::init(Suffies & suffies) {
	if (typeid(suffies)!=typeid(*_suffies_mvn)) {
		std::cout << "init mvn" << std::endl;
		std::cout << "current suffies: " << *_suffies_mvn << std::endl;
		std::cout << "incoming suffies: " << suffies << std::endl;
		assert (typeid(suffies)==typeid(*_suffies_mvn));
	}
	_suffies_mvn = &dynamic_cast<Suffies_MultivariateNormal&>(suffies);
	prepare();
}

Suffies & multivariate_normal_distribution::getSuffies() {
	return *_suffies_mvn;
}

void multivariate_normal_distribution::prepare() {
	_mean = _suffies_mvn->mu;
	_covar = _suffies_mvn->sigma;
	_D = _mean.size();
}

Suffies_Unity_MultivariateNormal* multivariate_normal_distribution::operator()(random_engine_t &generator) 
{
	_suffies_result = new Suffies_Unity_MultivariateNormal(_D);

	static std::normal_distribution<> dist;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(_covar);
	//	std::cout << "eigenSolver" << std::endl;
	_transform = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();

	_suffies_result->mu = _mean + _transform * 
		Eigen::VectorXd{ _mean.size() }.unaryExpr([&](auto x) { return dist(generator); });
	return _suffies_result;
}

double multivariate_normal_distribution::probability(data_t & data) const
{
	double* ptr = &data[0];
	size_t D;
	if (_regression) {
		D = data.size() - 1;
		double yd = data[D]; 
		Eigen::VectorXd y(1); y << yd;
		Eigen::Map< Eigen::VectorXd > value(ptr, D);
		auto y_proj = value.transpose() * _mean;
		const auto diff = y - y_proj;

		auto inverse = _covar.inverse();
		assert (inverse.cols() == diff.rows());
		auto exponent = -0.5 * diff.transpose() * inverse * diff;
		double det = _covar.determinant();
		auto constant = std::sqrt( std::pow( 2*M_PI, D ) * det); 
		return std::exp(exponent) / constant;
	} else {
		D = data.size();
		Eigen::Map< Eigen::VectorXd > value(ptr, D);
		const auto diff = value - _mean;

		auto inverse = _covar.inverse();
		assert (inverse.cols() == diff.rows());
		auto exponent = -0.5 * diff.transpose() * inverse * diff;
		double det = _covar.determinant();
		auto constant = std::sqrt( std::pow( 2*M_PI, D ) * det); 
		return std::exp(exponent) / constant;
	}
}

double multivariate_normal_distribution::probability(dataset_t & dataset) const
{
	assert (dataset.size() > 0);
	double result = 1.0; 
	for (int i = 0; i < (int)dataset.size(); ++i) {
		result *= probability(*dataset[i]);
	}
	return result;
}

double multivariate_normal_distribution::logprobability(data_t & data) const
{
	double* ptr = &data[0];
	size_t D;
	if (_regression) {
		D = data.size() - 1;
		double yd = data[D]; 
		Eigen::VectorXd y(1); y << yd;
		Eigen::Map< Eigen::VectorXd > value(ptr, D);
		auto y_proj = value.transpose() * _mean;
		const auto diff = y - y_proj;

		auto inverse = _covar.inverse();
		assert (inverse.cols() == diff.rows());
		auto exponent = -0.5 * diff.transpose() * inverse * diff;
		double det = _covar.determinant();
		auto constant = std::sqrt( std::pow( 2*M_PI, D ) * det); 
		return exponent - std::log(constant);
	} else {
		D = data.size();
		Eigen::Map< Eigen::VectorXd > value(ptr, D);
		const auto diff = value - _mean;

		auto inverse = _covar.inverse();
		assert (inverse.cols() == diff.rows());
		auto exponent = -0.5 * diff.transpose() * inverse * diff;
		double det = _covar.determinant();
		auto constant = std::sqrt( std::pow( 2*M_PI, D ) * det); 
		return exponent - std::log(constant);
	}
}

double multivariate_normal_distribution::logprobability(dataset_t & dataset) const
{
	assert (dataset.size() > 0);
	double result = 0.0; 
	for (int i = 0; i < (int)dataset.size(); ++i) {
		result += logprobability(*dataset[i]);
	}
	return result;
}
