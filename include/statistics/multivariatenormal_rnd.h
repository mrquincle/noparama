#include <Eigen/Dense>
#include <random>

class multivariate_normal_distribution {
	private:
		Eigen::VectorXd mean;
		Eigen::MatrixXd transform;

	public:

		multivariate_normal_distribution(Eigen::MatrixXd const& covar)
			: multivariate_normal_distribution(Eigen::VectorXd::Zero(covar.rows()), covar) 
		{
		}

		multivariate_normal_distribution(Eigen::VectorXd const& mean, Eigen::MatrixXd const& covar)
			: mean(mean)
		{
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
			transform = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
		}

		template<typename _UniformRandomNumberGenerator >
		Eigen::VectorXd operator()(_UniformRandomNumberGenerator & generator) const
		{
			static std::normal_distribution<> dist;

			return mean + transform * Eigen::VectorXd{ mean.size() }.unaryExpr([&](auto x) { return dist(generator); });
		}
};
