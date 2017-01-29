/*!
 * SufficientStatistics describes a probability distribution.
 */
class SufficientStatistics {

}

class SufficientStatistics_NormalInvWishart: SufficientStatistics {
	public:
		double mu;
		double kappa;
		double nu;
		double lambda;
}

class SufficientStatistics_NormalInvGamma: SufficientStatistics {
	public:
		double a;
		double b;
		double nu;
		Matrix2d Lambda;
}

class SufficientStatistics_Normal: SufficientStatistics {
	public:
		double mu;
		double sigma;
}

class SufficientStatistics_MultivariateNormal: SufficientStatistics {
	public:
		Eigen::VectorXd mu;
		Eigen::MatrixXd sigma;
}
