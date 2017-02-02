#pragma once

#include <Eigen/Dense>

/*!
 * Suffies describes a probability distribution. On the moment it does not contain any information on itself.
 *
 * Preferably the design should be such that there is not some kind of type field (to cast to derived classes). That
 * normally namely stinks of bad design. It might be okay to have a string or something like that for cosmetic
 * reasons (printing e.g.).
 */
class Suffies {

};

class Suffies_NormalInvWishart: public Suffies {
	public:
		Eigen::VectorXd mu;
		double kappa;
		double nu;
		Eigen::Matrix2d Lambda;
};

class Suffies_InvWishart: public Suffies {
	public:
		double nu;
		double lambda;
		Eigen::Matrix2d Lambda;
};

class Suffies_NormalInvGamma: public Suffies {
	public:
		Eigen::VectorXd mu;
		double a;
		double b;
		Eigen::Matrix2d Lambda;
};

class Suffies_Normal: public Suffies {
	public:
		double mu;
		double sigma;
};

/*!
 * A multivariate normal distribution with an identity matrix as covariance / information matrix.
 */
class Suffies_Unity_MultivariateNormal: public Suffies {
	public:
		Eigen::VectorXd mu;
};

/*!
 * A multivariate normal disribution centered at 0.
 */
class Suffies_ZeroCentered_MultivariateNormal: public Suffies {
	public:
		Eigen::MatrixXd sigma;
};

class Suffies_MultivariateNormal: public Suffies {
	public:
		Eigen::VectorXd mu;
		Eigen::MatrixXd sigma;
};

class Suffies_Dirichlet: public Suffies {
	public:
		double alpha;
		Suffies base_suffies;
};
