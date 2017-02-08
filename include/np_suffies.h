#pragma once

#include <Eigen/Dense>
#include <map>

enum distribution_type_t { 
	// normals
	Normal, Unity_Normal,
	// multivariate normals
	MultivariateNormal, Unity_MultivariateNormal, ZeroCentered_MultivariateNormal, 
	// inverse wisharts and alike, useful as priors for normals
	NormalInvWishart, InvWishart, NormalInvGamma, 
	// hierarchical priors
	Dirichlet
};

static std::map< distribution_type_t, const char * > distribution_type_str = {
	{Normal,                           "normal"},
	{Unity_Normal,                     "unity normal"},
	{MultivariateNormal,               "multivariate normal"},
	{Unity_MultivariateNormal,         "unity multivariate normal"},
	{ZeroCentered_MultivariateNormal,  "zero-centered multivariate normal"},
	{NormalInvWishart,                 "normal inverse-Wishart"},
	{InvWishart,                       "inverse-Wishart"},
	{NormalInvGamma,                   "normal inverse-Gamma"},
	{Dirichlet,                        "Dirichlet"}
};

/*!
 * Suffies describes a probability distribution. On the moment it does not contain any information on itself.
 *
 * Preferably the design should be such that there is not some kind of type field (to cast to derived classes). That
 * normally namely stinks of bad design. It might be okay to have a string or something like that for cosmetic
 * reasons (printing e.g.).
 */
class Suffies {
	public:
		int D;

		virtual ~Suffies() {
		};

		virtual void print(std::ostream& os) const {
			os << D << " (probably a copy by value rather than copy by reference somewhere?)";
		}

		friend std::ostream& operator<<(std::ostream& os, const Suffies& s);  
};

class Suffies_NormalInvWishart: public Suffies {
	public:
		Eigen::VectorXd mu;
		double kappa;
		double nu;
		Eigen::MatrixXd Lambda;
		
		Suffies_NormalInvWishart(int D): mu(D), Lambda(D,D) {
			Suffies::D = D;
		}
		
		void print(std::ostream& os) const {
			os << nu << " | " << kappa << " | " << nu << " | " << Lambda;
		}
};

class Suffies_InvWishart: public Suffies {
	public:
		double nu;
		Eigen::MatrixXd Lambda;
		
		Suffies_InvWishart(int D): Lambda(D,D) {
		}
		
		void print(std::ostream& os) const {
			os << nu << " | " << Lambda;
		}
};

class Suffies_NormalInvGamma: public Suffies {
	public:
		Eigen::VectorXd mu;
		double a;
		double b;
		Eigen::MatrixXd Lambda;
		
		Suffies_NormalInvGamma(int D): mu(D), Lambda(D,D) {
		}
		
		void print(std::ostream& os) const {
			os << mu << " | " << a << " | " << b << " | " << Lambda;
		}
};

class Suffies_Normal: public Suffies {
	public:
		double mu;
		double sigma;
		
		void print(std::ostream& os) const {
			os << mu << " | " << sigma;
		}
};

class Suffies_Unity_Normal: public Suffies {
	public:
		double mu;

		void print(std::ostream& os) const {
			os << mu;
		}
};

/*!
 * A multivariate normal distribution with an identity matrix as covariance / information matrix.
 */
class Suffies_Unity_MultivariateNormal: public Suffies {
	public:
		Eigen::VectorXd mu;

		Suffies_Unity_MultivariateNormal(int D): mu(D) {
		}
		
		void print(std::ostream& os) const {
			os << mu;
		}
};

/*!
 * A multivariate normal disribution centered at 0.
 */
class Suffies_ZeroCentered_MultivariateNormal: public Suffies {
	public:
		Eigen::MatrixXd sigma;
		
		Suffies_ZeroCentered_MultivariateNormal(int D): sigma(D, D) {
		}
		
		void print(std::ostream& os) const {
			os << sigma;
		}
};

static Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", "; ", "", "", "", "");

class Suffies_MultivariateNormal: public Suffies {
	public:
		Eigen::VectorXd mu;
		Eigen::MatrixXd sigma;
		
		Suffies_MultivariateNormal(int D): mu(D), sigma(D,D) {
		}

		void print(std::ostream& os) const {

			os << "[mu | sigma]: [" << mu.transpose() << " | " << sigma.format(CommaInitFmt) << "]";
		}
};

class Suffies_Dirichlet: public Suffies {
	public:
		double alpha;
//		Suffies base_suffies;
		
		void print(std::ostream& os) const {
			os << alpha;
		}
		
};
