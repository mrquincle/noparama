#pragma once

#include <Eigen/Dense>
#include <map>

/*!
 * The type of distribution ranges from a standard Normal to a Dirichlet distribution that can be used for 
 * nonparametric Bayesian models.
 */
enum distribution_type_t { 
	// gamma
	Gamma,
	// normals
	Normal, Unity_Normal,
	// multivariate normals
	MultivariateNormal, Unity_MultivariateNormal, ZeroCentered_MultivariateNormal, ScalarNoise_MultivariateNormal,
	// inverse wisharts and alike, useful as priors for normals
	NormalInvWishart, InvWishart, NormalInvGamma, 
	// hierarchical priors
	Dirichlet
};

static std::map< distribution_type_t, const char * > distribution_type_str = {
	{Gamma,                            "gamma"},
	{Normal,                           "normal"},
	{Unity_Normal,                     "unity normal"},
	{MultivariateNormal,               "multivariate normal"},
	{Unity_MultivariateNormal,         "unity multivariate normal"},
	{ZeroCentered_MultivariateNormal,  "zero-centered multivariate normal"},
	{ScalarNoise_MultivariateNormal,  "zero-centered multivariate normal"},
	{NormalInvWishart,                 "normal inverse-Wishart"},
	{InvWishart,                       "inverse-Wishart"},
	{NormalInvGamma,                   "normal inverse-Gamma"},
	{Dirichlet,                        "Dirichlet"}
};

/*! Concise way of representing a matrix on a single line.
 */
static Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", "; ", "", "", "", "");

/*!
 * Suffies describes a probability distribution. On the moment it does not contain any information on itself.
 *
 * Preferably the design should be such that there is not some kind of type field (to cast to derived classes). That
 * normally namely stinks of bad design. It might be okay to have a string or something like that for cosmetic
 * reasons (printing e.g.).
 */
class Suffies {
	public:
		//! dimension
		int D;

		virtual ~Suffies() {
		};

		virtual void print(std::ostream& os) const {
			os << D << " (probably a copy by value rather than copy by reference somewhere?)";
		}

		friend std::ostream& operator<<(std::ostream& os, const Suffies& s);  
};

class Suffies_Gamma: public Suffies {
	public:
		double alpha;
		double beta;
		
		void print(std::ostream& os) const {
			os << "g[alpha|beta]: " << alpha << " | " << beta;
		}
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
			os << "niw[nu|kappa|nu|Lambda]: " << nu << " | " << kappa << " | " << nu << " | " << Lambda.format(CommaInitFmt);
		}
};

class Suffies_InvWishart: public Suffies {
	public:
		double nu;
		Eigen::MatrixXd Lambda;
		
		Suffies_InvWishart(int D): Lambda(D,D) {
			Suffies::D = D;
		}
		
		void print(std::ostream& os) const {
			os << "iw[nu|Lambda]: " << nu << " | " << Lambda.format(CommaInitFmt);
		}
};

class Suffies_NormalInvGamma: public Suffies {
	public:
		Eigen::VectorXd mu;
		double alpha;
		double beta;
		//double gamma;
		Eigen::MatrixXd Lambda; // the multidimensional version
		
		Suffies_NormalInvGamma(int D): mu(D), Lambda(D,D) {
			Suffies::D = D;
		}
		
		void print(std::ostream& os) const {
			os << "nig[mu|alpha|beta|gamma]: " << mu << " | " << alpha << " | " << beta << " | " << Lambda.format(CommaInitFmt);
		}
};

class Suffies_Normal: public Suffies {
	public:
		double mu;
		double sigma;
		
		void print(std::ostream& os) const {
			os << "n[mu|sigma]: " << mu << " | " << sigma;
		}
};

class Suffies_Double: public Suffies {
	public:
		double val;

		void print(std::ostream& os) const {
			os << "[val]: " << val;
		}
};

/*!
 * A multivariate normal distribution with an identity matrix as covariance / information matrix.
 */
class Suffies_Unity_MultivariateNormal: public Suffies {
	public:
		Eigen::VectorXd mu;

		Suffies_Unity_MultivariateNormal(int D): mu(D) {
			Suffies::D = D;
		}
		
		void print(std::ostream& os) const {
			os << "mvn[mu]: " << mu;
		}
};

/*!
 * A multivariate normal disribution centered at 0.
 */
class Suffies_ZeroCentered_MultivariateNormal: public Suffies {
	public:
		Eigen::MatrixXd sigma;
		
		Suffies_ZeroCentered_MultivariateNormal(int D): sigma(D, D) {
			Suffies::D = D;
		}
		
		void print(std::ostream& os) const {
			os << "mvn[sigma]: " << sigma.format(CommaInitFmt);
		}
};

class Suffies_MultivariateNormal: public Suffies {
	public:
		Eigen::VectorXd mu;
		Eigen::MatrixXd sigma;
		
		Suffies_MultivariateNormal(int D): mu(D), sigma(D,D) {
			Suffies::D = D;
		}
		
		void print(std::ostream& os) const {
			os << "mvn[mu|sigma]: " << mu.transpose() << " | " << sigma.format(CommaInitFmt);
		}
};

class Suffies_ScalarNoise_MultivariateNormal: public Suffies {
	public:
		Eigen::VectorXd mu;
		double sigma;
		
		Suffies_ScalarNoise_MultivariateNormal(int D): mu(D)  {
			Suffies::D = D;
		}
		
		void print(std::ostream& os) const {
			os << "mvn[mu|sigma]: " << mu.transpose() << " | " << sigma;
		}
};

class Suffies_Dirichlet: public Suffies {
	public:
		double alpha;
//		Suffies base_suffies;
		
		void print(std::ostream& os) const {
			os << "d[alpha]: " << alpha;
		}
		
};
