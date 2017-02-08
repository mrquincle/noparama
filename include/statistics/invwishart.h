#pragma once

#include <Eigen/Dense>
#include <random>

#include <statistics/normal.h>

/*!
 * The base class distribution_t allows using this class, there were needs to be sampled from a probability 
 * distribution. The sufficient statistics are distribution dependent and are stored within this derived class.
 */
class inverse_wishart_distribution: public distribution_t {
	private:
		Suffies_InvWishart _suffies_iw;
		
		// Helper pdf
		Suffies_Normal _suffies_normal;
		int _D;

		// Result
		Suffies_ZeroCentered_MultivariateNormal *_suffies_result;
	public:
		inverse_wishart_distribution(Suffies_InvWishart & suffies_iw): 
			_suffies_iw(suffies_iw)
		{
			_distribution_type = InvWishart;

			_D = _suffies_iw.Lambda.cols();

			_suffies_normal.mu = _D;
			_suffies_normal.sigma = _suffies_iw.nu;
		}

		Suffies_ZeroCentered_MultivariateNormal *operator()(random_engine_t & generator) 
		{
			_suffies_result = new Suffies_ZeroCentered_MultivariateNormal(_D);

			normal_distribution normal_dist(_suffies_normal);

			Eigen::MatrixXd L( _suffies_iw.Lambda.llt().matrixL() );
			auto x = L.transpose() * normal_dist(generator)->mu;
			
			_suffies_result->sigma = x * x.transpose();
			
			return _suffies_result;
		}
};
