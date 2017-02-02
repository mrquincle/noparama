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
		
		// Result
		Suffies_ZeroCentered_MultivariateNormal _suffies_result;
	public:
		inverse_wishart_distribution(Suffies_InvWishart & suffies_iw): 
			_suffies_iw(suffies_iw)
		{
			int D = _suffies_iw.Lambda.cols();

			_suffies_normal.mu = D;
			_suffies_normal.sigma = _suffies_iw.nu;
		}

		template<typename _UniformRandomNumberGenerator >
		Suffies_ZeroCentered_MultivariateNormal & operator()(_UniformRandomNumberGenerator & generator) const
		{
			normal_distribution normal_dist(_suffies_normal);

			Eigen::MatrixXd L( _suffies_iw.Lambda.llt().matrixL() );
			auto x = L.transpose() * normal_dist(generator);
			
			_suffies_result.sigma = x * x.transpose();
			
			return _suffies_result;
		}
};
