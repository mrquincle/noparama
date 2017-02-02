#pragma once

#include <np_suffies.h>

#include <np_data.h>

class distribution_t {
	public:
		//! Constructor in case sufficient statistics are not used
		distribution_t();

		//! Constructor that takes not further defined sufficient statistics
		distribution_t(Suffies & suffies);

		virtual void init(Suffies & suffies);

		//! Generate random sample that can in turn be sufficient statistics for another distribution
		template<typename _UniformRandomNumberGenerator >
		Suffies & sample(_UniformRandomNumberGenerator generator);

		virtual double probability(Eigen::VectorXd & value) const;
		
		virtual double probability(dataset_t & dataset) const;
};
