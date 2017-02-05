#pragma once

#include <np_suffies.h>

#include <np_data.h>

class distribution_t {
	protected:
		Suffies _suffies;
	public:
		//! Constructor in case sufficient statistics are not used
		distribution_t() {};

		//! Constructor that takes not further defined sufficient statistics
		distribution_t(Suffies & suffies): _suffies(suffies) {
		};

		virtual void init(Suffies & suffies) {
			_suffies = suffies;
		};

		//! Generate random sample that can in turn be sufficient statistics for another distribution
		template<typename _UniformRandomNumberGenerator >
		Suffies & sample(_UniformRandomNumberGenerator generator) {
			return _suffies;
		};

		virtual double probability(dataset_t & dataset) const {
			return 0;
		};
		
		virtual double probability(data_t & data) const {
			return 0;
		};
};
