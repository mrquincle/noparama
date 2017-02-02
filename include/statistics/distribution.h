#pragma once

class distribution_t {


	public:
		//! Constructor in case sufficient statistics are not used
		distribution_t();

		//! Constructor that takes not further defined sufficient statistics
		distribution_t(Suffies suffies);

		//! Generate random sample that can in turn be sufficient statistics for another distribution
		template<typename _UniformRandomNumberGenerator >
		Suffies sample() const;
};
