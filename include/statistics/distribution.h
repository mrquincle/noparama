#pragma once

#include <iostream>
#include <random>
#include <memory>
#include <assert.h>

#include "np_suffies.h"
#include "np_data.h"

typedef std::default_random_engine random_engine_t;

#define TEST_IF_ALL_VIRTUALS_ARE_IMPLEMENTED     0

#define INIT_SHOULD_BE_IMPLEMENTED 1
#define GET_SUFFIES_SHOULD_BE_IMPLEMENTED 1
#define SAMPLE_OPERATOR_SHOULD_BE_IMPLEMENTED 1
#define PROBABILITY_FOR_DATASET_SHOULD_BE_IMPLEMENTED 1
#define PROBABILITY_FOR_DATAPOINT_SHOULD_BE_IMPLEMENTED 1
#define LOGPROBABILITY_FOR_DATASET_SHOULD_BE_IMPLEMENTED 1
#define LOGPROBABILITY_FOR_DATAPOINT_SHOULD_BE_IMPLEMENTED 1

class distribution_t {
	protected:
		distribution_type_t _distribution_type;
			
	public:
		//! Constructor in case the sufficient statistics are not used
		distribution_t() {};

		virtual ~distribution_t() {};

#if TEST_IF_ALL_VIRTUALS_ARE_IMPLEMENTED==1
		virtual void init(Suffies & suffies) = 0;

		virtual  Suffies* operator()(random_engine_t & generator) = 0;
		
		virtual double probability(dataset_t & dataset) const = 0;
		
		virtual double probability(data_t & data) const = 0;
		
		virtual double logprobability(dataset_t & dataset) const = 0;
		
		virtual double logprobability(data_t & data) const = 0;

#else
		virtual void init(Suffies & suffies) {
			std::cout << "Distribution: " << distribution_type_str[_distribution_type] << std::endl;
			assert(INIT_SHOULD_BE_IMPLEMENTED==0);
		};

		virtual Suffies & getSuffies() {
			std::cout << "Distribution: " << distribution_type_str[_distribution_type] << std::endl;
			assert(GET_SUFFIES_SHOULD_BE_IMPLEMENTED==0);
		}

		//! Generate random sample that can in turn be sufficient statistics for another distribution
		//virtual std::unique_ptr<Suffies> operator()(random_engine_t & generator) {
		virtual Suffies* operator()(random_engine_t & generator) {
			std::cout << "Distribution: " << distribution_type_str[_distribution_type] << std::endl;
			assert(SAMPLE_OPERATOR_SHOULD_BE_IMPLEMENTED==0);
			return NULL;
		}; 

		virtual double probability(dataset_t & dataset) const {
			std::cout << "Distribution: " << distribution_type_str[_distribution_type] << std::endl;
			assert(PROBABILITY_FOR_DATASET_SHOULD_BE_IMPLEMENTED==0);
			return 0;
		};
		
		virtual double probability(data_t & datum) const {
			std::cout << "Distribution: " << distribution_type_str[_distribution_type] << std::endl;
			assert(PROBABILITY_FOR_DATAPOINT_SHOULD_BE_IMPLEMENTED==0);
			return 0;
		}
		;
		virtual double logprobability(dataset_t & dataset) const {
			std::cout << "Distribution: " << distribution_type_str[_distribution_type] << std::endl;
			assert(LOGPROBABILITY_FOR_DATASET_SHOULD_BE_IMPLEMENTED==0);
			return 0;
		};
		
		virtual double logprobability(data_t & datum) const {
			std::cout << "Distribution: " << distribution_type_str[_distribution_type] << std::endl;
			assert(LOGPROBABILITY_FOR_DATAPOINT_SHOULD_BE_IMPLEMENTED==0);
			return 0;
		};
#endif
};
