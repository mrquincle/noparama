#include <data-driven/set_compare.h>

#include <data-driven/emd_meanshift.h>

#define SET_COMPARE_USE_PRIOR

set_compare::set_compare(
		Suffies_ScalarNoise_MultivariateNormal & suffies, 
		dataset_t & dataset_reference)
{
	_suffies_mvn = &suffies;
	_distribution_type = Unity_MultivariateNormal;
	_dataset_reference = dataset_reference;

	int size = 0, dim = 0;
	get_shape(dataset_reference, size, dim);
	_dataset_reference_raw = (float*)malloc(size * dim * sizeof(float));
	_dataset_reference_mean = (float*)malloc(dim * sizeof(float));
	// we assume the reference set has more points, so, then we do not need to allocate all the time
	_dataset_compare_raw = (float*)malloc(size * dim * sizeof(float));
	_dataset_compare_mean = (float*)malloc(dim * sizeof(float));

	_dataset_match = (float*)malloc(size * size * sizeof(float));

	get_raw(_dataset_reference, _dataset_reference_raw);
	calc_mean(_dataset_reference_raw, size, dim, _dataset_reference_mean);
}

void set_compare::calc_mean(float * data, int size, int dim, float * result) const {
	if (result == NULL) {
		std::cerr << "Error: result not allocated" << std::endl;
		assert(false);
	}
	double tmp[dim];
	for (int j = 0; j < dim; ++j) {
		tmp[j] = 0;
	}
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < dim; ++j) {
			tmp[j] += data[i*dim + j];
		}
	}
	for (int j = 0; j < dim; ++j) {
		result[j] = tmp[j]/(float)size;
	}
}

void set_compare::get_shape(const dataset_t & dataset, int & size, int & dim) const {
	size = (int)dataset.size();

	if (size < 1) {
		std::cerr << "Error: there must be a dataset as reference" << std::endl;
		return;
	}
	
	data_t *data = dataset[0];
	dim = (int)data->size();
}

void set_compare::get_raw(dataset_t & source, float* target) const {
	int N = 0, dim = 0;
	get_shape(source, N, dim);

	// fill matrix
	for (int i = 0; i < N; ++i) {
		data_t & data = *source[i];
		for (int j = 0; j < dim; ++j) {
			target[i*dim+j] = (float)data[j];
		}
	}
}

void set_compare::init(Suffies & suffies) {
	if (typeid(suffies)!=typeid(*_suffies_mvn)) {
		std::cout << "init mvn" << std::endl;
		std::cout << "current suffies: " << *_suffies_mvn << std::endl;
		std::cout << "incoming suffies: " << suffies << std::endl;
		assert (typeid(suffies)==typeid(*_suffies_mvn));
	}
	_suffies_mvn = &dynamic_cast<Suffies_ScalarNoise_MultivariateNormal&>(suffies);
}

Suffies & set_compare::getSuffies() {
	return *_suffies_mvn;
}

/*
 * Suffies_ScalarNoise_MultivariateNormal* set_compare::operator()(random_engine_t &generator) 
{

}
*/

double set_compare::probability(data_t & data) const
{
	// do not use
	assert(false);
}

/**
 * This should actually return a true probability. This class is however used in a likelihood setting. If all 
 * results are scaled by the same amount it doesn't matter for the MCMC algorithm.
 *
 * 1. Map original dataset to plain C 2D array and store as class variable (in prepare()).
 * 2. Calculate mean over incoming dataset 
 * 3. Map the given dataset to plain 2D array.
 * 4. Call the emd_meanshift function (or the approxmatch_cpu function (without offset calculations).
 * 5. Calculate the costs.
 */
double set_compare::probability(dataset_t & dataset) const
{
	assert (dataset.size() > 0);

	int b = 1;

	int n = 0, dim = 0;
	get_shape(_dataset_reference, n, dim);
	
	int m = 0;
	get_shape(dataset, m, dim);
	get_raw(dataset, _dataset_compare_raw);
#ifdef SET_COMPARE_USE_PRIOR
	double *mu = _suffies_mvn->mu.data();
	int mu_size = _suffies_mvn->mu.size();
	if (mu_size != dim) {
		std::cerr << "Error: dimensions do not match: " << mu_size << " vs " << dim << std::endl;
		assert (mu_size == dim);
	}
	for (int i = 0; i < mu_size; ++i) {
		_dataset_compare_mean[i] = (float)mu[i];
	}
#else
	calc_mean(_dataset_compare_raw, m, dim, _dataset_compare_mean);
#endif

	emd_global_offset(b, n, m, _dataset_reference_raw, _dataset_compare_raw, _dataset_match,
			_dataset_reference_mean, _dataset_compare_mean);

	float costs = 0.0;
	emd_mean_costs_global_offset(b, n, m, _dataset_reference_raw, _dataset_compare_raw, _dataset_match, 
			_dataset_reference_mean, _dataset_compare_mean, &costs);

	return (double)costs;
}

double set_compare::logprobability(data_t & data) const
{
	assert(false);
}

double set_compare::logprobability(dataset_t & dataset) const
{
	double prob = probability(dataset);
	assert(prob > 0);
	return std::log(prob);
}
