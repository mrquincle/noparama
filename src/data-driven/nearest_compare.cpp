#include <data-driven/nearest_compare.h>

#include <data-driven/emd_meanshift.h>
#include <helper/dim1algebra.hpp>
#include <math.h>

//#define EXCESSIVE_DEBUG

nearest_compare::nearest_compare(
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
	
	// map vector with <data> objects to an array
	get_raw(_dataset_reference, _dataset_reference_raw);
	// calculate the mean 
	calc_mean(_dataset_reference_raw, size, dim, _dataset_reference_mean);

	// create data compare object of right size
	int mu_size = _suffies_mvn->mu.size();
	_data_compare = new data_t(mu_size);
}

void nearest_compare::calc_mean(float * data, int size, int dim, float * result) const {
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

void nearest_compare::get_shape(const dataset_t & dataset, int & size, int & dim) const {
	size = (int)dataset.size();

	if (size < 1) {
		std::cerr << "Error: there must be a dataset as reference" << std::endl;
		return;
	}
	
	data_t *data = dataset[0];
	dim = (int)data->size();
}

/**
 * Create an array of size N * dim from N data items. Each data item occupies successive values in the array.
 */
void nearest_compare::get_raw(dataset_t & source, float* target) const {
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

void nearest_compare::init(Suffies & suffies) {
	if (typeid(suffies)!=typeid(*_suffies_mvn)) {
		std::cout << "init nearest compare" << std::endl;
		std::cout << "current suffies: " << *_suffies_mvn << std::endl;
		std::cout << "incoming suffies: " << suffies << std::endl;
		assert (typeid(suffies)==typeid(*_suffies_mvn));
	}
	_suffies_mvn = &dynamic_cast<Suffies_ScalarNoise_MultivariateNormal&>(suffies);
}

Suffies & nearest_compare::getSuffies() {
	return *_suffies_mvn;
}

/**
 * Get nearest neighbour.
 */
double nearest_compare::get_nearest_neighbour(data_t & test_data) const {

	assert(_dataset_reference.size() > 0);

	double min_dist;

	for (unsigned int i = 0; i < _dataset_reference.size(); ++i) {
		auto & data = *_dataset_reference[i];

		double dist = algebra::distance<double>(data.begin(), data.end(), test_data.begin(), test_data.end(), 
				algebra::DM_EUCLIDEAN); 

		if (i == 0 || dist < min_dist) {
			min_dist = dist;
		}
	}

	return min_dist;
}

/**
 * We can later use something like kd-tree, but for now just find the closest data point in the reference set and 
 * calculate the probability using a particular distance metric. For example, assume errors with respect to the
 * reference are distribution according to a normal distribution, or just calculate Euclidean distance.
 */
double nearest_compare::probability(data_t & data) const
{
	int D = data.size();
	double *mu_data = _suffies_mvn->mu.data();
	for (int i = 0; i < D; ++i) {
		(*_data_compare)[i] = data[i] - mu_data[i];
	}
	double dist = get_nearest_neighbour(*_data_compare);

	double prob = 1.0/sqrt(2.0*M_PI)*exp(-1/2.0*dist*dist);

#ifdef EXCESSIVE_DEBUG
	algebra::print(data.begin(), data.end());
	std::cout << "dist: " << dist << ", prob: " << prob << std::endl;
#endif
	if (prob < 0) {
		prob = 0.0;
	}
	return prob;
}

double nearest_compare::probability(dataset_t & dataset) const
{
	assert (dataset.size() > 0);

	double result = 1.0;
	for (auto data: dataset) {
		result *= probability(*data);
	}
	return result;
}

double nearest_compare::logprobability(data_t & data) const
{
	int D = data.size();
	double *mu_data = _suffies_mvn->mu.data();
	for (int i = 0; i < D; ++i) {
		(*_data_compare)[i] = data[i] - mu_data[i];
	}
	double dist = get_nearest_neighbour(*_data_compare);

	double prob = -1/2.0 * log(2.0*M_PI) + -1/2.0*dist*dist;
#ifdef EXCESSIVE_DEBUG
	algebra::print(data.begin(), data.end());
	std::cout << "dist: " << dist << ", prob: " << prob << std::endl;
#endif

	return prob;
}

double nearest_compare::logprobability(dataset_t & dataset) const
{
	assert (dataset.size() > 0);

	double result = 0.0;
	for (auto data: dataset) {
		result += logprobability(*data);
	}
#ifdef EXCESSIVE_DEBUG
	std::cout << "------------------------------------" << std::endl;
	std::cout << "total result: " << result << std::endl;
	std::cout << "------------------------------------" << std::endl;
#endif
	return result;
}
