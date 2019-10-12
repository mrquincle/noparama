#include <data-driven/set_compare.h>

#include <data-driven/emd_meanshift.h>

#include <math.h>

#define SET_COMPARE_USE_PRIOR

//#define EXCESSIVE_DEBUG

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
 * 4. Call the emd_meanshift function or the approxmatch_cpu function (without offset calculations).
 * 5. Calculate the costs.
 *
 * Mmmm.... It does seem to matter though... Yes, it does not return a probability. Moreover, it does not return
 * a unnormalized probability either. The probability of a dataset is not formed by a product of the probabilities of
 * individual data items. They are conditionally dependent. This means that a split action - even if the fit does not
 * change - leads to different results. Henceforth, there is an artifact for the number of items in a cluster and
 * the system will either converge to very few or diverge to very many clusters (depending on this artifact).
 *
 * If we use 1-nearest neighbour to calculate the probability and we make sure that the probability of a single data
 * item is strictly less than one. Then we can calculate probabilities of sets by just multiplication.
 *
 * It can happen that in this case individual lines (or groups of points) will be shifted around with mean parameters
 * that make sense for that subset. This is due to the fact that the conditional dependencies of points within a 
 * particular object is lost. The points can map to the same point on an object. The probability of a point matching
 * another point does not influence the probability for the "next" point fitting that object. That there will be only
 * a few objects (lines, features, whatever), used comes now only from the (sparse) Dirichlet prior.
 */
double set_compare::probability(dataset_t & dataset) const
{
	assert (dataset.size() > 0);

	int b = 1;

	int n = 0, tmp_dim = 0;
	get_shape(_dataset_reference, n, tmp_dim);
	
	int m = 0;
	get_shape(dataset, m, tmp_dim);
	assert(dim == tmp_dim);
	get_raw(dataset, _dataset_compare_raw);
#ifdef SET_COMPARE_USE_PRIOR
	double *mu = _suffies_mvn->mu.data();
	int mu_size = _suffies_mvn->mu.size();
	if (mu_size != tmp_dim) {
		std::cerr << "Error: tmp_dimensions do not match: " << mu_size << " vs " << tmp_dim << std::endl;
		assert (mu_size == tmp_dim);
	}
	for (int i = 0; i < mu_size; ++i) {
		_dataset_compare_mean[i] = ((float)mu[i]);
	}
	/*
	std::cout << "Mean: ";
	std::string sep = "";
	for (int i = 0; i < mu_size; ++i, sep = ", ") {
		std::cout << sep << _dataset_compare_mean[i];
	}
	std::cout << std::endl;
	*/
#else
	calc_mean(_dataset_compare_raw, m, tmp_dim, _dataset_compare_mean);
#endif
	
#ifdef EXCESSIVE_DEBUG
	std::cout << "Mean: " << _dataset_compare_mean[1] << std::endl;
#endif
	// if we pick the right mean the match matrix becomes almost singular, if very far off, it becomes uniformly distributed (it doesn't matter how to move)
//	_dataset_compare_mean[1] = 10;

	emd_global_offset(b, n, m, _dataset_reference_raw, _dataset_compare_raw, _dataset_match,
			_dataset_reference_mean, _dataset_compare_mean);

#ifdef EXCESSIVE_DEBUG
	std::cout << "Match: " << std::endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; ++j) {
			std::cout << _dataset_match[i*m+j] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
#endif
	float costs = 0.0;
	emd_mean_costs_global_offset(b, n, m, _dataset_reference_raw, _dataset_compare_raw, _dataset_match, 
			_dataset_reference_mean, _dataset_compare_mean, &costs);

#ifdef EXCESSIVE_DEBUG
	std::cout << "Costs: " << costs << std::endl;
#endif

	/*
	 * These are ad-hoc adjustments where we adjust the costs with the number of data points. We can tweak it a
	 * bit. For example, consider n the number of points in the reference dataset and m the number of points in
	 * the test set, then scale the costs by (m/n)^2. In other words, scale by 1 if m == n. However, if there
	 * are fewer points in the test set, say half of the reference set, scale by (1/2)^2. 
	 *
	 * This does not work however. Let us calculate the costs of X given R, where R is the reference set and
	 * X is half its cardinality. Then these costs are not just half of the original. If the points are 
	 * uniformly removed, what do we actually expect?
	 */
//	scale with number of datapoints, many data points are cons
	const float scale_factor = m/(float)n;
	costs = costs * pow(scale_factor, 2);
	//const float scale_factor = 1.0;
	//costs = costs * scale_factor;
#ifdef EXCESSIVE_DEBUG
	std::cout << "Scaled costs: " << costs << std::endl;
#endif
	/*
	 * Here we map the cost (which is a nonnegative scalar) to a value between 0 and 1. A cost of zero will be
	 * mapped to a likelihood of one. A cost of e.g. 10 will be mapped to a likelihood of exp(-10) which is 
	 * very small: 45e-6.
	 */
	double result = (exp(-costs) + 10e-9);
	return result;
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
