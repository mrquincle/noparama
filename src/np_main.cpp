#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <unordered_map>

// Include getopt
#include <unistd.h>
#include <experimental/filesystem>

#include <np_mcmc.h>
#include <np_data.h>
#include <np_results.h>

#include <membertrix.h>

#include <statistics/scalarnoise_multivariatenormal.h>
#include <statistics/multivariatenormal.h>
#include <statistics/dirichlet.h>
#include <statistics/normalinvwishart.h>
#include <statistics/normalinvgamma.h>

#include <data-driven/set_compare.h>
#include <data-driven/nearest_compare.h>

//#include <np_neal_algorithm2.h>
#include <np_neal_algorithm8.h>
#include <np_jain_neal_algorithm.h>
#include <np_triadic_algorithm.h>

#include <pretty_print.hpp>


#include <dim1algebra.hpp>

using namespace std;
using namespace algebra;

namespace fs = std::experimental::filesystem;

enum algorithm_t { /*algorithm2, */ algorithm8, jain_neal_split, triadic };

// forward declarations of test functions
void test_likelihood_function(std::default_random_engine &generator, dataset_t &dataset, ground_truth_t & ground_truth, dataset_t &ref_dataset, distribution_t & likelihood, dirichlet_process &hyper);

void disp_help(std::string appname) {
	std::string datafilename;
	std::cout << "noparama [version 0.1.72]" << std::endl;
	std::cout << std::endl;
	std::cout << "Usage: " << std::endl;
	std::cout << "  noparama <-d datafile> <-a algorithm8|jain_neal_split|triadic> <-T ticks> <-c regression|clustering|angular|points3d> <-r reference file>" << std::endl;
	std::cout << "For example:" << std::endl;
	datafilename = "datasets/twogaussians.data";
	std::cout << "  " << appname << " -d " << datafilename << " -a triadic -T 5000 -c regression" << std::endl;
}

void read_data(std::string & datafilename, representation_mode_t & representation_mode, ground_truth_t * ground_truth, dataset_t & dataset) {
	char _verbosity = Information;
	fout << "Load file: " << datafilename << std::endl;
	if (ground_truth) {
		ground_truth->clear();
	}
	dataset.clear();
	
	// potentially limit number of data items (for debug purposes)
	bool limit = false;
	int N = 16;
	if (limit) {
		fout << "Using limited dataset with only " << N << " items " << std::endl;
	}

	// The data file should have "a b c" on lines, separated by spaces (without quotes, each value of the type double).
	fout << "Read dataset" << std::endl;
	std::ifstream datafilehandle(datafilename);
	std::string line;
	double a, b, c, d;
	int n = 0;
	while (std::getline(datafilehandle, line)) {
		std::istringstream ss(line);
		ss >> a >> b >> c;
		switch(representation_mode) {
		case regression_mode:
			{
				data_t *data = new data_t(3);
				// prepend vector with constant for regression
				*data = { 1, a, b };
				dataset.push_back(data);
				if (ground_truth) {
					(*ground_truth)[n] = (int)c;
				}
				break;
			}
		case angular_mode: case clustering_mode: 
			{
				data_t *data = new data_t(2);
				*data = { a, b };
				dataset.push_back(data);
				if (ground_truth) {
					(*ground_truth)[n] = (int)c;
				}
				break;
			}
		case points3d_mode: 
			{
				data_t *data = new data_t(3);
				*data = { a, b, c};
				dataset.push_back(data);
				if (ground_truth) {
					if (ss >> d) {
						(*ground_truth)[n] = (int)d;
					} else {
						std::cerr << "Not enough columns" << std::endl;
					}
				}
			}
		}
		n++;
		if (limit && (n == N)) break;
	}

	if (n == 0 || dataset.empty()) {
		fout << "No data found... Check the file or the contents of the file." << std::endl;
		exit(7);
	}

	int I = 5;
	fout << "Display first and last " << I << " items of the dataset" << std::endl;
	if (ground_truth) {
		fout << "Include ground truth" << endl;
		for (int i = 0; i < I; ++i) {
			if (i == (int)ground_truth->size()) break;
			fout << "Data: " << *dataset[i] << " with ground truth " << (*ground_truth)[i] << std::endl;
		}
		for (int i = (int)dataset.size() - I; i < (int)dataset.size(); ++i) {
			if (i == (int)ground_truth->size()) break;
			fout << "Data: " << *dataset[i] << " with ground truth " << (*ground_truth)[i] << std::endl;
		}
	} else {
		fout << "No ground truth available" << endl;
		for (int i = 0; i < I; ++i) {
			fout << "Data: " << *dataset[i] << std::endl;
		}
		for (int i = dataset.size() - I; i < (int)dataset.size(); ++i) {
			fout << "Data: " << *dataset[i] << std::endl;
		}
	}

}

/**
 * With the dataset twogaussians it is important to realize that this did not originate from a Dirichlet process. Even
 * if alpha is chosen to be very low (say 0.0001), there will be multiple clusters generated, not just two. The purity
 * will be quite high (very few misclassifications) but the specificity is low (identification of multiple clusters
 * where there is only one).
 */
int main(int argc, char *argv[]) {
	char _verbosity = Debug;
	
	fout << "Welcome to noparama" << std::endl;

	// ---------------------------------------------------------------------------------------------------------------
	// Configuration parameters 
	// ---------------------------------------------------------------------------------------------------------------
	double alpha = 1;

	bool subsample = true;
	int subsample_size = 200;
//	subsample_size = 128;
//	subsample_size = 16;

	// If true, writes samples from (Dirichlet) prior to file fhyper.txt and quits 
	bool test_prior = false;

	bool test_likelihood = false;

	// ---------------------------------------------------------------------------------------------------------------

	fout << "The Dirichlet Process is run with alpha=" << alpha << std::endl;

	std::default_random_engine generator(std::random_device{}()); 

	// Configuration std::strings that can be filled from CLI arguments
	std::string datafilename;
	std::string reffilename;
	int token;
	std::unordered_map<int, std::string> flags;
	while( (token = getopt(argc, argv, "d:c:a:T:h?r:")) != EOF)
	{
		switch (token)
		{
			case 'd':
				flags['d'] = std::string(optarg);
				break;
			case 'c':
				flags['c'] = std::string(optarg);
				break;
			case 'a':
				flags['a'] = std::string(optarg);
				break;
			case 'T':
				flags['T'] = std::string(optarg);
				break;
			case 'h': case '?':
				disp_help(std::string(argv[0]));
				exit(1);
				break;
			case 'r':
				flags['r'] = std::string(optarg);
				break;
		}
	}

	// required, so show help and exit
	std::string algorithm_str = flags['a'];
	algorithm_t algorithm;
	if (flags.find('d') == flags.end() || flags.find('a') == flags.end()) {
		disp_help(std::string(argv[0]));
		exit(1);
	} else {
		datafilename = flags['d'];

		/*if (algorithm_str == "algorithm2") {
			algorithm = algorithm2;
		} else */
		if (algorithm_str == "algorithm8") {
			algorithm = algorithm8;
		} else if (algorithm_str == "jain_neal_split") {
			algorithm = jain_neal_split;
		} else if (algorithm_str == "triadic") {
			algorithm = triadic;
		} else {
			std::cerr << "Unknown algorithm: " << algorithm_str << std::endl;
			exit(1);
		}
	}

	if (flags.find('r') != flags.end()) {
		reffilename = flags['r'];
	}
	
	// optional T flag
	int T = 2000;
	if (flags.find('T') != flags.end()) {
		T = std::stoi(flags['T']);
	}

	// optional clustering/regression flag
	representation_mode_t representation_mode = regression_mode; 
	if (flags.find('c') != flags.end()) {
		std::string configuration = flags['c'];
		if (configuration == "regression") {
			fout << "Regression mode" << std::endl;
			representation_mode = regression_mode;
		} else if (configuration == "clustering") {
			fout << "Clustering mode" << std::endl;
			representation_mode = clustering_mode;
		} else if (configuration == "angular") {
			fout << "Angular mode (form of regression)" << std::endl;
			representation_mode = angular_mode;
		} else if (configuration == "points3d") {
			fout << "Likelihood comparison (using 3D point clouds)" << std::endl;
			representation_mode = points3d_mode;
		} else {
			fout << "Unknown mode: " << configuration << std::endl;
			fout << "Picked regression mode by default." << std::endl;
		}
	}

	// check if output directory exists and bail out if this is the case
	std::string justdatafilename = fs::path(datafilename).filename();
	std::string workspace = "output/" + algorithm_str + '/' + justdatafilename + '/';
	fout << "Have workspace: " << workspace << std::endl;
	if (fs::exists(workspace)) {
		std::cerr << "Directory already exists. Drop out. We don't want to overwrite it or do double work." << std::endl;
		exit(106);
	}
	
	fout << "Run MCMC sampler for " << T << " steps " << std::endl;

	// get dataset and ground truth
	dataset_t dataset;
	ground_truth_t ground_truth;
	if (subsample) {
		dataset_t complete_dataset;
		ground_truth_t complete_ground_truth;
		read_data(datafilename, representation_mode, &complete_ground_truth, complete_dataset);
		std::vector<int> indices;
		indices.resize(complete_dataset.size());
		fill_successively(indices.begin(), indices.end());
		random_order(indices.begin(), indices.end());
		assert (subsample_size <= (int)indices.size());
		for (int i = 0; i < subsample_size; ++i) {
			dataset.push_back(complete_dataset[indices[i]]);
			ground_truth[i] = complete_ground_truth[indices[i]];
		}
	} else {
		read_data(datafilename, representation_mode, &ground_truth, dataset);
	}

	// get reference dataset
	dataset_t ref_dataset;
	if (reffilename != "") {
		if (subsample) {
			dataset_t ref_complete_dataset;
			read_data(reffilename, representation_mode, NULL, ref_complete_dataset);
			std::vector<int> indices;
			indices.resize(ref_complete_dataset.size());
			fill_successively(indices.begin(), indices.end());
			random_order(indices.begin(), indices.end());
			assert (subsample_size <= (int)indices.size());
			for (int i = 0; i < subsample_size; ++i) {
				ref_dataset.push_back(ref_complete_dataset[indices[i]]);
			}
		} else {
			read_data(reffilename, representation_mode, NULL, ref_dataset);
			assert (dataset.size() == ref_dataset.size());
		}
	}

	// The likelihood function only is required w.r.t. type, parameters are normally set later
	distribution_t *likelihood;
	if (representation_mode == regression_mode || representation_mode == angular_mode) {
		fout << "Multivariate normal suffies (sigma scalar)" << std::endl;
		// note that we need to allocate the suffies or else it won't exists after the call to new multivar..
		Suffies_ScalarNoise_MultivariateNormal * suffies_mvn = new Suffies_ScalarNoise_MultivariateNormal(2);
		suffies_mvn->mu << 0, 0;
		suffies_mvn->sigma = 1;
		likelihood = new scalarnoise_multivariate_normal_distribution(*suffies_mvn, representation_mode);
	} else if (representation_mode == clustering_mode) {
		fout << "Multivariate normal suffies" << std::endl;
		Suffies_MultivariateNormal * suffies_mvn = new Suffies_MultivariateNormal(2);
		suffies_mvn->mu << 0, 0;
		suffies_mvn->sigma << 1, 0, 0, 1; // get rid of this
		likelihood = new multivariate_normal_distribution(*suffies_mvn);
	} else if (representation_mode == points3d_mode) {
		fout << "Set comparison 'suffies'" << std::endl;
		Suffies_ScalarNoise_MultivariateNormal * suffies_mvn = new Suffies_ScalarNoise_MultivariateNormal(3);
		suffies_mvn->mu << 0, 0, 0;
		suffies_mvn->sigma = 1;
		//Suffies_Unity_MultivariateNormal * suffies_umvn = new Suffies_Unity_MultivariateNormal(2);
		//suffies_umvn->mu << 0, 0;
		//likelihood = new set_compare(*suffies_mvn, ref_dataset);
		likelihood = new nearest_compare(*suffies_mvn, ref_dataset);
	} else {
		std::cerr << "Unknown likelihood" << std::endl;
		exit(107);
	}
		
	fout << "likelihood still: " << likelihood->getSuffies() << std::endl;

	// The hierarchical prior has an alpha of 1 and the base distribution is handed separately through the prior
	fout << "Dirichlet distribution" << std::endl;
	Suffies_Dirichlet suffies_dirichlet;
	suffies_dirichlet.alpha = alpha;
	distribution_t *prior;

	if (representation_mode == regression_mode || representation_mode == angular_mode) {
		fout << "Normal Inverse Gamma distribution" << std::endl;
		Suffies_NormalInvGamma * suffies_nig = new Suffies_NormalInvGamma(2);
		suffies_nig->mu << 0, 0;
		suffies_nig->alpha = 10;
		suffies_nig->beta = 0.1;
		suffies_nig->Lambda << 0.01, 0, 0, 0.01;
		prior = new normal_inverse_gamma_distribution(*suffies_nig);
	} else if (representation_mode == clustering_mode) {
		fout << "Normal Inverse Wishart distribution" << std::endl;
		Suffies_NormalInvWishart * suffies_niw = new Suffies_NormalInvWishart(2);
		suffies_niw->mu << 6, 6;
		suffies_niw->kappa = 1.0/500;
		suffies_niw->nu = 4;
		suffies_niw->Lambda << 0.01, 0, 0, 0.01;
		prior = new normal_inverse_wishart_distribution(*suffies_niw);
	} else if (representation_mode == points3d_mode) {
		// Check if indeed appropriate
		Suffies_NormalInvGamma * suffies_nig = new Suffies_NormalInvGamma(3);
		suffies_nig->mu << 0, 0, 0;
		// if alpha very small, dispersion gets really big
		suffies_nig->alpha = 10;
		// with beta at 0.1 most values are at around +/-5 till +/-15, with beta at 1, there's +/-0 up to +/- 5
		suffies_nig->beta = 1;
		suffies_nig->Lambda << 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01;
		prior = new normal_inverse_gamma_distribution(*suffies_nig);
	} else {
		std::cerr << "Unknown likelihood" << std::endl;
		exit(107);
	}
	fout << "Create Dirichlet Process" << std::endl;
	dirichlet_process hyper(suffies_dirichlet, *prior);

	fout << "Set up InitClusters object" << std::endl;
	InitClusters init_clusters(generator, hyper);

	// print hyper 
	if (test_prior) {
		fout << "Write everything to file fhyper.txt" << endl;
		int kTest = 1000;
		ofstream fhyper;
		fhyper.open("fhyper.txt");
		for (int k = 0; k < kTest; ++k) {
			// sample sufficient statistics from nonparametrics
			Suffies_MultivariateNormal* suffies = hyper.sample_base(generator);
			double * array = suffies->mu.data();
			string sep = "";
			for (int i = 0; i < suffies->mu.size(); ++i, sep = ", ") {
				fhyper << sep << array[i];
			}
			fhyper << endl;
		}
		fhyper.close();
		exit(0);
	}

	// To update cluster parameters we need prior (hyper parameters) and likelihood, the membership matrix is left 
	// invariant
	fout << "Set up UpdateClusters object" << std::endl;
	UpdateClusters update_clusters(generator, *likelihood, hyper);

	// To update the cluster population we need	to sample new clusters using hyper parameters and adjust existing ones
	// using prior and likelihood
	fout << "Set up UpdateClusterPopulation object" << std::endl;

	int subset_count;
	UpdateClusterPopulation *update_cluster_population;
	switch(algorithm) {
/*		case algorithm2: 
			{
				fout << "We will be using algorithm 2 by Neal" << std::endl;
				update_cluster_population = new NealAlgorithm2(generator, *likelihood, hyper);
				subset_count = 1;
				break;
			}
*/
		case algorithm8: 
			{
				fout << "We will be using algorithm 8 by Neal" << std::endl;
				update_cluster_population = new NealAlgorithm8(generator, *likelihood, hyper);
				subset_count = 1;
				break;
			}
		case jain_neal_split:
			{
				fout << "We will be using the Jain-Neal algorithm" << std::endl;
				update_cluster_population = new JainNealAlgorithm(generator, *likelihood, hyper);
				subset_count = 2;
				break;
			}
		case triadic:
			{
				fout << "We will be using the Triadic algorithm" << std::endl;
				update_cluster_population = new TriadicAlgorithm(generator, *likelihood, hyper);
//#define TEST_WITH_TWO
#ifdef TEST_WITH_TWO
				subset_count = 2;
#else
				subset_count = 3;
#endif
				break;
			}
	}

	if (test_likelihood) {
		test_likelihood_function(generator, dataset, ground_truth, ref_dataset, *likelihood, hyper);
		exit(0);
	}

	// create MCMC object
	fout << "Set up MCMC" << std::endl;
	MCMC & mcmc = *new MCMC(generator, init_clusters, update_clusters, *update_cluster_population, subset_count, *likelihood);

	fout << "Run MCMC for " << T << " steps" << std::endl;
	mcmc.run(dataset, T);

	fout << "Print statistics" << std::endl;
	update_cluster_population->printStatistics();
	
	// analyse and write out results
	std::chrono::time_point<std::chrono::system_clock> clock;
	clock = std::chrono::system_clock::now();
	std::time_t time = std::chrono::system_clock::to_time_t(clock);
	std::tm tm = *std::localtime(&time);

	std::stringstream tss; 
	tss << std::put_time(&tm, "%Y%m%d_%H:%M");
	std::string dirname = tss.str(); 
	std::string basename = "snapshot";
	
	fout << "Last item" << std::endl;
	const membertrix trix_snapshot = mcmc.getMembershipMatrix();
	Results results_snapshot(trix_snapshot, ground_truth);
	results_snapshot.write(workspace, dirname, basename);
	
	basename = "results";
	
	fout << "Write results" << std::endl;
	const membertrix trix = mcmc.getMaxLikelihoodMatrix();
	Results results(trix, ground_truth);
	results.write(workspace, dirname, basename);
	
	fout << "Memory deallocation" << std::endl;
	delete &mcmc;
	delete likelihood;
	delete prior;
	for (int i = 0; i < (int)dataset.size(); ++i) {
		delete dataset[i];
	}
	delete update_cluster_population;
}

void calc_mean(dataset_t &dataset, data_t &mean) {
	char _verbosity = Debug;
	int ldim = mean.size();
	assert (ldim > 0);
	for (int j = 0; j < ldim; ++j) {
		mean[j] = 0;
	}
	for (int i = 0; i < (int)dataset.size(); ++i) {
		for (int j = 0; j < ldim; ++j) {
			data_t *d = dataset[i];
			mean[j] += (*d)[j];
		}
	}
	for (int j = 0; j < ldim; ++j) {
		mean[j] /= dataset.size();
	}
	fout << "Mean ";
	string sep = "";
	for (int j = 0; j < ldim; ++j, sep = ", ") {
		cout << sep << mean[j];
	}
	cout << endl;
}

void test_likelihood_function(std::default_random_engine &generator, dataset_t &dataset, ground_truth_t & ground_truth, dataset_t &ref_dataset, distribution_t & likelihood, dirichlet_process &hyper) {
	char _verbosity = Debug;

	fout << "Test likelihood" << endl;
	int kTest = 10;
	int single_out_class = 0;
	dataset_t dataset_class0and1;
	dataset_t dataset_class0;
	dataset_t dataset_class1;
	dataset_class0and1.clear();
	dataset_class0.clear();
	dataset_class1.clear();
	data_t mean_class0and1(3);
	assert (ground_truth.size() == dataset.size());
	int imax = (int)dataset.size();
	for (int i = 0; i < imax; ++i) {
		if (ground_truth[i] == single_out_class) {
			dataset_class0.push_back(dataset[i]); 
		} else {
			dataset_class1.push_back(dataset[i]); 
		}
		dataset_class0and1.push_back(dataset[i]); 
	}
	if (dataset_class0and1.size() == 0) {
		fout << "No items in dataset found" << endl;
		exit(1);
	} else {
		fout << "Found " << dataset_class0and1.size() << " items" << endl;
	}
	assert(dataset_class0and1.size() > 0);
	assert(dataset_class0.size() > 0);
	assert(dataset_class1.size() > 0);
	calc_mean(dataset_class0and1, mean_class0and1);
	data_t mean(3);
	calc_mean(ref_dataset, mean);
	/*
	   fout << "Dataset: " << endl;
	   for (int i = 0; i < (int)ref_dataset.size(); ++i) {
	   fout << *ref_dataset[i] << endl;
	   }
	   fout << "Reference dataset: " << endl;
	   for (int i = 0; i < (int)dataset_class0and1.size(); ++i) {
	   fout << *dataset_class0and1[i] << endl;
	   }*/
	// see if the ones with mu closer to that of the sample mean of the ref_dataset gets actually higher likelihood values
	string sep = "";
	double probs[kTest];
	double probs0[kTest];
	double probs1[kTest];
//	double dists[kTest];
	double suffs[kTest];
	for (int k = 0; k < kTest; ++k) {
		// sample sufficient statistics from nonparametrics
		Suffies_MultivariateNormal* suffies = hyper.sample_base(generator);
		double * s_array = suffies->mu.data();
		cout << endl;
		fout << "Suffies: ";
		sep = "";
		for (int i = 0; i < 3; ++i, sep = ", ") {
			cout << sep << s_array[i];
		}
		cout << endl;
		likelihood.init(*suffies);

		double prob, prob0, prob1;
		prob0 = likelihood.probability(dataset_class0);
		fout << "probability class 0: " << prob0 << endl;
		
		prob1 = likelihood.probability(dataset_class1);
		fout << "probability class 1: " << prob1 << endl;
		prob = likelihood.probability(dataset_class0and1);
		fout << "probability class 0 and 1: " << prob << endl;
		//delete suffies;
		/*
		double dist = 0;
		for (int j = 0; j < 3; ++j) {
			dist += (s_array[j]-mean_class0and1[j]) * (s_array[j]-mean_class0and1[j]);
		}
		fout << "Distance: " << dist << endl;
		*/
		probs[k] = prob;
		probs0[k] = prob0;
		probs1[k] = prob1;
		//dists[k] = dist;
		suffs[k] = s_array[1];
	}
	fout << "Suffs: ";
	sep = "";
	for (int k = 0; k < kTest; ++k, sep = ", ") {
		cout << sep << std::fixed << setw(9) << setfill(' ')  << setprecision(5) << suffs[k];
	}
	cout << endl;
	fout << "Prob0: ";
	sep = "";
	for (int k = 0; k < kTest; ++k, sep = ", ") {
		cout << sep << std::fixed << setw(9) << setfill(' ')  << setprecision(5) << probs0[k];
	}
	cout << endl;
	fout << "Prob1: ";
	sep = "";
	for (int k = 0; k < kTest; ++k, sep = ", ") {
		cout << sep << std::fixed << setw(9) << setfill(' ')  << setprecision(5) << probs1[k];
	}
	cout << endl;
	fout << "Mults: ";
	sep = "";
	for (int k = 0; k < kTest; ++k, sep = ", ") {
		cout << sep << std::fixed << setw(9) << setfill(' ')  << setprecision(5) << probs0[k] * probs1[k];
	}
	cout << endl;
	fout << "Probs: ";
	sep = "";
	for (int k = 0; k < kTest; ++k, sep = ", ") {
		cout << sep << std::fixed << setw(9) << setfill(' ')  << setprecision(5) << probs[k];
	}
	cout << endl;
	fout << "Accep: ";
	sep = "";
	for (int k = 0; k < kTest; ++k, sep = ", ") {
		cout << sep << std::fixed << setw(9) << setfill(' ')  << setprecision(5) << ((probs0[k] * probs1[k] > probs[k]) ? "yes" : "no");
	}
	cout << endl;
	/*
	fout << "Dists: ";
	sep = "";
	for (int k = 0; k < kTest; ++k, sep = ", ") {
		cout << sep << std::fixed << setw(9) << setfill(' ') << setprecision(5) << dists[k];
	}
	cout << endl;
	*/
}
