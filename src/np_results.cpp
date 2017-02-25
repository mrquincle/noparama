#include <fstream>
#include <sys/stat.h>

#include <experimental/filesystem>

#include <np_results.h>
#include <pretty_print.hpp>

namespace fs = std::experimental::filesystem;

using namespace std;

Results::Results(const membertrix & membertrix, ground_truth_t & ground_truth): 
	_membertrix(membertrix), _ground_truth(ground_truth) {
}

/*
bool compare_assignments(assignments_t::value_type &item1, assignments_t::value_type &item2) {
	return item1.second < item2.second;
}

bool compare_clusters(clusters_t::value_type &item1, clusters_t::value_type &item2) {
	return item1.first < item2.first;
}*/

matrix_t & Results::calculateContingencyMatrix() {
	// make a copy (this also reshapes the matrices)
	membertrix trix_copy;
	trix_copy = _membertrix;

	std::vector<int> gt(_ground_truth.size());
	for (auto i = 0; i < (int)_ground_truth.size(); ++i) {
		gt[i] = _ground_truth.at(i);
	}
	fout << "Ground truth: " << gt << endl;

	std::vector<int> res;
	for (auto i = 0; i < (int)trix_copy.count(); ++i) {
		auto cluster_id = trix_copy.getCluster(i);
		res.push_back(cluster_id);
	}
	
	fout << "Results: " << res << endl;

	return _clustering_performance.calculateContingencyMatrix(gt, res);


	/*
	// assume id is representative for number of clusters (should not be very large)
	assignments_t::iterator max_cluster0 = std::max_element(_ground_truth.begin(), _ground_truth.end(), compare_assignments);
	assert ( max_cluster0 != _ground_truth.end());
	int max_cluster0_id = max_cluster0->second;
	int max_clusters0 = max_cluster0_id + 1;

	auto mClusters = trix_copy.getClusters();
	clusters_t::iterator max_cluster1 = std::max_element(mClusters.begin(), mClusters.end(), compare_clusters);
	assert ( max_cluster1 != mClusters.end());
	int max_cluster1_id = max_cluster1->first;
	int max_clusters1 = max_cluster1_id + 1;

	fout << "Create matrix of size " << max_clusters0 << " by " << max_clusters1 << endl;
	matrix_t &frequencies = *new matrix_t(max_clusters0, max_clusters1);
	
	// Fill Matrix
	
	for (auto truth: _ground_truth) {
		auto data_id = truth.first;
		auto cluster_id_truth = truth.second;
		auto cluster_id_result = trix_copy.getCluster(data_id);
		frequencies(cluster_id_truth, cluster_id_result)++;
	} 
	cout << frequencies << endl;

	return frequencies;
	*/
}

void Results::write(const string &workspace, const string & path, const string & basename) {

	matrix_t & mat = calculateContingencyMatrix();

	_clustering_performance.calculateSimilarity(mat);

	bool success;

	std::string ws_path = std::string(workspace + path);

	success = fs::create_directories(ws_path);
	if (!success) {
		fout << "Error creating directory" << endl;
		// for now, just continue
	}
	std::string latest = "LATEST";
	std::string ws_latest = std::string(workspace + latest);

	if (fs::exists(ws_latest)) {
		fout << "Remove symlink" << endl;
		fs::remove(ws_latest);
	} else {
		fout << "No previous symlink " << ws_latest << " found" << endl;
	}

	fs::create_symlink(path, ws_latest);

	_verbosity = 1;
	// write everything out
	int k = 0;
	const clusters_t &clusters = _membertrix.getClusters();
	
	int K = clusters.size();

	for (auto cluster_pair: clusters) {
		auto const &key = cluster_pair.first;
		fout << "Cluster ";
		_membertrix.print(key, cout);
		cout << endl;

		dataset_t *dataset = _membertrix.getData(key);
		/*
		for (auto data: *dataset) {
			fout << *data << endl;
		}
		*/

		stringstream sstream;
		sstream.str("");
		sstream << ws_path << '/' << basename << k << ".txt";
		string fname = sstream.str();
		fout << "Print to " << fname << endl;
		ofstream ofile;
		ofile.open(fname);
		for (auto data: *dataset) {
			ofile << (*data)[0] << " " << (*data)[1] << endl;
		}

		ofile.close();

		k++;
	}
	
	for (auto cluster_pair: clusters) {
		auto const &key = cluster_pair.first;
		fout << "Cluster ";
		_membertrix.print(key, cout);
		cout << endl;
	}

	stringstream sstream;
	sstream.str("");
	sstream << ws_path << '/' << basename << ".txt";
	string fname = sstream.str();

	ofstream ofile;
	ofile.open(fname);

	Eigen::IOFormat SpaceFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ", "", "", " ", "");

	ofile << "# name: mu" << endl;
	ofile << "# type: matrix" << endl;
	ofile << "# rows: " << K << endl;
	ofile << "# columns: 2" << endl;
	for (auto cluster_pair: clusters) {
		auto const &cluster = cluster_pair.second;
		
		Suffies_MultivariateNormal &suffies = (Suffies_MultivariateNormal&)cluster->getSuffies();
		ofile << suffies.mu.format(SpaceFmt) << endl;
	}

	Eigen::IOFormat NewlineFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " \n", " \n", "", "", " ", "");
	
	ofile << endl << endl;
	ofile << "# name: sigma" << endl;
	ofile << "# type: matrix" << endl;
	ofile << "# ndims: 3" << endl;
	ofile << " 2 2 " << K << endl;

	for (auto cluster_pair: clusters) {
		auto const &cluster = cluster_pair.second;
		
		Suffies_MultivariateNormal &suffies = (Suffies_MultivariateNormal&)cluster->getSuffies();
		ofile << suffies.sigma.format(NewlineFmt) << endl;
	}

	ofile.close();
	
	_verbosity = 4;


}
