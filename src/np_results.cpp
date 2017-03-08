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
		fout << "Remove (previous) symlink" << endl;
		fs::remove(ws_latest);
	} else {
		fout << "No previous symlink " << ws_latest << " found" << endl;
	}

	fs::create_symlink(path, ws_latest);

	_verbosity = 1;
	// write everything out
	int k = 0;
	const clusters_t &clusters = _membertrix.getClusters();
	
	for (auto cluster_pair: clusters) {
		auto const &key = cluster_pair.first;
		fout << "Write this cluster ";
		_membertrix.print(key, cout);
		cout << endl;

		dataset_t *dataset = _membertrix.getData(key);

		stringstream sstream;
		sstream.str("");
		sstream << ws_path << '/' << basename << k << ".txt";
		string fname = sstream.str();
		fout << "Write to file " << fname << endl;
		ofstream ofile;
		ofile.open(fname);
		for (auto data: *dataset) {
			for (auto d: *data) {
				ofile << d << " ";
			}			
			ofile << endl;
		}
		ofile.close();

		k++;
	}
	
	fout << "Write to file " << basename << ".txt" << endl;
	stringstream sstream;
	sstream.str("");
	sstream << ws_path << '/' << basename << ".txt";
	string fname = sstream.str();
	writeOctave(fname);
	_verbosity = 4;
}

void Results::writeOctave(const std::string & fname) {
	ofstream ofile;
	ofile.open(fname);

	Eigen::IOFormat SpaceFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ", "", "", " ", "");
	
	const clusters_t &clusters = _membertrix.getClusters();
	int K = clusters.size();
	if (K == 0) {
		return;
	}
	auto const example = (*clusters.begin()).second;
	distribution_type_t distribution_type = example->getSuffies().distribution_type;

	ofile << "# name: mu" << endl;
	ofile << "# type: matrix" << endl;
	ofile << "# rows: " << K << endl;
	ofile << "# columns: 2" << endl;
	for (auto cluster_pair: clusters) {
		auto const &cluster = cluster_pair.second;
		switch (distribution_type) {
			case MultivariateNormal: {
				Suffies_MultivariateNormal &suffies = (Suffies_MultivariateNormal&)cluster->getSuffies();
				ofile << suffies.mu.format(SpaceFmt) << endl;
				break;
			}
			case ScalarNoise_MultivariateNormal: {
				Suffies_ScalarNoise_MultivariateNormal &suffies = (Suffies_ScalarNoise_MultivariateNormal&)cluster->getSuffies();
				ofile << suffies.mu.format(SpaceFmt) << endl;
				break;
			}
			default: {
				fout << "Needs to be implemented!" << endl;
				assert(false);
			}
		}
	}

	Eigen::IOFormat NewlineFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " \n", " \n", "", "", " ", "");
	
	ofile << endl << endl;
	ofile << "# name: sigma" << endl;

	switch (distribution_type) {
		case MultivariateNormal: 
			ofile << "# type: matrix" << endl;
			ofile << "# ndims: 3" << endl;
			ofile << " 2 2 " << K << endl;
		break;
		case ScalarNoise_MultivariateNormal:
			ofile << "# type: matrix" << endl;
			ofile << "# rows: " << K << endl;
			ofile << "# columns: 1" << endl;
		break;
		default:
			fout << "Needs to be implemented!" << endl;
	}

	for (auto cluster_pair: clusters) {
		auto const &cluster = cluster_pair.second;
		switch (distribution_type) {
			case MultivariateNormal: {
				Suffies_MultivariateNormal &suffies = (Suffies_MultivariateNormal&)cluster->getSuffies();
				ofile << suffies.sigma.format(NewlineFmt) << endl;
				break;
			}
			case ScalarNoise_MultivariateNormal: {
				Suffies_ScalarNoise_MultivariateNormal &suffies = (Suffies_ScalarNoise_MultivariateNormal&)cluster->getSuffies();
				ofile << suffies.sigma << endl;
				break;
			}
			default: {
				fout << "Needs to be implemented!" << endl;
				assert(false);
			}
		}
	}

	ofile.close();
	
	fout << "Wrote all to file" << endl;
}
