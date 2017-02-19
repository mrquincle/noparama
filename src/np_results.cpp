#include <fstream>
#include <sys/stat.h>

#include <Eigen/Dense>

#include <np_results.h>

#include <pretty_print.hpp>

using namespace std;

Results::Results(const membertrix & membertrix, ground_truth_t & ground_truth): 
	_membertrix(membertrix), _ground_truth(ground_truth) {
}

void Results::write(const string & path, const string & basename) {
	const int err = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if (err != 0) {
		fout << "Error creating directory" << endl;
		// for now, just continue
	}

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
		sstream << path << '/' << basename << k << ".txt";
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
	sstream << path << '/' << basename << ".txt";
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
