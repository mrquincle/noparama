#include <membertrix.h>
#include <assert.h>

#include <pretty_print.hpp>
#include <iostream>
#include <iomanip>

using namespace std;

#define SEPARATE_STRUCTURE

#define ADDITIONAL_CHECKS

//#define NEVER_DELETING_CLUSTERS

membertrix::membertrix() {
	_verbosity = 4;
}

cluster_id_t membertrix::addCluster(cluster_t *cluster) {
	// use current number of columns as cluster index
	int cluster_index = _membership_matrix.cols();

#ifdef NEVER_DELETING_CLUSTERS
	assert(cluster_index == (int)_cluster_objects.size());
#endif

#ifdef ADDITIONAL_CHECKS
	for (auto cluster_pair: _cluster_objects) {
		auto const &prev = cluster_pair.second;
		fout << "Added cluster " << cluster_index << ": " << cluster->getSuffies() << endl;
		// can happen if accidentally the same cluster object is added twice
		assert(cluster != prev);

		fout << "Pointer to suffies object in cluster: " << &cluster->getSuffies() << endl;
		fout << "Pointer to prev.   object in cluster: " << &prev->getSuffies() << endl;
		// can happen if clusters do not have their own suffies object
		assert(&cluster->getSuffies() != &prev->getSuffies());
	}
#endif

	// add to column labels
	fout << "Add cluster with id " << cluster_index << endl;
	_cluster_objects.insert({cluster_index, cluster});

	// add to cluster datasets
	dataset_t *dataset = new dataset_t();
	_clusters_dataset.insert({cluster_index, dataset});

	// add to matrix
	_membership_matrix.conservativeResize(Eigen::NoChange, cluster_index + 1);
	_membership_matrix.col(cluster_index).setZero();

	return cluster_index;
}

data_id_t membertrix::addData(data_t & data) {
	// use current number of rows as data index
	int data_index = _membership_matrix.rows();
	
	assert(data_index == (int)_data_objects.size());

	// add to row labels
	_data_objects.push_back(&data);

	// add to matrix
	_membership_matrix.conservativeResize(data_index + 1, Eigen::NoChange);
	_membership_matrix.row(data_index).setZero();

	return data_index;
}

data_t * membertrix::getDatum(data_id_t data_id) {
//	fout << "get data " << data_id << endl;
//	fout << data_id << ": " << *_data_objects[data_id] << endl;
	return _data_objects[data_id];
}


np_error_t membertrix::assign(cluster_id_t cluster_id, data_id_t data_id) {
	if (_membership_matrix.row(data_id).any()) {
		fout << _membership_matrix.row(data_id) << endl;
		return error_already_assigned;
	}

	// write value in matrix
	_membership_matrix(data_id, cluster_id) = true;
	
#ifdef SEPARATE_STRUCTURE
	assert (exists(cluster_id));
	// add data to cluster to prepare for quick access through getData(cluster)
	dataset_t * cl = _clusters_dataset[cluster_id];
	cl->push_back(_data_objects[data_id]);
#endif

	return error_none;
}

bool membertrix::exists(cluster_id_t cluster_id) {
	auto it = _clusters_dataset.find(cluster_id);
	return (it != _clusters_dataset.end());
}

np_error_t membertrix::retract(cluster_id_t cluster_id, data_id_t data_id) {
//	fout << "Membership matrix: " << endl << _membership_matrix << endl;
	if (!_membership_matrix.row(data_id).any()) {
		return error_assignment_absent;
	}

	// remove value from matrix
	_membership_matrix(data_id, cluster_id) = false;

//	fout << *_clusters_dataset[cluster_id] << endl;

#ifdef SEPARATE_STRUCTURE
	// remove data from the dataset that belongs to cluster_id for the purpose to have an up-to-date getData(cluster)
	dataset_t * cl = _clusters_dataset[cluster_id];
	for (auto data: *cl) {
		fout << "Cluster " << cluster_id << " data: " << *data << endl;
	}
	auto elem = std::find(cl->begin(), cl->end(), _data_objects[data_id]);
	cl->erase(elem);
	fout << "After erasing: " << data_id << " from " << cluster_id << endl;
	for (auto data: *cl) {
		fout << *data << endl;
	}

	// if retract is using the iterator over getClusters, this will cause havoc
	if (empty(cluster_id)) {
		fout << "Delete cluster: " << cluster_id << endl;
		_cluster_objects.erase(cluster_id);
	}
#endif

	if (_membership_matrix.row(data_id).any()) {
		return error_assignment_remaining;
	}

	return error_none;
}

np_error_t membertrix::retract(data_id_t data_id) {
	cluster_id_t cluster_id = getCluster(data_id);
	return retract(cluster_id, data_id);
}

cluster_id_t membertrix::getCluster(data_id_t data_id) {
	// maxCoeff could be used, but would iterate till maximum is found, we would need find(...,'first');
	// we do not first create a col(), assuming directly accessing matrix(i,j) is faster
	for (int j = 0; j < _membership_matrix.cols(); ++j) {
		if (_membership_matrix(data_id, j)) {
			return j;
		}
	}
	return -1;
}

clusters_t & membertrix::getClusters() {
	fout << "Get clusters: " << _cluster_objects.size() << std::endl;
	for (auto cluster_pair: _cluster_objects) {
		auto const &key = cluster_pair.first;
		fout << "Cluster: " << key << std::endl;
		/*
		fout << "Cluster: ";
		print(key, std::cout);
		std::cout << std::endl;
		*/
	}
	return _cluster_objects;
}

void membertrix::print(cluster_id_t cluster_id, ostream &os) {
	cluster_t *cluster = _cluster_objects[cluster_id];
	os << setw(2) << cluster_id << \
		"[#" << count(cluster_id) << "]: " << \
		cluster->getSuffies();
}

void membertrix::print(std::ostream& os) const {
	os << _membership_matrix;
}

std::ostream& operator<<(std::ostream& os, const membertrix & m) {
	m.print(os);
	return os;
}

dataset_t* membertrix::getData(cluster_id_t cluster_id) {
#ifdef SEPARATE_STRUCTURE
	return _clusters_dataset[cluster_id];
#else
	// no separate structure, so we have to populate it from scratch
	auto cluster_data = _membership_matrix.col(cluster_id);

	dataset_t *dataset = new dataset_t();
	for (int i = 0; i < cluster_data.size(); ++i) {
		if (cluster_data(i)) {
			dataset->push_back(_data_objects[i]);
		}
	}
	return dataset; 
#endif
}

bool membertrix::empty(cluster_id_t cluster_id) {
	return _clusters_dataset[cluster_id]->size() == 0;
}

size_t membertrix::count(cluster_id_t cluster_id) {
	return _clusters_dataset[cluster_id]->size();
}

int membertrix::cleanup() {
	int clusters_removed = 0;
	clusters_t &clusters = getClusters();
	for (auto cluster_pair: clusters) {
		auto const &key = cluster_pair.first;
		fout << "list cluster: " << key << endl;
	}
	for (auto iter = clusters.begin(); iter != clusters.end(); ) {
		auto const &key = iter->first;
		fout << "check cluster: " << key << endl;

		if (empty(key)) {
			fout << "delete cluster " << key << endl;
			iter = clusters.erase(iter);
			fout << "cluster deleted: " << key << endl;
			clusters_removed++;
		} else {
			iter++;
		}
	}
	return clusters_removed;
}
