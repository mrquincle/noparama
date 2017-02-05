#include <membertrix.h>
#include <assert.h>

//using Eigen;

membertrix::membertrix() {
}

cluster_id_t membertrix::addCluster(cluster_t & cluster) {
	// use current number of columns as cluster index
	int cluster_index = _membership_matrix.cols();

	assert(cluster_index == (int)_cluster_objects.size());

	// add to column labels
	_cluster_objects.push_back(&cluster);

	// add to matrix
	_membership_matrix.conservativeResize(cluster_index + 1, Eigen::NoChange);

	return cluster_index;
}

data_id_t membertrix::addData(data_t & data) {
	// use current number of rows as data index
	int data_index = _membership_matrix.rows();

	// add to row labels
	_data_objects.push_back(&data);

	// add to matrix
	_membership_matrix.conservativeResize(Eigen::NoChange, data_index + 1);

	return data_index;
}

data_t & membertrix::getDatum(data_id_t data_id) {
	return *_data_objects[data_id];
}


np_error_t membertrix::assign(cluster_id_t cluster_id, data_id_t data_id) {
	if (_membership_matrix.row(data_id).any()) {
		return error_already_assigned;
	}

	// write value in matrix
	_membership_matrix(cluster_id, data_id) = true;

	// add data to cluster to prepare for quick access through getData(cluster)
	_clusters_dataset[cluster_id]->push_back(_data_objects[data_id]);

	return error_none;
}
		
np_error_t membertrix::retract(cluster_id_t cluster_id, data_id_t data_id) {
	if (!_membership_matrix.row(data_id).any()) {
		return error_assignment_absent;
	}

	// remove value from matrix
	_membership_matrix(cluster_id, data_id) = false;

	// remove data from cluster to update getData(cluster)
	dataset_t * cl = _clusters_dataset[cluster_id];
	auto elem = std::find(cl->begin(), cl->end(), _data_objects[data_id]);
	cl->erase(elem);

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
	return _cluster_objects;
}

dataset_t & membertrix::getData(cluster_id_t cluster_id) {
	return *_clusters_dataset[cluster_id];
}

size_t membertrix::count(cluster_id_t cluster_id) {
	return _clusters_dataset[cluster_id]->size();
}
