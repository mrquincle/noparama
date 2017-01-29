#include <membertrix.h>
#include <assert.h>

using Eigen;

membertrix::membertrix() {
}

id_cluster_t membertrix::addCluster(cluster_t & cluster) {
	// use current number of columns as cluster index
	int cluster_index = mat.columns();

	assert(cluster_index == _column_labels.size());

	// add to column labels
	_column_labels.push_back(cluster);

	// add to matrix
	_membership_matrix.conservativeResize(cluster_index + 1, NoChange);

	return cluster_index;
}

id_data_t membertrix::addData(data_t & data) {
	// use current number of rows as data index
	int data_index = mat.rows();

	// add to row labels
	_row_labels.push_back(data);

	// add to matrix
	_membership_matrix.conservativeResize(NoChange, data_index + 1);

	return data_index;
}

error_t membertrix::assign(id_cluster_t & id_cluster, id_data_t & id_data) {
	if (_membership_matrix.row(id_data).any()) {
		return error_already_assigned;
	}

	// write value in matrix
	_membership_matrix(id_cluster, id_data) = true;

	// add data to cluster set
	_cluster_sets(id_cluster).push_back(_row_labels[id_data]);

	return error_none;
}
		
error_t membertrix::retract(id_cluster_t & id_cluster, id_data_t & id_data) {
	if (!_membership_matrix.row(id_data).any()) {
		return error_assignment_absent;
	}

	// remove value from matrix
	_membership_matrix(id_cluster, id_data) = false;

	// remove value from cluster set
//	std::set<data_t>::iterator it;
//	it = _cluster_sets(id_cluster).find(_row_labels[id_data]);
//	_cluster_sets(id_cluster).erase(it);
	_cluster_sets(id_cluster).erase(_row_labels[id_data]);

	if (_membership_matrix.row(id_data).any()) {
		return error_assignment_remaining;
	}

	return error_none;
}

error_t membertrix::retract(id_data_t & id_data) {
	id_cluster = getCluster(id_data);
	retract(id_cluster, id_data);
}

id_cluster_t membertrix::getCluster(id_data_t & id_data) {
	// maxCoeff could be used, but would iterate till maximum is found, we would need find(...,'first');
	// we do not first create a col(), assuming directly accessing matrix(i,j) is faster
	for (size_t j = 0; j < _membership_matrix.cols(); ++j) {
		if (_membership_matrix(id_data, j)) {
			return j;
		}
	}
	return -1;
}

clusters_t & membertrix::getClusters() {
	return _column_labels;
}

dataset_t & membertrix::getData(id_cluster_t & id_cluster) {
	return _cluster_sets(id_cluster);
}

