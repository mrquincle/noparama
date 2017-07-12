#include <vector>
#include <iostream>
#include <fstream>

#include <Eigen/Dense>

#include <clustering_performance.h>

using namespace std;

clustering_performance::clustering_performance() {
}

matrix_t & clustering_performance::calculateContingencyMatrix(std::vector<int> & A, std::vector<int> & B) {
	auto maxA = std::max_element(A.begin(), A.end());
	assert ( maxA != A.end());

	auto maxB = std::max_element(B.begin(), B.end());
	assert ( maxB != B.end());

	int sizeA = *maxA + 1;
	int sizeB = *maxB + 1;

	cout << "Create matrix of size " << sizeA << "x" << sizeB << endl;
	matrix_t &frequencies = *new matrix_t(sizeA, sizeB);
	frequencies.setZero();
	
	for (int i = 0; i < (int)A.size(); ++i) { 
		int a = A[i];
		int b = B[i];
		frequencies(a,b) = frequencies(a,b) + 1;
	} 
	cout << frequencies << endl;

	return frequencies;
}

void clustering_performance::calculateSimilarity(matrix_t & frequencies) {

	// total number of pairs of clusters
	auto N = frequencies.sum();

	if (N == 0) {
		delete &frequencies;
		return;
	}

	auto A = frequencies;
	auto R = frequencies.rowwise().sum();
	auto C = frequencies.colwise().sum();
	
	_purity = frequencies.colwise().maxCoeff().sum() / (double)N;

	auto a = ((A.cwiseProduct(A) - A) / 2).sum();
	auto b = ((R.cwiseProduct(R) - R) / 2).sum();
	auto c = ((C.cwiseProduct(C) - C) / 2).sum();

	double S = ((N * N - N) / (double)2);
	
	if (S == 0) {
		delete &frequencies;
		return;
	}

	_rand_index = (2*a-b-c) / S + 1;

	// helper variables
	double b_PROD_c_DIV_S = b*c/S;
	double b_SUM_c_DIV_2 = (b+c)/(double)2;
	if (b_PROD_c_DIV_S == b_SUM_c_DIV_2) {
		delete &frequencies;
		return;
	}
	
	_adjusted_rand_index = (a - b_PROD_c_DIV_S) / (b_SUM_c_DIV_2 - b_PROD_c_DIV_S);

	cout << "Purity: " << _purity << endl;
	cout << "Rand Index: " << _rand_index << endl;
	cout << "Adjusted Rand Index: " << _adjusted_rand_index << endl;

	delete &frequencies;
}

void clustering_performance::write(std::string fname) {
	std::ofstream ofile;
	ofile.open(fname);

	ofile << "Purity: " << _purity << endl;
	ofile << "Rand Index: " << _rand_index << endl;
	ofile << "Adjusted Rand Index: " << _adjusted_rand_index << endl;

	ofile.close();
}
