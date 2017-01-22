#include <random>
#include <iostream>
#include <algorithm>

int main() {

	std::default_random_engine generator(std::random_device{}()); 
	std::uniform_real_distribution<double> distribution(0.0,1.0);        

	int K = 8;

	std::vector<double> weighted_likelihood(K);      
	for (int i = 0; i < K; ++i) {
		weighted_likelihood[i] = i*10;
	}
	std::cout << "Weighted likelihood: ";
	for (auto i: weighted_likelihood) std::cout << i << ' ';
	std::cout << std::endl;

	std::vector<double> cumsum_likelihood(K);      
	std::partial_sum(weighted_likelihood.begin(), weighted_likelihood.end(), cumsum_likelihood.begin()); 
	
	std::cout << "Cumulative sum of weighted likelihood: ";
	for (auto i: cumsum_likelihood) std::cout << i << ' ';
	std::cout << std::endl;

	std::vector<int> frequency(K); 
	
	int N = 280000;
	for (int i = 0; i < N; ++i) {
		double pick = distribution(generator) * cumsum_likelihood.back();

		auto lower = std::lower_bound(cumsum_likelihood.begin(), cumsum_likelihood.end(), pick);
		int index = std::distance(cumsum_likelihood.begin(), lower);
		frequency[index]++;
	}

	std::cout << "Frequencies: ";
	for (auto i: frequency) std::cout << i << ' ';
	std::cout << std::endl;
}
