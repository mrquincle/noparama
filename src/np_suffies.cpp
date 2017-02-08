#include <np_suffies.h>

std::ostream& operator<<(std::ostream& os, const Suffies & s) {
	s.print(os);
	return os;
}

/*
std::ostream& operator<<(std::ostream& os, const Suffies_MultivariateNormal& s) {
	os << s.mu << " | " << s.sigma;
	return os;
}

std::ostream& operator<<(std::ostream& os, const Suffies_Dirichlet& s) {
	os << s.alpha;
	return os;
}
*/
