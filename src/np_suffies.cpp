#include <np_suffies.h>

std::ostream& operator<<(std::ostream& os, const Suffies & s) {
	s.print(os);
	return os;
}

