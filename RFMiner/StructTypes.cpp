#include "StructTypes.h"

namespace PI {

	bool operator==(PatternInstance const& a, PatternInstance const& b) {
		return (a.sid == b.sid && a.l == b.l && a.r == b.r);
	}

	std::size_t PI::hash_value(PatternInstance const &a) {
		std::size_t seed = 0;
		boost::hash_combine(seed, a.sid);
		boost::hash_combine(seed, a.l);
		boost::hash_combine(seed, a.r);
		return seed;
	}
}