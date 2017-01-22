
/*
 * Data is a vector. The vector can exists out of dependent and independent variables.
 * A dataset is a set of data items. 
 */

typedef data_t std::vector<double>;

typedef dataset_t std::set<data_t*>;


/*
struct data {
	union {
//		double raw[];
//		double independent;
//		double dependent[];
	}
} data_t;
*/
