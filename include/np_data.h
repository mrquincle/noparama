#include <vector>
#include <set>

/**
 * Data is represented by a vector of doubles, data_t. 
 */
typedef std::vector<double> data_t;

/*!
 * A dataset is a set of data items. 
 */
typedef std::set<data_t*> dataset_t; 

/*!
 * The data items can be uniquely identified by an integer assuming the data set is not too big.
 */
typedef int data_id_t;

/*!
 * The data item exists of multiple variables, each of which can be dependent or independent.
 */
enum variable_type { vt_independent_variable, vt_dependent_variable };

/*!
 * The vector can exists out of dependent and independent variables. Preferably, information about the structure of 
 * the data vector should however be separate from that of the data itself. This is hence stored in a 
 */
typedef std::vector<variable_type> data_variable_type_t;

