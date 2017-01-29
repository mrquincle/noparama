
/*!
 * Data is represented by a vector of doubles, data_t. 
 */
typedef data_t std::vector<double>;

/*!
 * A dataset is a set of data items. 
 */
typedef dataset_t std::set<data_t*>;

/*!
 * The data items can be uniquely identified by an integer assuming the data set is not too big.
 */
typedef data_id_t int;

/*!
 * The data item exists of multiple variables, each of which can be dependent or independent.
 */
enum variable_type { vt_independent_variable, vt_dependent_variable };

/*!
 * The vector can exists out of dependent and independent variables. Preferably, information about the structure of 
 * the data vector should however be separate from that of the data itself. This is hence stored in a 
 */
typedef data_variable_type_t std::vector<variable_type>;

