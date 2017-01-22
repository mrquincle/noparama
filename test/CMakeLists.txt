include_directories ( "../include" )

set(SOURCE_FILES test_multivariate_normal_distribution.cpp)

add_executable(test_multivariate_normal_distribution ${SOURCE_FILES})

add_test(NAME test_multivariate_normal_distribution COMMAND test_multivariate_normal_distribution)

set(SOURCE_FILES test_weighted_vector.cpp)

add_executable(test_weighted_vector ${SOURCE_FILES})

add_test(NAME test_weighted_vector COMMAND test_weighted_vector)
