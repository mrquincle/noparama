include_directories ( "../include" )

set(TEST test_multivariate_normal_distribution)
set(TEST_SOURCE_FILES ../src/np_suffies.cpp ../src/statistics/multivariatenormal.cpp)
set(SOURCE_FILES ${TEST}.cpp ${TEST_SOURCE_FILES})
add_executable(${TEST} ${SOURCE_FILES})
add_test(NAME ${TEST} COMMAND ${TEST})

set(TEST test_weighted_vector)
set(SOURCE_FILES ${TEST}.cpp)
add_executable(${TEST} ${SOURCE_FILES})
add_test(NAME ${TEST} COMMAND ${TEST})

set(TEST test_mvn_likelihood)
set(TEST_SOURCE_FILES ../src/np_suffies.cpp ../src/statistics/multivariatenormal.cpp)
set(SOURCE_FILES ${TEST}.cpp ${TEST_SOURCE_FILES})
add_executable(${TEST} ${SOURCE_FILES})
add_test(NAME ${TEST} COMMAND ${TEST})

set(TEST test_membertrix)
set(TEST_SOURCE_FILES ../src/membertrix.cpp ../src/np_suffies.cpp)
set(SOURCE_FILES ${TEST}.cpp ${TEST_SOURCE_FILES})
add_executable(${TEST} ${SOURCE_FILES})
add_test(NAME ${TEST} COMMAND ${TEST})

set(TEST test_dirichlet)
set(TEST_SOURCE_FILES ../src/np_suffies.cpp ../src/statistics/multivariatenormal.cpp)
set(SOURCE_FILES ${TEST}.cpp ${TEST_SOURCE_FILES})
add_executable(${TEST} ${SOURCE_FILES})
add_test(NAME ${TEST} COMMAND ${TEST})
target_link_libraries( ${TEST} -lstdc++fs )
