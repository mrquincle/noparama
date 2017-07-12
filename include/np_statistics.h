#pragma once

#include <string.h>

typedef enum { simple_random_split, sams_prior, sams_random_walk, random_mixing } split_method_t;

typedef struct step {
	std::string type;
	int attempts;
	int cluster_events_accept;
	int cluster_events_reject;
	int target_likelihood_zero;
	int source_likelihood_zero;
	int likelihood_both_zero;
	int likelihood_both_nonzero;
	int likelihood_both_nonzero_accept;
	int likelihood_both_nonzero_reject;
} step_t;

struct statistics_t {
	step_t step[4];
};

