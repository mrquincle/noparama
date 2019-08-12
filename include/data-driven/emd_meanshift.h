#pragma once

// currently, experiment with objects in 2D
static const int dim = 3;

/**
 * Adaptation of EMD algorithm that calculates offsets and takes these into account.
 */
/*
void emd_meanshift(int b, int n, int m, const float * xy1, const float * xy2, float * match, 
    float * offset1, float * offset2);
*/
/**
 * Same offset for all points in the same cloud.
 */
void emd_global_offset(int b, int n, int m, const float * xy1, const float * xy2, float * match, 
    const float * offset1, const float * offset2);

void emd_mean_costs_global_offset(int b, int n, int m, float * xy1, float * xy2, float * match, 
    float * offset1, float * offset2, float * cost);
