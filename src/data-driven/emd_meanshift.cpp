#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <cassert>
#include <fstream>

#include <data-driven/emd_meanshift.h>

/*
extern "C" {
#include <sort_indices.h>
}

void emd_meanshift(int b, int n, int m, const float * cloud1, const float * cloud2, float * match, 
    float * offset1, float * offset2) {

  // here offset is just per point, not per pair
  calc_offset(n, cloud1, offset1);
  calc_offset(m, cloud2, offset2);

  emd_standard(b, n, m, cloud1, cloud2, match, offset1, offset2);
}
*/

double calc_dist(double dx, double dy, double dz) {
// return sqrt(dx*dx + dy*dy + dz*dz);
  double result = dx*dx + dy*dy + dz*dz;
  return sqrt(result);
}

void emd_global_offset(int b, int n, int m, const float * cloud1, const float * cloud2, float * match, 
    const float * offset1, const float * offset2) {

//  std::cout << "n = " << n << ", m = " << m << std::endl;

  /*
  if (m > n) {
    assert(m % n == 0);
  } else {
    assert(n % m == 0);
  }*/

  for (int i=0;i<b;i++){

    // decompose in such way that one factor is 1 and the other factor defines how often the cloud point fits
    int factorl=std::max(n,m)/n;
    int factorr=std::max(n,m)/m;
  
    // saturation says something about convergence, initialize at factor l and factor r.
    std::vector<double> saturatedl(n,double(factorl)), saturatedr(m,double(factorr));
    // weights for each pair of points
    std::vector<double> weight(n*m);
    // init match matrix to 0
    for (int j=0;j<n*m;j++)
      match[j]=0;
    // iterate over 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2
    for (int j=8;j>=-2;j--){
      // level is then -65536, -16384, ..., -4, -1, -1/4, -1/16, the latter of which is set to 0
      double level=-powf(4.0,j);
      if (j==-2)
	level=0;
      for (int k=0;k<n;k++){
	// note how we write offset1[0] rather than offset1[k*dim+0], offset is always the same
	double x1=cloud1[k*dim+0] - offset1[0];
	double y1=cloud1[k*dim+1] - offset1[1];
	double z1=cloud1[k*dim+2] - offset1[2];
	// this iterates over all points, they all have the same offset
	for (int l=0;l<m;l++){
	  double x2=cloud2[l*dim+0] - offset2[0];
	  double y2=cloud2[l*dim+1] - offset2[1];
	  double z2=cloud2[l*dim+2] - offset2[2];

	  // this is not sparse, that's why almost any point wants to contribute to other points
	  // even if there distance is not small
	  //double dist = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
	  double dist = calc_dist(x1-x2,y1-y2,z1-z2);
	  weight[k*m+l]=expf(level*dist)*saturatedr[l];
	  /*
	  if (j == 8) {
	    std::cout << "[" << k << "," << l << "]" << dist << ", ";
	    if (dist > 25) {
	      std::cout << std::endl;
//	      std::cout << "x1: " << x1 << ", ";
//	      std::cout << "y1: " << y1 << ", ";
//	      std::cout << "z1: " << z1 << std::endl;
	      std::cout << "x2: " << x2 << ", ";
	      std::cout << "y2: " << y2 << ", ";
	      std::cout << "z2: " << z2 << std::endl;
	      std::cout << "offset2 " << offset2[0] << ", ";
	      std::cout << offset2[1] << ", ";
	      std::cout << offset2[2] << std::endl;
	    }
//	    std::cout << saturatedr[l] << std::endl;
//	    std::cout << weight[k*m+l] << ',';
	  }*/
	}
      }
      // vector ss is sum for each l
      std::vector<double> ss(m,1e-9);
      for (int k=0;k<n;k++){
	double s=1e-9;
	// sum all weights
	for (int l=0;l<m;l++){
	  s+=weight[k*m+l];
	}
	// normalize with sum and multiply each point in k with saturation l
	for (int l=0;l<m;l++){
	  weight[k*m+l]=weight[k*m+l]/s*saturatedl[k];
	}
	// sum again for each point in l
	for (int l=0;l<m;l++) {
	  ss[l]+=weight[k*m+l];
	}
      }
      // normalize now over l
      for (int l=0;l<m;l++){
	double s=ss[l];
	double r=std::min(saturatedr[l]/s,1.0);
	ss[l]=r;
      }

      // bidding phase...
      // we should somehow find the biggest value v_i and then the second best entry w_i
      // vector ss2 is yet another sum
      std::vector<double> ss2(m,0);
      for (int k=0;k<n;k++){
	double s=0;
	for (int l=0;l<m;l++){
	  // we multiply the weights with ss
	  weight[k*m+l]*=ss[l];
	  // we add them to the sum s
	  s+=weight[k*m+l];
	  // we add them also to the sum ss2
	  ss2[l]+=weight[k*m+l];
	}
	// here we calculate saturated l as saturated l minus s
	saturatedl[k]=std::max(saturatedl[k]-s,0.0);
      }
      // write match matrix by adding weight, how is it only 0 or 1, it does not seem so.
      for (int kk=0;kk<n*m;kk++) {
	match[kk]+=weight[kk];
      }

      // saturation of r minus ss2
      for (int l=0;l<m;l++){
	saturatedr[l]=std::max(saturatedr[l]-ss2[l],0.0);
      }
    }
    cloud1+=n*dim;
    cloud2+=m*dim;
    match+=n*m;
  }
}

void emd_mean_costs_global_offset(int b, int n, int m, float * cloud1, float * cloud2, float * match, 
    float * offset1, float * offset2, float * cost) {

  *cost = 0.0;
  for (int i=0;i<b;i++){
    double s=0;
    for (int j=0;j<n;j++) {
      double x1=cloud1[j*dim+0] - offset1[0];
      double y1=cloud1[j*dim+1] - offset1[1];
      double z1=cloud1[j*dim+2] - offset1[2];
      for (int k=0;k<m;k++){
	double dx=(x1 - (cloud2[k*dim+0] - offset2[0]));
	double dy=(y1 - (cloud2[k*dim+1] - offset2[1]));
	double dz=(z1 - (cloud2[k*dim+2] - offset2[2]));
//	float d=sqrtf(dx*dx+dy*dy+dz*dz)*match[j*m+k];
	double d = calc_dist(dx,dy,dz) * match[j*m+k];
	s+=d;
      }
    }
    *cost+=(s/(n*m));
    cloud1+=n*dim;
    cloud2+=m*dim;
    match+=n*m;
  }
}
