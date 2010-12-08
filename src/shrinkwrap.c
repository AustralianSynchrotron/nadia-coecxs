
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <cmath>
#include "shrinkwrap.h"

using namespace std;

void apply_shrinkwrap(int nx, int ny, double *** recon,
		      double gauss_width, double threshold){

  //convolve
  convolve(nx,ny,recon,gauss_width);

  //threshold
  apply_threshold(nx,ny,recon,threshold);

}


//this is going to be very slow!!
void convolve(int nx, int ny, double *** array, double gauss_width){

  //only convolve up to 2 times the gaussian width
  const int cut_off = 2; 
  
  int half_range = cut_off*gauss_width;

  //make a temporary array
  double ** temp_array = new double*[nx];
  for(int i=0; i < nx; ++i){
    temp_array[i] = new double[ny];
    for(int j=0; j < ny; j++)
      temp_array[i][j]=0;
  }

  //now do the convolution
  //this is messy. First loop over the elements of
  //the array which was given as input
  for(int i=0; i < nx; i++){
    for(int j=0; j < ny; j++){
      
      //now loop over the colvoluted array we want to make.
      //Calculate the contribution to each element in it.
      
      for(int i2=i-half_range; i2 <= i+half_range; i2++){
	for(int j2=j-half_range; j2 <= j+half_range; j2++){
	  if(i2<nx && i2>=0 && j2 >=0 && j2<ny){
	    
	    temp_array[i2][j2]+=gauss_2d(i-i2,j-j2,
					 gauss_width,gauss_width, 
					 (*array)[i][j]);
	  }
	}
      }
      
    }
  }


  //now copy to the original array
  for(int i=0; i < nx; i++){
    for(int j=0; j < ny; j++)
      (*array)[i][j] = temp_array[i][j];
  }
  
  return;
}

/** threshold is a % of the maximum */
void apply_threshold(int nx, int ny, double *** array, 
	       double threshold){

  //find the maximum
  double max = 0;
  for(int i=0; i < nx; i++){
    for(int j=0; j < nx; j++){
      if( (*array)[i][j] > max)
	max = (*array)[i][j];
    }
  }

  cout << "max is:"<<max<<endl;

  for(int i=0; i < nx; i++){
    for(int j=0; j < nx; j++){
      if( (*array)[i][j] < (threshold*max) )
	(*array)[i][j] = 0.0;
    }
  }
}


double gauss_2d(double x, double y, 
		double sigma_x, double sigma_y, 
		double amp){  
  double x_part = (x*x)/(2.0*sigma_x*sigma_x);
  double y_part = (y*y)/(2.0*sigma_y*sigma_y);
  return amp*exp(-1*(x_part+y_part));
  
}
