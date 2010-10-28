#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <fftw3.h>
#include <cstdlib> 
#include "io.h"
#include "Complex_2D.h"
#include "Projection.h"
#include "FourierT.h"


using namespace std;


/**************************************/
int main(void){

  //get the data from file
  int n_x, n_y;
  double ** data;

  //read the data into an array
  int status = read_ppm("temp.ppm", &n_x, &n_y, &data);
  //int status = read_tiff("data/test_dat.tif", &n_x, &n_y, &data);
  if(!status){
    cout << "failed.. exiting"  << endl;
    return(1);
  }
  
  //fill the complex no. with image data
  Complex_2D in(n_x,n_y);
  for(int i=0; i<n_x; i++){
    for(int j=0; j<n_y; j++){
      in.set_value(i,j,REAL, 1.0/sqrt(2.0)*data[i][j]*pow(-1,i + j));
      in.set_value(i,j,IMAG, 1.0/sqrt(2.0)*data[i][j]*pow(-1,i + j));
    }
  }

  //fourier transform
  FourierT fft(n_x,n_y);
  fft.perform_forward_fft(&in);

  //write the fourier transform to file.
  double ** mag = new double*[n_x];
  for(int i=0; i < in.get_size_x(); i++)
    mag[i]= new double[n_y];

  in.get_2d(MAG, &mag);
  write_ppm("temp_f.ppm", n_x, n_y, mag);

  double ** support = new double*[n_x];
  for(int i=0; i<n_x; i++){
    support[i] = new double[n_y];
    for(int j=0; j<n_y; j++){
      if(i < n_x/4 || i > 3*n_x/4 || j < n_y/4 || j > 3*n_x/4)
	support[i][j]=0;
      else
	support[i][j]=1;
    }
  }

  //Make a random test image
  Complex_2D first_guess(n_x,n_y);
  for(int i=0; i<n_x; i++){
    for(int j=0; j<n_y; j++){
      double r = (255.0*rand()/(double) RAND_MAX)* pow(-1,i + j);
      double im = (255.0*rand()/(double) RAND_MAX)* pow(-1,i + j);
      first_guess.set_value(i,j,REAL,r); 
      first_guess.set_value(i,j,IMAG,im);
    }
  }

  Projection proj(&first_guess);
  proj.set_support(support);
  proj.set_intensity(mag);
  proj.set_algorithm(HIO);
  
  //apply the iterations
  for(int i=0; i<200; i++)
    proj.iterate();

  first_guess.get_2d(MAG,&mag);
  write_ppm("temp_rand.ppm", n_x, n_y, mag);

  //clean up
  for(int i=0; i< n_x; i++){
    delete[] mag[i];
  }

  delete[] mag;

  return 0;
}

