#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <fftw3.h>

#include "Complex_2D.h"
#include "FourierT.h"

using namespace std;


FourierT::FourierT(int x_size, int y_size){

  nx = x_size;
  ny = y_size;
  
  original = new fftw_complex[nx*ny];
  transformed = new fftw_complex[nx*ny];

  f_forward = fftw_plan_dft_2d(nx, ny, original, transformed, 
			       FFTW_FORWARD, FFTW_MEASURE);
  f_backward = fftw_plan_dft_2d(nx, ny, transformed, original, 
				FFTW_BACKWARD, FFTW_MEASURE);
  
};

FourierT::~FourierT(){

  //clean up
  fftw_destroy_plan(f_forward);
  fftw_destroy_plan(f_backward);

  delete[] original;
  delete[] transformed;

}



void FourierT::perform_forward_fft(Complex_2D * c_in, Complex_2D * c_out){

  copy_to_fftw_array(original, c_in);

  fftw_execute(f_forward);

  if(!c_out){
    copy_from_fftw_array(transformed, c_in);
  }
  else{
    copy_from_fftw_array(transformed, c_out); //true=invert
  }

}


void FourierT::perform_backward_fft(Complex_2D * c_in, Complex_2D *c_out){

  copy_to_fftw_array(transformed, c_in); //true=invert
  fftw_execute(f_backward);

  if(!c_out){
    copy_from_fftw_array(original, c_in);
  }
  else{
    copy_from_fftw_array(original, c_out);
  }
}

void FourierT::copy_to_fftw_array(fftw_complex * array , Complex_2D * c, bool invert){
  //check the dimensions:
  if(!c || c->get_size_x()!=nx || c->get_size_y()!=ny ){
    cout << "In FourierT::copy_to_fftw_array. Dimensions of "
	 << "the Complex_2D and fftw_complex do not "
	 << "match.. exiting" <<endl;
    exit(1);
  }
  int middle_x = nx/2;
  int middle_y = ny/2;

  if(nx%2==1 || ny%2==1)
    cout << "WARNING: The array dimensions are odd "
	 << "but we have assumed they are even when inverting an "
	 << "array after FFT. This will probably cause you issues..."<<endl;

  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
      if(invert){
	int j_new = j+middle_y; 
	int i_new = i+middle_x; 
	if(j >=  middle_y)
	  j_new = j_new - 2*middle_y; 
	if(i >=  middle_x)
	  i_new = i_new - 2*middle_x; 
	array[(i*ny) + j][REAL]=c->get_real(i_new,j_new);
	array[(i*ny) + j][IMAG]=c->get_imag(i_new,j_new);
      }
      else{
	array[(i*ny) + j][REAL]=c->get_real(i,j);
	array[(i*ny) + j][IMAG]=c->get_imag(i,j);
      }
    }
  } 
}

void FourierT::copy_from_fftw_array(fftw_complex * array, Complex_2D * c, bool invert){
  //check the dimensions:
  if(!c || c->get_size_x()!=nx || c->get_size_y()!=ny ){
    cout << "Dimensions of the Complex_2D and "
	 << "fftw_complex do not match.. exiting" <<endl;
    
    exit(1);
  }

  //always scale as FFTW doesn't normalise the result
  double scale_factor = 1.0/(sqrt(nx*ny));

  int middle_x = nx/2;
  int middle_y = ny/2;

  if(nx%2==1 || ny%2==1)
    cout << "WARNING: The array dimensions are odd "
	 << "but we have assumed they are even when inverting an "
	 << "array after FFT. This will probably cause you issues..."<<endl;

  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){

      //tricky code here which puts the corners back into the
      //center after FFT
      if(invert){
	int j_new = j+middle_y; 
	int i_new = i+middle_x; 
	if(j >=  middle_y)
	  j_new = j_new - 2*middle_y; 
	if(i >=  middle_x)
	  i_new = i_new - 2*middle_x; 

	c->set_real(i,j,(array[i_new*ny + j_new][REAL])*scale_factor);
	c->set_imag(i,j,(array[i_new*ny + j_new][IMAG])*scale_factor);
      }
      else{
	c->set_real(i,j,(array[i*ny + j][REAL])*scale_factor);
	c->set_imag(i,j,(array[i*ny + j][IMAG])*scale_factor);
      }
    }
  }
  
}

