#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <fftw3.h>
#include <cstdlib> 
#include "Complex_2D.h"
#include "FourierT.h"
#include "FCDI_IllumRecon.h"
#include "io.h" //
#include <sstream>

using namespace std;

FCDI_IllumRecon::FCDI_IllumRecon(Complex_2D * initial_guess,
				 double beam_wavelength,
				 double zone_focal_length,
				 double focal_detector_length,
				 double pixel_size){
  complex = initial_guess;
  wavelength = beam_wavelength;
  zone_to_focal_length = zone_focal_length;
  focal_to_detector_length = focal_detector_length;
  pixel_length = pixel_size;
  
  int nx = complex->get_size_x();
  int ny = complex->get_size_y();

  support = new double*[nx];
  intensity_sqrt = new double*[nx];
  for(int i=0; i< nx; i++){
    support[i]=new double[ny];
    intensity_sqrt[i]=new double[ny];
  }

  //make the fourier transform object
  fft = new FourierT(complex->get_size_x(), complex->get_size_y());


  forward_coefficients = new Complex_2D(nx,ny);
  backward_coefficients = new Complex_2D(nx,ny);
  forward_coefficients_const = new Complex_2D(nx,ny);
  backward_coefficients_const = new Complex_2D(nx,ny);


  //set-up the coefficients
  //it's easier to do it once and reuse the matrix.
  int x_mid = nx/2;
  int y_mid = ny/2;

  double scaling = beam_wavelength*focal_to_detector_length/(pixel_length*nx);

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      double phi = 0; //2*zone_to_focal_length + 2*focal_to_detector_length;

      phi += scaling*scaling*((x_mid-i)*(x_mid-i)+(y_mid-j)*(y_mid-j))/zone_to_focal_length;
      phi += scaling*scaling*((x_mid-i)*(x_mid-i)+(y_mid-j)*(y_mid-j))/focal_to_detector_length;
      //add in the rho part
      phi *= (M_PI/beam_wavelength);

      forward_coefficients->set_real(i,j,(1/beam_wavelength)*sin(phi+2*M_PI*zone_to_focal_length/beam_wavelength));
      forward_coefficients->set_imag(i,j,(-1/beam_wavelength)*cos(phi+2*M_PI*zone_to_focal_length/beam_wavelength));

      backward_coefficients->set_real(i,j,(-1/beam_wavelength)*sin(-1*phi-2*M_PI*focal_to_detector_length/beam_wavelength));
      backward_coefficients->set_imag(i,j,(1/beam_wavelength)*cos(-1*phi-2*M_PI*focal_to_detector_length/beam_wavelength));

      forward_coefficients_const->set_real(i,j,sin(2*M_PI*zone_to_focal_length/beam_wavelength));
      forward_coefficients_const->set_imag(i,j,-1*cos(2*M_PI*zone_to_focal_length/beam_wavelength));

      backward_coefficients_const->set_real(i,j,sin(-2*M_PI*focal_to_detector_length/beam_wavelength));
      backward_coefficients_const->set_imag(i,j,-1*cos(-2*M_PI*focal_to_detector_length/beam_wavelength));


    }
  }
  
}


FCDI_IllumRecon::~FCDI_IllumRecon(){
  delete fft;
  delete forward_coefficients;
  delete backward_coefficients;
  delete forward_coefficients_const;
  delete backward_coefficients_const;

  for(int i=0; i<complex->get_size_x() ; i++){
    delete[] support[i];
    delete[] intensity_sqrt[i];
  }
  delete[] support;
  delete[] intensity_sqrt;


  //  for(int i=0; i<NTERMS-1; i++)
  //  delete x[i];

}

double ** FCDI_IllumRecon::get_intensity_autocorrelation(){

  //make a temporary object
  int nx = complex->get_size_x();
  int ny = complex->get_size_y();

  Complex_2D temp_intensity(nx,ny);
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      //set the real and imaginary components using the magnitude.
      double component = (1.0/sqrt(2.0))*(intensity_sqrt[i][j]*intensity_sqrt[i][j]);
      //scale by "pow(-1,i + j)" to make the fourier transform centered.
      component*=pow(-1,i + j);
      temp_intensity.set_value(i,j,REAL, component);
      temp_intensity.set_value(i,j,IMAG, component);
    }
  }
  
  // fourier transform the intensity 
  fft->perform_backward_fft(&temp_intensity);  

  //allocate some memory for the output.
  double ** autoc = new double*[nx];
  for(int i=0; i < nx; i++)
    autoc[i]= new double[ny];

  //get the magnitude of the fourier transformed data.
  temp_intensity.get_2d(MAG, &autoc);

  return autoc;
}


void FCDI_IllumRecon::initialise_estimate(int seed){
  //initialise the random number generator
  srand(seed);

  int nx = complex->get_size_x();
  int ny = complex->get_size_y();

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      if(!support[i][j]){ //enforce the support condition on the inital guess
	complex->set_value(i,j,REAL,0); 
	complex->set_value(i,j,IMAG,0);
      }
      else{
	double r = (255.0*rand()/(double) RAND_MAX) * pow(-1,i + j);
	double im = (255.0*rand()/(double) RAND_MAX) * pow(-1,i + j);
	complex->set_value(i,j,REAL,r); 
	complex->set_value(i,j,IMAG,im);
      }
    }
  }
}

void FCDI_IllumRecon::set_support(double ** object_support){
  for(int i=0; i< complex->get_size_x(); i++){
    for(int j=0; j< complex->get_size_y(); j++){
      support[i][j] = object_support[i][j];
    }
  }
}


void FCDI_IllumRecon::set_intensity(double ** detector_intensity){
  for(int i=0; i< complex->get_size_x(); i++){
    for(int j=0; j< complex->get_size_y(); j++){
      intensity_sqrt[i][j] = detector_intensity[i][j]; //sqrt(detector_intensity[i][j]);
    }
  }
}

/**double FCDI_IllumRecon::calculate_error(bool error){
  do_error = error;
  }**/

double FCDI_IllumRecon::get_error(){
    
    return current_error;
}


int FCDI_IllumRecon::iterate(){
  
  //start with the guess in the zone plate plane.

  //assume the wavefield at the zone plate is initalised to Pexp(...)

  int nx = complex->get_size_x();
  int ny = complex->get_size_y();
  
  double ** result = new double*[nx];
  for(int i=0; i < nx; i++)
    result[i]= new double[ny];

  // complex->get_2d(MAG,&result);
  // write_ppm("before.ppm", nx, ny, result);

  //propogate to the focal plane.
  fft->perform_forward_fft(complex);    
  
  complex->get_2d(MAG,&result);
  write_ppm("1.ppm", nx, ny, result);
  cout << "1: "<< complex->get_norm() << endl;

  //multiply by A and B and constants (from nature paper 2006)
  complex->multiply(forward_coefficients,(1.0/zone_to_focal_length));

  //forward_coefficients->get_2d(MAG,&result);
  //write_ppm("1.forward.ppm", nx, ny, result);
  //cout << "forward_norm: "<< forward_coefficients->get_norm() << endl;

  //complex->get_2d(MAG,&result);
  //write_ppm("1.5.ppm", nx, ny, result);

  //propogate to the detector plane.
  fft->perform_forward_fft(complex);

  complex->multiply(forward_coefficients_const);

  complex->get_2d(MAG,&result);
  write_ppm("2.ppm", nx, ny, result);
  cout << "2: "<< complex->get_norm() << endl;

  //impose the amplitude constraint
  scale_intensity(complex);
  cout << "Current Error is: " << current_error<<endl;

  complex->get_2d(MAG,&result);
  write_ppm("3.ppm", nx, ny, result);
  cout << "3: "<< complex->get_norm() << endl;

  //go back to the focal plane
  fft->perform_backward_fft(complex); 

  complex->multiply(forward_coefficients_const);

  complex->get_2d(MAG,&result);
  write_ppm("4.ppm", nx, ny, result,log);
  cout << "4: "<< complex->get_norm() << endl;

  //multiply by constants
  complex->multiply(backward_coefficients,1.0/focal_to_detector_length);

  //complex->get_2d(MAG,&result);
  //write_ppm("4.5.ppm", nx, ny, result);

  //go back to zone plate plane. 
  fft->perform_backward_fft(complex); 
  
  complex->get_2d(MAG,&result);
  write_ppm("5.ppm", nx, ny, result);
  cout << "5: "<< complex->get_norm() << endl;

  //apply support constraint
  project_support(complex);

  complex->get_2d(MAG,&result);
  write_ppm("6.ppm", nx, ny, result);
  cout << "6: "<< complex->get_norm() << endl;

  for(int i=0; i < nx; i++)
    delete [] result[i];
  delete [] result;

  return SUCCESS;
}

void FCDI_IllumRecon::project_support(Complex_2D * c){
  int nx = c->get_size_x();
  int ny = c->get_size_y();

  for(int i=0; i< nx; ++i){
    for(int j=0; j< ny; ++j){
      if(support[i][j]==0){
	c->set_real(i,j,0);
	c->set_imag(i,j,0);
      }
    }
  }
}


void FCDI_IllumRecon::project_intensity(Complex_2D * c){
  fft->perform_forward_fft(c);    
  scale_intensity(c);
  fft->perform_backward_fft(c);  
}


void FCDI_IllumRecon::scale_intensity(Complex_2D * c){

  double norm2_mag=0;
  double norm2_diff=0;

  int nx = c->get_size_x();
  int ny = c->get_size_y();

  for(int i=0; i< nx; ++i){
    for(int j=0; j< ny; ++j){
      //scale
      double current_mag = c->get_mag(i,j);
      if(current_mag > 0.0)
	c->set_mag(i,j,intensity_sqrt[i][j]/current_mag);

      //calculate the error
      norm2_mag += intensity_sqrt[i][j]*intensity_sqrt[i][j];
      norm2_diff += pow(current_mag-intensity_sqrt[i][j],2);
    }
  }
  current_error = (norm2_diff/norm2_mag);
}

