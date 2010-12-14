#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <fftw3.h>
#include <cstdlib> 
#include <cmath>
#include "Complex_2D.h"
#include "Double_2D.h"
#include "FourierT.h"
#include "PlanarCDI.h"
//#include <time.h>
#include "io.h" //
#include <sstream>

using namespace std;

map<string,int> * PlanarCDI::algNameMap = PlanarCDI::set_up_algorithm_name_map();

map<string,int> * PlanarCDI::set_up_algorithm_name_map(){
  
  map<string,int> * temp_map = new map<string,int>;
  temp_map->insert(pair<string,int>("ER",ER));
  temp_map->insert(pair<string,int>("BIO",BIO));
  temp_map->insert(pair<string,int>("BOO",BOO));
  temp_map->insert(pair<string,int>("HIO",HIO));
  temp_map->insert(pair<string,int>("DM",DM));
  temp_map->insert(pair<string,int>("SF",SF));
  temp_map->insert(pair<string,int>("ASR",ASR));
  temp_map->insert(pair<string,int>("HPR",HPR));
  temp_map->insert(pair<string,int>("RAAR",RAAR));
  return temp_map;
}

PlanarCDI::PlanarCDI(Complex_2D * initial_guess){
  complex = initial_guess;
  
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

  temp_complex_PFS = new Complex_2D(nx,ny);
  temp_complex_PF = new Complex_2D(nx,ny);

  //set some defauts
  beta = 0.9;
  set_algorithm(HIO);

}


PlanarCDI::~PlanarCDI(){
  delete fft;

  for(int i=0; i<complex->get_size_x() ; i++){
    delete[] support[i];
    delete[] intensity_sqrt[i];
  }
  delete[] support;
  delete[] intensity_sqrt;

  //  for(int i=0; i<NTERMS-1; i++)
  //  delete x[i];

  delete temp_complex_PFS;
  delete temp_complex_PF;

  
  cout << "In the PlanarCDI destructor" << endl;

}

//double ** PlanarCDI::get_intensity_autocorrelation(){
void PlanarCDI::get_intensity_autocorrelation(Double_2D * autoc){

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
  //double ** autoc = new double*[nx];
  //for(int i=0; i < nx; i++)
  //  autoc[i]= new double[ny];

  //get the magnitude of the fourier transformed data.
  temp_intensity.get_2d(MAG, autoc);

}


void PlanarCDI::initialise_estimate(int seed){
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

void PlanarCDI::set_support(double ** object_support){
  for(int i=0; i< complex->get_size_x(); i++){
    for(int j=0; j< complex->get_size_y(); j++){
      support[i][j] = object_support[i][j];
    }
  }
}

void PlanarCDI::set_support(Double_2D object_support){
  for(int i=0; i< complex->get_size_x(); i++){
    for(int j=0; j< complex->get_size_y(); j++){
      support[i][j] = object_support.get(i,j);
    }
  }
}


void PlanarCDI::set_intensity(double ** detector_intensity){
  for(int i=0; i< complex->get_size_x(); i++){
    for(int j=0; j< complex->get_size_y(); j++){
      intensity_sqrt[i][j] = sqrt(detector_intensity[i][j]);
    }
  }
}

void PlanarCDI::set_intensity(Double_2D detector_intensity){
  for(int i=0; i< complex->get_size_x(); i++){
    for(int j=0; j< complex->get_size_y(); j++){
      intensity_sqrt[i][j] = sqrt(detector_intensity.get(i,j));
    }
  }
}


double PlanarCDI::get_error(){
  return current_error;
}

void PlanarCDI::apply_support(Complex_2D * c){
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


void PlanarCDI::project_intensity(Complex_2D * c){
  fft->perform_forward_fft(c);    
  scale_intensity(c);
  fft->perform_backward_fft(c);  
}


void PlanarCDI::scale_intensity(Complex_2D * c){

  double norm2_mag=0;
  double norm2_diff=0;

  int nx = c->get_size_x();
  int ny = c->get_size_y();

  for(int i=0; i< nx; ++i){
    for(int j=0; j< ny; ++j){
      //scale
      double current_mag = c->get_mag(i,j);
      if(current_mag > 0.0){
	c->set_mag(i,j,intensity_sqrt[i][j]/current_mag);
      }
      //calculate the error
      norm2_mag += intensity_sqrt[i][j]*intensity_sqrt[i][j];
      norm2_diff += pow(current_mag-intensity_sqrt[i][j],2);
    }
  }
  current_error = (norm2_diff/norm2_mag);

}


void PlanarCDI::set_algorithm(int alg){
  if(algorithm==alg){
    cout << "Warning you are trying to set the algorithm"
	 << " to the one already in use" << endl; 
    return;
  }

  algorithm = alg;
  
  switch(alg){

  case(ER): //checked
    set_custom_algorithm(0,0,0,1,0,0,0,0,0,0);
    break;
  case(BIO):
    set_custom_algorithm(0,beta,0,0,0,0,0,0,0,0);
    break;
  case(BOO):
    set_custom_algorithm(0,beta,0,0,0,0,0,0,1,0);
    break;
  case(HIO): //checked
    set_custom_algorithm(0,beta,1,0,0,0,0,0,0,0);
    break;
  case(DM): //checked
    set_custom_algorithm(beta,0,-1,0,-1,0,0,0,0,0);
    break;
  case(SF):
    set_custom_algorithm(0,1,0,1,0,0,0,0,0,0);
    break;
  case(ASR):
    set_custom_algorithm(0,1,1,0,0,0,0,0,0,0);
    break;
  case(HPR):
    set_custom_algorithm(0,beta,1,0,0,0,0,0,0,0);
    break;
  case(RAAR):
    set_custom_algorithm(0,1,beta,0,0,0,0,0,(1-beta),0);
    break;
  default:
    cout << "Algorithm unknown" <<endl;
  }

  print_algorithm();
}


void PlanarCDI::set_custom_algorithm(double m1, double m2, double m3, 
				      double m4, double m5, double m6, 
				      double m7, double m8,
				      double m9, double m10){
 
  algorithm_structure[PSF]=  m1 + m2 + m3 + m4;
  algorithm_structure[PFS]= -m1 + m5 + m6 + m7;
  algorithm_structure[PS] = -m3 - m6 - m8 + m10;
  algorithm_structure[PF] = -m2 - m5 + m8 + m9;
  algorithm_structure[PI] = -m4 - m7 - m9 - m10;
  
}

void PlanarCDI::print_algorithm(){
 
  cout << "x(k+1) = x(k) + ("
       << algorithm_structure[PSF]<<"*PsPf + "
       << algorithm_structure[PFS]<<"*PfPs + "
       << algorithm_structure[PS]<<"*Ps + "
       << algorithm_structure[PF]<<"*Pf + "
       << algorithm_structure[PI]<<"*I"
       << ")x(k)"<< endl;
  
}

int PlanarCDI::iterate(){

  //FS
  if(algorithm_structure[PFS]!=0){
    temp_complex_PFS->copy(complex);
    apply_support(temp_complex_PFS);
    project_intensity(temp_complex_PFS);
  }

  //SF & F
  if(algorithm_structure[PF]!=0||algorithm_structure[PSF]!=0){
    temp_complex_PF->copy(complex);
    project_intensity(temp_complex_PF);
  }

  
  //combine the result of the seperate operators

  double value_real, value_imag;
  int nx = complex->get_size_x();
  int ny = complex->get_size_y();

  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){

      //Add the identity
      value_real = (1+algorithm_structure[PI])*complex->get_real(i,j);
      value_imag = (1+algorithm_structure[PI])*complex->get_imag(i,j);

      //Add the component from the PfPs operator
      if(algorithm_structure[PFS]!=0){
	value_real+=algorithm_structure[PFS]*temp_complex_PFS->get_real(i,j);
	value_imag+=algorithm_structure[PFS]*temp_complex_PFS->get_imag(i,j);
      }

      //Add the component from the Pf operator
      if(algorithm_structure[PF]!=0){
	value_real+=algorithm_structure[PF]*temp_complex_PF->get_real(i,j);
	value_imag+=algorithm_structure[PF]*temp_complex_PF->get_imag(i,j);
      }

      //Add the support
      if(support[i][j]!=0){
	//PS
	if(algorithm_structure[PS]!=0){
	  value_real+=algorithm_structure[PS]*complex->get_real(i,j);
	  value_imag+=algorithm_structure[PS]*complex->get_imag(i,j);
	}
	//PSF
	if(algorithm_structure[PSF]!=0){
	  value_real+=algorithm_structure[PSF]*temp_complex_PF->get_real(i,j);
	  value_imag+=algorithm_structure[PSF]*temp_complex_PF->get_imag(i,j);
	}
      }
      
      complex->set_real(i,j,value_real);
      complex->set_imag(i,j,value_imag);
   
    }
  }

  return SUCCESS;
}


void PlanarCDI::apply_shrinkwrap(double gauss_width, double threshold){
  
  int nx = complex->get_size_x();
  int ny = complex->get_size_y();
  double ** recon = new double*[nx];
  for(int i=0; i < nx; ++i)
    recon[i] = new double[ny];
  complex->get_2d(MAG,&recon);
  
  //convolve
  convolve(nx,ny,&recon,gauss_width);
  
  //threshold
  apply_threshold(nx,ny,&recon,threshold);

  set_support(recon);

  for(int i=0; i < nx; ++i)
    delete [] recon[i];
  delete recon;
  
}


void PlanarCDI::convolve(int nx, int ny, double *** array, double gauss_width){

  //to speed up computation we only convolve 
  //up to 2 times the gaussian width
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

  //delete the temporary array
  for(int i=0; i < nx; ++i)
    delete [] temp_array[i];
  delete [] temp_array;

  return;
}

/** threshold is a % of the maximum */
void PlanarCDI::apply_threshold(int nx, int ny, double *** array, 
			     double threshold){
  
  //find the maximum
  double max = 0;
  for(int i=0; i < nx; i++){
    for(int j=0; j < nx; j++){
      if( (*array)[i][j] > max)
	max = (*array)[i][j];
    }
  }

  //apply the threshold
  for(int i=0; i < nx; i++){
    for(int j=0; j < nx; j++){
      if( (*array)[i][j] < (threshold*max) )
	(*array)[i][j] = 0.0;
    }
  }
}


double PlanarCDI::gauss_2d(double x, double y, 
			double sigma_x, double sigma_y, 
			double amp){  
  double x_part = (x*x)/(2.0*sigma_x*sigma_x);
  double y_part = (y*y)/(2.0*sigma_y*sigma_y);
  return amp*exp(-1*(x_part+y_part));
  
}
