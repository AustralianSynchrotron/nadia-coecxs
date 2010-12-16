#include <iostream>
#include <string>
#include <fftw3.h>
#include <cstdlib> 
#include <cmath>
#include "Complex_2D.h"
#include "PlanarCDI.h"
#include "io.h" 

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



/***************************************************************/
PlanarCDI::PlanarCDI(Complex_2D & initial_guess)
  : complex(initial_guess),
    nx(initial_guess.get_size_x()),
    ny(initial_guess.get_size_y()),
    fft(nx,ny),
    temp_complex_PFS(nx,ny),
    temp_complex_PF(nx,ny),
    beta(0.9),
    support(nx,ny),
    intensity_sqrt(nx,ny){
  
  /**  support = new double*[nx];
  intensity_sqrt = new double*[nx];
  for(int i=0; i< nx; i++){
    support[i]=new double[ny];
    intensity_sqrt[i]=new double[ny];
    }**/

  set_algorithm(HIO);

}


PlanarCDI::~PlanarCDI(){
  //delete fft;

  /**  for(int i=0; i<nx ; i++){
    delete[] support[i];
    delete[] intensity_sqrt[i];
  }
  delete[] support;
  delete[] intensity_sqrt;**/

  //  for(int i=0; i<NTERMS-1; i++)
  //  delete x[i];

  //delete temp_complex_PFS;
  //delete temp_complex_PF;

}

//double ** PlanarCDI::get_intensity_autocorrelation(){
void PlanarCDI::get_intensity_autocorrelation(Double_2D & autoc){

  //make a temporary object
  Complex_2D temp_intensity(nx,ny);
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      //set the real and imaginary components using the magnitude.
      double component = (1.0/sqrt(2.0))*(intensity_sqrt.get(i,j)*intensity_sqrt.get(i,j));
      //scale by "pow(-1,i + j)" to make the fourier transform centered.
      component*=pow(-1,i + j);
      temp_intensity.set_value(i,j,REAL, component);
      temp_intensity.set_value(i,j,IMAG, component);
    }
  }
  
  // fourier transform the intensity 
  fft.perform_backward_fft(temp_intensity);  

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

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      if(!support.get(i,j)){ //enforce the support condition on the inital guess
	complex.set_value(i,j,REAL,0); 
	complex.set_value(i,j,IMAG,0);
      }
      else{
	double r = (255.0*rand()/(double) RAND_MAX) * pow(-1,i + j);
	double im = (255.0*rand()/(double) RAND_MAX) * pow(-1,i + j);
	complex.set_value(i,j,REAL,r); 
	complex.set_value(i,j,IMAG,im);
      }
    }
  }
}

void PlanarCDI::set_support(const Double_2D & object_support){
  for(int i=0; i< nx; i++){
    for(int j=0; j< ny; j++){
      support.set(i,j,object_support.get(i,j));
    }
  }
}


void PlanarCDI::set_intensity(const Double_2D &detector_intensity){

  for(int i=0; i< nx; i++){
    for(int j=0; j< ny; j++){
      intensity_sqrt.set(i,j,sqrt(detector_intensity.get(i,j)));
    }
  }

}


double PlanarCDI::get_error(){
  return current_error;
}

void PlanarCDI::apply_support(Complex_2D & c){
  for(int i=0; i< nx; ++i){
    for(int j=0; j< ny; ++j){
      if(support.get(i,j)==0){
	c.set_real(i,j,0);
	c.set_imag(i,j,0);
      }
    }
  }
}


void PlanarCDI::project_intensity(Complex_2D & c){
  fft.perform_forward_fft(c);    
  scale_intensity(c);
  fft.perform_backward_fft(c);  
}


void PlanarCDI::scale_intensity(Complex_2D & c){

  double norm2_mag=0;
  double norm2_diff=0;

  for(int i=0; i< nx; ++i){
    for(int j=0; j< ny; ++j){
      //scale
      double current_mag = c.get_mag(i,j);
      if(current_mag > 0.0){
	c.set_mag(i,j,intensity_sqrt.get(i,j)/current_mag);
      }
      //calculate the error
      norm2_mag += intensity_sqrt.get(i,j)*intensity_sqrt.get(i,j);
      norm2_diff += pow(current_mag-intensity_sqrt.get(i,j),2);
    }
  }
  current_error = (norm2_diff/norm2_mag);

}


void PlanarCDI::set_algorithm(int alg){
  /**  if(algorithm==alg){
    cout << "Warning you are trying to set the algorithm"
	 << " to the one already in use" << endl; 
    return;
    }**/

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
 
  cout << "Currently using the algorithm: "
       << "x(k+1) = x(k) + ("
       << algorithm_structure[PSF]<<"*PsPf + "
       << algorithm_structure[PFS]<<"*PfPs + "
       << algorithm_structure[PS]<<"*Ps + "
       << algorithm_structure[PF]<<"*Pf + "
       << algorithm_structure[PI]<<"*I"
       << ")x(k)"<< endl;
  cout << "Ps - support constraint, Pf - modulus constraint" << endl;
  
}

int PlanarCDI::iterate(){

  //FS
  if(algorithm_structure[PFS]!=0){
    temp_complex_PFS.copy(complex);
    apply_support(temp_complex_PFS);
    project_intensity(temp_complex_PFS);
  }

  //SF & F
  if(algorithm_structure[PF]!=0||algorithm_structure[PSF]!=0){
    temp_complex_PF.copy(complex);
    project_intensity(temp_complex_PF);
  }

  
  //combine the result of the seperate operators

  double value_real, value_imag;

  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){

      //Add the identity
      value_real = (1+algorithm_structure[PI])*complex.get_real(i,j);
      value_imag = (1+algorithm_structure[PI])*complex.get_imag(i,j);

      //Add the component from the PfPs operator
      if(algorithm_structure[PFS]!=0){
	value_real+=algorithm_structure[PFS]*temp_complex_PFS.get_real(i,j);
	value_imag+=algorithm_structure[PFS]*temp_complex_PFS.get_imag(i,j);
      }

      //Add the component from the Pf operator
      if(algorithm_structure[PF]!=0){
	value_real+=algorithm_structure[PF]*temp_complex_PF.get_real(i,j);
	value_imag+=algorithm_structure[PF]*temp_complex_PF.get_imag(i,j);
      }

      //Add the support
      if(support.get(i,j)!=0){
	//PS
	if(algorithm_structure[PS]!=0){
	  value_real+=algorithm_structure[PS]*complex.get_real(i,j);
	  value_imag+=algorithm_structure[PS]*complex.get_imag(i,j);
	}
	//PSF
	if(algorithm_structure[PSF]!=0){
	  value_real+=algorithm_structure[PSF]*temp_complex_PF.get_real(i,j);
	  value_imag+=algorithm_structure[PSF]*temp_complex_PF.get_imag(i,j);
	}
      }
      
      complex.set_real(i,j,value_real);
      complex.set_imag(i,j,value_imag);
   
    }
  }

  return SUCCESS;
}

void PlanarCDI::get_support(Double_2D & object_support){
  for(int i=0; i < nx; i++)
    for(int j=0; j < ny; j++)
      object_support.set(i,j,support.get(i,j));
}


void PlanarCDI::apply_shrinkwrap(double gauss_width, double threshold){
  
  Double_2D recon(nx,ny);
  complex.get_2d(MAG,recon);

  //write_ppm("shrink_1.ppm",recon);

  //convolve
  convolve(recon,gauss_width);
  
  //write_ppm("shrink_2.ppm",recon);

  //threshold
  apply_threshold(recon,threshold);

  //write_ppm("shrink_3.ppm",recon);

  set_support(recon);

}


void PlanarCDI::convolve(Double_2D & array, double gauss_width){

  //to speed up computation we only convolve 
  //up to 4 pixels away from the gaussian peak
  int half_range = 4; 
  
  //make a temporary array to hold the smeared image
  Double_2D temp_array(nx,ny);
  
  //make a temporary array to hold the gaussian distribution.
  Double_2D gauss_dist(half_range+1, half_range+1);
  for(int i=0; i <= half_range; i++){
    for(int j=0; j <= half_range; j++){
      double denom = 2.0*gauss_width*gauss_width;
      gauss_dist.set(i,j,exp(-1*(i*i+j*j)/denom ) );
    }
  }      

  //now do the convolution
  //this is messy. First loop over the elements of
  //the array which was given as input
  for(int i=0; i < nx; i++){
    for(int j=0; j < ny; j++){
      
      //now loop over the colvoluted array (the one we want to make).
      //Calculate the contribution to each element in it.
      
      for(int i2=i-half_range; i2 <= i+half_range; i2++){
	for(int j2=j-half_range; j2 <= j+half_range; j2++){
	  if(i2<nx && i2>=0 && j2 >=0 && j2<ny){
	    double smeared_value = temp_array.get(i2,j2);
	    smeared_value += array.get(i,j)*gauss_dist.get(fabs(i-i2),fabs(j-j2));
	    temp_array.set(i2,j2,smeared_value); 
	  }
	}
      }
      
    }
  }

  //now copy to the original array
  for(int i=0; i < nx; i++){
    for(int j=0; j < ny; j++)
      array.set(i,j,temp_array.get(i,j));
  }
  
}

/** threshold is a % of the maximum */
void PlanarCDI::apply_threshold(Double_2D & array, 
				double threshold){
  
  //find the maximum
  double max = 0;
  for(int i=0; i < nx; i++){
    for(int j=0; j < nx; j++){
      if( array.get(i,j) > max)
	max = array.get(i,j);
    }
  }
  
  //apply the threshold
  for(int i=0; i < nx; i++){
    for(int j=0; j < nx; j++){
      if( array.get(i,j) < (threshold*max) )
	array.set(i,j,0.0);
    }
  }
}


