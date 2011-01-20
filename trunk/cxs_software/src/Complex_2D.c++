#include <iostream>  
#include "Complex_2D.h"
#include "Double_2D.h"
#include <stdlib.h>
#include <string.h>

using namespace std;

Complex_2D::Complex_2D(int x_size, int y_size){

  nx = x_size;
  ny = y_size;

  array = new fftw_complex[nx*ny];

}

Complex_2D::~Complex_2D(){

  delete[] array;

}

void Complex_2D::set_value(int x, int y, int component, double value){
  
  if(check_bounds(x,y)==FAILURE){
    cout << "can not set value out of array bounds" << endl;
    exit(1);
  }
  
  switch(component){

  case REAL :
    set_real(x,y,value);
    break;
  case IMAG :
    set_imag(x,y,value);
    break;
  default:
    cout << "Value type in Complex_2D::set_value is unknown: " 
	 << component << ". Must be REAL or IMAG" << endl;
    exit(1);
  }
}
 
double Complex_2D::get_value(int x, int y, int type) const {
  //by default we check that the value is within the bounds of the
  //array, but this can be turned off for optimisation.
  if(check_bounds(x,y)==FAILURE){
    cout << "can not get value out of array bounds" << endl;
    exit(1);
  }
  switch(type){
  case MAG:
    return get_mag(x,y);
  case REAL:
    return get_real(x,y);
  case IMAG:
    return get_imag(x,y);
  case PHASE: //goes between 0 and 2pi i.e. always positive.
    if( atan2(get_imag(x,y),get_real(x,y)) <0 )
      return atan2(get_imag(x,y),get_real(x,y)) + 2*M_PI;
    return atan2(get_imag(x,y),get_real(x,y));
  case MAG_SQ:
    return pow(get_mag(x,y),2);
  default:
    cout << "value type in Complex_2D::get_value is unknown" << endl;
    exit(1);
  }
}

void Complex_2D::get_2d(int type, Double_2D & result) const {
  
  for(int i=0; i < nx; i++)
    for(int j=0; j < ny; j++){
      result.set(i,j,get_value(i,j,type));
    }
}



void Complex_2D::scale(double scale_factor){
  
  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){ 
      array[i*ny + j][REAL]*=scale_factor;
      array[i*ny + j][IMAG]*=scale_factor;

    }
  }

}

void Complex_2D::add(Complex_2D & c2, double scale){

  if(nx!=c2.get_size_x() || ny!=c2.get_size_y()){
    cout << "in Complex_2D::add, the dimensions of the "
      "input Complex_2D do not match the dimensions of "
      "this Complex_2D object" << endl;
    exit(1);
  }

  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
      array[i*ny + j][REAL]+=scale*c2.get_real(i,j);
      array[i*ny + j][IMAG]+=scale*c2.get_imag(i,j);
    }
  }
}

void Complex_2D::multiply(Complex_2D & c2, double scale){

  if(nx!=c2.get_size_x() || ny!=c2.get_size_y()){
    cout << "in Complex_2D::add, the dimensions of the "
      "input Complex_2D do not match the dimensions of "
      "this Complex_2D object" << endl;
    exit(1);
  }
  
  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
      double new_real = c2.get_real(i,j)*this->get_real(i,j)
	- c2.get_imag(i,j)*this->get_imag(i,j);
      double new_imag = c2.get_real(i,j)*this->get_imag(i,j)
	+ c2.get_imag(i,j)*this->get_real(i,j);
      
      array[i*ny + j][REAL]=scale*new_real;
      array[i*ny + j][IMAG]=scale*new_imag;
    }
  }
}

double Complex_2D::get_norm() const {

  double norm_squared=0;

  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
      norm_squared += pow(get_imag(i,j),2)+pow(get_real(i,j),2);
    }
  }
  return sqrt(norm_squared);
}

Complex_2D * Complex_2D::clone() const {

  Complex_2D * new_complex = new Complex_2D(nx,ny);
  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
      new_complex->set_real(i,j, get_real(i,j));
      new_complex->set_imag(i,j, get_imag(i,j));
    }
  } 
 
  return new_complex;
}

void Complex_2D::copy(Complex_2D & c){
  //todo: check the bounds......
  
  memcpy(array,c.array,sizeof(fftw_complex)*nx*ny);
}


void Complex_2D::invert(){

  int middle_x = nx/2;
  int middle_y = ny/2;

  if(nx%2==1 || ny%2==1)
    cout << "WARNING: The array dimensions are odd "
	 << "but we have assumed they are even when inverting an "
	 << "array after FFT. This will probably cause you issues..."<<endl;

  for(int i=0; i < nx; ++i){
    for(int j=0; j < middle_y; ++j){
      
	int j_new = j+middle_y; 
	int i_new = i+middle_x; 

	if(i >=  middle_x)
	  i_new = i_new - 2*middle_x;

	double temp_rl = get_real(i_new,j_new);
	double temp_im = get_imag(i_new,j_new);

	set_real(i_new,j_new,get_real(i,j));
	set_imag(i_new,j_new,get_imag(i,j));

	set_real(i,j,temp_rl);
	set_imag(i,j,temp_im);
    }
  }

}


int Complex_2D::check_bounds(int x, int y) const{
 
  if(x < 0 || x >= nx || y < 0 || y >=ny )
      return FAILURE;
 
   return SUCCESS;
}
     
