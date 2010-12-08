#include <iostream>  
#include "Complex_2D.h"
#include <stdlib.h>


using namespace std;

Complex_2D::Complex_2D(int x_size, int y_size){

  nx = x_size;
  ny = y_size;

  array = new double**[nx];
  
  for(int i=0; i < nx; ++i){
    array[i] = new double*[ny];
    for(int j=0; j < ny; ++j)
      array[i][j] = new double[2];
  }  
}

Complex_2D::~Complex_2D(){

  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
      delete[] array[i][j];
    }  
    delete[] array[i];
  }
  delete[] array;

}

void Complex_2D::set_value(int x, int y, int component, double value){
  
  if(check_bounds(x,y)==FAILURE){
    cout << "can not set value out of array bounds" << endl;
    exit(1);
  }
  
  switch(component){

  case REAL :
    array[x][y][REAL]=value;
    break;
  case IMAG :
    array[x][y][IMAG]=value;
    break;
  default:
    cout << "Value type in Complex_2D::set_value is unknown: " 
	 << component << ". Must be REAL or IMAG" << endl;
    exit(1);
  }
}
 
double Complex_2D::get_value(int x, int y, int type){
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
    /**    if(get_real(x,y)==0){
      if(get_imag(x,y)==0)
    	return 0;
      if(get_imag(x,y)>0)
	return M_PI/2.0;
      else
	return 3*M_PI/2.0;
	}**/
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


void Complex_2D::get_2d(int type, double *** result){
  
  if(result==0){
    result = new double**;
    *result = new double*[nx];  
    for(int i=0; i < nx; ++i)
      (*result)[i] = new double[ny];
  }

  for(int i=0; i < nx; i++)
    for(int j=0; j < ny; j++){
      (*result)[i][j] = get_value(i,j,type);
    }
}



void Complex_2D::scale(double scale_factor){
  
  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
      array[i][j][REAL]*=scale_factor;
      array[i][j][IMAG]*=scale_factor;
    }
  }
  //  cout << "scaled array[0][0][0]="<<array[0][0][0]<<endl;

}

void Complex_2D::add(Complex_2D * c2, double scale){

  if(nx!=c2->get_size_x() || ny!=c2->get_size_y()){
    cout << "in Complex_2D::add, the dimensions of the "
      "input Complex_2D do not match the dimensions of "
      "this Complex_2D object" << endl;
    exit(1);
  }
  
  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
      array[i][j][REAL]+=scale*c2->get_real(i,j);
      array[i][j][IMAG]+=scale*c2->get_imag(i,j);
    }
  }

}

void Complex_2D::multiply(Complex_2D * c2, double scale){

  if(nx!=c2->get_size_x() || ny!=c2->get_size_y()){
    cout << "in Complex_2D::add, the dimensions of the "
      "input Complex_2D do not match the dimensions of "
      "this Complex_2D object" << endl;
    exit(1);
  }
  
  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
      double new_real = c2->get_real(i,j)*this->get_real(i,j)
	- c2->get_imag(i,j)*this->get_imag(i,j);
      double new_imag = c2->get_real(i,j)*this->get_imag(i,j)
	+ c2->get_imag(i,j)*this->get_real(i,j);
      array[i][j][REAL]=scale*new_real;
      array[i][j][IMAG]=scale*new_imag;
    }
  }
}

double Complex_2D::get_norm(){

  double norm_squared=0;

  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
      norm_squared += pow(array[i][j][IMAG],2) + pow(array[i][j][REAL],2);
    }
  }
  return sqrt(norm_squared);
}

Complex_2D * Complex_2D::clone(){

  Complex_2D * new_complex = new Complex_2D(nx,ny);
  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
      new_complex->set_real(i,j, array[i][j][REAL]);
      new_complex->set_imag(i,j, array[i][j][IMAG]);
    }
  } 
 
  return new_complex;
}

/** Fill this Complex_2D with the values from "c". */
void Complex_2D::copy(Complex_2D * c){

  //todo: check the bounds......

  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
      array[i][j][REAL]=c->get_real(i,j);
      array[i][j][IMAG]=c->get_imag(i,j);
    }
  }  

}

void Complex_2D::invert(){

  //Complex_2D * temp = clone();
    
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

	//	if(j >=  middle_y)
	//	  j_new = j_new - 2*middle_y; 
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


int Complex_2D::check_bounds(int x, int y){
 
  if(x < 0 || x >= nx || y < 0 || y >=ny )
      return FAILURE;
 
   return SUCCESS;
}
     
