// Copyright 2011 Nadia Davidson for The ARC Centre of Excellence in
// Coherent X-ray Science. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

#include <iostream>
#include <Complex_2D.h>
#include <stdlib.h>
#include <cstring>
#include <fftw3.h>
#include <types.h>
#include <complex>

using namespace std;

template<class T>
ComplexR_2D<T>::ComplexR_2D(int x_size, int y_size){

  //set the array size
  nx = x_size;
  ny = y_size;
  malloc_size=sizeof(FFTW_COMPLEX)*nx*ny;
  array = (FFTW_COMPLEX*) FFTW_MALLOC(malloc_size);
  //allocate memory for the array


  //initalise the fftw plans to null (not created yet. We will
  //create them when needed to avoid unnecessary time overhead).
  f_forward = 0;
  f_backward = 0;
  fftw_type = FFTW_MEASURE;
}

/*
 * Constructor which doesn't allocate any memory, use with caution.
 */
template<class T>
ComplexR_2D<T>::ComplexR_2D(){
	f_forward = 0;
	f_backward = 0;
	fftw_type = FFTW_MEASURE;
};

/*Copy Constructor*/
template<class T>
ComplexR_2D<T>::ComplexR_2D(const ComplexR_2D & object){

  //set the array size
  nx = object.get_size_x();
  ny = object.get_size_y();

  array = (FFTW_COMPLEX*) FFTW_MALLOC(sizeof(FFTW_COMPLEX)*nx*ny);
  memcpy(array, object.array, sizeof(FFTW_COMPLEX)*nx*ny);

  //initalise the fftw plans to null (not created yet. We will
  //create them when needed to avoid unnecessary time overhead).
  f_forward = 0;
  f_backward = 0;
  fftw_type = FFTW_MEASURE;

  /*
     for(int i=0; i < object.get_size_x(); i++){
     for(int j=0; j < object.get_size_y(); j++){

     set_real(i, j, object.get_real(i, j));
     set_imag(i, j, object.get_imag(i, j));
     }
     }
   */
}

/*Copy Constructor from a Double_2D*/
template<class T>
ComplexR_2D<T>::ComplexR_2D(const Double_2D& object){

  //set the array size
  nx = object.get_size_x();
  ny = object.get_size_y();

  array=(FFTW_COMPLEX*) FFTW_MALLOC(sizeof(FFTW_COMPLEX)*nx*ny);


  //initalise the fftw plans to null (not created yet. We will
  //create them when needed to avoid unnecessary time overhead).
  f_forward = 0;
  f_backward = 0;
  fftw_type = FFTW_MEASURE;

  for(int i=0; i < object.get_size_x(); i++){
    for(int j=0; j < object.get_size_y(); j++){

      set_real(i, j, object.get(i, j));
      set_imag(i, j, 0);
    }
  }
}

/*Destructor*/
template<class T>
ComplexR_2D<T>::~ComplexR_2D(){

if (nx>0){
	FFTW_FREE(array);
}
if(f_forward)
  FFTW_DESTROY_PLAN(f_forward);
if(f_backward)
  FFTW_DESTROY_PLAN(f_backward);

}


//set the value at positions x,y. See Complex_2D.h for more info.
template<class T>
void ComplexR_2D<T>::set_value(int x, int y, int component, T value){

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
  case MAG :
    set_mag(x,y,value);
    break;
  case PHASE :
    set_phase(x,y,value);
    break;
  default:
    cout << "Value type in Complex_2D::set_value is unknown: "
      << component << ". Must be REAL, IMAG, MAG or PHASE" << endl;
    exit(1);
  }
}

//get the value at positions x,y. See Complex_2D.h for more info.
template<class T>
T ComplexR_2D<T>::get_value(int x, int y, int type) const {
  //by default we check that the value is within the bounds of the
  //array, but this can be turned off for optimisation.
 // if(check_bounds(x,y)==FAILURE){
 //   cout << "can not get value out of array bounds" << endl;
 //   exit(1);
  //}
  switch(type){
  case MAG:
    return get_mag(x,y);
  case REAL:
    return get_real(x,y);
  case IMAG:
    return get_imag(x,y);
  case PHASE: //goes between -pi and pi
    return get_phase(x,y);
  case MAG_SQ:
    return pow(get_mag(x,y),2); //the square of the magnitude
  default:
    cout << "value type in Complex_2D::get_value is unknown" << endl;
    exit(1);
  }
}

//like get() but we do it for the entire array not just a single value.
template<class T>
void ComplexR_2D<T>::get_2d(int type, Double_2D & result) const {

  for(int i=0; i < nx; i++)
    for(int j=0; j < ny; j++){
      result.set(i,j,get_value(i,j,type));
    }
}


//scale all the values in the array by the given factor.
template<class T>
void ComplexR_2D<T>::scale(T scale_factor){

  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
      array[i*ny + j][REAL]*=scale_factor;
      array[i*ny + j][IMAG]*=scale_factor;

    }
  }

}

//add another Complex_2D to this one.
template<class T>
void ComplexR_2D<T>::add(ComplexR_2D<T> & c2, T scale){

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

//multiply another Complex_2D with this one.
template<class T>
void ComplexR_2D<T>::multiply(ComplexR_2D<T> & c2, T scale){

  if(nx!=c2.get_size_x() || ny!=c2.get_size_y()){
    cout << "in Complex_2D::multiply, the dimensions of the "
      "input Complex_2D do not match the dimensions of "
      "this Complex_2D object" << endl;
    exit(1);
  }

  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
      // values are multiplied in the usual way
      // if c1 = a + ib and c2 = d + ie
      // then the new c1 is:
      // c1 = (a*d - b*e) + i(a*e + b*d)
      T new_real = c2.get_real(i,j)*this->get_real(i,j)
	- c2.get_imag(i,j)*this->get_imag(i,j);
      T new_imag = c2.get_real(i,j)*this->get_imag(i,j)
	+ c2.get_imag(i,j)*this->get_real(i,j);

      //and set the values
      set_real(i,j,scale*new_real);
      set_imag(i,j,scale*new_imag);
    }
  }
}

//multiply another Complex_2D with this one.
template<class T>
void ComplexR_2D<T>::multiply(Double_2D & c2, T scale){

  if(nx!=c2.get_size_x() || ny!=c2.get_size_y()){
    cout << "in Complex_2D::multiply, the dimensions of the "
      "input Double_2D do not match the dimensions of "
      "this Complex_2D object" << endl;
    exit(1);
  }

  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
      // values are multiplied in the usual way
      // if c1 = a + ib and c2 = d + ie
      // then the new c1 is:
      // c1 = (a*d - b*e) + i(a*e + b*d)
      T new_real = c2.get(i,j)*this->get_real(i,j);
      //- c2.get_imag(i,j)*this->get_imag(i,j);
      T new_imag = c2.get(i,j)*this->get_imag(i,j);
      //+ c2.get_imag(i,j)*this->get_real(i,j);

      //and set the values
      set_real(i,j,scale*new_real);
      set_imag(i,j,scale*new_imag);
    }
  }
}

template<class T>
T ComplexR_2D<T>::get_norm() const {

  T norm_squared=0;

  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
      norm_squared += pow(get_imag(i,j),2)+pow(get_real(i,j),2);
    }
  }
  return sqrt(norm_squared);
}

template<class T>
void ComplexR_2D<T>::conjugate() {

  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
      set_imag(i,j,-1*get_imag(i,j));
    }
  }

}

//make a new complex 2d that has idential values to this one
template<class T>
ComplexR_2D<T> * ComplexR_2D<T>::clone() const {

  ComplexR_2D<T>  * new_complex = new ComplexR_2D<T>(nx,ny);
  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
      new_complex->set_real(i,j, get_real(i,j));
      new_complex->set_imag(i,j, get_imag(i,j));
    }
  }

  return new_complex;
}


//copy another array
template<class T>
void ComplexR_2D<T>::copy(const ComplexR_2D<T> & c){

  //check the bounds
  if(c.get_size_x()!=get_size_x()||
      c.get_size_y()!=get_size_y()){
    cout << "Trying to copy an array with different dimensions... "
      << "exiting"<<endl;
    exit(1);
  }
  //copy
 std::memcpy(array,c.array,sizeof(FFTW_COMPLEX)*nx*ny);

}


//invert (and scale if we want to).
template<class T>
void ComplexR_2D<T>::invert(bool scale){

  int middle_x = nx/2;
  int middle_y = ny/2;

  T scale_factor = 1;
  if(scale)
    scale_factor = 1.0/(sqrt(nx*ny));

  if(nx%2==1 || ny%2==1)
    cout << "WARNING: The array dimensions are odd "
      << "but we have assumed they are even when inverting an "
      << "array after FFT. This will probably cause you issues..."
      << endl;

  for(int i=0; i < nx; ++i){
    for(int j=0; j < middle_y; ++j){

      int j_new = j+middle_y;
      int i_new = i+middle_x;

      if(i >=  middle_x)
	i_new = i_new - 2*middle_x;

      T temp_rl = get_real(i_new,j_new);
      T temp_im = get_imag(i_new,j_new);

      set_real(i_new,j_new,get_real(i,j)*scale_factor);
      set_imag(i_new,j_new,get_imag(i,j)*scale_factor);

      set_real(i,j,temp_rl*scale_factor);
      set_imag(i,j,temp_im*scale_factor);
    }
  }

}

template<class T>
void ComplexR_2D<T>::mirror(){

  int middle_x = nx/2;
  int middle_y = ny/2;

  //double scale_factor = 1;
  //if(scale)
  //scale_factor = 1.0/(sqrt(nx*ny));

  if(nx%2==1 || ny%2==1)
    cout << "WARNING: The array dimensions are odd "
      << "but we have assumed they are even when flipping an "
      << "array. This will probably cause you issues..."
      << endl;

  for(int i=0; i < nx; ++i){
    for(int j=0; j < middle_y; ++j){

      int j_new = ny-j-1;
      int i_new = nx-i-1;//i+middle_x;

      //if(i >=  middle_x)
      //i_new = i_new - 2*middle_x;

      T temp_rl = get_real(i,j_new);
      T temp_im = get_imag(i,j_new);

      set_real(i,j_new,get_real(i,j));
      set_imag(i,j_new,get_imag(i,j));

      set_real(i,j,temp_rl);
      set_imag(i,j,temp_im);
    }
  }
}


template<class T>
int ComplexR_2D<T>::check_bounds(int x, int y) const{

  if(x < 0 || x >= nx || y < 0 || y >=ny )
    return FAILURE;

  return SUCCESS;
}

template<class T>
void * ComplexR_2D<T>::get_array_ptr(){
	return &array;
}

template<class T>
void ComplexR_2D<T>::set_array_ptr(void * input_array, int x_size, int y_size){
	array = static_cast<FFTW_COMPLEX*>(input_array);
	nx=x_size;
	ny=y_size;
}

template<class T>
void ComplexR_2D<T>::unset_sizes(){
	nx=0;
	ny=0;
}

template<class T>
ComplexR_2D<T> ComplexR_2D<T>::get_padded(int x_add, int y_add){

  ComplexR_2D<T> padded(nx+2*x_add, ny+2*y_add);

  for(int i=0; i<x_add; i++){
    for(int j=0; j<ny+2*y_add; j++){
      padded.set_real(i, j, 0);
      padded.set_imag(i, j, 0);
    }
  }

  for(int i=nx+x_add; i<nx+2*x_add; i++){
    for(int j=0; j<ny+2*y_add; j++){
      padded.set_real(i, j, 0);
      padded.set_imag(i, j, 0);
    }
  }

  for(int i=x_add; i<nx; i++){
    for(int j=0; j<y_add; j++){
      padded.set_real(i, j, 0);
      padded.set_imag(i, j, 0);
    }
  }

  for(int i=x_add; i<nx; i++){
    for(int j=ny+y_add; j<ny+2*y_add; j++){
      padded.set_real(i, j, 0);
      padded.set_imag(i, j, 0);
    }
  }

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      padded.set_real(i+x_add, j+y_add, get_real(i, j));
      padded.set_imag(i+x_add, j+y_add, get_imag(i, j));
      //std::cout<<i<<", "<<j<<", "<<x_add<<"\n";
    }
  }

  return padded;
}

template<class T>
ComplexR_2D<T> ComplexR_2D<T>::get_unpadded(int x_add, int y_add){

  ComplexR_2D<T> unpadded(nx-2*x_add, ny-2*y_add);

  for(int i=0; i<nx-2*x_add; i++){
    for(int j=0; j<ny-2*y_add; j++){
      unpadded.set_real(i, j, get_real(i+x_add, j+y_add));
      unpadded.set_imag(i, j, get_imag(i+x_add, j+y_add));
      //      std::cout<<"Goodbye\n";

    }
  }
  return unpadded;
}

template<class T>
void ComplexR_2D<T>::initialise_fft(){

#if defined(MULTI_THREADED)
  int numthreads=FFTW_INIT_THREADS();
  FFTW_PLAN_WITH_NTHREADS(NUM_THREADS);
#endif
  //creating the plan will erase the content of the array
  //so we need to be a bit tricky here.
  FFTW_COMPLEX *tmp_array;
  tmp_array = (FFTW_COMPLEX*) FFTW_MALLOC(sizeof(FFTW_COMPLEX)*nx*ny);
  std::memcpy(array,tmp_array,sizeof(FFTW_COMPLEX)*nx*ny);
  //tmp_array = (FFTW_COMPLEX*) FFTW_MALLOC(sizeof(FFTW_COMPLEX)*nx*ny);
  //make the plans
  f_backward = FFTW_PLAN_DFT_2D(nx, ny, array, array,
				FFTW_BACKWARD, fftw_type);
  f_forward = FFTW_PLAN_DFT_2D(nx, ny,array,array,
			       FFTW_FORWARD, fftw_type);

  std::memcpy(array,tmp_array,sizeof(FFTW_COMPLEX)*nx*ny);
  FFTW_FREE(tmp_array);



}

template<class T>
void ComplexR_2D<T>::perform_forward_fft(){

  //make a new forward fft plan if we haven't made one already.
  if(f_forward==0 )
    initialise_fft();
  FFTW_EXECUTE(f_forward);

}

template<class T>
void ComplexR_2D<T>::perform_backward_fft(){

  //make a new backward fft plan if we haven't made one already.
  if(f_backward==0)
    initialise_fft();
  FFTW_EXECUTE(f_backward);
}


//this object is fourier transformed and the result placed in 'result'
template<class T>
void ComplexR_2D<T>::perform_backward_fft_real(Double_2D & result){

  FFTW_PLAN fftw;
  fftw=FFTW_PLAN_DFT_C2R_2D(nx,ny,array,result.array,FFTW_ESTIMATE);
  FFTW_EXECUTE(fftw);
  FFTW_DESTROY_PLAN(fftw);
}

//'result' is fourier transformed and the result placed in this object
template<class T>
void  ComplexR_2D<T>::perform_forward_fft_real(Double_2D & input){
  FFTW_PLAN fftw;
  fftw=FFTW_PLAN_DFT_R2C_2D(nx,ny,input.array,array,FFTW_ESTIMATE);
  FFTW_EXECUTE(fftw);
  FFTW_DESTROY_PLAN(fftw);
}

#ifdef DOUBLE_PRECISION
template class ComplexR_2D<double>;
#else
template class ComplexR_2D<float>;
#endif



