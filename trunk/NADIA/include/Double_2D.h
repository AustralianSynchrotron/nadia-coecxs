// Copyright 2011 Nadia Davidson for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

/**
 * @file Double_2D.h
 * @class Double_2D
 * @author  Nadia Davidson 
 * @date Last modified on 4/5/2012 by T'Mir D. Julius
 * 
 *
 * @brief A 2-dimensional array of doubles
 *
 * This class represents a 2D field of doubles. Setter and getter
 * methods are provided along with some other useful functions. This
 * is a convenient form to store and pass image data in the
 * reconstruction software.  A better implementation of this class
 * would use templates, but restricting the type to double is
 * sufficient our needs.
 */

#ifndef DOUBLE_2D_H
#define DOUBLE_2D_H
//#define DOUBLE_PRECISION 

#include <math.h>
#include <cstring>
#include <iostream>

template <class T>
class Real_2D{

 protected:

  /** the underlying 2-D array */
  T * array;
  
  /** the size in x */
  int nx;

  /** the size in y */
  int ny;

 public:

  template <class TT>friend class ComplexR_2D;
  /**
   * A constructor which creates an empty array (of no size).  Note
   * that memory has not been allocated if this method is used.
   */
  Real_2D():nx(0),ny(0){};
  
  /**
   * Constructor that creates a 2D object with the given dimensions.
   * 
   * @param x_size The number of samplings in the horizontal direction
   * @param y_size The number of samplings in the vertical direction
   */
  Real_2D(int x_size, int y_size){
    allocate_memory(x_size,y_size);
  };
  
  /**
   * Destructor. Memory is deallocated here.
   */
  ~Real_2D(){
    if(nx > 0 ){
      delete [] array;
    }  
    
  };
     
  /**
   * Allocate memory for the array. This should only be used if
   * the constructor was called with no parameters!
   * 
   * @param x_size The number of samplings in the horizontal direction
   * @param y_size The number of samplings in the vertical direction
   */ 
  void allocate_memory(int x_size, int y_size){
    nx = x_size;
    ny = y_size;
    array = new T[nx*ny];
    memset(array, 0, sizeof(T)*nx*ny);

    return;
    /**    if(!array){
      std::cout << "Could not allocated memory. Exiting.."<<std::endl;
      exit();
      }**/
    //std::cout<<"Sizeof(T) is "<<sizeof(T)*nx*ny<<"\n\n";

//    memset(array, 0, sizeof(T)*nx*ny);
    /*for(int i=0; i<nx; i++)
      for(int j=0; j<ny; j++)
	array[i*ny+j]=0;
*/
  };


  /**
   * Copy the contents of another Real_2D array to this one.  Note
   * that this is a quick copy method and no bounds checking is done.
   * 
   * @param other_array The array to copy from
   */
  void copy(const Real_2D<T> & other_array){
    memcpy(array,other_array.array, sizeof(T)*nx*ny);
  };

  /**
   * Set the value at positions (x,y) WARNING: no bound checking is
   * done!
   *
   * @param x The horizontal position 
   * @param y The vertical position
   * @param value The value to set
   */
  inline void set(int x, int y, T value){
    array[x*ny+y]=value;
  };


  /**
   * Get the value at positions (x,y) WARNING: no bound checking is
   * done!
   *
   * @param x The horizontal position 
   * @param y The vertical position
   * @return The value at (x,y)
   */
  inline T get(int x, int y) const {
    if(x < 0 || x > nx){
      return(0);
    }
    if(y < 0 || y > ny){
      return(0);
    } 

    return array[x*ny+y];
  };

  /**
   * Get the size in x;
   * 
   * @return The number of horizontal points.
   *  
   */
  inline int get_size_x() const {
    return nx;
  };

  /**
   * Get the size in y;
   * 
   * @return The number of vertical points.
   *  
   */
  inline int get_size_y() const {
    return ny;
  };

  /**
    *Find the mirror image of the image
   */
  void mirror(){

    int middle_x = nx/2;
    int middle_y = ny/2;

    //double scale_factor = 1;
    //if(scale)
    //scale_factor = 1.0/(sqrt(nx*ny));

    if(nx%2==1 || ny%2==1)
      std::cout << "WARNING: The array dimensions are odd "
	<< "but we have assumed they are even when flipping an "
	<< "array. This will probably cause you issues..."
	<< std::endl;

    for(int i=0; i < nx; ++i){
      for(int j=0; j < middle_y; ++j){

	int j_new = ny-j-1;
	int i_new = nx-i-1;//i+middle_x;

	//if(i >=  middle_x)
	//i_new = i_new - 2*middle_x;

	double temp_rl = get(i,j_new);

	set(i,j_new,get(i,j));

	set(i,j,temp_rl);
      }
    }
  }



  /**
   * Get the sum of all values in the array. This is useful to
   * determine normalisation values.
   * 
   * @return The sum of all values in the array
   *  
   */
  T get_sum() const{
    T total = 0;
    for(int i=0; i<nx; i++)
      for(int j=0; j<ny; j++)
	total+=array[i*ny+j];
    return total;
  };

  /**
   * Get the sum of the abolute value of all elements in the array.
   * 
   * @return The sum of the absolute values of all elements in the array
   *  
   */
  T get_abs_sum() const{
    T total = 0;
    for(int i=0; i<nx; i++)
      for(int j=0; j<ny; j++)
	total+=fabs(array[i*ny+j]);
    return total;
  };


  /**
   * Get the average of the values of the elements in the array.
   *
   *@return The average of the values of the elements in the array.
   *
   *
   */
  T get_ave() const{
    T sum = get_sum();
    sum = sum/(get_size_x()*get_size_y()); 
    return sum;
  }

  /**
   * Get the maximum of all values in the array. 
   * 
   * @return The maximum value in the array
   *  
   */
  T get_max() const{
    if(nx==0||ny==0)
      return 0;

    T max = array[0];
    for(int i=0; i<nx; i++){
      for(int j=0; j<ny; j++){
	if(array[i*ny+j]>max)
	  max = array[i*ny+j];
      }
    }
    return max;

  };

  /**
   * Get the minimum of all values in the array. 
   * 
   * @return The minimum value in the array 
   */
  T get_min() const {
    if(nx==0||ny==0)
      return 0;

    T min = array[0];
    for(int i=0; i<nx; i++){
      for(int j=0; j<ny; j++){
	if(array[i*ny+j]<min)
	  min = array[i*ny+j];
      }
    }
    return min;

  };

  void add(const Real_2D<T> & other_array, double norm=1.0){
    for(int i=0; i< nx; i++)
      for(int j=0; j< ny; j++)
	array[i*ny+j]+=norm*other_array.get(i,j);
  };

  void scale(double scale_by){
    for(int i=0; i< nx; i++)
      for(int j=0; j< ny; j++)
	array[i*ny+j]*=scale_by;
  };


  /**
   * Square all array values in-place. 
   */
  void square(){
    for(int i=0; i< nx; i++)
      for(int j=0; j< ny; j++)
	array[i*ny+j] *= array[i*ny+j];
  };

  /**
   * Square-root all array values in-place. 
   */
  void sq_root(){
    for(int i=0; i< nx; i++)
      for(int j=0; j< ny; j++)
	array[i*ny+j] = sqrt(array[i*ny+j]);
  };

  /**
   * Set the copy constructor
   */

  Real_2D(const Real_2D& object){

    allocate_memory(object.get_size_x(),object.get_size_y());
    memcpy(array,object.array, sizeof(T)*nx*ny);

  };

  /**
   * Set the assignment operator
   */
  Real_2D& operator=(const Real_2D& rhs){

    //Clean up
    if(nx > 0 ){
      delete [] array;
    }

    //Construct again
    allocate_memory(rhs.get_size_x(),rhs.get_size_y());
    memcpy(array, rhs.array, sizeof(T)*nx*ny);
    /*
       for(int i =0; i<rhs.get_size_x(); i++){
       for(int j=0; j<rhs.get_size_y(); j++){
       set(i, j, rhs.get(i, j));
       }
       }*/
  };

};

#ifndef DOUBLE_PRECISION
typedef Real_2D<float> Double_2D;
#else
typedef Real_2D<double> Double_2D;
#endif


#endif
