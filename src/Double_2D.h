#ifndef DOUBLE_2D_H
#define DOUBLE_2D_H


/**
 * @file Double_2D.h
 * @class Double_2D
 * @author  Nadia Davidson 
 *
 * @breif A 2-dimensional array of complex numbers 
 *
 * This class represents a 2D complex field. Setter and getter methods
 * are provided along with some other useful functions. Double_2D
 * objects are used in the CDI reconstruction to represent the ESW
 * in a single plane.
 */
class Double_2D{
  
  double * array;
  int nx, ny;

 public:

  /**
   * Constructor that creates a 2D object with the given dimensions.
   * 
   * @param x_size The number of samplings in the horizontal direction
   * @param y_size The number of samplings in the vertical direction
   */
  
  Double_2D():nx(0),ny(0){};

  Double_2D(int x_size, int y_size){
    allocate_memory(x_size,y_size);
  }

  /**
   * Destructor
   */
  ~Double_2D(){
    if(nx > 0 ){
      /**      for(int i=0; i < nx; ++i)
	delete [] array[i];
      delete [] array;
      }**/
      delete [] array;}
  };


  void allocate_memory(int x_size, int y_size){
    nx = x_size;
    ny = y_size;
    /**    array = new double*[nx];
    for(int i=0; i < nx; ++i)
    array[i] = new double[ny];**/
    array = new double[nx*ny];
    for(int i=0; i<nx; i++)
      for(int j=0; j<ny; j++)
	array[i*ny+j]=0;
  }

  /**
   * Set the value at point x,y. Note that this is
   * the slow method. For fast, but unsafe, methods 
   * use set_real or set_imag.
   * 
   * @param x The x position
   * @param y The y position
   * @param type Which component to set. The options are either: 
   * "REAL" or "IMAG"
   * @value The value which it will be set to
   *  
   */

  /** WARNING: no bound checking is done! */
  inline void set(int x, int y, double value){
    //     array[x][y]=value;
    array[x*ny+y]=value;
  };

  inline double get(int x, int y) const {
    //    return array[x][y]; 
    return array[x*ny+y];
  };

  /**
   * Get the size in x;
   * 
   * @return The number of horizontal points.
   *  
   */
  inline double get_size_x() const {
    return nx;
  };

  /**
   * Get the size in y;
   * 
   * @return The number of vertical points.
   *  
   */
  inline double get_size_y() const {
    return ny;
  };

  double get_sum() const {
    double total = 0;
    for(int i=0; i<nx; i++)
      for(int j=0; j<ny; j++)
	total+=array[i*ny+j];
    return total;
  };
  
  
};

#endif
