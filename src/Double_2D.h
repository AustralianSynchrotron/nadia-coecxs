#ifndef REAL_2D_H
#define REAL_2D_H


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
  
  double ** array;
  int nx, ny;

 public:

  /**
   * Constructor that creates a 2D object with the given dimensions.
   * 
   * @param x_size The number of samplings in the horizontal direction
   * @param y_size The number of samplings in the vertical direction
   */
  Double_2D(int x_size, int y_size){
    nx = x_size;
    ny = y_size;
    array = new double*[nx];
    for(int i=0; i < nx; ++i)
      array[i] = new double[ny];
  }

  /**
   * Destructor
   */
  ~Double_2D(){
    for(int i=0; i < nx; ++i)
      delete [] array[i];
    delete array;
  };


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
     array[x][y]=value;
  };

  inline double get(int x, int y){
    return array[x][y];
  };

  /**
   * Get the size in x;
   * 
   * @return The number of horizontal points.
   *  
   */
  inline double get_size_x(){
    return nx;
  };

  /**
   * Get the size in y;
   * 
   * @return The number of vertical points.
   *  
   */
  inline double get_size_y(){
    return ny;
  };

  
};

#endif
