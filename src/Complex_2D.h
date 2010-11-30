#ifndef COMPLEX_2D_H
#define COMPLEX_2D_H


#include <math.h>

#define FAILURE 0
#define SUCCESS 1

#define REAL 0
#define IMAG 1
#define MAG 2
#define PHASE 3
#define MAG_SQ 4


/**
 * @file Complex_2D.h
 * @class Complex_2D
 * @author  Nadia Davidson 
 *
 * @breif A 2-dimensional array of complex numbers 
 *
 * This class represents a 2D complex field. Setter and getter methods
 * are provided along with some other useful functions. Complex_2D
 * objects are used in the CDI reconstruction to represent the ESW
 * in a single plane.
 */
class Complex_2D{
  
  double *** array;
  int nx, ny;

 public:

  /**
   * Constructor that creates a 2D object with the given dimensions.
   * 
   * @param x_size The number of samplings in the horizontal direction
   * @param y_size The number of samplings in the vertical direction
   */
  Complex_2D(int x_size, int y_size);

  /**
   * Destructor
   */
  ~Complex_2D();


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
  void set_value(int x, int y, int type, double value);


  /**
   * Set the real components at point x,y. Note that this is
   * an unsafe method as no bounds checking is performed.
   * 
   * @param x The x position
   * @param y The y position
   * @value The value which it will be set to
   *  
   */
  inline void set_real(int x, int y, double value){
    array[x][y][REAL]=value;
  };

  /**
   * Set the imaginary components at point x,y. Note that this is
   * an unsafe method as no bounds checking is performed.
   * 
   * @param x The x position
   * @param y The y position
   * @value The value which it will be set to
   *  
   */
  inline void set_imag(int x, int y, double value){
    array[x][y][IMAG]=value;
  };

  /**
   * Set the magnitude at point x,y. Note that this is
   * an unsafe method as no bounds checking is performed.
   * 
   * @param x The x position
   * @param y The y position
   * @value The value which it will be set to
   *  
   */
  inline void set_mag(int x, int y, double value){
    array[x][y][REAL]*=value;
    array[x][y][IMAG]*=value;
  };

  /**
   * Get the real components at point x,y. Note that this is
   * an unsafe method as no bounds checking is performed.
   * 
   * @param x The horizontal position
   * @param y The vertical position
   * @return The value at (x,y)  
   */
  inline double get_real(int x, int y){
    return array[x][y][REAL];
  };

  /**
   * Get the imaginary components at point x,y. Note that this is
   * an unsafe method as no bounds checking is performed.
   * 
   * @param x The horizontal position
   * @param y The vertical position
   * @return The value at (x,y)  
   */
  inline double get_imag(int x, int y){
    return array[x][y][IMAG];
  };
  
  /**
   * Get the magnitude at point x,y, @f$ \sqrt{\mathrm{real}^2 + \mathrm{imag}^2} @f$
   * Note that this is an unsafe method as no bounds checking is performed.
   * 
   * @param x The horizontal position
   * @param y The vertical position
   * @return The value at (x,y)  
   */
  inline double get_mag(int x, int y){
    return sqrt(array[x][y][REAL]*array[x][y][REAL]+
		array[x][y][IMAG]*array[x][y][IMAG]);
  };


  /**
   * Get the value at point x,y. Note that this is
   * the slow method. For fast, but unsafe, methods 
   * use get_real, get_imag or get_mag.
   * 
   * @param x The x position
   * @param y The y position
   * @param type Which component to set. The options are either: 
   * "REAL","IMAG","MAG","MAG_SQ" or "PHASE"
   * @return The value
   *  
   */
  double get_value(int x, int y, int type); 


  /**
   * Get the size in x;
   * 
   * @return The number of horizontal points.
   *  
   */
  double get_size_x(){
    return nx;
  };

  /**
   * Get the size in y;
   * 
   * @return The number of vertical points.
   *  
   */
  double get_size_y(){
    return ny;
  };

  /**
   * Get a mapping of the 2D array onto real space. 
   * 
   * @param type Which type of value is wanted. The options are either: 
   * "REAL","IMAG","MAG","MAG_SQ" or "PHASE"
   * @param result A pointer to a 2D array. The array will be filled with the
   * result. Note: This method does not allocated memory for the array,
   * so this should be done before making the call.
   */
  void get_2d(int type, double *** result=0);

  /**
   * Scale the real and imaginary components of the array by a factor. 
   * 
   * @param scale_factor Number which is multiplied by all components
   * of the field.
   */
  void scale(double scale_factor);

  /**
   * Add another Complex_2D to this Complex_2D. The values in this object
   * will be modified.
   * 
   * @param c2 The Complex_2D to add.
   * @param scale If this value is non-empty, or not 1, c2 will be
   * scaled before being added to the Complex_2D,@f $\mathrm{this} = \mathrm{this} +
   * \mathrm{scale} \times \mathrm{c2} $ @f . Using this function is more efficient than
   * calling Complex_2D::scale() followed by Complex_2D:add()
   * separately.
   */
  void add(Complex_2D *c2, double scale=1);


  /**
   * Multiple another Complex_2D to this Complex_2D. The values in this object
   * will be modified.
   * 
   * @param c2 The Complex_2D to add.
   * @param scale If this value is non-empty, or not 1, c2 will be
   * scaled before being added to the Complex_2D,@f $\mathrm{this} = \mathrm{this} \times
   * \mathrm{scale} \times \mathrm{c2} $ @f . Using this function is more efficient than
   * calling Complex_2D::scale() followed by Complex_2D:multiply()
   * separately.
   */
  void multiply(Complex_2D *c2, double scale=1);



  /**
   * Get the norm of this Complex_2D: @f $     $ @f.
   * 
   * @return @f $ \sqrt{ \sum_{x,y}{ \mathrm{|C(x,y)|^2} }  } $ @f.
   * @todo Check that this still works. This is not really useful, maybe I
   * should get ride of it.
   */
  double get_norm();

  /**
   * Create a new Complex_2D with the same values as this one.
   * 
   * @return The new Complex_2D
   */
  Complex_2D * clone();

  /**
   * Copy the values from another Complex_2D to this one.
   * 
   * @param c The Complex_2D which will be copied from.
   */
  void copy(Complex_2D * c);


  void invert();


 private:
    
  /**
   * Check that an (x,y) position is within the bounds of the array
   * 
   * @param x The horizontal position to check
   * @param y The vertical position to check
   */
  int check_bounds(int x, int y);
  
};

#endif
