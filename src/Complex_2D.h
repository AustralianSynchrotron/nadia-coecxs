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

class Complex_2D{
  
  double *** array;
  int nx, ny;

 public:

  Complex_2D(int x_size, int y_size);
  ~Complex_2D();

  void set_value(int x, int y, int type, double value);

  inline void set_real(int x, int y, double value){
    array[x][y][REAL]=value;
  };

  inline void set_imag(int x, int y, double value){
    array[x][y][IMAG]=value;
  };

  inline void set_mag(int x, int y, double value){
    array[x][y][REAL]*=value;
    array[x][y][IMAG]*=value;
  };

  inline double get_real(int x, int y){
    return array[x][y][REAL];
  };

  inline double get_imag(int x, int y){
    return array[x][y][IMAG];
  };

  inline double get_mag(int x, int y){
    return sqrt(array[x][y][REAL]*array[x][y][REAL]+
		array[x][y][IMAG]*array[x][y][IMAG]);
  };

  double get_value(int x, int y, int type); 

  double get_size_x(){
    return nx;
  };

  double get_size_y(){
    return ny;
  };

  void get_2d(int type, double *** results=0);
  void scale(double scale_factor);
  void add(Complex_2D *c2, double scale=1);
  double get_norm();
  
  Complex_2D * clone();
  void copy(Complex_2D * c);

 private:
  
  int check_bounds(int x, int y);
  
};

#endif
