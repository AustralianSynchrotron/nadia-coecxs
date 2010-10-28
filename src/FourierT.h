#ifndef FOURIERT_H
#define FOURIERT_H


//#include <fftw3.h>


#define FAILURE 0
#define SUCCESS 1

#define REAL 0
#define IMAG 1
#define MAG 2
#define PHASE 3

class Complex_2D;

class FourierT{

  //initialise some arrays which will 
  //hold the data and solutions.
 private:
  fftw_complex * original;
  fftw_complex * transformed;
  
  fftw_plan f_forward;
  fftw_plan f_backward;

  int nx;
  int ny;

 public:

  FourierT(int x_size, int y_size);
  ~FourierT();

  void perform_forward_fft(Complex_2D * c_in, Complex_2D * c_out=NULL);
  void perform_backward_fft(Complex_2D * c_in, Complex_2D * c_out=NULL);

  
 private:
  void copy_to_fftw_array(fftw_complex * array , Complex_2D * c);
  void copy_from_fftw_array(fftw_complex * array , Complex_2D * c);
    
};

#endif
