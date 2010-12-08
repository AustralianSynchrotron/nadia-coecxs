#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <fftw3.h>
#include <cstdlib> 
#include "Complex_2D.h"
#include "FourierT.h"
#include "PlanarCDI.h"
//#include <time.h>
#include "io.h" //
#include <sstream>

using namespace std;

PlanarCDI::PlanarCDI(Complex_2D * initial_guess)
  : PhaseRetrievalBase(initial_guess){};

int PlanarCDI::iterate(){
  
  //FS
  if(algorithm_structure[PFS]!=0){
    temp_complex_PFS->copy(complex);
    project_support(temp_complex_PFS);
    project_intensity(temp_complex_PFS);
  }

  //SF & F
  if(algorithm_structure[PF]!=0||algorithm_structure[PSF]!=0){
    temp_complex_PF->copy(complex);
    project_intensity(temp_complex_PF);
  }
  
  //combine the result of the seperate operators

  double value_real, value_imag;
  int nx = complex->get_size_x();
  int ny = complex->get_size_y();

  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){

      //Add the identity
      value_real = (1+algorithm_structure[PI])*complex->get_real(i,j);
      value_imag = (1+algorithm_structure[PI])*complex->get_imag(i,j);

      //Add the component from the PfPs operator
      if(algorithm_structure[PFS]!=0){
	value_real+=algorithm_structure[PFS]*temp_complex_PFS->get_real(i,j);
	value_imag+=algorithm_structure[PFS]*temp_complex_PFS->get_imag(i,j);
      }

      //Add the component from the Pf operator
      if(algorithm_structure[PF]!=0){
	value_real+=algorithm_structure[PF]*temp_complex_PF->get_real(i,j);
	value_imag+=algorithm_structure[PF]*temp_complex_PF->get_imag(i,j);
      }

      //Add the support
      if(support[i][j]!=0){
	//PS
	if(algorithm_structure[PS]!=0){
	  value_real+=algorithm_structure[PS]*complex->get_real(i,j);
	  value_imag+=algorithm_structure[PS]*complex->get_imag(i,j);
	}
	//PSF
	if(algorithm_structure[PSF]!=0){
	  value_real+=algorithm_structure[PSF]*temp_complex_PF->get_real(i,j);
	  value_imag+=algorithm_structure[PSF]*temp_complex_PF->get_imag(i,j);
	}
      }
      
      complex->set_real(i,j,value_real);
      complex->set_imag(i,j,value_imag);
   
    }
  }

  return SUCCESS;
}

