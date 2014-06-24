// Copyright 2011 Nadia Davidson for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

/**
 * @file cplx2tiff.c
 * 
 * \a cplx2tiff.exe - Extract part of a complex binary file (2D fftw
 * format) and save as a tiff (grey-scale, 16 bit file).
 *  The real, imaginary, magnitude, phase and magnitude squared
 * can be extracted.
 * 
 * \par Usage: cplx2tiff.exe \<input cplx file\> \<output tiff file\> \<component type\> \<size in x\> \<size in y\> 
 * \par
 * where component type is one of:
 * - 0: REAL 
 * - 1: IMAG 
 * - 2: MAG 
 * - 3: PHASE 
 * - 4: MAG_SQ 
 *
 * \par Example:
 * \verbatim cplx2tiff.exe my_white_field.cplx my_white_field_illum.pgm 4 1024 1024 \endverbatim
 * Extract the magnitude squared from a reconstucted white-field file. 
 * 
 **/

#include <iostream>
#include <stdlib.h>
#include <io.h>
#include <Complex_2D.h>
#include <Double_2D.h>

using namespace std;

/**************************************/
int main(int argc, char * argv[]){

  //check for 5 arguements
  if(argc!=6 ){
    cout << "Wrong number of arguments." << endl;
    cout << "Usage: cplx2tiff <input cplx file> <output tiff file> "
	 << "<component type> <dim x> <dim y>" << endl;
    cout << "  where component type is one of:" <<endl;
    cout << "     0 - REAL" << endl;
    cout << "     1 - IMAG" << endl;
    cout << "     2 - MAG" << endl;
    cout << "     3 - PHASE" << endl;
    cout << "     4 - MAG_SQ" << endl;
    return 1;
  }

  //read the data block in the file
  int nx = atoi(argv[4]);
  int ny = atoi(argv[5]);
  Complex_2D complex(nx,ny);
  int status;

  //read the data into an array
  status = read_cplx(argv[1], complex);
  
  if(!status){
    cout << "failed.. exiting"  << endl;
    return(1);
  }

  Double_2D real(nx,ny);

  complex.get_2d(atoi(argv[3]),real);
  
  //write the data to a file
  status = write_tiff(argv[2], real);

      return 0;
}
