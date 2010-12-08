#include <iostream>
#include <stdlib.h>
#include "io.h"
#include "Complex_2D.h"

using namespace std;

/**************************************/
int main(int argc, char * argv[]){

  //check for 5 arguements
  if(argc!=6 ){
    cout << "Wrong number of arguments. Usage: " 
	 << "dbin2ppm <input dbin file> <output real part ppm file> "
	 << "<output imag part ppm file> <dim x> <dim y>" << endl;
    return 1;
  }

  //read the data block in the file
  int nx = atoi(argv[4]);
  int ny = atoi(argv[5]);
  Complex_2D complex(nx,ny);
  int status;

  //read the data into an array
  status = read_cplx(argv[1], nx, ny, &complex);
  
  if(!status){
    cout << "failed.. exiting"  << endl;
    return(1);
  }

  double ** mag = new double*[nx];
  double ** phase = new double*[nx];
  for(int i=0; i < nx; i++){
    mag[i]= new double[ny];
    phase[i]= new double[ny];
  }
  complex.get_2d(MAG,&mag);
  complex.get_2d(PHASE,&phase);
  
  //write the data to a file
  write_ppm(argv[2], nx, ny, mag);
  write_ppm(argv[3], nx, ny, phase);
      
  for(int i=0; i < nx; i++){
    delete [] mag[i];
    delete [] phase[i];
  }
  delete [] mag;
  delete [] phase;
  
  return 0;
}
