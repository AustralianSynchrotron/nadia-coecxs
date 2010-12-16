#include <iostream>
#include <stdlib.h>
#include "io.h"
#include "Complex_2D.h"
#include "Double_2D.h"

using namespace std;

/**************************************/
int main(int argc, char * argv[]){

  //check for 5 arguements
  if(argc!=6 ){
    cout << "Wrong number of arguments. Usage: " 
	 << "dbin2ppm <input dbin file> <name of ppm file> "
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
  write_ppm(argv[2], real);
      
  return 0;
}
