#include <iostream>
#include <stdlib.h>
#include "io.h"

using namespace std;

/**************************************/
int main(int argc, char * argv[]){

  //check for 2 or 3 arguements
  if(argc<3 || argc>4 ){
    cout << "Wrong number of arguments. Usage: " 
	 << "hdf2ppm <input hdf file> <output ppm file>" << endl;
    return 1;
  }

  //read the data block in the file
  double ** data;
  int nx, ny;
  int status;

  //read the data into an array
  if(argc==4)
    status = read_hdf4(argv[1], &nx, &ny, &data, argv[3]);
  else
    status = read_hdf4(argv[1], &nx, &ny, &data);
  
  if(!status){
    cout << "failed.. exiting"  << endl;
    return(1);
  }
  
  //write the data to a file
  write_ppm(argv[2], nx, ny, data);
      
  return 0;
}
