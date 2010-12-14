#include <iostream>
#include <stdlib.h>
#include "io.h"

using namespace std;

/**************************************/
int main(int argc, char * argv[]){

  //check for 5 arguements
  if(argc!=5 ){
    cout << "Wrong number of arguments. Usage: " 
	 << "dbin2ppm <input dbin file> <output ppm file> <dim x> <dim y>" << endl;
    return 1;
  }

  //read the data block in the file
  double ** data;
  int nx = atoi(argv[3]);
  int ny = atoi(argv[4]);
  int status;

  //read the data into an array

  status = read_dbin(argv[1], nx, ny, &data);
  
  if(!status){
    cout << "failed.. exiting"  << endl;
    return(1);
  }
  
  //write the data to a file
  write_ppm(argv[2], nx, ny, data, true);
      
  return 0;
}