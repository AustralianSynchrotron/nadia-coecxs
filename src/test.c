#include <iostream>
#include <stdlib.h>
#include "io.h"

using namespace std;

/**************************************/
int main(int argc, char * argv[]){

  //check for 2 or 3 arguements
  if(argc<3 || argc>5 ){
    cout << "Wrong number of arguments" << endl;
    return 1;
  }

  //read the data block in the file
  double ** data;
  int nx, ny;

  //read the data into an array
  int status = read_ppm(argv[1], &nx, &ny, &data);
  if(!status){
    cout << "failed.. exiting"  << endl;
    return(1);
  }
  
  //write the data to a file
  write_ppm(argv[2], nx, ny, data);
      
  return 0;
}
