// Copyright 2011 Nadia Davidson for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

/**
 * @file hdf2dbin.c
 * 
 * \a hdf2dbin.exe - Extract the data from a HDF4 file and save as 
 * a dbin Convert a binary file (2D double / 64 bit, format). 
 * 
 * Usage: hdf2dbin.exe \<input hdf file\> \<output dbin file\> \<optional: data block name\> \par 
 *
 * \par Example:
 * \verbatim  hdf2dbin.exe my_data.hdf my_data.dbin  \endverbatim
 * 
 **/

#include <iostream>
#include <stdlib.h>
#include <Double_2D.h>
#include <io.h>

using namespace std;

/**************************************/
int main(int argc, char * argv[]){

  //check for 2 or 3 arguements
  if(argc<3 || argc>4 ){
    cout << "Wrong number of arguments. Usage: " 
	 << "hdf2dbin <input hdf file> <output dbin file>" << endl;
    return 1;
  }

  //read the data block in the file
  Double_2D data;
  int status;

  //read the data into an array
  if(argc==4)
    status = read_hdf4(argv[1], data, argv[3]);
  else
    status = read_hdf4(argv[1], data);
  
  if(!status){
    cout << "failed.. exiting"  << endl;
    return(1);
  }
  
  //write the data to a file
  write_dbin(argv[2], data);
      
  return 0;
}
