#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <cmath>
#include "tiffio.h"
#include "io.h"
#include "Double_2D.h"

using namespace std;

#define FAILURE 0
#define SUCCESS 1

#define ONE_BYTE 1
#define TWO_BYTE 2
#define FOUR_BYTE 4

/***************************************************************/
//bad code which allows the array type to be 
//decided at run-time.
/***************************************************************/
class anonomous_array{

  int type_;
  uint8  * a_uint8;
  uint16 * a_uint16;
  uint32 * a_uint32;
  
 public:

  anonomous_array(int type, int size){
    type_ = type;

     a_uint8=0;
     a_uint16=0;
     a_uint32=0;

    switch( type ){
    case ONE_BYTE:
      a_uint8 = new uint8[size];
      break;
    case TWO_BYTE:
      a_uint16 = new uint16[size];
      break;
    case FOUR_BYTE:
      a_uint32 = new uint32[size];
      break;
      
    default:
      cout << "Not familiar with type.. exiting.." << endl;
    }    
  }
  
  ~anonomous_array(){
    switch( type_ ){
    case ONE_BYTE:
      delete[] a_uint8;
      break;
    case TWO_BYTE:
      delete[] a_uint16;
      break;
    case FOUR_BYTE:
      delete[] a_uint32;
      break;

    default:
      cout << "Not familiar with the HDF4 type.. exiting.." << endl;
    }    
  }

  void * return_array(){
    switch( type_ ){
    case ONE_BYTE:
      return a_uint8;
    case TWO_BYTE:
      return a_uint16;
    case FOUR_BYTE:
      return a_uint32;
     default:
      cout << "Not familiar with the HDF4 type.. exiting.." << endl;
      return FAILURE;
    }  
  }

  double return_array_value(int i){

    //coverting to double. Not ideal.
    switch( type_ ){
    case ONE_BYTE:
      return a_uint8[i];
    case TWO_BYTE:
      return a_uint16[i];
    case FOUR_BYTE:
      return a_uint32[i];
    default:
      cout << "Not familiar with the HDF4 type.. exiting.." << endl;
      return FAILURE;
    }          
  }


  
};

/***************************************************************/

/***************************************************************/

int read_tiff(string file_name, Double_2D & data){
 
  //open the input file:
  TIFF* tif = TIFFOpen(file_name.c_str(), "r");
  
  if (tif) {
    int dircount = 0;
    for( ; TIFFReadDirectory(tif); dircount++);
    if(dircount > 1)
      cout << "Multiple directories in the file "<<file_name
	   << "... only using the first" <<endl;
  }
  else{
    cout << "Could not open the file "<<file_name<<endl;
    return FAILURE;
  }
  
  uint32 w, h;
  uint16 bits_per_sample;
  uint16 samples_per_pixel;

  TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
  TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
  TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bits_per_sample);
  TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &samples_per_pixel);

  int strip_size = TIFFStripSize(tif);
  int bytes_per_pixel = bits_per_sample*samples_per_pixel/8;
  int pixels_per_strip = strip_size/bytes_per_pixel;
  
  double * grey_image = new double[w*h]; 
  uint32 * colour_image = new uint32[w*h];
  
  if(samples_per_pixel>1){ //see if the image is colour

    cout << "Processing colour image" << endl;
    TIFFReadRGBAImage(tif, w, h, colour_image, 0);
  }
  else{ //otherwise if the image is grey scale

    cout << "Processing grey scale image" << endl;

    anonomous_array * buffer;

    switch(bytes_per_pixel){ //otherwise set up the buffer for grey scale
    case(ONE_BYTE):
      buffer = new anonomous_array(ONE_BYTE,pixels_per_strip);
    break;
    case(TWO_BYTE):
      buffer = new anonomous_array(TWO_BYTE,pixels_per_strip);
    break;
    case(FOUR_BYTE):
      buffer = new anonomous_array(FOUR_BYTE,pixels_per_strip);
      break;
    default:
      cout << "Confused about the tiff image.." <<endl;
      return FAILURE;
    }
  

    for ( tstrip_t strip = 0; strip < TIFFNumberOfStrips(tif); strip++){
      int bytes_read = TIFFReadEncodedStrip(tif, strip, buffer->return_array(), (tsize_t) - 1);
      for( int i=0; i< bytes_read/bytes_per_pixel; i++ ){
	grey_image[i+strip*pixels_per_strip] = buffer->return_array_value(i);
      }
    }
    
  }

  //copy to the image array
  if(data.get_size_x()==0){
    cout << "allocating memory"<<endl;
    data.allocate_memory(w,h);
  }

  for(int i=0; i < w; ++i){
    for(int j=0; j< h; ++j){
      if(samples_per_pixel>1){//if the image is colour we take the sum of colour value
	int pixel = colour_image[(h-j)*w+i];
	data.set(i,j,TIFFGetR(pixel)+TIFFGetG(pixel)+TIFFGetB(pixel));
      }
      else //grey scale:
	data.set(i,j,grey_image[j*w+i]);
    }
  }
  
  delete[] grey_image;
  delete[] colour_image;

  TIFFClose(tif);
  
  return SUCCESS; //success
    
}
