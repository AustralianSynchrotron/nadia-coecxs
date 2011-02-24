#include <string>

class Complex_2D;
class Double_2D;

/**
 * Read a ppm file. Returns a 2D of the data. 
 *
 * @param file_name The name of the file to read from
 * @param data The array to be filled with data
 */
int read_ppm(std::string file_name, Double_2D & data);

/**
 * Read a tiff file. Returns a 2D of the data.
 * 
 * @param file_name The name of the file to read from
 * @param data The array to be filled with data
 */
int read_tiff(std::string file_name, Double_2D & data);

/**
 * Read a HDF4 file. Returns a 2D of the data. 
 *
 * @param file_name The name of the file to read from
 * @param data The array to be filled with data
 * @param data_name The name of the block in the HDF4 where the data is.
 * By default it looks for the "data" block.
 */
int read_hdf4(std::string file_name, Double_2D & data, 
	      char * data_name="data");



int read_dbin(std::string file_name, int nx, int ny, Double_2D & data);

int read_cplx(std::string file_name, Complex_2D & complex);



/** 
 * Write a 2D array to a ppm file
 *
 * @param file_name The name of the file to write to
 * @param data The array to be written to file
 * @param log_scale Output on log scale? true/false. Default is false.
 */ 
int write_ppm(std::string file_name, const Double_2D & data, 
	      bool log_scale=false);


int write_dbin(std::string file_name, const Double_2D & data);


int write_cplx(std::string file_name, const Complex_2D & complex);


int write_tiff(std::string file_name, const Double_2D & data);
