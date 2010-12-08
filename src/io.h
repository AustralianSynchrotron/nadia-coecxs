#include <string>

class Complex_2D;

/** Write a 2D array to a ppm file */
int write_ppm(std::string file_name, int nx, int ny, 
	      double ** data, bool log_scale=false);

/** Read a ppm file. Returns a 2D of the data. */
int read_ppm(std::string file_name, int * nx, int  * ny, double *** data);

/** Read a tiff file. Returns a 2D of the data. */
int read_tiff(std::string file_name, int * nx, int  * ny, double *** data);

/** Read a HDF4 file. Returns a 2D of the data. */
int read_hdf4(std::string file_name, int * nx, int  * ny, double *** data, 
	      char * data_name="data");


int read_dbin(std::string file_name, int nx, int ny, double *** data);


int read_cplx(std::string file_name, int nx, int ny, 
	      Complex_2D * complex);
