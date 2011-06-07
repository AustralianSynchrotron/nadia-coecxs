#include "utils.h"
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include "io.h"
#include "Complex_2D.h"
#include "Double_2D.h"
#include "FresnelCDI.h"
#include <cstdlib>
#include <math.h>
#include <string>

using namespace std;

#define BINS 100
#define sgn(x) (( x > 0 ) - ( x < 0 ))

double sgrad(Double_2D & image, int i, int j, double & dir,
	     double x=100, double y=100){

  if(i<1 || j<1 
     || i > (image.get_size_x()-2)
     || j > (image.get_size_y()-2))
    return 0;

      //do the x direction
  double value_x  = -1*image.get(i-1,j-1);
  value_x -=  2*image.get(i-1,j);
  value_x -=  1*image.get(i-1,j+1);

  value_x +=  1*image.get(i+1,j-1);
  value_x +=  2*image.get(i+1,j);
  value_x +=  1*image.get(i+1,j+1);
  
  //do the y direction
  double value_y  = -1*image.get(i-1,j-1);
  value_y -=  2*image.get(i,j-1);
  value_y -=  1*image.get(i+1,j-1);

  value_y +=  1*image.get(i-1,j+1);
  value_y +=  2*image.get(i,j+1);
  value_y +=  1*image.get(i+1,j+1);
  
  dir = atan2(value_y,value_x);

  if(x==100&&y==100)
    return sqrt(value_x*value_x + value_y*value_y);
  else 
    return x*value_x + y*value_y;
}

void crop(Double_2D & image, Double_2D & new_image, int x_start, int y_start){
  int nx = new_image.get_size_x();
  int ny = new_image.get_size_y();

  if( image.get_size_x() - x_start < nx || image.get_size_y() - y_start < ny){
    cout << "In crop(), the new image given is too small "
	 << "to hold the content of the cropped image .. exiting"
	 << endl;
    exit(1);
  }
  
  for(int i=0; i < nx ; i++){
    for(int j=0; j < ny ; j++){
      new_image.set(i,j,image.get(x_start+i,y_start+j));
    }
  }

}

void rescale(Double_2D & image, double scale){

  int nx = image.get_size_x();
  int ny = image.get_size_y();

  int middle_x = nx/2;
  int middle_y = ny/2;

  Double_2D image_temp(nx,ny);

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      
      double i_ = scale*(i-middle_x) + middle_x;
      double j_ = scale*(j-middle_y) + middle_y;

      //this is the tricky bit.
      //use a weighted average
      double total_weights=0;
      double weighted_value=0;

      for(double x = (int) i_; x < i_ + scale; x++ ){
	for(double y = (int) j_; y < j_ + scale; y++){

	  //check that we are still in range
	  if(x<nx&&x>=0&&y<ny&&y>=0){

	    //subtract the remainders from one.
	    double x_fraction = 1;
	    double y_fraction = 1;

	    if(x < i_)
	      x_fraction = 1 - fmod(i_, 1.0);

	    if(y < j_)
	      y_fraction = 1 - fmod(j_, 1.0);

	    if(x > i_ + scale -1)
	      x_fraction = fmod((i_+scale), 1.0);

	    if(y > j_ + scale - 1)
	      y_fraction = fmod((j_+scale), 1.0);

	    //add to the weight
	    double this_weight = x_fraction*y_fraction;

	    total_weights += this_weight;
	    weighted_value += this_weight*image.get(x,y);

	  }
	}
      }

      if(total_weights==0)
	image_temp.set(i,j,0);
      else
	image_temp.set(i,j,weighted_value/total_weights);      
    }
  }

  image.copy(image_temp);

}


//no good
double calculate_high_frequency_ratio(Double_2D & image){

  double value = 0;
  double total = 0;

  int nx = image.get_size_x();
  int ny = image.get_size_y();

  double * in = new double[nx*ny];
  double * out = new double[nx*ny];

  fftw_plan fftw = fftw_plan_r2r_2d(nx,ny, in, out,
				    FFTW_REDFT00, FFTW_REDFT00,
				    FFTW_ESTIMATE);
  

  //copy image to the in array
  for(int i=0; i < nx ; i++){
    for(int j=0; j < ny ; j++){
      in[i*ny + j] = image.get(i,j);
    }
  }

  fftw_execute(fftw);

  Double_2D output(nx,ny);
  //copy image to the in array
  for(int i=0; i < nx ; i++){
    for(int j=0; j < ny ; j++){
      output.set(i,j,fabs(out[i*ny + j]));
      total+=fabs(out[i*ny + j]);
    }
  } 
  write_image("fftw.ppm",output,true);

  value = output.get(0,0);
  //value += output.get(0,1);
  //value = output.get(1,0);
  //value += output.get(1,1);

  fftw_destroy_plan(fftw);  
  delete[] in;
  delete[] out;

  return value;

}

//calculate the chi2 between different images
double diff_of_squares(Double_2D & image1, Double_2D & image2){
  int nx = image1.get_size_x();
  int ny = image1.get_size_y();
  
  int sum = 0;
  double total = 0;

  for(int i=1; i < nx ; i++){
    for(int j=1; j < ny ; j++){
      sum+=pow(image1.get(i,j)-image2.get(i,j),2);
    }
  }

  return sqrt(sum);
  
}

double count_pixels(Double_2D & image, double threshold){
  int nx = image.get_size_x();
  int ny = image.get_size_y();
  
  int sum = 0;
  double total = 0;

  for(int i=1; i < nx ; i++){
    for(int j=1; j < ny ; j++){
      total+=image.get(i,j);
    }
  }

  double scale = (nx*ny)/total;
  
  for(int i=1; i < nx ; i++){
    for(int j=1; j < ny ; j++){
      if(image.get(i,j)>threshold)
	sum++;
    }
  }

  return sum;
}

//same as calculate_average_energy_density
double deviation_from_zero(Double_2D & image){
  int nx = image.get_size_x();
  int ny = image.get_size_y();

  double sum = 0;
  double total = 0;

  for(int i=1; i < nx ; i++){
    for(int j=1; j < ny ; j++){
      total+=image.get(i,j);
    }
  }

  double scale = (nx*ny)/total;


  for(int i=1; i < nx ; i++){
    for(int j=1; j < ny ; j++){
      
      double temp = (scale*image.get(i,j)+1);
      sum+=temp*temp;
      //      sum += log10(temp);
    }
  }

  return sum; //log10(nx*ny);

}


//no good for focal-sample distance optimisation
//good for normalisation optimisation
double calculate_average_energy_density(Double_2D & image){
  int nx = image.get_size_x();
  int ny = image.get_size_y();

  double sum = 0;
  double total = 0;
  double value = 0;

  for(int i=1; i < nx ; i++){
    for(int j=1; j < ny ; j++){
      total+=image.get(i,j);
    }
  }
  
  double scale = (nx*ny)/total;
  total=0;

  for(int i=1; i < nx-1 ; i++){
    for(int j=1; j < ny-1 ; j++){

      sum=scale*image.get(i,j);
      sum*=scale*image.get(i-1,j);
      sum*=scale*image.get(i+1,j);
      sum*=scale*image.get(i,j-1);
      sum*=scale*image.get(i,j+1);
      
      total += scale*image.get(i,j);
      value += sum;
    }
  }

  return value;

}


//
double simple(Double_2D & image, double scale){
  int nx = image.get_size_x();
  int ny = image.get_size_y();
  
  Double_2D image_copy(nx,ny);
  image_copy.copy(image);
  image_copy.scale(scale);

  double total = 0;

  //loop over the array once to sort into bins
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      total+=pow(image_copy.get(i,j),2);
    }
  }
  
  return sqrt(total)/(nx*ny);
}

//no good
double calculate_image_entropy(Double_2D & image){
  
  double max = image.get_max();
  double min = image.get_min();
  double counts[BINS]={0};
  int bin;

  int nx = image.get_size_x();
  int ny = image.get_size_y();

  double total = nx*ny;
  //cout << "here" <<endl;

  //loop over the array once to sort into bins
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      //cout << "i,j ";
      if(image.get(i,j)==max)
	bin = BINS-1;
      else
	bin = BINS*((image.get(i,j)-min)/(max-min));
      //cout << "bin="<<bin;
      counts[bin]+=1.0/total;
      //cout << " count="<<counts[bin]<<endl;
    }
  }

  //cout << "here2" <<endl;

  //loop over the bins to get the entropy
  double entropy = 0;

  //  double max = 

  //  cout << endl << "a = [";
  for(int c=1; c < BINS; c++){

    /**    cout << counts[c];
    if(c!=BINS-1)
      cout << ",";
    else
    cout << "]"<<endl;**/
      
    if(counts[c]!=0){
      //      cout << " log2(counts[c]) = " <<log2(counts[c])<<endl;
      entropy-= counts[c]*log2(counts[c]);
    }
  }
  
  //  cout << "Zeros: "<< counts[0] << endl;
  return counts[2]; //entropy;

}


double calculate_image_entropy_2(Double_2D & image){
  
  double max = image.get_max();
  double min = image.get_min();
  double counts[BINS][BINS]={0};
  int bin_x;
  int bin_y;

  int nx = image.get_size_x();
  int ny = image.get_size_y();

  double total = 2*(nx*ny-1);

  //loop over the array once to sort into bins
  for(int i=1; i<nx; i++){
    for(int j=1; j<ny; j++){

      if(image.get(i,j)==max)
	bin_x = BINS-1;
      else
	bin_x = BINS*((image.get(i,j)-min)/(max-min));

      if(image.get(i-1,j)==max)
	bin_y = BINS-1;
      else
	bin_y = BINS*((image.get(i-1,j)-min)/(max-min));

      counts[bin_x][bin_y]+=1.0/total;

      if(image.get(i,j-1)==max)
	bin_y = BINS-1;
      else
	bin_y = BINS*((image.get(i,j-1)-min)/(max-min));

      counts[bin_x][bin_y]+=1.0/total;

      //      cout << "x,y="<<bin_x<<","<<bin_y<<" counts="<<counts[bin_x][bin_y]<< endl;

    }
  }

  //loop over the bins to get the entropy
  double entropy = 0;

  for(int c_x=0; c_x < BINS; c_x++){
    for(int c_y=0; c_y < c_x+1; c_y++){
      
      double c = counts[c_x][c_y] + counts[c_y][c_x];


      
      if(c!=0.0&&c_x!=0&&c_y!=0){
	//cout << "cx="<<c_x<<" cy="<<c_y<<" c="<<c<<endl;
	entropy-= c*log2(c);
      }
      
    }
  }

  return entropy;
}


double laplace_gradient(Double_2D & image){

  double total = 0;
  double value ;

  int nx = image.get_size_x();
  int ny = image.get_size_y();

  // double max = image.get_max()-image.get_min();
  //  threshold = pow(threshold*max,2);

  //  double min_allowed = threshold*4*max*max;

  Double_2D output(nx,ny);

  //copy image to the in array
  for(int i=2; i < nx-2 ; i++){
    for(int j=2; j < ny-2 ; j++){

      value = 0;
      
      for(int i_=i-2; i_ < i+3 ; i_++){
	for(int j_=j-2; j_ < j+3 ; j_++){
	  
	  if(i_==i&&j_==j)
	    value  += 24*image.get(i_,j_);
	  else
	    value  -= image.get(i_,j_);      
	}
      }

      output.set(i,j,fabs(value));
      total += fabs(value) ;
    }
  }
  
  char buf[50];
  static int counter = 0;
  cout << output.get_max() << endl;
  sprintf(buf,"laplace_%i.tiff",counter);
  write_image(buf,output);

  counter++;

  /**Double_2D output_2(nx,ny);

  value = 0;
  //  total = 0;
  double total_t = 0;

  //now see how isolated they are.
  for(int i=1; i < nx-1 ; i++){
    for(int j=1; j < ny-1 ; j++){
      value=0;
      //  total+=output.get(i,j);
      if(output.get(i,j)>min_allowed){
	//	total++;
	value+=output.get(i-1,j)/(max*max);
	value+=output.get(i+1,j)/(max*max);
	value+=output.get(i,j-1)/(max*max);
	value+=output.get(i,j+1)/(max*max);

	total_t+=value;
      }

      output_2.set(i,j,value);

    }
  }**/

  return total;

}


void convolve(Double_2D & array, double gauss_width, 
	      int pixel_cut_off){
  //to speed up computation we only convolve 
  //up to 4 pixels away from the gaussian peak

  //make a temporary array to hold the smeared image
  double nx = array.get_size_x();
  double ny = array.get_size_y();

  Double_2D temp_array(nx,ny);
  
  //make a temporary array to hold the gaussian distribution.
  Double_2D gauss_dist(pixel_cut_off+1, pixel_cut_off+1);
  for(int i=0; i <= pixel_cut_off; i++){
    for(int j=0; j <= pixel_cut_off; j++){
      double denom = 2.0*gauss_width*gauss_width;
      gauss_dist.set(i,j,exp(-1*(i*i+j*j)/denom ) );
    }
  }      

  //now do the convolution
  //this is messy. First loop over the elements of
  //the array which was given as input
  double new_value;
  for(int i=0; i < nx; i++){
    for(int j=0; j < ny; j++){
      
      //now loop over the colvoluted array (the one we want to make).
      //Calculate the contribution to each element in it.
      
      new_value = 0;
      
      for(int i2=i-pixel_cut_off; i2 <= i+pixel_cut_off; i2++){
	for(int j2=j-pixel_cut_off; j2 <= j+pixel_cut_off; j2++){
	  if(i2<nx && i2>=0 && j2>=0 && j2<ny){
	    new_value += array.get(i2,j2)*gauss_dist.get(fabs(i-i2),fabs(j-j2));
	  }
	}
      }
      temp_array.set(i,j,new_value); 
    }
  }

  array.copy(temp_array);

}


double sobel_gradient(Double_2D & image){

  int nx = image.get_size_x();
  int ny = image.get_size_y();

  Double_2D mask(nx,ny);
  mask.copy(image);
  mask.scale(1.0/image.get_sum());
  
  convolve(mask,3.0,9);

  char buf[50];
  static int counter = 0;
  //sprintf(buf,"image_conv_%i.tiff",counter);
  //write_image(buf,image); 

  double value_x ;
  double value_y;

  // double max = image.get_max()-image.get_min();
  //  threshold = pow(threshold*max,2);

  //  double min_allowed = threshold*4*max*max;

  Double_2D output(nx,ny);
  Double_2D direction(nx,ny);

  double smallest_max = 0;
  double no_of_max = 2;
  double max_total = 0;

  for(int i=2; i < nx-2 ; i++){
    for(int j=2; j < ny-2 ; j++){

      //do the x direction
      value_x  = -1*mask.get(i-1,j-1);
      value_x -=  2*mask.get(i-1,j);
      value_x -=  1*mask.get(i-1,j+1);

      value_x +=  1*mask.get(i,j-1);
      value_x +=  2*mask.get(i,j);
      value_x +=  1*mask.get(i,j+1);

      //do the y direction
      value_y  = -1*mask.get(i-1,j-1);
      value_y -=  2*mask.get(i,j-1);
      value_y -=  1*mask.get(i+1,j-1);

      value_y +=  1*mask.get(i-1,j);
      value_y +=  2*mask.get(i,j);
      value_y +=  1*mask.get(i+1,j);

      double value_total = sqrt(value_x*value_x + value_y*value_y);
      output.set(i,j,value_total);
    }

  }

  
  double max = output.get_max();
  double cut = 0.2*max;
  double total = 0;

  for(int i=2; i < nx-2 ; i++){
    for(int j=2; j < ny-2 ; j++){

	//do the x direction
	value_x  = -1*image.get(i-1,j-1);
	value_x -=  2*image.get(i-1,j);
	value_x -=  1*image.get(i-1,j+1);
      
	value_x +=  1*image.get(i+1,j-1);
	value_x +=  2*image.get(i+1,j);
	value_x +=  1*image.get(i+1,j+1);

	//do the y direction
	value_y  = -1*image.get(i-1,j-1);
	value_y -=  2*image.get(i,j-1);
	value_y -=  1*image.get(i+1,j-1);
	
	value_y +=  1*image.get(i-1,j+1);
	value_y +=  2*image.get(i,j+1);
	value_y +=  1*image.get(i+1,j+1);
	  
	//value_x =  image.get(i-1,j-1)-image.get(i,j);
	//value_y =  image.get(i,j-1)-image.get(i-1,j);

	double value_total = sqrt(value_x*value_x + value_y*value_y);
	
	if(output.get(i,j)>cut){
	  mask.set(i,j,value_total);
	  direction.set(i,j,atan2(value_y,value_x));
	}
	else{
	  mask.set(i,j,0);
	  direction.set(i,j,0);
	}
	
	total+=value_total;
    }
  }



  //  cout << "count is: " << number << endl;


  //  cout << "max is: " << max << endl;
  //  cout << "normalised max: " << max_total/image.get_sum() << endl;

  sprintf(buf,"sobel_direction_%i.tiff",counter);
  write_image(buf,direction);
  counter++;

  //  cout << "Entropy of gradient is :"<< calculate_image_entropy(output)<<endl;

  return total; //total_t; // output.get_max();

}



//good for focal-sample distance!!
//not so good for normalisation.
double calculate_gradients(Double_2D & image, double threshold){

  double value = 0;
  double total = 0;

  //  image.scale(image.get_sum());

  int nx = image.get_size_x();
  int ny = image.get_size_y();

  double max = image.get_max()-image.get_min();
  //threshold = pow(threshold*max,2);

  double min_allowed = threshold*2*max*max;

  max = 0;

  Double_2D output(nx,ny);
  //copy image to the in array
  for(int i=1; i < nx-1 ; i++){
    for(int j=1; j < ny-1 ; j++){
      /**      value =pow(fabs(image.get(i,j)-image.get(i-1,j)),2);      
      value+=pow(fabs(image.get(i,j)-image.get(i,j-1)),2);
      value+=pow(fabs(image.get(i,j)-image.get(i+1,j)),2);      
      value+=pow(fabs(image.get(i,j)-image.get(i,j+1)),2);**/

      value =fabs( pow(image.get(i,j)-image.get(i-1,j),2 ) );     
      //		   pow(image.get(i,j)-image.get(i+1,j),3) );
      
      value+= fabs( pow(image.get(i,j)-image.get(i,j-1),2 ) );     
		    //		    pow(image.get(i,j)-image.get(i,j+1),3) );

      output.set(i,j,fabs(value));
      total+=value;

      if(value> max)
	max = value;
    }
  }

  Double_2D output_2(nx,ny);

  value = 0;
  total = 0;
  double total_t = 0;

  //now see how isolated they are.
  for(int i=1; i < nx-1 ; i++){
    for(int j=1; j < ny-1 ; j++){
      value=0;
      //  total+=output.get(i,j);
      if(output.get(i,j)>min_allowed){
	total++;
	value+=output.get(i-1,j);///(max);
	value+=output.get(i+1,j);///(max);
	value+=output.get(i,j-1);///(max);
	value+=output.get(i,j+1);///(max);

	/**	if(output.get(i-1,j)>min_allowed)
	  value++;
	if(output.get(i+1,j)>min_allowed)
	  value++;
	if(output.get(i,j-1)>min_allowed)
	  value++;
	if(output.get(i,j+1)>min_allowed)
	value++; **/

	total_t+=value;
      }
      
      output_2.set(i,j,value);
      
    }
  }

  char buf[50];
  static int counter = 0;
  sprintf(buf,"grad_%i.tiff",counter);
  write_image(buf,output_2);
  counter++;

  return total; //total_t;///total;

}

//no good
double vollaths_4(Double_2D & image){
  int nx = image.get_size_x();
  int ny = image.get_size_y();

  double sum = 0;

  for(int i=2; i < nx ; i++){
    for(int j=0; j < ny ; j++){
      sum += image.get(i,j)*image.get(i-1,j)
	-image.get(i,j)*image.get(i-2,j);
    }
  }

  return sum;
}

//no good
double vollaths_5(Double_2D & image){

  int nx = image.get_size_x();
  int ny = image.get_size_y();

  double sum_x = 0;
  double sum_y = 0;
  double total = 0;

  for(int i=1; i < nx ; i++){
    for(int j=1; j < ny ; j++){
      sum_x += image.get(i,j)*image.get(i-1,j);
      sum_y += image.get(i,j)*image.get(i,j-1);
    }
  }

  return sum_x*sum_y/(nx*nx*ny*ny); //sum_x*sum_y*total/(pow(nx*ny,3));
}

double line_out(Double_2D & image){

  int nx = image.get_size_x();
  int ny = image.get_size_y();

  int y = ny/2.0;
  int x_start = nx/4.0;
  int x_end = 3*nx/4.0;

  double value;
  
  for(int x=x_start; x < x_end+1 ; x++){
    value = 0;
    if(x > x_start)
      value = fabs(image.get(x,y)-image.get(x-1,y));

    //    if(x>(x_start+1))
    //  value = fabs(image.get(x,y)*image.get(x-1,y)
    //		   -image.get(x,y)*image.get(x-2,y));
    value*=value;
    cout << ", " << value;
  }
  cout << endl; 
  
  return 0;  
}

double edge_grad(Double_2D & image, Double_2D & mask){

  int nx = image.get_size_x();
  int ny = image.get_size_y();
  
  double value = 0;
  double edge_points = 0;

  Double_2D direction(nx,ny);

  for(int x=1 ; x < nx ; x++){
    for(int y=1 ; y < ny ; y++){
      
      if(mask.get(x,y)!=0){
	double dir;
	double this_value = sgrad(image,x,y,dir);
	double dir_m = mask.get(x,y);
	double projection = cos(dir)*cos(dir_m)+sin(dir)*sin(dir_m);

	//	cout << "value="<<this_value<<" dir_m="<<dir_m<<"projection="<< projection<<endl;

	value += this_value; //projection*this_value;
	edge_points++;  //+= projection;
	direction.set(x,y,dir);
      }
    }
  }

  static int counter = 0;
  char buf[80];
  sprintf(buf,"edges_grad_%i.tiff",counter);
  write_image(buf,direction);
  counter++;

  //  cout << "value="<<value<<" edge_points="<<edge_points<<endl;

  return value/((double) edge_points);
}


double edges(Double_2D & image){

  int nx = image.get_size_x();
  int ny = image.get_size_y();
  
  Double_2D temp(nx,ny);
  
  double threshold = 0.05*image.get_max();
  double value = 0;
  double edge_points = 0;
  double direction = 0;

  for(int x=1 ; x < nx ; x++){

    bool edge_found = false;
    
    for(int y=1 ; y < ny ; y++){

      if(edge_found==false && image.get(x,y)>threshold ){
	edge_found = true;
	//image.set(x,y,60000);
	if(temp.get(x,y)==0.0){
	  value += sgrad(image,x,y,direction); //fabs(image.get(x,y)-image.get(x-1,y));
	  temp.set(x,y,M_PI/2.0);
	  edge_points++;
	}
      }
    }

    edge_found = false;
    
    for(int y=ny-1 ; y > 0 ; y--){
      
      //if(edge_found==true && image.get(x,y)<threshold){ 
      //edge_found = false;
      if(edge_found==false && image.get(x,y)>threshold ){
      	edge_found = true;
	//image.set(x,y,60000);
	if(temp.get(x,y)==0.0){
	  value += sgrad(image,x,y,direction); //fabs(image.get(x,y)-image.get(x-1,y));
	  temp.set(x,y,-M_PI/2.0);
	  edge_points++;
	}
      }
      
    }
  }
  
  
  for(int y=1 ; y < ny ; y++){
    
    bool edge_found = false;
    
    for(int x=1 ; x < nx ; x++){

      if(edge_found==false && image.get(x,y)>threshold ){
	edge_found = true;
	//image.set(x,y,60000);
	if(fabs(cos(temp.get(x,y)))<0.01){
	  value += sgrad(image,x,y,direction); //fabs(image.get(x,y)-image.get(x,y-1));
	  if(temp.get(x,y)==0)
	    temp.set(x,y,2.0*M_PI);
	  else
	    temp.set(x,y,temp.get(x,y)/2.0);
	  edge_points++;
	}
      }
    }
    
    edge_found = false;
    
    for(int x=nx-1 ; x > 0 ; x--){
      
      if(edge_found==false && image.get(x,y)>threshold ){
      	edge_found = true;
	//if(edge_found==true && image.get(x,y)<threshold ){
	//edge_found = false;
	//image.set(x,y,60000);
	if(fabs(cos(temp.get(x,y)))<0.01){
	  value += sgrad(image,x,y,direction); //fabs(image.get(x,y)-image.get(x,y-1));
	  edge_points++;
	  if(temp.get(x,y)==0)
	    temp.set(x,y,M_PI);
	  else
	    temp.set(x,y,M_PI-(temp.get(x,y)/2.0));
	}
      }
      
    } 
  }

  /**  static int counter = 0;
  char buf[80];
  sprintf(buf,"edges_%i.tiff",counter);
  write_image(buf,temp);
  sprintf(buf,"esw_%i.tiff",counter);
  write_image(buf,image);
  counter++;**/

  //  image.copy(temp);
  
  return value/edge_points;  
}

//no good
double calculate_mean_difference(Double_2D & image){

  double mean = 0;
  double value = 0;

  int nx = image.get_size_x();
  int ny = image.get_size_y();

  Double_2D output(nx,ny);

  //loop once to get the mean
  for(int i=0; i < nx ; i++){
    for(int j=0; j < ny ; j++){
      mean +=image.get(i,j);
    }
  }
  //normalise the mean
  mean=mean/((double)nx*ny);

  //loop again to get the average difference
  for(int i=0; i < nx ; i++){
    for(int j=0; j < ny ; j++){
      value += pow(image.get(i,j)-mean,2);
    }
  }
  return value/((double)nx*ny);
  
}


void align(Double_2D & image1, Double_2D & image2, 
	   int & offset_x, int & offset_y, 
	   int step_size, int min_x, int max_x,
	   int min_y, int max_y){
  
  //divide by 8.... until step_size is 1.

  if(step_size<1.0)
    return;

  //make a low res. version of the image.
  int nx = image1.get_size_x();
  int ny = image1.get_size_y();

  int nx_small = image2.get_size_x();
  int ny_small = image2.get_size_y();

  //make a low res. version of the image.
  int nx_low = nx/step_size;
  int ny_low = ny/step_size;

  int nx_small_low = nx_small/step_size;
  int ny_small_low = ny_small/step_size;

  if(min_x==max_x&&min_y==max_y){
    min_x = 1-nx_low; //1-nx_small_low;
    max_x = nx_low;
    min_y = 1-ny_low; //1-ny_small_low;
    max_y = ny_low;
  }

  Double_2D temp1(nx_low,ny_low);
  Double_2D temp2(nx_small_low,ny_small_low);
  
  for(int x=0; x<nx_low; x++){
    for(int y=0; y<ny_low; y++){
      double value =0;
      for(int i=0; i<step_size ; i++){
	for(int j=0; j<step_size ; j++){
	  value+=image1.get(step_size*x+i,step_size*y+j);
	}
      }
      temp1.set(x,y,value/(step_size*step_size));
    }
  }

  for(int x=0; x<nx_small_low; x++){
    for(int y=0; y<ny_small_low; y++){
      double value =0;
      for(int i=0; i<step_size ; i++){
	for(int j=0; j<step_size ; j++){
	  value+=image2.get(step_size*x+i,step_size*y+j);
	}
      }
      temp2.set(x,y,value/(step_size*step_size));
    }
  }

  //now calculate the correlation for the low res. image
  double max_correlation = 0;

  
  //loop over possible alignments
  for(int x=min_x; x<max_x; x++){
    for(int y=min_y; y<max_y; y++){
      
      double correlation = 0;
      double pixels = 0;

      //loop over the image to calculate the correlation metric
      for(int i=0; i<nx_low; i++){
	for(int j=0; j<ny_low; j++){
      
	  int i_ = i - x/step_size;
	  int j_ = j - y/step_size;
	  
	  if(i_>=0&&i_<nx_small_low&&j_>=0&&j_<ny_small_low){
	    correlation+=temp1.get(i,j)*temp2.get(i_,j_);
	    if(temp1.get(i,j)!=0.0 && temp2.get(i_,j_)!=0.0)
	      pixels++;
	  }
	  
	}
      }
      
      correlation=correlation/((double)pixels);

      if(correlation > max_correlation){
	max_correlation = correlation;
	offset_x = x;
	offset_y = y;
      }
	

    }
  }

  cout << "min="<<min_x<< "max="<<max_x<<endl;
  cout << "max_correlation is :"<<max_correlation <<endl;
  static int counter = 0;
  char buff[90];
  sprintf(buff,"temp1_%i.tiff",counter);
  write_image(buff,temp1);
  sprintf(buff,"temp2_%i.tiff",counter);
  write_image(buff,temp2);
  counter++;

  //now decrease the step size and the search area:
  align(image1, image2, 
	offset_x, offset_y, 
	step_size/2.0, 
	offset_x-2*step_size, offset_x+2*step_size,
	offset_y-2*step_size, offset_y+2*step_size);

  return;
  
}


void align_even_better(Double_2D & first_image, Double_2D & second_image,
		       int & offset_x, int & offset_y,
		       int min_x, int max_x,
		       int min_y, int max_y,
		       Double_2D * first_image_weights,
		       Double_2D * second_image_weights,
		       double overlap_fraction){

  const int MAX_PIXELS = 1024;

  int nx_1 = first_image.get_size_x();
  int ny_1 = first_image.get_size_y();
  int nx_2 = second_image.get_size_x();
  int ny_2 = second_image.get_size_y();

  int nx, ny;
  if(nx_1>nx_2)
    nx = nx_1*2;
  else
    nx = nx_2*2;
  if(ny_1>ny_2)
    ny = ny_1*2;
  else
    ny = ny_2*2;

  //if the image is too big for this to work we'll
  //shrink is down, align the shrunken images using this method,
  //and then apply the shifting alignment method.
  if(nx > MAX_PIXELS*2 || ny > MAX_PIXELS*2){

    //get smaller versions of the images
    Double_2D * small_first = new Double_2D(MAX_PIXELS, MAX_PIXELS);
    Double_2D * small_second = new Double_2D(MAX_PIXELS,MAX_PIXELS);

    Double_2D * small_weight_first = new Double_2D(MAX_PIXELS, MAX_PIXELS);
    Double_2D * small_weight_second = new Double_2D(MAX_PIXELS,MAX_PIXELS);

    shrink(first_image, *small_first);
    shrink(second_image, *small_second);
    shrink(*first_image_weights, *small_weight_first);
    shrink(*second_image_weights, *small_weight_second);

    //apply this alignment method with the smaller images
    align_even_better(*small_first, *small_second,
		      offset_x, offset_y,
		      min_x, max_x,
		      min_y, max_y,
		      small_weight_first,
		      small_weight_second,
		      overlap_fraction);   

    delete small_first;
    delete small_second;

    delete small_weight_first;
    delete small_weight_second;

    cout << "done small alignment" <<endl;

    //apply this shifting alignment method to get the final precision.
    align(first_image, second_image, 
	  offset_x, offset_y, 1 , 
	  -nx/MAX_PIXELS,nx/MAX_PIXELS,
	  -ny/MAX_PIXELS,ny/MAX_PIXELS);

    cout << "done fine adjusting: -" << (nx/MAX_PIXELS) << endl;
    
    return;
  }
  
  cout << "here 0"<<endl;

  //check these boundaries later
  if(min_x==max_x || min_y==max_y){
    min_x = 1-nx/2.0; 
    max_x = nx/2.0;
    min_y = 1-ny/2.0;
    max_y = ny/2.0;
  }

  /**Double_2D image1(nx,ny);
  Double_2D image2(nx,ny);

  Double_2D weights1(nx,ny);
  Double_2D weights2(nx,ny);

  cout << "here 1"<<endl; **/

  /**  for(int i=0.5*(nx-nx_1); i < 0.5*(nx+nx_1); i++){
    for(int j=0.5*(ny-ny_1); j < 0.5*(ny+ny_1); j++){
      int i_ = i-0.5*(nx-nx_1);
      int j_ = j-0.5*(ny-ny_1);
      
      image1.set(i,j,first_image.get(i_,j_));
      if(!first_image_weights || first_image_weights->get(i_,j_)>0)
	weights1.set(i,j,1.0);

      if(i_<nx_2&&j_<ny_2){
	
	image2.set(i,j,second_image.get(i_,j_));
	
	if(!second_image_weights || second_image_weights->get(i_,j_)>0)
	  weights2.set(i,j,1.0);
      }
    }
    }**/

  cout << "here 2"<<endl;

  double * temp_img_1 = new double[nx*ny];
  double * temp_img_2 = new double[nx*ny];

  double * temp_img_1_weight = new double[nx*ny];
  double * temp_img_2_weight = new double[nx*ny];

  fftw_complex* temp_fft_1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny); 
  fftw_complex* temp_fft_2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny); 


  fftw_plan fftw1 = fftw_plan_dft_r2c_2d(nx,ny, temp_img_1, temp_fft_1,
					 FFTW_ESTIMATE);
  fftw_plan fftw2 = fftw_plan_dft_r2c_2d(nx,ny, temp_img_2, temp_fft_2,
					 FFTW_ESTIMATE);

  fftw_plan fftw1_w = fftw_plan_dft_r2c_2d(nx,ny, temp_img_1_weight, temp_fft_1,
					 FFTW_ESTIMATE);
  fftw_plan fftw2_w = fftw_plan_dft_r2c_2d(nx,ny, temp_img_2_weight, temp_fft_2,
					 FFTW_ESTIMATE);

  //backwards
  fftw_complex* fft_total = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny); 
  fftw_plan fftw_inv = fftw_plan_dft_c2r_2d(nx,ny, fft_total, temp_img_1,
					    FFTW_ESTIMATE);

  double weight_sum = 0 ;
  
  //copy the image information into the temporary arrays
  for(int i=0; i < nx ; i++){
    for(int j=0; j < ny ; j++){
   
      temp_img_1[i*ny + j] = 0.0;
      temp_img_2[i*ny + j] = 0.0;

      temp_img_1_weight[i*ny + j] = 0.0;
      temp_img_2_weight[i*ny + j] = 0.0;

      if(i>=0.5*(nx-nx_1) && i < 0.5*(nx+nx_1) &&
	 j>=0.5*(ny-ny_1) && j < 0.5*(ny+ny_1)){
	    
	    int i_ = i-0.5*(nx-nx_1);
	    int j_ = j-0.5*(ny-ny_1);
	    
	    //image1.set(i,j,first_image.get(i_,j_));
	    if(!first_image_weights || first_image_weights->get(i_,j_)>0){
	      temp_img_1[i*ny + j] = first_image.get(i_,j_);
	      temp_img_1_weight[i*ny + j] = 1.0;
	    }
	    
	    if(i_<nx_2&&j_<ny_2 && (!second_image_weights || second_image_weights->get(i_,j_)>0)){
	      temp_img_2[i*ny + j] = second_image.get(i_,j_);
	      temp_img_2_weight[i*ny + j] = 1.0;
	    }

	    if(nx_1*ny_1 < nx_2*ny_2)
	      weight_sum+=temp_img_1_weight[i*ny + j];
	    else
	      weight_sum+=temp_img_2_weight[i*ny + j];
	    
	 }
    }
  }

  cout << "here 4"<<endl;

  fftw_execute(fftw1);
  fftw_execute(fftw2);

  //  Complex_2D my_fft(nx,ny);

  cout << "here 5"<<endl;

  //copy image to the in array
  for(int i=0; i < nx ; i++){
    for(int j=0; j < ny ; j++){

      double rl = temp_fft_1[i*ny + j][0]*temp_fft_2[i*ny + j][0]
	+temp_fft_1[i*ny + j][1]*temp_fft_2[i*ny + j][1];
      double im = temp_fft_1[i*ny + j][0]*temp_fft_2[i*ny + j][1]
	-temp_fft_1[i*ny + j][1]*temp_fft_2[i*ny + j][0];

      fft_total[i*ny + j][0] = rl; 
      fft_total[i*ny + j][1] = im; 

      //      temp_img_1[i*ny + j] = weights1.get(i,j);
      //   temp_img_2[i*ny + j] = weights2.get(i,j);

    }
  }

  cout << "here 6"<<endl;

  fftw_execute(fftw1_w);
  fftw_execute(fftw2_w);

  fftw_execute(fftw_inv);

  for(int i=0; i < nx ; i++){
    for(int j=0; j < ny ; j++){

      double rl = temp_fft_1[i*ny + j][0]*temp_fft_2[i*ny + j][0]
	+temp_fft_1[i*ny + j][1]*temp_fft_2[i*ny + j][1];
      double im = temp_fft_1[i*ny + j][0]*temp_fft_2[i*ny + j][1]
	-temp_fft_1[i*ny + j][1]*temp_fft_2[i*ny + j][0];
     
      fft_total[i*ny + j][0] = rl; 
      fft_total[i*ny + j][1] = im; 

      temp_img_2[i*ny + j] = temp_img_1[i*ny + j]/(nx*ny);

    }
  }

  fftw_execute(fftw_inv);

  double max = 0;

  double overlap_factor = weight_sum*nx*ny;

  //copy image to an output  array
  //  Double_2D temp_image(nx,ny); //<-

  for(int i=0; i < nx ; i++){
    for(int j=0; j < ny ; j++){


      if(temp_img_1[i*ny + j]/overlap_factor > overlap_fraction){

	//	temp_image.set(i,j, temp_img_2[i*ny + j]/temp_img_1[i*ny + j]);

	if(temp_img_2[i*ny + j]/temp_img_1[i*ny + j] > max){

	  int temp_i = - i;
	  int temp_j = - j;

	  if(i >=nx/2.0)
	    temp_i = - i + nx;

	  if(j >=  ny/2.0)
	    temp_j = - j + ny;

	  if(temp_i>=min_x && temp_i<max_x &&
	     temp_j>=min_y && temp_j<max_y){
	    
	    max = temp_img_2[i*ny + j]/temp_img_1[i*ny + j];
	    offset_x = temp_i;
	    offset_y = temp_j;
	  }

	  //cout <<temp_img_2[i*ny + j]<<" "<<temp_img_1[i*ny + j]/overlap_factor<<endl;

	}
      }

    }
  }
  
  //  write_image("temp_image.tiff",temp_image,false);
  //cout <<temp_image.get_max()<<endl;

  delete [] temp_img_1;
  delete [] temp_img_2;
  delete [] temp_img_1_weight;
  delete [] temp_img_2_weight;

  fftw_free(temp_fft_1);
  fftw_free(temp_fft_2);
  fftw_free(fft_total);

  fftw_destroy_plan(fftw1);
  fftw_destroy_plan(fftw2);
  fftw_destroy_plan(fftw_inv);
  

}

////////////////////////////////

void interpolate( Double_2D & original, Double_2D & big){

  int nx = original.get_size_x();
  int ny = original.get_size_y();
 
  int bnx = big.get_size_x();
  int bny = big.get_size_y();

  if(bnx % nx !=0 || bny % ny !=0 ){

    cerr<< "In interpolate(Complex_2D & original, Complex_2D & small) " 
	<< "the size of 'big' must be an integer "
	<< "multiple of the size of 'original'" << endl;
    exit(0);
  }

  int scale_x = bnx/nx;
  int scale_y = bny/ny;

  //if we are close to the edge of the image there's no enough
  //information to interpolate, so just get the pixel value to
  //the non-interpolated value.
  //if(lower_i<0 || lower_j<0 || lower_i >= (nx-1) || lower_j >= (ny-1) ){
  //  double value = original.get(i,j);
  //  big.set(new_i+di,new_j+dj,value);
  //  sum+=value;
  // }

  double ** x_pos = new double*[scale_x];
  double ** y_pos = new double*[scale_y];

  //pretabulate the spacings 
  for(int di=0; di < scale_x; di++){
    x_pos[di] = new double[scale_y];
    y_pos[di] = new double[scale_x];
    
    for(int dj=0; dj < scale_y; dj++){
      x_pos[di][dj] = (di + 0.5*((scale_x+1) % 2)) /((double) scale_x);
      y_pos[di][dj] = (dj + 0.5*((scale_y+1) % 2)) /((double) scale_y);
    }
  }

  for(int i=1; i<nx-1; i++){
    for(int j=1; j<ny-1; j++){

      double sum = 0;

      //      int new_i = i*scale_x;
      //      int new_j = j*scale_y;

      double f00=original.get(i,j);
      double f10=original.get(i+1,j);
      double f01=original.get(i,j+1);
      double f11=original.get(i+1,j+1);

      double new_i = (i+0.5)*scale_x ;
      double new_j = (j+0.5)*scale_y ;

      for(int di=0; di < scale_x; di++){
	for(int dj=0; dj < scale_y; dj++){

	  //double x = (di + 0.5*((scale_x+1) % 2)) /((double) scale_x);
	  //double y = (dj + 0.5*((scale_y+1) % 2)) /((double) scale_y);

	  double x = x_pos[di][dj];
	  double y = y_pos[di][dj];

	  double x_1 = 1-x;
	  double y_1 = 1-y;

	  double value = f00*(x_1)*(y_1) + f10*x*(y_1) + f01*(x_1)*y + f11*x*y;
	  big.set(new_i + di , new_j + dj , value);

	}
      }
    }
  }

  for(int di=0; di < scale_x; di++){
    delete [] x_pos[di];
    delete [] y_pos[di];
  }
  delete [] x_pos;
  delete [] y_pos;
  
}




void interpolate( Complex_2D & original, Complex_2D & big){

  int nx = original.get_size_x();
  int ny = original.get_size_y();
 
  int bnx = big.get_size_x();
  int bny = big.get_size_y();

  if(bnx % nx !=0 || bny % ny !=0 ){

    cerr<< "In interpolate(Complex_2D & original, Complex_2D & small) " 
	<< "the size of 'big' must be an integer "
	<< "multiple of the size of 'original'" << endl;
    exit(0);
  }

  int scale_x = bnx/nx;
  int scale_y = bny/ny;


  //if we are close to the edge of the image there's no enough
  //information to interpolate, so just get the pixel value to
  //the non-interpolated value.
  //if(lower_i<0 || lower_j<0 || lower_i >= (nx-1) || lower_j >= (ny-1) ){
  //  double value = original.get(i,j);
  //  big.set(new_i+di,new_j+dj,value);
  //  sum+=value;
  // }

  double ** x_pos = new double*[scale_x];
  double ** y_pos = new double*[scale_x];

  //pretabulate the spacings 
  for(int di=0; di < scale_x; di++){
    x_pos[di] = new double[scale_y];
    y_pos[di] = new double[scale_y];
    
    for(int dj=0; dj < scale_y; dj++){
      x_pos[di][dj] = (di + 0.5*((scale_x+1) % 2)) /((double) scale_x);
      y_pos[di][dj] = (dj + 0.5*((scale_y+1) % 2)) /((double) scale_y);
    }
  }

  for(int i=1; i<nx-1; i++){
    for(int j=1; j<ny-1; j++){

      double sum = 0;

      //      int new_i = i*scale_x;
      //      int new_j = j*scale_y;

      double f00r=original.get_real(i,j);
      double f10r=original.get_real(i+1,j);
      double f01r=original.get_real(i,j+1);
      double f11r=original.get_real(i+1,j+1);

      double f00i=original.get_imag(i,j);
      double f10i=original.get_imag(i+1,j);
      double f01i=original.get_imag(i,j+1);
      double f11i=original.get_imag(i+1,j+1);

      double new_i = (i+0.5)*scale_x ;
      double new_j = (j+0.5)*scale_y ;

      for(int di=0; di < scale_x; di++){
	for(int dj=0; dj < scale_y; dj++){

	  //double x = (di + 0.5*((scale_x+1) % 2)) /((double) scale_x);
	  //double y = (dj + 0.5*((scale_y+1) % 2)) /((double) scale_y);

	  double x = x_pos[di][dj];
	  double y = y_pos[di][dj];

	  double x_1 = 1-x;
	  double y_1 = 1-y;

	  double value_r = f00r*(x_1)*(y_1) + f10r*x*(y_1) + f01r*(x_1)*y + f11r*x*y;
	  big.set_real(new_i + di , new_j + dj , value_r);

	  double value_i = f00i*(x_1)*(y_1) + f10i*x*(y_1) + f01i*(x_1)*y + f11i*x*y;
	  big.set_imag(new_i + di, new_j + dj, value_i);

	}
      }
	  /**	  double x = 0.5 + (di+0.5)/((double) scale_x);
	  double y = 0.5 + (dj+0.5)/((double) scale_y);

	  int lower_i;// = i-1;
	  int lower_j;// = j-1;
	  
	  if(x>1){
	    x -= 1;
	    lower_i = i;
	  }
	  else{
	    lower_i = i-1;
	  }

	  if(y>1){
	    y -= 1;
	    lower_j = j;
	  }
	  else{
	    lower_j = j-1;
	  }

	  //	  else{
	  double f00=original.get(lower_i,lower_j);
	  double f10=original.get(lower_i+1,lower_j);
	  double f01=original.get(lower_i,lower_j+1);
	  double f11=original.get(lower_i+1,lower_j+1);
	  
	  //do the actual interpolation.
	  double value = f00*(1-x)*(1-y) + f10*x*(1-y) + f01*(1-x)*y + f11*x*y;
	  big.set(new_i+di,new_j+dj,value);
	  sum+=value; **/
	  
	    //	  }
      //	}
      // }
      
      /**      if(sum!=0){
	double norm = (original.get(i,j)*scale_x*scale_y)/sum;
	
	//loop again to normalise
	for(int di=0; di < scale_x; di++){
	  for(int dj=0; dj < scale_y; dj++){
	    double new_value = big.get(new_i+di,new_j+dj)*norm;
	    big.set(new_i+di,new_j+dj,new_value);
	  }
	}
	} **/
    }
  }

  for(int di=0; di < scale_x; di++){
    delete [] x_pos[di];
    delete [] y_pos[di];
  }
  delete [] x_pos;
  delete [] y_pos;
  
}


void shrink( Complex_2D & original, Complex_2D & small){

  int nx = original.get_size_x();
  int ny = original.get_size_y();
  int snx = small.get_size_x();
  int sny = small.get_size_y();

  //small.copy(original);
  if( (nx % snx) !=0 || (ny % sny) !=0 ){

    cerr<< "In shrink(Complex_2D & original, Complex_2D & small) " 
	<< "the size of 'original' must be an integer "
	<< "multiple of the size of 'small'" << endl;
    exit(0);
  }

  int scale_x = nx/snx;
  int scale_y = ny/sny;

  for(int i=0; i < snx; i++){
    for(int j=0; j < sny; j++){

      double real_value=0;
      double imag_value=0;

      for(int di=0; di < scale_x; di++){
	for(int dj=0; dj < scale_y; dj++){
	  int new_i = i*scale_x+di;
	  int new_j = j*scale_y+dj;
	  real_value+=original.get_real(new_i,new_j);
	  imag_value+=original.get_imag(new_i,new_j);
	}
      }
      
      double scale = scale_x*scale_y;
      small.set_real(i,j,real_value/scale);
      small.set_imag(i,j,imag_value/scale);
    }
  }

}

void shrink( Double_2D & original, Double_2D & small){

  int nx = original.get_size_x();
  int ny = original.get_size_y();
  int snx = small.get_size_x();
  int sny = small.get_size_y();

  //small.copy(original);
  if( (nx % snx) !=0 || (ny % sny) !=0 ){

    cerr<< "In shrink(Complex_2D & original, Complex_2D & small) " 
	<< "the size of 'original' must be an integer "
	<< "multiple of the size of 'small'" << endl;
    exit(0);
  }

  int scale_x = nx/snx;
  int scale_y = ny/sny;

  for(int i=0; i < snx; i++){
    for(int j=0; j < sny; j++){

      double value=0;

      for(int di=0; di < scale_x; di++){
	for(int dj=0; dj < scale_y; dj++){
	  int new_i = i*scale_x+di;
	  int new_j = j*scale_y+dj;
	  value+=original.get(new_i,new_j);
	}
      }
      
      double scale = scale_x*scale_y;
      small.set(i,j,value/scale);
    }
  }

}


////////////////////////////////

