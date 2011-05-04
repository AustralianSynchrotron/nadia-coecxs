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

      for(double x = (int) i_; x <= i_ + scale; x++ ){
	for(double y = (int) j_; y <= j_ + scale; y++){

	  //check that we are still in range
	  if(i_<nx&&i_>=0&&j_<ny&&j_>=0){

	    //subtract the remainders from one.
	    double x_fraction = 1;
	    double y_fraction = 1;

	    if(x<i_)
	      x_fraction = 1 - fmod(i_, 1.0);
	    if(y<j_)
	      y_fraction = 1 - fmod(j_, 1.0);

	    if(x > i_ + scale)
	      x_fraction = fmod(i_, 1.0);

	    if(y > j_ + scale)
	      y_fraction = fmod(j_, 1.0);

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
  
  double total = 0;

  //loop over the array once to sort into bins
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      total+=pow(image.get(i,j),4);
    }
  }
  
  return total;
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

  cout << endl << "a = [";
  for(int c=0; c < BINS; c++){

    cout << counts[c];
    if(c!=BINS-1)
      cout << ",";
    else
      cout << "]"<<endl;
      
    if(counts[c]!=0){
      //      cout << " log2(counts[c]) = " <<log2(counts[c])<<endl;
      entropy-= counts[c]*log2(counts[c]);
    }
  }
  
  //cout << "here3 "<<entropy << endl;
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

  image.scale(image.get_sum());

  convolve(image,3.0,4);

  char buf[50];
  static int counter = 0;
  sprintf(buf,"image_conv_%i.tiff",counter);
  write_image(buf,image);

  double total = 0;
  double value_x ;
  double value_y;

  int nx = image.get_size_x();
  int ny = image.get_size_y();

  // double max = image.get_max()-image.get_min();
  //  threshold = pow(threshold*max,2);

  //  double min_allowed = threshold*4*max*max;

  Double_2D output(nx,ny);

  double max = 0;
  double smallest_max = 0;
  double no_of_max = 2;
  double max_total = 0;

  //copy image to the in array
  for(int i=2; i < nx-2 ; i++){
    for(int j=2; j < ny-2 ; j++){

      //do the x direction
      value_x  = -1*image.get(i-1,j-1);
      value_x -=  2*image.get(i-1,j);
      value_x -=  1*image.get(i-1,j+1);

      //    value_x -=  1*image.get(i-2,j);

      value_x +=  1*image.get(i+1,j-1);
      value_x +=  2*image.get(i+1,j);
      value_x +=  1*image.get(i+1,j+1);

      //      value_x -=  1*image.get(i+2,j);

      //do the y direction
      value_y  = -1*image.get(i-1,j-1);
      value_y -=  2*image.get(i,j-1);
      value_y -=  1*image.get(i+1,j-1);

      //      value_y -=  1*image.get(i,j-2);

      value_y +=  1*image.get(i-1,j+1);
      value_y +=  2*image.get(i,j+1);
      value_y +=  1*image.get(i+1,j+1);

      // value_x = image.get(i-1,j-1) - image.get(i,j); 
      //value_y = image.get(i-1,j) - image.get(i,j-1); 

      //      value_y +=  1*image.get(i,j+2);

      double value_total = sqrt(value_x*value_x + value_y*value_y);

      output.set(i,j,value_total);
      total += value_total;

      if(max<value_total){
	smallest_max = max;
	max = value_total;
      }
      else if(smallest_max<value_total){
	smallest_max = value_total;
      }
    }
  }

  max_total = (max + smallest_max)/2.0;
  
  double threshold = 0.5;
  double mean = total/(double)(nx*ny);
  cout << "mean is: " << mean << endl;
  cout << "max_combined is: " << max_total << endl;
  //double max = output.get_max()<<endl;

  //now threshold 
  /**  double total_t=0.0;
  double number=0.0;
  for(int i=0; i < nx ; i++){
    for(int j=0; j < ny ; j++){
      if(output.get(i,j) < mean){ //max*threshold)
	output.set(i,j,0);
      }
      else{
	total_t+=(output.get(i,j)-mean)/(max-mean);
	number++;
      }
    }
  }

  cout << "count is: " << number << endl;


  cout << "max is: " << max << endl;
  cout << "normalised max: " << max_total/image.get_sum() << endl; **/

  sprintf(buf,"sobel_%i.tiff",counter);
  write_image(buf,output);
  counter++;

  /**
  cout << "Entropy of gradient is :"<< calculate_image_entropy(output)<<endl;**/

  return output.get_max();

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

  double min_allowed = threshold*4*max*max;

  max = 0;

  Double_2D output(nx,ny);
  //copy image to the in array
  for(int i=1; i < nx-1 ; i++){
    for(int j=1; j < ny-1 ; j++){
      value = pow(fabs(image.get(i,j)-image.get(i-1,j)),2);      
      value+=pow(fabs(image.get(i,j)-image.get(i,j-1)),2);
      value+= pow(fabs(image.get(i,j)-image.get(i+1,j)),2);      
      value+=pow(fabs(image.get(i,j)-image.get(i,j+1)),2);

      output.set(i,j,value);
      total+=value;

      if(value> max)
	max = value;
    }
  }

  Double_2D output_2(nx,ny);

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

  return total_t/total;

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

  double sum = 0;
  double total = 0;

  for(int i=1; i < nx ; i++){
    for(int j=0; j < ny ; j++){
      sum += image.get(i,j)*image.get(i-1,j);
      total += image.get(i,j)*image.get(i,j);
    }
  }
  return sum-(total*total);
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


