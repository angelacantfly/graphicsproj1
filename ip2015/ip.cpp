#include "ip.h"
#include "main.h"
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <time.h>



/*
* convolve with a box filter
*/
Image* ip_blur_box (Image* src, int size)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}


/*
* convolve with a gaussian filter
*/
Image* ip_blur_gaussian (Image* src, int size, double sigma)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}


/*
* convolve with a triangle filter
*/
Image* ip_blur_triangle (Image* src, int size)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}


/*
* interpolate with a black image
*/
Image* ip_brighten (Image* src, double alpha)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}


/*
* shift colors
*/
Image* ip_color_shift(Image* src)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}


/*
* use a mask image for a per-pixel alpha value to perform
* interpolation with a second image
*/
Image* ip_composite (Image* src1, Image* src2, 
					 Image* mask)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}


/*
* interpolate with the average intensity of the src image
*/
Image* ip_contrast (Image* src, double alpha)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}



/*
* convolve an image with a kernel
*/
Image* ip_convolve (Image* src, int size, double* kernel )
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}



/*
*  create cropped version of image
*/
Image* ip_crop (Image* src, int x0, int y0, int x1, int y1)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}
/*
* convolve with an edge detection kernel
*/
Image* ip_edge_detect (Image* src)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}


/*
* extract channel of input image
*/
Image* ip_extract (Image* src, int channel)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}


/*
* create your own fun warp
*/
Image* ip_fun_warp (Image* src)
{

	//  ask user for input parameters here including resampling method and, 
	//  if gaussian resampling is used, its filtersize and sigma
	//  if you implement more than one warp, you should ask the 
	//  user to chose the one to perform here too!

	return NULL;
}
/*
* create a new image with values equal to the psychosomatic intensities
* of the source image
*/
Image* ip_grey (Image* src)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}


/*
*  shift image by dx and dy (modulo width & height)
*/

Image* ip_image_shift (Image* src, double dx, double dy)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}
/*
* interpolate an image with another image
*/
Image* ip_interpolate (Image* src1, Image* src2, double alpha)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}
/*
* invert input image
*/
Image* ip_invert (Image* src)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}


/*
* define your own filter
* you need to request any input paraters here, not in control.cpp
*/

Image* ip_misc(Image* src)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}




/*
* round each pixel to the nearest value in the new number of bits
*/
Image* ip_quantize_simple (Image* src, int bitsPerChannel)
{

	cerr << "This function is not implemented." << endl;
	return NULL;
}


/*
* dither each pixel to the nearest value in the new number of bits
* using a static 4x4 matrix
*/
Image* ip_quantize_ordered (Image* src, int bitsPerChannel)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}


/*
* dither each pixel to the nearest value in the new number of bits
* using error diffusion
*/
Image* ip_quantize_fs (Image* src, int bitsPerChannel)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}

/*
* nearest neighbor sample
*/
Pixel ip_resample_nearest(Image* src, double x, double y) {
	cerr << "This function is not implemented." << endl;
	Pixel myPixel(0,0,0);

	return myPixel;
}

/*
* bilinear sample
*/

Pixel ip_resample_bilinear(Image* src, double x, double y) {
	cerr << "This function is not implemented." << endl;
	Pixel myPixel(0,0,0);
	return myPixel;
}

/*
* gausian sample
*/
Pixel ip_resample_gaussian(Image* src, double x, double y, int size, double sigma)
{
	cerr << "This function is not implemented." << endl;
	Pixel myPixel(0,0,0);
	return myPixel;
}

/*
* rotate image using one of three sampling techniques
*/
Image* ip_rotate (Image* src, double theta, int x, int y, int samplingMode, 
				  int gaussianFilterSize, double gaussianSigma)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}


/*
* change saturation
*/
Image* ip_saturate (Image* src, double alpha)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}


/*
* scale image using one of three sampling techniques
*/
Image* ip_scale (Image* src, double xFac, double yFac, int samplingMode, 
				 int gaussianFilterSize, double gaussianSigma)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}


/*
* threshold image
*/
Image* ip_threshold (Image* src, double cutoff)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}




