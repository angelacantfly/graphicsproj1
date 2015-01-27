#include "ip.h"
#include "main.h"
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define RED	0
#define GREEN	1
#define BLUE	2

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
    // get width and height
    int width = src->getWidth();
    int height = src->getHeight();
    Image* blackImage = new Image(width,height);
    
    if (alpha > 0 && alpha <= 1)
        return ip_interpolate(src, blackImage, alpha);
    
    if (alpha > 1)
        return ip_interpolate(src, blackImage, alpha);

    return NULL;
    
}


/*
* shift colors
*/
Image* ip_color_shift(Image* src)
{
    // get width and height
    int width = src->getWidth();
    int height = src->getHeight();
    
    Image* newImage =  new Image(width, height);
    
    for (int w = 0 ; w < width; ++w) {
        for (int h = 0; h < height; ++h) {
            newImage->setPixel(w, h, GREEN, src->getPixel(w, h, RED));
            newImage->setPixel(w, h, BLUE, src->getPixel(w, h, GREEN));
            newImage->setPixel(w, h, RED, src->getPixel(w, h, BLUE));
        }
    }
    
	cerr << "Done!" << endl;
	return newImage;
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
    int width = src->getWidth();
    int height = src->getHeight();
    
    Image* grey = new Image(width, height);
    

    for (int w = 0; w <width; ++w)
        for (int h = 0; h<height; ++h )
        {
            grey->setPixel(w, h, RED, 0.5);
            grey->setPixel(w, h, GREEN, 0.5);
            grey->setPixel(w, h, BLUE, 0.5);
        }

    
    return ip_interpolate(src, grey, alpha);

}

double correctChannel(double value)
{
    if (value > 1) return 1;
    if (value < 0) return 0;
    return value;
}

/*
* convolve an image with a kernel
*/
Image* ip_convolve (Image* src, int size, double* kernel )
{
    int width = src->getWidth();
    int height = src->getHeight();
    assert(size%2 ==1); // makes sure that the size is odd number
    
    double currentRed;
    double currentGreen;
    double currentBlue;
    int countNeightbor = 0;
    
    int currentX;
    int currentY;
    
    Image* newImage =  new Image(width, height);
    
    for (int w = 0 ; w < width; ++w) {
        for (int h = 0; h < height; ++h) {
            currentRed = 0;
            currentGreen = 0;
            currentBlue = 0;
            countNeightbor = 0;
            
            for (int nx = -(size-1)/2; nx < (size-1)/2 + 1; ++nx )
                for (int ny = -(size-1)/2; ny < (size-1)/2 + 1; ++ny)
                {
                    currentX = w + nx;
                    currentY = h + ny;
                    if (currentX < 0 || currentX >= width || currentY < 0 || currentY >= height) {
                        continue;
                    }
                    currentRed += src->getPixel(currentX, currentY, RED) * kernel[countNeightbor];
                    currentGreen += src->getPixel(currentX, currentY, GREEN) * kernel[countNeightbor];
                    currentBlue += src->getPixel(currentX, currentY, BLUE) * kernel[countNeightbor];
                    countNeightbor++;
                }
            
            newImage->setPixel(w, h, RED, correctChannel(currentRed) );
            newImage->setPixel(w, h, GREEN, correctChannel(currentGreen));
            newImage->setPixel(w, h, BLUE, correctChannel(currentBlue));
            
        }
    }
    
    cerr << "Done!" << endl;
    return newImage;
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
    
    double* kernel = new double[9];
    
    for(int i = 0; i < 9; i++){
        kernel[i] = -1;
    }
    
    kernel[4] = 8;
    
    return ip_convolve(src, 3, kernel);
}


/*
* extract channel of input image
*/
Image* ip_extract (Image* src, int channel)
{
    // get width and height
    int width = src->getWidth();
    int height = src->getHeight();
    
    Image* newImage =  new Image(width, height);
    
    for (int w = 0 ; w < width; ++w) {
        for (int h = 0; h < height; ++h) {
            for (int c = 0; c < 3; ++c)
                if (c == channel) newImage->setPixel(w, h, c, src->getPixel(w, h, c));
        }
    }
    
    cerr << "Done!" << endl;
    return newImage;
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
    // get width and height
    int width = src->getWidth();
    int height = src->getHeight();
    double grey;
    double RED_CONSTANT = .2126;
    double GREEN_CONSTANT = .7152;
    double BLUE_CONSTANT = .0722;
    
    
    Image* newImage =  new Image(width, height);
    
    for (int w = 0 ; w < width; ++w) {
        for (int h = 0; h < height; ++h) {
            grey = RED_CONSTANT * src->getPixel(w, h, RED) +
                    GREEN_CONSTANT * src->getPixel(w, h, GREEN) +
                    BLUE_CONSTANT * src->getPixel(w, h, BLUE);
            newImage->setPixel(w, h, GREEN, grey);
            newImage->setPixel(w, h, BLUE, grey);
            newImage->setPixel(w, h, RED, grey);
        }
    }
    
    cerr << "Done!" << endl;
    return newImage;
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
    //FIXME: negative alpha value freaks out and doesn't end
    // Only seems to be a problem with brighten function
    // works for saturate and contrast
    
    // get width and height
    int width = src1->getWidth();
    int height = src1->getHeight();
    Pixel src1Pixel;
    Pixel src2Pixel;
    double currentChannel;
    
    
    Image* newImage =  new Image(width, height);
    
    for (int w = 0 ; w < width; ++w) {
        for (int h = 0; h < height; ++h) {
            src1Pixel = src1->getPixel(w, h);
            src2Pixel = src2->getPixel(w, h);
            for (int c = 0; c < 3; ++c)
            {
                currentChannel = alpha * src1Pixel.getColor(c) + (1.0-alpha)* src2Pixel.getColor(c);
                newImage->setPixel(w, h, c, correctChannel(currentChannel));

            }
        }
    }
    
    cerr << "Done!" << endl;
    return newImage;
}
/*
* invert input image
*/
Image* ip_invert (Image* src)
{
    int width = src->getWidth();
    int height = src->getHeight();
    

    Image* grey = new Image(width,height);
    
    for (int w =0; w<width; ++w)
        for (int h=0; h<height; ++h)
        {
            grey->setPixel(w, h, RED, 0.5);
            grey->setPixel(w, h, GREEN, 0.5);
            grey->setPixel(w, h, BLUE, 0.5);
        }
    
    return ip_interpolate(src, grey, -1);
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
    Image* greyImage = ip_grey(src);

    return ip_interpolate(src, greyImage, alpha);

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
    // get width and height
    int width = src->getWidth();
    int height = src->getHeight();
    
    
    Image* newImage =  new Image(width, height);
    
    for (int w = 0 ; w < width; ++w) {
        for (int h = 0; h < height; ++h) {
            for (int c = 0; c < 3; ++ c)
                if (src->getPixel(w, h, c) > cutoff) {
                    newImage->setPixel(w, h, c, 1);
                }
        }
    }
    
    cerr << "Done!" << endl;
    return newImage;
}




