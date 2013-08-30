#ifndef __quality_mesure_H_
#define __quality_mesure_H_
#include "defineall.h"

double PSNR(IplImage *original,IplImage *distorted);
double SSIM(IplImage *original,IplImage *distorted);
#endif
