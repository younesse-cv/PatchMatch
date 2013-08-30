#ifndef __quality_mesure_H__
#define __quality_mesure_H__
#include "defineall.h"

double PSNR(IplImage *original,IplImage *distorted);
double SSIM(IplImage *original,IplImage *distorted);
#endif
