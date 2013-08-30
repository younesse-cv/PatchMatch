#ifndef __MASKED_IMAGE_H_
#define __MASKED_IMAGE_H_
#include "defineall.h"

void initSimilarity();
MaskedImage_P initMaskedImage(IplImage* image, int** mask);
MaskedImage_P initNewMaskedImage(int width, int height);
void freeMaskedImage(MaskedImage_P mIm);
int isMasked(MaskedImage_P mIm, int x, int y);
void setMask(MaskedImage_P mIm, int x, int y, int value);
int constainsMasked(MaskedImage_P mIm, int x, int y, int S);
int distanceMaskedImage(MaskedImage_P source,int xs,int ys, MaskedImage_P target,int xt,int yt, int S);
MaskedImage_P downsample(MaskedImage_P source);
MaskedImage_P upscale(MaskedImage_P source, int newW,int newH);
MaskedImage_P copyMaskedImage(MaskedImage_P mIm);

int getSampleMaskedImage(MaskedImage_P mIm, int x, int y, int band);////////////
void setSampleMaskedImage(MaskedImage_P mIm, int x, int y, int band, int value);/////////////


// Variables globales
extern double* G_globalSimilarity;
extern int G_initSim;

#endif
