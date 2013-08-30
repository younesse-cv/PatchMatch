#ifndef __INPAINT_H_
#define __INPAINT_H_
#include "defineall.h"

Inpaint_P initInpaint();

void addEltInpaintingPyramid(Inpaint_P imp, MaskedImage_P elt);
MaskedImage_P ExpectationMaximization(Inpaint_P imp, int level);
void ExpectationStep(NNF_P nnf, int sourceToTarget, double*** vote, MaskedImage_P source, int upscale);
void weightedCopy(MaskedImage_P src, int xs, int ys, double*** vote, int xd, int yd, double w);
void MaximizationStep(MaskedImage_P target, double*** vote);
IplImage* inpaint(Inpaint_P imp, IplImage* input, int ** mask, int radius);
void freeInpaintingPyramid(Inpaint_P imp);

// Variables globales
extern double* G_globalSimilarity;
extern int G_initSim;
#endif
