#ifndef __NNF_H_
#define __NNF_H_
#include "defineall.h"

NNF_P initNNF(MaskedImage_P input, MaskedImage_P output, int patchsize);
void allocNNFField(NNF_P nnf);
void freeNNFField(NNF_P nnf);
void freeNNF(NNF_P nnf);
void randomize(NNF_P nnf);
void initializeNNFFromOtherNNF(NNF_P nnf, NNF_P otherNnf);
void initializeNNF(NNF_P nnf);
void minimizeNNF(NNF_P nnf, int pass);
void minimizeLinkNNF(NNF_P nnf, int x, int y, int dir);
int distanceNNF(NNF_P nnf, int x, int y, int xp,int yp);

#endif
