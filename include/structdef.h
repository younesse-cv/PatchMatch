#ifndef __STRUCT_PATCHMATCH_H_
#define __STRUCT_PATCHMATCH_H_

// the maximum value returned by MaskedImage.distance() 
#define DSCALE 65535
#define INCREASE_PYRAMID_SIZE_RATE 2


typedef struct{
	// image data
	int ** mask;
	IplImage* image;

	// array for converting distance to similarity
    double * similarity;

	int isNew;
} MaskedImage_T;

typedef MaskedImage_T* MaskedImage_P;

typedef struct{
	// image 
	MaskedImage_P input;
	MaskedImage_P output;

	//  patch size
	int S;

	// Nearest-Neighbor Field 1 pixel = { x_target, y_target, distance_scaled } 
	int *** field;
	int fieldH;
	int fieldW;
} NNF_T;
typedef NNF_T* NNF_P;

typedef struct{
	//initial image
	MaskedImage_P initial;
	
	// Nearest-Neighbor Fields
	NNF_P nnf_SourceToTarget;
	NNF_P nnf_TargetToSource;
	
	// patch radius
	int radius;
	
	// Pyramid of downsampled initial images
	MaskedImage_P* pyramid;
	int nbEltPyramid;
	int nbEltMaxPyramid;
} Inpaint_T;

typedef Inpaint_T* Inpaint_P;

#endif
