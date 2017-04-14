#include "defineall.h"

// Disable to measure the PSNR or the SSIM measure 

#define MESURES 0

double* G_globalSimilarity = NULL;
int G_initSim = 0;

double max1(double a, double b)
{
	return (a + b + fabs(a-b) ) / 2;
}

double min1(double a, double b)
{
	return (a + b - fabs(a-b) ) / 2;
}

int main(int argc, char** argv) 

{

	printf("------------------------------------------------------------------ \n");
	printf("                  Inpainting based PatchMatch                      \n");
	printf(" ------------------------------------------------------------------\n");
	printf("      PatchMatch : A Randomized Correspondence Algorithm           \n");
	printf("               for Structural Image Editing                        \n");
	printf("   by C.Barnes,E.Shechtman,A.Finkelstein and Dan B.Goldman         \n");
	printf("  ACM Transactions on Graphics (Proc. SIGGRAPH), vol.28, aug-2009  \n");
	printf("              For more information please refer to:                \n");
	printf(" http://www.cs.princeton.edu/gfx/pubs/Barnes_2009_PAR/index.php    \n");
	printf("\n");
	printf("                   @Author: Younesse ANDAM                         \n"); 
	printf("            Contact: younesse.andam@gmail.com                      \n");
	printf("                Copyright (c) 2010-2011                            \n");


	clock_t tic,toc;
	float cpu_time;     /* Total CPU time in minutes */ 


	uchar* data;
	CvSize size;
	double height,width;

	//Put your file name here
	char fileNameInput[150] = "C:\\Users\\PVHP9868\\Desktop\\PATCH_MATCH_APPLICATION_CONSOLE\\image_files\\forest\\forest.bmp";
	char fileNameMasked[150] = "C:\\Users\\PVHP9868\\Desktop\\PATCH_MATCH_APPLICATION_CONSOLE\\image_files\\forest\\forest_pruned.bmp";

	char fileNameOutput[150];
	char fileNameOutputMasked[150];
	strcpy(fileNameOutput, fileNameInput);
	strcpy(fileNameOutputMasked, fileNameMasked);
	strcat(fileNameOutput, "_RESULT.png");
	strcat(fileNameOutputMasked, "res.png");

	IplImage* input_img=NULL , *maskimage=NULL,*distorted=NULL , *output_img=NULL,
		*inpaint_mask_gray=NULL,*input_gray=NULL, *distorted_gray=NULL;

	input_img=cvLoadImage(fileNameInput, 1);

    if (!input_img) {
		printf("Could not load image file: %s\n",fileNameInput); 
		exit(1);
	}

	size = cvGetSize(input_img);
	height=size.height;
	width=size.width;

	distorted=cvLoadImage(fileNameMasked, 1);

    if (!distorted) {
        printf("Could not load image file: %s\n",fileNameMasked);
        exit(1);
	}

	//Mask computation
	maskimage=cvCreateImage(size,8,1);
	cvZero(maskimage);

	input_gray=cvCreateImage(size,8,1);
	cvZero(input_gray);

	distorted_gray=cvCreateImage(size,8,1);
	cvZero(distorted_gray);

	cvCvtColor(distorted,maskimage,CV_BGR2GRAY);

    for ( int i=0 ; i < height ; ++i )
		for ( int j=0 ; j<width;  ++j )
            if (cvGet2D(maskimage,i,j).val[0] !=255 )
				cvSet2D(maskimage,i,j,cvScalar(0,0,0));

	//display images
	cvNamedWindow("Original image", CV_WINDOW_AUTOSIZE);
	cvShowImage("Original image",input_img);

	cvNamedWindow("MASK", CV_WINDOW_AUTOSIZE);
	cvShowImage("MASK",distorted);

	// generate mask array from mask image
	int channels=maskimage->nChannels;
	int step = maskimage->widthStep/sizeof(uchar);
	int ** mask;
	mask = (int **) calloc(int(height),sizeof(int*));
    for ( int i=0 ; i<height ; i++)
		mask[i] = (int *) calloc(int(width),sizeof(int));

	printf("----------------------------------------------------------------------\n");
	printf("\n");
	printf("Computing, please wait, this operation may take several minutes...\n");

	data = (uchar *)maskimage->imageData;
	//Timer: tic, toc
	tic = clock ();
    for ( int i = 0 ; i < height ; ++i )
        for ( int j = 0 ; j < width ; ++j )
            if ( data[i*step+j*channels]==255 )
				mask[i][j]=1;	   


		Inpaint_P inp = initInpaint();
		output_img = inpaint(inp, input_img, (int**)mask, 2);
        if (!cvSaveImage(fileNameOutput,output_img))

			printf("/!\\/!\\/!\\/!\\/!\\/!\\Could not save the resultant image. Check the path of saving.../!\\/!\\/!\\/!\\/!\\/!\\\n");

		//cvSaveImage("mask.bmp",maskimage);
		
/*---------------------------------------------- Quality Mesure ----------------------------------------------*/ 
// It will be useful so far to measure the quality of the recontructed areas against the original ones 

#if MESURES
		double ssim;
		double psnr;

		psnr=PSNR(input_img,output_img );
		ssim=SSIM(input_img,output_img);

		printf("Quality mesures: \n");

		printf("--> PSNR=%f dB\n",psnr);
		printf("\n");
		printf("--> SSIM= %f \n",ssim);
#endif
		/*------------------------------------------------------------------------------------------------------------*/
		printf("\n");
		printf("DONE\n");
		printf("The result is saved in: %s\n",fileNameOutput);
		toc = clock ();
		cpu_time = float((toc - tic)/CLOCKS_PER_SEC) ;

		int secondes;
		int time=(int)(cpu_time/60);
		secondes=int(cpu_time-time*60);

		printf("The CPU time is %i minutes and %i secondes ",time,secondes);

		cvWaitKey();

		//free memory

		cvDestroyAllWindows();

		for( int i = 0 ; i < height ; ++i)
			free(mask[i]);
		free(mask);

		cvReleaseImage(&input_img);
		cvReleaseImage(&maskimage);
		cvReleaseImage(&output_img);
		cvReleaseImage(&distorted);
		cvReleaseImage(&inpaint_mask_gray);
		cvReleaseImage(&input_gray);
		cvReleaseImage(&distorted_gray);

		return 0;
}
