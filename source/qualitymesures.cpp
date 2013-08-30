#include "defineall.h"

double PSNR(IplImage *Original , IplImage *inpainted)
{
    //First let's convert the images to Y component only
    double sumR = 0., sumG = 0., sumB = 0., sumRGB;
    CvSize size;
    double psnr, EQM = 0.;
    size = cvGetSize(Original);
    double h=size.height;
    double w=size.width;

    for (int i=0 ; i<h ; i++ )
        for ( int j=0 ; j<w; j++ ) {
            sumR+= pow( (cvGet2D(Original,i,j).val[0] - cvGet2D(inpainted,i,j).val[0]),2 );
            sumG+= pow( (cvGet2D(Original,i,j).val[1] - cvGet2D(inpainted,i,j).val[1]),2 );
            sumB+= pow( (cvGet2D(Original,i,j).val[2] - cvGet2D(inpainted,i,j).val[2]),2 );
        }
    sumRGB=sumR+sumG+sumB;
    EQM= sumRGB /(3*h*w);
    psnr=10.0*log10( 65025/EQM );
    return psnr;
}

/*This function uses the algorithme of Rabah Mehdi
    *The equivalent of Zhou Wang's SSIM matlab code using OpenCV.
    * from http://www.cns.nyu.edu/~zwang/files/research/ssim/index.html
    * The measure is described in :
    * "Image quality assessment: From error measurement to structural similarity"*/

double SSIM(IplImage *original,IplImage *distorted)
{
    double C1 = 6.100210, C2 = 108.102210;
    double SSIM;

    IplImage
            *img1=NULL, *img2=NULL, *img1_img2=NULL,
            *img1_temp=NULL, *img2_temp=NULL,
            *img1_sq=NULL, *img2_sq=NULL,
            *mu1=NULL, *mu2=NULL,
            *mu1_sq=NULL, *mu2_sq=NULL, *mu1_mu2=NULL,
            *sigma1_sq=NULL, *sigma2_sq=NULL, *sigma12=NULL,
            *ssim_map=NULL, *temp1=NULL, *temp2=NULL, *temp3=NULL;

    /***************************** INITS **********************************/
    if (original==NULL || distorted==NULL)
        return -1;

    int x=original->width, y=original->height;
    int nChan=original->nChannels, d=IPL_DEPTH_32F;
    CvSize size = cvSize(x, y);

    img1 = cvCreateImage( size, d, nChan);
    img2 = cvCreateImage( size, d, nChan);

    cvConvert(original, img1);
    cvConvert(distorted, img2);
    cvReleaseImage(&original);
    cvReleaseImage(&distorted);


    img1_sq = cvCreateImage( size, d, nChan);
    img2_sq = cvCreateImage( size, d, nChan);
    img1_img2 = cvCreateImage( size, d, nChan);

    cvPow( img1, img1_sq, 2 );
    cvPow( img2, img2_sq, 2 );
    cvMul( img1, img2, img1_img2, 1 );

    mu1 = cvCreateImage( size, d, nChan);
    mu2 = cvCreateImage( size, d, nChan);

    mu1_sq = cvCreateImage( size, d, nChan);
    mu2_sq = cvCreateImage( size, d, nChan);
    mu1_mu2 = cvCreateImage( size, d, nChan);


    sigma1_sq = cvCreateImage( size, d, nChan);
    sigma2_sq = cvCreateImage( size, d, nChan);
    sigma12 = cvCreateImage( size, d, nChan);

    temp1 = cvCreateImage( size, d, nChan);
    temp2 = cvCreateImage( size, d, nChan);
    temp3 = cvCreateImage( size, d, nChan);

    ssim_map = cvCreateImage( size, d, nChan);

    /*************************** END INITS **********************************/
    //////////////////////////////////////////////////////////////////////////
    // PRELIMINARY COMPUTING
    cvSmooth( img1, mu1, CV_GAUSSIAN, 11, 11, 1.10 );
    cvSmooth( img2, mu2, CV_GAUSSIAN, 11, 11, 1.10 );

    cvPow( mu1, mu1_sq, 2 );
    cvPow( mu2, mu2_sq, 2 );
    cvMul( mu1, mu2, mu1_mu2, 1 );

    cvSmooth( img1_sq, sigma1_sq, CV_GAUSSIAN, 11, 11, 1.10 );
    cvAddWeighted( sigma1_sq, 1, mu1_sq, -1, 0, sigma1_sq );

    cvSmooth( img2_sq, sigma2_sq, CV_GAUSSIAN, 11, 11, 1.10 );
    cvAddWeighted( sigma2_sq, 1, mu2_sq, -1, 0, sigma2_sq );

    cvSmooth( img1_img2, sigma12, CV_GAUSSIAN, 11, 11, 1.10 );
    cvAddWeighted( sigma12, 1, mu1_mu2, -1, 0, sigma12 );

    //////////////////////////////////////////////////////////////////////////
    // FORMULA
    // (2*mu1_mu2 + C1)
    cvScale( mu1_mu2, temp1, 2 );
    cvAddS( temp1, cvScalarAll(C1), temp1 );

    // (2*sigma12 + C2)
    cvScale( sigma12, temp2, 2 );
    cvAddS( temp2, cvScalarAll(C2), temp2 );

    // ((2*mu1_mu2 + C1).*(2*sigma12 + C2))
    cvMul( temp1, temp2, temp3, 1 );

    // (mu1_sq + mu2_sq + C1)
    cvAdd( mu1_sq, mu2_sq, temp1 );
    cvAddS( temp1, cvScalarAll(C1), temp1 );

    // (sigma1_sq + sigma2_sq + C2)
    cvAdd( sigma1_sq, sigma2_sq, temp2 );
    cvAddS( temp2, cvScalarAll(C2), temp2 );

    // ((mu1_sq + mu2_sq + C1).*(sigma1_sq + sigma2_sq + C2))
    cvMul( temp1, temp2, temp1, 1 );

    // ((2*mu1_mu2 + C1).*(2*sigma12 + C2))./((mu1_sq + mu2_sq + C1).*(sigma1_sq + sigma2_sq + C2))
    cvDiv( temp3, temp1, ssim_map, 1 );
    CvScalar index_scalar = cvAvg( ssim_map );
    // through observation, there is approximately
    //  error max with the original matlab program
    SSIM=(index_scalar.val[2]+index_scalar.val[1]+index_scalar.val[0])/3;

    return SSIM;
}
