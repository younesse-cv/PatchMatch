#include "maskedimage.h"

void initSimilarity()
{
	int i, j, k, length;
	double base[11] = {1.0, 0.99, 0.96, 0.83, 0.38, 0.11, 0.02, 0.005, 0.0006, 0.0001, 0 };
	double t, vj, vk;
	length = (DSCALE+1);
    if (!G_initSim) {
		G_globalSimilarity = (double *) calloc(length, sizeof(double));
        for ( i=0 ; i<length ; ++i) {
			t = (double)i/length;
			j = (int)(100*t);
			k=j+1;
			vj = (j<11)?base[j]:0;
			vk = (k<11)?base[k]:0;
			G_globalSimilarity[i] = vj + (100*t-j)*(vk-vj);
		}
	}
	G_initSim = 1;
}

void initSimilarity2()
{
	int i, length;
	double t_halfmax=0.04, t, coef;
	length = (DSCALE+1);
	if (!G_initSim){
		G_globalSimilarity = (double *) calloc(length, sizeof(double));
		coef = -log(0.5)/pow(t_halfmax,2);
        for (i=0;i<length;i++) {
			t = (double)i/length;
			G_globalSimilarity[i] = exp(-(t*t)*coef);
		}
	}
	G_initSim = 1;
}

MaskedImage_P initMaskedImage(IplImage* image, int** mask)
{
	MaskedImage_P mIm = (MaskedImage_P)malloc(sizeof(MaskedImage_T));
	// image data
	mIm->mask = mask;
	mIm->image = image;

	initSimilarity();
	mIm->similarity = G_globalSimilarity;

	mIm->isNew=0;

	return mIm;
}

MaskedImage_P initMaskedImageFromImage(IplImage* image, int** mask)
{
	MaskedImage_P mIm = (MaskedImage_P)malloc(sizeof(MaskedImage_T));

	// image data
	mIm->mask = mask;
	mIm->image = image;

	initSimilarity();
	mIm->similarity = G_globalSimilarity;

	mIm->isNew=0;

	return mIm;
}

MaskedImage_P initNewMaskedImage(int width, int height)
{
	int j;
	MaskedImage_P mIm = (MaskedImage_P)malloc(sizeof(MaskedImage_T));

	// image data
	mIm->mask = (int**)malloc(sizeof(int*)*height);
    for (j=0;j<height;j++)
		mIm->mask[j] = (int*)calloc(width, sizeof(int));

	mIm->image = cvCreateImage(cvSize(width,height),IPL_DEPTH_8U,3);

	initSimilarity();
	mIm->similarity = G_globalSimilarity;

	mIm->isNew=1;

	return mIm;
}

void freeMaskedImage(MaskedImage_P mIm)
{
	int j;
    if (mIm!=NULL) {
        if (mIm->isNew) {
            if (mIm->mask!=NULL) {
                for (j=0;j<mIm->image->height;j++)
					free(mIm->mask[j]);

				free(mIm->mask);
				mIm->mask=NULL;
			}
            if (mIm->image!=NULL) {
				cvReleaseImage(&mIm->image);
				mIm->image=NULL;
			}
		}
		free(mIm);
		mIm = NULL;
	}
}


int getSampleMaskedImage(MaskedImage_P mIm, int x, int y, int band)
{
	uchar* data;
	int channels=mIm->image->nChannels;
	int step = mIm->image->widthStep/sizeof(uchar);
	data = (uchar *)mIm->image->imageData;
	return data[x*step+y*channels+band];
}

void setSampleMaskedImage(MaskedImage_P mIm, int x, int y, int band, int value)
{
	uchar* data;
	int channels=mIm->image->nChannels;
	int step = mIm->image->widthStep/sizeof(uchar);
	data = (uchar *)mIm->image->imageData;
	data[x*step+y*channels+band] = value;

}

int isMasked(MaskedImage_P mIm, int x, int y)
{
	if (mIm==NULL || mIm->mask==NULL)
		return 0; 
	return mIm->mask[x][y];
}

void setMask(MaskedImage_P mIm, int x, int y, int value) {
	if (mIm==NULL || mIm->mask==NULL)
		return;
	mIm->mask[x][y]=0<value;
}

// return true if the patch contains one (or more) masked pixel
int constainsMasked(MaskedImage_P mIm, int x, int y, int S)
{
	int dy, dx;
	int xs, ys;
    for (dy=-S;dy<=S;dy++) {
        for (dx=-S;dx<=S;dx++) {
			xs=x+dx; 
			ys=y+dy;
            if (xs<0 || xs>=mIm->image->height)
                continue;
            if (ys<0 || ys>=mIm->image->width)
                continue;
            if (mIm->mask[xs][ys])
                return 1;
		}
	}
	return 0;
}

// distance between two patches in two images
int distanceMaskedImage(MaskedImage_P source,int xs,int ys, MaskedImage_P target,int xt,int yt, int S)
{
	long double distance=0;
	long double wsum=0, ssdmax = 9*255*255;
	int dy, dx, band;
	int xks, yks;
	int xkt, ykt;
	long double ssd; 
	long res;
	int s_value, t_value, s_gx, t_gx, s_gy, t_gy;

	// for each pixel in the source patch
    for ( dy=-S ; dy<=S ; ++dy ) {
        for ( dx=-S ; dx<=S ; ++dx ) {

			xks = xs+dx;
			yks = ys+dy;
			xkt=xt+dx;
			ykt=yt+dy;
			wsum++;

			if ( xks<1 || xks>=source->image->height-1 ) {distance++; continue;}
			if ( yks<1 || yks>=source->image->width-1 ) {distance++; continue;}

			// cannot use masked pixels as a valid source of information
			if (isMasked(source, xks, yks)) {distance++; continue;}

			// corresponding pixel in the target patch
			if (xkt<1 || xkt>=target->image->height-1) {distance++; continue;}
			if (ykt<1 || ykt>=target->image->width-1) {distance++; continue;}

			// cannot use masked pixels as a valid source of information
			if (isMasked(target, xkt, ykt)) {distance++; continue;}

			ssd=0;
            for (band=0; band<3; ++band) {
				// pixel values
				s_value = getSampleMaskedImage(source, xks, yks, band);
				t_value = getSampleMaskedImage(source, xkt, ykt, band);

				// pixel horizontal gradients (Gx)
				s_gx = 128+(getSampleMaskedImage(source, xks+1, yks, band) - getSampleMaskedImage(source, xks-1, yks, band))/2;
				t_gx = 128+(getSampleMaskedImage(target, xkt+1, ykt, band) - getSampleMaskedImage(target, xkt-1, ykt, band))/2;

				// pixel vertical gradients (Gy)
				s_gy = 128+(getSampleMaskedImage(source, xks, yks+1, band) - getSampleMaskedImage(source, xks, yks-1, band))/2;
				t_gy = 128+(getSampleMaskedImage(target, xkt, ykt+1, band) - getSampleMaskedImage(target, xkt, ykt-1, band))/2;

				ssd += pow((long double)s_value-t_value , 2); // distance between values in [0,255^2]
				ssd += pow((long double)s_gx-t_gx , 2); // distance between Gx in [0,255^2]
				ssd += pow((long double)s_gy-t_gy , 2); // distance between Gy in [0,255^2]
			}

			// add pixel distance to global patch distance
			distance += ssd/ssdmax;
		}
	}

	res = (int)(DSCALE*distance/wsum);
	if (res < 0 || res > DSCALE) return DSCALE;
	return res;
}

// return a copy of the image
MaskedImage_P copyMaskedImage(MaskedImage_P mIm)
{
	int j, x;
	int W = mIm->image->width;
	int H = mIm->image->height;
	int ** newmask;
	MaskedImage_P copy;
	newmask = (int**)malloc(sizeof(int*)*H);
    for (j=0;j<H;j++) {
		newmask[j] = (int*)calloc(W, sizeof(int));
		for(x=0;x<W;x++)
			newmask[j][x] = mIm->mask[j][x];
	}

	IplImage* newimage = cvCloneImage(mIm->image);

	copy = initMaskedImage(newimage,newmask);
	copy->isNew=1;

	return copy;
}

// return a downsampled image (factor 1/2)
MaskedImage_P downsample2(MaskedImage_P source) {
	int kernel[6] = {1,2,4,4,2,1};
	int H, W;
	int x, y;
	int xs, ys;
	int dx, dy;
	int xk, yk, ky, k;
	int r=0,g=0,b=0,m=0,ksum=0;
	H=source->image->height;
	W=source->image->width;
	int newW=W/2, newH=H/2;

	MaskedImage_P newimage = initNewMaskedImage(newW, newH);
	xs=0;
	for(x=0;x<newH;++x) {
		ys=0;
		for(y=0;y<newW;++y) {
			r=0; g=0; b=0; m=0; ksum=0;

			for(dy=-2;dy<=3;dy++) {
				yk=ys+dy;
				if (yk<0 || yk>=W) continue;
				ky = kernel[2+dy];
				for(dx=-2;dx<=3;dx++) {
					xk = xs+dx;
					if (xk<0 || xk>=H) continue;

					if (source->mask[xk][yk]) continue;

					k = kernel[2+dx]*ky;
					r+= k*getSampleMaskedImage(source, xk, yk, 0);
					g+= k*getSampleMaskedImage(source, xk, yk, 1);
					b+= k*getSampleMaskedImage(source, xk, yk, 2);
					ksum+=k;
					m++;
				}
			}
			if (ksum>0) {r/=ksum; g/=ksum; b/=ksum;}

			if (m!=0) {
				setSampleMaskedImage(newimage, x, y, 0, r); 
				setSampleMaskedImage(newimage, x, y, 1, g); 
				setSampleMaskedImage(newimage, x, y, 2, b);
				setMask(newimage, x, y, 0);
			} else {
				setMask(newimage, x, y, 1);
				setSampleMaskedImage(newimage, x, y, 0, 0); 
				setSampleMaskedImage(newimage, x, y, 1, 0); 
				setSampleMaskedImage(newimage, x, y, 2, 0); 
			}
			ys+=2;
		}
		xs+=2;
	}

	return newimage;
}

// return a downsampled image (factor 1/2)
MaskedImage_P downsample(MaskedImage_P source)
{
	int kernel[6] = {1,5,10,10,5,1};
	int H, W;
	int x, y;
	int dx, dy;
	int xk, yk, ky, k;
	int r=0,g=0,b=0,m=0,ksum=0;
	H=source->image->height;
	W=source->image->width;
	int newW=W/2, newH=H/2;

	MaskedImage_P newimage = initNewMaskedImage(newW, newH);
    for (x=0;x<H-1;x+=2) {
        for (y=0;y<W-1;y+=2) {
			r=0; g=0; b=0; m=0; ksum=0;

            for (dy=-2;dy<=3;++dy) {
				yk=y+dy;
                if (yk<0 || yk>=W)
                    continue;
				ky = kernel[2+dy];
                for (dx=-2;dx<=3;++dx) {
					xk = x+dx;
                    if (xk<0 || xk>=H)
                        continue;

                    if (source->mask[xk][yk])
                        continue;

					k = kernel[2+dx]*ky;
					r+= k*getSampleMaskedImage(source, xk, yk, 0);
					g+= k*getSampleMaskedImage(source, xk, yk, 1);
					b+= k*getSampleMaskedImage(source, xk, yk, 2);
					ksum+=k;
					m++;
				}
			}
            if (ksum>0) {
                r/=ksum;
                g/=ksum;
                b/=ksum;
            }

			if (m!=0) {
				setSampleMaskedImage(newimage, x/2, y/2, 0, r); 
				setSampleMaskedImage(newimage, x/2, y/2, 1, g); 
				setSampleMaskedImage(newimage, x/2, y/2, 2, b); 
				setMask(newimage, x/2, y/2, 0);
			} else {
				setMask(newimage, x/2, y/2, 1);
			}
		}
	}

	return newimage;
}

// return an upscaled image
MaskedImage_P upscale(MaskedImage_P source, int newW,int newH)
{
	int x, y;
	int xs, ys;
	int H, W;
	H=source->image->height;
	W=source->image->width;
	MaskedImage_P newimage = initNewMaskedImage(newW, newH);

    for (x=0;x<newH;x++) {
        for (y=0;y<newW;y++) {

			// original pixel
			ys = (y*W)/newW;
			xs = (x*H)/newH;

			// copy to new image
			if (!source->mask[xs][ys]) {
				setSampleMaskedImage(newimage, x, y, 0, getSampleMaskedImage(source, xs, ys, 0)); 
				setSampleMaskedImage(newimage, x, y, 1, getSampleMaskedImage(source, xs, ys, 1)); 
				setSampleMaskedImage(newimage, x, y, 2, getSampleMaskedImage(source, xs, ys, 2)); 
				setMask(newimage, x, y, 0);
			} else {
				setSampleMaskedImage(newimage, x, y, 0, 0); 
				setSampleMaskedImage(newimage, x, y, 1, 0); 
				setSampleMaskedImage(newimage, x, y, 2, 0); 
				setMask(newimage, x, y, 1);
			}
		}
	}

	return newimage;
}
