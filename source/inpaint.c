#include "inpaint.h"

/** This algorithme uses a version proposed by Xavier Philippeau
*/

Inpaint_P initInpaint()
{
    Inpaint_P inp = (Inpaint_P)malloc(sizeof(Inpaint_T));
    //initial image
    inp->initial = NULL;

    // Nearest-Neighbor Fields
    inp->nnf_SourceToTarget = NULL;
    inp->nnf_TargetToSource = NULL;

    // Pyramid of downsampled initial images
    inp->pyramid = NULL;
    inp->nbEltPyramid = 0;
    inp->nbEltMaxPyramid = 0;

    return inp;
}

void addEltInpaintingPyramid(Inpaint_P imp, MaskedImage_P elt)
{
    int inc = INCREASE_PYRAMID_SIZE_RATE;
    if (inc<2)
        inc = 2;

    if (imp->pyramid == NULL || imp->nbEltMaxPyramid == 0) {
        imp->nbEltMaxPyramid = inc;
        imp->pyramid = (MaskedImage_P*)malloc(sizeof(MaskedImage_P)*imp->nbEltMaxPyramid);
    } else if (imp->nbEltPyramid == imp->nbEltMaxPyramid) {
        imp->nbEltMaxPyramid = imp->nbEltMaxPyramid*inc;
        imp->pyramid = (MaskedImage_P*)realloc(imp->pyramid, sizeof(MaskedImage_P)*imp->nbEltMaxPyramid);
    }

    imp->pyramid[imp->nbEltPyramid] = elt;
    imp->nbEltPyramid++;
}

// EM-Like algorithm (see "PatchMatch" - page 6)
// Returns a double sized target image
MaskedImage_P ExpectationMaximization(Inpaint_P imp, int level)
{
    int emloop, x, y, H, W, i, j;
    double*** vote;

    int iterEM = 1+2*level;
    int iterNNF = int(min1(7,1+level));

    int upscaled;
    MaskedImage_P newsource;
    MaskedImage_P source = imp->nnf_SourceToTarget->input;
    MaskedImage_P target = imp->nnf_SourceToTarget->output;
    MaskedImage_P newtarget = NULL;

    CvSize size;

    printf("EM loop (em=%d,nnf=%d) : ", iterEM, iterNNF);

    // EM Loop
    for (emloop=1;emloop<=iterEM;emloop++) {
        printf(" %d", 1+iterEM-emloop);
        // set the new target as current target
        if (newtarget!=NULL) {
            imp->nnf_SourceToTarget->output = newtarget;
            imp->nnf_TargetToSource->input = newtarget;
            target = newtarget;
            newtarget = NULL;
        }
        // -- add constraint to the NNF
        H = source->image->height;
        W = source->image->width;
        // we force the link between unmasked patch in source/target
        for ( x=0 ; x<H; ++x)
            for ( y=0 ; y<W ; ++y)
                if (!constainsMasked(source, x, y, imp->radius)) {
                    imp->nnf_SourceToTarget->field[x][y][0] = x;
                    imp->nnf_SourceToTarget->field[x][y][1] = y;
                    imp->nnf_SourceToTarget->field[x][y][2] = 0;
                }

        H = target->image->height;
        W = target->image->width;
        for ( x=0 ; x<H ; ++x)
            for ( y=0 ; y<W ; ++y)
                if (!constainsMasked(source, x, y, imp->radius)) {
                    imp->nnf_TargetToSource->field[x][y][0] = x;
                    imp->nnf_TargetToSource->field[x][y][1] = y;
                    imp->nnf_TargetToSource->field[x][y][2] = 0;
                }
        // -- minimize the NNF
        minimizeNNF(imp->nnf_SourceToTarget, iterNNF);
        minimizeNNF(imp->nnf_TargetToSource, iterNNF);

        // -- Now we rebuild the target using best patches from source

        upscaled = 0;

        // Instead of upsizing the final target, we build the last target from the next level source image
        // So the final target is less blurry (see "Space-Time Video Completion" - page 5)
        if (level>=1 && (emloop==iterEM)) {
            newsource = imp->pyramid[level-1];
            newtarget = upscale(target, newsource->image->width,newsource->image->height);
            upscaled = 1;
        } else {
            newsource = imp->pyramid[level];
            newtarget = copyMaskedImage(target);
            upscaled = 0;
        }

        // --- EXPECTATION STEP ---

        // votes for best patch from NNF Source->Target (completeness) and Target->Source (coherence)

        vote = (double ***)malloc(newtarget->image->height*sizeof(double **));
        for ( i=0 ; i<newtarget->image->height ; ++i ){
            vote[i] = (double **)malloc(newtarget->image->width*sizeof(double *));
            for  ( j=0 ; j<newtarget->image->width ; ++j) {
                vote[i][j] = (double *)calloc(4, sizeof(double));
            }
        }

        ExpectationStep(imp->nnf_SourceToTarget, 1, vote, newsource, upscaled);
        ExpectationStep(imp->nnf_TargetToSource, 0, vote, newsource, upscaled);

        // --- MAXIMIZATION STEP ---

        // compile votes and update pixel values
        MaximizationStep(newtarget, vote);

        size = cvSize(imp->initial->image->width,imp->initial->image->height);
        IplImage* result=cvCreateImage(size,IPL_DEPTH_8U,3);

        cvResize(newtarget->image,result,CV_INTER_LINEAR);

        for (i=0;i<newtarget->image->height;i++) {
            for (j=0;j<newtarget->image->width;j++)
                free(vote[i][j]);

            free(vote[i]);
        }
        free(vote);
    }

    printf("\n");

    return newtarget;
}


// Expectation Step : vote for best estimations of each pixel
void ExpectationStep(NNF_P nnf, int sourceToTarget, double*** vote, MaskedImage_P source, int upscale)
{
    int y, x, H, W, Ho, Wo, xp, yp, dp, dy, dx;
    int xs,ys,xt,yt;
    int*** field = nnf->field;
    int R = nnf->S; /////////////int R = nnf->PatchSize;
    double w;

    H = nnf->input->image->height;
    W = nnf->input->image->width;
    Ho = nnf->output->image->height;
    Wo = nnf->output->image->width;
    for ( x=0 ; x<H ; ++x) {
        for ( y=0 ; y<W; ++y) {
            // x,y = center pixel of patch in input

            // xp,yp = center pixel of best corresponding patch in output
            xp=field[x][y][0];
            yp=field[x][y][1];
            dp=field[x][y][2];

            // similarity measure between the two patches
            w = G_globalSimilarity[dp];

            // vote for each pixel inside the input patch
            for ( dy=-R ; dy<=R ; ++dy) {
                for ( dx=-R ; dx<=R; ++dx) {

                    // get corresponding pixel in output patch
                    if (sourceToTarget)
                    { xs=x+dx; ys=y+dy;	xt=xp+dx; yt=yp+dy;}
                    else
                    { xs=xp+dx; ys=yp+dy; xt=x+dx; yt=y+dy; }

                    if (xs<0 || xs>=H) continue;
                    if (ys<0 || ys>=W) continue;
                    if (xt<0 || xt>=Ho) continue;
                    if (yt<0 || yt>=Wo) continue;

                    // add vote for the value
                    if (upscale) {
                        weightedCopy(source, 2*xs,   2*ys,   vote, 2*xt,   2*yt,   w);
                        weightedCopy(source, 2*xs+1, 2*ys,   vote, 2*xt+1, 2*yt,   w);
                        weightedCopy(source, 2*xs,   2*ys+1, vote, 2*xt,   2*yt+1, w);
                        weightedCopy(source, 2*xs+1, 2*ys+1, vote, 2*xt+1, 2*yt+1, w);
                    } else {
                        weightedCopy(source, xs, ys, vote, xt, yt, w);
                    }
                }
            }
        }
    }
}

void weightedCopy(MaskedImage_P src, int xs, int ys, double*** vote, int xd,int yd, double w)
{
    if (isMasked(src, xs, ys))
        return;

    vote[xd][yd][0] += w*getSampleMaskedImage(src, xs, ys, 0);
    vote[xd][yd][1] += w*getSampleMaskedImage(src, xs, ys, 1);
    vote[xd][yd][2] += w*getSampleMaskedImage(src, xs, ys, 2);
    vote[xd][yd][3] += w;
}


// Maximization Step : Maximum likelihood of target pixel
void MaximizationStep(MaskedImage_P target, double*** vote)
{
    int y, x, H, W, r, g, b;
    H = target->image->height;
    W = target->image->width;
    for( x=0 ; x<H ; ++x){
        for( y=0 ; y<W ; ++y){
            if (vote[x][y][3]>0) {
                r = (int) (vote[x][y][0]/vote[x][y][3]);
                g = (int) (vote[x][y][1]/vote[x][y][3]);
                b = (int) (vote[x][y][2]/vote[x][y][3]);

                setSampleMaskedImage(target, x, y, 0, r );
                setSampleMaskedImage(target, x, y, 1, g );
                setSampleMaskedImage(target, x, y, 2, b );
                setMask(target, x, y, 0);
            }


        }
    }
}
IplImage* inpaint(Inpaint_P imp, IplImage* input, int ** mask, int radius)
{
    int level, y, x;
    NNF_P new_nnf, new_nnf_rev;
    // initial image
    imp->initial = initMaskedImage(input, mask);
    IplImage* tmp = NULL;
    CvSize size;
    size = cvSize(imp->initial->image->width,imp->initial->image->height);

    // patch radius
    imp->radius = radius;

    // working copies
    MaskedImage_P source = imp->initial;
    MaskedImage_P target = NULL;

    printf("Build pyramid of images...\n");

    // build pyramid of downscaled images
    addEltInpaintingPyramid(imp, source);
    while (source->image->width>radius && source->image->height>radius) {
        source = downsample(source);
        addEltInpaintingPyramid(imp, source);
    }
    int maxlevel=imp->nbEltPyramid;

    // for each level of the pyramid
    cvNamedWindow("Progression", CV_WINDOW_AUTOSIZE);
    for (level=maxlevel-1 ; level>0 ; level--) {
        printf( "\n*** Processing -  Zoom 1:%d ***" , 1<<level );

        // create Nearest-Neighbor Fields (direct and reverse)
        source = imp->pyramid[level];

        printf("initialize NNF...");
        if (level==maxlevel-1) {
            // at first, we use the same image for target and source
            // and use random data as initial guess
            target = copyMaskedImage(source);

            // we consider that the target contains no masked pixels in the firt time

            for( x = 0 ; x < target -> image -> height ; x++ )
                for( y = 0 ; y < target -> image -> width ; y++ )
                    setMask(target, x, y, 0);

            imp->nnf_SourceToTarget = initNNF(source, target, radius);
            randomize(imp->nnf_SourceToTarget);

            imp->nnf_TargetToSource = initNNF(target, source, radius);
            randomize(imp->nnf_TargetToSource);

        } else {
            // then, we use the rebuilt (upscaled) target
            // and re-use the previous NNF as initial guess
            new_nnf = initNNF(source, target, radius);
            initializeNNFFromOtherNNF(new_nnf, imp->nnf_SourceToTarget);
            imp->nnf_SourceToTarget = new_nnf;

            new_nnf_rev = initNNF(target, source, radius);
            initializeNNFFromOtherNNF(new_nnf_rev, imp->nnf_TargetToSource);
            imp->nnf_TargetToSource = new_nnf_rev;
        }

        target = ExpectationMaximization(imp, level);
        if (tmp != NULL )
            cvReleaseImage(&tmp);

        tmp=cvCreateImage(size,IPL_DEPTH_8U,3);
        cvResize(target->image,tmp,CV_INTER_LINEAR);

        cvShowImage("Progression",tmp);
        cvMoveWindow("Progression",750, 100);
        cvWaitKey(1);
    }


    return target->image;
    cvReleaseImage(&tmp);
}

void freeInpaintingPyramid(Inpaint_P inp)
{
    int i;
    if (inp->pyramid != NULL) {
        for ( i=0 ; i<inp->nbEltPyramid ; ++i)
            freeMaskedImage(inp->pyramid[i]);

        free(inp->pyramid);
        inp->pyramid = NULL;
    }
}
