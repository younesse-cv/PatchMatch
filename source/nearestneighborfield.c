#include "nearestneighborfield.h"

/**
* Nearest-Neighbor Field (see PatchMatch algorithm) 
*  This algorithme uses a version proposed by Xavier Philippeau
*
*/

NNF_P initNNF(MaskedImage_P input, MaskedImage_P output, int patchsize)
{
	NNF_P nnf = (NNF_P)malloc(sizeof(NNF_T));
	nnf->input = input;
	nnf->output= output;
	nnf->S = patchsize;
	nnf->field=NULL;

	return nnf;
}

void allocNNFField(NNF_P nnf)
{
	int i, j;
	if (nnf!=NULL){
		nnf->fieldH=nnf->input->image->height;
		nnf->fieldW=nnf->input->image->width;
		nnf->field = (int ***) malloc(sizeof(int**)*nnf->fieldH);

        for ( i=0 ; i < nnf->fieldH ; i++ ) {
			nnf->field[i] = (int **) malloc(sizeof(int*)*nnf->fieldW);
            for (j=0 ; j<nnf->fieldW ; j++ ) {
				nnf->field[i][j] = (int *) calloc(3,sizeof(int));
			}
		}
	}
}

void freeNNFField(NNF_P nnf)
{
	int i, j;
	if ( nnf->field != NULL ){
		for ( i=0 ; i < nnf->fieldH ; ++i ){
			for ( j=0 ; j < nnf->fieldW ; ++j ){
				free( nnf->field[i][j] );
			}
			free(nnf->field[i]);
		}
		free(nnf->field);
		nnf->field=NULL;
	}
}

void freeNNF(NNF_P nnf)
{
    if (nnf!=NULL) {
		freeNNFField(nnf);
		free(nnf);
		nnf=NULL;
	}
}

// initialize field with random values
void randomize(NNF_P nnf)
{
	int i, j;
	// field
	allocNNFField(nnf);
    for (i=0; i<nnf->input->image->height; ++i){
        for (j=0; j<nnf->input->image->width; ++j){
			nnf->field[i][j][0] = rand() % nnf->output->image->height +1;
			nnf->field[i][j][1] = rand() % nnf->output->image->width +1;
			nnf->field[i][j][2] = DSCALE;
		}
	}
	initializeNNF(nnf);
}

// initialize field from an existing (possibily smaller) NNF
void initializeNNFFromOtherNNF(NNF_P nnf, NNF_P otherNnf)
{
	int fx, fy, x, y, xlow, ylow;
	// field
	allocNNFField(nnf);
	fy = nnf->fieldW/otherNnf->fieldW;
	fx = nnf->fieldH/otherNnf->fieldH;
    for (x=0;x<nnf->fieldH;++x) {
        for (y=0;y<nnf->fieldW;++y) {
			xlow = int(min1(x/fx, otherNnf->input->image->height-1));
			ylow = int(min1(y/fy, otherNnf->input->image->width-1));
			nnf->field[x][y][0] = otherNnf->field[xlow][ylow][0]*fx;  
			nnf->field[x][y][1] = otherNnf->field[xlow][ylow][1]*fy;
			nnf->field[x][y][2] = DSCALE;
		}
	}
	initializeNNF(nnf);
}

// compute initial value of the distance term
void initializeNNF(NNF_P nnf)
{
	int y, x;
	int iter=0, maxretry=20;
    for (x=0;x<nnf->fieldH;++x) {
        for (y=0;y<nnf->fieldW;++y) {

			nnf->field[x][y][2] = distanceNNF(nnf, x,y,  nnf->field[x][y][0],nnf->field[x][y][1]);
			// if the distance is INFINITY (all pixels masked ?), try to find a better link
			iter=0;
            while ( nnf->field[x][y][2] == DSCALE && iter<maxretry) {
				nnf->field[x][y][0] = rand() % nnf->output->image->height +1;
				nnf->field[x][y][1] = rand() % nnf->output->image->width +1;
				nnf->field[x][y][2] = distanceNNF(nnf, x,y,  nnf->field[x][y][0],nnf->field[x][y][1]);
				iter++;
			}
		}
	}
}

// multi-pass NN-field minimization (see "PatchMatch" - page 4)
void minimizeNNF(NNF_P nnf, int pass)
{
	int i, y, x;
	int min_x=0, min_y=0, max_y=nnf->input->image->width-1, max_x=nnf->input->image->height-1;
	// multi-pass minimization
    for (i=0;i<pass;i++) {
		// scanline order
        for (x=min_x;x<max_x;++x)
            for (y=min_y;y<=max_y;++y)
                if (nnf->field[x][y][2]>0)
                    minimizeLinkNNF(nnf, x,y,+1);

		// reverse scanline order
        for (x=max_x;x>=min_x;x--)
            for (y=max_y;y>=min_y;y--)
                if (nnf->field[x][y][2]>0)
                    minimizeLinkNNF(nnf, x,y,-1);
	}
}

// minimize a single link (see "PatchMatch" - page 4)
void minimizeLinkNNF(NNF_P nnf, int x, int y, int dir)
{
	int xp,yp,dp,wi, xpi, ypi;
	//Propagation Up/Down
	if (x-dir>0 && x-dir<nnf->input->image->height) {
		xp = nnf->field[x-dir][y][0]+dir;
		yp = nnf->field[x-dir][y][1];
		dp = distanceNNF(nnf,x,y, xp,yp);
		if (dp<nnf->field[x][y][2]) {
			nnf->field[x][y][0] = xp;
			nnf->field[x][y][1] = yp;
			nnf->field[x][y][2] = dp;
		}
	}
	//Propagation Left/Right
	if (y-dir>0 && y-dir<nnf->input->image->width) {
		xp = nnf->field[x][y-dir][0];
		yp = nnf->field[x][y-dir][1]+dir;
		dp = distanceNNF(nnf,x,y, xp,yp);
		if (dp<nnf->field[x][y][2]) {
			nnf->field[x][y][0] = xp;
			nnf->field[x][y][1] = yp;
			nnf->field[x][y][2] = dp;
		}
	}
	//Random search
	wi=nnf->output->image->width;
	xpi=nnf->field[x][y][0];
	ypi=nnf->field[x][y][1];
	int r=0;
    while (wi>0) {
		r=(rand() % (2*wi)) -wi;
		xp = xpi + r;
		r=(rand() % (2*wi)) -wi;
		yp = ypi + r;
		xp = int(max1(0, min1(nnf->output->image->height-1, xp )));
		yp = int(max1(0, min1(nnf->output->image->width-1, yp )));

		dp = distanceNNF(nnf,x,y, xp,yp);
		if (dp<nnf->field[x][y][2]) {
			nnf->field[x][y][0] = xp;
			nnf->field[x][y][1] = yp;
			nnf->field[x][y][2] = dp;
		}
		wi/=2;
	}
}

// compute distance between two patch 
int distanceNNF(NNF_P nnf, int x,int y, int xp,int yp)
{
	return distanceMaskedImage(nnf->input,x,y, nnf->output,xp,yp, nnf->S);
}

