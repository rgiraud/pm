//rebuilt_mex.c //
//
//Rebuilts the image a with the stack b and the NNF map returned by the
//patchmatch algorithm. Block wise algorithm
//
//
//INPUTS: - b: Image B to extract patches
//        - nnf: NNF map returned by the PatchMatch algorithm
//        - patch_w: patch size
//
//OUTPUTS: - a_rebuilt: Image A rebuilt from patches in B


#include "math.h"
#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>



//Rebuilt function
void rebuilt_2D(float *b, int*nnf, float *a_rebuilt, int ha, int wa,
        int hb, int d, int sizea, int sizeb, int patch_w, int patch_bw){
    
    
    float *table_sum = (float*) calloc(sizea,sizeof(float));
    
    //Every voxel treatment
    int i, j, bx, by, pos_current, dx, dy, dc;
    
    for (j=patch_w; j<wa-patch_w; j++) {
        for (i=patch_w; i<ha-patch_w; i++) {
            
            pos_current = i + j*ha;
            
            bx = nnf[pos_current];
            by = nnf[pos_current + sizea];
            
            for (dx = -patch_bw; dx <= patch_bw; dx++) {
                for (dy =-patch_bw; dy <= patch_bw; dy++) {
                    
                    for (dc = 0; dc < d; dc++) {
                        a_rebuilt[pos_current+dy+dx*ha+sizea*dc]
                                += b[by+dy+(bx+dx)*hb +sizeb*dc];
                    }
                    
                    table_sum[pos_current+dy+dx*ha] += 1;
                    
                }
            }
        }
        
    }
    
    //Normalization
    for (i=0; i<sizea; i++) {
        for (dc=0; dc<d; dc++) {
            a_rebuilt[i+dc*sizea] = (float) a_rebuilt[i+dc*sizea]/table_sum[i];
        }
    }
    
    free(table_sum);
    
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[]) {
    
    //INPUTS
    /*Image loading*/
    float* b = (float*)mxGetPr(prhs[0]);
    int* nnf = (int*)mxGetPr(prhs[1]);
    
    /*Parameters*/
    int patch_w = mxGetScalar(prhs[2]);
    int patch_bw = patch_w;
    
    
    //B dimensions
    int ndimb = mxGetNumberOfDimensions(prhs[0]);
    const int *dimsb = mxGetDimensions(prhs[0]);
    
    //A dimensions
    const int *dimsnnf = mxGetDimensions(prhs[1]);
    
    int ha, wa, hb, wb, d, i;
    int dimsa[3];
    int sizea, sizeb;
    
    float* a_rebuilt = NULL;
    float* pta_rebuilt = NULL;
    
    
    d = 1; //channel number
    if (ndimb == 3) {
        d = dimsb[2];
    }
    
    dimsa[0] = dimsnnf[0];
    dimsa[1] = dimsnnf[1];
    dimsa[2] = d;
    
    
    //Sizes
    ha = dimsa[0];
    wa = dimsa[1];
    hb = dimsb[0];
    wb = dimsb[1];
    
    
    printf("%d, %d, %d, %d, %d, %d\n", ha, wa, hb, wb, d, ndimb);
    mexEvalString("drawnow");
    
    
    sizea = ha*wa;
    sizeb = hb*wb;
    
    // A rebuilt
    plhs[0] = mxCreateNumericArray(ndimb, dimsa, mxGetClassID(prhs[0]),
            mxREAL);
    a_rebuilt = (float*)mxGetPr(plhs[0]);
    pta_rebuilt = a_rebuilt;
    for(i=0; i<sizea*d; i++)
        *pta_rebuilt++ = 0.0;
    
    
    rebuilt_2D(b, nnf, a_rebuilt, ha, wa, hb, d, sizea, sizeb, patch_w, patch_bw);
    
    
    
    
}