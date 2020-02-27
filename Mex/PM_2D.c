//PM_2D.c //
//
//Computes for all patches in an image A, correspondences in an image B.
//
//
//INPUTS: - a: image A
//        - b: image B
//        - patch_w: patch size
//        - iter: iteration number of PatchMatch inner algorithm
//
//OUTPUTS: - nnf: nearest neighbor field map storing positions of matches
//         - nnfd: distance associated to the selected matches


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mex.h"
#include "matrix.h"

#ifndef MAX
#define MAX(a, b) ((a)>(b)?(a):(b))
#endif

#ifndef MIN
#define MIN(a, b) ((a)<(b)?(a):(b))
#endif

#define EPSILON 0.0000001

#ifndef enc_type
#define enc_type unsigned char
#endif



// Measure distance between entire 2 patches,
// terminating early if we exceed a cutoff distance
float dist_2D(float *a, float *b, int ax, int ay, int bx, int by,
        int ha, int  d, int hb, int size_a, int size_b, int patch_w,
        float cutoff) {
    
    float ans = 0;
    int dx, dy, dz;
    float ac = 0, dc = 0, bc =0;
    
    for (dx = -patch_w; dx <= patch_w; dx++) {
        for (dy = -patch_w; dy <= patch_w; dy++) {
            for (dz=0; dz<d; dz++) {
                
                ac = a[(ax+dx)*ha+ay+dy+dz*size_a];
                bc = b[(bx+dx)*hb+by+dy+dz*size_b];
                dc = ac - bc;
                
                ans += dc*dc;
            }
        }
        if (ans >= cutoff) { return cutoff; }
    }
    
    
    return ans;
}



// Compare two distance and change the optimal ptr if the new distance is
// lower than the previous best 
void improve_guess_2D(float *a, float *b, int ax, int ay,
        int *xbest_ptr, int *ybest_ptr, float *distbest_ptr,
        int bx, int by, int ha, int d, int hb, int patch_w,
        int size_a, int size_b) {
    
    //Distance computation
    float dist = dist_2D(a, b, ax, ay, bx, by, ha, d, hb, size_a, size_b,
            patch_w, *distbest_ptr);
    
    //Comparison
    if (dist < *distbest_ptr) {
        *distbest_ptr = dist;
        *xbest_ptr = bx;
        *ybest_ptr = by;
    }
    
}



void  pm_hist_rgb(int*nnf_out, float* nnfd_out, float * a, float * b ,
        int ha, int wa, int hb, int wb,  int d, int patch_w, int limit) {
    
    // Initialize with random nearest neighbor field (NNF)
    // Effective width and height (possible upper left corners of patches)
    int aew = wa-patch_w, aeh = ha-patch_w;
    int bew = wb-patch_w, beh = hb-patch_w;
    
    int size_a = ha*wa;
    int size_b = hb*wb;
    
    int pos, pos_shift;
    int ax, ay, bx, by, ax_off;
    int iter = 0;
    int ystart, yend, ychange, xstart, xend, xchange;
    int xbest, ybest;
    float distbest;
    int distbest_temp;
    float *distbest_ptr;
    int *xbest_ptr;
    int *ybest_ptr;
    int mag = 100;
    int rand_val;
    int xp, yp, xmin, xmax, ymin, ymax;
    int i;
    int * nnf = (int *) calloc(ha*wa*2, sizeof(int));
    float * nnfd = (float *) calloc(ha*wa, sizeof(float));
    
    
    // RANDOM INITIALIZATION
    for (ax = patch_w; ax < aew; ax++) {
        ax_off =ha*ax;
        for (ay = patch_w; ay < aeh; ay++) {
            
            pos = ax_off+ay;
            
            xmin = MAX(ax-mag, patch_w);
            xmax = MIN(ax+mag+1, bew-1);
            
            ymin = MAX(ay-mag, patch_w);
            ymax = MIN(ay+mag+1, beh-1);
            
            
            //random match
            rand_val = rand();
            bx = xmin+rand_val%(xmax-xmin);
            rand_val = rand();
            by = ymin+rand_val%(ymax-ymin);
            
            //Map init
            nnf[pos] = bx;
            nnf[pos + size_a] = by;
            nnfd[pos] = (float) dist_2D(a, b, ax, ay, bx, by, ha, d, hb, size_a, size_b, patch_w, FLT_MAX);
            
        }
    }
    
    //Recopy out
    for (i=0; i<ha*wa; i++) {
        nnfd_out[i] = nnfd[i];
        nnf_out[i] = nnf[i];
    }
    for (i=ha*wa; i<ha*wa*2; i++) {
        nnf_out[i] = nnf[i];
    }
    
    while (iter < limit) {
        
        // In each iteration, improve the NNF, by looping in scanline or
        // reverse-scanline order
        ystart = patch_w+1; yend = aeh; ychange = 1;
        xstart = patch_w+1; xend = aew; xchange = 1;
        
        if (iter % 2 == 1) {
            xstart = xend-1; xend = patch_w-1; xchange = -1;
            ystart = yend-1; yend = patch_w-1; ychange = -1;
        }
        
        for (ay = ystart; ay != yend; ay += ychange) {
            for (ax = xstart; ax != xend; ax += xchange) {
                
                //Current position in a
                pos = ha*ax+ay;
                
                // Current (best) guess
                xbest = nnf[pos];
                ybest = nnf[pos + size_a];
                distbest = nnfd[pos];
                distbest_temp = distbest;
                distbest_ptr = &distbest;
                xbest_ptr = &xbest;
                ybest_ptr = &ybest;
                
                // PROPAGATION: Improve current guess by trying instead
                // correspondences from left and above (below and right on
                // odd iterations)
                
                //XSHIFT
                if ((unsigned) (ax - xchange-patch_w) < (unsigned) aew-patch_w) {
                    pos_shift = ha*(ax-xchange)+ay;
                    xp = nnf[pos_shift]  + xchange;
                    yp = nnf[pos_shift + size_a];
                    if ((unsigned) xp-patch_w < (unsigned) bew-patch_w) {
                        if ((xp !=  *xbest_ptr)  || (yp !=  *ybest_ptr))
                            improve_guess_2D(a, b, ax, ay, xbest_ptr, ybest_ptr,
                                    distbest_ptr, xp, yp, ha, d,  hb, patch_w, size_a, size_b);
                    }
                }
                
                
                //YSHIFT
                if ((unsigned) (ay - ychange-patch_w) < (unsigned) aeh-patch_w) {
                    pos_shift = ha*ax+ay -ychange;
                    xp = nnf[pos_shift];
                    yp = nnf[pos_shift + size_a] + ychange;
                    if ((unsigned) yp-patch_w < (unsigned) beh-patch_w) {
                        if ((xp !=  *xbest_ptr)  || (yp !=  *ybest_ptr))
                            
                            improve_guess_2D(a, b, ax, ay, xbest_ptr, ybest_ptr,
                                    distbest_ptr, xp, yp, ha, d,  hb, patch_w, size_a, size_b);
                        
                    }
                }
                
                // RANDOM SEARCH: Improve current guess by searching in
                // boxes of exponentially decreasing size around the
                // current position
                
                // Sampling window
                for (mag = MAX(hb,wb); mag >= 1; mag /= 2) {
                    
                    xbest = *xbest_ptr;
                    ybest = *ybest_ptr;
                    
                    xmin = MAX(xbest-mag, patch_w);
                    xmax = MIN(xbest+mag+1, bew-1);
                    if(xmin == xmax) continue;
                    
                    ymin = MAX(ybest-mag, patch_w);
                    ymax = MIN(ybest+mag+1, beh-1);
                    if(ymin == ymax) continue;
                    
                    //Random match
                    xp = (int) xmin+rand()%MAX(xmax-xmin, 1);
                    yp = (int) ymin+rand()%MAX(ymax-ymin, 1);
                    
                    if (((unsigned) yp-patch_w < (unsigned) beh)
                    && ((unsigned) xp-patch_w < (unsigned) bew)) {
                        improve_guess_2D(a, b, ax, ay, xbest_ptr, ybest_ptr,
                                distbest_ptr, xp, yp, ha, d,  hb, patch_w, size_a, size_b);
                    }
                }
                
                //Map updates
                nnf[pos] = *xbest_ptr;
                nnf[pos + size_a] =*ybest_ptr;
                nnfd[pos] = *distbest_ptr;
                
            }
            
        }
        
        for (i=0; i<ha*wa; i++) {
            nnfd_out[i + (iter+1)*size_a] = nnfd[i];
            nnf_out[i + (iter+1)*size_a*2] = nnf[i];
        }
        for (i=ha*wa; i<ha*wa*2; i++) {
            nnf_out[i + (iter+1)*size_a*2] = nnf[i];
        }
        iter++;
    }
    
}










//////////////////////////////////////////////////////////////////////////
/////////////////////////////////// MAIN /////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[]){
    
    float * a;
    float * b;
    const int * a_dims;
    const int * b_dims;
    int ha, wa, hb, wb, d;
    int dims[4];
    int * nnf;
    float * nnfd;
    int idx = 2;
    int patch_w, iter;
    int img_rgb_dims;
    
    a = (float*) mxGetPr(prhs[0]);
    b = (float*) mxGetPr(prhs[1]);
    a_dims = mxGetDimensions(prhs[0]);
    ha = a_dims[0];
    wa = a_dims[1];
    b_dims = mxGetDimensions(prhs[1]);
    hb = b_dims[0];
    wb = b_dims[1];
    
    img_rgb_dims = mxGetNumberOfDimensions(prhs[1]);
    
    patch_w = (int)  mxGetScalar(prhs[idx++]);
    iter = (int) mxGetScalar(prhs[idx++]);
    
    dims[0] = ha;
    dims[1] = wa;
    dims[2] = 2;
    dims[3] = iter+1;
    plhs[0] = mxCreateNumericArray(4, dims, mxINT32_CLASS, mxREAL);
    nnf = (int*)mxGetPr(plhs[0]);
    
    dims[2] = iter+1;
    plhs[1] = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL);
    nnfd = (float*)mxGetPr(plhs[1]);
    
    if (img_rgb_dims == 3)
        d = a_dims[2];
    else
        d = 1;
    
    pm_hist_rgb(nnf, nnfd, a, b, ha, wa, hb, wb, d, patch_w, iter);
    
    
}