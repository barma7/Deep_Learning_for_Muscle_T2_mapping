#include <math.h>
#include "mex.h"
#include "matrix.h"

void epg(double *Sig, double T1, double T2, double esp, double alpha[], mwSize etl, mwSize runs) {
    
    /*  Define basis states  */
    mwSize nk = 2*etl+1;
    
    /*  Define indexing variables  */
    mwSize h,i,j,k=0;
    
    /* mwSize jn[nk]; */
    
    /*  Define computation variables (scalars)  */
    double E1 = exp(-0.5*esp/T1);
    double E2 = exp(-0.5*esp/T2);
    double ca, sa, cahp, cahm, sah;
    
    /*  Define computation variables (matrices)  */
    double *Z  = (double*)mxCalloc(nk,sizeof(double));
    double *F  = (double*)mxCalloc(nk,sizeof(double));
    double *Zk = (double*)mxCalloc(nk,sizeof(double));
    double *Fk = (double*)mxCalloc(nk,sizeof(double));
    mwSize *jn = (mwSize*)mxCalloc(2*etl+1,sizeof(mwSize));
    
    /*  Define a k-state index  */
    
    for (i = 0; i < nk; i++) {
        jn[i] = nk-1 - i;
    }
    
    /*  Loop through each run (could parallelize)  */
    for (h=0; h < runs; h++) {
        
        /*  Reset computational arrays  */
        for (j=0; j < nk; j++)
            Z[j] = 0.0;
        for (j=0; j < nk; j++)
            F[j] = 0.0;
        for (j=0; j < nk; j++)
            Zk[j] = 0.0;
        for (j=0; j < nk; j++)
            Fk[j] = 0.0;
        
        /*  Set transverse magnetization following 90deg. excitation  */
        /*  then one precession/decay period                          */
        F[nk/2-1] = E2;
        
        /*  Loop through echo train  */
        for (i=0; i < etl; i++) {
             
            /*  Apply nutation to all k-states  */
            ca   = cos(alpha[k]);
            cahp = 0.5+0.5*ca;
            cahm = 0.5-0.5*ca;
            sa   = sin(alpha[k]);
            sah  = 0.5*sa;
            for (j=0; j < nk; j++)
                Fk[j] = cahp*F[j] + cahm*F[jn[j]] + sa*Z[j];
            for (j=0; j < nk; j++)
                Zk[j] = -sah*F[j] + sah*F[jn[j]] + ca*Z[j];
            
            /*  Apply first precession and relaxation  */
            /*  (pre-echo)                             */
            for (j=0; j < nk-1; j++)
                F[j] = Fk[j+1] * E2;
            for (j=0; j < nk-1; j++)
                Z[j] = Zk[j] * E1;
            
            /*  Get F(k==0) echo amplitudes  */
            Sig[k++] = F[nk/2];
            
            /*  Apply second precession and relaxation  */
            /*  (post-echo)                             */
            for (j=0; j < nk-1; j++)
                F[j] = F[j+1] * E2;
            for (j=0; j < nk-1; j++)
                Z[j] = Z[j] * E1;
        }
    }
    return;
}

void mexFunction( int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[] )
{
    /*  Initialize basic pointers  */
    mwSize mrows,ncols;
    
    /*  Assign inputs  */
    double T1     = mxGetScalar(prhs[0]);
    double T2     = mxGetScalar(prhs[1]);
    double esp    = mxGetScalar(prhs[2]);
    double *alpha = mxGetPr(prhs[3]);
    
    /* Declare output pointer */
    double *Sig;
    
    /*  Check for proper number of arguments.  */
    if(nrhs!=4) {
        mexErrMsgTxt("Four inputs required.");
    } else if(nlhs>1) {
        mexErrMsgTxt("Too many output arguments");
    }
    
    /*  Get size of flip angles  */
    mrows = mxGetM(prhs[3]);
    ncols = mxGetN(prhs[3]);
    
    /*  Create matrix for the return argument  */
    plhs[0] = mxCreateDoubleMatrix(mrows,ncols,mxREAL);
    
    /*  Assign pointers for signal output  */
    Sig = mxGetPr(plhs[0]);
    
    /*  Call epg computation engine  */
    epg(Sig,T1,T2,esp,alpha,mrows,ncols);
}
