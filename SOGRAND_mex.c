/*
* Guessing Random Additive Noise Decoding (GRAND)
* All code is subject to license:
* GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf
*/

/* SOGRAND: 1-line ORBGRAND with list output and block soft output
*  Optional argument allows basic ORBGRAND or a fixed offset.
*  
* Install: mex -O SOGRAND.c
*
* [chat, score, NT, curL, pNL, APP] = 1lineORB(llr, H, IC, L, Tmax, ConfMin, even);
*
* n        :                 codelength
*
* chat     :(n,L) matrix,    L best candidates c
* score    :(4,L) matrix,    LW, PM, #guess, listConf of the candidates
* NT       :scalar,          Number of guesses
* curL     :scalar,          real list size
* pNL      :scalar,          probability of the transmitted codeword is not in the list
* APP      :scalar,          probability of each decoding being correct
*
* LLR      :(n,1) matrix,    log likelihood ratio of the channel outputs
* H        :(n*s,1) matrix,  Parity-check matrix
* IC       :scalar,          intercept for 1-line ORBGRAND (IC<0 for automatic, 0 for basic ORBGRAND)
* L        :scalar,          list size of the list decoder
* Tmax     :scalar,          Maximum number of guesses
* ConfMin  :scalar,          Required listConf
* even     :scalar,          1 if code is even; 0 otherwise
*
* References to be cited when using this code:
* [1] K. R. Duffy, J. Li, M. Médard. "Capacity-Achieving Guessing Random Additive Noise Decoding". IEEE Transactions on Information Theory 65 (2019): 4023-4040.
* [2] K. R. Duffy, W. An, M. Médard. "Ordered reliability bits guessing random additive noise decoding." IEEE Transactions on Signal Processing 70 (2022): 4528-4542.
* [3] K. Galligan, P. Yuan, M. Médard, K. R. Duffy. "Upgrade error detection to prediction with GRAND". IEEE Globecom (2023).
* [4] P. Yuan, M. Medard, K. Galligan, K. R. Duffy. "Soft-output (SO) GRAND and Iterative Decoding to Outperform LDPC Codes". IEEE Transactions on Wireless Communications (2025).
*/

/* C Functions */
#include <mex.h>
#include <math.h>
#include <stdint.h>
#include <stddef.h>
#define Inf 0x7fffffff

/* Parity_Check */
uint8_t ParityCheck(uint8_t *c, uint8_t *H, uint64_t n, uint64_t s) {
    uint8_t syndrome;
    for (size_t j = 0; j < s; j++){
        syndrome = 0;
        for (size_t i = 0; i < n; i++)
            syndrome ^= ( c[i] * H[j*n + i] );
        if (syndrome == 1)
            return 0;
    }
    return 1;
}

/* Find Max */
int32_t findMax(int32_t a, int32_t b) {
    return !(b > a) ? a : b;
}

/* Hard Decision */
void HardDec(uint8_t *c, double *llr, uint64_t n) {
    for (size_t i = 0; i < n; i++){
        if (llr[i] > 0.0)
            c[i] = 0;
        else
            c[i] = 1;
    }
}

/* Determine parity of binary sequence */
int parity(uint8_t array[], uint64_t n)
{
    int i;
    int sum = 0;
    for (i = 0; i < n; i++)
    {
        sum = sum + array[i];
    }
    sum = sum % 2;
    return sum;
}

/* Determine probability Bernoulli sequence have given parity */
double prob_parity(int parity_cHD, double *absL, uint64_t n)
{
    int i;
    int sum = 0;
    double prob_even = 1;

    for (i = 0; i < n; i++)
    {
        prob_even = prob_even*(1-2*exp(-absL[i])/(1+exp(-absL[i]))); 
    }
    prob_even = 0.5*(1.0+prob_even);
    if (parity_cHD==0)
         return prob_even;
    prob_even = 1-prob_even;
    return prob_even;
}

/* Add TEP */
void AddTEP(uint8_t *c, uint8_t *cHD, uint8_t *TEP, size_t *perm, uint64_t n) {
    for (size_t i = 0; i < n; i++)
        c[perm[i]] = cHD[perm[i]] ^ TEP[i];
}

/* Jacobian logarithm (approximated) */
double JacLog(double x) {
    if(x > 50)
        return x;
    else if(x < -50)
        return 0.0;
    else
        return log(1+exp(x));
}

/* Quick Sort (ascend) */
void QuickSort (double *a, size_t *perm, uint64_t n) {
    if (n < 2) return;
    uint64_t i, j;
    double t, p;
    size_t tt;
    p = a[n / 2];
    for (i = 0, j = n - 1;; i++, j--) {
        while (a[i] < p) i++;
        while (p < a[j]) j--;
        if (i >= j) break;
        t = a[i];
        a[i] = a[j];
        a[j] = t;
        tt = perm[i];
        perm[i] = perm[j];
        perm[j] = tt;
    }
    QuickSort(a, perm, i);
    QuickSort(a + i, perm + i, n - i);
}

/* compute PM */
double getPM_HD(double *absL, uint64_t n) {
    double pm = 0;
    for(size_t i = 0; i < n; i++){
        pm += JacLog( -absL[i] );
    }
    return pm;
}

/* compute PM based on PM_HD */
double getPM(uint8_t *TEP, double *absL, double PM_HD, uint64_t n) {
    double pm = PM_HD;
    for(size_t i = 0; i < n; i++){
        if (TEP[i] == 1)
            pm += (JacLog(absL[i]) - JacLog( -absL[i] ));
    }
    return pm;
}

/* compute List Confidence */
double getLConf(double *pNL, double P_notGuess, uint64_t cur_L, double *score, uint64_t s, uint8_t even) {
    double P_positive = 0.0;
    for(size_t i = 0; i <= cur_L; i++)
        P_positive += exp(-score[4*i+1]);
    if(even==1) {
        s--;
    }
    double P_negative = pow(2.0, -(double)s) * P_notGuess;
    pNL[0] = P_negative;
    return (P_positive / (P_positive + P_negative));
}

/* Mountain Build */
void mountain_build(int32_t *u, int32_t k, int32_t w, int32_t W1, int32_t n1){
    for(size_t i = k + 1; i < w; i++)
        u[i] = u[k];
    uint64_t W2 = W1;
    for(size_t i = 0; i < w; i++)
        W2 -= u[i];
    uint64_t q = (uint64_t)floor( W2 / (n1 - u[k]) );
    uint64_t r = W2 - q*(n1 - u[k]);
    if (q != 0){
        for(size_t i = w-q; i < w; i++)
            u[i] = n1;
    }
    if (w > q)
        u[w-q-1] = u[w-q-1] + r;
}

/* Recover probability that each decoding is correct*/
void getAPP(uint64_t cur_L, double *score, double *APP) {
    
    double P_positive = 0.0;
    for(size_t i = 0; i < cur_L; i++)
        P_positive += exp(-score[4*i+1]);
    
    double denominator = score[4*cur_L-4+3]/P_positive;
    for(size_t i = 0; i < cur_L; i++)
        APP[i]=exp(-score[4*i+1])*denominator;    
}

/* Main Function */
void SOGRAND(double *chat, double *score, double *T, double *curL, double *pNL, double *APP, double *llr, uint8_t *H, uint64_t n, uint64_t s, int32_t IC, uint64_t L, uint64_t Tmax, double thres, uint8_t even){
    /* Create vectors */
    size_t *perm = calloc(n, sizeof(size_t));
    uint8_t *cHD = calloc(n, sizeof(uint8_t));
    uint8_t parity_cHD;
    uint8_t *TEP = calloc(n, sizeof(uint8_t));
    uint8_t *c   = calloc(n, sizeof(uint8_t));
    double *absL = calloc(n, sizeof(double));
    int32_t *u   = calloc(n, sizeof(int32_t));
    int32_t *d   = calloc(n, sizeof(int32_t));
    int32_t *D   = calloc(n, sizeof(int32_t));
    for(size_t i = 0; i < n; i++)
        perm[i] = i;
    for(size_t i = 0; i < 4*L; i++)
        score[i] = 0;
    for(size_t i = 0; i < L; i++)
        APP[i] = 0;
    uint64_t cur_L = 0;
    /* Initialize */
    HardDec(cHD, llr, n);
    parity_cHD = parity(cHD,n);
    pNL[0] = 0.0;
    

    if (Tmax==0){
        Tmax=Inf;
    }    

    for (size_t i = 0; i < n; i++){
        TEP[i] = 0;
        absL[i] = fabs(llr[i]);
    }

    double P_notGuess = 1.0; // If an odd code
    if (even==1) {           // If an even code, there's a minor change
        P_notGuess = prob_parity(parity_cHD,absL,n);
    }
    double current_score=0;

    double PM_HD = getPM_HD(absL, n);
    QuickSort(absL, perm, n);
    /* Intercept */
    if (IC < 0){
        double beta = ( absL[ (uint64_t)round((double)n/2) - 1] - absL[0] ) / ( round((double)n/2) - 1 ); // slope
        IC = findMax( (int32_t)round( absL[0]/beta -1 ), 0 );
    }
    /* Main GRAND */
    AddTEP(c, cHD, TEP, perm, n);
    T[0] = 1;
    if (parity_cHD==0||even==0)  
        P_notGuess -= exp(-getPM(TEP, absL, PM_HD, n)); // Only remove the query's likelihood if parity is 0
    if (ParityCheck(c, H, n, s) == 1){
        score[0] = 0;
        score[1] = getPM(TEP, absL, PM_HD, n);
        for(size_t i = 0; i < n; i++)
            chat[i] = c[i];
        score[2] = 1;
        score[3] = getLConf(pNL, P_notGuess, cur_L, score, s, even);
        cur_L++;
        curL[0] = (double)cur_L;
        if ( (score[3] > thres) || (cur_L == L) ){
            getAPP(cur_L, score, APP);
            free(perm);
            free(cHD);
            free(TEP);
            free(c);
            free(absL);
            free(u);
            free(D);
            free(d);
            return;
        }
    }
    /* main loop */
    int32_t w = 0;
    int parity_w;
    int32_t W = 0;
    int32_t wt = IC + 1;
    int32_t wt_max = IC*n+n*(n+1)/2;
    int32_t W1;
    int32_t n1;
    int32_t k;
    double temp = 1+2*( (double)n + (double)IC );
    while ( (cur_L < L) && (T[0] < Tmax) && (wt <= wt_max) ){
        w = findMax(1, (int32_t)ceil((temp-sqrt(pow((double)temp, 2.0)-8*wt))/2) );
        parity_w = w%2;
        if (even==1) {
          if (parity_w!=parity_cHD){
            w ++;
          }
        }
        while ( w <= n ){
            W = wt - IC*w;
            if(W < w*(w+1)/2)
                break;
            /* landslide */
            W1 = W - w*(w+1)/2;
            n1 = n - w;
            for (size_t i = 0; i < w; i++)
                u[i] = 0;
            mountain_build(u,0,w,W1,n1);
            for (size_t i = 0; i < n; i++)
                TEP[i] = 0;
            for (size_t i = 0; i < w; i++)
                TEP[i+u[i]] = 1;
            AddTEP(c, cHD, TEP, perm, n);
            T[0]++;
            P_notGuess -= exp(-getPM(TEP, absL, PM_HD, n));

            if (ParityCheck(c, H, n, s) == 1){
                score[4*cur_L] = wt;
                score[4*cur_L+1] = getPM(TEP, absL, PM_HD, n);
                for(size_t i = 0; i < n; i++)
                    chat[cur_L*n + i] = c[i];
                score[4*cur_L+2] = T[0];
                score[4*cur_L+3] = getLConf(pNL, P_notGuess, cur_L, score, s, even);
                cur_L++;
                curL[0] = (double)cur_L;
                if ( (score[4*(cur_L-1)+3] > thres) || (cur_L == L)){
                    getAPP(cur_L, score, APP);
                    free(perm);
                    free(cHD);
                    free(TEP);
                    free(c);
                    free(absL);
                    free(u);
                    free(D);
                    free(d);
                    return;
                }
            }
            /* end */
            for (size_t i = 0; i < w - 1; i++)
                d[i] = u[i+1] - u[i];
            d[w-1] = 0;
            D[w-1] = d[w-1];
            for (size_t i = 1; i < w; i++)
                D[w-i-1] = D[w-i] + d[w-i-1];
            while( D[0] >= 2 ){
                k = 0;
                for (size_t i = w-1; i > 0; i--){
                    if (D[i] >= 2){
                        k = i;
                        break;
                    }
                }
                u[k] ++;
                mountain_build(u,k,w,W1,n1);
                for (size_t i = 0; i < n; i++)
                    TEP[i] = 0;
                for (size_t i = 0; i < w; i++)
                    TEP[i+u[i]] = 1;
                AddTEP(c, cHD, TEP, perm, n);
                T[0]++;
                P_notGuess -= exp(-getPM(TEP, absL, PM_HD, n));
                if (ParityCheck(c, H, n, s) == 1){
                    score[4*cur_L] = wt;
                    score[4*cur_L+1] = getPM(TEP, absL, PM_HD, n);
                    for(size_t i = 0; i < n; i++)
                        chat[cur_L*n + i] = c[i];
                    score[4*cur_L+2] = T[0];
                    score[4*cur_L+3] = getLConf(pNL, P_notGuess, cur_L, score, s, even);
                    cur_L++;
                    curL[0] = (double)cur_L;
                    if ( (score[4*(cur_L-1)+3] > thres) || (cur_L == L) ){
                        getAPP(cur_L, score, APP);
                        free(perm);
                        free(cHD);
                        free(TEP);
                        free(c);
                        free(absL);
                        free(u);
                        free(D);
                        free(d);
                        return;
                    }
                }
                for (size_t i = 0; i < w - 1; i++)
                    d[i] = u[i+1] - u[i];
                d[w-1] = 0;
                D[w-1] = d[w-1];
                for (size_t i = 1; i < w; i++)
                    D[w-i-1] = D[w-i] + d[w-i-1];
            }
            w ++;
            parity_w = w%2;
            if (even==1) {
                if (parity_w!=parity_cHD){
                w ++;
                }
            }
        }
        wt ++;
    }
    /* Clean up allocated memory */
    free(perm);
    free(cHD);
    free(TEP);
    free(c);
    free(absL);
    free(u);
    free(D);
    free(d);
}

/* Mexfunction Interface */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    /* Check for proper number of arguments */
    if( nrhs != 7 )
        mexErrMsgTxt("Must have 7 input arguments.");
    if( nlhs != 6 )
        mexErrMsgTxt("Must have 6 output arguments.");
    if( mxGetN(prhs[0]) != 1 || !mxIsClass(prhs[0], "double"))
        mexErrMsgTxt("First Input (LLR) must be a column-vector of type double.");
    uint64_t n = mxGetNumberOfElements(prhs[0]); // code length n
    if( mxGetN(prhs[1]) != 1 || !mxIsClass(prhs[1], "uint8"))
        mexErrMsgTxt("Second Input (Parity check matrix) must be a column-vector of length s*n and type uint8.");
    uint64_t N = mxGetNumberOfElements(prhs[1]); // s*n
    uint64_t s = (uint64_t) N/n; // number of constraints
    if( mxGetNumberOfElements(prhs[2]) != 1 || !mxIsClass(prhs[2], "int32"))
        mexErrMsgTxt("Third Input (Intercept) must be an int32 scalar.");
    if( mxGetNumberOfElements(prhs[3]) != 1 || !mxIsClass(prhs[3], "uint64"))
        mexErrMsgTxt("Third Input (List size) must be an uint64 scalar.");
    if( mxGetNumberOfElements(prhs[4]) != 1 || !mxIsClass(prhs[4], "uint64"))
        mexErrMsgTxt("Forth Input (Maximum number of guesses) must be an uint64 scalar.");
    if( mxGetNumberOfElements(prhs[5]) != 1 || !mxIsClass(prhs[5], "double"))
        mexErrMsgTxt("Fifth Input (Required listConf) must be an double scalar.");
    /* input */
    double *llr = mxGetPr(prhs[0]);
    uint8_t *H = (uint8_t *)mxGetData(prhs[1]);
    int32_t IC = (int32_t)mxGetScalar(prhs[2]);
    uint64_t L = (uint64_t)mxGetScalar(prhs[3]);
    uint64_t Tmax = (uint64_t)mxGetScalar(prhs[4]);
    double thres = (double)mxGetScalar(prhs[5]);
    uint8_t even = (uint64_t)mxGetScalar(prhs[6]);
    /* create the output vector */
    plhs[0] = mxCreateDoubleMatrix(n, L, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(4, L, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[5] = mxCreateDoubleMatrix(1, L, mxREAL);
    double *chat  = mxGetData(plhs[0]);
    double *score = mxGetData(plhs[1]);
    double *T     = mxGetData(plhs[2]);
    double *curL  = mxGetData(plhs[3]);
    double *pNL    = mxGetData(plhs[4]);
    double *APP    = mxGetData(plhs[5]);
    /* use C functions in mexfunction */
    SOGRAND(chat, score, T, curL, pNL, APP, llr, H, n, s, IC, L, Tmax, thres, even);
}
