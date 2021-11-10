#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

#if 0
/* Elements
0  H
1  C
2  O
3  N
*/
#endif

/* Species */
#define H_ID 0
#define H2_ID 1
#define O_ID 2
#define O2_ID 3
#define OH_ID 4
#define H2O_ID 5
#define CO_ID 6
#define CO2_ID 7
#define CH3_ID 8
#define CH4_ID 9
#define HO2_ID 10
#define H2O2_ID 11
#define CH2O_ID 12
#define C2H6_ID 13
#define C2H4_ID 14
#define C2H5_ID 15
#define C2H_ID 16
#define C2H2_ID 17
#define CH3OH_ID 18
#define CH2CO_ID 19
#define HCCO_ID 20
#define CH3CHO_ID 21
#define C3H4A_ID 22
#define C3H4P_ID 23
#define C3H6_ID 24
#define C4H6_ID 25
#define C4H7_ID 26
#define C4H81_ID 27
#define PC4H9_ID 28
#define CH3COCH2_ID 29
#define C2H5CHO_ID 30
#define C5H9_ID 31
#define C5H101_ID 32
#define C5H111_ID 33
#define CH3O2_ID 34
#define CH3O2H_ID 35
#define C2H3CO_ID 36
#define C2H3CHO_ID 37
#define C3H5A_ID 38
#define C3H3_ID 39
#define NC3H7CHO_ID 40
#define C2H5COCH2_ID 41
#define CH3CHCO_ID 42
#define NC3H7COCH3_ID 43
#define NC3H7COCH2_ID 44
#define NC7H16_ID 45
#define C7H151_ID 46
#define C7H15O2_ID 47
#define C7H14OOHO2_ID 48
#define C7H14O_ID 49
#define NC7KET_ID 50
#define N2_ID 51

#define NUM_ELEMENTS 4
#define NUM_SPECIES 52
#define NUM_REACTIONS 0

#define NUM_FIT 4


/* ALWAYS on CPU stuff -- can have different def depending on if we are CPU or GPU based. Defined in mechanism.cpp  */
void atomicWeight(amrex::Real *  awt);
/* MISC  */
void CKAWT(amrex::Real *  awt);
void CKNCF(int * ncf);
void CKSYME_STR(amrex::Vector<std::string>& ename);
void CKSYMS_STR(amrex::Vector<std::string>& kname);
/* SPARSE INFORMATION  */
void SPARSITY_INFO(int * nJdata, int * consP, int NCELLS);
void SPARSITY_INFO_SYST(int * nJdata, int * consP, int NCELLS);
void SPARSITY_INFO_SYST_SIMPLIFIED(int * nJdata, int * consP);
void SPARSITY_PREPROC_CSC(int * rowVals, int * colPtrs, int * consP, int NCELLS);
void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, int * consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtrs, int * consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, int * consP);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, int * consP, int base);



/*A few mechanism parameters */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKINDX(int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 4;
    *kk = 52;
    *ii = 0;
    *nfit = -1; /*Why do you need this anyway ?  */
}


/* inverse molecular weights  */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void get_imw(amrex::Real *imw_new){
    imw_new[0] = 1.0/1.007970;/*H */
    imw_new[1] = 1.0/2.015940;/*H2 */
    imw_new[2] = 1.0/15.999400;/*O */
    imw_new[3] = 1.0/31.998800;/*O2 */
    imw_new[4] = 1.0/17.007370;/*OH */
    imw_new[5] = 1.0/18.015340;/*H2O */
    imw_new[6] = 1.0/28.010550;/*CO */
    imw_new[7] = 1.0/44.009950;/*CO2 */
    imw_new[8] = 1.0/15.035060;/*CH3 */
    imw_new[9] = 1.0/16.043030;/*CH4 */
    imw_new[10] = 1.0/33.006770;/*HO2 */
    imw_new[11] = 1.0/34.014740;/*H2O2 */
    imw_new[12] = 1.0/30.026490;/*CH2O */
    imw_new[13] = 1.0/30.070120;/*C2H6 */
    imw_new[14] = 1.0/28.054180;/*C2H4 */
    imw_new[15] = 1.0/29.062150;/*C2H5 */
    imw_new[16] = 1.0/25.030270;/*C2H */
    imw_new[17] = 1.0/26.038240;/*C2H2 */
    imw_new[18] = 1.0/32.042430;/*CH3OH */
    imw_new[19] = 1.0/42.037640;/*CH2CO */
    imw_new[20] = 1.0/41.029670;/*HCCO */
    imw_new[21] = 1.0/44.053580;/*CH3CHO */
    imw_new[22] = 1.0/40.065330;/*C3H4-A */
    imw_new[23] = 1.0/40.065330;/*C3H4-P */
    imw_new[24] = 1.0/42.081270;/*C3H6 */
    imw_new[25] = 1.0/54.092420;/*C4H6 */
    imw_new[26] = 1.0/55.100390;/*C4H7 */
    imw_new[27] = 1.0/56.108360;/*C4H8-1 */
    imw_new[28] = 1.0/57.116330;/*PC4H9 */
    imw_new[29] = 1.0/57.072700;/*CH3COCH2 */
    imw_new[30] = 1.0/58.080670;/*C2H5CHO */
    imw_new[31] = 1.0/69.127480;/*C5H9 */
    imw_new[32] = 1.0/70.135450;/*C5H10-1 */
    imw_new[33] = 1.0/71.143420;/*C5H11-1 */
    imw_new[34] = 1.0/47.033860;/*CH3O2 */
    imw_new[35] = 1.0/48.041830;/*CH3O2H */
    imw_new[36] = 1.0/55.056760;/*C2H3CO */
    imw_new[37] = 1.0/56.064730;/*C2H3CHO */
    imw_new[38] = 1.0/41.073300;/*C3H5-A */
    imw_new[39] = 1.0/39.057360;/*C3H3 */
    imw_new[40] = 1.0/72.107760;/*NC3H7CHO */
    imw_new[41] = 1.0/71.099790;/*C2H5COCH2 */
    imw_new[42] = 1.0/56.064730;/*CH3CHCO */
    imw_new[43] = 1.0/86.134850;/*NC3H7COCH3 */
    imw_new[44] = 1.0/85.126880;/*NC3H7COCH2 */
    imw_new[45] = 1.0/100.205570;/*NC7H16 */
    imw_new[46] = 1.0/99.197600;/*C7H15-1 */
    imw_new[47] = 1.0/131.196400;/*C7H15O2 */
    imw_new[48] = 1.0/163.195200;/*C7H14OOHO2 */
    imw_new[49] = 1.0/114.189030;/*C7H14O */
    imw_new[50] = 1.0/146.187830;/*NC7KET */
    imw_new[51] = 1.0/28.013400;/*N2 */
}

/* molecular weights  */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void get_mw(amrex::Real *mw_new){
    mw_new[0] = 1.007970;/*H */
    mw_new[1] = 2.015940;/*H2 */
    mw_new[2] = 15.999400;/*O */
    mw_new[3] = 31.998800;/*O2 */
    mw_new[4] = 17.007370;/*OH */
    mw_new[5] = 18.015340;/*H2O */
    mw_new[6] = 28.010550;/*CO */
    mw_new[7] = 44.009950;/*CO2 */
    mw_new[8] = 15.035060;/*CH3 */
    mw_new[9] = 16.043030;/*CH4 */
    mw_new[10] = 33.006770;/*HO2 */
    mw_new[11] = 34.014740;/*H2O2 */
    mw_new[12] = 30.026490;/*CH2O */
    mw_new[13] = 30.070120;/*C2H6 */
    mw_new[14] = 28.054180;/*C2H4 */
    mw_new[15] = 29.062150;/*C2H5 */
    mw_new[16] = 25.030270;/*C2H */
    mw_new[17] = 26.038240;/*C2H2 */
    mw_new[18] = 32.042430;/*CH3OH */
    mw_new[19] = 42.037640;/*CH2CO */
    mw_new[20] = 41.029670;/*HCCO */
    mw_new[21] = 44.053580;/*CH3CHO */
    mw_new[22] = 40.065330;/*C3H4-A */
    mw_new[23] = 40.065330;/*C3H4-P */
    mw_new[24] = 42.081270;/*C3H6 */
    mw_new[25] = 54.092420;/*C4H6 */
    mw_new[26] = 55.100390;/*C4H7 */
    mw_new[27] = 56.108360;/*C4H8-1 */
    mw_new[28] = 57.116330;/*PC4H9 */
    mw_new[29] = 57.072700;/*CH3COCH2 */
    mw_new[30] = 58.080670;/*C2H5CHO */
    mw_new[31] = 69.127480;/*C5H9 */
    mw_new[32] = 70.135450;/*C5H10-1 */
    mw_new[33] = 71.143420;/*C5H11-1 */
    mw_new[34] = 47.033860;/*CH3O2 */
    mw_new[35] = 48.041830;/*CH3O2H */
    mw_new[36] = 55.056760;/*C2H3CO */
    mw_new[37] = 56.064730;/*C2H3CHO */
    mw_new[38] = 41.073300;/*C3H5-A */
    mw_new[39] = 39.057360;/*C3H3 */
    mw_new[40] = 72.107760;/*NC3H7CHO */
    mw_new[41] = 71.099790;/*C2H5COCH2 */
    mw_new[42] = 56.064730;/*CH3CHCO */
    mw_new[43] = 86.134850;/*NC3H7COCH3 */
    mw_new[44] = 85.126880;/*NC3H7COCH2 */
    mw_new[45] = 100.205570;/*NC7H16 */
    mw_new[46] = 99.197600;/*C7H15-1 */
    mw_new[47] = 131.196400;/*C7H15O2 */
    mw_new[48] = 163.195200;/*C7H14OOHO2 */
    mw_new[49] = 114.189030;/*C7H14O */
    mw_new[50] = 146.187830;/*NC7KET */
    mw_new[51] = 28.013400;/*N2 */
}



/* Returns R, Rc, Patm */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRP(amrex::Real *  ru, amrex::Real *  ruc, amrex::Real *  pa)
{
     *ru  = 8.31446261815324e+07; 
     *ruc = 1.98721558317399615845; 
     *pa  = 1.01325e+06; 
}


/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void cv_R(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +1.50104422e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            +2.29812400e+00
            +8.24944200e-04 * tc[1]
            -8.14301500e-07 * tc[2]
            -9.47543400e-11 * tc[3]
            +4.13487200e-13 * tc[4];
        /*species 2: O */
        species[2] =
            +1.94642900e+00
            -1.63816600e-03 * tc[1]
            +2.42103200e-06 * tc[2]
            -1.60284300e-09 * tc[3]
            +3.89069600e-13 * tc[4];
        /*species 3: O2 */
        species[3] =
            +2.21293600e+00
            +1.12748600e-03 * tc[1]
            -5.75615000e-07 * tc[2]
            +1.31387700e-09 * tc[3]
            -8.76855400e-13 * tc[4];
        /*species 5: H2O */
        species[5] =
            +2.38684200e+00
            +3.47498200e-03 * tc[1]
            -6.35469600e-06 * tc[2]
            +6.96858100e-09 * tc[3]
            -2.50658800e-12 * tc[4];
        /*species 6: CO */
        species[6] =
            +2.26245200e+00
            +1.51194100e-03 * tc[1]
            -3.88175500e-06 * tc[2]
            +5.58194400e-09 * tc[3]
            -2.47495100e-12 * tc[4];
        /*species 7: CO2 */
        species[7] =
            +1.27572500e+00
            +9.92207200e-03 * tc[1]
            -1.04091100e-05 * tc[2]
            +6.86668700e-09 * tc[3]
            -2.11728000e-12 * tc[4];
        /*species 8: CH3 */
        species[8] =
            +1.43044300e+00
            +1.11241000e-02 * tc[1]
            -1.68022000e-05 * tc[2]
            +1.62182900e-08 * tc[3]
            -5.86495300e-12 * tc[4];
        /*species 9: CH4 */
        species[9] =
            -2.21258500e-01
            +1.74766800e-02 * tc[1]
            -2.78340900e-05 * tc[2]
            +3.04970800e-08 * tc[3]
            -1.22393100e-11 * tc[4];
        /*species 11: H2O2 */
        species[11] =
            +2.38875400e+00
            +6.56922600e-03 * tc[1]
            -1.48501300e-07 * tc[2]
            -4.62580600e-09 * tc[3]
            +2.47151500e-12 * tc[4];
        /*species 12: CH2O */
        species[12] =
            +6.52731000e-01
            +1.26314400e-02 * tc[1]
            -1.88816800e-05 * tc[2]
            +2.05003100e-08 * tc[3]
            -8.41323700e-12 * tc[4];
        /*species 14: C2H4 */
        species[14] =
            -1.86148800e+00
            +2.79616300e-02 * tc[1]
            -3.38867700e-05 * tc[2]
            +2.78515200e-08 * tc[3]
            -9.73787900e-12 * tc[4];
        /*species 15: C2H5 */
        species[15] =
            +1.69070200e+00
            +8.71913300e-03 * tc[1]
            +4.41983900e-06 * tc[2]
            +9.33870300e-10 * tc[3]
            -3.92777300e-12 * tc[4];
        /*species 16: C2H */
        species[16] =
            +1.73770400e+00
            +8.04844600e-03 * tc[1]
            -9.24431000e-06 * tc[2]
            +6.52525900e-09 * tc[3]
            -1.93958000e-12 * tc[4];
        /*species 17: C2H2 */
        species[17] =
            +1.01356200e+00
            +1.51904500e-02 * tc[1]
            -1.61631900e-05 * tc[2]
            +9.07899200e-09 * tc[3]
            -1.91274600e-12 * tc[4];
        /*species 18: CH3OH */
        species[18] =
            +1.66011500e+00
            +7.34150800e-03 * tc[1]
            +7.17005100e-06 * tc[2]
            -8.79319400e-09 * tc[3]
            +2.39057000e-12 * tc[4];
        /*species 19: CH2CO */
        species[19] =
            +1.97497100e+00
            +1.21187100e-02 * tc[1]
            -2.34504600e-06 * tc[2]
            -6.46668500e-09 * tc[3]
            +3.90564900e-12 * tc[4];
        /*species 20: HCCO */
        species[20] =
            +4.04796500e+00
            +4.45347800e-03 * tc[1]
            +2.26828300e-07 * tc[2]
            -1.48209500e-09 * tc[3]
            +2.25074200e-13 * tc[4];
        /*species 21: CH3CHO */
        species[21] =
            +1.50569500e+00
            +1.33699100e-02 * tc[1]
            +4.67195300e-06 * tc[2]
            -1.12814000e-08 * tc[3]
            +4.26356600e-12 * tc[4];
        /*species 39: C3H3 */
        species[39] =
            +3.75420000e+00
            +1.10802800e-02 * tc[1]
            +2.79332300e-07 * tc[2]
            -5.47921200e-09 * tc[3]
            +1.94962900e-12 * tc[4];
        /*species 51: N2 */
        species[51] =
            +2.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +1.50104422e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            +1.99142300e+00
            +7.00064400e-04 * tc[1]
            -5.63382900e-08 * tc[2]
            -9.23157800e-12 * tc[3]
            +1.58275200e-15 * tc[4];
        /*species 2: O */
        species[2] =
            +1.54206000e+00
            -2.75506200e-05 * tc[1]
            -3.10280300e-09 * tc[2]
            +4.55106700e-12 * tc[3]
            -4.36805200e-16 * tc[4];
        /*species 3: O2 */
        species[3] =
            +2.69757800e+00
            +6.13519700e-04 * tc[1]
            -1.25884200e-07 * tc[2]
            +1.77528100e-11 * tc[3]
            -1.13643500e-15 * tc[4];
        /*species 5: H2O */
        species[5] =
            +1.67214600e+00
            +3.05629300e-03 * tc[1]
            -8.73026000e-07 * tc[2]
            +1.20099600e-10 * tc[3]
            -6.39161800e-15 * tc[4];
        /*species 6: CO */
        species[6] =
            +2.02507800e+00
            +1.44268900e-03 * tc[1]
            -5.63082800e-07 * tc[2]
            +1.01858100e-10 * tc[3]
            -6.91095200e-15 * tc[4];
        /*species 7: CO2 */
        species[7] =
            +3.45362300e+00
            +3.14016900e-03 * tc[1]
            -1.27841100e-06 * tc[2]
            +2.39399700e-10 * tc[3]
            -1.66903300e-14 * tc[4];
        /*species 8: CH3 */
        species[8] =
            +1.84405200e+00
            +6.13797400e-03 * tc[1]
            -2.23034500e-06 * tc[2]
            +3.78516100e-10 * tc[3]
            -2.45215900e-14 * tc[4];
        /*species 9: CH4 */
        species[9] =
            +6.83479000e-01
            +1.02372400e-02 * tc[1]
            -3.87512900e-06 * tc[2]
            +6.78558500e-10 * tc[3]
            -4.50342300e-14 * tc[4];
        /*species 11: H2O2 */
        species[11] =
            +3.57316700e+00
            +4.33613600e-03 * tc[1]
            -1.47468900e-06 * tc[2]
            +2.34890400e-10 * tc[3]
            -1.43165400e-14 * tc[4];
        /*species 12: CH2O */
        species[12] =
            +1.99560600e+00
            +6.68132100e-03 * tc[1]
            -2.62895500e-06 * tc[2]
            +4.73715300e-10 * tc[3]
            -3.21251700e-14 * tc[4];
        /*species 14: C2H4 */
        species[14] =
            +2.52841900e+00
            +1.14851800e-02 * tc[1]
            -4.41838500e-06 * tc[2]
            +7.84460100e-10 * tc[3]
            -5.26684800e-14 * tc[4];
        /*species 15: C2H5 */
        species[15] =
            +6.19048000e+00
            +6.48407700e-03 * tc[1]
            -6.42806500e-07 * tc[2]
            -2.34787900e-10 * tc[3]
            +3.88087700e-14 * tc[4];
        /*species 16: C2H */
        species[16] =
            +2.98636700e+00
            +3.14312300e-03 * tc[1]
            -1.26724300e-06 * tc[2]
            +2.92436300e-10 * tc[3]
            -2.71632000e-14 * tc[4];
        /*species 17: C2H2 */
        species[17] =
            +3.43677000e+00
            +5.37603900e-03 * tc[1]
            -1.91281700e-06 * tc[2]
            +3.28637900e-10 * tc[3]
            -2.15671000e-14 * tc[4];
        /*species 18: CH3OH */
        species[18] =
            +3.02906100e+00
            +9.37659300e-03 * tc[1]
            -3.05025400e-06 * tc[2]
            +4.35879300e-10 * tc[3]
            -2.22472300e-14 * tc[4];
        /*species 19: CH2CO */
        species[19] =
            +5.03881700e+00
            +5.80484000e-03 * tc[1]
            -1.92095400e-06 * tc[2]
            +2.79448500e-10 * tc[3]
            -1.45886800e-14 * tc[4];
        /*species 20: HCCO */
        species[20] =
            +5.75807300e+00
            +2.00040000e-03 * tc[1]
            -2.02760700e-07 * tc[2]
            -1.04113200e-10 * tc[3]
            +1.96516500e-14 * tc[4];
        /*species 21: CH3CHO */
        species[21] =
            +4.86865000e+00
            +1.07942400e-02 * tc[1]
            -3.64553000e-06 * tc[2]
            +5.41291200e-10 * tc[3]
            -2.89684400e-14 * tc[4];
        /*species 39: C3H3 */
        species[39] =
            +7.83104700e+00
            +4.35719500e-03 * tc[1]
            -4.10906700e-07 * tc[2]
            -2.36872300e-10 * tc[3]
            +4.37652000e-14 * tc[4];
        /*species 51: N2 */
        species[51] =
            +1.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }

    /*species with midpoint at T=1357 kelvin */
    if (T < 1357) {
        /*species 4: OH */
        species[4] =
            +2.43586219e+00
            +2.02235804e-04 * tc[1]
            -1.13546412e-07 * tc[2]
            +2.42445149e-10 * tc[3]
            -7.43651031e-14 * tc[4];
    } else {
        /*species 4: OH */
        species[4] =
            +1.62599754e+00
            +1.31992406e-03 * tc[1]
            -3.59724670e-07 * tc[2]
            +4.25630800e-11 * tc[3]
            -1.82048016e-15 * tc[4];
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 10: HO2 */
        species[10] =
            +2.18310656e+00
            +3.66767950e-03 * tc[1]
            -9.32385122e-07 * tc[2]
            -3.25852919e-10 * tc[3]
            +1.51139912e-13 * tc[4];
        /*species 35: CH3O2H */
        species[35] =
            +2.23442817e+00
            +1.90129767e-02 * tc[1]
            -1.13386287e-05 * tc[2]
            +3.40306653e-09 * tc[3]
            -4.11830222e-13 * tc[4];
        /*species 46: C7H15-1 */
        species[46] =
            -1.49957041e+00
            +8.08826467e-02 * tc[1]
            -5.00532754e-05 * tc[2]
            +1.56549308e-08 * tc[3]
            -1.96616227e-12 * tc[4];
        /*species 47: C7H15O2 */
        species[47] =
            -6.44747083e-01
            +9.42381007e-02 * tc[1]
            -6.66955106e-05 * tc[2]
            +2.54795583e-08 * tc[3]
            -4.13211105e-12 * tc[4];
        /*species 48: C7H14OOHO2 */
        species[48] =
            +4.09124790e-01
            +1.10720834e-01 * tc[1]
            -8.71178075e-05 * tc[2]
            +3.63662749e-08 * tc[3]
            -6.28350689e-12 * tc[4];
    } else {
        /*species 10: HO2 */
        species[10] =
            +3.10547423e+00
            +2.38452835e-03 * tc[1]
            -8.06347989e-07 * tc[2]
            +1.24191723e-10 * tc[3]
            -7.16400108e-15 * tc[4];
        /*species 35: CH3O2H */
        species[35] =
            +7.43117091e+00
            +8.06817909e-03 * tc[1]
            -2.77094921e-06 * tc[2]
            +4.31332243e-10 * tc[3]
            -2.50692146e-14 * tc[4];
        /*species 46: C7H15-1 */
        species[46] =
            +2.07940709e+01
            +3.26280243e-02 * tc[1]
            -1.11138244e-05 * tc[2]
            +1.72067148e-09 * tc[3]
            -9.96366999e-14 * tc[4];
        /*species 47: C7H15O2 */
        species[47] =
            +2.50535640e+01
            +3.43831688e-02 * tc[1]
            -1.18713247e-05 * tc[2]
            +1.85450241e-09 * tc[3]
            -1.08052399e-13 * tc[4];
        /*species 48: C7H14OOHO2 */
        species[48] =
            +3.18358296e+01
            +3.32767285e-02 * tc[1]
            -1.15833285e-05 * tc[2]
            +1.81958638e-09 * tc[3]
            -1.06434419e-13 * tc[4];
    }

    /*species with midpoint at T=1384 kelvin */
    if (T < 1384) {
        /*species 13: C2H6 */
        species[13] =
            -1.02528543e+00
            +2.40764754e-02 * tc[1]
            -1.11893472e-05 * tc[2]
            +2.08340901e-09 * tc[3]
            -5.29868616e-14 * tc[4];
    } else {
        /*species 13: C2H6 */
        species[13] =
            +5.10683385e+00
            +1.29236361e-02 * tc[1]
            -4.42527196e-06 * tc[2]
            +6.87391726e-10 * tc[3]
            -3.98901732e-14 * tc[4];
    }

    /*species with midpoint at T=1400 kelvin */
    if (T < 1400) {
        /*species 22: C3H4-A */
        species[22] =
            +1.53983100e+00
            +1.63343700e-02 * tc[1]
            -1.76495000e-06 * tc[2]
            -4.64736500e-09 * tc[3]
            +1.72913100e-12 * tc[4];
        /*species 23: C3H4-P */
        species[23] =
            +2.02973000e+00
            +1.49896100e-02 * tc[1]
            -1.39850000e-06 * tc[2]
            -3.96961900e-09 * tc[3]
            +1.38821700e-12 * tc[4];
        /*species 42: CH3CHCO */
        species[42] =
            +4.83801190e-01
            +3.22203013e-02 * tc[1]
            -2.70250033e-05 * tc[2]
            +1.20499164e-08 * tc[3]
            -2.18365931e-12 * tc[4];
    } else {
        /*species 22: C3H4-A */
        species[22] =
            +8.77625600e+00
            +5.30213800e-03 * tc[1]
            -3.70111800e-07 * tc[2]
            -3.02638600e-10 * tc[3]
            +5.08958100e-14 * tc[4];
        /*species 23: C3H4-P */
        species[23] =
            +8.76810200e+00
            +5.21915100e-03 * tc[1]
            -3.75314000e-07 * tc[2]
            -2.99219100e-10 * tc[3]
            +5.10787800e-14 * tc[4];
        /*species 42: CH3CHCO */
        species[42] =
            +9.02191230e+00
            +9.56966300e-03 * tc[1]
            -3.26221644e-06 * tc[2]
            +5.05231706e-10 * tc[3]
            -2.92593257e-14 * tc[4];
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 24: C3H6 */
        species[24] =
            -6.05384556e-01
            +2.89107662e-02 * tc[1]
            -1.54886808e-05 * tc[2]
            +3.88814209e-09 * tc[3]
            -3.37890352e-13 * tc[4];
    } else {
        /*species 24: C3H6 */
        species[24] =
            +7.01595958e+00
            +1.37023634e-02 * tc[1]
            -4.66249733e-06 * tc[2]
            +7.21254402e-10 * tc[3]
            -4.17370126e-14 * tc[4];
    }

    /*species with midpoint at T=1398 kelvin */
    if (T < 1398) {
        /*species 25: C4H6 */
        species[25] =
            -2.43095121e+00
            +4.78706062e-02 * tc[1]
            -4.15446800e-05 * tc[2]
            +1.91549552e-08 * tc[3]
            -3.57158507e-12 * tc[4];
    } else {
        /*species 25: C4H6 */
        species[25] =
            +1.01633789e+01
            +1.37163965e-02 * tc[1]
            -4.69715783e-06 * tc[2]
            +7.29693836e-10 * tc[3]
            -4.23486203e-14 * tc[4];
    }

    /*species with midpoint at T=1392 kelvin */
    if (T < 1392) {
        /*species 26: C4H7 */
        species[26] =
            -1.35050835e+00
            +4.26511243e-02 * tc[1]
            -2.90979373e-05 * tc[2]
            +1.05403914e-08 * tc[3]
            -1.60059854e-12 * tc[4];
        /*species 27: C4H8-1 */
        species[27] =
            -1.83137209e+00
            +4.52580978e-02 * tc[1]
            -2.93658559e-05 * tc[2]
            +1.00220436e-08 * tc[3]
            -1.43191680e-12 * tc[4];
        /*species 31: C5H9 */
        species[31] =
            -2.38013950e+00
            +5.57608487e-02 * tc[1]
            -3.70143928e-05 * tc[2]
            +1.26883901e-08 * tc[3]
            -1.78538835e-12 * tc[4];
        /*species 32: C5H10-1 */
        species[32] =
            -2.06223481e+00
            +5.74218294e-02 * tc[1]
            -3.74486890e-05 * tc[2]
            +1.27364989e-08 * tc[3]
            -1.79609789e-12 * tc[4];
    } else {
        /*species 26: C4H7 */
        species[26] =
            +1.02103578e+01
            +1.60483196e-02 * tc[1]
            -5.46502292e-06 * tc[2]
            +8.45941053e-10 * tc[3]
            -4.89772739e-14 * tc[4];
        /*species 27: C4H8-1 */
        species[27] =
            +1.03508668e+01
            +1.80617877e-02 * tc[1]
            -6.16093029e-06 * tc[2]
            +9.54652959e-10 * tc[3]
            -5.53089641e-14 * tc[4];
        /*species 31: C5H9 */
        species[31] =
            +1.31860454e+01
            +2.07128899e-02 * tc[1]
            -7.06960617e-06 * tc[2]
            +1.09607133e-09 * tc[3]
            -6.35322208e-14 * tc[4];
        /*species 32: C5H10-1 */
        species[32] =
            +1.35851539e+01
            +2.24072471e-02 * tc[1]
            -7.63348025e-06 * tc[2]
            +1.18188966e-09 * tc[3]
            -6.84385139e-14 * tc[4];
    }

    /*species with midpoint at T=1395 kelvin */
    if (T < 1395) {
        /*species 28: PC4H9 */
        species[28] =
            -1.43777972e+00
            +4.78972364e-02 * tc[1]
            -3.14023159e-05 * tc[2]
            +1.09786472e-08 * tc[3]
            -1.62010664e-12 * tc[4];
    } else {
        /*species 28: PC4H9 */
        species[28] =
            +1.11510082e+01
            +1.94310717e-02 * tc[1]
            -6.61577950e-06 * tc[2]
            +1.02375136e-09 * tc[3]
            -5.92529707e-14 * tc[4];
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 29: CH3COCH2 */
        species[29] =
            +8.03391870e-01
            +3.01407085e-02 * tc[1]
            -1.93505552e-05 * tc[2]
            +6.38199034e-09 * tc[3]
            -8.66103180e-13 * tc[4];
        /*species 45: NC7H16 */
        species[45] =
            -2.26836187e+00
            +8.54355820e-02 * tc[1]
            -5.25346786e-05 * tc[2]
            +1.62945721e-08 * tc[3]
            -2.02394925e-12 * tc[4];
    } else {
        /*species 29: CH3COCH2 */
        species[29] =
            +9.23039750e+00
            +1.16494161e-02 * tc[1]
            -4.01005537e-06 * tc[2]
            +6.25205246e-10 * tc[3]
            -3.63784362e-14 * tc[4];
        /*species 45: NC7H16 */
        species[45] =
            +2.12148969e+01
            +3.47675750e-02 * tc[1]
            -1.18407129e-05 * tc[2]
            +1.83298478e-09 * tc[3]
            -1.06130266e-13 * tc[4];
    }

    /*species with midpoint at T=1378 kelvin */
    if (T < 1378) {
        /*species 30: C2H5CHO */
        species[30] =
            +1.16308444e+00
            +2.95501264e-02 * tc[1]
            -1.52446252e-05 * tc[2]
            +3.49503947e-09 * tc[3]
            -2.38896627e-13 * tc[4];
        /*species 40: NC3H7CHO */
        species[40] =
            +8.74159590e-01
            +4.19240315e-02 * tc[1]
            -2.35148779e-05 * tc[2]
            +6.26913673e-09 * tc[3]
            -6.09443908e-13 * tc[4];
    } else {
        /*species 30: C2H5CHO */
        species[30] =
            +9.24276950e+00
            +1.39641989e-02 * tc[1]
            -4.76248001e-06 * tc[2]
            +7.38105706e-10 * tc[3]
            -4.27759503e-14 * tc[4];
        /*species 40: NC3H7CHO */
        species[40] =
            +1.25988068e+01
            +1.81652474e-02 * tc[1]
            -6.17844458e-06 * tc[2]
            +9.55980208e-10 * tc[3]
            -5.53442958e-14 * tc[4];
    }

    /*species with midpoint at T=1397 kelvin */
    if (T < 1397) {
        /*species 33: C5H11-1 */
        species[33] =
            -1.90525591e+00
            +6.10632852e-02 * tc[1]
            -4.09491825e-05 * tc[2]
            +1.46093470e-08 * tc[3]
            -2.18859615e-12 * tc[4];
        /*species 38: C3H5-A */
        species[38] =
            -1.52913196e+00
            +3.34559100e-02 * tc[1]
            -2.53401027e-05 * tc[2]
            +1.02865754e-08 * tc[3]
            -1.73258340e-12 * tc[4];
    } else {
        /*species 33: C5H11-1 */
        species[33] =
            +1.43234740e+01
            +2.39041200e-02 * tc[1]
            -8.14771619e-06 * tc[2]
            +1.26176236e-09 * tc[3]
            -7.30677335e-14 * tc[4];
        /*species 38: C3H5-A */
        species[38] =
            +7.45883958e+00
            +1.12695483e-02 * tc[1]
            -3.83792864e-06 * tc[2]
            +5.94059119e-10 * tc[3]
            -3.43918030e-14 * tc[4];
    }

    /*species with midpoint at T=1385 kelvin */
    if (T < 1385) {
        /*species 34: CH3O2 */
        species[34] =
            +3.26146906e+00
            +1.00873599e-02 * tc[1]
            -3.21506184e-06 * tc[2]
            +2.09409267e-10 * tc[3]
            +4.18339103e-14 * tc[4];
    } else {
        /*species 34: CH3O2 */
        species[34] =
            +4.95787891e+00
            +7.90728626e-03 * tc[1]
            -2.68246234e-06 * tc[2]
            +4.13891337e-10 * tc[3]
            -2.39007330e-14 * tc[4];
    }

    /*species with midpoint at T=1402 kelvin */
    if (T < 1402) {
        /*species 36: C2H3CO */
        species[36] =
            +3.62420130e-01
            +3.15273972e-02 * tc[1]
            -3.00218935e-05 * tc[2]
            +1.48167112e-08 * tc[3]
            -2.87971530e-12 * tc[4];
    } else {
        /*species 36: C2H3CO */
        species[36] =
            +8.37467676e+00
            +7.91296900e-03 * tc[1]
            -2.67198280e-06 * tc[2]
            +4.11115430e-10 * tc[3]
            -2.36978981e-14 * tc[4];
    }

    /*species with midpoint at T=1393 kelvin */
    if (T < 1393) {
        /*species 37: C2H3CHO */
        species[37] =
            -7.07644838e-01
            +3.54321417e-02 * tc[1]
            -2.94936324e-05 * tc[2]
            +1.28100124e-08 * tc[3]
            -2.26144108e-12 * tc[4];
    } else {
        /*species 37: C2H3CHO */
        species[37] =
            +9.41849590e+00
            +9.48963321e-03 * tc[1]
            -3.29310529e-06 * tc[2]
            +5.16279203e-10 * tc[3]
            -3.01587291e-14 * tc[4];
    }

    /*species with midpoint at T=1383 kelvin */
    if (T < 1383) {
        /*species 41: C2H5COCH2 */
        species[41] =
            +5.40138560e-01
            +4.39486258e-02 * tc[1]
            -2.97002421e-05 * tc[2]
            +1.05495313e-08 * tc[3]
            -1.58598769e-12 * tc[4];
        /*species 43: NC3H7COCH3 */
        species[43] =
            +2.40373610e-01
            +5.42494047e-02 * tc[1]
            -3.14825332e-05 * tc[2]
            +9.14492321e-09 * tc[3]
            -1.07665566e-12 * tc[4];
    } else {
        /*species 41: C2H5COCH2 */
        species[41] =
            +1.32098738e+01
            +1.57866459e-02 * tc[1]
            -5.50529183e-06 * tc[2]
            +8.65870540e-10 * tc[3]
            -5.06913329e-14 * tc[4];
        /*species 43: NC3H7COCH3 */
        species[43] =
            +1.51841564e+01
            +2.35245881e-02 * tc[1]
            -8.09675055e-06 * tc[2]
            +1.26230690e-09 * tc[3]
            -7.34487839e-14 * tc[4];
    }

    /*species with midpoint at T=1389 kelvin */
    if (T < 1389) {
        /*species 44: NC3H7COCH2 */
        species[44] =
            -4.17007290e-02
            +5.68162532e-02 * tc[1]
            -3.99112781e-05 * tc[2]
            +1.52671514e-08 * tc[3]
            -2.49221047e-12 * tc[4];
    } else {
        /*species 44: NC3H7COCH2 */
        species[44] =
            +1.51502419e+01
            +2.14093466e-02 * tc[1]
            -7.36059614e-06 * tc[2]
            +1.14656816e-09 * tc[3]
            -6.66712722e-14 * tc[4];
    }

    /*species with midpoint at T=1399 kelvin */
    if (T < 1399) {
        /*species 49: C7H14O */
        species[49] =
            -6.39477206e+00
            +1.01846577e-01 * tc[1]
            -7.60046236e-05 * tc[2]
            +2.96538020e-08 * tc[3]
            -4.74854370e-12 * tc[4];
    } else {
        /*species 49: C7H14O */
        species[49] =
            +2.22692731e+01
            +3.25585391e-02 * tc[1]
            -1.11624906e-05 * tc[2]
            +1.73573658e-09 * tc[3]
            -1.00813215e-13 * tc[4];
    }

    /*species with midpoint at T=1387 kelvin */
    if (T < 1387) {
        /*species 50: NC7KET */
        species[50] =
            +5.29366920e-01
            +9.58173466e-02 * tc[1]
            -6.96688520e-05 * tc[2]
            +2.69540382e-08 * tc[3]
            -4.38728126e-12 * tc[4];
    } else {
        /*species 50: NC7KET */
        species[50] =
            +2.78330388e+01
            +3.20168096e-02 * tc[1]
            -1.11508456e-05 * tc[2]
            +1.75226159e-09 * tc[3]
            -1.02520451e-13 * tc[4];
    }
    return;
}


/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void cp_R(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.50104422e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            +3.29812400e+00
            +8.24944200e-04 * tc[1]
            -8.14301500e-07 * tc[2]
            -9.47543400e-11 * tc[3]
            +4.13487200e-13 * tc[4];
        /*species 2: O */
        species[2] =
            +2.94642900e+00
            -1.63816600e-03 * tc[1]
            +2.42103200e-06 * tc[2]
            -1.60284300e-09 * tc[3]
            +3.89069600e-13 * tc[4];
        /*species 3: O2 */
        species[3] =
            +3.21293600e+00
            +1.12748600e-03 * tc[1]
            -5.75615000e-07 * tc[2]
            +1.31387700e-09 * tc[3]
            -8.76855400e-13 * tc[4];
        /*species 5: H2O */
        species[5] =
            +3.38684200e+00
            +3.47498200e-03 * tc[1]
            -6.35469600e-06 * tc[2]
            +6.96858100e-09 * tc[3]
            -2.50658800e-12 * tc[4];
        /*species 6: CO */
        species[6] =
            +3.26245200e+00
            +1.51194100e-03 * tc[1]
            -3.88175500e-06 * tc[2]
            +5.58194400e-09 * tc[3]
            -2.47495100e-12 * tc[4];
        /*species 7: CO2 */
        species[7] =
            +2.27572500e+00
            +9.92207200e-03 * tc[1]
            -1.04091100e-05 * tc[2]
            +6.86668700e-09 * tc[3]
            -2.11728000e-12 * tc[4];
        /*species 8: CH3 */
        species[8] =
            +2.43044300e+00
            +1.11241000e-02 * tc[1]
            -1.68022000e-05 * tc[2]
            +1.62182900e-08 * tc[3]
            -5.86495300e-12 * tc[4];
        /*species 9: CH4 */
        species[9] =
            +7.78741500e-01
            +1.74766800e-02 * tc[1]
            -2.78340900e-05 * tc[2]
            +3.04970800e-08 * tc[3]
            -1.22393100e-11 * tc[4];
        /*species 11: H2O2 */
        species[11] =
            +3.38875400e+00
            +6.56922600e-03 * tc[1]
            -1.48501300e-07 * tc[2]
            -4.62580600e-09 * tc[3]
            +2.47151500e-12 * tc[4];
        /*species 12: CH2O */
        species[12] =
            +1.65273100e+00
            +1.26314400e-02 * tc[1]
            -1.88816800e-05 * tc[2]
            +2.05003100e-08 * tc[3]
            -8.41323700e-12 * tc[4];
        /*species 14: C2H4 */
        species[14] =
            -8.61488000e-01
            +2.79616300e-02 * tc[1]
            -3.38867700e-05 * tc[2]
            +2.78515200e-08 * tc[3]
            -9.73787900e-12 * tc[4];
        /*species 15: C2H5 */
        species[15] =
            +2.69070200e+00
            +8.71913300e-03 * tc[1]
            +4.41983900e-06 * tc[2]
            +9.33870300e-10 * tc[3]
            -3.92777300e-12 * tc[4];
        /*species 16: C2H */
        species[16] =
            +2.73770400e+00
            +8.04844600e-03 * tc[1]
            -9.24431000e-06 * tc[2]
            +6.52525900e-09 * tc[3]
            -1.93958000e-12 * tc[4];
        /*species 17: C2H2 */
        species[17] =
            +2.01356200e+00
            +1.51904500e-02 * tc[1]
            -1.61631900e-05 * tc[2]
            +9.07899200e-09 * tc[3]
            -1.91274600e-12 * tc[4];
        /*species 18: CH3OH */
        species[18] =
            +2.66011500e+00
            +7.34150800e-03 * tc[1]
            +7.17005100e-06 * tc[2]
            -8.79319400e-09 * tc[3]
            +2.39057000e-12 * tc[4];
        /*species 19: CH2CO */
        species[19] =
            +2.97497100e+00
            +1.21187100e-02 * tc[1]
            -2.34504600e-06 * tc[2]
            -6.46668500e-09 * tc[3]
            +3.90564900e-12 * tc[4];
        /*species 20: HCCO */
        species[20] =
            +5.04796500e+00
            +4.45347800e-03 * tc[1]
            +2.26828300e-07 * tc[2]
            -1.48209500e-09 * tc[3]
            +2.25074200e-13 * tc[4];
        /*species 21: CH3CHO */
        species[21] =
            +2.50569500e+00
            +1.33699100e-02 * tc[1]
            +4.67195300e-06 * tc[2]
            -1.12814000e-08 * tc[3]
            +4.26356600e-12 * tc[4];
        /*species 39: C3H3 */
        species[39] =
            +4.75420000e+00
            +1.10802800e-02 * tc[1]
            +2.79332300e-07 * tc[2]
            -5.47921200e-09 * tc[3]
            +1.94962900e-12 * tc[4];
        /*species 51: N2 */
        species[51] =
            +3.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.50104422e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            +2.99142300e+00
            +7.00064400e-04 * tc[1]
            -5.63382900e-08 * tc[2]
            -9.23157800e-12 * tc[3]
            +1.58275200e-15 * tc[4];
        /*species 2: O */
        species[2] =
            +2.54206000e+00
            -2.75506200e-05 * tc[1]
            -3.10280300e-09 * tc[2]
            +4.55106700e-12 * tc[3]
            -4.36805200e-16 * tc[4];
        /*species 3: O2 */
        species[3] =
            +3.69757800e+00
            +6.13519700e-04 * tc[1]
            -1.25884200e-07 * tc[2]
            +1.77528100e-11 * tc[3]
            -1.13643500e-15 * tc[4];
        /*species 5: H2O */
        species[5] =
            +2.67214600e+00
            +3.05629300e-03 * tc[1]
            -8.73026000e-07 * tc[2]
            +1.20099600e-10 * tc[3]
            -6.39161800e-15 * tc[4];
        /*species 6: CO */
        species[6] =
            +3.02507800e+00
            +1.44268900e-03 * tc[1]
            -5.63082800e-07 * tc[2]
            +1.01858100e-10 * tc[3]
            -6.91095200e-15 * tc[4];
        /*species 7: CO2 */
        species[7] =
            +4.45362300e+00
            +3.14016900e-03 * tc[1]
            -1.27841100e-06 * tc[2]
            +2.39399700e-10 * tc[3]
            -1.66903300e-14 * tc[4];
        /*species 8: CH3 */
        species[8] =
            +2.84405200e+00
            +6.13797400e-03 * tc[1]
            -2.23034500e-06 * tc[2]
            +3.78516100e-10 * tc[3]
            -2.45215900e-14 * tc[4];
        /*species 9: CH4 */
        species[9] =
            +1.68347900e+00
            +1.02372400e-02 * tc[1]
            -3.87512900e-06 * tc[2]
            +6.78558500e-10 * tc[3]
            -4.50342300e-14 * tc[4];
        /*species 11: H2O2 */
        species[11] =
            +4.57316700e+00
            +4.33613600e-03 * tc[1]
            -1.47468900e-06 * tc[2]
            +2.34890400e-10 * tc[3]
            -1.43165400e-14 * tc[4];
        /*species 12: CH2O */
        species[12] =
            +2.99560600e+00
            +6.68132100e-03 * tc[1]
            -2.62895500e-06 * tc[2]
            +4.73715300e-10 * tc[3]
            -3.21251700e-14 * tc[4];
        /*species 14: C2H4 */
        species[14] =
            +3.52841900e+00
            +1.14851800e-02 * tc[1]
            -4.41838500e-06 * tc[2]
            +7.84460100e-10 * tc[3]
            -5.26684800e-14 * tc[4];
        /*species 15: C2H5 */
        species[15] =
            +7.19048000e+00
            +6.48407700e-03 * tc[1]
            -6.42806500e-07 * tc[2]
            -2.34787900e-10 * tc[3]
            +3.88087700e-14 * tc[4];
        /*species 16: C2H */
        species[16] =
            +3.98636700e+00
            +3.14312300e-03 * tc[1]
            -1.26724300e-06 * tc[2]
            +2.92436300e-10 * tc[3]
            -2.71632000e-14 * tc[4];
        /*species 17: C2H2 */
        species[17] =
            +4.43677000e+00
            +5.37603900e-03 * tc[1]
            -1.91281700e-06 * tc[2]
            +3.28637900e-10 * tc[3]
            -2.15671000e-14 * tc[4];
        /*species 18: CH3OH */
        species[18] =
            +4.02906100e+00
            +9.37659300e-03 * tc[1]
            -3.05025400e-06 * tc[2]
            +4.35879300e-10 * tc[3]
            -2.22472300e-14 * tc[4];
        /*species 19: CH2CO */
        species[19] =
            +6.03881700e+00
            +5.80484000e-03 * tc[1]
            -1.92095400e-06 * tc[2]
            +2.79448500e-10 * tc[3]
            -1.45886800e-14 * tc[4];
        /*species 20: HCCO */
        species[20] =
            +6.75807300e+00
            +2.00040000e-03 * tc[1]
            -2.02760700e-07 * tc[2]
            -1.04113200e-10 * tc[3]
            +1.96516500e-14 * tc[4];
        /*species 21: CH3CHO */
        species[21] =
            +5.86865000e+00
            +1.07942400e-02 * tc[1]
            -3.64553000e-06 * tc[2]
            +5.41291200e-10 * tc[3]
            -2.89684400e-14 * tc[4];
        /*species 39: C3H3 */
        species[39] =
            +8.83104700e+00
            +4.35719500e-03 * tc[1]
            -4.10906700e-07 * tc[2]
            -2.36872300e-10 * tc[3]
            +4.37652000e-14 * tc[4];
        /*species 51: N2 */
        species[51] =
            +2.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }

    /*species with midpoint at T=1357 kelvin */
    if (T < 1357) {
        /*species 4: OH */
        species[4] =
            +3.43586219e+00
            +2.02235804e-04 * tc[1]
            -1.13546412e-07 * tc[2]
            +2.42445149e-10 * tc[3]
            -7.43651031e-14 * tc[4];
    } else {
        /*species 4: OH */
        species[4] =
            +2.62599754e+00
            +1.31992406e-03 * tc[1]
            -3.59724670e-07 * tc[2]
            +4.25630800e-11 * tc[3]
            -1.82048016e-15 * tc[4];
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 10: HO2 */
        species[10] =
            +3.18310656e+00
            +3.66767950e-03 * tc[1]
            -9.32385122e-07 * tc[2]
            -3.25852919e-10 * tc[3]
            +1.51139912e-13 * tc[4];
        /*species 35: CH3O2H */
        species[35] =
            +3.23442817e+00
            +1.90129767e-02 * tc[1]
            -1.13386287e-05 * tc[2]
            +3.40306653e-09 * tc[3]
            -4.11830222e-13 * tc[4];
        /*species 46: C7H15-1 */
        species[46] =
            -4.99570406e-01
            +8.08826467e-02 * tc[1]
            -5.00532754e-05 * tc[2]
            +1.56549308e-08 * tc[3]
            -1.96616227e-12 * tc[4];
        /*species 47: C7H15O2 */
        species[47] =
            +3.55252917e-01
            +9.42381007e-02 * tc[1]
            -6.66955106e-05 * tc[2]
            +2.54795583e-08 * tc[3]
            -4.13211105e-12 * tc[4];
        /*species 48: C7H14OOHO2 */
        species[48] =
            +1.40912479e+00
            +1.10720834e-01 * tc[1]
            -8.71178075e-05 * tc[2]
            +3.63662749e-08 * tc[3]
            -6.28350689e-12 * tc[4];
    } else {
        /*species 10: HO2 */
        species[10] =
            +4.10547423e+00
            +2.38452835e-03 * tc[1]
            -8.06347989e-07 * tc[2]
            +1.24191723e-10 * tc[3]
            -7.16400108e-15 * tc[4];
        /*species 35: CH3O2H */
        species[35] =
            +8.43117091e+00
            +8.06817909e-03 * tc[1]
            -2.77094921e-06 * tc[2]
            +4.31332243e-10 * tc[3]
            -2.50692146e-14 * tc[4];
        /*species 46: C7H15-1 */
        species[46] =
            +2.17940709e+01
            +3.26280243e-02 * tc[1]
            -1.11138244e-05 * tc[2]
            +1.72067148e-09 * tc[3]
            -9.96366999e-14 * tc[4];
        /*species 47: C7H15O2 */
        species[47] =
            +2.60535640e+01
            +3.43831688e-02 * tc[1]
            -1.18713247e-05 * tc[2]
            +1.85450241e-09 * tc[3]
            -1.08052399e-13 * tc[4];
        /*species 48: C7H14OOHO2 */
        species[48] =
            +3.28358296e+01
            +3.32767285e-02 * tc[1]
            -1.15833285e-05 * tc[2]
            +1.81958638e-09 * tc[3]
            -1.06434419e-13 * tc[4];
    }

    /*species with midpoint at T=1384 kelvin */
    if (T < 1384) {
        /*species 13: C2H6 */
        species[13] =
            -2.52854344e-02
            +2.40764754e-02 * tc[1]
            -1.11893472e-05 * tc[2]
            +2.08340901e-09 * tc[3]
            -5.29868616e-14 * tc[4];
    } else {
        /*species 13: C2H6 */
        species[13] =
            +6.10683385e+00
            +1.29236361e-02 * tc[1]
            -4.42527196e-06 * tc[2]
            +6.87391726e-10 * tc[3]
            -3.98901732e-14 * tc[4];
    }

    /*species with midpoint at T=1400 kelvin */
    if (T < 1400) {
        /*species 22: C3H4-A */
        species[22] =
            +2.53983100e+00
            +1.63343700e-02 * tc[1]
            -1.76495000e-06 * tc[2]
            -4.64736500e-09 * tc[3]
            +1.72913100e-12 * tc[4];
        /*species 23: C3H4-P */
        species[23] =
            +3.02973000e+00
            +1.49896100e-02 * tc[1]
            -1.39850000e-06 * tc[2]
            -3.96961900e-09 * tc[3]
            +1.38821700e-12 * tc[4];
        /*species 42: CH3CHCO */
        species[42] =
            +1.48380119e+00
            +3.22203013e-02 * tc[1]
            -2.70250033e-05 * tc[2]
            +1.20499164e-08 * tc[3]
            -2.18365931e-12 * tc[4];
    } else {
        /*species 22: C3H4-A */
        species[22] =
            +9.77625600e+00
            +5.30213800e-03 * tc[1]
            -3.70111800e-07 * tc[2]
            -3.02638600e-10 * tc[3]
            +5.08958100e-14 * tc[4];
        /*species 23: C3H4-P */
        species[23] =
            +9.76810200e+00
            +5.21915100e-03 * tc[1]
            -3.75314000e-07 * tc[2]
            -2.99219100e-10 * tc[3]
            +5.10787800e-14 * tc[4];
        /*species 42: CH3CHCO */
        species[42] =
            +1.00219123e+01
            +9.56966300e-03 * tc[1]
            -3.26221644e-06 * tc[2]
            +5.05231706e-10 * tc[3]
            -2.92593257e-14 * tc[4];
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 24: C3H6 */
        species[24] =
            +3.94615444e-01
            +2.89107662e-02 * tc[1]
            -1.54886808e-05 * tc[2]
            +3.88814209e-09 * tc[3]
            -3.37890352e-13 * tc[4];
    } else {
        /*species 24: C3H6 */
        species[24] =
            +8.01595958e+00
            +1.37023634e-02 * tc[1]
            -4.66249733e-06 * tc[2]
            +7.21254402e-10 * tc[3]
            -4.17370126e-14 * tc[4];
    }

    /*species with midpoint at T=1398 kelvin */
    if (T < 1398) {
        /*species 25: C4H6 */
        species[25] =
            -1.43095121e+00
            +4.78706062e-02 * tc[1]
            -4.15446800e-05 * tc[2]
            +1.91549552e-08 * tc[3]
            -3.57158507e-12 * tc[4];
    } else {
        /*species 25: C4H6 */
        species[25] =
            +1.11633789e+01
            +1.37163965e-02 * tc[1]
            -4.69715783e-06 * tc[2]
            +7.29693836e-10 * tc[3]
            -4.23486203e-14 * tc[4];
    }

    /*species with midpoint at T=1392 kelvin */
    if (T < 1392) {
        /*species 26: C4H7 */
        species[26] =
            -3.50508352e-01
            +4.26511243e-02 * tc[1]
            -2.90979373e-05 * tc[2]
            +1.05403914e-08 * tc[3]
            -1.60059854e-12 * tc[4];
        /*species 27: C4H8-1 */
        species[27] =
            -8.31372089e-01
            +4.52580978e-02 * tc[1]
            -2.93658559e-05 * tc[2]
            +1.00220436e-08 * tc[3]
            -1.43191680e-12 * tc[4];
        /*species 31: C5H9 */
        species[31] =
            -1.38013950e+00
            +5.57608487e-02 * tc[1]
            -3.70143928e-05 * tc[2]
            +1.26883901e-08 * tc[3]
            -1.78538835e-12 * tc[4];
        /*species 32: C5H10-1 */
        species[32] =
            -1.06223481e+00
            +5.74218294e-02 * tc[1]
            -3.74486890e-05 * tc[2]
            +1.27364989e-08 * tc[3]
            -1.79609789e-12 * tc[4];
    } else {
        /*species 26: C4H7 */
        species[26] =
            +1.12103578e+01
            +1.60483196e-02 * tc[1]
            -5.46502292e-06 * tc[2]
            +8.45941053e-10 * tc[3]
            -4.89772739e-14 * tc[4];
        /*species 27: C4H8-1 */
        species[27] =
            +1.13508668e+01
            +1.80617877e-02 * tc[1]
            -6.16093029e-06 * tc[2]
            +9.54652959e-10 * tc[3]
            -5.53089641e-14 * tc[4];
        /*species 31: C5H9 */
        species[31] =
            +1.41860454e+01
            +2.07128899e-02 * tc[1]
            -7.06960617e-06 * tc[2]
            +1.09607133e-09 * tc[3]
            -6.35322208e-14 * tc[4];
        /*species 32: C5H10-1 */
        species[32] =
            +1.45851539e+01
            +2.24072471e-02 * tc[1]
            -7.63348025e-06 * tc[2]
            +1.18188966e-09 * tc[3]
            -6.84385139e-14 * tc[4];
    }

    /*species with midpoint at T=1395 kelvin */
    if (T < 1395) {
        /*species 28: PC4H9 */
        species[28] =
            -4.37779725e-01
            +4.78972364e-02 * tc[1]
            -3.14023159e-05 * tc[2]
            +1.09786472e-08 * tc[3]
            -1.62010664e-12 * tc[4];
    } else {
        /*species 28: PC4H9 */
        species[28] =
            +1.21510082e+01
            +1.94310717e-02 * tc[1]
            -6.61577950e-06 * tc[2]
            +1.02375136e-09 * tc[3]
            -5.92529707e-14 * tc[4];
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 29: CH3COCH2 */
        species[29] =
            +1.80339187e+00
            +3.01407085e-02 * tc[1]
            -1.93505552e-05 * tc[2]
            +6.38199034e-09 * tc[3]
            -8.66103180e-13 * tc[4];
        /*species 45: NC7H16 */
        species[45] =
            -1.26836187e+00
            +8.54355820e-02 * tc[1]
            -5.25346786e-05 * tc[2]
            +1.62945721e-08 * tc[3]
            -2.02394925e-12 * tc[4];
    } else {
        /*species 29: CH3COCH2 */
        species[29] =
            +1.02303975e+01
            +1.16494161e-02 * tc[1]
            -4.01005537e-06 * tc[2]
            +6.25205246e-10 * tc[3]
            -3.63784362e-14 * tc[4];
        /*species 45: NC7H16 */
        species[45] =
            +2.22148969e+01
            +3.47675750e-02 * tc[1]
            -1.18407129e-05 * tc[2]
            +1.83298478e-09 * tc[3]
            -1.06130266e-13 * tc[4];
    }

    /*species with midpoint at T=1378 kelvin */
    if (T < 1378) {
        /*species 30: C2H5CHO */
        species[30] =
            +2.16308444e+00
            +2.95501264e-02 * tc[1]
            -1.52446252e-05 * tc[2]
            +3.49503947e-09 * tc[3]
            -2.38896627e-13 * tc[4];
        /*species 40: NC3H7CHO */
        species[40] =
            +1.87415959e+00
            +4.19240315e-02 * tc[1]
            -2.35148779e-05 * tc[2]
            +6.26913673e-09 * tc[3]
            -6.09443908e-13 * tc[4];
    } else {
        /*species 30: C2H5CHO */
        species[30] =
            +1.02427695e+01
            +1.39641989e-02 * tc[1]
            -4.76248001e-06 * tc[2]
            +7.38105706e-10 * tc[3]
            -4.27759503e-14 * tc[4];
        /*species 40: NC3H7CHO */
        species[40] =
            +1.35988068e+01
            +1.81652474e-02 * tc[1]
            -6.17844458e-06 * tc[2]
            +9.55980208e-10 * tc[3]
            -5.53442958e-14 * tc[4];
    }

    /*species with midpoint at T=1397 kelvin */
    if (T < 1397) {
        /*species 33: C5H11-1 */
        species[33] =
            -9.05255912e-01
            +6.10632852e-02 * tc[1]
            -4.09491825e-05 * tc[2]
            +1.46093470e-08 * tc[3]
            -2.18859615e-12 * tc[4];
        /*species 38: C3H5-A */
        species[38] =
            -5.29131958e-01
            +3.34559100e-02 * tc[1]
            -2.53401027e-05 * tc[2]
            +1.02865754e-08 * tc[3]
            -1.73258340e-12 * tc[4];
    } else {
        /*species 33: C5H11-1 */
        species[33] =
            +1.53234740e+01
            +2.39041200e-02 * tc[1]
            -8.14771619e-06 * tc[2]
            +1.26176236e-09 * tc[3]
            -7.30677335e-14 * tc[4];
        /*species 38: C3H5-A */
        species[38] =
            +8.45883958e+00
            +1.12695483e-02 * tc[1]
            -3.83792864e-06 * tc[2]
            +5.94059119e-10 * tc[3]
            -3.43918030e-14 * tc[4];
    }

    /*species with midpoint at T=1385 kelvin */
    if (T < 1385) {
        /*species 34: CH3O2 */
        species[34] =
            +4.26146906e+00
            +1.00873599e-02 * tc[1]
            -3.21506184e-06 * tc[2]
            +2.09409267e-10 * tc[3]
            +4.18339103e-14 * tc[4];
    } else {
        /*species 34: CH3O2 */
        species[34] =
            +5.95787891e+00
            +7.90728626e-03 * tc[1]
            -2.68246234e-06 * tc[2]
            +4.13891337e-10 * tc[3]
            -2.39007330e-14 * tc[4];
    }

    /*species with midpoint at T=1402 kelvin */
    if (T < 1402) {
        /*species 36: C2H3CO */
        species[36] =
            +1.36242013e+00
            +3.15273972e-02 * tc[1]
            -3.00218935e-05 * tc[2]
            +1.48167112e-08 * tc[3]
            -2.87971530e-12 * tc[4];
    } else {
        /*species 36: C2H3CO */
        species[36] =
            +9.37467676e+00
            +7.91296900e-03 * tc[1]
            -2.67198280e-06 * tc[2]
            +4.11115430e-10 * tc[3]
            -2.36978981e-14 * tc[4];
    }

    /*species with midpoint at T=1393 kelvin */
    if (T < 1393) {
        /*species 37: C2H3CHO */
        species[37] =
            +2.92355162e-01
            +3.54321417e-02 * tc[1]
            -2.94936324e-05 * tc[2]
            +1.28100124e-08 * tc[3]
            -2.26144108e-12 * tc[4];
    } else {
        /*species 37: C2H3CHO */
        species[37] =
            +1.04184959e+01
            +9.48963321e-03 * tc[1]
            -3.29310529e-06 * tc[2]
            +5.16279203e-10 * tc[3]
            -3.01587291e-14 * tc[4];
    }

    /*species with midpoint at T=1383 kelvin */
    if (T < 1383) {
        /*species 41: C2H5COCH2 */
        species[41] =
            +1.54013856e+00
            +4.39486258e-02 * tc[1]
            -2.97002421e-05 * tc[2]
            +1.05495313e-08 * tc[3]
            -1.58598769e-12 * tc[4];
        /*species 43: NC3H7COCH3 */
        species[43] =
            +1.24037361e+00
            +5.42494047e-02 * tc[1]
            -3.14825332e-05 * tc[2]
            +9.14492321e-09 * tc[3]
            -1.07665566e-12 * tc[4];
    } else {
        /*species 41: C2H5COCH2 */
        species[41] =
            +1.42098738e+01
            +1.57866459e-02 * tc[1]
            -5.50529183e-06 * tc[2]
            +8.65870540e-10 * tc[3]
            -5.06913329e-14 * tc[4];
        /*species 43: NC3H7COCH3 */
        species[43] =
            +1.61841564e+01
            +2.35245881e-02 * tc[1]
            -8.09675055e-06 * tc[2]
            +1.26230690e-09 * tc[3]
            -7.34487839e-14 * tc[4];
    }

    /*species with midpoint at T=1389 kelvin */
    if (T < 1389) {
        /*species 44: NC3H7COCH2 */
        species[44] =
            +9.58299271e-01
            +5.68162532e-02 * tc[1]
            -3.99112781e-05 * tc[2]
            +1.52671514e-08 * tc[3]
            -2.49221047e-12 * tc[4];
    } else {
        /*species 44: NC3H7COCH2 */
        species[44] =
            +1.61502419e+01
            +2.14093466e-02 * tc[1]
            -7.36059614e-06 * tc[2]
            +1.14656816e-09 * tc[3]
            -6.66712722e-14 * tc[4];
    }

    /*species with midpoint at T=1399 kelvin */
    if (T < 1399) {
        /*species 49: C7H14O */
        species[49] =
            -5.39477206e+00
            +1.01846577e-01 * tc[1]
            -7.60046236e-05 * tc[2]
            +2.96538020e-08 * tc[3]
            -4.74854370e-12 * tc[4];
    } else {
        /*species 49: C7H14O */
        species[49] =
            +2.32692731e+01
            +3.25585391e-02 * tc[1]
            -1.11624906e-05 * tc[2]
            +1.73573658e-09 * tc[3]
            -1.00813215e-13 * tc[4];
    }

    /*species with midpoint at T=1387 kelvin */
    if (T < 1387) {
        /*species 50: NC7KET */
        species[50] =
            +1.52936692e+00
            +9.58173466e-02 * tc[1]
            -6.96688520e-05 * tc[2]
            +2.69540382e-08 * tc[3]
            -4.38728126e-12 * tc[4];
    } else {
        /*species 50: NC7KET */
        species[50] =
            +2.88330388e+01
            +3.20168096e-02 * tc[1]
            -1.11508456e-05 * tc[2]
            +1.75226159e-09 * tc[3]
            -1.02520451e-13 * tc[4];
    }
    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void gibbs(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];
    amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.547474660000000e+04 * invT
            +2.966385537000000e+00
            -2.501044220000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -1.012521000000000e+03 * invT
            +6.592218000000000e+00
            -3.298124000000000e+00 * tc[0]
            -4.124721000000000e-04 * tc[1]
            +1.357169166666667e-07 * tc[2]
            +7.896194999999999e-12 * tc[3]
            -2.067436000000000e-14 * tc[4];
        /*species 2: O */
        species[2] =
            +2.914764000000000e+04 * invT
            -1.756599999999997e-02
            -2.946429000000000e+00 * tc[0]
            +8.190830000000000e-04 * tc[1]
            -4.035053333333333e-07 * tc[2]
            +1.335702500000000e-10 * tc[3]
            -1.945348000000000e-14 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.005249000000000e+03 * invT
            -2.821802000000000e+00
            -3.212936000000000e+00 * tc[0]
            -5.637430000000000e-04 * tc[1]
            +9.593583333333333e-08 * tc[2]
            -1.094897500000000e-10 * tc[3]
            +4.384277000000000e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.020811000000000e+04 * invT
            +7.966090000000001e-01
            -3.386842000000000e+00 * tc[0]
            -1.737491000000000e-03 * tc[1]
            +1.059116000000000e-06 * tc[2]
            -5.807150833333333e-10 * tc[3]
            +1.253294000000000e-13 * tc[4];
        /*species 6: CO */
        species[6] =
            -1.431054000000000e+04 * invT
            -1.586445000000000e+00
            -3.262452000000000e+00 * tc[0]
            -7.559705000000000e-04 * tc[1]
            +6.469591666666667e-07 * tc[2]
            -4.651620000000000e-10 * tc[3]
            +1.237475500000000e-13 * tc[4];
        /*species 7: CO2 */
        species[7] =
            -4.837314000000000e+04 * invT
            -7.912765000000000e+00
            -2.275725000000000e+00 * tc[0]
            -4.961036000000000e-03 * tc[1]
            +1.734851666666667e-06 * tc[2]
            -5.722239166666667e-10 * tc[3]
            +1.058640000000000e-13 * tc[4];
        /*species 8: CH3 */
        species[8] =
            +1.642378000000000e+04 * invT
            -4.359351000000000e+00
            -2.430443000000000e+00 * tc[0]
            -5.562050000000000e-03 * tc[1]
            +2.800366666666666e-06 * tc[2]
            -1.351524166666667e-09 * tc[3]
            +2.932476500000000e-13 * tc[4];
        /*species 9: CH4 */
        species[9] =
            -9.825228999999999e+03 * invT
            -1.294344850000000e+01
            -7.787415000000000e-01 * tc[0]
            -8.738340000000001e-03 * tc[1]
            +4.639015000000000e-06 * tc[2]
            -2.541423333333333e-09 * tc[3]
            +6.119655000000000e-13 * tc[4];
        /*species 11: H2O2 */
        species[11] =
            -1.766315000000000e+04 * invT
            -3.396609000000000e+00
            -3.388754000000000e+00 * tc[0]
            -3.284613000000000e-03 * tc[1]
            +2.475021666666666e-08 * tc[2]
            +3.854838333333333e-10 * tc[3]
            -1.235757500000000e-13 * tc[4];
        /*species 12: CH2O */
        species[12] =
            -1.486540000000000e+04 * invT
            -1.213208900000000e+01
            -1.652731000000000e+00 * tc[0]
            -6.315720000000000e-03 * tc[1]
            +3.146946666666667e-06 * tc[2]
            -1.708359166666667e-09 * tc[3]
            +4.206618500000000e-13 * tc[4];
        /*species 14: C2H4 */
        species[14] =
            +5.573046000000000e+03 * invT
            -2.507297800000000e+01
            +8.614880000000000e-01 * tc[0]
            -1.398081500000000e-02 * tc[1]
            +5.647795000000000e-06 * tc[2]
            -2.320960000000000e-09 * tc[3]
            +4.868939500000000e-13 * tc[4];
        /*species 15: C2H5 */
        species[15] =
            +1.287040000000000e+04 * invT
            -9.447498000000000e+00
            -2.690702000000000e+00 * tc[0]
            -4.359566500000000e-03 * tc[1]
            -7.366398333333332e-07 * tc[2]
            -7.782252500000001e-11 * tc[3]
            +1.963886500000000e-13 * tc[4];
        /*species 16: C2H */
        species[16] =
            +6.683813000000000e+04 * invT
            -4.562516000000000e+00
            -2.737704000000000e+00 * tc[0]
            -4.024223000000000e-03 * tc[1]
            +1.540718333333333e-06 * tc[2]
            -5.437715833333333e-10 * tc[3]
            +9.697900000000001e-14 * tc[4];
        /*species 17: C2H2 */
        species[17] =
            +2.612444000000000e+04 * invT
            -6.791815999999999e+00
            -2.013562000000000e+00 * tc[0]
            -7.595225000000000e-03 * tc[1]
            +2.693865000000000e-06 * tc[2]
            -7.565826666666667e-10 * tc[3]
            +9.563730000000000e-14 * tc[4];
        /*species 18: CH3OH */
        species[18] =
            -2.535348000000000e+04 * invT
            -8.572515000000001e+00
            -2.660115000000000e+00 * tc[0]
            -3.670754000000000e-03 * tc[1]
            -1.195008500000000e-06 * tc[2]
            +7.327661666666667e-10 * tc[3]
            -1.195285000000000e-13 * tc[4];
        /*species 19: CH2CO */
        species[19] =
            -7.632637000000000e+03 * invT
            -5.698582000000000e+00
            -2.974971000000000e+00 * tc[0]
            -6.059355000000000e-03 * tc[1]
            +3.908410000000000e-07 * tc[2]
            +5.388904166666666e-10 * tc[3]
            -1.952824500000000e-13 * tc[4];
        /*species 20: HCCO */
        species[20] =
            +1.965892000000000e+04 * invT
            +4.566121099999999e+00
            -5.047965000000000e+00 * tc[0]
            -2.226739000000000e-03 * tc[1]
            -3.780471666666667e-08 * tc[2]
            +1.235079166666667e-10 * tc[3]
            -1.125371000000000e-14 * tc[4];
        /*species 21: CH3CHO */
        species[21] =
            -2.124589000000000e+04 * invT
            -1.084519500000000e+01
            -2.505695000000000e+00 * tc[0]
            -6.684955000000000e-03 * tc[1]
            -7.786588333333333e-07 * tc[2]
            +9.401166666666668e-10 * tc[3]
            -2.131783000000000e-13 * tc[4];
        /*species 39: C3H3 */
        species[39] =
            +3.988883000000000e+04 * invT
            +4.168745100000000e+00
            -4.754200000000000e+00 * tc[0]
            -5.540140000000000e-03 * tc[1]
            -4.655538333333333e-08 * tc[2]
            +4.566010000000000e-10 * tc[3]
            -9.748145000000000e-14 * tc[4];
        /*species 51: N2 */
        species[51] =
            -1.020900000000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041200000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427500000000e-13 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.547474660000000e+04 * invT
            +2.966385537000000e+00
            -2.501044220000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -8.350340000000000e+02 * invT
            +4.346533000000000e+00
            -2.991423000000000e+00 * tc[0]
            -3.500322000000000e-04 * tc[1]
            +9.389715000000000e-09 * tc[2]
            +7.692981666666667e-13 * tc[3]
            -7.913760000000000e-17 * tc[4];
        /*species 2: O */
        species[2] =
            +2.923080000000000e+04 * invT
            -2.378248000000000e+00
            -2.542060000000000e+00 * tc[0]
            +1.377531000000000e-05 * tc[1]
            +5.171338333333333e-10 * tc[2]
            -3.792555833333334e-13 * tc[3]
            +2.184026000000000e-17 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.233930000000000e+03 * invT
            +5.084119999999999e-01
            -3.697578000000000e+00 * tc[0]
            -3.067598500000000e-04 * tc[1]
            +2.098070000000000e-08 * tc[2]
            -1.479400833333333e-12 * tc[3]
            +5.682175000000001e-17 * tc[4];
        /*species 5: H2O */
        species[5] =
            -2.989921000000000e+04 * invT
            -4.190671000000000e+00
            -2.672146000000000e+00 * tc[0]
            -1.528146500000000e-03 * tc[1]
            +1.455043333333333e-07 * tc[2]
            -1.000830000000000e-11 * tc[3]
            +3.195809000000000e-16 * tc[4];
        /*species 6: CO */
        species[6] =
            -1.426835000000000e+04 * invT
            -3.083140000000000e+00
            -3.025078000000000e+00 * tc[0]
            -7.213445000000000e-04 * tc[1]
            +9.384713333333334e-08 * tc[2]
            -8.488174999999999e-12 * tc[3]
            +3.455476000000000e-16 * tc[4];
        /*species 7: CO2 */
        species[7] =
            -4.896696000000000e+04 * invT
            +5.409018900000000e+00
            -4.453623000000000e+00 * tc[0]
            -1.570084500000000e-03 * tc[1]
            +2.130685000000000e-07 * tc[2]
            -1.994997500000000e-11 * tc[3]
            +8.345165000000000e-16 * tc[4];
        /*species 8: CH3 */
        species[8] =
            +1.643781000000000e+04 * invT
            -2.608645000000000e+00
            -2.844052000000000e+00 * tc[0]
            -3.068987000000000e-03 * tc[1]
            +3.717241666666666e-07 * tc[2]
            -3.154300833333333e-11 * tc[3]
            +1.226079500000000e-15 * tc[4];
        /*species 9: CH4 */
        species[9] =
            -1.008079000000000e+04 * invT
            -7.939916000000000e+00
            -1.683479000000000e+00 * tc[0]
            -5.118620000000000e-03 * tc[1]
            +6.458548333333333e-07 * tc[2]
            -5.654654166666667e-11 * tc[3]
            +2.251711500000000e-15 * tc[4];
        /*species 11: H2O2 */
        species[11] =
            -1.800696000000000e+04 * invT
            +4.072030000000000e+00
            -4.573167000000000e+00 * tc[0]
            -2.168068000000000e-03 * tc[1]
            +2.457815000000000e-07 * tc[2]
            -1.957420000000000e-11 * tc[3]
            +7.158270000000000e-16 * tc[4];
        /*species 12: CH2O */
        species[12] =
            -1.532037000000000e+04 * invT
            -3.916966000000000e+00
            -2.995606000000000e+00 * tc[0]
            -3.340660500000000e-03 * tc[1]
            +4.381591666666666e-07 * tc[2]
            -3.947627500000000e-11 * tc[3]
            +1.606258500000000e-15 * tc[4];
        /*species 14: C2H4 */
        species[14] =
            +4.428289000000000e+03 * invT
            +1.298030000000000e+00
            -3.528419000000000e+00 * tc[0]
            -5.742590000000000e-03 * tc[1]
            +7.363975000000000e-07 * tc[2]
            -6.537167500000001e-11 * tc[3]
            +2.633424000000000e-15 * tc[4];
        /*species 15: C2H5 */
        species[15] =
            +1.067455000000000e+04 * invT
            +2.197137000000000e+01
            -7.190480000000000e+00 * tc[0]
            -3.242038500000000e-03 * tc[1]
            +1.071344166666667e-07 * tc[2]
            +1.956565833333333e-11 * tc[3]
            -1.940438500000000e-15 * tc[4];
        /*species 16: C2H */
        species[16] =
            +6.655884000000000e+04 * invT
            +2.795304000000000e+00
            -3.986367000000000e+00 * tc[0]
            -1.571561500000000e-03 * tc[1]
            +2.112071666666667e-07 * tc[2]
            -2.436969166666667e-11 * tc[3]
            +1.358160000000000e-15 * tc[4];
        /*species 17: C2H2 */
        species[17] =
            +2.566766000000000e+04 * invT
            +7.237108000000000e+00
            -4.436770000000000e+00 * tc[0]
            -2.688019500000000e-03 * tc[1]
            +3.188028333333333e-07 * tc[2]
            -2.738649166666667e-11 * tc[3]
            +1.078355000000000e-15 * tc[4];
        /*species 18: CH3OH */
        species[18] =
            -2.615791000000000e+04 * invT
            +1.650865000000000e+00
            -4.029061000000000e+00 * tc[0]
            -4.688296500000000e-03 * tc[1]
            +5.083756666666666e-07 * tc[2]
            -3.632327500000000e-11 * tc[3]
            +1.112361500000000e-15 * tc[4];
        /*species 19: CH2CO */
        species[19] =
            -8.583402000000000e+03 * invT
            +1.369639800000000e+01
            -6.038817000000000e+00 * tc[0]
            -2.902420000000000e-03 * tc[1]
            +3.201590000000000e-07 * tc[2]
            -2.328737500000000e-11 * tc[3]
            +7.294340000000000e-16 * tc[4];
        /*species 20: HCCO */
        species[20] =
            +1.901513000000000e+04 * invT
            +1.582933500000000e+01
            -6.758073000000000e+00 * tc[0]
            -1.000200000000000e-03 * tc[1]
            +3.379345000000000e-08 * tc[2]
            +8.676100000000000e-12 * tc[3]
            -9.825825000000000e-16 * tc[4];
        /*species 21: CH3CHO */
        species[21] =
            -2.264569000000000e+04 * invT
            +1.188159600000000e+01
            -5.868650000000000e+00 * tc[0]
            -5.397120000000000e-03 * tc[1]
            +6.075883333333333e-07 * tc[2]
            -4.510760000000000e-11 * tc[3]
            +1.448422000000000e-15 * tc[4];
        /*species 39: C3H3 */
        species[39] =
            +3.847420000000000e+04 * invT
            +3.061023700000000e+01
            -8.831047000000000e+00 * tc[0]
            -2.178597500000000e-03 * tc[1]
            +6.848445000000000e-08 * tc[2]
            +1.973935833333333e-11 * tc[3]
            -2.188260000000000e-15 * tc[4];
        /*species 51: N2 */
        species[51] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439885000000000e-04 * tc[1]
            +9.474601666666666e-08 * tc[2]
            -8.414199999999999e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
    }

    /*species with midpoint at T=1357 kelvin */
    if (T < 1357) {
        /*species 4: OH */
        species[4] =
            +3.743212520000000e+03 * invT
            +9.857209199999999e-01
            -3.435862190000000e+00 * tc[0]
            -1.011179020000000e-04 * tc[1]
            +1.892440200000000e-08 * tc[2]
            -2.020376241666666e-11 * tc[3]
            +3.718255155000000e-15 * tc[4];
    } else {
        /*species 4: OH */
        species[4] =
            +4.120853740000000e+03 * invT
            -4.480675530000000e+00
            -2.625997540000000e+00 * tc[0]
            -6.599620300000000e-04 * tc[1]
            +5.995411166666667e-08 * tc[2]
            -3.546923333333333e-12 * tc[3]
            +9.102400800000000e-17 * tc[4];
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 10: HO2 */
        species[10] =
            +8.091810130000000e+02 * invT
            -5.210604430000000e+00
            -3.183106560000000e+00 * tc[0]
            -1.833839750000000e-03 * tc[1]
            +1.553975203333333e-07 * tc[2]
            +2.715440991666667e-11 * tc[3]
            -7.556995600000000e-15 * tc[4];
        /*species 35: CH3O2H */
        species[35] =
            -1.771979260000000e+04 * invT
            -6.021811320000000e+00
            -3.234428170000000e+00 * tc[0]
            -9.506488350000000e-03 * tc[1]
            +1.889771450000000e-06 * tc[2]
            -2.835888775000000e-10 * tc[3]
            +2.059151110000000e-14 * tc[4];
        /*species 46: C7H15-1 */
        species[46] =
            -1.045902230000000e+03 * invT
            -3.515597150599999e+01
            +4.995704060000000e-01 * tc[0]
            -4.044132335000000e-02 * tc[1]
            +8.342212566666667e-06 * tc[2]
            -1.304577566666667e-09 * tc[3]
            +9.830811350000000e-14 * tc[4];
        /*species 47: C7H15O2 */
        species[47] =
            -2.157952540000000e+04 * invT
            -3.258809058300000e+01
            -3.552529170000000e-01 * tc[0]
            -4.711905035000000e-02 * tc[1]
            +1.111591843333333e-05 * tc[2]
            -2.123296525000000e-09 * tc[3]
            +2.066055525000000e-13 * tc[4];
        /*species 48: C7H14OOHO2 */
        species[48] =
            -3.463047880000000e+04 * invT
            -3.021321351000000e+01
            -1.409124790000000e+00 * tc[0]
            -5.536041700000000e-02 * tc[1]
            +1.451963458333333e-05 * tc[2]
            -3.030522908333333e-09 * tc[3]
            +3.141753445000000e-13 * tc[4];
    } else {
        /*species 10: HO2 */
        species[10] =
            +3.981276890000000e+02 * invT
            +9.803158699999996e-01
            -4.105474230000000e+00 * tc[0]
            -1.192264175000000e-03 * tc[1]
            +1.343913315000000e-07 * tc[2]
            -1.034931025000000e-11 * tc[3]
            +3.582000540000000e-16 * tc[4];
        /*species 35: CH3O2H */
        species[35] =
            -1.966787710000000e+04 * invT
            +2.754823381000000e+01
            -8.431170910000001e+00 * tc[0]
            -4.034089545000000e-03 * tc[1]
            +4.618248683333333e-07 * tc[2]
            -3.594435358333333e-11 * tc[3]
            +1.253460730000000e-15 * tc[4];
        /*species 46: C7H15-1 */
        species[46] =
            -9.209382210000000e+03 * invT
            +1.082895020000000e+02
            -2.179407090000000e+01 * tc[0]
            -1.631401215000000e-02 * tc[1]
            +1.852304066666667e-06 * tc[2]
            -1.433892900000000e-10 * tc[3]
            +4.981834995000000e-15 * tc[4];
        /*species 47: C7H15O2 */
        species[47] =
            -3.068420950000000e+04 * invT
            +1.314611270000000e+02
            -2.605356400000000e+01 * tc[0]
            -1.719158440000000e-02 * tc[1]
            +1.978554116666667e-06 * tc[2]
            -1.545418675000000e-10 * tc[3]
            +5.402619949999999e-15 * tc[4];
        /*species 48: C7H14OOHO2 */
        species[48] =
            -4.528730030000000e+04 * invT
            +1.688115396000000e+02
            -3.283582960000000e+01 * tc[0]
            -1.663836425000000e-02 * tc[1]
            +1.930554750000000e-06 * tc[2]
            -1.516321983333333e-10 * tc[3]
            +5.321720950000000e-15 * tc[4];
    }

    /*species with midpoint at T=1384 kelvin */
    if (T < 1384) {
        /*species 13: C2H6 */
        species[13] =
            -1.123455340000000e+04 * invT
            -2.119016043440000e+01
            +2.528543440000000e-02 * tc[0]
            -1.203823770000000e-02 * tc[1]
            +1.864891200000000e-06 * tc[2]
            -1.736174175000000e-10 * tc[3]
            +2.649343080000000e-15 * tc[4];
    } else {
        /*species 13: C2H6 */
        species[13] =
            -1.375000140000000e+04 * invT
            +1.911495885000000e+01
            -6.106833850000000e+00 * tc[0]
            -6.461818050000000e-03 * tc[1]
            +7.375453266666666e-07 * tc[2]
            -5.728264383333333e-11 * tc[3]
            +1.994508660000000e-15 * tc[4];
    }

    /*species with midpoint at T=1400 kelvin */
    if (T < 1400) {
        /*species 22: C3H4-A */
        species[22] =
            +2.251243000000000e+04 * invT
            -7.395871000000000e+00
            -2.539831000000000e+00 * tc[0]
            -8.167185000000000e-03 * tc[1]
            +2.941583333333333e-07 * tc[2]
            +3.872804166666666e-10 * tc[3]
            -8.645655000000001e-14 * tc[4];
        /*species 23: C3H4-P */
        species[23] =
            +2.148408000000000e+04 * invT
            -4.974864000000001e+00
            -3.029730000000000e+00 * tc[0]
            -7.494805000000000e-03 * tc[1]
            +2.330833333333333e-07 * tc[2]
            +3.308015833333333e-10 * tc[3]
            -6.941085000000000e-14 * tc[4];
        /*species 42: CH3CHCO */
        species[42] =
            -1.152765400000000e+04 * invT
            -1.567140561000000e+01
            -1.483801190000000e+00 * tc[0]
            -1.611015065000000e-02 * tc[1]
            +4.504167216666666e-06 * tc[2]
            -1.004159700000000e-09 * tc[3]
            +1.091829655000000e-13 * tc[4];
    } else {
        /*species 22: C3H4-A */
        species[22] =
            +1.954972000000000e+04 * invT
            +4.054686600000000e+01
            -9.776256000000000e+00 * tc[0]
            -2.651069000000000e-03 * tc[1]
            +6.168530000000000e-08 * tc[2]
            +2.521988333333334e-11 * tc[3]
            -2.544790500000000e-15 * tc[4];
        /*species 23: C3H4-P */
        species[23] =
            +1.860277000000000e+04 * invT
            +3.997488200000000e+01
            -9.768102000000001e+00 * tc[0]
            -2.609575500000000e-03 * tc[1]
            +6.255233333333333e-08 * tc[2]
            +2.493492500000000e-11 * tc[3]
            -2.553939000000000e-15 * tc[4];
        /*species 42: CH3CHCO */
        species[42] =
            -1.424827380000000e+04 * invT
            +3.780490960000000e+01
            -1.002191230000000e+01 * tc[0]
            -4.784831500000000e-03 * tc[1]
            +5.437027400000000e-07 * tc[2]
            -4.210264216666667e-11 * tc[3]
            +1.462966285000000e-15 * tc[4];
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 24: C3H6 */
        species[24] =
            +1.066881640000000e+03 * invT
            -2.150575815600000e+01
            -3.946154440000000e-01 * tc[0]
            -1.445538310000000e-02 * tc[1]
            +2.581446800000000e-06 * tc[2]
            -3.240118408333333e-10 * tc[3]
            +1.689451760000000e-14 * tc[4];
    } else {
        /*species 24: C3H6 */
        species[24] =
            -1.878212710000000e+03 * invT
            +2.803202638000000e+01
            -8.015959580000001e+00 * tc[0]
            -6.851181700000000e-03 * tc[1]
            +7.770828883333333e-07 * tc[2]
            -6.010453350000000e-11 * tc[3]
            +2.086850630000000e-15 * tc[4];
    }

    /*species with midpoint at T=1398 kelvin */
    if (T < 1398) {
        /*species 25: C4H6 */
        species[25] =
            +1.175513140000000e+04 * invT
            -3.051353451000000e+01
            +1.430951210000000e+00 * tc[0]
            -2.393530310000000e-02 * tc[1]
            +6.924113333333333e-06 * tc[2]
            -1.596246266666667e-09 * tc[3]
            +1.785792535000000e-13 * tc[4];
    } else {
        /*species 25: C4H6 */
        species[25] =
            +7.790397700000000e+03 * invT
            +4.814817380000000e+01
            -1.116337890000000e+01 * tc[0]
            -6.858198250000000e-03 * tc[1]
            +7.828596383333334e-07 * tc[2]
            -6.080781966666666e-11 * tc[3]
            +2.117431015000000e-15 * tc[4];
    }

    /*species with midpoint at T=1392 kelvin */
    if (T < 1392) {
        /*species 26: C4H7 */
        species[26] =
            +1.499335910000000e+04 * invT
            -2.708007795200000e+01
            +3.505083520000000e-01 * tc[0]
            -2.132556215000000e-02 * tc[1]
            +4.849656216666667e-06 * tc[2]
            -8.783659500000000e-10 * tc[3]
            +8.002992699999999e-14 * tc[4];
        /*species 27: C4H8-1 */
        species[27] =
            -1.578750350000000e+03 * invT
            -3.033979568900000e+01
            +8.313720890000000e-01 * tc[0]
            -2.262904890000000e-02 * tc[1]
            +4.894309316666667e-06 * tc[2]
            -8.351702999999999e-10 * tc[3]
            +7.159584000000000e-14 * tc[4];
        /*species 31: C5H9 */
        species[31] =
            +1.255898240000000e+04 * invT
            -3.402426990000000e+01
            +1.380139500000000e+00 * tc[0]
            -2.788042435000000e-02 * tc[1]
            +6.169065466666666e-06 * tc[2]
            -1.057365841666667e-09 * tc[3]
            +8.926941750000000e-14 * tc[4];
        /*species 32: C5H10-1 */
        species[32] =
            -4.465466660000000e+03 * invT
            -3.333621381000000e+01
            +1.062234810000000e+00 * tc[0]
            -2.871091470000000e-02 * tc[1]
            +6.241448166666667e-06 * tc[2]
            -1.061374908333333e-09 * tc[3]
            +8.980489449999999e-14 * tc[4];
    } else {
        /*species 26: C4H7 */
        species[26] =
            +1.090419370000000e+04 * invT
            +4.676965930000000e+01
            -1.121035780000000e+01 * tc[0]
            -8.024159800000000e-03 * tc[1]
            +9.108371533333334e-07 * tc[2]
            -7.049508775000000e-11 * tc[3]
            +2.448863695000000e-15 * tc[4];
        /*species 27: C4H8-1 */
        species[27] =
            -5.978710380000000e+03 * invT
            +4.778781060000000e+01
            -1.135086680000000e+01 * tc[0]
            -9.030893850000001e-03 * tc[1]
            +1.026821715000000e-06 * tc[2]
            -7.955441325000001e-11 * tc[3]
            +2.765448205000000e-15 * tc[4];
        /*species 31: C5H9 */
        species[31] =
            +7.004961350000000e+03 * invT
            +6.563622270000000e+01
            -1.418604540000000e+01 * tc[0]
            -1.035644495000000e-02 * tc[1]
            +1.178267695000000e-06 * tc[2]
            -9.133927750000000e-11 * tc[3]
            +3.176611040000000e-15 * tc[4];
        /*species 32: C5H10-1 */
        species[32] =
            -1.008982050000000e+04 * invT
            +6.695354750000000e+01
            -1.458515390000000e+01 * tc[0]
            -1.120362355000000e-02 * tc[1]
            +1.272246708333333e-06 * tc[2]
            -9.849080500000001e-11 * tc[3]
            +3.421925695000000e-15 * tc[4];
    }

    /*species with midpoint at T=1395 kelvin */
    if (T < 1395) {
        /*species 28: PC4H9 */
        species[28] =
            +7.689452480000000e+03 * invT
            -2.912305292500000e+01
            +4.377797250000000e-01 * tc[0]
            -2.394861820000000e-02 * tc[1]
            +5.233719316666667e-06 * tc[2]
            -9.148872666666667e-10 * tc[3]
            +8.100533200000001e-14 * tc[4];
    } else {
        /*species 28: PC4H9 */
        species[28] =
            +3.172319420000000e+03 * invT
            +5.149359040000000e+01
            -1.215100820000000e+01 * tc[0]
            -9.715535850000000e-03 * tc[1]
            +1.102629916666667e-06 * tc[2]
            -8.531261333333333e-11 * tc[3]
            +2.962648535000000e-15 * tc[4];
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 29: CH3COCH2 */
        species[29] =
            -5.372332610000000e+03 * invT
            -1.600124893000000e+01
            -1.803391870000000e+00 * tc[0]
            -1.507035425000000e-02 * tc[1]
            +3.225092533333333e-06 * tc[2]
            -5.318325283333333e-10 * tc[3]
            +4.330515900000000e-14 * tc[4];
        /*species 45: NC7H16 */
        species[45] =
            -2.565865650000000e+04 * invT
            -3.664165307000000e+01
            +1.268361870000000e+00 * tc[0]
            -4.271779100000000e-02 * tc[1]
            +8.755779766666667e-06 * tc[2]
            -1.357881008333333e-09 * tc[3]
            +1.011974625000000e-13 * tc[4];
    } else {
        /*species 29: CH3COCH2 */
        species[29] =
            -8.443804620000001e+03 * invT
            +3.815011950000000e+01
            -1.023039750000000e+01 * tc[0]
            -5.824708050000000e-03 * tc[1]
            +6.683425616666666e-07 * tc[2]
            -5.210043716666667e-11 * tc[3]
            +1.818921810000000e-15 * tc[4];
        /*species 45: NC7H16 */
        species[45] =
            -3.427600810000000e+04 * invT
            +1.145189165000000e+02
            -2.221489690000000e+01 * tc[0]
            -1.738378750000000e-02 * tc[1]
            +1.973452150000000e-06 * tc[2]
            -1.527487316666667e-10 * tc[3]
            +5.306513300000000e-15 * tc[4];
    }

    /*species with midpoint at T=1378 kelvin */
    if (T < 1378) {
        /*species 30: C2H5CHO */
        species[30] =
            -2.422601370000000e+04 * invT
            -1.395225036000000e+01
            -2.163084440000000e+00 * tc[0]
            -1.477506320000000e-02 * tc[1]
            +2.540770866666666e-06 * tc[2]
            -2.912532891666666e-10 * tc[3]
            +1.194483135000000e-14 * tc[4];
        /*species 40: NC3H7CHO */
        species[40] =
            -2.710321940000000e+04 * invT
            -1.728269781000000e+01
            -1.874159590000000e+00 * tc[0]
            -2.096201575000000e-02 * tc[1]
            +3.919146316666667e-06 * tc[2]
            -5.224280608333334e-10 * tc[3]
            +3.047219540000000e-14 * tc[4];
    } else {
        /*species 30: C2H5CHO */
        species[30] =
            -2.741431350000000e+04 * invT
            +3.877745380000000e+01
            -1.024276950000000e+01 * tc[0]
            -6.982099450000000e-03 * tc[1]
            +7.937466683333333e-07 * tc[2]
            -6.150880883333334e-11 * tc[3]
            +2.138797515000000e-15 * tc[4];
        /*species 40: NC3H7CHO */
        species[40] =
            -3.158453480000000e+04 * invT
            +5.877782960000000e+01
            -1.359880680000000e+01 * tc[0]
            -9.082623700000000e-03 * tc[1]
            +1.029740763333333e-06 * tc[2]
            -7.966501733333333e-11 * tc[3]
            +2.767214790000000e-15 * tc[4];
    }

    /*species with midpoint at T=1397 kelvin */
    if (T < 1397) {
        /*species 33: C5H11-1 */
        species[33] =
            +4.839953030000000e+03 * invT
            -3.346275221200000e+01
            +9.052559120000000e-01 * tc[0]
            -3.053164260000000e-02 * tc[1]
            +6.824863750000000e-06 * tc[2]
            -1.217445583333333e-09 * tc[3]
            +1.094298075000000e-13 * tc[4];
        /*species 38: C3H5-A */
        species[38] =
            +1.938342260000000e+04 * invT
            -2.583584505800000e+01
            +5.291319580000000e-01 * tc[0]
            -1.672795500000000e-02 * tc[1]
            +4.223350450000001e-06 * tc[2]
            -8.572146166666666e-10 * tc[3]
            +8.662917000000000e-14 * tc[4];
    } else {
        /*species 33: C5H11-1 */
        species[33] =
            -9.232416370000000e+02 * invT
            +7.027635990000000e+01
            -1.532347400000000e+01 * tc[0]
            -1.195206000000000e-02 * tc[1]
            +1.357952698333333e-06 * tc[2]
            -1.051468633333333e-10 * tc[3]
            +3.653386675000001e-15 * tc[4];
        /*species 38: C3H5-A */
        species[38] =
            +1.635760920000000e+04 * invT
            +3.103978458000000e+01
            -8.458839579999999e+00 * tc[0]
            -5.634774150000000e-03 * tc[1]
            +6.396547733333333e-07 * tc[2]
            -4.950492658333333e-11 * tc[3]
            +1.719590150000000e-15 * tc[4];
    }

    /*species with midpoint at T=1385 kelvin */
    if (T < 1385) {
        /*species 34: CH3O2 */
        species[34] =
            -6.843942590000000e+02 * invT
            -9.018341400000001e-01
            -4.261469060000000e+00 * tc[0]
            -5.043679950000000e-03 * tc[1]
            +5.358436400000000e-07 * tc[2]
            -1.745077225000000e-11 * tc[3]
            -2.091695515000000e-15 * tc[4];
    } else {
        /*species 34: CH3O2 */
        species[34] =
            -1.535748380000000e+03 * invT
            +1.067751777000000e+01
            -5.957878910000000e+00 * tc[0]
            -3.953643130000000e-03 * tc[1]
            +4.470770566666667e-07 * tc[2]
            -3.449094475000000e-11 * tc[3]
            +1.195036650000000e-15 * tc[4];
    }

    /*species with midpoint at T=1402 kelvin */
    if (T < 1402) {
        /*species 36: C2H3CO */
        species[36] =
            +4.257702150000000e+03 * invT
            -1.590023447000000e+01
            -1.362420130000000e+00 * tc[0]
            -1.576369860000000e-02 * tc[1]
            +5.003648916666666e-06 * tc[2]
            -1.234725933333333e-09 * tc[3]
            +1.439857650000000e-13 * tc[4];
    } else {
        /*species 36: C2H3CO */
        species[36] =
            +1.929695140000000e+03 * invT
            +3.346394636000000e+01
            -9.374676760000000e+00 * tc[0]
            -3.956484500000000e-03 * tc[1]
            +4.453304666666666e-07 * tc[2]
            -3.425961916666666e-11 * tc[3]
            +1.184894905000000e-15 * tc[4];
    }

    /*species with midpoint at T=1393 kelvin */
    if (T < 1393) {
        /*species 37: C2H3CHO */
        species[37] =
            -1.165215840000000e+04 * invT
            -2.259547283800000e+01
            -2.923551620000000e-01 * tc[0]
            -1.771607085000000e-02 * tc[1]
            +4.915605400000000e-06 * tc[2]
            -1.067501033333333e-09 * tc[3]
            +1.130720540000000e-13 * tc[4];
    } else {
        /*species 37: C2H3CHO */
        species[37] =
            -1.496302810000000e+04 * invT
            +4.114200200000001e+01
            -1.041849590000000e+01 * tc[0]
            -4.744816605000000e-03 * tc[1]
            +5.488508816666667e-07 * tc[2]
            -4.302326691666667e-11 * tc[3]
            +1.507936455000000e-15 * tc[4];
    }

    /*species with midpoint at T=1383 kelvin */
    if (T < 1383) {
        /*species 41: C2H5COCH2 */
        species[41] =
            -9.507965050000001e+03 * invT
            -1.843052554000000e+01
            -1.540138560000000e+00 * tc[0]
            -2.197431290000000e-02 * tc[1]
            +4.950040350000000e-06 * tc[2]
            -8.791276083333333e-10 * tc[3]
            +7.929938449999999e-14 * tc[4];
        /*species 43: NC3H7COCH3 */
        species[43] =
            -3.362764720000000e+04 * invT
            -2.237463319000000e+01
            -1.240373610000000e+00 * tc[0]
            -2.712470235000000e-02 * tc[1]
            +5.247088866666667e-06 * tc[2]
            -7.620769341666667e-10 * tc[3]
            +5.383278300000000e-14 * tc[4];
    } else {
        /*species 41: C2H5COCH2 */
        species[41] =
            -1.412849510000000e+04 * invT
            +6.292316490000000e+01
            -1.420987380000000e+01 * tc[0]
            -7.893322950000000e-03 * tc[1]
            +9.175486383333333e-07 * tc[2]
            -7.215587833333333e-11 * tc[3]
            +2.534566645000000e-15 * tc[4];
        /*species 43: NC3H7COCH3 */
        species[43] =
            -3.933754740000000e+04 * invT
            +7.448221789999999e+01
            -1.618415640000000e+01 * tc[0]
            -1.176229405000000e-02 * tc[1]
            +1.349458425000000e-06 * tc[2]
            -1.051922416666667e-10 * tc[3]
            +3.672439195000000e-15 * tc[4];
    }

    /*species with midpoint at T=1389 kelvin */
    if (T < 1389) {
        /*species 44: NC3H7COCH2 */
        species[44] =
            -1.230622380000000e+04 * invT
            -2.245298912900000e+01
            -9.582992710000000e-01 * tc[0]
            -2.840812660000000e-02 * tc[1]
            +6.651879683333334e-06 * tc[2]
            -1.272262616666667e-09 * tc[3]
            +1.246105235000000e-13 * tc[4];
    } else {
        /*species 44: NC3H7COCH2 */
        species[44] =
            -1.769685610000000e+04 * invT
            +7.453673030000000e+01
            -1.615024190000000e+01 * tc[0]
            -1.070467330000000e-02 * tc[1]
            +1.226766023333333e-06 * tc[2]
            -9.554734666666667e-11 * tc[3]
            +3.333563610000000e-15 * tc[4];
    }

    /*species with midpoint at T=1399 kelvin */
    if (T < 1399) {
        /*species 49: C7H14O */
        species[49] =
            -2.671182950000000e+04 * invT
            -6.056790396000000e+01
            +5.394772060000000e+00 * tc[0]
            -5.092328850000000e-02 * tc[1]
            +1.266743726666667e-05 * tc[2]
            -2.471150166666667e-09 * tc[3]
            +2.374271850000000e-13 * tc[4];
    } else {
        /*species 49: C7H14O */
        species[49] =
            -3.643496640000000e+04 * invT
            +1.212150009000000e+02
            -2.326927310000000e+01 * tc[0]
            -1.627926955000000e-02 * tc[1]
            +1.860415100000000e-06 * tc[2]
            -1.446447150000000e-10 * tc[3]
            +5.040660750000000e-15 * tc[4];
    }

    /*species with midpoint at T=1387 kelvin */
    if (T < 1387) {
        /*species 50: NC7KET */
        species[50] =
            -5.238395420000000e+04 * invT
            -2.785576648000000e+01
            -1.529366920000000e+00 * tc[0]
            -4.790867330000000e-02 * tc[1]
            +1.161147533333333e-05 * tc[2]
            -2.246169850000000e-09 * tc[3]
            +2.193640630000000e-13 * tc[4];
    } else {
        /*species 50: NC7KET */
        species[50] =
            -6.201424740000000e+04 * invT
            +1.463327428000000e+02
            -2.883303880000000e+01 * tc[0]
            -1.600840480000000e-02 * tc[1]
            +1.858474266666667e-06 * tc[2]
            -1.460217991666667e-10 * tc[3]
            +5.126022550000000e-15 * tc[4];
    }
    return;
}


/*compute the a/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void helmholtz(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];
    amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.54747466e+04 * invT
            +1.96638554e+00
            -2.50104422e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -1.01252100e+03 * invT
            +5.59221800e+00
            -3.29812400e+00 * tc[0]
            -4.12472100e-04 * tc[1]
            +1.35716917e-07 * tc[2]
            +7.89619500e-12 * tc[3]
            -2.06743600e-14 * tc[4];
        /*species 2: O */
        species[2] =
            +2.91476400e+04 * invT
            -1.01756600e+00
            -2.94642900e+00 * tc[0]
            +8.19083000e-04 * tc[1]
            -4.03505333e-07 * tc[2]
            +1.33570250e-10 * tc[3]
            -1.94534800e-14 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.00524900e+03 * invT
            -3.82180200e+00
            -3.21293600e+00 * tc[0]
            -5.63743000e-04 * tc[1]
            +9.59358333e-08 * tc[2]
            -1.09489750e-10 * tc[3]
            +4.38427700e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.02081100e+04 * invT
            -2.03391000e-01
            -3.38684200e+00 * tc[0]
            -1.73749100e-03 * tc[1]
            +1.05911600e-06 * tc[2]
            -5.80715083e-10 * tc[3]
            +1.25329400e-13 * tc[4];
        /*species 6: CO */
        species[6] =
            -1.43105400e+04 * invT
            -2.58644500e+00
            -3.26245200e+00 * tc[0]
            -7.55970500e-04 * tc[1]
            +6.46959167e-07 * tc[2]
            -4.65162000e-10 * tc[3]
            +1.23747550e-13 * tc[4];
        /*species 7: CO2 */
        species[7] =
            -4.83731400e+04 * invT
            -8.91276500e+00
            -2.27572500e+00 * tc[0]
            -4.96103600e-03 * tc[1]
            +1.73485167e-06 * tc[2]
            -5.72223917e-10 * tc[3]
            +1.05864000e-13 * tc[4];
        /*species 8: CH3 */
        species[8] =
            +1.64237800e+04 * invT
            -5.35935100e+00
            -2.43044300e+00 * tc[0]
            -5.56205000e-03 * tc[1]
            +2.80036667e-06 * tc[2]
            -1.35152417e-09 * tc[3]
            +2.93247650e-13 * tc[4];
        /*species 9: CH4 */
        species[9] =
            -9.82522900e+03 * invT
            -1.39434485e+01
            -7.78741500e-01 * tc[0]
            -8.73834000e-03 * tc[1]
            +4.63901500e-06 * tc[2]
            -2.54142333e-09 * tc[3]
            +6.11965500e-13 * tc[4];
        /*species 11: H2O2 */
        species[11] =
            -1.76631500e+04 * invT
            -4.39660900e+00
            -3.38875400e+00 * tc[0]
            -3.28461300e-03 * tc[1]
            +2.47502167e-08 * tc[2]
            +3.85483833e-10 * tc[3]
            -1.23575750e-13 * tc[4];
        /*species 12: CH2O */
        species[12] =
            -1.48654000e+04 * invT
            -1.31320890e+01
            -1.65273100e+00 * tc[0]
            -6.31572000e-03 * tc[1]
            +3.14694667e-06 * tc[2]
            -1.70835917e-09 * tc[3]
            +4.20661850e-13 * tc[4];
        /*species 14: C2H4 */
        species[14] =
            +5.57304600e+03 * invT
            -2.60729780e+01
            +8.61488000e-01 * tc[0]
            -1.39808150e-02 * tc[1]
            +5.64779500e-06 * tc[2]
            -2.32096000e-09 * tc[3]
            +4.86893950e-13 * tc[4];
        /*species 15: C2H5 */
        species[15] =
            +1.28704000e+04 * invT
            -1.04474980e+01
            -2.69070200e+00 * tc[0]
            -4.35956650e-03 * tc[1]
            -7.36639833e-07 * tc[2]
            -7.78225250e-11 * tc[3]
            +1.96388650e-13 * tc[4];
        /*species 16: C2H */
        species[16] =
            +6.68381300e+04 * invT
            -5.56251600e+00
            -2.73770400e+00 * tc[0]
            -4.02422300e-03 * tc[1]
            +1.54071833e-06 * tc[2]
            -5.43771583e-10 * tc[3]
            +9.69790000e-14 * tc[4];
        /*species 17: C2H2 */
        species[17] =
            +2.61244400e+04 * invT
            -7.79181600e+00
            -2.01356200e+00 * tc[0]
            -7.59522500e-03 * tc[1]
            +2.69386500e-06 * tc[2]
            -7.56582667e-10 * tc[3]
            +9.56373000e-14 * tc[4];
        /*species 18: CH3OH */
        species[18] =
            -2.53534800e+04 * invT
            -9.57251500e+00
            -2.66011500e+00 * tc[0]
            -3.67075400e-03 * tc[1]
            -1.19500850e-06 * tc[2]
            +7.32766167e-10 * tc[3]
            -1.19528500e-13 * tc[4];
        /*species 19: CH2CO */
        species[19] =
            -7.63263700e+03 * invT
            -6.69858200e+00
            -2.97497100e+00 * tc[0]
            -6.05935500e-03 * tc[1]
            +3.90841000e-07 * tc[2]
            +5.38890417e-10 * tc[3]
            -1.95282450e-13 * tc[4];
        /*species 20: HCCO */
        species[20] =
            +1.96589200e+04 * invT
            +3.56612110e+00
            -5.04796500e+00 * tc[0]
            -2.22673900e-03 * tc[1]
            -3.78047167e-08 * tc[2]
            +1.23507917e-10 * tc[3]
            -1.12537100e-14 * tc[4];
        /*species 21: CH3CHO */
        species[21] =
            -2.12458900e+04 * invT
            -1.18451950e+01
            -2.50569500e+00 * tc[0]
            -6.68495500e-03 * tc[1]
            -7.78658833e-07 * tc[2]
            +9.40116667e-10 * tc[3]
            -2.13178300e-13 * tc[4];
        /*species 39: C3H3 */
        species[39] =
            +3.98888300e+04 * invT
            +3.16874510e+00
            -4.75420000e+00 * tc[0]
            -5.54014000e-03 * tc[1]
            -4.65553833e-08 * tc[2]
            +4.56601000e-10 * tc[3]
            -9.74814500e-14 * tc[4];
        /*species 51: N2 */
        species[51] =
            -1.02090000e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120000e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242750e-13 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.54747466e+04 * invT
            +1.96638554e+00
            -2.50104422e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -8.35034000e+02 * invT
            +3.34653300e+00
            -2.99142300e+00 * tc[0]
            -3.50032200e-04 * tc[1]
            +9.38971500e-09 * tc[2]
            +7.69298167e-13 * tc[3]
            -7.91376000e-17 * tc[4];
        /*species 2: O */
        species[2] =
            +2.92308000e+04 * invT
            -3.37824800e+00
            -2.54206000e+00 * tc[0]
            +1.37753100e-05 * tc[1]
            +5.17133833e-10 * tc[2]
            -3.79255583e-13 * tc[3]
            +2.18402600e-17 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.23393000e+03 * invT
            -4.91588000e-01
            -3.69757800e+00 * tc[0]
            -3.06759850e-04 * tc[1]
            +2.09807000e-08 * tc[2]
            -1.47940083e-12 * tc[3]
            +5.68217500e-17 * tc[4];
        /*species 5: H2O */
        species[5] =
            -2.98992100e+04 * invT
            -5.19067100e+00
            -2.67214600e+00 * tc[0]
            -1.52814650e-03 * tc[1]
            +1.45504333e-07 * tc[2]
            -1.00083000e-11 * tc[3]
            +3.19580900e-16 * tc[4];
        /*species 6: CO */
        species[6] =
            -1.42683500e+04 * invT
            -4.08314000e+00
            -3.02507800e+00 * tc[0]
            -7.21344500e-04 * tc[1]
            +9.38471333e-08 * tc[2]
            -8.48817500e-12 * tc[3]
            +3.45547600e-16 * tc[4];
        /*species 7: CO2 */
        species[7] =
            -4.89669600e+04 * invT
            +4.40901890e+00
            -4.45362300e+00 * tc[0]
            -1.57008450e-03 * tc[1]
            +2.13068500e-07 * tc[2]
            -1.99499750e-11 * tc[3]
            +8.34516500e-16 * tc[4];
        /*species 8: CH3 */
        species[8] =
            +1.64378100e+04 * invT
            -3.60864500e+00
            -2.84405200e+00 * tc[0]
            -3.06898700e-03 * tc[1]
            +3.71724167e-07 * tc[2]
            -3.15430083e-11 * tc[3]
            +1.22607950e-15 * tc[4];
        /*species 9: CH4 */
        species[9] =
            -1.00807900e+04 * invT
            -8.93991600e+00
            -1.68347900e+00 * tc[0]
            -5.11862000e-03 * tc[1]
            +6.45854833e-07 * tc[2]
            -5.65465417e-11 * tc[3]
            +2.25171150e-15 * tc[4];
        /*species 11: H2O2 */
        species[11] =
            -1.80069600e+04 * invT
            +3.07203000e+00
            -4.57316700e+00 * tc[0]
            -2.16806800e-03 * tc[1]
            +2.45781500e-07 * tc[2]
            -1.95742000e-11 * tc[3]
            +7.15827000e-16 * tc[4];
        /*species 12: CH2O */
        species[12] =
            -1.53203700e+04 * invT
            -4.91696600e+00
            -2.99560600e+00 * tc[0]
            -3.34066050e-03 * tc[1]
            +4.38159167e-07 * tc[2]
            -3.94762750e-11 * tc[3]
            +1.60625850e-15 * tc[4];
        /*species 14: C2H4 */
        species[14] =
            +4.42828900e+03 * invT
            +2.98030000e-01
            -3.52841900e+00 * tc[0]
            -5.74259000e-03 * tc[1]
            +7.36397500e-07 * tc[2]
            -6.53716750e-11 * tc[3]
            +2.63342400e-15 * tc[4];
        /*species 15: C2H5 */
        species[15] =
            +1.06745500e+04 * invT
            +2.09713700e+01
            -7.19048000e+00 * tc[0]
            -3.24203850e-03 * tc[1]
            +1.07134417e-07 * tc[2]
            +1.95656583e-11 * tc[3]
            -1.94043850e-15 * tc[4];
        /*species 16: C2H */
        species[16] =
            +6.65588400e+04 * invT
            +1.79530400e+00
            -3.98636700e+00 * tc[0]
            -1.57156150e-03 * tc[1]
            +2.11207167e-07 * tc[2]
            -2.43696917e-11 * tc[3]
            +1.35816000e-15 * tc[4];
        /*species 17: C2H2 */
        species[17] =
            +2.56676600e+04 * invT
            +6.23710800e+00
            -4.43677000e+00 * tc[0]
            -2.68801950e-03 * tc[1]
            +3.18802833e-07 * tc[2]
            -2.73864917e-11 * tc[3]
            +1.07835500e-15 * tc[4];
        /*species 18: CH3OH */
        species[18] =
            -2.61579100e+04 * invT
            +6.50865000e-01
            -4.02906100e+00 * tc[0]
            -4.68829650e-03 * tc[1]
            +5.08375667e-07 * tc[2]
            -3.63232750e-11 * tc[3]
            +1.11236150e-15 * tc[4];
        /*species 19: CH2CO */
        species[19] =
            -8.58340200e+03 * invT
            +1.26963980e+01
            -6.03881700e+00 * tc[0]
            -2.90242000e-03 * tc[1]
            +3.20159000e-07 * tc[2]
            -2.32873750e-11 * tc[3]
            +7.29434000e-16 * tc[4];
        /*species 20: HCCO */
        species[20] =
            +1.90151300e+04 * invT
            +1.48293350e+01
            -6.75807300e+00 * tc[0]
            -1.00020000e-03 * tc[1]
            +3.37934500e-08 * tc[2]
            +8.67610000e-12 * tc[3]
            -9.82582500e-16 * tc[4];
        /*species 21: CH3CHO */
        species[21] =
            -2.26456900e+04 * invT
            +1.08815960e+01
            -5.86865000e+00 * tc[0]
            -5.39712000e-03 * tc[1]
            +6.07588333e-07 * tc[2]
            -4.51076000e-11 * tc[3]
            +1.44842200e-15 * tc[4];
        /*species 39: C3H3 */
        species[39] =
            +3.84742000e+04 * invT
            +2.96102370e+01
            -8.83104700e+00 * tc[0]
            -2.17859750e-03 * tc[1]
            +6.84844500e-08 * tc[2]
            +1.97393583e-11 * tc[3]
            -2.18826000e-15 * tc[4];
        /*species 51: N2 */
        species[51] =
            -9.22797700e+02 * invT
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988500e-04 * tc[1]
            +9.47460167e-08 * tc[2]
            -8.41420000e-12 * tc[3]
            +3.37667550e-16 * tc[4];
    }

    /*species with midpoint at T=1357 kelvin */
    if (T < 1357) {
        /*species 4: OH */
        species[4] =
            +3.74321252e+03 * invT
            -1.42790800e-02
            -3.43586219e+00 * tc[0]
            -1.01117902e-04 * tc[1]
            +1.89244020e-08 * tc[2]
            -2.02037624e-11 * tc[3]
            +3.71825515e-15 * tc[4];
    } else {
        /*species 4: OH */
        species[4] =
            +4.12085374e+03 * invT
            -5.48067553e+00
            -2.62599754e+00 * tc[0]
            -6.59962030e-04 * tc[1]
            +5.99541117e-08 * tc[2]
            -3.54692333e-12 * tc[3]
            +9.10240080e-17 * tc[4];
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 10: HO2 */
        species[10] =
            +8.09181013e+02 * invT
            -6.21060443e+00
            -3.18310656e+00 * tc[0]
            -1.83383975e-03 * tc[1]
            +1.55397520e-07 * tc[2]
            +2.71544099e-11 * tc[3]
            -7.55699560e-15 * tc[4];
        /*species 35: CH3O2H */
        species[35] =
            -1.77197926e+04 * invT
            -7.02181132e+00
            -3.23442817e+00 * tc[0]
            -9.50648835e-03 * tc[1]
            +1.88977145e-06 * tc[2]
            -2.83588878e-10 * tc[3]
            +2.05915111e-14 * tc[4];
        /*species 46: C7H15-1 */
        species[46] =
            -1.04590223e+03 * invT
            -3.61559715e+01
            +4.99570406e-01 * tc[0]
            -4.04413234e-02 * tc[1]
            +8.34221257e-06 * tc[2]
            -1.30457757e-09 * tc[3]
            +9.83081135e-14 * tc[4];
        /*species 47: C7H15O2 */
        species[47] =
            -2.15795254e+04 * invT
            -3.35880906e+01
            -3.55252917e-01 * tc[0]
            -4.71190503e-02 * tc[1]
            +1.11159184e-05 * tc[2]
            -2.12329653e-09 * tc[3]
            +2.06605552e-13 * tc[4];
        /*species 48: C7H14OOHO2 */
        species[48] =
            -3.46304788e+04 * invT
            -3.12132135e+01
            -1.40912479e+00 * tc[0]
            -5.53604170e-02 * tc[1]
            +1.45196346e-05 * tc[2]
            -3.03052291e-09 * tc[3]
            +3.14175344e-13 * tc[4];
    } else {
        /*species 10: HO2 */
        species[10] =
            +3.98127689e+02 * invT
            -1.96841300e-02
            -4.10547423e+00 * tc[0]
            -1.19226417e-03 * tc[1]
            +1.34391332e-07 * tc[2]
            -1.03493103e-11 * tc[3]
            +3.58200054e-16 * tc[4];
        /*species 35: CH3O2H */
        species[35] =
            -1.96678771e+04 * invT
            +2.65482338e+01
            -8.43117091e+00 * tc[0]
            -4.03408955e-03 * tc[1]
            +4.61824868e-07 * tc[2]
            -3.59443536e-11 * tc[3]
            +1.25346073e-15 * tc[4];
        /*species 46: C7H15-1 */
        species[46] =
            -9.20938221e+03 * invT
            +1.07289502e+02
            -2.17940709e+01 * tc[0]
            -1.63140122e-02 * tc[1]
            +1.85230407e-06 * tc[2]
            -1.43389290e-10 * tc[3]
            +4.98183499e-15 * tc[4];
        /*species 47: C7H15O2 */
        species[47] =
            -3.06842095e+04 * invT
            +1.30461127e+02
            -2.60535640e+01 * tc[0]
            -1.71915844e-02 * tc[1]
            +1.97855412e-06 * tc[2]
            -1.54541867e-10 * tc[3]
            +5.40261995e-15 * tc[4];
        /*species 48: C7H14OOHO2 */
        species[48] =
            -4.52873003e+04 * invT
            +1.67811540e+02
            -3.28358296e+01 * tc[0]
            -1.66383642e-02 * tc[1]
            +1.93055475e-06 * tc[2]
            -1.51632198e-10 * tc[3]
            +5.32172095e-15 * tc[4];
    }

    /*species with midpoint at T=1384 kelvin */
    if (T < 1384) {
        /*species 13: C2H6 */
        species[13] =
            -1.12345534e+04 * invT
            -2.21901604e+01
            +2.52854344e-02 * tc[0]
            -1.20382377e-02 * tc[1]
            +1.86489120e-06 * tc[2]
            -1.73617417e-10 * tc[3]
            +2.64934308e-15 * tc[4];
    } else {
        /*species 13: C2H6 */
        species[13] =
            -1.37500014e+04 * invT
            +1.81149589e+01
            -6.10683385e+00 * tc[0]
            -6.46181805e-03 * tc[1]
            +7.37545327e-07 * tc[2]
            -5.72826438e-11 * tc[3]
            +1.99450866e-15 * tc[4];
    }

    /*species with midpoint at T=1400 kelvin */
    if (T < 1400) {
        /*species 22: C3H4-A */
        species[22] =
            +2.25124300e+04 * invT
            -8.39587100e+00
            -2.53983100e+00 * tc[0]
            -8.16718500e-03 * tc[1]
            +2.94158333e-07 * tc[2]
            +3.87280417e-10 * tc[3]
            -8.64565500e-14 * tc[4];
        /*species 23: C3H4-P */
        species[23] =
            +2.14840800e+04 * invT
            -5.97486400e+00
            -3.02973000e+00 * tc[0]
            -7.49480500e-03 * tc[1]
            +2.33083333e-07 * tc[2]
            +3.30801583e-10 * tc[3]
            -6.94108500e-14 * tc[4];
        /*species 42: CH3CHCO */
        species[42] =
            -1.15276540e+04 * invT
            -1.66714056e+01
            -1.48380119e+00 * tc[0]
            -1.61101506e-02 * tc[1]
            +4.50416722e-06 * tc[2]
            -1.00415970e-09 * tc[3]
            +1.09182966e-13 * tc[4];
    } else {
        /*species 22: C3H4-A */
        species[22] =
            +1.95497200e+04 * invT
            +3.95468660e+01
            -9.77625600e+00 * tc[0]
            -2.65106900e-03 * tc[1]
            +6.16853000e-08 * tc[2]
            +2.52198833e-11 * tc[3]
            -2.54479050e-15 * tc[4];
        /*species 23: C3H4-P */
        species[23] =
            +1.86027700e+04 * invT
            +3.89748820e+01
            -9.76810200e+00 * tc[0]
            -2.60957550e-03 * tc[1]
            +6.25523333e-08 * tc[2]
            +2.49349250e-11 * tc[3]
            -2.55393900e-15 * tc[4];
        /*species 42: CH3CHCO */
        species[42] =
            -1.42482738e+04 * invT
            +3.68049096e+01
            -1.00219123e+01 * tc[0]
            -4.78483150e-03 * tc[1]
            +5.43702740e-07 * tc[2]
            -4.21026422e-11 * tc[3]
            +1.46296628e-15 * tc[4];
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 24: C3H6 */
        species[24] =
            +1.06688164e+03 * invT
            -2.25057582e+01
            -3.94615444e-01 * tc[0]
            -1.44553831e-02 * tc[1]
            +2.58144680e-06 * tc[2]
            -3.24011841e-10 * tc[3]
            +1.68945176e-14 * tc[4];
    } else {
        /*species 24: C3H6 */
        species[24] =
            -1.87821271e+03 * invT
            +2.70320264e+01
            -8.01595958e+00 * tc[0]
            -6.85118170e-03 * tc[1]
            +7.77082888e-07 * tc[2]
            -6.01045335e-11 * tc[3]
            +2.08685063e-15 * tc[4];
    }

    /*species with midpoint at T=1398 kelvin */
    if (T < 1398) {
        /*species 25: C4H6 */
        species[25] =
            +1.17551314e+04 * invT
            -3.15135345e+01
            +1.43095121e+00 * tc[0]
            -2.39353031e-02 * tc[1]
            +6.92411333e-06 * tc[2]
            -1.59624627e-09 * tc[3]
            +1.78579253e-13 * tc[4];
    } else {
        /*species 25: C4H6 */
        species[25] =
            +7.79039770e+03 * invT
            +4.71481738e+01
            -1.11633789e+01 * tc[0]
            -6.85819825e-03 * tc[1]
            +7.82859638e-07 * tc[2]
            -6.08078197e-11 * tc[3]
            +2.11743101e-15 * tc[4];
    }

    /*species with midpoint at T=1392 kelvin */
    if (T < 1392) {
        /*species 26: C4H7 */
        species[26] =
            +1.49933591e+04 * invT
            -2.80800780e+01
            +3.50508352e-01 * tc[0]
            -2.13255622e-02 * tc[1]
            +4.84965622e-06 * tc[2]
            -8.78365950e-10 * tc[3]
            +8.00299270e-14 * tc[4];
        /*species 27: C4H8-1 */
        species[27] =
            -1.57875035e+03 * invT
            -3.13397957e+01
            +8.31372089e-01 * tc[0]
            -2.26290489e-02 * tc[1]
            +4.89430932e-06 * tc[2]
            -8.35170300e-10 * tc[3]
            +7.15958400e-14 * tc[4];
        /*species 31: C5H9 */
        species[31] =
            +1.25589824e+04 * invT
            -3.50242699e+01
            +1.38013950e+00 * tc[0]
            -2.78804243e-02 * tc[1]
            +6.16906547e-06 * tc[2]
            -1.05736584e-09 * tc[3]
            +8.92694175e-14 * tc[4];
        /*species 32: C5H10-1 */
        species[32] =
            -4.46546666e+03 * invT
            -3.43362138e+01
            +1.06223481e+00 * tc[0]
            -2.87109147e-02 * tc[1]
            +6.24144817e-06 * tc[2]
            -1.06137491e-09 * tc[3]
            +8.98048945e-14 * tc[4];
    } else {
        /*species 26: C4H7 */
        species[26] =
            +1.09041937e+04 * invT
            +4.57696593e+01
            -1.12103578e+01 * tc[0]
            -8.02415980e-03 * tc[1]
            +9.10837153e-07 * tc[2]
            -7.04950877e-11 * tc[3]
            +2.44886369e-15 * tc[4];
        /*species 27: C4H8-1 */
        species[27] =
            -5.97871038e+03 * invT
            +4.67878106e+01
            -1.13508668e+01 * tc[0]
            -9.03089385e-03 * tc[1]
            +1.02682171e-06 * tc[2]
            -7.95544133e-11 * tc[3]
            +2.76544820e-15 * tc[4];
        /*species 31: C5H9 */
        species[31] =
            +7.00496135e+03 * invT
            +6.46362227e+01
            -1.41860454e+01 * tc[0]
            -1.03564449e-02 * tc[1]
            +1.17826770e-06 * tc[2]
            -9.13392775e-11 * tc[3]
            +3.17661104e-15 * tc[4];
        /*species 32: C5H10-1 */
        species[32] =
            -1.00898205e+04 * invT
            +6.59535475e+01
            -1.45851539e+01 * tc[0]
            -1.12036235e-02 * tc[1]
            +1.27224671e-06 * tc[2]
            -9.84908050e-11 * tc[3]
            +3.42192569e-15 * tc[4];
    }

    /*species with midpoint at T=1395 kelvin */
    if (T < 1395) {
        /*species 28: PC4H9 */
        species[28] =
            +7.68945248e+03 * invT
            -3.01230529e+01
            +4.37779725e-01 * tc[0]
            -2.39486182e-02 * tc[1]
            +5.23371932e-06 * tc[2]
            -9.14887267e-10 * tc[3]
            +8.10053320e-14 * tc[4];
    } else {
        /*species 28: PC4H9 */
        species[28] =
            +3.17231942e+03 * invT
            +5.04935904e+01
            -1.21510082e+01 * tc[0]
            -9.71553585e-03 * tc[1]
            +1.10262992e-06 * tc[2]
            -8.53126133e-11 * tc[3]
            +2.96264853e-15 * tc[4];
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 29: CH3COCH2 */
        species[29] =
            -5.37233261e+03 * invT
            -1.70012489e+01
            -1.80339187e+00 * tc[0]
            -1.50703542e-02 * tc[1]
            +3.22509253e-06 * tc[2]
            -5.31832528e-10 * tc[3]
            +4.33051590e-14 * tc[4];
        /*species 45: NC7H16 */
        species[45] =
            -2.56586565e+04 * invT
            -3.76416531e+01
            +1.26836187e+00 * tc[0]
            -4.27177910e-02 * tc[1]
            +8.75577977e-06 * tc[2]
            -1.35788101e-09 * tc[3]
            +1.01197462e-13 * tc[4];
    } else {
        /*species 29: CH3COCH2 */
        species[29] =
            -8.44380462e+03 * invT
            +3.71501195e+01
            -1.02303975e+01 * tc[0]
            -5.82470805e-03 * tc[1]
            +6.68342562e-07 * tc[2]
            -5.21004372e-11 * tc[3]
            +1.81892181e-15 * tc[4];
        /*species 45: NC7H16 */
        species[45] =
            -3.42760081e+04 * invT
            +1.13518917e+02
            -2.22148969e+01 * tc[0]
            -1.73837875e-02 * tc[1]
            +1.97345215e-06 * tc[2]
            -1.52748732e-10 * tc[3]
            +5.30651330e-15 * tc[4];
    }

    /*species with midpoint at T=1378 kelvin */
    if (T < 1378) {
        /*species 30: C2H5CHO */
        species[30] =
            -2.42260137e+04 * invT
            -1.49522504e+01
            -2.16308444e+00 * tc[0]
            -1.47750632e-02 * tc[1]
            +2.54077087e-06 * tc[2]
            -2.91253289e-10 * tc[3]
            +1.19448313e-14 * tc[4];
        /*species 40: NC3H7CHO */
        species[40] =
            -2.71032194e+04 * invT
            -1.82826978e+01
            -1.87415959e+00 * tc[0]
            -2.09620158e-02 * tc[1]
            +3.91914632e-06 * tc[2]
            -5.22428061e-10 * tc[3]
            +3.04721954e-14 * tc[4];
    } else {
        /*species 30: C2H5CHO */
        species[30] =
            -2.74143135e+04 * invT
            +3.77774538e+01
            -1.02427695e+01 * tc[0]
            -6.98209945e-03 * tc[1]
            +7.93746668e-07 * tc[2]
            -6.15088088e-11 * tc[3]
            +2.13879752e-15 * tc[4];
        /*species 40: NC3H7CHO */
        species[40] =
            -3.15845348e+04 * invT
            +5.77778296e+01
            -1.35988068e+01 * tc[0]
            -9.08262370e-03 * tc[1]
            +1.02974076e-06 * tc[2]
            -7.96650173e-11 * tc[3]
            +2.76721479e-15 * tc[4];
    }

    /*species with midpoint at T=1397 kelvin */
    if (T < 1397) {
        /*species 33: C5H11-1 */
        species[33] =
            +4.83995303e+03 * invT
            -3.44627522e+01
            +9.05255912e-01 * tc[0]
            -3.05316426e-02 * tc[1]
            +6.82486375e-06 * tc[2]
            -1.21744558e-09 * tc[3]
            +1.09429807e-13 * tc[4];
        /*species 38: C3H5-A */
        species[38] =
            +1.93834226e+04 * invT
            -2.68358451e+01
            +5.29131958e-01 * tc[0]
            -1.67279550e-02 * tc[1]
            +4.22335045e-06 * tc[2]
            -8.57214617e-10 * tc[3]
            +8.66291700e-14 * tc[4];
    } else {
        /*species 33: C5H11-1 */
        species[33] =
            -9.23241637e+02 * invT
            +6.92763599e+01
            -1.53234740e+01 * tc[0]
            -1.19520600e-02 * tc[1]
            +1.35795270e-06 * tc[2]
            -1.05146863e-10 * tc[3]
            +3.65338668e-15 * tc[4];
        /*species 38: C3H5-A */
        species[38] =
            +1.63576092e+04 * invT
            +3.00397846e+01
            -8.45883958e+00 * tc[0]
            -5.63477415e-03 * tc[1]
            +6.39654773e-07 * tc[2]
            -4.95049266e-11 * tc[3]
            +1.71959015e-15 * tc[4];
    }

    /*species with midpoint at T=1385 kelvin */
    if (T < 1385) {
        /*species 34: CH3O2 */
        species[34] =
            -6.84394259e+02 * invT
            -1.90183414e+00
            -4.26146906e+00 * tc[0]
            -5.04367995e-03 * tc[1]
            +5.35843640e-07 * tc[2]
            -1.74507723e-11 * tc[3]
            -2.09169552e-15 * tc[4];
    } else {
        /*species 34: CH3O2 */
        species[34] =
            -1.53574838e+03 * invT
            +9.67751777e+00
            -5.95787891e+00 * tc[0]
            -3.95364313e-03 * tc[1]
            +4.47077057e-07 * tc[2]
            -3.44909447e-11 * tc[3]
            +1.19503665e-15 * tc[4];
    }

    /*species with midpoint at T=1402 kelvin */
    if (T < 1402) {
        /*species 36: C2H3CO */
        species[36] =
            +4.25770215e+03 * invT
            -1.69002345e+01
            -1.36242013e+00 * tc[0]
            -1.57636986e-02 * tc[1]
            +5.00364892e-06 * tc[2]
            -1.23472593e-09 * tc[3]
            +1.43985765e-13 * tc[4];
    } else {
        /*species 36: C2H3CO */
        species[36] =
            +1.92969514e+03 * invT
            +3.24639464e+01
            -9.37467676e+00 * tc[0]
            -3.95648450e-03 * tc[1]
            +4.45330467e-07 * tc[2]
            -3.42596192e-11 * tc[3]
            +1.18489490e-15 * tc[4];
    }

    /*species with midpoint at T=1393 kelvin */
    if (T < 1393) {
        /*species 37: C2H3CHO */
        species[37] =
            -1.16521584e+04 * invT
            -2.35954728e+01
            -2.92355162e-01 * tc[0]
            -1.77160709e-02 * tc[1]
            +4.91560540e-06 * tc[2]
            -1.06750103e-09 * tc[3]
            +1.13072054e-13 * tc[4];
    } else {
        /*species 37: C2H3CHO */
        species[37] =
            -1.49630281e+04 * invT
            +4.01420020e+01
            -1.04184959e+01 * tc[0]
            -4.74481661e-03 * tc[1]
            +5.48850882e-07 * tc[2]
            -4.30232669e-11 * tc[3]
            +1.50793646e-15 * tc[4];
    }

    /*species with midpoint at T=1383 kelvin */
    if (T < 1383) {
        /*species 41: C2H5COCH2 */
        species[41] =
            -9.50796505e+03 * invT
            -1.94305255e+01
            -1.54013856e+00 * tc[0]
            -2.19743129e-02 * tc[1]
            +4.95004035e-06 * tc[2]
            -8.79127608e-10 * tc[3]
            +7.92993845e-14 * tc[4];
        /*species 43: NC3H7COCH3 */
        species[43] =
            -3.36276472e+04 * invT
            -2.33746332e+01
            -1.24037361e+00 * tc[0]
            -2.71247023e-02 * tc[1]
            +5.24708887e-06 * tc[2]
            -7.62076934e-10 * tc[3]
            +5.38327830e-14 * tc[4];
    } else {
        /*species 41: C2H5COCH2 */
        species[41] =
            -1.41284951e+04 * invT
            +6.19231649e+01
            -1.42098738e+01 * tc[0]
            -7.89332295e-03 * tc[1]
            +9.17548638e-07 * tc[2]
            -7.21558783e-11 * tc[3]
            +2.53456664e-15 * tc[4];
        /*species 43: NC3H7COCH3 */
        species[43] =
            -3.93375474e+04 * invT
            +7.34822179e+01
            -1.61841564e+01 * tc[0]
            -1.17622941e-02 * tc[1]
            +1.34945843e-06 * tc[2]
            -1.05192242e-10 * tc[3]
            +3.67243920e-15 * tc[4];
    }

    /*species with midpoint at T=1389 kelvin */
    if (T < 1389) {
        /*species 44: NC3H7COCH2 */
        species[44] =
            -1.23062238e+04 * invT
            -2.34529891e+01
            -9.58299271e-01 * tc[0]
            -2.84081266e-02 * tc[1]
            +6.65187968e-06 * tc[2]
            -1.27226262e-09 * tc[3]
            +1.24610524e-13 * tc[4];
    } else {
        /*species 44: NC3H7COCH2 */
        species[44] =
            -1.76968561e+04 * invT
            +7.35367303e+01
            -1.61502419e+01 * tc[0]
            -1.07046733e-02 * tc[1]
            +1.22676602e-06 * tc[2]
            -9.55473467e-11 * tc[3]
            +3.33356361e-15 * tc[4];
    }

    /*species with midpoint at T=1399 kelvin */
    if (T < 1399) {
        /*species 49: C7H14O */
        species[49] =
            -2.67118295e+04 * invT
            -6.15679040e+01
            +5.39477206e+00 * tc[0]
            -5.09232885e-02 * tc[1]
            +1.26674373e-05 * tc[2]
            -2.47115017e-09 * tc[3]
            +2.37427185e-13 * tc[4];
    } else {
        /*species 49: C7H14O */
        species[49] =
            -3.64349664e+04 * invT
            +1.20215001e+02
            -2.32692731e+01 * tc[0]
            -1.62792695e-02 * tc[1]
            +1.86041510e-06 * tc[2]
            -1.44644715e-10 * tc[3]
            +5.04066075e-15 * tc[4];
    }

    /*species with midpoint at T=1387 kelvin */
    if (T < 1387) {
        /*species 50: NC7KET */
        species[50] =
            -5.23839542e+04 * invT
            -2.88557665e+01
            -1.52936692e+00 * tc[0]
            -4.79086733e-02 * tc[1]
            +1.16114753e-05 * tc[2]
            -2.24616985e-09 * tc[3]
            +2.19364063e-13 * tc[4];
    } else {
        /*species 50: NC7KET */
        species[50] =
            -6.20142474e+04 * invT
            +1.45332743e+02
            -2.88330388e+01 * tc[0]
            -1.60084048e-02 * tc[1]
            +1.85847427e-06 * tc[2]
            -1.46021799e-10 * tc[3]
            +5.12602255e-15 * tc[4];
    }
    return;
}


/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void speciesInternalEnergy(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];
    amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +1.50104422e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54747466e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +2.29812400e+00
            +4.12472100e-04 * tc[1]
            -2.71433833e-07 * tc[2]
            -2.36885850e-11 * tc[3]
            +8.26974400e-14 * tc[4]
            -1.01252100e+03 * invT;
        /*species 2: O */
        species[2] =
            +1.94642900e+00
            -8.19083000e-04 * tc[1]
            +8.07010667e-07 * tc[2]
            -4.00710750e-10 * tc[3]
            +7.78139200e-14 * tc[4]
            +2.91476400e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +2.21293600e+00
            +5.63743000e-04 * tc[1]
            -1.91871667e-07 * tc[2]
            +3.28469250e-10 * tc[3]
            -1.75371080e-13 * tc[4]
            -1.00524900e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +2.38684200e+00
            +1.73749100e-03 * tc[1]
            -2.11823200e-06 * tc[2]
            +1.74214525e-09 * tc[3]
            -5.01317600e-13 * tc[4]
            -3.02081100e+04 * invT;
        /*species 6: CO */
        species[6] =
            +2.26245200e+00
            +7.55970500e-04 * tc[1]
            -1.29391833e-06 * tc[2]
            +1.39548600e-09 * tc[3]
            -4.94990200e-13 * tc[4]
            -1.43105400e+04 * invT;
        /*species 7: CO2 */
        species[7] =
            +1.27572500e+00
            +4.96103600e-03 * tc[1]
            -3.46970333e-06 * tc[2]
            +1.71667175e-09 * tc[3]
            -4.23456000e-13 * tc[4]
            -4.83731400e+04 * invT;
        /*species 8: CH3 */
        species[8] =
            +1.43044300e+00
            +5.56205000e-03 * tc[1]
            -5.60073333e-06 * tc[2]
            +4.05457250e-09 * tc[3]
            -1.17299060e-12 * tc[4]
            +1.64237800e+04 * invT;
        /*species 9: CH4 */
        species[9] =
            -2.21258500e-01
            +8.73834000e-03 * tc[1]
            -9.27803000e-06 * tc[2]
            +7.62427000e-09 * tc[3]
            -2.44786200e-12 * tc[4]
            -9.82522900e+03 * invT;
        /*species 11: H2O2 */
        species[11] =
            +2.38875400e+00
            +3.28461300e-03 * tc[1]
            -4.95004333e-08 * tc[2]
            -1.15645150e-09 * tc[3]
            +4.94303000e-13 * tc[4]
            -1.76631500e+04 * invT;
        /*species 12: CH2O */
        species[12] =
            +6.52731000e-01
            +6.31572000e-03 * tc[1]
            -6.29389333e-06 * tc[2]
            +5.12507750e-09 * tc[3]
            -1.68264740e-12 * tc[4]
            -1.48654000e+04 * invT;
        /*species 14: C2H4 */
        species[14] =
            -1.86148800e+00
            +1.39808150e-02 * tc[1]
            -1.12955900e-05 * tc[2]
            +6.96288000e-09 * tc[3]
            -1.94757580e-12 * tc[4]
            +5.57304600e+03 * invT;
        /*species 15: C2H5 */
        species[15] =
            +1.69070200e+00
            +4.35956650e-03 * tc[1]
            +1.47327967e-06 * tc[2]
            +2.33467575e-10 * tc[3]
            -7.85554600e-13 * tc[4]
            +1.28704000e+04 * invT;
        /*species 16: C2H */
        species[16] =
            +1.73770400e+00
            +4.02422300e-03 * tc[1]
            -3.08143667e-06 * tc[2]
            +1.63131475e-09 * tc[3]
            -3.87916000e-13 * tc[4]
            +6.68381300e+04 * invT;
        /*species 17: C2H2 */
        species[17] =
            +1.01356200e+00
            +7.59522500e-03 * tc[1]
            -5.38773000e-06 * tc[2]
            +2.26974800e-09 * tc[3]
            -3.82549200e-13 * tc[4]
            +2.61244400e+04 * invT;
        /*species 18: CH3OH */
        species[18] =
            +1.66011500e+00
            +3.67075400e-03 * tc[1]
            +2.39001700e-06 * tc[2]
            -2.19829850e-09 * tc[3]
            +4.78114000e-13 * tc[4]
            -2.53534800e+04 * invT;
        /*species 19: CH2CO */
        species[19] =
            +1.97497100e+00
            +6.05935500e-03 * tc[1]
            -7.81682000e-07 * tc[2]
            -1.61667125e-09 * tc[3]
            +7.81129800e-13 * tc[4]
            -7.63263700e+03 * invT;
        /*species 20: HCCO */
        species[20] =
            +4.04796500e+00
            +2.22673900e-03 * tc[1]
            +7.56094333e-08 * tc[2]
            -3.70523750e-10 * tc[3]
            +4.50148400e-14 * tc[4]
            +1.96589200e+04 * invT;
        /*species 21: CH3CHO */
        species[21] =
            +1.50569500e+00
            +6.68495500e-03 * tc[1]
            +1.55731767e-06 * tc[2]
            -2.82035000e-09 * tc[3]
            +8.52713200e-13 * tc[4]
            -2.12458900e+04 * invT;
        /*species 39: C3H3 */
        species[39] =
            +3.75420000e+00
            +5.54014000e-03 * tc[1]
            +9.31107667e-08 * tc[2]
            -1.36980300e-09 * tc[3]
            +3.89925800e-13 * tc[4]
            +3.98888300e+04 * invT;
        /*species 51: N2 */
        species[51] =
            +2.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
    } else {
        /*species 0: H */
        species[0] =
            +1.50104422e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54747466e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +1.99142300e+00
            +3.50032200e-04 * tc[1]
            -1.87794300e-08 * tc[2]
            -2.30789450e-12 * tc[3]
            +3.16550400e-16 * tc[4]
            -8.35034000e+02 * invT;
        /*species 2: O */
        species[2] =
            +1.54206000e+00
            -1.37753100e-05 * tc[1]
            -1.03426767e-09 * tc[2]
            +1.13776675e-12 * tc[3]
            -8.73610400e-17 * tc[4]
            +2.92308000e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +2.69757800e+00
            +3.06759850e-04 * tc[1]
            -4.19614000e-08 * tc[2]
            +4.43820250e-12 * tc[3]
            -2.27287000e-16 * tc[4]
            -1.23393000e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +1.67214600e+00
            +1.52814650e-03 * tc[1]
            -2.91008667e-07 * tc[2]
            +3.00249000e-11 * tc[3]
            -1.27832360e-15 * tc[4]
            -2.98992100e+04 * invT;
        /*species 6: CO */
        species[6] =
            +2.02507800e+00
            +7.21344500e-04 * tc[1]
            -1.87694267e-07 * tc[2]
            +2.54645250e-11 * tc[3]
            -1.38219040e-15 * tc[4]
            -1.42683500e+04 * invT;
        /*species 7: CO2 */
        species[7] =
            +3.45362300e+00
            +1.57008450e-03 * tc[1]
            -4.26137000e-07 * tc[2]
            +5.98499250e-11 * tc[3]
            -3.33806600e-15 * tc[4]
            -4.89669600e+04 * invT;
        /*species 8: CH3 */
        species[8] =
            +1.84405200e+00
            +3.06898700e-03 * tc[1]
            -7.43448333e-07 * tc[2]
            +9.46290250e-11 * tc[3]
            -4.90431800e-15 * tc[4]
            +1.64378100e+04 * invT;
        /*species 9: CH4 */
        species[9] =
            +6.83479000e-01
            +5.11862000e-03 * tc[1]
            -1.29170967e-06 * tc[2]
            +1.69639625e-10 * tc[3]
            -9.00684600e-15 * tc[4]
            -1.00807900e+04 * invT;
        /*species 11: H2O2 */
        species[11] =
            +3.57316700e+00
            +2.16806800e-03 * tc[1]
            -4.91563000e-07 * tc[2]
            +5.87226000e-11 * tc[3]
            -2.86330800e-15 * tc[4]
            -1.80069600e+04 * invT;
        /*species 12: CH2O */
        species[12] =
            +1.99560600e+00
            +3.34066050e-03 * tc[1]
            -8.76318333e-07 * tc[2]
            +1.18428825e-10 * tc[3]
            -6.42503400e-15 * tc[4]
            -1.53203700e+04 * invT;
        /*species 14: C2H4 */
        species[14] =
            +2.52841900e+00
            +5.74259000e-03 * tc[1]
            -1.47279500e-06 * tc[2]
            +1.96115025e-10 * tc[3]
            -1.05336960e-14 * tc[4]
            +4.42828900e+03 * invT;
        /*species 15: C2H5 */
        species[15] =
            +6.19048000e+00
            +3.24203850e-03 * tc[1]
            -2.14268833e-07 * tc[2]
            -5.86969750e-11 * tc[3]
            +7.76175400e-15 * tc[4]
            +1.06745500e+04 * invT;
        /*species 16: C2H */
        species[16] =
            +2.98636700e+00
            +1.57156150e-03 * tc[1]
            -4.22414333e-07 * tc[2]
            +7.31090750e-11 * tc[3]
            -5.43264000e-15 * tc[4]
            +6.65588400e+04 * invT;
        /*species 17: C2H2 */
        species[17] =
            +3.43677000e+00
            +2.68801950e-03 * tc[1]
            -6.37605667e-07 * tc[2]
            +8.21594750e-11 * tc[3]
            -4.31342000e-15 * tc[4]
            +2.56676600e+04 * invT;
        /*species 18: CH3OH */
        species[18] =
            +3.02906100e+00
            +4.68829650e-03 * tc[1]
            -1.01675133e-06 * tc[2]
            +1.08969825e-10 * tc[3]
            -4.44944600e-15 * tc[4]
            -2.61579100e+04 * invT;
        /*species 19: CH2CO */
        species[19] =
            +5.03881700e+00
            +2.90242000e-03 * tc[1]
            -6.40318000e-07 * tc[2]
            +6.98621250e-11 * tc[3]
            -2.91773600e-15 * tc[4]
            -8.58340200e+03 * invT;
        /*species 20: HCCO */
        species[20] =
            +5.75807300e+00
            +1.00020000e-03 * tc[1]
            -6.75869000e-08 * tc[2]
            -2.60283000e-11 * tc[3]
            +3.93033000e-15 * tc[4]
            +1.90151300e+04 * invT;
        /*species 21: CH3CHO */
        species[21] =
            +4.86865000e+00
            +5.39712000e-03 * tc[1]
            -1.21517667e-06 * tc[2]
            +1.35322800e-10 * tc[3]
            -5.79368800e-15 * tc[4]
            -2.26456900e+04 * invT;
        /*species 39: C3H3 */
        species[39] =
            +7.83104700e+00
            +2.17859750e-03 * tc[1]
            -1.36968900e-07 * tc[2]
            -5.92180750e-11 * tc[3]
            +8.75304000e-15 * tc[4]
            +3.84742000e+04 * invT;
        /*species 51: N2 */
        species[51] =
            +1.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }

    /*species with midpoint at T=1357 kelvin */
    if (T < 1357) {
        /*species 4: OH */
        species[4] =
            +2.43586219e+00
            +1.01117902e-04 * tc[1]
            -3.78488040e-08 * tc[2]
            +6.06112872e-11 * tc[3]
            -1.48730206e-14 * tc[4]
            +3.74321252e+03 * invT;
    } else {
        /*species 4: OH */
        species[4] =
            +1.62599754e+00
            +6.59962030e-04 * tc[1]
            -1.19908223e-07 * tc[2]
            +1.06407700e-11 * tc[3]
            -3.64096032e-16 * tc[4]
            +4.12085374e+03 * invT;
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 10: HO2 */
        species[10] =
            +2.18310656e+00
            +1.83383975e-03 * tc[1]
            -3.10795041e-07 * tc[2]
            -8.14632298e-11 * tc[3]
            +3.02279824e-14 * tc[4]
            +8.09181013e+02 * invT;
        /*species 35: CH3O2H */
        species[35] =
            +2.23442817e+00
            +9.50648835e-03 * tc[1]
            -3.77954290e-06 * tc[2]
            +8.50766632e-10 * tc[3]
            -8.23660444e-14 * tc[4]
            -1.77197926e+04 * invT;
        /*species 46: C7H15-1 */
        species[46] =
            -1.49957041e+00
            +4.04413234e-02 * tc[1]
            -1.66844251e-05 * tc[2]
            +3.91373270e-09 * tc[3]
            -3.93232454e-13 * tc[4]
            -1.04590223e+03 * invT;
        /*species 47: C7H15O2 */
        species[47] =
            -6.44747083e-01
            +4.71190503e-02 * tc[1]
            -2.22318369e-05 * tc[2]
            +6.36988958e-09 * tc[3]
            -8.26422210e-13 * tc[4]
            -2.15795254e+04 * invT;
        /*species 48: C7H14OOHO2 */
        species[48] =
            +4.09124790e-01
            +5.53604170e-02 * tc[1]
            -2.90392692e-05 * tc[2]
            +9.09156872e-09 * tc[3]
            -1.25670138e-12 * tc[4]
            -3.46304788e+04 * invT;
    } else {
        /*species 10: HO2 */
        species[10] =
            +3.10547423e+00
            +1.19226417e-03 * tc[1]
            -2.68782663e-07 * tc[2]
            +3.10479308e-11 * tc[3]
            -1.43280022e-15 * tc[4]
            +3.98127689e+02 * invT;
        /*species 35: CH3O2H */
        species[35] =
            +7.43117091e+00
            +4.03408955e-03 * tc[1]
            -9.23649737e-07 * tc[2]
            +1.07833061e-10 * tc[3]
            -5.01384292e-15 * tc[4]
            -1.96678771e+04 * invT;
        /*species 46: C7H15-1 */
        species[46] =
            +2.07940709e+01
            +1.63140122e-02 * tc[1]
            -3.70460813e-06 * tc[2]
            +4.30167870e-10 * tc[3]
            -1.99273400e-14 * tc[4]
            -9.20938221e+03 * invT;
        /*species 47: C7H15O2 */
        species[47] =
            +2.50535640e+01
            +1.71915844e-02 * tc[1]
            -3.95710823e-06 * tc[2]
            +4.63625602e-10 * tc[3]
            -2.16104798e-14 * tc[4]
            -3.06842095e+04 * invT;
        /*species 48: C7H14OOHO2 */
        species[48] =
            +3.18358296e+01
            +1.66383642e-02 * tc[1]
            -3.86110950e-06 * tc[2]
            +4.54896595e-10 * tc[3]
            -2.12868838e-14 * tc[4]
            -4.52873003e+04 * invT;
    }

    /*species with midpoint at T=1384 kelvin */
    if (T < 1384) {
        /*species 13: C2H6 */
        species[13] =
            -1.02528543e+00
            +1.20382377e-02 * tc[1]
            -3.72978240e-06 * tc[2]
            +5.20852252e-10 * tc[3]
            -1.05973723e-14 * tc[4]
            -1.12345534e+04 * invT;
    } else {
        /*species 13: C2H6 */
        species[13] =
            +5.10683385e+00
            +6.46181805e-03 * tc[1]
            -1.47509065e-06 * tc[2]
            +1.71847932e-10 * tc[3]
            -7.97803464e-15 * tc[4]
            -1.37500014e+04 * invT;
    }

    /*species with midpoint at T=1400 kelvin */
    if (T < 1400) {
        /*species 22: C3H4-A */
        species[22] =
            +1.53983100e+00
            +8.16718500e-03 * tc[1]
            -5.88316667e-07 * tc[2]
            -1.16184125e-09 * tc[3]
            +3.45826200e-13 * tc[4]
            +2.25124300e+04 * invT;
        /*species 23: C3H4-P */
        species[23] =
            +2.02973000e+00
            +7.49480500e-03 * tc[1]
            -4.66166667e-07 * tc[2]
            -9.92404750e-10 * tc[3]
            +2.77643400e-13 * tc[4]
            +2.14840800e+04 * invT;
        /*species 42: CH3CHCO */
        species[42] =
            +4.83801190e-01
            +1.61101506e-02 * tc[1]
            -9.00833443e-06 * tc[2]
            +3.01247910e-09 * tc[3]
            -4.36731862e-13 * tc[4]
            -1.15276540e+04 * invT;
    } else {
        /*species 22: C3H4-A */
        species[22] =
            +8.77625600e+00
            +2.65106900e-03 * tc[1]
            -1.23370600e-07 * tc[2]
            -7.56596500e-11 * tc[3]
            +1.01791620e-14 * tc[4]
            +1.95497200e+04 * invT;
        /*species 23: C3H4-P */
        species[23] =
            +8.76810200e+00
            +2.60957550e-03 * tc[1]
            -1.25104667e-07 * tc[2]
            -7.48047750e-11 * tc[3]
            +1.02157560e-14 * tc[4]
            +1.86027700e+04 * invT;
        /*species 42: CH3CHCO */
        species[42] =
            +9.02191230e+00
            +4.78483150e-03 * tc[1]
            -1.08740548e-06 * tc[2]
            +1.26307927e-10 * tc[3]
            -5.85186514e-15 * tc[4]
            -1.42482738e+04 * invT;
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 24: C3H6 */
        species[24] =
            -6.05384556e-01
            +1.44553831e-02 * tc[1]
            -5.16289360e-06 * tc[2]
            +9.72035522e-10 * tc[3]
            -6.75780704e-14 * tc[4]
            +1.06688164e+03 * invT;
    } else {
        /*species 24: C3H6 */
        species[24] =
            +7.01595958e+00
            +6.85118170e-03 * tc[1]
            -1.55416578e-06 * tc[2]
            +1.80313601e-10 * tc[3]
            -8.34740252e-15 * tc[4]
            -1.87821271e+03 * invT;
    }

    /*species with midpoint at T=1398 kelvin */
    if (T < 1398) {
        /*species 25: C4H6 */
        species[25] =
            -2.43095121e+00
            +2.39353031e-02 * tc[1]
            -1.38482267e-05 * tc[2]
            +4.78873880e-09 * tc[3]
            -7.14317014e-13 * tc[4]
            +1.17551314e+04 * invT;
    } else {
        /*species 25: C4H6 */
        species[25] =
            +1.01633789e+01
            +6.85819825e-03 * tc[1]
            -1.56571928e-06 * tc[2]
            +1.82423459e-10 * tc[3]
            -8.46972406e-15 * tc[4]
            +7.79039770e+03 * invT;
    }

    /*species with midpoint at T=1392 kelvin */
    if (T < 1392) {
        /*species 26: C4H7 */
        species[26] =
            -1.35050835e+00
            +2.13255622e-02 * tc[1]
            -9.69931243e-06 * tc[2]
            +2.63509785e-09 * tc[3]
            -3.20119708e-13 * tc[4]
            +1.49933591e+04 * invT;
        /*species 27: C4H8-1 */
        species[27] =
            -1.83137209e+00
            +2.26290489e-02 * tc[1]
            -9.78861863e-06 * tc[2]
            +2.50551090e-09 * tc[3]
            -2.86383360e-13 * tc[4]
            -1.57875035e+03 * invT;
        /*species 31: C5H9 */
        species[31] =
            -2.38013950e+00
            +2.78804243e-02 * tc[1]
            -1.23381309e-05 * tc[2]
            +3.17209752e-09 * tc[3]
            -3.57077670e-13 * tc[4]
            +1.25589824e+04 * invT;
        /*species 32: C5H10-1 */
        species[32] =
            -2.06223481e+00
            +2.87109147e-02 * tc[1]
            -1.24828963e-05 * tc[2]
            +3.18412472e-09 * tc[3]
            -3.59219578e-13 * tc[4]
            -4.46546666e+03 * invT;
    } else {
        /*species 26: C4H7 */
        species[26] =
            +1.02103578e+01
            +8.02415980e-03 * tc[1]
            -1.82167431e-06 * tc[2]
            +2.11485263e-10 * tc[3]
            -9.79545478e-15 * tc[4]
            +1.09041937e+04 * invT;
        /*species 27: C4H8-1 */
        species[27] =
            +1.03508668e+01
            +9.03089385e-03 * tc[1]
            -2.05364343e-06 * tc[2]
            +2.38663240e-10 * tc[3]
            -1.10617928e-14 * tc[4]
            -5.97871038e+03 * invT;
        /*species 31: C5H9 */
        species[31] =
            +1.31860454e+01
            +1.03564449e-02 * tc[1]
            -2.35653539e-06 * tc[2]
            +2.74017833e-10 * tc[3]
            -1.27064442e-14 * tc[4]
            +7.00496135e+03 * invT;
        /*species 32: C5H10-1 */
        species[32] =
            +1.35851539e+01
            +1.12036235e-02 * tc[1]
            -2.54449342e-06 * tc[2]
            +2.95472415e-10 * tc[3]
            -1.36877028e-14 * tc[4]
            -1.00898205e+04 * invT;
    }

    /*species with midpoint at T=1395 kelvin */
    if (T < 1395) {
        /*species 28: PC4H9 */
        species[28] =
            -1.43777972e+00
            +2.39486182e-02 * tc[1]
            -1.04674386e-05 * tc[2]
            +2.74466180e-09 * tc[3]
            -3.24021328e-13 * tc[4]
            +7.68945248e+03 * invT;
    } else {
        /*species 28: PC4H9 */
        species[28] =
            +1.11510082e+01
            +9.71553585e-03 * tc[1]
            -2.20525983e-06 * tc[2]
            +2.55937840e-10 * tc[3]
            -1.18505941e-14 * tc[4]
            +3.17231942e+03 * invT;
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 29: CH3COCH2 */
        species[29] =
            +8.03391870e-01
            +1.50703542e-02 * tc[1]
            -6.45018507e-06 * tc[2]
            +1.59549759e-09 * tc[3]
            -1.73220636e-13 * tc[4]
            -5.37233261e+03 * invT;
        /*species 45: NC7H16 */
        species[45] =
            -2.26836187e+00
            +4.27177910e-02 * tc[1]
            -1.75115595e-05 * tc[2]
            +4.07364302e-09 * tc[3]
            -4.04789850e-13 * tc[4]
            -2.56586565e+04 * invT;
    } else {
        /*species 29: CH3COCH2 */
        species[29] =
            +9.23039750e+00
            +5.82470805e-03 * tc[1]
            -1.33668512e-06 * tc[2]
            +1.56301312e-10 * tc[3]
            -7.27568724e-15 * tc[4]
            -8.44380462e+03 * invT;
        /*species 45: NC7H16 */
        species[45] =
            +2.12148969e+01
            +1.73837875e-02 * tc[1]
            -3.94690430e-06 * tc[2]
            +4.58246195e-10 * tc[3]
            -2.12260532e-14 * tc[4]
            -3.42760081e+04 * invT;
    }

    /*species with midpoint at T=1378 kelvin */
    if (T < 1378) {
        /*species 30: C2H5CHO */
        species[30] =
            +1.16308444e+00
            +1.47750632e-02 * tc[1]
            -5.08154173e-06 * tc[2]
            +8.73759867e-10 * tc[3]
            -4.77793254e-14 * tc[4]
            -2.42260137e+04 * invT;
        /*species 40: NC3H7CHO */
        species[40] =
            +8.74159590e-01
            +2.09620158e-02 * tc[1]
            -7.83829263e-06 * tc[2]
            +1.56728418e-09 * tc[3]
            -1.21888782e-13 * tc[4]
            -2.71032194e+04 * invT;
    } else {
        /*species 30: C2H5CHO */
        species[30] =
            +9.24276950e+00
            +6.98209945e-03 * tc[1]
            -1.58749334e-06 * tc[2]
            +1.84526426e-10 * tc[3]
            -8.55519006e-15 * tc[4]
            -2.74143135e+04 * invT;
        /*species 40: NC3H7CHO */
        species[40] =
            +1.25988068e+01
            +9.08262370e-03 * tc[1]
            -2.05948153e-06 * tc[2]
            +2.38995052e-10 * tc[3]
            -1.10688592e-14 * tc[4]
            -3.15845348e+04 * invT;
    }

    /*species with midpoint at T=1397 kelvin */
    if (T < 1397) {
        /*species 33: C5H11-1 */
        species[33] =
            -1.90525591e+00
            +3.05316426e-02 * tc[1]
            -1.36497275e-05 * tc[2]
            +3.65233675e-09 * tc[3]
            -4.37719230e-13 * tc[4]
            +4.83995303e+03 * invT;
        /*species 38: C3H5-A */
        species[38] =
            -1.52913196e+00
            +1.67279550e-02 * tc[1]
            -8.44670090e-06 * tc[2]
            +2.57164385e-09 * tc[3]
            -3.46516680e-13 * tc[4]
            +1.93834226e+04 * invT;
    } else {
        /*species 33: C5H11-1 */
        species[33] =
            +1.43234740e+01
            +1.19520600e-02 * tc[1]
            -2.71590540e-06 * tc[2]
            +3.15440590e-10 * tc[3]
            -1.46135467e-14 * tc[4]
            -9.23241637e+02 * invT;
        /*species 38: C3H5-A */
        species[38] =
            +7.45883958e+00
            +5.63477415e-03 * tc[1]
            -1.27930955e-06 * tc[2]
            +1.48514780e-10 * tc[3]
            -6.87836060e-15 * tc[4]
            +1.63576092e+04 * invT;
    }

    /*species with midpoint at T=1385 kelvin */
    if (T < 1385) {
        /*species 34: CH3O2 */
        species[34] =
            +3.26146906e+00
            +5.04367995e-03 * tc[1]
            -1.07168728e-06 * tc[2]
            +5.23523168e-11 * tc[3]
            +8.36678206e-15 * tc[4]
            -6.84394259e+02 * invT;
    } else {
        /*species 34: CH3O2 */
        species[34] =
            +4.95787891e+00
            +3.95364313e-03 * tc[1]
            -8.94154113e-07 * tc[2]
            +1.03472834e-10 * tc[3]
            -4.78014660e-15 * tc[4]
            -1.53574838e+03 * invT;
    }

    /*species with midpoint at T=1402 kelvin */
    if (T < 1402) {
        /*species 36: C2H3CO */
        species[36] =
            +3.62420130e-01
            +1.57636986e-02 * tc[1]
            -1.00072978e-05 * tc[2]
            +3.70417780e-09 * tc[3]
            -5.75943060e-13 * tc[4]
            +4.25770215e+03 * invT;
    } else {
        /*species 36: C2H3CO */
        species[36] =
            +8.37467676e+00
            +3.95648450e-03 * tc[1]
            -8.90660933e-07 * tc[2]
            +1.02778857e-10 * tc[3]
            -4.73957962e-15 * tc[4]
            +1.92969514e+03 * invT;
    }

    /*species with midpoint at T=1393 kelvin */
    if (T < 1393) {
        /*species 37: C2H3CHO */
        species[37] =
            -7.07644838e-01
            +1.77160709e-02 * tc[1]
            -9.83121080e-06 * tc[2]
            +3.20250310e-09 * tc[3]
            -4.52288216e-13 * tc[4]
            -1.16521584e+04 * invT;
    } else {
        /*species 37: C2H3CHO */
        species[37] =
            +9.41849590e+00
            +4.74481661e-03 * tc[1]
            -1.09770176e-06 * tc[2]
            +1.29069801e-10 * tc[3]
            -6.03174582e-15 * tc[4]
            -1.49630281e+04 * invT;
    }

    /*species with midpoint at T=1383 kelvin */
    if (T < 1383) {
        /*species 41: C2H5COCH2 */
        species[41] =
            +5.40138560e-01
            +2.19743129e-02 * tc[1]
            -9.90008070e-06 * tc[2]
            +2.63738282e-09 * tc[3]
            -3.17197538e-13 * tc[4]
            -9.50796505e+03 * invT;
        /*species 43: NC3H7COCH3 */
        species[43] =
            +2.40373610e-01
            +2.71247023e-02 * tc[1]
            -1.04941777e-05 * tc[2]
            +2.28623080e-09 * tc[3]
            -2.15331132e-13 * tc[4]
            -3.36276472e+04 * invT;
    } else {
        /*species 41: C2H5COCH2 */
        species[41] =
            +1.32098738e+01
            +7.89332295e-03 * tc[1]
            -1.83509728e-06 * tc[2]
            +2.16467635e-10 * tc[3]
            -1.01382666e-14 * tc[4]
            -1.41284951e+04 * invT;
        /*species 43: NC3H7COCH3 */
        species[43] =
            +1.51841564e+01
            +1.17622941e-02 * tc[1]
            -2.69891685e-06 * tc[2]
            +3.15576725e-10 * tc[3]
            -1.46897568e-14 * tc[4]
            -3.93375474e+04 * invT;
    }

    /*species with midpoint at T=1389 kelvin */
    if (T < 1389) {
        /*species 44: NC3H7COCH2 */
        species[44] =
            -4.17007290e-02
            +2.84081266e-02 * tc[1]
            -1.33037594e-05 * tc[2]
            +3.81678785e-09 * tc[3]
            -4.98442094e-13 * tc[4]
            -1.23062238e+04 * invT;
    } else {
        /*species 44: NC3H7COCH2 */
        species[44] =
            +1.51502419e+01
            +1.07046733e-02 * tc[1]
            -2.45353205e-06 * tc[2]
            +2.86642040e-10 * tc[3]
            -1.33342544e-14 * tc[4]
            -1.76968561e+04 * invT;
    }

    /*species with midpoint at T=1399 kelvin */
    if (T < 1399) {
        /*species 49: C7H14O */
        species[49] =
            -6.39477206e+00
            +5.09232885e-02 * tc[1]
            -2.53348745e-05 * tc[2]
            +7.41345050e-09 * tc[3]
            -9.49708740e-13 * tc[4]
            -2.67118295e+04 * invT;
    } else {
        /*species 49: C7H14O */
        species[49] =
            +2.22692731e+01
            +1.62792695e-02 * tc[1]
            -3.72083020e-06 * tc[2]
            +4.33934145e-10 * tc[3]
            -2.01626430e-14 * tc[4]
            -3.64349664e+04 * invT;
    }

    /*species with midpoint at T=1387 kelvin */
    if (T < 1387) {
        /*species 50: NC7KET */
        species[50] =
            +5.29366920e-01
            +4.79086733e-02 * tc[1]
            -2.32229507e-05 * tc[2]
            +6.73850955e-09 * tc[3]
            -8.77456252e-13 * tc[4]
            -5.23839542e+04 * invT;
    } else {
        /*species 50: NC7KET */
        species[50] =
            +2.78330388e+01
            +1.60084048e-02 * tc[1]
            -3.71694853e-06 * tc[2]
            +4.38065397e-10 * tc[3]
            -2.05040902e-14 * tc[4]
            -6.20142474e+04 * invT;
    }
    return;
}


/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void speciesEnthalpy(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];
    amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.50104422e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54747466e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +3.29812400e+00
            +4.12472100e-04 * tc[1]
            -2.71433833e-07 * tc[2]
            -2.36885850e-11 * tc[3]
            +8.26974400e-14 * tc[4]
            -1.01252100e+03 * invT;
        /*species 2: O */
        species[2] =
            +2.94642900e+00
            -8.19083000e-04 * tc[1]
            +8.07010667e-07 * tc[2]
            -4.00710750e-10 * tc[3]
            +7.78139200e-14 * tc[4]
            +2.91476400e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +3.21293600e+00
            +5.63743000e-04 * tc[1]
            -1.91871667e-07 * tc[2]
            +3.28469250e-10 * tc[3]
            -1.75371080e-13 * tc[4]
            -1.00524900e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +3.38684200e+00
            +1.73749100e-03 * tc[1]
            -2.11823200e-06 * tc[2]
            +1.74214525e-09 * tc[3]
            -5.01317600e-13 * tc[4]
            -3.02081100e+04 * invT;
        /*species 6: CO */
        species[6] =
            +3.26245200e+00
            +7.55970500e-04 * tc[1]
            -1.29391833e-06 * tc[2]
            +1.39548600e-09 * tc[3]
            -4.94990200e-13 * tc[4]
            -1.43105400e+04 * invT;
        /*species 7: CO2 */
        species[7] =
            +2.27572500e+00
            +4.96103600e-03 * tc[1]
            -3.46970333e-06 * tc[2]
            +1.71667175e-09 * tc[3]
            -4.23456000e-13 * tc[4]
            -4.83731400e+04 * invT;
        /*species 8: CH3 */
        species[8] =
            +2.43044300e+00
            +5.56205000e-03 * tc[1]
            -5.60073333e-06 * tc[2]
            +4.05457250e-09 * tc[3]
            -1.17299060e-12 * tc[4]
            +1.64237800e+04 * invT;
        /*species 9: CH4 */
        species[9] =
            +7.78741500e-01
            +8.73834000e-03 * tc[1]
            -9.27803000e-06 * tc[2]
            +7.62427000e-09 * tc[3]
            -2.44786200e-12 * tc[4]
            -9.82522900e+03 * invT;
        /*species 11: H2O2 */
        species[11] =
            +3.38875400e+00
            +3.28461300e-03 * tc[1]
            -4.95004333e-08 * tc[2]
            -1.15645150e-09 * tc[3]
            +4.94303000e-13 * tc[4]
            -1.76631500e+04 * invT;
        /*species 12: CH2O */
        species[12] =
            +1.65273100e+00
            +6.31572000e-03 * tc[1]
            -6.29389333e-06 * tc[2]
            +5.12507750e-09 * tc[3]
            -1.68264740e-12 * tc[4]
            -1.48654000e+04 * invT;
        /*species 14: C2H4 */
        species[14] =
            -8.61488000e-01
            +1.39808150e-02 * tc[1]
            -1.12955900e-05 * tc[2]
            +6.96288000e-09 * tc[3]
            -1.94757580e-12 * tc[4]
            +5.57304600e+03 * invT;
        /*species 15: C2H5 */
        species[15] =
            +2.69070200e+00
            +4.35956650e-03 * tc[1]
            +1.47327967e-06 * tc[2]
            +2.33467575e-10 * tc[3]
            -7.85554600e-13 * tc[4]
            +1.28704000e+04 * invT;
        /*species 16: C2H */
        species[16] =
            +2.73770400e+00
            +4.02422300e-03 * tc[1]
            -3.08143667e-06 * tc[2]
            +1.63131475e-09 * tc[3]
            -3.87916000e-13 * tc[4]
            +6.68381300e+04 * invT;
        /*species 17: C2H2 */
        species[17] =
            +2.01356200e+00
            +7.59522500e-03 * tc[1]
            -5.38773000e-06 * tc[2]
            +2.26974800e-09 * tc[3]
            -3.82549200e-13 * tc[4]
            +2.61244400e+04 * invT;
        /*species 18: CH3OH */
        species[18] =
            +2.66011500e+00
            +3.67075400e-03 * tc[1]
            +2.39001700e-06 * tc[2]
            -2.19829850e-09 * tc[3]
            +4.78114000e-13 * tc[4]
            -2.53534800e+04 * invT;
        /*species 19: CH2CO */
        species[19] =
            +2.97497100e+00
            +6.05935500e-03 * tc[1]
            -7.81682000e-07 * tc[2]
            -1.61667125e-09 * tc[3]
            +7.81129800e-13 * tc[4]
            -7.63263700e+03 * invT;
        /*species 20: HCCO */
        species[20] =
            +5.04796500e+00
            +2.22673900e-03 * tc[1]
            +7.56094333e-08 * tc[2]
            -3.70523750e-10 * tc[3]
            +4.50148400e-14 * tc[4]
            +1.96589200e+04 * invT;
        /*species 21: CH3CHO */
        species[21] =
            +2.50569500e+00
            +6.68495500e-03 * tc[1]
            +1.55731767e-06 * tc[2]
            -2.82035000e-09 * tc[3]
            +8.52713200e-13 * tc[4]
            -2.12458900e+04 * invT;
        /*species 39: C3H3 */
        species[39] =
            +4.75420000e+00
            +5.54014000e-03 * tc[1]
            +9.31107667e-08 * tc[2]
            -1.36980300e-09 * tc[3]
            +3.89925800e-13 * tc[4]
            +3.98888300e+04 * invT;
        /*species 51: N2 */
        species[51] =
            +3.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
    } else {
        /*species 0: H */
        species[0] =
            +2.50104422e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54747466e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +2.99142300e+00
            +3.50032200e-04 * tc[1]
            -1.87794300e-08 * tc[2]
            -2.30789450e-12 * tc[3]
            +3.16550400e-16 * tc[4]
            -8.35034000e+02 * invT;
        /*species 2: O */
        species[2] =
            +2.54206000e+00
            -1.37753100e-05 * tc[1]
            -1.03426767e-09 * tc[2]
            +1.13776675e-12 * tc[3]
            -8.73610400e-17 * tc[4]
            +2.92308000e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +3.69757800e+00
            +3.06759850e-04 * tc[1]
            -4.19614000e-08 * tc[2]
            +4.43820250e-12 * tc[3]
            -2.27287000e-16 * tc[4]
            -1.23393000e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +2.67214600e+00
            +1.52814650e-03 * tc[1]
            -2.91008667e-07 * tc[2]
            +3.00249000e-11 * tc[3]
            -1.27832360e-15 * tc[4]
            -2.98992100e+04 * invT;
        /*species 6: CO */
        species[6] =
            +3.02507800e+00
            +7.21344500e-04 * tc[1]
            -1.87694267e-07 * tc[2]
            +2.54645250e-11 * tc[3]
            -1.38219040e-15 * tc[4]
            -1.42683500e+04 * invT;
        /*species 7: CO2 */
        species[7] =
            +4.45362300e+00
            +1.57008450e-03 * tc[1]
            -4.26137000e-07 * tc[2]
            +5.98499250e-11 * tc[3]
            -3.33806600e-15 * tc[4]
            -4.89669600e+04 * invT;
        /*species 8: CH3 */
        species[8] =
            +2.84405200e+00
            +3.06898700e-03 * tc[1]
            -7.43448333e-07 * tc[2]
            +9.46290250e-11 * tc[3]
            -4.90431800e-15 * tc[4]
            +1.64378100e+04 * invT;
        /*species 9: CH4 */
        species[9] =
            +1.68347900e+00
            +5.11862000e-03 * tc[1]
            -1.29170967e-06 * tc[2]
            +1.69639625e-10 * tc[3]
            -9.00684600e-15 * tc[4]
            -1.00807900e+04 * invT;
        /*species 11: H2O2 */
        species[11] =
            +4.57316700e+00
            +2.16806800e-03 * tc[1]
            -4.91563000e-07 * tc[2]
            +5.87226000e-11 * tc[3]
            -2.86330800e-15 * tc[4]
            -1.80069600e+04 * invT;
        /*species 12: CH2O */
        species[12] =
            +2.99560600e+00
            +3.34066050e-03 * tc[1]
            -8.76318333e-07 * tc[2]
            +1.18428825e-10 * tc[3]
            -6.42503400e-15 * tc[4]
            -1.53203700e+04 * invT;
        /*species 14: C2H4 */
        species[14] =
            +3.52841900e+00
            +5.74259000e-03 * tc[1]
            -1.47279500e-06 * tc[2]
            +1.96115025e-10 * tc[3]
            -1.05336960e-14 * tc[4]
            +4.42828900e+03 * invT;
        /*species 15: C2H5 */
        species[15] =
            +7.19048000e+00
            +3.24203850e-03 * tc[1]
            -2.14268833e-07 * tc[2]
            -5.86969750e-11 * tc[3]
            +7.76175400e-15 * tc[4]
            +1.06745500e+04 * invT;
        /*species 16: C2H */
        species[16] =
            +3.98636700e+00
            +1.57156150e-03 * tc[1]
            -4.22414333e-07 * tc[2]
            +7.31090750e-11 * tc[3]
            -5.43264000e-15 * tc[4]
            +6.65588400e+04 * invT;
        /*species 17: C2H2 */
        species[17] =
            +4.43677000e+00
            +2.68801950e-03 * tc[1]
            -6.37605667e-07 * tc[2]
            +8.21594750e-11 * tc[3]
            -4.31342000e-15 * tc[4]
            +2.56676600e+04 * invT;
        /*species 18: CH3OH */
        species[18] =
            +4.02906100e+00
            +4.68829650e-03 * tc[1]
            -1.01675133e-06 * tc[2]
            +1.08969825e-10 * tc[3]
            -4.44944600e-15 * tc[4]
            -2.61579100e+04 * invT;
        /*species 19: CH2CO */
        species[19] =
            +6.03881700e+00
            +2.90242000e-03 * tc[1]
            -6.40318000e-07 * tc[2]
            +6.98621250e-11 * tc[3]
            -2.91773600e-15 * tc[4]
            -8.58340200e+03 * invT;
        /*species 20: HCCO */
        species[20] =
            +6.75807300e+00
            +1.00020000e-03 * tc[1]
            -6.75869000e-08 * tc[2]
            -2.60283000e-11 * tc[3]
            +3.93033000e-15 * tc[4]
            +1.90151300e+04 * invT;
        /*species 21: CH3CHO */
        species[21] =
            +5.86865000e+00
            +5.39712000e-03 * tc[1]
            -1.21517667e-06 * tc[2]
            +1.35322800e-10 * tc[3]
            -5.79368800e-15 * tc[4]
            -2.26456900e+04 * invT;
        /*species 39: C3H3 */
        species[39] =
            +8.83104700e+00
            +2.17859750e-03 * tc[1]
            -1.36968900e-07 * tc[2]
            -5.92180750e-11 * tc[3]
            +8.75304000e-15 * tc[4]
            +3.84742000e+04 * invT;
        /*species 51: N2 */
        species[51] =
            +2.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }

    /*species with midpoint at T=1357 kelvin */
    if (T < 1357) {
        /*species 4: OH */
        species[4] =
            +3.43586219e+00
            +1.01117902e-04 * tc[1]
            -3.78488040e-08 * tc[2]
            +6.06112872e-11 * tc[3]
            -1.48730206e-14 * tc[4]
            +3.74321252e+03 * invT;
    } else {
        /*species 4: OH */
        species[4] =
            +2.62599754e+00
            +6.59962030e-04 * tc[1]
            -1.19908223e-07 * tc[2]
            +1.06407700e-11 * tc[3]
            -3.64096032e-16 * tc[4]
            +4.12085374e+03 * invT;
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 10: HO2 */
        species[10] =
            +3.18310656e+00
            +1.83383975e-03 * tc[1]
            -3.10795041e-07 * tc[2]
            -8.14632298e-11 * tc[3]
            +3.02279824e-14 * tc[4]
            +8.09181013e+02 * invT;
        /*species 35: CH3O2H */
        species[35] =
            +3.23442817e+00
            +9.50648835e-03 * tc[1]
            -3.77954290e-06 * tc[2]
            +8.50766632e-10 * tc[3]
            -8.23660444e-14 * tc[4]
            -1.77197926e+04 * invT;
        /*species 46: C7H15-1 */
        species[46] =
            -4.99570406e-01
            +4.04413234e-02 * tc[1]
            -1.66844251e-05 * tc[2]
            +3.91373270e-09 * tc[3]
            -3.93232454e-13 * tc[4]
            -1.04590223e+03 * invT;
        /*species 47: C7H15O2 */
        species[47] =
            +3.55252917e-01
            +4.71190503e-02 * tc[1]
            -2.22318369e-05 * tc[2]
            +6.36988958e-09 * tc[3]
            -8.26422210e-13 * tc[4]
            -2.15795254e+04 * invT;
        /*species 48: C7H14OOHO2 */
        species[48] =
            +1.40912479e+00
            +5.53604170e-02 * tc[1]
            -2.90392692e-05 * tc[2]
            +9.09156872e-09 * tc[3]
            -1.25670138e-12 * tc[4]
            -3.46304788e+04 * invT;
    } else {
        /*species 10: HO2 */
        species[10] =
            +4.10547423e+00
            +1.19226417e-03 * tc[1]
            -2.68782663e-07 * tc[2]
            +3.10479308e-11 * tc[3]
            -1.43280022e-15 * tc[4]
            +3.98127689e+02 * invT;
        /*species 35: CH3O2H */
        species[35] =
            +8.43117091e+00
            +4.03408955e-03 * tc[1]
            -9.23649737e-07 * tc[2]
            +1.07833061e-10 * tc[3]
            -5.01384292e-15 * tc[4]
            -1.96678771e+04 * invT;
        /*species 46: C7H15-1 */
        species[46] =
            +2.17940709e+01
            +1.63140122e-02 * tc[1]
            -3.70460813e-06 * tc[2]
            +4.30167870e-10 * tc[3]
            -1.99273400e-14 * tc[4]
            -9.20938221e+03 * invT;
        /*species 47: C7H15O2 */
        species[47] =
            +2.60535640e+01
            +1.71915844e-02 * tc[1]
            -3.95710823e-06 * tc[2]
            +4.63625602e-10 * tc[3]
            -2.16104798e-14 * tc[4]
            -3.06842095e+04 * invT;
        /*species 48: C7H14OOHO2 */
        species[48] =
            +3.28358296e+01
            +1.66383642e-02 * tc[1]
            -3.86110950e-06 * tc[2]
            +4.54896595e-10 * tc[3]
            -2.12868838e-14 * tc[4]
            -4.52873003e+04 * invT;
    }

    /*species with midpoint at T=1384 kelvin */
    if (T < 1384) {
        /*species 13: C2H6 */
        species[13] =
            -2.52854344e-02
            +1.20382377e-02 * tc[1]
            -3.72978240e-06 * tc[2]
            +5.20852252e-10 * tc[3]
            -1.05973723e-14 * tc[4]
            -1.12345534e+04 * invT;
    } else {
        /*species 13: C2H6 */
        species[13] =
            +6.10683385e+00
            +6.46181805e-03 * tc[1]
            -1.47509065e-06 * tc[2]
            +1.71847932e-10 * tc[3]
            -7.97803464e-15 * tc[4]
            -1.37500014e+04 * invT;
    }

    /*species with midpoint at T=1400 kelvin */
    if (T < 1400) {
        /*species 22: C3H4-A */
        species[22] =
            +2.53983100e+00
            +8.16718500e-03 * tc[1]
            -5.88316667e-07 * tc[2]
            -1.16184125e-09 * tc[3]
            +3.45826200e-13 * tc[4]
            +2.25124300e+04 * invT;
        /*species 23: C3H4-P */
        species[23] =
            +3.02973000e+00
            +7.49480500e-03 * tc[1]
            -4.66166667e-07 * tc[2]
            -9.92404750e-10 * tc[3]
            +2.77643400e-13 * tc[4]
            +2.14840800e+04 * invT;
        /*species 42: CH3CHCO */
        species[42] =
            +1.48380119e+00
            +1.61101506e-02 * tc[1]
            -9.00833443e-06 * tc[2]
            +3.01247910e-09 * tc[3]
            -4.36731862e-13 * tc[4]
            -1.15276540e+04 * invT;
    } else {
        /*species 22: C3H4-A */
        species[22] =
            +9.77625600e+00
            +2.65106900e-03 * tc[1]
            -1.23370600e-07 * tc[2]
            -7.56596500e-11 * tc[3]
            +1.01791620e-14 * tc[4]
            +1.95497200e+04 * invT;
        /*species 23: C3H4-P */
        species[23] =
            +9.76810200e+00
            +2.60957550e-03 * tc[1]
            -1.25104667e-07 * tc[2]
            -7.48047750e-11 * tc[3]
            +1.02157560e-14 * tc[4]
            +1.86027700e+04 * invT;
        /*species 42: CH3CHCO */
        species[42] =
            +1.00219123e+01
            +4.78483150e-03 * tc[1]
            -1.08740548e-06 * tc[2]
            +1.26307927e-10 * tc[3]
            -5.85186514e-15 * tc[4]
            -1.42482738e+04 * invT;
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 24: C3H6 */
        species[24] =
            +3.94615444e-01
            +1.44553831e-02 * tc[1]
            -5.16289360e-06 * tc[2]
            +9.72035522e-10 * tc[3]
            -6.75780704e-14 * tc[4]
            +1.06688164e+03 * invT;
    } else {
        /*species 24: C3H6 */
        species[24] =
            +8.01595958e+00
            +6.85118170e-03 * tc[1]
            -1.55416578e-06 * tc[2]
            +1.80313601e-10 * tc[3]
            -8.34740252e-15 * tc[4]
            -1.87821271e+03 * invT;
    }

    /*species with midpoint at T=1398 kelvin */
    if (T < 1398) {
        /*species 25: C4H6 */
        species[25] =
            -1.43095121e+00
            +2.39353031e-02 * tc[1]
            -1.38482267e-05 * tc[2]
            +4.78873880e-09 * tc[3]
            -7.14317014e-13 * tc[4]
            +1.17551314e+04 * invT;
    } else {
        /*species 25: C4H6 */
        species[25] =
            +1.11633789e+01
            +6.85819825e-03 * tc[1]
            -1.56571928e-06 * tc[2]
            +1.82423459e-10 * tc[3]
            -8.46972406e-15 * tc[4]
            +7.79039770e+03 * invT;
    }

    /*species with midpoint at T=1392 kelvin */
    if (T < 1392) {
        /*species 26: C4H7 */
        species[26] =
            -3.50508352e-01
            +2.13255622e-02 * tc[1]
            -9.69931243e-06 * tc[2]
            +2.63509785e-09 * tc[3]
            -3.20119708e-13 * tc[4]
            +1.49933591e+04 * invT;
        /*species 27: C4H8-1 */
        species[27] =
            -8.31372089e-01
            +2.26290489e-02 * tc[1]
            -9.78861863e-06 * tc[2]
            +2.50551090e-09 * tc[3]
            -2.86383360e-13 * tc[4]
            -1.57875035e+03 * invT;
        /*species 31: C5H9 */
        species[31] =
            -1.38013950e+00
            +2.78804243e-02 * tc[1]
            -1.23381309e-05 * tc[2]
            +3.17209752e-09 * tc[3]
            -3.57077670e-13 * tc[4]
            +1.25589824e+04 * invT;
        /*species 32: C5H10-1 */
        species[32] =
            -1.06223481e+00
            +2.87109147e-02 * tc[1]
            -1.24828963e-05 * tc[2]
            +3.18412472e-09 * tc[3]
            -3.59219578e-13 * tc[4]
            -4.46546666e+03 * invT;
    } else {
        /*species 26: C4H7 */
        species[26] =
            +1.12103578e+01
            +8.02415980e-03 * tc[1]
            -1.82167431e-06 * tc[2]
            +2.11485263e-10 * tc[3]
            -9.79545478e-15 * tc[4]
            +1.09041937e+04 * invT;
        /*species 27: C4H8-1 */
        species[27] =
            +1.13508668e+01
            +9.03089385e-03 * tc[1]
            -2.05364343e-06 * tc[2]
            +2.38663240e-10 * tc[3]
            -1.10617928e-14 * tc[4]
            -5.97871038e+03 * invT;
        /*species 31: C5H9 */
        species[31] =
            +1.41860454e+01
            +1.03564449e-02 * tc[1]
            -2.35653539e-06 * tc[2]
            +2.74017833e-10 * tc[3]
            -1.27064442e-14 * tc[4]
            +7.00496135e+03 * invT;
        /*species 32: C5H10-1 */
        species[32] =
            +1.45851539e+01
            +1.12036235e-02 * tc[1]
            -2.54449342e-06 * tc[2]
            +2.95472415e-10 * tc[3]
            -1.36877028e-14 * tc[4]
            -1.00898205e+04 * invT;
    }

    /*species with midpoint at T=1395 kelvin */
    if (T < 1395) {
        /*species 28: PC4H9 */
        species[28] =
            -4.37779725e-01
            +2.39486182e-02 * tc[1]
            -1.04674386e-05 * tc[2]
            +2.74466180e-09 * tc[3]
            -3.24021328e-13 * tc[4]
            +7.68945248e+03 * invT;
    } else {
        /*species 28: PC4H9 */
        species[28] =
            +1.21510082e+01
            +9.71553585e-03 * tc[1]
            -2.20525983e-06 * tc[2]
            +2.55937840e-10 * tc[3]
            -1.18505941e-14 * tc[4]
            +3.17231942e+03 * invT;
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 29: CH3COCH2 */
        species[29] =
            +1.80339187e+00
            +1.50703542e-02 * tc[1]
            -6.45018507e-06 * tc[2]
            +1.59549759e-09 * tc[3]
            -1.73220636e-13 * tc[4]
            -5.37233261e+03 * invT;
        /*species 45: NC7H16 */
        species[45] =
            -1.26836187e+00
            +4.27177910e-02 * tc[1]
            -1.75115595e-05 * tc[2]
            +4.07364302e-09 * tc[3]
            -4.04789850e-13 * tc[4]
            -2.56586565e+04 * invT;
    } else {
        /*species 29: CH3COCH2 */
        species[29] =
            +1.02303975e+01
            +5.82470805e-03 * tc[1]
            -1.33668512e-06 * tc[2]
            +1.56301312e-10 * tc[3]
            -7.27568724e-15 * tc[4]
            -8.44380462e+03 * invT;
        /*species 45: NC7H16 */
        species[45] =
            +2.22148969e+01
            +1.73837875e-02 * tc[1]
            -3.94690430e-06 * tc[2]
            +4.58246195e-10 * tc[3]
            -2.12260532e-14 * tc[4]
            -3.42760081e+04 * invT;
    }

    /*species with midpoint at T=1378 kelvin */
    if (T < 1378) {
        /*species 30: C2H5CHO */
        species[30] =
            +2.16308444e+00
            +1.47750632e-02 * tc[1]
            -5.08154173e-06 * tc[2]
            +8.73759867e-10 * tc[3]
            -4.77793254e-14 * tc[4]
            -2.42260137e+04 * invT;
        /*species 40: NC3H7CHO */
        species[40] =
            +1.87415959e+00
            +2.09620158e-02 * tc[1]
            -7.83829263e-06 * tc[2]
            +1.56728418e-09 * tc[3]
            -1.21888782e-13 * tc[4]
            -2.71032194e+04 * invT;
    } else {
        /*species 30: C2H5CHO */
        species[30] =
            +1.02427695e+01
            +6.98209945e-03 * tc[1]
            -1.58749334e-06 * tc[2]
            +1.84526426e-10 * tc[3]
            -8.55519006e-15 * tc[4]
            -2.74143135e+04 * invT;
        /*species 40: NC3H7CHO */
        species[40] =
            +1.35988068e+01
            +9.08262370e-03 * tc[1]
            -2.05948153e-06 * tc[2]
            +2.38995052e-10 * tc[3]
            -1.10688592e-14 * tc[4]
            -3.15845348e+04 * invT;
    }

    /*species with midpoint at T=1397 kelvin */
    if (T < 1397) {
        /*species 33: C5H11-1 */
        species[33] =
            -9.05255912e-01
            +3.05316426e-02 * tc[1]
            -1.36497275e-05 * tc[2]
            +3.65233675e-09 * tc[3]
            -4.37719230e-13 * tc[4]
            +4.83995303e+03 * invT;
        /*species 38: C3H5-A */
        species[38] =
            -5.29131958e-01
            +1.67279550e-02 * tc[1]
            -8.44670090e-06 * tc[2]
            +2.57164385e-09 * tc[3]
            -3.46516680e-13 * tc[4]
            +1.93834226e+04 * invT;
    } else {
        /*species 33: C5H11-1 */
        species[33] =
            +1.53234740e+01
            +1.19520600e-02 * tc[1]
            -2.71590540e-06 * tc[2]
            +3.15440590e-10 * tc[3]
            -1.46135467e-14 * tc[4]
            -9.23241637e+02 * invT;
        /*species 38: C3H5-A */
        species[38] =
            +8.45883958e+00
            +5.63477415e-03 * tc[1]
            -1.27930955e-06 * tc[2]
            +1.48514780e-10 * tc[3]
            -6.87836060e-15 * tc[4]
            +1.63576092e+04 * invT;
    }

    /*species with midpoint at T=1385 kelvin */
    if (T < 1385) {
        /*species 34: CH3O2 */
        species[34] =
            +4.26146906e+00
            +5.04367995e-03 * tc[1]
            -1.07168728e-06 * tc[2]
            +5.23523168e-11 * tc[3]
            +8.36678206e-15 * tc[4]
            -6.84394259e+02 * invT;
    } else {
        /*species 34: CH3O2 */
        species[34] =
            +5.95787891e+00
            +3.95364313e-03 * tc[1]
            -8.94154113e-07 * tc[2]
            +1.03472834e-10 * tc[3]
            -4.78014660e-15 * tc[4]
            -1.53574838e+03 * invT;
    }

    /*species with midpoint at T=1402 kelvin */
    if (T < 1402) {
        /*species 36: C2H3CO */
        species[36] =
            +1.36242013e+00
            +1.57636986e-02 * tc[1]
            -1.00072978e-05 * tc[2]
            +3.70417780e-09 * tc[3]
            -5.75943060e-13 * tc[4]
            +4.25770215e+03 * invT;
    } else {
        /*species 36: C2H3CO */
        species[36] =
            +9.37467676e+00
            +3.95648450e-03 * tc[1]
            -8.90660933e-07 * tc[2]
            +1.02778857e-10 * tc[3]
            -4.73957962e-15 * tc[4]
            +1.92969514e+03 * invT;
    }

    /*species with midpoint at T=1393 kelvin */
    if (T < 1393) {
        /*species 37: C2H3CHO */
        species[37] =
            +2.92355162e-01
            +1.77160709e-02 * tc[1]
            -9.83121080e-06 * tc[2]
            +3.20250310e-09 * tc[3]
            -4.52288216e-13 * tc[4]
            -1.16521584e+04 * invT;
    } else {
        /*species 37: C2H3CHO */
        species[37] =
            +1.04184959e+01
            +4.74481661e-03 * tc[1]
            -1.09770176e-06 * tc[2]
            +1.29069801e-10 * tc[3]
            -6.03174582e-15 * tc[4]
            -1.49630281e+04 * invT;
    }

    /*species with midpoint at T=1383 kelvin */
    if (T < 1383) {
        /*species 41: C2H5COCH2 */
        species[41] =
            +1.54013856e+00
            +2.19743129e-02 * tc[1]
            -9.90008070e-06 * tc[2]
            +2.63738282e-09 * tc[3]
            -3.17197538e-13 * tc[4]
            -9.50796505e+03 * invT;
        /*species 43: NC3H7COCH3 */
        species[43] =
            +1.24037361e+00
            +2.71247023e-02 * tc[1]
            -1.04941777e-05 * tc[2]
            +2.28623080e-09 * tc[3]
            -2.15331132e-13 * tc[4]
            -3.36276472e+04 * invT;
    } else {
        /*species 41: C2H5COCH2 */
        species[41] =
            +1.42098738e+01
            +7.89332295e-03 * tc[1]
            -1.83509728e-06 * tc[2]
            +2.16467635e-10 * tc[3]
            -1.01382666e-14 * tc[4]
            -1.41284951e+04 * invT;
        /*species 43: NC3H7COCH3 */
        species[43] =
            +1.61841564e+01
            +1.17622941e-02 * tc[1]
            -2.69891685e-06 * tc[2]
            +3.15576725e-10 * tc[3]
            -1.46897568e-14 * tc[4]
            -3.93375474e+04 * invT;
    }

    /*species with midpoint at T=1389 kelvin */
    if (T < 1389) {
        /*species 44: NC3H7COCH2 */
        species[44] =
            +9.58299271e-01
            +2.84081266e-02 * tc[1]
            -1.33037594e-05 * tc[2]
            +3.81678785e-09 * tc[3]
            -4.98442094e-13 * tc[4]
            -1.23062238e+04 * invT;
    } else {
        /*species 44: NC3H7COCH2 */
        species[44] =
            +1.61502419e+01
            +1.07046733e-02 * tc[1]
            -2.45353205e-06 * tc[2]
            +2.86642040e-10 * tc[3]
            -1.33342544e-14 * tc[4]
            -1.76968561e+04 * invT;
    }

    /*species with midpoint at T=1399 kelvin */
    if (T < 1399) {
        /*species 49: C7H14O */
        species[49] =
            -5.39477206e+00
            +5.09232885e-02 * tc[1]
            -2.53348745e-05 * tc[2]
            +7.41345050e-09 * tc[3]
            -9.49708740e-13 * tc[4]
            -2.67118295e+04 * invT;
    } else {
        /*species 49: C7H14O */
        species[49] =
            +2.32692731e+01
            +1.62792695e-02 * tc[1]
            -3.72083020e-06 * tc[2]
            +4.33934145e-10 * tc[3]
            -2.01626430e-14 * tc[4]
            -3.64349664e+04 * invT;
    }

    /*species with midpoint at T=1387 kelvin */
    if (T < 1387) {
        /*species 50: NC7KET */
        species[50] =
            +1.52936692e+00
            +4.79086733e-02 * tc[1]
            -2.32229507e-05 * tc[2]
            +6.73850955e-09 * tc[3]
            -8.77456252e-13 * tc[4]
            -5.23839542e+04 * invT;
    } else {
        /*species 50: NC7KET */
        species[50] =
            +2.88330388e+01
            +1.60084048e-02 * tc[1]
            -3.71694853e-06 * tc[2]
            +4.38065397e-10 * tc[3]
            -2.05040902e-14 * tc[4]
            -6.20142474e+04 * invT;
    }
    return;
}


/*compute the S/R at the given temperature (Eq 21) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void speciesEntropy(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.50104422e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.65341317e-01 ;
        /*species 1: H2 */
        species[1] =
            +3.29812400e+00 * tc[0]
            +8.24944200e-04 * tc[1]
            -4.07150750e-07 * tc[2]
            -3.15847800e-11 * tc[3]
            +1.03371800e-13 * tc[4]
            -3.29409400e+00 ;
        /*species 2: O */
        species[2] =
            +2.94642900e+00 * tc[0]
            -1.63816600e-03 * tc[1]
            +1.21051600e-06 * tc[2]
            -5.34281000e-10 * tc[3]
            +9.72674000e-14 * tc[4]
            +2.96399500e+00 ;
        /*species 3: O2 */
        species[3] =
            +3.21293600e+00 * tc[0]
            +1.12748600e-03 * tc[1]
            -2.87807500e-07 * tc[2]
            +4.37959000e-10 * tc[3]
            -2.19213850e-13 * tc[4]
            +6.03473800e+00 ;
        /*species 5: H2O */
        species[5] =
            +3.38684200e+00 * tc[0]
            +3.47498200e-03 * tc[1]
            -3.17734800e-06 * tc[2]
            +2.32286033e-09 * tc[3]
            -6.26647000e-13 * tc[4]
            +2.59023300e+00 ;
        /*species 6: CO */
        species[6] =
            +3.26245200e+00 * tc[0]
            +1.51194100e-03 * tc[1]
            -1.94087750e-06 * tc[2]
            +1.86064800e-09 * tc[3]
            -6.18737750e-13 * tc[4]
            +4.84889700e+00 ;
        /*species 7: CO2 */
        species[7] =
            +2.27572500e+00 * tc[0]
            +9.92207200e-03 * tc[1]
            -5.20455500e-06 * tc[2]
            +2.28889567e-09 * tc[3]
            -5.29320000e-13 * tc[4]
            +1.01884900e+01 ;
        /*species 8: CH3 */
        species[8] =
            +2.43044300e+00 * tc[0]
            +1.11241000e-02 * tc[1]
            -8.40110000e-06 * tc[2]
            +5.40609667e-09 * tc[3]
            -1.46623825e-12 * tc[4]
            +6.78979400e+00 ;
        /*species 9: CH4 */
        species[9] =
            +7.78741500e-01 * tc[0]
            +1.74766800e-02 * tc[1]
            -1.39170450e-05 * tc[2]
            +1.01656933e-08 * tc[3]
            -3.05982750e-12 * tc[4]
            +1.37221900e+01 ;
        /*species 11: H2O2 */
        species[11] =
            +3.38875400e+00 * tc[0]
            +6.56922600e-03 * tc[1]
            -7.42506500e-08 * tc[2]
            -1.54193533e-09 * tc[3]
            +6.17878750e-13 * tc[4]
            +6.78536300e+00 ;
        /*species 12: CH2O */
        species[12] =
            +1.65273100e+00 * tc[0]
            +1.26314400e-02 * tc[1]
            -9.44084000e-06 * tc[2]
            +6.83343667e-09 * tc[3]
            -2.10330925e-12 * tc[4]
            +1.37848200e+01 ;
        /*species 14: C2H4 */
        species[14] =
            -8.61488000e-01 * tc[0]
            +2.79616300e-02 * tc[1]
            -1.69433850e-05 * tc[2]
            +9.28384000e-09 * tc[3]
            -2.43446975e-12 * tc[4]
            +2.42114900e+01 ;
        /*species 15: C2H5 */
        species[15] =
            +2.69070200e+00 * tc[0]
            +8.71913300e-03 * tc[1]
            +2.20991950e-06 * tc[2]
            +3.11290100e-10 * tc[3]
            -9.81943250e-13 * tc[4]
            +1.21382000e+01 ;
        /*species 16: C2H */
        species[16] =
            +2.73770400e+00 * tc[0]
            +8.04844600e-03 * tc[1]
            -4.62215500e-06 * tc[2]
            +2.17508633e-09 * tc[3]
            -4.84895000e-13 * tc[4]
            +7.30022000e+00 ;
        /*species 17: C2H2 */
        species[17] =
            +2.01356200e+00 * tc[0]
            +1.51904500e-02 * tc[1]
            -8.08159500e-06 * tc[2]
            +3.02633067e-09 * tc[3]
            -4.78186500e-13 * tc[4]
            +8.80537800e+00 ;
        /*species 18: CH3OH */
        species[18] =
            +2.66011500e+00 * tc[0]
            +7.34150800e-03 * tc[1]
            +3.58502550e-06 * tc[2]
            -2.93106467e-09 * tc[3]
            +5.97642500e-13 * tc[4]
            +1.12326300e+01 ;
        /*species 19: CH2CO */
        species[19] =
            +2.97497100e+00 * tc[0]
            +1.21187100e-02 * tc[1]
            -1.17252300e-06 * tc[2]
            -2.15556167e-09 * tc[3]
            +9.76412250e-13 * tc[4]
            +8.67355300e+00 ;
        /*species 20: HCCO */
        species[20] =
            +5.04796500e+00 * tc[0]
            +4.45347800e-03 * tc[1]
            +1.13414150e-07 * tc[2]
            -4.94031667e-10 * tc[3]
            +5.62685500e-14 * tc[4]
            +4.81843900e-01 ;
        /*species 21: CH3CHO */
        species[21] =
            +2.50569500e+00 * tc[0]
            +1.33699100e-02 * tc[1]
            +2.33597650e-06 * tc[2]
            -3.76046667e-09 * tc[3]
            +1.06589150e-12 * tc[4]
            +1.33508900e+01 ;
        /*species 39: C3H3 */
        species[39] =
            +4.75420000e+00 * tc[0]
            +1.10802800e-02 * tc[1]
            +1.39666150e-07 * tc[2]
            -1.82640400e-09 * tc[3]
            +4.87407250e-13 * tc[4]
            +5.85454900e-01 ;
        /*species 51: N2 */
        species[51] =
            +3.29867700e+00 * tc[0]
            +1.40824000e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213750e-13 * tc[4]
            +3.95037200e+00 ;
    } else {
        /*species 0: H */
        species[0] =
            +2.50104422e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.65341317e-01 ;
        /*species 1: H2 */
        species[1] =
            +2.99142300e+00 * tc[0]
            +7.00064400e-04 * tc[1]
            -2.81691450e-08 * tc[2]
            -3.07719267e-12 * tc[3]
            +3.95688000e-16 * tc[4]
            -1.35511000e+00 ;
        /*species 2: O */
        species[2] =
            +2.54206000e+00 * tc[0]
            -2.75506200e-05 * tc[1]
            -1.55140150e-09 * tc[2]
            +1.51702233e-12 * tc[3]
            -1.09201300e-16 * tc[4]
            +4.92030800e+00 ;
        /*species 3: O2 */
        species[3] =
            +3.69757800e+00 * tc[0]
            +6.13519700e-04 * tc[1]
            -6.29421000e-08 * tc[2]
            +5.91760333e-12 * tc[3]
            -2.84108750e-16 * tc[4]
            +3.18916600e+00 ;
        /*species 5: H2O */
        species[5] =
            +2.67214600e+00 * tc[0]
            +3.05629300e-03 * tc[1]
            -4.36513000e-07 * tc[2]
            +4.00332000e-11 * tc[3]
            -1.59790450e-15 * tc[4]
            +6.86281700e+00 ;
        /*species 6: CO */
        species[6] =
            +3.02507800e+00 * tc[0]
            +1.44268900e-03 * tc[1]
            -2.81541400e-07 * tc[2]
            +3.39527000e-11 * tc[3]
            -1.72773800e-15 * tc[4]
            +6.10821800e+00 ;
        /*species 7: CO2 */
        species[7] =
            +4.45362300e+00 * tc[0]
            +3.14016900e-03 * tc[1]
            -6.39205500e-07 * tc[2]
            +7.97999000e-11 * tc[3]
            -4.17258250e-15 * tc[4]
            -9.55395900e-01 ;
        /*species 8: CH3 */
        species[8] =
            +2.84405200e+00 * tc[0]
            +6.13797400e-03 * tc[1]
            -1.11517250e-06 * tc[2]
            +1.26172033e-10 * tc[3]
            -6.13039750e-15 * tc[4]
            +5.45269700e+00 ;
        /*species 9: CH4 */
        species[9] =
            +1.68347900e+00 * tc[0]
            +1.02372400e-02 * tc[1]
            -1.93756450e-06 * tc[2]
            +2.26186167e-10 * tc[3]
            -1.12585575e-14 * tc[4]
            +9.62339500e+00 ;
        /*species 11: H2O2 */
        species[11] =
            +4.57316700e+00 * tc[0]
            +4.33613600e-03 * tc[1]
            -7.37344500e-07 * tc[2]
            +7.82968000e-11 * tc[3]
            -3.57913500e-15 * tc[4]
            +5.01137000e-01 ;
        /*species 12: CH2O */
        species[12] =
            +2.99560600e+00 * tc[0]
            +6.68132100e-03 * tc[1]
            -1.31447750e-06 * tc[2]
            +1.57905100e-10 * tc[3]
            -8.03129250e-15 * tc[4]
            +6.91257200e+00 ;
        /*species 14: C2H4 */
        species[14] =
            +3.52841900e+00 * tc[0]
            +1.14851800e-02 * tc[1]
            -2.20919250e-06 * tc[2]
            +2.61486700e-10 * tc[3]
            -1.31671200e-14 * tc[4]
            +2.23038900e+00 ;
        /*species 15: C2H5 */
        species[15] =
            +7.19048000e+00 * tc[0]
            +6.48407700e-03 * tc[1]
            -3.21403250e-07 * tc[2]
            -7.82626333e-11 * tc[3]
            +9.70219250e-15 * tc[4]
            -1.47808900e+01 ;
        /*species 16: C2H */
        species[16] =
            +3.98636700e+00 * tc[0]
            +3.14312300e-03 * tc[1]
            -6.33621500e-07 * tc[2]
            +9.74787667e-11 * tc[3]
            -6.79080000e-15 * tc[4]
            +1.19106300e+00 ;
        /*species 17: C2H2 */
        species[17] =
            +4.43677000e+00 * tc[0]
            +5.37603900e-03 * tc[1]
            -9.56408500e-07 * tc[2]
            +1.09545967e-10 * tc[3]
            -5.39177500e-15 * tc[4]
            -2.80033800e+00 ;
        /*species 18: CH3OH */
        species[18] =
            +4.02906100e+00 * tc[0]
            +9.37659300e-03 * tc[1]
            -1.52512700e-06 * tc[2]
            +1.45293100e-10 * tc[3]
            -5.56180750e-15 * tc[4]
            +2.37819600e+00 ;
        /*species 19: CH2CO */
        species[19] =
            +6.03881700e+00 * tc[0]
            +5.80484000e-03 * tc[1]
            -9.60477000e-07 * tc[2]
            +9.31495000e-11 * tc[3]
            -3.64717000e-15 * tc[4]
            -7.65758100e+00 ;
        /*species 20: HCCO */
        species[20] =
            +6.75807300e+00 * tc[0]
            +2.00040000e-03 * tc[1]
            -1.01380350e-07 * tc[2]
            -3.47044000e-11 * tc[3]
            +4.91291250e-15 * tc[4]
            -9.07126200e+00 ;
        /*species 21: CH3CHO */
        species[21] =
            +5.86865000e+00 * tc[0]
            +1.07942400e-02 * tc[1]
            -1.82276500e-06 * tc[2]
            +1.80430400e-10 * tc[3]
            -7.24211000e-15 * tc[4]
            -6.01294600e+00 ;
        /*species 39: C3H3 */
        species[39] =
            +8.83104700e+00 * tc[0]
            +4.35719500e-03 * tc[1]
            -2.05453350e-07 * tc[2]
            -7.89574333e-11 * tc[3]
            +1.09413000e-14 * tc[4]
            -2.17791900e+01 ;
        /*species 51: N2 */
        species[51] =
            +2.92664000e+00 * tc[0]
            +1.48797700e-03 * tc[1]
            -2.84238050e-07 * tc[2]
            +3.36568000e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
    }

    /*species with midpoint at T=1357 kelvin */
    if (T < 1357) {
        /*species 4: OH */
        species[4] =
            +3.43586219e+00 * tc[0]
            +2.02235804e-04 * tc[1]
            -5.67732060e-08 * tc[2]
            +8.08150497e-11 * tc[3]
            -1.85912758e-14 * tc[4]
            +2.45014127e+00 ;
    } else {
        /*species 4: OH */
        species[4] =
            +2.62599754e+00 * tc[0]
            +1.31992406e-03 * tc[1]
            -1.79862335e-07 * tc[2]
            +1.41876933e-11 * tc[3]
            -4.55120040e-16 * tc[4]
            +7.10667307e+00 ;
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 10: HO2 */
        species[10] =
            +3.18310656e+00 * tc[0]
            +3.66767950e-03 * tc[1]
            -4.66192561e-07 * tc[2]
            -1.08617640e-10 * tc[3]
            +3.77849780e-14 * tc[4]
            +8.39371099e+00 ;
        /*species 35: CH3O2H */
        species[35] =
            +3.23442817e+00 * tc[0]
            +1.90129767e-02 * tc[1]
            -5.66931435e-06 * tc[2]
            +1.13435551e-09 * tc[3]
            -1.02957555e-13 * tc[4]
            +9.25623949e+00 ;
        /*species 46: C7H15-1 */
        species[46] =
            -4.99570406e-01 * tc[0]
            +8.08826467e-02 * tc[1]
            -2.50266377e-05 * tc[2]
            +5.21831027e-09 * tc[3]
            -4.91540568e-13 * tc[4]
            +3.46564011e+01 ;
        /*species 47: C7H15O2 */
        species[47] =
            +3.55252917e-01 * tc[0]
            +9.42381007e-02 * tc[1]
            -3.33477553e-05 * tc[2]
            +8.49318610e-09 * tc[3]
            -1.03302776e-12 * tc[4]
            +3.29433435e+01 ;
        /*species 48: C7H14OOHO2 */
        species[48] =
            +1.40912479e+00 * tc[0]
            +1.10720834e-01 * tc[1]
            -4.35589037e-05 * tc[2]
            +1.21220916e-08 * tc[3]
            -1.57087672e-12 * tc[4]
            +3.16223383e+01 ;
    } else {
        /*species 10: HO2 */
        species[10] =
            +4.10547423e+00 * tc[0]
            +2.38452835e-03 * tc[1]
            -4.03173995e-07 * tc[2]
            +4.13972410e-11 * tc[3]
            -1.79100027e-15 * tc[4]
            +3.12515836e+00 ;
        /*species 35: CH3O2H */
        species[35] =
            +8.43117091e+00 * tc[0]
            +8.06817909e-03 * tc[1]
            -1.38547461e-06 * tc[2]
            +1.43777414e-10 * tc[3]
            -6.26730365e-15 * tc[4]
            -1.91170629e+01 ;
        /*species 46: C7H15-1 */
        species[46] =
            +2.17940709e+01 * tc[0]
            +3.26280243e-02 * tc[1]
            -5.55691220e-06 * tc[2]
            +5.73557160e-10 * tc[3]
            -2.49091750e-14 * tc[4]
            -8.64954311e+01 ;
        /*species 47: C7H15O2 */
        species[47] =
            +2.60535640e+01 * tc[0]
            +3.43831688e-02 * tc[1]
            -5.93566235e-06 * tc[2]
            +6.18167470e-10 * tc[3]
            -2.70130997e-14 * tc[4]
            -1.05407563e+02 ;
        /*species 48: C7H14OOHO2 */
        species[48] =
            +3.28358296e+01 * tc[0]
            +3.32767285e-02 * tc[1]
            -5.79166425e-06 * tc[2]
            +6.06528793e-10 * tc[3]
            -2.66086048e-14 * tc[4]
            -1.35975710e+02 ;
    }

    /*species with midpoint at T=1384 kelvin */
    if (T < 1384) {
        /*species 13: C2H6 */
        species[13] =
            -2.52854344e-02 * tc[0]
            +2.40764754e-02 * tc[1]
            -5.59467360e-06 * tc[2]
            +6.94469670e-10 * tc[3]
            -1.32467154e-14 * tc[4]
            +2.11648750e+01 ;
    } else {
        /*species 13: C2H6 */
        species[13] =
            +6.10683385e+00 * tc[0]
            +1.29236361e-02 * tc[1]
            -2.21263598e-06 * tc[2]
            +2.29130575e-10 * tc[3]
            -9.97254330e-15 * tc[4]
            -1.30081250e+01 ;
    }

    /*species with midpoint at T=1400 kelvin */
    if (T < 1400) {
        /*species 22: C3H4-A */
        species[22] =
            +2.53983100e+00 * tc[0]
            +1.63343700e-02 * tc[1]
            -8.82475000e-07 * tc[2]
            -1.54912167e-09 * tc[3]
            +4.32282750e-13 * tc[4]
            +9.93570200e+00 ;
        /*species 23: C3H4-P */
        species[23] =
            +3.02973000e+00 * tc[0]
            +1.49896100e-02 * tc[1]
            -6.99250000e-07 * tc[2]
            -1.32320633e-09 * tc[3]
            +3.47054250e-13 * tc[4]
            +8.00459400e+00 ;
        /*species 42: CH3CHCO */
        species[42] =
            +1.48380119e+00 * tc[0]
            +3.22203013e-02 * tc[1]
            -1.35125016e-05 * tc[2]
            +4.01663880e-09 * tc[3]
            -5.45914828e-13 * tc[4]
            +1.71552068e+01 ;
    } else {
        /*species 22: C3H4-A */
        species[22] =
            +9.77625600e+00 * tc[0]
            +5.30213800e-03 * tc[1]
            -1.85055900e-07 * tc[2]
            -1.00879533e-10 * tc[3]
            +1.27239525e-14 * tc[4]
            -3.07706100e+01 ;
        /*species 23: C3H4-P */
        species[23] =
            +9.76810200e+00 * tc[0]
            +5.21915100e-03 * tc[1]
            -1.87657000e-07 * tc[2]
            -9.97397000e-11 * tc[3]
            +1.27696950e-14 * tc[4]
            -3.02067800e+01 ;
        /*species 42: CH3CHCO */
        species[42] =
            +1.00219123e+01 * tc[0]
            +9.56966300e-03 * tc[1]
            -1.63110822e-06 * tc[2]
            +1.68410569e-10 * tc[3]
            -7.31483142e-15 * tc[4]
            -2.77829973e+01 ;
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 24: C3H6 */
        species[24] =
            +3.94615444e-01 * tc[0]
            +2.89107662e-02 * tc[1]
            -7.74434040e-06 * tc[2]
            +1.29604736e-09 * tc[3]
            -8.44725880e-14 * tc[4]
            +2.19003736e+01 ;
    } else {
        /*species 24: C3H6 */
        species[24] =
            +8.01595958e+00 * tc[0]
            +1.37023634e-02 * tc[1]
            -2.33124867e-06 * tc[2]
            +2.40418134e-10 * tc[3]
            -1.04342532e-14 * tc[4]
            -2.00160668e+01 ;
    }

    /*species with midpoint at T=1398 kelvin */
    if (T < 1398) {
        /*species 25: C4H6 */
        species[25] =
            -1.43095121e+00 * tc[0]
            +4.78706062e-02 * tc[1]
            -2.07723400e-05 * tc[2]
            +6.38498507e-09 * tc[3]
            -8.92896267e-13 * tc[4]
            +2.90825833e+01 ;
    } else {
        /*species 25: C4H6 */
        species[25] =
            +1.11633789e+01 * tc[0]
            +1.37163965e-02 * tc[1]
            -2.34857892e-06 * tc[2]
            +2.43231279e-10 * tc[3]
            -1.05871551e-14 * tc[4]
            -3.69847949e+01 ;
    }

    /*species with midpoint at T=1392 kelvin */
    if (T < 1392) {
        /*species 26: C4H7 */
        species[26] =
            -3.50508352e-01 * tc[0]
            +4.26511243e-02 * tc[1]
            -1.45489687e-05 * tc[2]
            +3.51346380e-09 * tc[3]
            -4.00149635e-13 * tc[4]
            +2.67295696e+01 ;
        /*species 27: C4H8-1 */
        species[27] =
            -8.31372089e-01 * tc[0]
            +4.52580978e-02 * tc[1]
            -1.46829280e-05 * tc[2]
            +3.34068120e-09 * tc[3]
            -3.57979200e-13 * tc[4]
            +2.95084236e+01 ;
        /*species 31: C5H9 */
        species[31] =
            -1.38013950e+00 * tc[0]
            +5.57608487e-02 * tc[1]
            -1.85071964e-05 * tc[2]
            +4.22946337e-09 * tc[3]
            -4.46347087e-13 * tc[4]
            +3.26441304e+01 ;
        /*species 32: C5H10-1 */
        species[32] =
            -1.06223481e+00 * tc[0]
            +5.74218294e-02 * tc[1]
            -1.87243445e-05 * tc[2]
            +4.24549963e-09 * tc[3]
            -4.49024472e-13 * tc[4]
            +3.22739790e+01 ;
    } else {
        /*species 26: C4H7 */
        species[26] =
            +1.12103578e+01 * tc[0]
            +1.60483196e-02 * tc[1]
            -2.73251146e-06 * tc[2]
            +2.81980351e-10 * tc[3]
            -1.22443185e-14 * tc[4]
            -3.55593015e+01 ;
        /*species 27: C4H8-1 */
        species[27] =
            +1.13508668e+01 * tc[0]
            +1.80617877e-02 * tc[1]
            -3.08046515e-06 * tc[2]
            +3.18217653e-10 * tc[3]
            -1.38272410e-14 * tc[4]
            -3.64369438e+01 ;
        /*species 31: C5H9 */
        species[31] =
            +1.41860454e+01 * tc[0]
            +2.07128899e-02 * tc[1]
            -3.53480309e-06 * tc[2]
            +3.65357110e-10 * tc[3]
            -1.58830552e-14 * tc[4]
            -5.14501773e+01 ;
        /*species 32: C5H10-1 */
        species[32] =
            +1.45851539e+01 * tc[0]
            +2.24072471e-02 * tc[1]
            -3.81674012e-06 * tc[2]
            +3.93963220e-10 * tc[3]
            -1.71096285e-14 * tc[4]
            -5.23683936e+01 ;
    }

    /*species with midpoint at T=1395 kelvin */
    if (T < 1395) {
        /*species 28: PC4H9 */
        species[28] =
            -4.37779725e-01 * tc[0]
            +4.78972364e-02 * tc[1]
            -1.57011580e-05 * tc[2]
            +3.65954907e-09 * tc[3]
            -4.05026660e-13 * tc[4]
            +2.86852732e+01 ;
    } else {
        /*species 28: PC4H9 */
        species[28] =
            +1.21510082e+01 * tc[0]
            +1.94310717e-02 * tc[1]
            -3.30788975e-06 * tc[2]
            +3.41250453e-10 * tc[3]
            -1.48132427e-14 * tc[4]
            -3.93425822e+01 ;
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 29: CH3COCH2 */
        species[29] =
            +1.80339187e+00 * tc[0]
            +3.01407085e-02 * tc[1]
            -9.67527760e-06 * tc[2]
            +2.12733011e-09 * tc[3]
            -2.16525795e-13 * tc[4]
            +1.78046408e+01 ;
        /*species 45: NC7H16 */
        species[45] =
            -1.26836187e+00 * tc[0]
            +8.54355820e-02 * tc[1]
            -2.62673393e-05 * tc[2]
            +5.43152403e-09 * tc[3]
            -5.05987313e-13 * tc[4]
            +3.53732912e+01 ;
    } else {
        /*species 29: CH3COCH2 */
        species[29] =
            +1.02303975e+01 * tc[0]
            +1.16494161e-02 * tc[1]
            -2.00502768e-06 * tc[2]
            +2.08401749e-10 * tc[3]
            -9.09460905e-15 * tc[4]
            -2.79197220e+01 ;
        /*species 45: NC7H16 */
        species[45] =
            +2.22148969e+01 * tc[0]
            +3.47675750e-02 * tc[1]
            -5.92035645e-06 * tc[2]
            +6.10994927e-10 * tc[3]
            -2.65325665e-14 * tc[4]
            -9.23040196e+01 ;
    }

    /*species with midpoint at T=1378 kelvin */
    if (T < 1378) {
        /*species 30: C2H5CHO */
        species[30] =
            +2.16308444e+00 * tc[0]
            +2.95501264e-02 * tc[1]
            -7.62231260e-06 * tc[2]
            +1.16501316e-09 * tc[3]
            -5.97241567e-14 * tc[4]
            +1.61153348e+01 ;
        /*species 40: NC3H7CHO */
        species[40] =
            +1.87415959e+00 * tc[0]
            +4.19240315e-02 * tc[1]
            -1.17574389e-05 * tc[2]
            +2.08971224e-09 * tc[3]
            -1.52360977e-13 * tc[4]
            +1.91568574e+01 ;
    } else {
        /*species 30: C2H5CHO */
        species[30] =
            +1.02427695e+01 * tc[0]
            +1.39641989e-02 * tc[1]
            -2.38124000e-06 * tc[2]
            +2.46035235e-10 * tc[3]
            -1.06939876e-14 * tc[4]
            -2.85346843e+01 ;
        /*species 40: NC3H7CHO */
        species[40] =
            +1.35988068e+01 * tc[0]
            +1.81652474e-02 * tc[1]
            -3.08922229e-06 * tc[2]
            +3.18660069e-10 * tc[3]
            -1.38360739e-14 * tc[4]
            -4.51790228e+01 ;
    }

    /*species with midpoint at T=1397 kelvin */
    if (T < 1397) {
        /*species 33: C5H11-1 */
        species[33] =
            -9.05255912e-01 * tc[0]
            +6.10632852e-02 * tc[1]
            -2.04745912e-05 * tc[2]
            +4.86978233e-09 * tc[3]
            -5.47149037e-13 * tc[4]
            +3.25574963e+01 ;
        /*species 38: C3H5-A */
        species[38] =
            -5.29131958e-01 * tc[0]
            +3.34559100e-02 * tc[1]
            -1.26700514e-05 * tc[2]
            +3.42885847e-09 * tc[3]
            -4.33145850e-13 * tc[4]
            +2.53067131e+01 ;
    } else {
        /*species 33: C5H11-1 */
        species[33] =
            +1.53234740e+01 * tc[0]
            +2.39041200e-02 * tc[1]
            -4.07385810e-06 * tc[2]
            +4.20587453e-10 * tc[3]
            -1.82669334e-14 * tc[4]
            -5.49528859e+01 ;
        /*species 38: C3H5-A */
        species[38] =
            +8.45883958e+00 * tc[0]
            +1.12695483e-02 * tc[1]
            -1.91896432e-06 * tc[2]
            +1.98019706e-10 * tc[3]
            -8.59795075e-15 * tc[4]
            -2.25809450e+01 ;
    }

    /*species with midpoint at T=1385 kelvin */
    if (T < 1385) {
        /*species 34: CH3O2 */
        species[34] =
            +4.26146906e+00 * tc[0]
            +1.00873599e-02 * tc[1]
            -1.60753092e-06 * tc[2]
            +6.98030890e-11 * tc[3]
            +1.04584776e-14 * tc[4]
            +5.16330320e+00 ;
    } else {
        /*species 34: CH3O2 */
        species[34] =
            +5.95787891e+00 * tc[0]
            +7.90728626e-03 * tc[1]
            -1.34123117e-06 * tc[2]
            +1.37963779e-10 * tc[3]
            -5.97518325e-15 * tc[4]
            -4.71963886e+00 ;
    }

    /*species with midpoint at T=1402 kelvin */
    if (T < 1402) {
        /*species 36: C2H3CO */
        species[36] =
            +1.36242013e+00 * tc[0]
            +3.15273972e-02 * tc[1]
            -1.50109468e-05 * tc[2]
            +4.93890373e-09 * tc[3]
            -7.19928825e-13 * tc[4]
            +1.72626546e+01 ;
    } else {
        /*species 36: C2H3CO */
        species[36] =
            +9.37467676e+00 * tc[0]
            +7.91296900e-03 * tc[1]
            -1.33599140e-06 * tc[2]
            +1.37038477e-10 * tc[3]
            -5.92447452e-15 * tc[4]
            -2.40892696e+01 ;
    }

    /*species with midpoint at T=1393 kelvin */
    if (T < 1393) {
        /*species 37: C2H3CHO */
        species[37] =
            +2.92355162e-01 * tc[0]
            +3.54321417e-02 * tc[1]
            -1.47468162e-05 * tc[2]
            +4.27000413e-09 * tc[3]
            -5.65360270e-13 * tc[4]
            +2.28878280e+01 ;
    } else {
        /*species 37: C2H3CHO */
        species[37] =
            +1.04184959e+01 * tc[0]
            +9.48963321e-03 * tc[1]
            -1.64655265e-06 * tc[2]
            +1.72093068e-10 * tc[3]
            -7.53968227e-15 * tc[4]
            -3.07235061e+01 ;
    }

    /*species with midpoint at T=1383 kelvin */
    if (T < 1383) {
        /*species 41: C2H5COCH2 */
        species[41] =
            +1.54013856e+00 * tc[0]
            +4.39486258e-02 * tc[1]
            -1.48501211e-05 * tc[2]
            +3.51651043e-09 * tc[3]
            -3.96496922e-13 * tc[4]
            +1.99706641e+01 ;
        /*species 43: NC3H7COCH3 */
        species[43] =
            +1.24037361e+00 * tc[0]
            +5.42494047e-02 * tc[1]
            -1.57412666e-05 * tc[2]
            +3.04830774e-09 * tc[3]
            -2.69163915e-13 * tc[4]
            +2.36150068e+01 ;
    } else {
        /*species 41: C2H5COCH2 */
        species[41] =
            +1.42098738e+01 * tc[0]
            +1.57866459e-02 * tc[1]
            -2.75264591e-06 * tc[2]
            +2.88623513e-10 * tc[3]
            -1.26728332e-14 * tc[4]
            -4.87132911e+01 ;
        /*species 43: NC3H7COCH3 */
        species[43] =
            +1.61841564e+01 * tc[0]
            +2.35245881e-02 * tc[1]
            -4.04837528e-06 * tc[2]
            +4.20768967e-10 * tc[3]
            -1.83621960e-14 * tc[4]
            -5.82980615e+01 ;
    }

    /*species with midpoint at T=1389 kelvin */
    if (T < 1389) {
        /*species 44: NC3H7COCH2 */
        species[44] =
            +9.58299271e-01 * tc[0]
            +5.68162532e-02 * tc[1]
            -1.99556391e-05 * tc[2]
            +5.08905047e-09 * tc[3]
            -6.23052617e-13 * tc[4]
            +2.34112884e+01 ;
    } else {
        /*species 44: NC3H7COCH2 */
        species[44] =
            +1.61502419e+01 * tc[0]
            +2.14093466e-02 * tc[1]
            -3.68029807e-06 * tc[2]
            +3.82189387e-10 * tc[3]
            -1.66678180e-14 * tc[4]
            -5.83864884e+01 ;
    }

    /*species with midpoint at T=1399 kelvin */
    if (T < 1399) {
        /*species 49: C7H14O */
        species[49] =
            -5.39477206e+00 * tc[0]
            +1.01846577e-01 * tc[1]
            -3.80023118e-05 * tc[2]
            +9.88460067e-09 * tc[3]
            -1.18713592e-12 * tc[4]
            +5.51731319e+01 ;
    } else {
        /*species 49: C7H14O */
        species[49] =
            +2.32692731e+01 * tc[0]
            +3.25585391e-02 * tc[1]
            -5.58124530e-06 * tc[2]
            +5.78578860e-10 * tc[3]
            -2.52033037e-14 * tc[4]
            -9.79457278e+01 ;
    }

    /*species with midpoint at T=1387 kelvin */
    if (T < 1387) {
        /*species 50: NC7KET */
        species[50] =
            +1.52936692e+00 * tc[0]
            +9.58173466e-02 * tc[1]
            -3.48344260e-05 * tc[2]
            +8.98467940e-09 * tc[3]
            -1.09682032e-12 * tc[4]
            +2.93851334e+01 ;
    } else {
        /*species 50: NC7KET */
        species[50] =
            +2.88330388e+01 * tc[0]
            +3.20168096e-02 * tc[1]
            -5.57542280e-06 * tc[2]
            +5.84087197e-10 * tc[3]
            -2.56301128e-14 * tc[4]
            -1.17499704e+02 ;
    }
    return;
}


/*Returns the mean specific heat at CP (Eq. 33) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPBL(amrex::Real *  T, amrex::Real *  x,  amrex::Real *  cpbl)
{
    int id; /*loop counter */
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real cpor[52]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 52; ++id) {
        result += x[id]*cpor[id];
    }

    *cpbl = result * 8.31446261815324e+07;
}


/*Returns the mean specific heat at CP (Eq. 34) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPBS(amrex::Real *  T, amrex::Real *  y,  amrex::Real *  cpbs)
{
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real cpor[52], tresult[52]; /* temporary storage */
    amrex::Real imw[52];

    get_imw(imw);

    cp_R(cpor, tc);

    for (int i = 0; i < 52; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 52; i++)
    {
        result += tresult[i];
    }

    *cpbs = result * 8.31446261815324e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVBL(amrex::Real *  T, amrex::Real *  x,  amrex::Real *  cvbl)
{
    int id; /*loop counter */
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real cvor[52]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 52; ++id) {
        result += x[id]*cvor[id];
    }

    *cvbl = result * 8.31446261815324e+07;
}


/*Returns the mean specific heat at CV (Eq. 36) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVBS(amrex::Real *  T, amrex::Real *  y,  amrex::Real *  cvbs)
{
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real cvor[52]; /* temporary storage */
    amrex::Real imw[52];

    get_imw(imw);

    cv_R(cvor, tc);

    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*H */
    result += cvor[1]*y[1]*imw[1]; /*H2 */
    result += cvor[2]*y[2]*imw[2]; /*O */
    result += cvor[3]*y[3]*imw[3]; /*O2 */
    result += cvor[4]*y[4]*imw[4]; /*OH */
    result += cvor[5]*y[5]*imw[5]; /*H2O */
    result += cvor[6]*y[6]*imw[6]; /*CO */
    result += cvor[7]*y[7]*imw[7]; /*CO2 */
    result += cvor[8]*y[8]*imw[8]; /*CH3 */
    result += cvor[9]*y[9]*imw[9]; /*CH4 */
    result += cvor[10]*y[10]*imw[10]; /*HO2 */
    result += cvor[11]*y[11]*imw[11]; /*H2O2 */
    result += cvor[12]*y[12]*imw[12]; /*CH2O */
    result += cvor[13]*y[13]*imw[13]; /*C2H6 */
    result += cvor[14]*y[14]*imw[14]; /*C2H4 */
    result += cvor[15]*y[15]*imw[15]; /*C2H5 */
    result += cvor[16]*y[16]*imw[16]; /*C2H */
    result += cvor[17]*y[17]*imw[17]; /*C2H2 */
    result += cvor[18]*y[18]*imw[18]; /*CH3OH */
    result += cvor[19]*y[19]*imw[19]; /*CH2CO */
    result += cvor[20]*y[20]*imw[20]; /*HCCO */
    result += cvor[21]*y[21]*imw[21]; /*CH3CHO */
    result += cvor[22]*y[22]*imw[22]; /*C3H4-A */
    result += cvor[23]*y[23]*imw[23]; /*C3H4-P */
    result += cvor[24]*y[24]*imw[24]; /*C3H6 */
    result += cvor[25]*y[25]*imw[25]; /*C4H6 */
    result += cvor[26]*y[26]*imw[26]; /*C4H7 */
    result += cvor[27]*y[27]*imw[27]; /*C4H8-1 */
    result += cvor[28]*y[28]*imw[28]; /*PC4H9 */
    result += cvor[29]*y[29]*imw[29]; /*CH3COCH2 */
    result += cvor[30]*y[30]*imw[30]; /*C2H5CHO */
    result += cvor[31]*y[31]*imw[31]; /*C5H9 */
    result += cvor[32]*y[32]*imw[32]; /*C5H10-1 */
    result += cvor[33]*y[33]*imw[33]; /*C5H11-1 */
    result += cvor[34]*y[34]*imw[34]; /*CH3O2 */
    result += cvor[35]*y[35]*imw[35]; /*CH3O2H */
    result += cvor[36]*y[36]*imw[36]; /*C2H3CO */
    result += cvor[37]*y[37]*imw[37]; /*C2H3CHO */
    result += cvor[38]*y[38]*imw[38]; /*C3H5-A */
    result += cvor[39]*y[39]*imw[39]; /*C3H3 */
    result += cvor[40]*y[40]*imw[40]; /*NC3H7CHO */
    result += cvor[41]*y[41]*imw[41]; /*C2H5COCH2 */
    result += cvor[42]*y[42]*imw[42]; /*CH3CHCO */
    result += cvor[43]*y[43]*imw[43]; /*NC3H7COCH3 */
    result += cvor[44]*y[44]*imw[44]; /*NC3H7COCH2 */
    result += cvor[45]*y[45]*imw[45]; /*NC7H16 */
    result += cvor[46]*y[46]*imw[46]; /*C7H15-1 */
    result += cvor[47]*y[47]*imw[47]; /*C7H15O2 */
    result += cvor[48]*y[48]*imw[48]; /*C7H14OOHO2 */
    result += cvor[49]*y[49]*imw[49]; /*C7H14O */
    result += cvor[50]*y[50]*imw[50]; /*NC7KET */
    result += cvor[51]*y[51]*imw[51]; /*N2 */

    *cvbs = result * 8.31446261815324e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHBML(amrex::Real *  T, amrex::Real *  x,  amrex::Real *  hbml)
{
    int id; /*loop counter */
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real hml[52]; /* temporary storage */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 52; ++id) {
        result += x[id]*hml[id];
    }

    *hbml = result * RT;
}


/*Returns mean enthalpy of mixture in mass units */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHBMS(amrex::Real *  T, amrex::Real *  y,  amrex::Real *  hbms)
{
    amrex::Real result = 0;
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real hml[52], tmp[52]; /* temporary storage */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[52];

    get_imw(imw);

    speciesEnthalpy(hml, tc);

    int id;
    for (id = 0; id < 52; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 52; ++id) {
        result += tmp[id];
    }

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUBML(amrex::Real *  T, amrex::Real *  x,  amrex::Real *  ubml)
{
    int id; /*loop counter */
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real uml[52]; /* temporary energy array */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 52; ++id) {
        result += x[id]*uml[id];
    }

    *ubml = result * RT;
}


/*get mean internal energy in mass units */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUBMS(amrex::Real *  T, amrex::Real *  y,  amrex::Real *  ubms)
{
    amrex::Real result = 0;
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real ums[52]; /* temporary energy array */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[52];

    get_imw(imw);

    speciesInternalEnergy(ums, tc);

    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*H */
    result += y[1]*ums[1]*imw[1]; /*H2 */
    result += y[2]*ums[2]*imw[2]; /*O */
    result += y[3]*ums[3]*imw[3]; /*O2 */
    result += y[4]*ums[4]*imw[4]; /*OH */
    result += y[5]*ums[5]*imw[5]; /*H2O */
    result += y[6]*ums[6]*imw[6]; /*CO */
    result += y[7]*ums[7]*imw[7]; /*CO2 */
    result += y[8]*ums[8]*imw[8]; /*CH3 */
    result += y[9]*ums[9]*imw[9]; /*CH4 */
    result += y[10]*ums[10]*imw[10]; /*HO2 */
    result += y[11]*ums[11]*imw[11]; /*H2O2 */
    result += y[12]*ums[12]*imw[12]; /*CH2O */
    result += y[13]*ums[13]*imw[13]; /*C2H6 */
    result += y[14]*ums[14]*imw[14]; /*C2H4 */
    result += y[15]*ums[15]*imw[15]; /*C2H5 */
    result += y[16]*ums[16]*imw[16]; /*C2H */
    result += y[17]*ums[17]*imw[17]; /*C2H2 */
    result += y[18]*ums[18]*imw[18]; /*CH3OH */
    result += y[19]*ums[19]*imw[19]; /*CH2CO */
    result += y[20]*ums[20]*imw[20]; /*HCCO */
    result += y[21]*ums[21]*imw[21]; /*CH3CHO */
    result += y[22]*ums[22]*imw[22]; /*C3H4-A */
    result += y[23]*ums[23]*imw[23]; /*C3H4-P */
    result += y[24]*ums[24]*imw[24]; /*C3H6 */
    result += y[25]*ums[25]*imw[25]; /*C4H6 */
    result += y[26]*ums[26]*imw[26]; /*C4H7 */
    result += y[27]*ums[27]*imw[27]; /*C4H8-1 */
    result += y[28]*ums[28]*imw[28]; /*PC4H9 */
    result += y[29]*ums[29]*imw[29]; /*CH3COCH2 */
    result += y[30]*ums[30]*imw[30]; /*C2H5CHO */
    result += y[31]*ums[31]*imw[31]; /*C5H9 */
    result += y[32]*ums[32]*imw[32]; /*C5H10-1 */
    result += y[33]*ums[33]*imw[33]; /*C5H11-1 */
    result += y[34]*ums[34]*imw[34]; /*CH3O2 */
    result += y[35]*ums[35]*imw[35]; /*CH3O2H */
    result += y[36]*ums[36]*imw[36]; /*C2H3CO */
    result += y[37]*ums[37]*imw[37]; /*C2H3CHO */
    result += y[38]*ums[38]*imw[38]; /*C3H5-A */
    result += y[39]*ums[39]*imw[39]; /*C3H3 */
    result += y[40]*ums[40]*imw[40]; /*NC3H7CHO */
    result += y[41]*ums[41]*imw[41]; /*C2H5COCH2 */
    result += y[42]*ums[42]*imw[42]; /*CH3CHCO */
    result += y[43]*ums[43]*imw[43]; /*NC3H7COCH3 */
    result += y[44]*ums[44]*imw[44]; /*NC3H7COCH2 */
    result += y[45]*ums[45]*imw[45]; /*NC7H16 */
    result += y[46]*ums[46]*imw[46]; /*C7H15-1 */
    result += y[47]*ums[47]*imw[47]; /*C7H15O2 */
    result += y[48]*ums[48]*imw[48]; /*C7H14OOHO2 */
    result += y[49]*ums[49]*imw[49]; /*C7H14O */
    result += y[50]*ums[50]*imw[50]; /*NC7KET */
    result += y[51]*ums[51]*imw[51]; /*N2 */

    *ubms = result * RT;
}


/*get mixture entropy in molar units */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSBML(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  sbml)
{
    int id; /*loop counter */
    amrex::Real result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    amrex::Real logPratio = log ( *P / 1013250.0 ); 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real sor[52]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 52; ++id) {
        result += x[id]*(sor[id]-log((x[id]+1e-100))-logPratio);
    }

    *sbml = result * 8.31446261815324e+07;
}


/*get mixture entropy in mass units */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSBMS(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  sbms)
{
    amrex::Real result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    amrex::Real logPratio = log ( *P / 1013250.0 ); 
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real sor[52]; /* temporary storage */
    amrex::Real x[52]; /* need a ytx conversion */
    amrex::Real YOW = 0; /*See Eq 4, 6 in CK Manual */
    amrex::Real imw[52];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*O2 */
    YOW += y[4]*imw[4]; /*OH */
    YOW += y[5]*imw[5]; /*H2O */
    YOW += y[6]*imw[6]; /*CO */
    YOW += y[7]*imw[7]; /*CO2 */
    YOW += y[8]*imw[8]; /*CH3 */
    YOW += y[9]*imw[9]; /*CH4 */
    YOW += y[10]*imw[10]; /*HO2 */
    YOW += y[11]*imw[11]; /*H2O2 */
    YOW += y[12]*imw[12]; /*CH2O */
    YOW += y[13]*imw[13]; /*C2H6 */
    YOW += y[14]*imw[14]; /*C2H4 */
    YOW += y[15]*imw[15]; /*C2H5 */
    YOW += y[16]*imw[16]; /*C2H */
    YOW += y[17]*imw[17]; /*C2H2 */
    YOW += y[18]*imw[18]; /*CH3OH */
    YOW += y[19]*imw[19]; /*CH2CO */
    YOW += y[20]*imw[20]; /*HCCO */
    YOW += y[21]*imw[21]; /*CH3CHO */
    YOW += y[22]*imw[22]; /*C3H4-A */
    YOW += y[23]*imw[23]; /*C3H4-P */
    YOW += y[24]*imw[24]; /*C3H6 */
    YOW += y[25]*imw[25]; /*C4H6 */
    YOW += y[26]*imw[26]; /*C4H7 */
    YOW += y[27]*imw[27]; /*C4H8-1 */
    YOW += y[28]*imw[28]; /*PC4H9 */
    YOW += y[29]*imw[29]; /*CH3COCH2 */
    YOW += y[30]*imw[30]; /*C2H5CHO */
    YOW += y[31]*imw[31]; /*C5H9 */
    YOW += y[32]*imw[32]; /*C5H10-1 */
    YOW += y[33]*imw[33]; /*C5H11-1 */
    YOW += y[34]*imw[34]; /*CH3O2 */
    YOW += y[35]*imw[35]; /*CH3O2H */
    YOW += y[36]*imw[36]; /*C2H3CO */
    YOW += y[37]*imw[37]; /*C2H3CHO */
    YOW += y[38]*imw[38]; /*C3H5-A */
    YOW += y[39]*imw[39]; /*C3H3 */
    YOW += y[40]*imw[40]; /*NC3H7CHO */
    YOW += y[41]*imw[41]; /*C2H5COCH2 */
    YOW += y[42]*imw[42]; /*CH3CHCO */
    YOW += y[43]*imw[43]; /*NC3H7COCH3 */
    YOW += y[44]*imw[44]; /*NC3H7COCH2 */
    YOW += y[45]*imw[45]; /*NC7H16 */
    YOW += y[46]*imw[46]; /*C7H15-1 */
    YOW += y[47]*imw[47]; /*C7H15O2 */
    YOW += y[48]*imw[48]; /*C7H14OOHO2 */
    YOW += y[49]*imw[49]; /*C7H14O */
    YOW += y[50]*imw[50]; /*NC7KET */
    YOW += y[51]*imw[51]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(1.007970*YOW); 
    x[1] = y[1]/(2.015940*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(28.010550*YOW); 
    x[7] = y[7]/(44.009950*YOW); 
    x[8] = y[8]/(15.035060*YOW); 
    x[9] = y[9]/(16.043030*YOW); 
    x[10] = y[10]/(33.006770*YOW); 
    x[11] = y[11]/(34.014740*YOW); 
    x[12] = y[12]/(30.026490*YOW); 
    x[13] = y[13]/(30.070120*YOW); 
    x[14] = y[14]/(28.054180*YOW); 
    x[15] = y[15]/(29.062150*YOW); 
    x[16] = y[16]/(25.030270*YOW); 
    x[17] = y[17]/(26.038240*YOW); 
    x[18] = y[18]/(32.042430*YOW); 
    x[19] = y[19]/(42.037640*YOW); 
    x[20] = y[20]/(41.029670*YOW); 
    x[21] = y[21]/(44.053580*YOW); 
    x[22] = y[22]/(40.065330*YOW); 
    x[23] = y[23]/(40.065330*YOW); 
    x[24] = y[24]/(42.081270*YOW); 
    x[25] = y[25]/(54.092420*YOW); 
    x[26] = y[26]/(55.100390*YOW); 
    x[27] = y[27]/(56.108360*YOW); 
    x[28] = y[28]/(57.116330*YOW); 
    x[29] = y[29]/(57.072700*YOW); 
    x[30] = y[30]/(58.080670*YOW); 
    x[31] = y[31]/(69.127480*YOW); 
    x[32] = y[32]/(70.135450*YOW); 
    x[33] = y[33]/(71.143420*YOW); 
    x[34] = y[34]/(47.033860*YOW); 
    x[35] = y[35]/(48.041830*YOW); 
    x[36] = y[36]/(55.056760*YOW); 
    x[37] = y[37]/(56.064730*YOW); 
    x[38] = y[38]/(41.073300*YOW); 
    x[39] = y[39]/(39.057360*YOW); 
    x[40] = y[40]/(72.107760*YOW); 
    x[41] = y[41]/(71.099790*YOW); 
    x[42] = y[42]/(56.064730*YOW); 
    x[43] = y[43]/(86.134850*YOW); 
    x[44] = y[44]/(85.126880*YOW); 
    x[45] = y[45]/(100.205570*YOW); 
    x[46] = y[46]/(99.197600*YOW); 
    x[47] = y[47]/(131.196400*YOW); 
    x[48] = y[48]/(163.195200*YOW); 
    x[49] = y[49]/(114.189030*YOW); 
    x[50] = y[50]/(146.187830*YOW); 
    x[51] = y[51]/(28.013400*YOW); 
    speciesEntropy(sor, tc);
    /*Perform computation in Eq 42 and 43 */
    result += x[0]*(sor[0]-log((x[0]+1e-100))-logPratio);
    result += x[1]*(sor[1]-log((x[1]+1e-100))-logPratio);
    result += x[2]*(sor[2]-log((x[2]+1e-100))-logPratio);
    result += x[3]*(sor[3]-log((x[3]+1e-100))-logPratio);
    result += x[4]*(sor[4]-log((x[4]+1e-100))-logPratio);
    result += x[5]*(sor[5]-log((x[5]+1e-100))-logPratio);
    result += x[6]*(sor[6]-log((x[6]+1e-100))-logPratio);
    result += x[7]*(sor[7]-log((x[7]+1e-100))-logPratio);
    result += x[8]*(sor[8]-log((x[8]+1e-100))-logPratio);
    result += x[9]*(sor[9]-log((x[9]+1e-100))-logPratio);
    result += x[10]*(sor[10]-log((x[10]+1e-100))-logPratio);
    result += x[11]*(sor[11]-log((x[11]+1e-100))-logPratio);
    result += x[12]*(sor[12]-log((x[12]+1e-100))-logPratio);
    result += x[13]*(sor[13]-log((x[13]+1e-100))-logPratio);
    result += x[14]*(sor[14]-log((x[14]+1e-100))-logPratio);
    result += x[15]*(sor[15]-log((x[15]+1e-100))-logPratio);
    result += x[16]*(sor[16]-log((x[16]+1e-100))-logPratio);
    result += x[17]*(sor[17]-log((x[17]+1e-100))-logPratio);
    result += x[18]*(sor[18]-log((x[18]+1e-100))-logPratio);
    result += x[19]*(sor[19]-log((x[19]+1e-100))-logPratio);
    result += x[20]*(sor[20]-log((x[20]+1e-100))-logPratio);
    result += x[21]*(sor[21]-log((x[21]+1e-100))-logPratio);
    result += x[22]*(sor[22]-log((x[22]+1e-100))-logPratio);
    result += x[23]*(sor[23]-log((x[23]+1e-100))-logPratio);
    result += x[24]*(sor[24]-log((x[24]+1e-100))-logPratio);
    result += x[25]*(sor[25]-log((x[25]+1e-100))-logPratio);
    result += x[26]*(sor[26]-log((x[26]+1e-100))-logPratio);
    result += x[27]*(sor[27]-log((x[27]+1e-100))-logPratio);
    result += x[28]*(sor[28]-log((x[28]+1e-100))-logPratio);
    result += x[29]*(sor[29]-log((x[29]+1e-100))-logPratio);
    result += x[30]*(sor[30]-log((x[30]+1e-100))-logPratio);
    result += x[31]*(sor[31]-log((x[31]+1e-100))-logPratio);
    result += x[32]*(sor[32]-log((x[32]+1e-100))-logPratio);
    result += x[33]*(sor[33]-log((x[33]+1e-100))-logPratio);
    result += x[34]*(sor[34]-log((x[34]+1e-100))-logPratio);
    result += x[35]*(sor[35]-log((x[35]+1e-100))-logPratio);
    result += x[36]*(sor[36]-log((x[36]+1e-100))-logPratio);
    result += x[37]*(sor[37]-log((x[37]+1e-100))-logPratio);
    result += x[38]*(sor[38]-log((x[38]+1e-100))-logPratio);
    result += x[39]*(sor[39]-log((x[39]+1e-100))-logPratio);
    result += x[40]*(sor[40]-log((x[40]+1e-100))-logPratio);
    result += x[41]*(sor[41]-log((x[41]+1e-100))-logPratio);
    result += x[42]*(sor[42]-log((x[42]+1e-100))-logPratio);
    result += x[43]*(sor[43]-log((x[43]+1e-100))-logPratio);
    result += x[44]*(sor[44]-log((x[44]+1e-100))-logPratio);
    result += x[45]*(sor[45]-log((x[45]+1e-100))-logPratio);
    result += x[46]*(sor[46]-log((x[46]+1e-100))-logPratio);
    result += x[47]*(sor[47]-log((x[47]+1e-100))-logPratio);
    result += x[48]*(sor[48]-log((x[48]+1e-100))-logPratio);
    result += x[49]*(sor[49]-log((x[49]+1e-100))-logPratio);
    result += x[50]*(sor[50]-log((x[50]+1e-100))-logPratio);
    result += x[51]*(sor[51]-log((x[51]+1e-100))-logPratio);
    /*Scale by R/W */
    *sbms = result * 8.31446261815324e+07 * YOW;
}

/* get temperature given internal energy in mass units and mass fracs */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void GET_T_GIVEN_EY(amrex::Real *  e, amrex::Real *  y, amrex::Real *  t, int * ierr)
{
#ifdef CONVERGENCE
    const int maxiter = 5000;
    const amrex::Real tol  = 1.e-12;
#else
    const int maxiter = 200;
    const amrex::Real tol  = 1.e-6;
#endif
    amrex::Real ein  = *e;
    amrex::Real tmin = 90;/*max lower bound for thermo def */
    amrex::Real tmax = 4000;/*min upper bound for thermo def */
    amrex::Real e1,emin,emax,cv,t1,dt;
    int i;/* loop counter */
    CKUBMS(&tmin, y, &emin);
    CKUBMS(&tmax, y, &emax);
    if (ein < emin) {
        /*Linear Extrapolation below tmin */
        CKCVBS(&tmin, y, &cv);
        *t = tmin - (emin-ein)/cv;
        *ierr = 1;
        return;
    }
    if (ein > emax) {
        /*Linear Extrapolation above tmax */
        CKCVBS(&tmax, y, &cv);
        *t = tmax - (emax-ein)/cv;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin);
    }
    for (i = 0; i < maxiter; ++i) {
        CKUBMS(&t1,y,&e1);
        CKCVBS(&t1,y,&cv);
        dt = (ein - e1) / cv;
        if (dt > 100.) { dt = 100.; }
        else if (dt < -100.) { dt = -100.; }
        else if (fabs(dt) < tol) break;
        else if (t1+dt == t1) break;
        t1 += dt;
    }
    *t = t1;
    *ierr = 0;
    return;
}

/* get temperature given enthalpy in mass units and mass fracs */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void GET_T_GIVEN_HY(amrex::Real *  h, amrex::Real *  y, amrex::Real *  t, int * ierr)
{
#ifdef CONVERGENCE
    const int maxiter = 5000;
    const amrex::Real tol  = 1.e-12;
#else
    const int maxiter = 200;
    const amrex::Real tol  = 1.e-6;
#endif
    amrex::Real hin  = *h;
    amrex::Real tmin = 90;/*max lower bound for thermo def */
    amrex::Real tmax = 4000;/*min upper bound for thermo def */
    amrex::Real h1,hmin,hmax,cp,t1,dt;
    int i;/* loop counter */
    CKHBMS(&tmin, y, &hmin);
    CKHBMS(&tmax, y, &hmax);
    if (hin < hmin) {
        /*Linear Extrapolation below tmin */
        CKCPBS(&tmin, y, &cp);
        *t = tmin - (hmin-hin)/cp;
        *ierr = 1;
        return;
    }
    if (hin > hmax) {
        /*Linear Extrapolation above tmax */
        CKCPBS(&tmax, y, &cp);
        *t = tmax - (hmax-hin)/cp;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(hmax-hmin)*(hin-hmin);
    }
    for (i = 0; i < maxiter; ++i) {
        CKHBMS(&t1,y,&h1);
        CKCPBS(&t1,y,&cp);
        dt = (hin - h1) / cp;
        if (dt > 100.) { dt = 100.; }
        else if (dt < -100.) { dt = -100.; }
        else if (fabs(dt) < tol) break;
        else if (t1+dt == t1) break;
        t1 += dt;
    }
    *t = t1;
    *ierr = 0;
    return;
}


/*Compute P = rhoRT/W(x) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPX(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x, amrex::Real *  P)
{
    amrex::Real XW = 0;/* To hold mean molecular wt */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*28.010550; /*CO */
    XW += x[7]*44.009950; /*CO2 */
    XW += x[8]*15.035060; /*CH3 */
    XW += x[9]*16.043030; /*CH4 */
    XW += x[10]*33.006770; /*HO2 */
    XW += x[11]*34.014740; /*H2O2 */
    XW += x[12]*30.026490; /*CH2O */
    XW += x[13]*30.070120; /*C2H6 */
    XW += x[14]*28.054180; /*C2H4 */
    XW += x[15]*29.062150; /*C2H5 */
    XW += x[16]*25.030270; /*C2H */
    XW += x[17]*26.038240; /*C2H2 */
    XW += x[18]*32.042430; /*CH3OH */
    XW += x[19]*42.037640; /*CH2CO */
    XW += x[20]*41.029670; /*HCCO */
    XW += x[21]*44.053580; /*CH3CHO */
    XW += x[22]*40.065330; /*C3H4-A */
    XW += x[23]*40.065330; /*C3H4-P */
    XW += x[24]*42.081270; /*C3H6 */
    XW += x[25]*54.092420; /*C4H6 */
    XW += x[26]*55.100390; /*C4H7 */
    XW += x[27]*56.108360; /*C4H8-1 */
    XW += x[28]*57.116330; /*PC4H9 */
    XW += x[29]*57.072700; /*CH3COCH2 */
    XW += x[30]*58.080670; /*C2H5CHO */
    XW += x[31]*69.127480; /*C5H9 */
    XW += x[32]*70.135450; /*C5H10-1 */
    XW += x[33]*71.143420; /*C5H11-1 */
    XW += x[34]*47.033860; /*CH3O2 */
    XW += x[35]*48.041830; /*CH3O2H */
    XW += x[36]*55.056760; /*C2H3CO */
    XW += x[37]*56.064730; /*C2H3CHO */
    XW += x[38]*41.073300; /*C3H5-A */
    XW += x[39]*39.057360; /*C3H3 */
    XW += x[40]*72.107760; /*NC3H7CHO */
    XW += x[41]*71.099790; /*C2H5COCH2 */
    XW += x[42]*56.064730; /*CH3CHCO */
    XW += x[43]*86.134850; /*NC3H7COCH3 */
    XW += x[44]*85.126880; /*NC3H7COCH2 */
    XW += x[45]*100.205570; /*NC7H16 */
    XW += x[46]*99.197600; /*C7H15-1 */
    XW += x[47]*131.196400; /*C7H15O2 */
    XW += x[48]*163.195200; /*C7H14OOHO2 */
    XW += x[49]*114.189030; /*C7H14O */
    XW += x[50]*146.187830; /*NC7KET */
    XW += x[51]*28.013400; /*N2 */
    *P = *rho * 8.31446261815324e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPY(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  P)
{
    amrex::Real YOW = 0;/* for computing mean MW */
    amrex::Real imw[52];

    get_imw(imw);

    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*O2 */
    YOW += y[4]*imw[4]; /*OH */
    YOW += y[5]*imw[5]; /*H2O */
    YOW += y[6]*imw[6]; /*CO */
    YOW += y[7]*imw[7]; /*CO2 */
    YOW += y[8]*imw[8]; /*CH3 */
    YOW += y[9]*imw[9]; /*CH4 */
    YOW += y[10]*imw[10]; /*HO2 */
    YOW += y[11]*imw[11]; /*H2O2 */
    YOW += y[12]*imw[12]; /*CH2O */
    YOW += y[13]*imw[13]; /*C2H6 */
    YOW += y[14]*imw[14]; /*C2H4 */
    YOW += y[15]*imw[15]; /*C2H5 */
    YOW += y[16]*imw[16]; /*C2H */
    YOW += y[17]*imw[17]; /*C2H2 */
    YOW += y[18]*imw[18]; /*CH3OH */
    YOW += y[19]*imw[19]; /*CH2CO */
    YOW += y[20]*imw[20]; /*HCCO */
    YOW += y[21]*imw[21]; /*CH3CHO */
    YOW += y[22]*imw[22]; /*C3H4-A */
    YOW += y[23]*imw[23]; /*C3H4-P */
    YOW += y[24]*imw[24]; /*C3H6 */
    YOW += y[25]*imw[25]; /*C4H6 */
    YOW += y[26]*imw[26]; /*C4H7 */
    YOW += y[27]*imw[27]; /*C4H8-1 */
    YOW += y[28]*imw[28]; /*PC4H9 */
    YOW += y[29]*imw[29]; /*CH3COCH2 */
    YOW += y[30]*imw[30]; /*C2H5CHO */
    YOW += y[31]*imw[31]; /*C5H9 */
    YOW += y[32]*imw[32]; /*C5H10-1 */
    YOW += y[33]*imw[33]; /*C5H11-1 */
    YOW += y[34]*imw[34]; /*CH3O2 */
    YOW += y[35]*imw[35]; /*CH3O2H */
    YOW += y[36]*imw[36]; /*C2H3CO */
    YOW += y[37]*imw[37]; /*C2H3CHO */
    YOW += y[38]*imw[38]; /*C3H5-A */
    YOW += y[39]*imw[39]; /*C3H3 */
    YOW += y[40]*imw[40]; /*NC3H7CHO */
    YOW += y[41]*imw[41]; /*C2H5COCH2 */
    YOW += y[42]*imw[42]; /*CH3CHCO */
    YOW += y[43]*imw[43]; /*NC3H7COCH3 */
    YOW += y[44]*imw[44]; /*NC3H7COCH2 */
    YOW += y[45]*imw[45]; /*NC7H16 */
    YOW += y[46]*imw[46]; /*C7H15-1 */
    YOW += y[47]*imw[47]; /*C7H15O2 */
    YOW += y[48]*imw[48]; /*C7H14OOHO2 */
    YOW += y[49]*imw[49]; /*C7H14O */
    YOW += y[50]*imw[50]; /*NC7KET */
    YOW += y[51]*imw[51]; /*N2 */
    *P = *rho * 8.31446261815324e+07 * (*T) * YOW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(c) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPC(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  c,  amrex::Real *  P)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    amrex::Real W = 0;
    amrex::Real sumC = 0;
    W += c[0]*1.007970; /*H */
    W += c[1]*2.015940; /*H2 */
    W += c[2]*15.999400; /*O */
    W += c[3]*31.998800; /*O2 */
    W += c[4]*17.007370; /*OH */
    W += c[5]*18.015340; /*H2O */
    W += c[6]*28.010550; /*CO */
    W += c[7]*44.009950; /*CO2 */
    W += c[8]*15.035060; /*CH3 */
    W += c[9]*16.043030; /*CH4 */
    W += c[10]*33.006770; /*HO2 */
    W += c[11]*34.014740; /*H2O2 */
    W += c[12]*30.026490; /*CH2O */
    W += c[13]*30.070120; /*C2H6 */
    W += c[14]*28.054180; /*C2H4 */
    W += c[15]*29.062150; /*C2H5 */
    W += c[16]*25.030270; /*C2H */
    W += c[17]*26.038240; /*C2H2 */
    W += c[18]*32.042430; /*CH3OH */
    W += c[19]*42.037640; /*CH2CO */
    W += c[20]*41.029670; /*HCCO */
    W += c[21]*44.053580; /*CH3CHO */
    W += c[22]*40.065330; /*C3H4-A */
    W += c[23]*40.065330; /*C3H4-P */
    W += c[24]*42.081270; /*C3H6 */
    W += c[25]*54.092420; /*C4H6 */
    W += c[26]*55.100390; /*C4H7 */
    W += c[27]*56.108360; /*C4H8-1 */
    W += c[28]*57.116330; /*PC4H9 */
    W += c[29]*57.072700; /*CH3COCH2 */
    W += c[30]*58.080670; /*C2H5CHO */
    W += c[31]*69.127480; /*C5H9 */
    W += c[32]*70.135450; /*C5H10-1 */
    W += c[33]*71.143420; /*C5H11-1 */
    W += c[34]*47.033860; /*CH3O2 */
    W += c[35]*48.041830; /*CH3O2H */
    W += c[36]*55.056760; /*C2H3CO */
    W += c[37]*56.064730; /*C2H3CHO */
    W += c[38]*41.073300; /*C3H5-A */
    W += c[39]*39.057360; /*C3H3 */
    W += c[40]*72.107760; /*NC3H7CHO */
    W += c[41]*71.099790; /*C2H5COCH2 */
    W += c[42]*56.064730; /*CH3CHCO */
    W += c[43]*86.134850; /*NC3H7COCH3 */
    W += c[44]*85.126880; /*NC3H7COCH2 */
    W += c[45]*100.205570; /*NC7H16 */
    W += c[46]*99.197600; /*C7H15-1 */
    W += c[47]*131.196400; /*C7H15O2 */
    W += c[48]*163.195200; /*C7H14OOHO2 */
    W += c[49]*114.189030; /*C7H14O */
    W += c[50]*146.187830; /*NC7KET */
    W += c[51]*28.013400; /*N2 */

    for (id = 0; id < 52; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31446261815324e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOX(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  rho)
{
    amrex::Real XW = 0;/* To hold mean molecular wt */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*28.010550; /*CO */
    XW += x[7]*44.009950; /*CO2 */
    XW += x[8]*15.035060; /*CH3 */
    XW += x[9]*16.043030; /*CH4 */
    XW += x[10]*33.006770; /*HO2 */
    XW += x[11]*34.014740; /*H2O2 */
    XW += x[12]*30.026490; /*CH2O */
    XW += x[13]*30.070120; /*C2H6 */
    XW += x[14]*28.054180; /*C2H4 */
    XW += x[15]*29.062150; /*C2H5 */
    XW += x[16]*25.030270; /*C2H */
    XW += x[17]*26.038240; /*C2H2 */
    XW += x[18]*32.042430; /*CH3OH */
    XW += x[19]*42.037640; /*CH2CO */
    XW += x[20]*41.029670; /*HCCO */
    XW += x[21]*44.053580; /*CH3CHO */
    XW += x[22]*40.065330; /*C3H4-A */
    XW += x[23]*40.065330; /*C3H4-P */
    XW += x[24]*42.081270; /*C3H6 */
    XW += x[25]*54.092420; /*C4H6 */
    XW += x[26]*55.100390; /*C4H7 */
    XW += x[27]*56.108360; /*C4H8-1 */
    XW += x[28]*57.116330; /*PC4H9 */
    XW += x[29]*57.072700; /*CH3COCH2 */
    XW += x[30]*58.080670; /*C2H5CHO */
    XW += x[31]*69.127480; /*C5H9 */
    XW += x[32]*70.135450; /*C5H10-1 */
    XW += x[33]*71.143420; /*C5H11-1 */
    XW += x[34]*47.033860; /*CH3O2 */
    XW += x[35]*48.041830; /*CH3O2H */
    XW += x[36]*55.056760; /*C2H3CO */
    XW += x[37]*56.064730; /*C2H3CHO */
    XW += x[38]*41.073300; /*C3H5-A */
    XW += x[39]*39.057360; /*C3H3 */
    XW += x[40]*72.107760; /*NC3H7CHO */
    XW += x[41]*71.099790; /*C2H5COCH2 */
    XW += x[42]*56.064730; /*CH3CHCO */
    XW += x[43]*86.134850; /*NC3H7COCH3 */
    XW += x[44]*85.126880; /*NC3H7COCH2 */
    XW += x[45]*100.205570; /*NC7H16 */
    XW += x[46]*99.197600; /*C7H15-1 */
    XW += x[47]*131.196400; /*C7H15O2 */
    XW += x[48]*163.195200; /*C7H14OOHO2 */
    XW += x[49]*114.189030; /*C7H14O */
    XW += x[50]*146.187830; /*NC7KET */
    XW += x[51]*28.013400; /*N2 */
    *rho = *P * XW / (8.31446261815324e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOY(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  rho)
{
    amrex::Real YOW = 0;
    amrex::Real tmp[52];
    amrex::Real imw[52];

    get_imw(imw);

    for (int i = 0; i < 52; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 52; i++)
    {
        YOW += tmp[i];
    }

    *rho = *P / (8.31446261815324e+07 * (*T) * YOW);/*rho = P*W/(R*T) */
    return;
}


/*Compute rho = P*W(c)/(R*T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOC(amrex::Real *  P, amrex::Real *  T, amrex::Real *  c,  amrex::Real *  rho)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    amrex::Real W = 0;
    amrex::Real sumC = 0;
    W += c[0]*1.007970; /*H */
    W += c[1]*2.015940; /*H2 */
    W += c[2]*15.999400; /*O */
    W += c[3]*31.998800; /*O2 */
    W += c[4]*17.007370; /*OH */
    W += c[5]*18.015340; /*H2O */
    W += c[6]*28.010550; /*CO */
    W += c[7]*44.009950; /*CO2 */
    W += c[8]*15.035060; /*CH3 */
    W += c[9]*16.043030; /*CH4 */
    W += c[10]*33.006770; /*HO2 */
    W += c[11]*34.014740; /*H2O2 */
    W += c[12]*30.026490; /*CH2O */
    W += c[13]*30.070120; /*C2H6 */
    W += c[14]*28.054180; /*C2H4 */
    W += c[15]*29.062150; /*C2H5 */
    W += c[16]*25.030270; /*C2H */
    W += c[17]*26.038240; /*C2H2 */
    W += c[18]*32.042430; /*CH3OH */
    W += c[19]*42.037640; /*CH2CO */
    W += c[20]*41.029670; /*HCCO */
    W += c[21]*44.053580; /*CH3CHO */
    W += c[22]*40.065330; /*C3H4-A */
    W += c[23]*40.065330; /*C3H4-P */
    W += c[24]*42.081270; /*C3H6 */
    W += c[25]*54.092420; /*C4H6 */
    W += c[26]*55.100390; /*C4H7 */
    W += c[27]*56.108360; /*C4H8-1 */
    W += c[28]*57.116330; /*PC4H9 */
    W += c[29]*57.072700; /*CH3COCH2 */
    W += c[30]*58.080670; /*C2H5CHO */
    W += c[31]*69.127480; /*C5H9 */
    W += c[32]*70.135450; /*C5H10-1 */
    W += c[33]*71.143420; /*C5H11-1 */
    W += c[34]*47.033860; /*CH3O2 */
    W += c[35]*48.041830; /*CH3O2H */
    W += c[36]*55.056760; /*C2H3CO */
    W += c[37]*56.064730; /*C2H3CHO */
    W += c[38]*41.073300; /*C3H5-A */
    W += c[39]*39.057360; /*C3H3 */
    W += c[40]*72.107760; /*NC3H7CHO */
    W += c[41]*71.099790; /*C2H5COCH2 */
    W += c[42]*56.064730; /*CH3CHCO */
    W += c[43]*86.134850; /*NC3H7COCH3 */
    W += c[44]*85.126880; /*NC3H7COCH2 */
    W += c[45]*100.205570; /*NC7H16 */
    W += c[46]*99.197600; /*C7H15-1 */
    W += c[47]*131.196400; /*C7H15O2 */
    W += c[48]*163.195200; /*C7H14OOHO2 */
    W += c[49]*114.189030; /*C7H14O */
    W += c[50]*146.187830; /*NC7KET */
    W += c[51]*28.013400; /*N2 */

    for (id = 0; id < 52; ++id) {
        sumC += c[id];
    }
    *rho = *P * W / (sumC * (*T) * 8.31446261815324e+07); /*rho = PW/(R*T) */

    return;
}


/*get molecular weight for all species */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWT( amrex::Real *  wt)
{
    get_mw(wt);
}

/*given y[species]: mass fractions */
/*returns mean molecular weight (gm/mole) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKMMWY(amrex::Real *  y,  amrex::Real *  wtm)
{
    amrex::Real YOW = 0;
    amrex::Real tmp[52];
    amrex::Real imw[52];

    get_imw(imw);

    for (int i = 0; i < 52; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 52; i++)
    {
        YOW += tmp[i];
    }

    *wtm = 1.0 / YOW;
    return;
}


/*given x[species]: mole fractions */
/*returns mean molecular weight (gm/mole) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKMMWX(amrex::Real *  x,  amrex::Real *  wtm)
{
    amrex::Real XW = 0;/* see Eq 4 in CK Manual */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*28.010550; /*CO */
    XW += x[7]*44.009950; /*CO2 */
    XW += x[8]*15.035060; /*CH3 */
    XW += x[9]*16.043030; /*CH4 */
    XW += x[10]*33.006770; /*HO2 */
    XW += x[11]*34.014740; /*H2O2 */
    XW += x[12]*30.026490; /*CH2O */
    XW += x[13]*30.070120; /*C2H6 */
    XW += x[14]*28.054180; /*C2H4 */
    XW += x[15]*29.062150; /*C2H5 */
    XW += x[16]*25.030270; /*C2H */
    XW += x[17]*26.038240; /*C2H2 */
    XW += x[18]*32.042430; /*CH3OH */
    XW += x[19]*42.037640; /*CH2CO */
    XW += x[20]*41.029670; /*HCCO */
    XW += x[21]*44.053580; /*CH3CHO */
    XW += x[22]*40.065330; /*C3H4-A */
    XW += x[23]*40.065330; /*C3H4-P */
    XW += x[24]*42.081270; /*C3H6 */
    XW += x[25]*54.092420; /*C4H6 */
    XW += x[26]*55.100390; /*C4H7 */
    XW += x[27]*56.108360; /*C4H8-1 */
    XW += x[28]*57.116330; /*PC4H9 */
    XW += x[29]*57.072700; /*CH3COCH2 */
    XW += x[30]*58.080670; /*C2H5CHO */
    XW += x[31]*69.127480; /*C5H9 */
    XW += x[32]*70.135450; /*C5H10-1 */
    XW += x[33]*71.143420; /*C5H11-1 */
    XW += x[34]*47.033860; /*CH3O2 */
    XW += x[35]*48.041830; /*CH3O2H */
    XW += x[36]*55.056760; /*C2H3CO */
    XW += x[37]*56.064730; /*C2H3CHO */
    XW += x[38]*41.073300; /*C3H5-A */
    XW += x[39]*39.057360; /*C3H3 */
    XW += x[40]*72.107760; /*NC3H7CHO */
    XW += x[41]*71.099790; /*C2H5COCH2 */
    XW += x[42]*56.064730; /*CH3CHCO */
    XW += x[43]*86.134850; /*NC3H7COCH3 */
    XW += x[44]*85.126880; /*NC3H7COCH2 */
    XW += x[45]*100.205570; /*NC7H16 */
    XW += x[46]*99.197600; /*C7H15-1 */
    XW += x[47]*131.196400; /*C7H15O2 */
    XW += x[48]*163.195200; /*C7H14OOHO2 */
    XW += x[49]*114.189030; /*C7H14O */
    XW += x[50]*146.187830; /*NC7KET */
    XW += x[51]*28.013400; /*N2 */
    *wtm = XW;

    return;
}

/*given c[species]: molar concentration */
/*returns mean molecular weight (gm/mole) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKMMWC(amrex::Real *  c,  amrex::Real *  wtm)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    amrex::Real W = 0;
    amrex::Real sumC = 0;
    W += c[0]*1.007970; /*H */
    W += c[1]*2.015940; /*H2 */
    W += c[2]*15.999400; /*O */
    W += c[3]*31.998800; /*O2 */
    W += c[4]*17.007370; /*OH */
    W += c[5]*18.015340; /*H2O */
    W += c[6]*28.010550; /*CO */
    W += c[7]*44.009950; /*CO2 */
    W += c[8]*15.035060; /*CH3 */
    W += c[9]*16.043030; /*CH4 */
    W += c[10]*33.006770; /*HO2 */
    W += c[11]*34.014740; /*H2O2 */
    W += c[12]*30.026490; /*CH2O */
    W += c[13]*30.070120; /*C2H6 */
    W += c[14]*28.054180; /*C2H4 */
    W += c[15]*29.062150; /*C2H5 */
    W += c[16]*25.030270; /*C2H */
    W += c[17]*26.038240; /*C2H2 */
    W += c[18]*32.042430; /*CH3OH */
    W += c[19]*42.037640; /*CH2CO */
    W += c[20]*41.029670; /*HCCO */
    W += c[21]*44.053580; /*CH3CHO */
    W += c[22]*40.065330; /*C3H4-A */
    W += c[23]*40.065330; /*C3H4-P */
    W += c[24]*42.081270; /*C3H6 */
    W += c[25]*54.092420; /*C4H6 */
    W += c[26]*55.100390; /*C4H7 */
    W += c[27]*56.108360; /*C4H8-1 */
    W += c[28]*57.116330; /*PC4H9 */
    W += c[29]*57.072700; /*CH3COCH2 */
    W += c[30]*58.080670; /*C2H5CHO */
    W += c[31]*69.127480; /*C5H9 */
    W += c[32]*70.135450; /*C5H10-1 */
    W += c[33]*71.143420; /*C5H11-1 */
    W += c[34]*47.033860; /*CH3O2 */
    W += c[35]*48.041830; /*CH3O2H */
    W += c[36]*55.056760; /*C2H3CO */
    W += c[37]*56.064730; /*C2H3CHO */
    W += c[38]*41.073300; /*C3H5-A */
    W += c[39]*39.057360; /*C3H3 */
    W += c[40]*72.107760; /*NC3H7CHO */
    W += c[41]*71.099790; /*C2H5COCH2 */
    W += c[42]*56.064730; /*CH3CHCO */
    W += c[43]*86.134850; /*NC3H7COCH3 */
    W += c[44]*85.126880; /*NC3H7COCH2 */
    W += c[45]*100.205570; /*NC7H16 */
    W += c[46]*99.197600; /*C7H15-1 */
    W += c[47]*131.196400; /*C7H15O2 */
    W += c[48]*163.195200; /*C7H14OOHO2 */
    W += c[49]*114.189030; /*C7H14O */
    W += c[50]*146.187830; /*NC7KET */
    W += c[51]*28.013400; /*N2 */

    for (id = 0; id < 52; ++id) {
        sumC += c[id];
    }
    /* CK provides no guard against divison by zero */
    *wtm = W/sumC;

    return;
}


/*get Cp/R as a function of T  */
/*for all species (Eq 19) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPOR(amrex::Real *  T, amrex::Real *  cpor)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpor, tc);
}


/*get H/RT as a function of T  */
/*for all species (Eq 20) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHORT(amrex::Real *  T, amrex::Real *  hort)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEnthalpy(hort, tc);
}


/*get S/R as a function of T  */
/*for all species (Eq 21) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSOR(amrex::Real *  T, amrex::Real *  sor)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sor, tc);
}

/*convert y[species] (mass fracs) to x[species] (mole fracs) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKYTX(amrex::Real *  y,  amrex::Real *  x)
{
    amrex::Real YOW = 0;
    amrex::Real tmp[52];
    amrex::Real imw[52];

    get_imw(imw);

    for (int i = 0; i < 52; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 52; i++)
    {
        YOW += tmp[i];
    }

    amrex::Real YOWINV = 1.0/YOW;

    for (int i = 0; i < 52; i++)
    {
        x[i] = y[i]*imw[i]*YOWINV;
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKYTCP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  c)
{
    amrex::Real YOW = 0;
    amrex::Real PWORT;
    amrex::Real imw[52];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    for (int i = 0; i < 52; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 52; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31446261815324e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 52; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKYTCR(amrex::Real *  rho, amrex::Real * /*T*/, amrex::Real * y,  amrex::Real * c)
{
    amrex::Real imw[52];

    get_imw(imw);

    for (int i = 0; i < 52; i++)
    {
        c[i] = (*rho)  * y[i] * imw[i];
    }
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTY(amrex::Real *  x,  amrex::Real *  y)
{
    amrex::Real XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*28.010550; /*CO */
    XW += x[7]*44.009950; /*CO2 */
    XW += x[8]*15.035060; /*CH3 */
    XW += x[9]*16.043030; /*CH4 */
    XW += x[10]*33.006770; /*HO2 */
    XW += x[11]*34.014740; /*H2O2 */
    XW += x[12]*30.026490; /*CH2O */
    XW += x[13]*30.070120; /*C2H6 */
    XW += x[14]*28.054180; /*C2H4 */
    XW += x[15]*29.062150; /*C2H5 */
    XW += x[16]*25.030270; /*C2H */
    XW += x[17]*26.038240; /*C2H2 */
    XW += x[18]*32.042430; /*CH3OH */
    XW += x[19]*42.037640; /*CH2CO */
    XW += x[20]*41.029670; /*HCCO */
    XW += x[21]*44.053580; /*CH3CHO */
    XW += x[22]*40.065330; /*C3H4-A */
    XW += x[23]*40.065330; /*C3H4-P */
    XW += x[24]*42.081270; /*C3H6 */
    XW += x[25]*54.092420; /*C4H6 */
    XW += x[26]*55.100390; /*C4H7 */
    XW += x[27]*56.108360; /*C4H8-1 */
    XW += x[28]*57.116330; /*PC4H9 */
    XW += x[29]*57.072700; /*CH3COCH2 */
    XW += x[30]*58.080670; /*C2H5CHO */
    XW += x[31]*69.127480; /*C5H9 */
    XW += x[32]*70.135450; /*C5H10-1 */
    XW += x[33]*71.143420; /*C5H11-1 */
    XW += x[34]*47.033860; /*CH3O2 */
    XW += x[35]*48.041830; /*CH3O2H */
    XW += x[36]*55.056760; /*C2H3CO */
    XW += x[37]*56.064730; /*C2H3CHO */
    XW += x[38]*41.073300; /*C3H5-A */
    XW += x[39]*39.057360; /*C3H3 */
    XW += x[40]*72.107760; /*NC3H7CHO */
    XW += x[41]*71.099790; /*C2H5COCH2 */
    XW += x[42]*56.064730; /*CH3CHCO */
    XW += x[43]*86.134850; /*NC3H7COCH3 */
    XW += x[44]*85.126880; /*NC3H7COCH2 */
    XW += x[45]*100.205570; /*NC7H16 */
    XW += x[46]*99.197600; /*C7H15-1 */
    XW += x[47]*131.196400; /*C7H15O2 */
    XW += x[48]*163.195200; /*C7H14OOHO2 */
    XW += x[49]*114.189030; /*C7H14O */
    XW += x[50]*146.187830; /*NC7KET */
    XW += x[51]*28.013400; /*N2 */
    /*Now compute conversion */
    amrex::Real XWinv = 1.0/XW;
    y[0] = x[0]*1.007970*XWinv; 
    y[1] = x[1]*2.015940*XWinv; 
    y[2] = x[2]*15.999400*XWinv; 
    y[3] = x[3]*31.998800*XWinv; 
    y[4] = x[4]*17.007370*XWinv; 
    y[5] = x[5]*18.015340*XWinv; 
    y[6] = x[6]*28.010550*XWinv; 
    y[7] = x[7]*44.009950*XWinv; 
    y[8] = x[8]*15.035060*XWinv; 
    y[9] = x[9]*16.043030*XWinv; 
    y[10] = x[10]*33.006770*XWinv; 
    y[11] = x[11]*34.014740*XWinv; 
    y[12] = x[12]*30.026490*XWinv; 
    y[13] = x[13]*30.070120*XWinv; 
    y[14] = x[14]*28.054180*XWinv; 
    y[15] = x[15]*29.062150*XWinv; 
    y[16] = x[16]*25.030270*XWinv; 
    y[17] = x[17]*26.038240*XWinv; 
    y[18] = x[18]*32.042430*XWinv; 
    y[19] = x[19]*42.037640*XWinv; 
    y[20] = x[20]*41.029670*XWinv; 
    y[21] = x[21]*44.053580*XWinv; 
    y[22] = x[22]*40.065330*XWinv; 
    y[23] = x[23]*40.065330*XWinv; 
    y[24] = x[24]*42.081270*XWinv; 
    y[25] = x[25]*54.092420*XWinv; 
    y[26] = x[26]*55.100390*XWinv; 
    y[27] = x[27]*56.108360*XWinv; 
    y[28] = x[28]*57.116330*XWinv; 
    y[29] = x[29]*57.072700*XWinv; 
    y[30] = x[30]*58.080670*XWinv; 
    y[31] = x[31]*69.127480*XWinv; 
    y[32] = x[32]*70.135450*XWinv; 
    y[33] = x[33]*71.143420*XWinv; 
    y[34] = x[34]*47.033860*XWinv; 
    y[35] = x[35]*48.041830*XWinv; 
    y[36] = x[36]*55.056760*XWinv; 
    y[37] = x[37]*56.064730*XWinv; 
    y[38] = x[38]*41.073300*XWinv; 
    y[39] = x[39]*39.057360*XWinv; 
    y[40] = x[40]*72.107760*XWinv; 
    y[41] = x[41]*71.099790*XWinv; 
    y[42] = x[42]*56.064730*XWinv; 
    y[43] = x[43]*86.134850*XWinv; 
    y[44] = x[44]*85.126880*XWinv; 
    y[45] = x[45]*100.205570*XWinv; 
    y[46] = x[46]*99.197600*XWinv; 
    y[47] = x[47]*131.196400*XWinv; 
    y[48] = x[48]*163.195200*XWinv; 
    y[49] = x[49]*114.189030*XWinv; 
    y[50] = x[50]*146.187830*XWinv; 
    y[51] = x[51]*28.013400*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTCP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  c)
{
    int id; /*loop counter */
    amrex::Real PORT = (*P)/(8.31446261815324e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 52; ++id) {
        c[id] = x[id]*PORT;
    }

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTCR(amrex::Real *  rho, amrex::Real * /*T*/, amrex::Real *  x, amrex::Real *  c)
{
    int id; /*loop counter */
    amrex::Real XW = 0; /*See Eq 4, 11 in CK Manual */
    amrex::Real ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*28.010550; /*CO */
    XW += x[7]*44.009950; /*CO2 */
    XW += x[8]*15.035060; /*CH3 */
    XW += x[9]*16.043030; /*CH4 */
    XW += x[10]*33.006770; /*HO2 */
    XW += x[11]*34.014740; /*H2O2 */
    XW += x[12]*30.026490; /*CH2O */
    XW += x[13]*30.070120; /*C2H6 */
    XW += x[14]*28.054180; /*C2H4 */
    XW += x[15]*29.062150; /*C2H5 */
    XW += x[16]*25.030270; /*C2H */
    XW += x[17]*26.038240; /*C2H2 */
    XW += x[18]*32.042430; /*CH3OH */
    XW += x[19]*42.037640; /*CH2CO */
    XW += x[20]*41.029670; /*HCCO */
    XW += x[21]*44.053580; /*CH3CHO */
    XW += x[22]*40.065330; /*C3H4-A */
    XW += x[23]*40.065330; /*C3H4-P */
    XW += x[24]*42.081270; /*C3H6 */
    XW += x[25]*54.092420; /*C4H6 */
    XW += x[26]*55.100390; /*C4H7 */
    XW += x[27]*56.108360; /*C4H8-1 */
    XW += x[28]*57.116330; /*PC4H9 */
    XW += x[29]*57.072700; /*CH3COCH2 */
    XW += x[30]*58.080670; /*C2H5CHO */
    XW += x[31]*69.127480; /*C5H9 */
    XW += x[32]*70.135450; /*C5H10-1 */
    XW += x[33]*71.143420; /*C5H11-1 */
    XW += x[34]*47.033860; /*CH3O2 */
    XW += x[35]*48.041830; /*CH3O2H */
    XW += x[36]*55.056760; /*C2H3CO */
    XW += x[37]*56.064730; /*C2H3CHO */
    XW += x[38]*41.073300; /*C3H5-A */
    XW += x[39]*39.057360; /*C3H3 */
    XW += x[40]*72.107760; /*NC3H7CHO */
    XW += x[41]*71.099790; /*C2H5COCH2 */
    XW += x[42]*56.064730; /*CH3CHCO */
    XW += x[43]*86.134850; /*NC3H7COCH3 */
    XW += x[44]*85.126880; /*NC3H7COCH2 */
    XW += x[45]*100.205570; /*NC7H16 */
    XW += x[46]*99.197600; /*C7H15-1 */
    XW += x[47]*131.196400; /*C7H15O2 */
    XW += x[48]*163.195200; /*C7H14OOHO2 */
    XW += x[49]*114.189030; /*C7H14O */
    XW += x[50]*146.187830; /*NC7KET */
    XW += x[51]*28.013400; /*N2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 52; ++id) {
        c[id] = x[id]*ROW;
    }

    return;
}


/*convert c[species] (molar conc) to x[species] (mole fracs) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCTX(amrex::Real *  c, amrex::Real *  x)
{
    int id; /*loop counter */
    amrex::Real sumC = 0; 

    /*compute sum of c  */
    for (id = 0; id < 52; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    amrex::Real sumCinv = 1.0/sumC;
    for (id = 0; id < 52; ++id) {
        x[id] = c[id]*sumCinv;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCTY(amrex::Real *  c, amrex::Real *  y)
{
    amrex::Real CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*1.007970; /*H */
    CW += c[1]*2.015940; /*H2 */
    CW += c[2]*15.999400; /*O */
    CW += c[3]*31.998800; /*O2 */
    CW += c[4]*17.007370; /*OH */
    CW += c[5]*18.015340; /*H2O */
    CW += c[6]*28.010550; /*CO */
    CW += c[7]*44.009950; /*CO2 */
    CW += c[8]*15.035060; /*CH3 */
    CW += c[9]*16.043030; /*CH4 */
    CW += c[10]*33.006770; /*HO2 */
    CW += c[11]*34.014740; /*H2O2 */
    CW += c[12]*30.026490; /*CH2O */
    CW += c[13]*30.070120; /*C2H6 */
    CW += c[14]*28.054180; /*C2H4 */
    CW += c[15]*29.062150; /*C2H5 */
    CW += c[16]*25.030270; /*C2H */
    CW += c[17]*26.038240; /*C2H2 */
    CW += c[18]*32.042430; /*CH3OH */
    CW += c[19]*42.037640; /*CH2CO */
    CW += c[20]*41.029670; /*HCCO */
    CW += c[21]*44.053580; /*CH3CHO */
    CW += c[22]*40.065330; /*C3H4-A */
    CW += c[23]*40.065330; /*C3H4-P */
    CW += c[24]*42.081270; /*C3H6 */
    CW += c[25]*54.092420; /*C4H6 */
    CW += c[26]*55.100390; /*C4H7 */
    CW += c[27]*56.108360; /*C4H8-1 */
    CW += c[28]*57.116330; /*PC4H9 */
    CW += c[29]*57.072700; /*CH3COCH2 */
    CW += c[30]*58.080670; /*C2H5CHO */
    CW += c[31]*69.127480; /*C5H9 */
    CW += c[32]*70.135450; /*C5H10-1 */
    CW += c[33]*71.143420; /*C5H11-1 */
    CW += c[34]*47.033860; /*CH3O2 */
    CW += c[35]*48.041830; /*CH3O2H */
    CW += c[36]*55.056760; /*C2H3CO */
    CW += c[37]*56.064730; /*C2H3CHO */
    CW += c[38]*41.073300; /*C3H5-A */
    CW += c[39]*39.057360; /*C3H3 */
    CW += c[40]*72.107760; /*NC3H7CHO */
    CW += c[41]*71.099790; /*C2H5COCH2 */
    CW += c[42]*56.064730; /*CH3CHCO */
    CW += c[43]*86.134850; /*NC3H7COCH3 */
    CW += c[44]*85.126880; /*NC3H7COCH2 */
    CW += c[45]*100.205570; /*NC7H16 */
    CW += c[46]*99.197600; /*C7H15-1 */
    CW += c[47]*131.196400; /*C7H15O2 */
    CW += c[48]*163.195200; /*C7H14OOHO2 */
    CW += c[49]*114.189030; /*C7H14O */
    CW += c[50]*146.187830; /*NC7KET */
    CW += c[51]*28.013400; /*N2 */
    /*Now compute conversion */
    amrex::Real CWinv = 1.0/CW;
    y[0] = c[0]*1.007970*CWinv; 
    y[1] = c[1]*2.015940*CWinv; 
    y[2] = c[2]*15.999400*CWinv; 
    y[3] = c[3]*31.998800*CWinv; 
    y[4] = c[4]*17.007370*CWinv; 
    y[5] = c[5]*18.015340*CWinv; 
    y[6] = c[6]*28.010550*CWinv; 
    y[7] = c[7]*44.009950*CWinv; 
    y[8] = c[8]*15.035060*CWinv; 
    y[9] = c[9]*16.043030*CWinv; 
    y[10] = c[10]*33.006770*CWinv; 
    y[11] = c[11]*34.014740*CWinv; 
    y[12] = c[12]*30.026490*CWinv; 
    y[13] = c[13]*30.070120*CWinv; 
    y[14] = c[14]*28.054180*CWinv; 
    y[15] = c[15]*29.062150*CWinv; 
    y[16] = c[16]*25.030270*CWinv; 
    y[17] = c[17]*26.038240*CWinv; 
    y[18] = c[18]*32.042430*CWinv; 
    y[19] = c[19]*42.037640*CWinv; 
    y[20] = c[20]*41.029670*CWinv; 
    y[21] = c[21]*44.053580*CWinv; 
    y[22] = c[22]*40.065330*CWinv; 
    y[23] = c[23]*40.065330*CWinv; 
    y[24] = c[24]*42.081270*CWinv; 
    y[25] = c[25]*54.092420*CWinv; 
    y[26] = c[26]*55.100390*CWinv; 
    y[27] = c[27]*56.108360*CWinv; 
    y[28] = c[28]*57.116330*CWinv; 
    y[29] = c[29]*57.072700*CWinv; 
    y[30] = c[30]*58.080670*CWinv; 
    y[31] = c[31]*69.127480*CWinv; 
    y[32] = c[32]*70.135450*CWinv; 
    y[33] = c[33]*71.143420*CWinv; 
    y[34] = c[34]*47.033860*CWinv; 
    y[35] = c[35]*48.041830*CWinv; 
    y[36] = c[36]*55.056760*CWinv; 
    y[37] = c[37]*56.064730*CWinv; 
    y[38] = c[38]*41.073300*CWinv; 
    y[39] = c[39]*39.057360*CWinv; 
    y[40] = c[40]*72.107760*CWinv; 
    y[41] = c[41]*71.099790*CWinv; 
    y[42] = c[42]*56.064730*CWinv; 
    y[43] = c[43]*86.134850*CWinv; 
    y[44] = c[44]*85.126880*CWinv; 
    y[45] = c[45]*100.205570*CWinv; 
    y[46] = c[46]*99.197600*CWinv; 
    y[47] = c[47]*131.196400*CWinv; 
    y[48] = c[48]*163.195200*CWinv; 
    y[49] = c[49]*114.189030*CWinv; 
    y[50] = c[50]*146.187830*CWinv; 
    y[51] = c[51]*28.013400*CWinv; 

    return;
}


/*get specific heat at constant volume as a function  */
/*of T for all species (molar units) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVML(amrex::Real *  T,  amrex::Real *  cvml)
{
    int id; /*loop counter */
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 52; ++id) {
        cvml[id] *= 8.31446261815324e+07;
    }
}


/*get specific heat at constant pressure as a  */
/*function of T for all species (molar units) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPML(amrex::Real *  T,  amrex::Real *  cpml)
{
    int id; /*loop counter */
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 52; ++id) {
        cpml[id] *= 8.31446261815324e+07;
    }
}


/*get internal energy as a function  */
/*of T for all species (molar units) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUML(amrex::Real *  T,  amrex::Real *  uml)
{
    int id; /*loop counter */
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 52; ++id) {
        uml[id] *= RT;
    }
}


/*get enthalpy as a function  */
/*of T for all species (molar units) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHML(amrex::Real *  T,  amrex::Real *  hml)
{
    int id; /*loop counter */
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 52; ++id) {
        hml[id] *= RT;
    }
}


/*Returns the standard-state entropies in molar units */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSML(amrex::Real *  T,  amrex::Real *  sml)
{
    int id; /*loop counter */
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 52; ++id) {
        sml[id] *= 8.31446261815324e+07;
    }
}


/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVMS(amrex::Real *  T,  amrex::Real *  cvms)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvms, tc);
    /*multiply by R/molecularweight */
    cvms[0] *= 8.248720317224957e+07; /*H */
    cvms[1] *= 4.124360158612479e+07; /*H2 */
    cvms[2] *= 5.196734013871295e+06; /*O */
    cvms[3] *= 2.598367006935648e+06; /*O2 */
    cvms[4] *= 4.888740950630956e+06; /*OH */
    cvms[5] *= 4.615212712140454e+06; /*H2O */
    cvms[6] *= 2.968332509769797e+06; /*CO */
    cvms[7] *= 1.889223372931176e+06; /*CO2 */
    cvms[8] *= 5.530049509714786e+06; /*CH3 */
    cvms[9] *= 5.182601178301878e+06; /*CH4 */
    cvms[10] *= 2.519017346487778e+06; /*HO2 */
    cvms[11] *= 2.444370475315478e+06; /*H2O2 */
    cvms[12] *= 2.769042474879095e+06; /*CH2O */
    cvms[13] *= 2.765024754857393e+06; /*C2H6 */
    cvms[14] *= 2.963716144315478e+06; /*C2H4 */
    cvms[15] *= 2.860924817383862e+06; /*C2H5 */
    cvms[16] *= 3.321763056552422e+06; /*C2H */
    cvms[17] *= 3.193173815954242e+06; /*C2H2 */
    cvms[18] *= 2.594828987112788e+06; /*CH3OH */
    cvms[19] *= 1.977861416138784e+06; /*CH2CO */
    cvms[20] *= 2.026451252996488e+06; /*HCCO */
    cvms[21] *= 1.887352314648035e+06; /*CH3CHO */
    cvms[22] *= 2.075226291198210e+06; /*C3H4-A */
    cvms[23] *= 2.075226291198210e+06; /*C3H4-P */
    cvms[24] *= 1.975810762876985e+06; /*C3H6 */
    cvms[25] *= 1.537084607816260e+06; /*C4H6 */
    cvms[26] *= 1.508966201174482e+06; /*C4H7 */
    cvms[27] *= 1.481858072157739e+06; /*C4H8-1 */
    cvms[28] *= 1.455706733635239e+06; /*PC4H9 */
    cvms[29] *= 1.456819568401922e+06; /*CH3COCH2 */
    cvms[30] *= 1.431536967144704e+06; /*C2H5CHO */
    cvms[31] *= 1.202772416722444e+06; /*C5H9 */
    cvms[32] *= 1.185486457726191e+06; /*C5H10-1 */
    cvms[33] *= 1.168690318535887e+06; /*C5H11-1 */
    cvms[34] *= 1.767761059405552e+06; /*CH3O2 */
    cvms[35] *= 1.730671504010825e+06; /*CH3O2H */
    cvms[36] *= 1.510161988855363e+06; /*C2H3CO */
    cvms[37] *= 1.483011265398627e+06; /*C2H3CHO */
    cvms[38] *= 2.024298660724422e+06; /*C3H5-A */
    cvms[39] *= 2.128782543969495e+06; /*C3H3 */
    cvms[40] *= 1.153060727188480e+06; /*NC3H7CHO */
    cvms[41] *= 1.169407479002855e+06; /*C2H5COCH2 */
    cvms[42] *= 1.483011265398627e+06; /*CH3CHCO */
    cvms[43] *= 9.652843904822775e+05; /*NC3H7COCH3 */
    cvms[44] *= 9.767141258029473e+05; /*NC3H7COCH2 */
    cvms[45] *= 8.297405641376262e+05; /*NC7H16 */
    cvms[46] *= 8.381717519529947e+05; /*C7H15-1 */
    cvms[47] *= 6.337416741734712e+05; /*C7H15O2 */
    cvms[48] *= 5.094796059046614e+05; /*C7H14OOHO2 */
    cvms[49] *= 7.281314692097165e+05; /*C7H14O */
    cvms[50] *= 5.687520375774945e+05; /*NC7KET */
    cvms[51] *= 2.968030520448514e+06; /*N2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPMS(amrex::Real *  T,  amrex::Real *  cpms)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 8.248720317224957e+07; /*H */
    cpms[1] *= 4.124360158612479e+07; /*H2 */
    cpms[2] *= 5.196734013871295e+06; /*O */
    cpms[3] *= 2.598367006935648e+06; /*O2 */
    cpms[4] *= 4.888740950630956e+06; /*OH */
    cpms[5] *= 4.615212712140454e+06; /*H2O */
    cpms[6] *= 2.968332509769797e+06; /*CO */
    cpms[7] *= 1.889223372931176e+06; /*CO2 */
    cpms[8] *= 5.530049509714786e+06; /*CH3 */
    cpms[9] *= 5.182601178301878e+06; /*CH4 */
    cpms[10] *= 2.519017346487778e+06; /*HO2 */
    cpms[11] *= 2.444370475315478e+06; /*H2O2 */
    cpms[12] *= 2.769042474879095e+06; /*CH2O */
    cpms[13] *= 2.765024754857393e+06; /*C2H6 */
    cpms[14] *= 2.963716144315478e+06; /*C2H4 */
    cpms[15] *= 2.860924817383862e+06; /*C2H5 */
    cpms[16] *= 3.321763056552422e+06; /*C2H */
    cpms[17] *= 3.193173815954242e+06; /*C2H2 */
    cpms[18] *= 2.594828987112788e+06; /*CH3OH */
    cpms[19] *= 1.977861416138784e+06; /*CH2CO */
    cpms[20] *= 2.026451252996488e+06; /*HCCO */
    cpms[21] *= 1.887352314648035e+06; /*CH3CHO */
    cpms[22] *= 2.075226291198210e+06; /*C3H4-A */
    cpms[23] *= 2.075226291198210e+06; /*C3H4-P */
    cpms[24] *= 1.975810762876985e+06; /*C3H6 */
    cpms[25] *= 1.537084607816260e+06; /*C4H6 */
    cpms[26] *= 1.508966201174482e+06; /*C4H7 */
    cpms[27] *= 1.481858072157739e+06; /*C4H8-1 */
    cpms[28] *= 1.455706733635239e+06; /*PC4H9 */
    cpms[29] *= 1.456819568401922e+06; /*CH3COCH2 */
    cpms[30] *= 1.431536967144704e+06; /*C2H5CHO */
    cpms[31] *= 1.202772416722444e+06; /*C5H9 */
    cpms[32] *= 1.185486457726191e+06; /*C5H10-1 */
    cpms[33] *= 1.168690318535887e+06; /*C5H11-1 */
    cpms[34] *= 1.767761059405552e+06; /*CH3O2 */
    cpms[35] *= 1.730671504010825e+06; /*CH3O2H */
    cpms[36] *= 1.510161988855363e+06; /*C2H3CO */
    cpms[37] *= 1.483011265398627e+06; /*C2H3CHO */
    cpms[38] *= 2.024298660724422e+06; /*C3H5-A */
    cpms[39] *= 2.128782543969495e+06; /*C3H3 */
    cpms[40] *= 1.153060727188480e+06; /*NC3H7CHO */
    cpms[41] *= 1.169407479002855e+06; /*C2H5COCH2 */
    cpms[42] *= 1.483011265398627e+06; /*CH3CHCO */
    cpms[43] *= 9.652843904822775e+05; /*NC3H7COCH3 */
    cpms[44] *= 9.767141258029473e+05; /*NC3H7COCH2 */
    cpms[45] *= 8.297405641376262e+05; /*NC7H16 */
    cpms[46] *= 8.381717519529947e+05; /*C7H15-1 */
    cpms[47] *= 6.337416741734712e+05; /*C7H15O2 */
    cpms[48] *= 5.094796059046614e+05; /*C7H14OOHO2 */
    cpms[49] *= 7.281314692097165e+05; /*C7H14O */
    cpms[50] *= 5.687520375774945e+05; /*NC7KET */
    cpms[51] *= 2.968030520448514e+06; /*N2 */
}


/*Returns internal energy in mass units (Eq 30.) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUMS(amrex::Real *  T,  amrex::Real *  ums)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[52];

    get_imw(imw);

    speciesInternalEnergy(ums, tc);

    for (int i = 0; i < 52; i++)
    {
        ums[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHMS(amrex::Real *  T,  amrex::Real *  hms)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[52];

    get_imw(imw);

    speciesEnthalpy(hms, tc);

    for (int i = 0; i < 52; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


/*Returns the entropies in mass units (Eq 28.) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSMS(amrex::Real *  T,  amrex::Real *  sms)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sms, tc);
    /*multiply by R/molecularweight */
    sms[0] *= 8.248720317224957e+07; /*H */
    sms[1] *= 4.124360158612479e+07; /*H2 */
    sms[2] *= 5.196734013871295e+06; /*O */
    sms[3] *= 2.598367006935648e+06; /*O2 */
    sms[4] *= 4.888740950630956e+06; /*OH */
    sms[5] *= 4.615212712140454e+06; /*H2O */
    sms[6] *= 2.968332509769797e+06; /*CO */
    sms[7] *= 1.889223372931176e+06; /*CO2 */
    sms[8] *= 5.530049509714786e+06; /*CH3 */
    sms[9] *= 5.182601178301878e+06; /*CH4 */
    sms[10] *= 2.519017346487778e+06; /*HO2 */
    sms[11] *= 2.444370475315478e+06; /*H2O2 */
    sms[12] *= 2.769042474879095e+06; /*CH2O */
    sms[13] *= 2.765024754857393e+06; /*C2H6 */
    sms[14] *= 2.963716144315478e+06; /*C2H4 */
    sms[15] *= 2.860924817383862e+06; /*C2H5 */
    sms[16] *= 3.321763056552422e+06; /*C2H */
    sms[17] *= 3.193173815954242e+06; /*C2H2 */
    sms[18] *= 2.594828987112788e+06; /*CH3OH */
    sms[19] *= 1.977861416138784e+06; /*CH2CO */
    sms[20] *= 2.026451252996488e+06; /*HCCO */
    sms[21] *= 1.887352314648035e+06; /*CH3CHO */
    sms[22] *= 2.075226291198210e+06; /*C3H4-A */
    sms[23] *= 2.075226291198210e+06; /*C3H4-P */
    sms[24] *= 1.975810762876985e+06; /*C3H6 */
    sms[25] *= 1.537084607816260e+06; /*C4H6 */
    sms[26] *= 1.508966201174482e+06; /*C4H7 */
    sms[27] *= 1.481858072157739e+06; /*C4H8-1 */
    sms[28] *= 1.455706733635239e+06; /*PC4H9 */
    sms[29] *= 1.456819568401922e+06; /*CH3COCH2 */
    sms[30] *= 1.431536967144704e+06; /*C2H5CHO */
    sms[31] *= 1.202772416722444e+06; /*C5H9 */
    sms[32] *= 1.185486457726191e+06; /*C5H10-1 */
    sms[33] *= 1.168690318535887e+06; /*C5H11-1 */
    sms[34] *= 1.767761059405552e+06; /*CH3O2 */
    sms[35] *= 1.730671504010825e+06; /*CH3O2H */
    sms[36] *= 1.510161988855363e+06; /*C2H3CO */
    sms[37] *= 1.483011265398627e+06; /*C2H3CHO */
    sms[38] *= 2.024298660724422e+06; /*C3H5-A */
    sms[39] *= 2.128782543969495e+06; /*C3H3 */
    sms[40] *= 1.153060727188480e+06; /*NC3H7CHO */
    sms[41] *= 1.169407479002855e+06; /*C2H5COCH2 */
    sms[42] *= 1.483011265398627e+06; /*CH3CHCO */
    sms[43] *= 9.652843904822775e+05; /*NC3H7COCH3 */
    sms[44] *= 9.767141258029473e+05; /*NC3H7COCH2 */
    sms[45] *= 8.297405641376262e+05; /*NC7H16 */
    sms[46] *= 8.381717519529947e+05; /*C7H15-1 */
    sms[47] *= 6.337416741734712e+05; /*C7H15O2 */
    sms[48] *= 5.094796059046614e+05; /*C7H14OOHO2 */
    sms[49] *= 7.281314692097165e+05; /*C7H14O */
    sms[50] *= 5.687520375774945e+05; /*NC7KET */
    sms[51] *= 2.968030520448514e+06; /*N2 */
}


/* GPU version of productionRate: no more use of thermo namespace vectors  */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void comp_qfqr(amrex::Real *  qf, amrex::Real * qr, amrex::Real * sc, amrex::Real * sc_qss, amrex::Real * tc, amrex::Real invT)
{

    return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void  productionRate(amrex::Real * wdot, amrex::Real * sc, amrex::Real T)
{
    amrex::Real tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    amrex::Real invT = 1.0 / tc[1];


    for (int i = 0; i < 52; ++i) {
        wdot[i] = 0.0;
    }

    return;
}



/*compute the production rate for each species */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWC(amrex::Real *  T, amrex::Real *  C,  amrex::Real *  wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 52; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 52; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[52]; /*temporary storage */
    amrex::Real YOW = 0; 
    amrex::Real PWORT; 
    amrex::Real imw[52];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*O2 */
    YOW += y[4]*imw[4]; /*OH */
    YOW += y[5]*imw[5]; /*H2O */
    YOW += y[6]*imw[6]; /*CO */
    YOW += y[7]*imw[7]; /*CO2 */
    YOW += y[8]*imw[8]; /*CH3 */
    YOW += y[9]*imw[9]; /*CH4 */
    YOW += y[10]*imw[10]; /*HO2 */
    YOW += y[11]*imw[11]; /*H2O2 */
    YOW += y[12]*imw[12]; /*CH2O */
    YOW += y[13]*imw[13]; /*C2H6 */
    YOW += y[14]*imw[14]; /*C2H4 */
    YOW += y[15]*imw[15]; /*C2H5 */
    YOW += y[16]*imw[16]; /*C2H */
    YOW += y[17]*imw[17]; /*C2H2 */
    YOW += y[18]*imw[18]; /*CH3OH */
    YOW += y[19]*imw[19]; /*CH2CO */
    YOW += y[20]*imw[20]; /*HCCO */
    YOW += y[21]*imw[21]; /*CH3CHO */
    YOW += y[22]*imw[22]; /*C3H4-A */
    YOW += y[23]*imw[23]; /*C3H4-P */
    YOW += y[24]*imw[24]; /*C3H6 */
    YOW += y[25]*imw[25]; /*C4H6 */
    YOW += y[26]*imw[26]; /*C4H7 */
    YOW += y[27]*imw[27]; /*C4H8-1 */
    YOW += y[28]*imw[28]; /*PC4H9 */
    YOW += y[29]*imw[29]; /*CH3COCH2 */
    YOW += y[30]*imw[30]; /*C2H5CHO */
    YOW += y[31]*imw[31]; /*C5H9 */
    YOW += y[32]*imw[32]; /*C5H10-1 */
    YOW += y[33]*imw[33]; /*C5H11-1 */
    YOW += y[34]*imw[34]; /*CH3O2 */
    YOW += y[35]*imw[35]; /*CH3O2H */
    YOW += y[36]*imw[36]; /*C2H3CO */
    YOW += y[37]*imw[37]; /*C2H3CHO */
    YOW += y[38]*imw[38]; /*C3H5-A */
    YOW += y[39]*imw[39]; /*C3H3 */
    YOW += y[40]*imw[40]; /*NC3H7CHO */
    YOW += y[41]*imw[41]; /*C2H5COCH2 */
    YOW += y[42]*imw[42]; /*CH3CHCO */
    YOW += y[43]*imw[43]; /*NC3H7COCH3 */
    YOW += y[44]*imw[44]; /*NC3H7COCH2 */
    YOW += y[45]*imw[45]; /*NC7H16 */
    YOW += y[46]*imw[46]; /*C7H15-1 */
    YOW += y[47]*imw[47]; /*C7H15O2 */
    YOW += y[48]*imw[48]; /*C7H14OOHO2 */
    YOW += y[49]*imw[49]; /*C7H14O */
    YOW += y[50]*imw[50]; /*NC7KET */
    YOW += y[51]*imw[51]; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31446261815324e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]*imw[0]; 
    c[1] = PWORT * y[1]*imw[1]; 
    c[2] = PWORT * y[2]*imw[2]; 
    c[3] = PWORT * y[3]*imw[3]; 
    c[4] = PWORT * y[4]*imw[4]; 
    c[5] = PWORT * y[5]*imw[5]; 
    c[6] = PWORT * y[6]*imw[6]; 
    c[7] = PWORT * y[7]*imw[7]; 
    c[8] = PWORT * y[8]*imw[8]; 
    c[9] = PWORT * y[9]*imw[9]; 
    c[10] = PWORT * y[10]*imw[10]; 
    c[11] = PWORT * y[11]*imw[11]; 
    c[12] = PWORT * y[12]*imw[12]; 
    c[13] = PWORT * y[13]*imw[13]; 
    c[14] = PWORT * y[14]*imw[14]; 
    c[15] = PWORT * y[15]*imw[15]; 
    c[16] = PWORT * y[16]*imw[16]; 
    c[17] = PWORT * y[17]*imw[17]; 
    c[18] = PWORT * y[18]*imw[18]; 
    c[19] = PWORT * y[19]*imw[19]; 
    c[20] = PWORT * y[20]*imw[20]; 
    c[21] = PWORT * y[21]*imw[21]; 
    c[22] = PWORT * y[22]*imw[22]; 
    c[23] = PWORT * y[23]*imw[23]; 
    c[24] = PWORT * y[24]*imw[24]; 
    c[25] = PWORT * y[25]*imw[25]; 
    c[26] = PWORT * y[26]*imw[26]; 
    c[27] = PWORT * y[27]*imw[27]; 
    c[28] = PWORT * y[28]*imw[28]; 
    c[29] = PWORT * y[29]*imw[29]; 
    c[30] = PWORT * y[30]*imw[30]; 
    c[31] = PWORT * y[31]*imw[31]; 
    c[32] = PWORT * y[32]*imw[32]; 
    c[33] = PWORT * y[33]*imw[33]; 
    c[34] = PWORT * y[34]*imw[34]; 
    c[35] = PWORT * y[35]*imw[35]; 
    c[36] = PWORT * y[36]*imw[36]; 
    c[37] = PWORT * y[37]*imw[37]; 
    c[38] = PWORT * y[38]*imw[38]; 
    c[39] = PWORT * y[39]*imw[39]; 
    c[40] = PWORT * y[40]*imw[40]; 
    c[41] = PWORT * y[41]*imw[41]; 
    c[42] = PWORT * y[42]*imw[42]; 
    c[43] = PWORT * y[43]*imw[43]; 
    c[44] = PWORT * y[44]*imw[44]; 
    c[45] = PWORT * y[45]*imw[45]; 
    c[46] = PWORT * y[46]*imw[46]; 
    c[47] = PWORT * y[47]*imw[47]; 
    c[48] = PWORT * y[48]*imw[48]; 
    c[49] = PWORT * y[49]*imw[49]; 
    c[50] = PWORT * y[50]*imw[50]; 
    c[51] = PWORT * y[51]*imw[51]; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 52; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[52]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 52; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 52; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[52]; /*temporary storage */
    amrex::Real imw[52];

    get_imw(imw);

    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]*imw[0]; 
    c[1] = 1e6 * (*rho) * y[1]*imw[1]; 
    c[2] = 1e6 * (*rho) * y[2]*imw[2]; 
    c[3] = 1e6 * (*rho) * y[3]*imw[3]; 
    c[4] = 1e6 * (*rho) * y[4]*imw[4]; 
    c[5] = 1e6 * (*rho) * y[5]*imw[5]; 
    c[6] = 1e6 * (*rho) * y[6]*imw[6]; 
    c[7] = 1e6 * (*rho) * y[7]*imw[7]; 
    c[8] = 1e6 * (*rho) * y[8]*imw[8]; 
    c[9] = 1e6 * (*rho) * y[9]*imw[9]; 
    c[10] = 1e6 * (*rho) * y[10]*imw[10]; 
    c[11] = 1e6 * (*rho) * y[11]*imw[11]; 
    c[12] = 1e6 * (*rho) * y[12]*imw[12]; 
    c[13] = 1e6 * (*rho) * y[13]*imw[13]; 
    c[14] = 1e6 * (*rho) * y[14]*imw[14]; 
    c[15] = 1e6 * (*rho) * y[15]*imw[15]; 
    c[16] = 1e6 * (*rho) * y[16]*imw[16]; 
    c[17] = 1e6 * (*rho) * y[17]*imw[17]; 
    c[18] = 1e6 * (*rho) * y[18]*imw[18]; 
    c[19] = 1e6 * (*rho) * y[19]*imw[19]; 
    c[20] = 1e6 * (*rho) * y[20]*imw[20]; 
    c[21] = 1e6 * (*rho) * y[21]*imw[21]; 
    c[22] = 1e6 * (*rho) * y[22]*imw[22]; 
    c[23] = 1e6 * (*rho) * y[23]*imw[23]; 
    c[24] = 1e6 * (*rho) * y[24]*imw[24]; 
    c[25] = 1e6 * (*rho) * y[25]*imw[25]; 
    c[26] = 1e6 * (*rho) * y[26]*imw[26]; 
    c[27] = 1e6 * (*rho) * y[27]*imw[27]; 
    c[28] = 1e6 * (*rho) * y[28]*imw[28]; 
    c[29] = 1e6 * (*rho) * y[29]*imw[29]; 
    c[30] = 1e6 * (*rho) * y[30]*imw[30]; 
    c[31] = 1e6 * (*rho) * y[31]*imw[31]; 
    c[32] = 1e6 * (*rho) * y[32]*imw[32]; 
    c[33] = 1e6 * (*rho) * y[33]*imw[33]; 
    c[34] = 1e6 * (*rho) * y[34]*imw[34]; 
    c[35] = 1e6 * (*rho) * y[35]*imw[35]; 
    c[36] = 1e6 * (*rho) * y[36]*imw[36]; 
    c[37] = 1e6 * (*rho) * y[37]*imw[37]; 
    c[38] = 1e6 * (*rho) * y[38]*imw[38]; 
    c[39] = 1e6 * (*rho) * y[39]*imw[39]; 
    c[40] = 1e6 * (*rho) * y[40]*imw[40]; 
    c[41] = 1e6 * (*rho) * y[41]*imw[41]; 
    c[42] = 1e6 * (*rho) * y[42]*imw[42]; 
    c[43] = 1e6 * (*rho) * y[43]*imw[43]; 
    c[44] = 1e6 * (*rho) * y[44]*imw[44]; 
    c[45] = 1e6 * (*rho) * y[45]*imw[45]; 
    c[46] = 1e6 * (*rho) * y[46]*imw[46]; 
    c[47] = 1e6 * (*rho) * y[47]*imw[47]; 
    c[48] = 1e6 * (*rho) * y[48]*imw[48]; 
    c[49] = 1e6 * (*rho) * y[49]*imw[49]; 
    c[50] = 1e6 * (*rho) * y[50]*imw[50]; 
    c[51] = 1e6 * (*rho) * y[51]*imw[51]; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 52; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[52]; /*temporary storage */
    amrex::Real XW = 0; /*See Eq 4, 11 in CK Manual */
    amrex::Real ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*28.010550; /*CO */
    XW += x[7]*44.009950; /*CO2 */
    XW += x[8]*15.035060; /*CH3 */
    XW += x[9]*16.043030; /*CH4 */
    XW += x[10]*33.006770; /*HO2 */
    XW += x[11]*34.014740; /*H2O2 */
    XW += x[12]*30.026490; /*CH2O */
    XW += x[13]*30.070120; /*C2H6 */
    XW += x[14]*28.054180; /*C2H4 */
    XW += x[15]*29.062150; /*C2H5 */
    XW += x[16]*25.030270; /*C2H */
    XW += x[17]*26.038240; /*C2H2 */
    XW += x[18]*32.042430; /*CH3OH */
    XW += x[19]*42.037640; /*CH2CO */
    XW += x[20]*41.029670; /*HCCO */
    XW += x[21]*44.053580; /*CH3CHO */
    XW += x[22]*40.065330; /*C3H4-A */
    XW += x[23]*40.065330; /*C3H4-P */
    XW += x[24]*42.081270; /*C3H6 */
    XW += x[25]*54.092420; /*C4H6 */
    XW += x[26]*55.100390; /*C4H7 */
    XW += x[27]*56.108360; /*C4H8-1 */
    XW += x[28]*57.116330; /*PC4H9 */
    XW += x[29]*57.072700; /*CH3COCH2 */
    XW += x[30]*58.080670; /*C2H5CHO */
    XW += x[31]*69.127480; /*C5H9 */
    XW += x[32]*70.135450; /*C5H10-1 */
    XW += x[33]*71.143420; /*C5H11-1 */
    XW += x[34]*47.033860; /*CH3O2 */
    XW += x[35]*48.041830; /*CH3O2H */
    XW += x[36]*55.056760; /*C2H3CO */
    XW += x[37]*56.064730; /*C2H3CHO */
    XW += x[38]*41.073300; /*C3H5-A */
    XW += x[39]*39.057360; /*C3H3 */
    XW += x[40]*72.107760; /*NC3H7CHO */
    XW += x[41]*71.099790; /*C2H5COCH2 */
    XW += x[42]*56.064730; /*CH3CHCO */
    XW += x[43]*86.134850; /*NC3H7COCH3 */
    XW += x[44]*85.126880; /*NC3H7COCH2 */
    XW += x[45]*100.205570; /*NC7H16 */
    XW += x[46]*99.197600; /*C7H15-1 */
    XW += x[47]*131.196400; /*C7H15O2 */
    XW += x[48]*163.195200; /*C7H14OOHO2 */
    XW += x[49]*114.189030; /*C7H14O */
    XW += x[50]*146.187830; /*NC7KET */
    XW += x[51]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 52; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 52; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void dcvpRdT(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 1: H2 */
        species[1] =
            +8.24944200e-04
            -1.62860300e-06 * tc[1]
            -2.84263020e-10 * tc[2]
            +1.65394880e-12 * tc[3];
        /*species 2: O */
        species[2] =
            -1.63816600e-03
            +4.84206400e-06 * tc[1]
            -4.80852900e-09 * tc[2]
            +1.55627840e-12 * tc[3];
        /*species 3: O2 */
        species[3] =
            +1.12748600e-03
            -1.15123000e-06 * tc[1]
            +3.94163100e-09 * tc[2]
            -3.50742160e-12 * tc[3];
        /*species 5: H2O */
        species[5] =
            +3.47498200e-03
            -1.27093920e-05 * tc[1]
            +2.09057430e-08 * tc[2]
            -1.00263520e-11 * tc[3];
        /*species 6: CO */
        species[6] =
            +1.51194100e-03
            -7.76351000e-06 * tc[1]
            +1.67458320e-08 * tc[2]
            -9.89980400e-12 * tc[3];
        /*species 7: CO2 */
        species[7] =
            +9.92207200e-03
            -2.08182200e-05 * tc[1]
            +2.06000610e-08 * tc[2]
            -8.46912000e-12 * tc[3];
        /*species 8: CH3 */
        species[8] =
            +1.11241000e-02
            -3.36044000e-05 * tc[1]
            +4.86548700e-08 * tc[2]
            -2.34598120e-11 * tc[3];
        /*species 9: CH4 */
        species[9] =
            +1.74766800e-02
            -5.56681800e-05 * tc[1]
            +9.14912400e-08 * tc[2]
            -4.89572400e-11 * tc[3];
        /*species 11: H2O2 */
        species[11] =
            +6.56922600e-03
            -2.97002600e-07 * tc[1]
            -1.38774180e-08 * tc[2]
            +9.88606000e-12 * tc[3];
        /*species 12: CH2O */
        species[12] =
            +1.26314400e-02
            -3.77633600e-05 * tc[1]
            +6.15009300e-08 * tc[2]
            -3.36529480e-11 * tc[3];
        /*species 14: C2H4 */
        species[14] =
            +2.79616300e-02
            -6.77735400e-05 * tc[1]
            +8.35545600e-08 * tc[2]
            -3.89515160e-11 * tc[3];
        /*species 15: C2H5 */
        species[15] =
            +8.71913300e-03
            +8.83967800e-06 * tc[1]
            +2.80161090e-09 * tc[2]
            -1.57110920e-11 * tc[3];
        /*species 16: C2H */
        species[16] =
            +8.04844600e-03
            -1.84886200e-05 * tc[1]
            +1.95757770e-08 * tc[2]
            -7.75832000e-12 * tc[3];
        /*species 17: C2H2 */
        species[17] =
            +1.51904500e-02
            -3.23263800e-05 * tc[1]
            +2.72369760e-08 * tc[2]
            -7.65098400e-12 * tc[3];
        /*species 18: CH3OH */
        species[18] =
            +7.34150800e-03
            +1.43401020e-05 * tc[1]
            -2.63795820e-08 * tc[2]
            +9.56228000e-12 * tc[3];
        /*species 19: CH2CO */
        species[19] =
            +1.21187100e-02
            -4.69009200e-06 * tc[1]
            -1.94000550e-08 * tc[2]
            +1.56225960e-11 * tc[3];
        /*species 20: HCCO */
        species[20] =
            +4.45347800e-03
            +4.53656600e-07 * tc[1]
            -4.44628500e-09 * tc[2]
            +9.00296800e-13 * tc[3];
        /*species 21: CH3CHO */
        species[21] =
            +1.33699100e-02
            +9.34390600e-06 * tc[1]
            -3.38442000e-08 * tc[2]
            +1.70542640e-11 * tc[3];
        /*species 39: C3H3 */
        species[39] =
            +1.10802800e-02
            +5.58664600e-07 * tc[1]
            -1.64376360e-08 * tc[2]
            +7.79851600e-12 * tc[3];
        /*species 51: N2 */
        species[51] =
            +1.40824000e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77942000e-12 * tc[3];
    } else {
        /*species 0: H */
        species[0] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 1: H2 */
        species[1] =
            +7.00064400e-04
            -1.12676580e-07 * tc[1]
            -2.76947340e-11 * tc[2]
            +6.33100800e-15 * tc[3];
        /*species 2: O */
        species[2] =
            -2.75506200e-05
            -6.20560600e-09 * tc[1]
            +1.36532010e-11 * tc[2]
            -1.74722080e-15 * tc[3];
        /*species 3: O2 */
        species[3] =
            +6.13519700e-04
            -2.51768400e-07 * tc[1]
            +5.32584300e-11 * tc[2]
            -4.54574000e-15 * tc[3];
        /*species 5: H2O */
        species[5] =
            +3.05629300e-03
            -1.74605200e-06 * tc[1]
            +3.60298800e-10 * tc[2]
            -2.55664720e-14 * tc[3];
        /*species 6: CO */
        species[6] =
            +1.44268900e-03
            -1.12616560e-06 * tc[1]
            +3.05574300e-10 * tc[2]
            -2.76438080e-14 * tc[3];
        /*species 7: CO2 */
        species[7] =
            +3.14016900e-03
            -2.55682200e-06 * tc[1]
            +7.18199100e-10 * tc[2]
            -6.67613200e-14 * tc[3];
        /*species 8: CH3 */
        species[8] =
            +6.13797400e-03
            -4.46069000e-06 * tc[1]
            +1.13554830e-09 * tc[2]
            -9.80863600e-14 * tc[3];
        /*species 9: CH4 */
        species[9] =
            +1.02372400e-02
            -7.75025800e-06 * tc[1]
            +2.03567550e-09 * tc[2]
            -1.80136920e-13 * tc[3];
        /*species 11: H2O2 */
        species[11] =
            +4.33613600e-03
            -2.94937800e-06 * tc[1]
            +7.04671200e-10 * tc[2]
            -5.72661600e-14 * tc[3];
        /*species 12: CH2O */
        species[12] =
            +6.68132100e-03
            -5.25791000e-06 * tc[1]
            +1.42114590e-09 * tc[2]
            -1.28500680e-13 * tc[3];
        /*species 14: C2H4 */
        species[14] =
            +1.14851800e-02
            -8.83677000e-06 * tc[1]
            +2.35338030e-09 * tc[2]
            -2.10673920e-13 * tc[3];
        /*species 15: C2H5 */
        species[15] =
            +6.48407700e-03
            -1.28561300e-06 * tc[1]
            -7.04363700e-10 * tc[2]
            +1.55235080e-13 * tc[3];
        /*species 16: C2H */
        species[16] =
            +3.14312300e-03
            -2.53448600e-06 * tc[1]
            +8.77308900e-10 * tc[2]
            -1.08652800e-13 * tc[3];
        /*species 17: C2H2 */
        species[17] =
            +5.37603900e-03
            -3.82563400e-06 * tc[1]
            +9.85913700e-10 * tc[2]
            -8.62684000e-14 * tc[3];
        /*species 18: CH3OH */
        species[18] =
            +9.37659300e-03
            -6.10050800e-06 * tc[1]
            +1.30763790e-09 * tc[2]
            -8.89889200e-14 * tc[3];
        /*species 19: CH2CO */
        species[19] =
            +5.80484000e-03
            -3.84190800e-06 * tc[1]
            +8.38345500e-10 * tc[2]
            -5.83547200e-14 * tc[3];
        /*species 20: HCCO */
        species[20] =
            +2.00040000e-03
            -4.05521400e-07 * tc[1]
            -3.12339600e-10 * tc[2]
            +7.86066000e-14 * tc[3];
        /*species 21: CH3CHO */
        species[21] =
            +1.07942400e-02
            -7.29106000e-06 * tc[1]
            +1.62387360e-09 * tc[2]
            -1.15873760e-13 * tc[3];
        /*species 39: C3H3 */
        species[39] =
            +4.35719500e-03
            -8.21813400e-07 * tc[1]
            -7.10616900e-10 * tc[2]
            +1.75060800e-13 * tc[3];
        /*species 51: N2 */
        species[51] =
            +1.48797700e-03
            -1.13695220e-06 * tc[1]
            +3.02911200e-10 * tc[2]
            -2.70134040e-14 * tc[3];
    }

    /*species with midpoint at T=1357 kelvin */
    if (T < 1357) {
        /*species 4: OH */
        species[4] =
            +2.02235804e-04
            -2.27092824e-07 * tc[1]
            +7.27335447e-10 * tc[2]
            -2.97460412e-13 * tc[3];
    } else {
        /*species 4: OH */
        species[4] =
            +1.31992406e-03
            -7.19449340e-07 * tc[1]
            +1.27689240e-10 * tc[2]
            -7.28192064e-15 * tc[3];
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 10: HO2 */
        species[10] =
            +3.66767950e-03
            -1.86477024e-06 * tc[1]
            -9.77558757e-10 * tc[2]
            +6.04559648e-13 * tc[3];
        /*species 35: CH3O2H */
        species[35] =
            +1.90129767e-02
            -2.26772574e-05 * tc[1]
            +1.02091996e-08 * tc[2]
            -1.64732089e-12 * tc[3];
        /*species 46: C7H15-1 */
        species[46] =
            +8.08826467e-02
            -1.00106551e-04 * tc[1]
            +4.69647924e-08 * tc[2]
            -7.86464908e-12 * tc[3];
        /*species 47: C7H15O2 */
        species[47] =
            +9.42381007e-02
            -1.33391021e-04 * tc[1]
            +7.64386749e-08 * tc[2]
            -1.65284442e-11 * tc[3];
        /*species 48: C7H14OOHO2 */
        species[48] =
            +1.10720834e-01
            -1.74235615e-04 * tc[1]
            +1.09098825e-07 * tc[2]
            -2.51340276e-11 * tc[3];
    } else {
        /*species 10: HO2 */
        species[10] =
            +2.38452835e-03
            -1.61269598e-06 * tc[1]
            +3.72575169e-10 * tc[2]
            -2.86560043e-14 * tc[3];
        /*species 35: CH3O2H */
        species[35] =
            +8.06817909e-03
            -5.54189842e-06 * tc[1]
            +1.29399673e-09 * tc[2]
            -1.00276858e-13 * tc[3];
        /*species 46: C7H15-1 */
        species[46] =
            +3.26280243e-02
            -2.22276488e-05 * tc[1]
            +5.16201444e-09 * tc[2]
            -3.98546800e-13 * tc[3];
        /*species 47: C7H15O2 */
        species[47] =
            +3.43831688e-02
            -2.37426494e-05 * tc[1]
            +5.56350723e-09 * tc[2]
            -4.32209596e-13 * tc[3];
        /*species 48: C7H14OOHO2 */
        species[48] =
            +3.32767285e-02
            -2.31666570e-05 * tc[1]
            +5.45875914e-09 * tc[2]
            -4.25737676e-13 * tc[3];
    }

    /*species with midpoint at T=1384 kelvin */
    if (T < 1384) {
        /*species 13: C2H6 */
        species[13] =
            +2.40764754e-02
            -2.23786944e-05 * tc[1]
            +6.25022703e-09 * tc[2]
            -2.11947446e-13 * tc[3];
    } else {
        /*species 13: C2H6 */
        species[13] =
            +1.29236361e-02
            -8.85054392e-06 * tc[1]
            +2.06217518e-09 * tc[2]
            -1.59560693e-13 * tc[3];
    }

    /*species with midpoint at T=1400 kelvin */
    if (T < 1400) {
        /*species 22: C3H4-A */
        species[22] =
            +1.63343700e-02
            -3.52990000e-06 * tc[1]
            -1.39420950e-08 * tc[2]
            +6.91652400e-12 * tc[3];
        /*species 23: C3H4-P */
        species[23] =
            +1.49896100e-02
            -2.79700000e-06 * tc[1]
            -1.19088570e-08 * tc[2]
            +5.55286800e-12 * tc[3];
        /*species 42: CH3CHCO */
        species[42] =
            +3.22203013e-02
            -5.40500066e-05 * tc[1]
            +3.61497492e-08 * tc[2]
            -8.73463724e-12 * tc[3];
    } else {
        /*species 22: C3H4-A */
        species[22] =
            +5.30213800e-03
            -7.40223600e-07 * tc[1]
            -9.07915800e-10 * tc[2]
            +2.03583240e-13 * tc[3];
        /*species 23: C3H4-P */
        species[23] =
            +5.21915100e-03
            -7.50628000e-07 * tc[1]
            -8.97657300e-10 * tc[2]
            +2.04315120e-13 * tc[3];
        /*species 42: CH3CHCO */
        species[42] =
            +9.56966300e-03
            -6.52443288e-06 * tc[1]
            +1.51569512e-09 * tc[2]
            -1.17037303e-13 * tc[3];
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 24: C3H6 */
        species[24] =
            +2.89107662e-02
            -3.09773616e-05 * tc[1]
            +1.16644263e-08 * tc[2]
            -1.35156141e-12 * tc[3];
    } else {
        /*species 24: C3H6 */
        species[24] =
            +1.37023634e-02
            -9.32499466e-06 * tc[1]
            +2.16376321e-09 * tc[2]
            -1.66948050e-13 * tc[3];
    }

    /*species with midpoint at T=1398 kelvin */
    if (T < 1398) {
        /*species 25: C4H6 */
        species[25] =
            +4.78706062e-02
            -8.30893600e-05 * tc[1]
            +5.74648656e-08 * tc[2]
            -1.42863403e-11 * tc[3];
    } else {
        /*species 25: C4H6 */
        species[25] =
            +1.37163965e-02
            -9.39431566e-06 * tc[1]
            +2.18908151e-09 * tc[2]
            -1.69394481e-13 * tc[3];
    }

    /*species with midpoint at T=1392 kelvin */
    if (T < 1392) {
        /*species 26: C4H7 */
        species[26] =
            +4.26511243e-02
            -5.81958746e-05 * tc[1]
            +3.16211742e-08 * tc[2]
            -6.40239416e-12 * tc[3];
        /*species 27: C4H8-1 */
        species[27] =
            +4.52580978e-02
            -5.87317118e-05 * tc[1]
            +3.00661308e-08 * tc[2]
            -5.72766720e-12 * tc[3];
        /*species 31: C5H9 */
        species[31] =
            +5.57608487e-02
            -7.40287856e-05 * tc[1]
            +3.80651703e-08 * tc[2]
            -7.14155340e-12 * tc[3];
        /*species 32: C5H10-1 */
        species[32] =
            +5.74218294e-02
            -7.48973780e-05 * tc[1]
            +3.82094967e-08 * tc[2]
            -7.18439156e-12 * tc[3];
    } else {
        /*species 26: C4H7 */
        species[26] =
            +1.60483196e-02
            -1.09300458e-05 * tc[1]
            +2.53782316e-09 * tc[2]
            -1.95909096e-13 * tc[3];
        /*species 27: C4H8-1 */
        species[27] =
            +1.80617877e-02
            -1.23218606e-05 * tc[1]
            +2.86395888e-09 * tc[2]
            -2.21235856e-13 * tc[3];
        /*species 31: C5H9 */
        species[31] =
            +2.07128899e-02
            -1.41392123e-05 * tc[1]
            +3.28821399e-09 * tc[2]
            -2.54128883e-13 * tc[3];
        /*species 32: C5H10-1 */
        species[32] =
            +2.24072471e-02
            -1.52669605e-05 * tc[1]
            +3.54566898e-09 * tc[2]
            -2.73754056e-13 * tc[3];
    }

    /*species with midpoint at T=1395 kelvin */
    if (T < 1395) {
        /*species 28: PC4H9 */
        species[28] =
            +4.78972364e-02
            -6.28046318e-05 * tc[1]
            +3.29359416e-08 * tc[2]
            -6.48042656e-12 * tc[3];
    } else {
        /*species 28: PC4H9 */
        species[28] =
            +1.94310717e-02
            -1.32315590e-05 * tc[1]
            +3.07125408e-09 * tc[2]
            -2.37011883e-13 * tc[3];
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 29: CH3COCH2 */
        species[29] =
            +3.01407085e-02
            -3.87011104e-05 * tc[1]
            +1.91459710e-08 * tc[2]
            -3.46441272e-12 * tc[3];
        /*species 45: NC7H16 */
        species[45] =
            +8.54355820e-02
            -1.05069357e-04 * tc[1]
            +4.88837163e-08 * tc[2]
            -8.09579700e-12 * tc[3];
    } else {
        /*species 29: CH3COCH2 */
        species[29] =
            +1.16494161e-02
            -8.02011074e-06 * tc[1]
            +1.87561574e-09 * tc[2]
            -1.45513745e-13 * tc[3];
        /*species 45: NC7H16 */
        species[45] =
            +3.47675750e-02
            -2.36814258e-05 * tc[1]
            +5.49895434e-09 * tc[2]
            -4.24521064e-13 * tc[3];
    }

    /*species with midpoint at T=1378 kelvin */
    if (T < 1378) {
        /*species 30: C2H5CHO */
        species[30] =
            +2.95501264e-02
            -3.04892504e-05 * tc[1]
            +1.04851184e-08 * tc[2]
            -9.55586508e-13 * tc[3];
        /*species 40: NC3H7CHO */
        species[40] =
            +4.19240315e-02
            -4.70297558e-05 * tc[1]
            +1.88074102e-08 * tc[2]
            -2.43777563e-12 * tc[3];
    } else {
        /*species 30: C2H5CHO */
        species[30] =
            +1.39641989e-02
            -9.52496002e-06 * tc[1]
            +2.21431712e-09 * tc[2]
            -1.71103801e-13 * tc[3];
        /*species 40: NC3H7CHO */
        species[40] =
            +1.81652474e-02
            -1.23568892e-05 * tc[1]
            +2.86794062e-09 * tc[2]
            -2.21377183e-13 * tc[3];
    }

    /*species with midpoint at T=1397 kelvin */
    if (T < 1397) {
        /*species 33: C5H11-1 */
        species[33] =
            +6.10632852e-02
            -8.18983650e-05 * tc[1]
            +4.38280410e-08 * tc[2]
            -8.75438460e-12 * tc[3];
        /*species 38: C3H5-A */
        species[38] =
            +3.34559100e-02
            -5.06802054e-05 * tc[1]
            +3.08597262e-08 * tc[2]
            -6.93033360e-12 * tc[3];
    } else {
        /*species 33: C5H11-1 */
        species[33] =
            +2.39041200e-02
            -1.62954324e-05 * tc[1]
            +3.78528708e-09 * tc[2]
            -2.92270934e-13 * tc[3];
        /*species 38: C3H5-A */
        species[38] =
            +1.12695483e-02
            -7.67585728e-06 * tc[1]
            +1.78217736e-09 * tc[2]
            -1.37567212e-13 * tc[3];
    }

    /*species with midpoint at T=1385 kelvin */
    if (T < 1385) {
        /*species 34: CH3O2 */
        species[34] =
            +1.00873599e-02
            -6.43012368e-06 * tc[1]
            +6.28227801e-10 * tc[2]
            +1.67335641e-13 * tc[3];
    } else {
        /*species 34: CH3O2 */
        species[34] =
            +7.90728626e-03
            -5.36492468e-06 * tc[1]
            +1.24167401e-09 * tc[2]
            -9.56029320e-14 * tc[3];
    }

    /*species with midpoint at T=1402 kelvin */
    if (T < 1402) {
        /*species 36: C2H3CO */
        species[36] =
            +3.15273972e-02
            -6.00437870e-05 * tc[1]
            +4.44501336e-08 * tc[2]
            -1.15188612e-11 * tc[3];
    } else {
        /*species 36: C2H3CO */
        species[36] =
            +7.91296900e-03
            -5.34396560e-06 * tc[1]
            +1.23334629e-09 * tc[2]
            -9.47915924e-14 * tc[3];
    }

    /*species with midpoint at T=1393 kelvin */
    if (T < 1393) {
        /*species 37: C2H3CHO */
        species[37] =
            +3.54321417e-02
            -5.89872648e-05 * tc[1]
            +3.84300372e-08 * tc[2]
            -9.04576432e-12 * tc[3];
    } else {
        /*species 37: C2H3CHO */
        species[37] =
            +9.48963321e-03
            -6.58621058e-06 * tc[1]
            +1.54883761e-09 * tc[2]
            -1.20634916e-13 * tc[3];
    }

    /*species with midpoint at T=1383 kelvin */
    if (T < 1383) {
        /*species 41: C2H5COCH2 */
        species[41] =
            +4.39486258e-02
            -5.94004842e-05 * tc[1]
            +3.16485939e-08 * tc[2]
            -6.34395076e-12 * tc[3];
        /*species 43: NC3H7COCH3 */
        species[43] =
            +5.42494047e-02
            -6.29650664e-05 * tc[1]
            +2.74347696e-08 * tc[2]
            -4.30662264e-12 * tc[3];
    } else {
        /*species 41: C2H5COCH2 */
        species[41] =
            +1.57866459e-02
            -1.10105837e-05 * tc[1]
            +2.59761162e-09 * tc[2]
            -2.02765332e-13 * tc[3];
        /*species 43: NC3H7COCH3 */
        species[43] =
            +2.35245881e-02
            -1.61935011e-05 * tc[1]
            +3.78692070e-09 * tc[2]
            -2.93795136e-13 * tc[3];
    }

    /*species with midpoint at T=1389 kelvin */
    if (T < 1389) {
        /*species 44: NC3H7COCH2 */
        species[44] =
            +5.68162532e-02
            -7.98225562e-05 * tc[1]
            +4.58014542e-08 * tc[2]
            -9.96884188e-12 * tc[3];
    } else {
        /*species 44: NC3H7COCH2 */
        species[44] =
            +2.14093466e-02
            -1.47211923e-05 * tc[1]
            +3.43970448e-09 * tc[2]
            -2.66685089e-13 * tc[3];
    }

    /*species with midpoint at T=1399 kelvin */
    if (T < 1399) {
        /*species 49: C7H14O */
        species[49] =
            +1.01846577e-01
            -1.52009247e-04 * tc[1]
            +8.89614060e-08 * tc[2]
            -1.89941748e-11 * tc[3];
    } else {
        /*species 49: C7H14O */
        species[49] =
            +3.25585391e-02
            -2.23249812e-05 * tc[1]
            +5.20720974e-09 * tc[2]
            -4.03252860e-13 * tc[3];
    }

    /*species with midpoint at T=1387 kelvin */
    if (T < 1387) {
        /*species 50: NC7KET */
        species[50] =
            +9.58173466e-02
            -1.39337704e-04 * tc[1]
            +8.08621146e-08 * tc[2]
            -1.75491250e-11 * tc[3];
    } else {
        /*species 50: NC7KET */
        species[50] =
            +3.20168096e-02
            -2.23016912e-05 * tc[1]
            +5.25678477e-09 * tc[2]
            -4.10081804e-13 * tc[3];
    }
    return;
}

/*compute an approx to the reaction Jacobian */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void aJacobian_precond(amrex::Real *  J, amrex::Real *  sc, amrex::Real T, int HP)
{
    for (int i=0; i<2809; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[52];
    for (int k=0; k<52; k++) {
        wdot[k] = 0.0;
    }

    amrex::Real tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    amrex::Real invT = 1.0 / tc[1];
    amrex::Real invT2 = invT * invT;

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    amrex::Real refC = 101325 / 8.31446 / T;
    amrex::Real refCinv = 1.0 / refC;

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int k = 0; k < 52; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[52];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[52];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[52];
    amrex::Real Pr, fPr, F, k_0, logPr;
    amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
    amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
    amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const amrex::Real ln10 = log(10.0);
    const amrex::Real log10e = 1.0/log(10.0);
    amrex::Real c_R[52], dcRdT[52], e_RT[52];
    amrex::Real * eh_RT;
    if (HP) {
        cp_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        eh_RT = &h_RT[0];
    }
    else {
        cv_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        speciesInternalEnergy(e_RT, tc);
        eh_RT = &e_RT[0];
    }

    amrex::Real cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;
    for (int k = 0; k < 52; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[2756+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 52; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 52; ++m) {
            dehmixdc += eh_RT[m]*J[k*53+m];
        }
        J[k*53+52] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[2808] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}

/*compute an approx to the reaction Jacobian (for preconditioning) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void DWDOT_SIMPLIFIED(amrex::Real *  J, amrex::Real *  sc, amrex::Real *  Tp, int * HP)
{
    amrex::Real c[52];

    for (int k=0; k<52; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian_precond(J, c, *Tp, *HP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<52; k++) {
        J[2756+k] *= 1.e-6;
        J[k*53+52] *= 1.e6;
    }

    return;
}

/*compute the reaction Jacobian on GPU */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void aJacobian(amrex::Real * J, amrex::Real * sc, amrex::Real T, int consP)
{

    for (int i=0; i<2809; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[52];
    for (int k=0; k<52; k++) {
        wdot[k] = 0.0;
    }

    amrex::Real tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    amrex::Real invT = 1.0 / tc[1];
    amrex::Real invT2 = invT * invT;


    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    amrex::Real refC = 101325 / 8.31446 / T;
    amrex::Real refCinv = 1.0 / refC;

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int k = 0; k < 52; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[52];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[52];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[52];
    amrex::Real Pr, fPr, F, k_0, logPr;
    amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
    amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
    amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const amrex::Real ln10 = log(10.0);
    const amrex::Real log10e = 1.0/log(10.0);
    amrex::Real c_R[52], dcRdT[52], e_RT[52];
    amrex::Real * eh_RT;
    if (consP) {
        cp_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        eh_RT = &h_RT[0];
    }
    else {
        cv_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        speciesInternalEnergy(e_RT, tc);
        eh_RT = &e_RT[0];
    }

    amrex::Real cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;
    for (int k = 0; k < 52; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[2756+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 52; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 52; ++m) {
            dehmixdc += eh_RT[m]*J[k*53+m];
        }
        J[k*53+52] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[2808] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;

return;
}


/*compute the reaction Jacobian */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void DWDOT(amrex::Real *  J, amrex::Real *  sc, amrex::Real *  Tp, int * consP)
{
    amrex::Real c[52];

    for (int k=0; k<52; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<52; k++) {
        J[2756+k] *= 1.e-6;
        J[k*53+52] *= 1.e6;
    }

    return;
}


/* Transport function declarations  */


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetLENIMC(int* LENIMC ) {
    *LENIMC = 210;}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetLENRMC(int* LENRMC ) {
    *LENRMC = 53196;}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetNO(int* NO ) {
    *NO = 4;}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetKK(int* KK ) {
    *KK = 52;}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetNLITE(int* NLITE ) {
    *NLITE = 2;}


/*Patm in ergs/cm3 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetPATM(amrex::Real* PATM) {
    *PATM =   0.1013250000000000E+07;}


/*the molecular weights in g/mol */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetWT(amrex::Real* WT ) {
    WT[0] = 1.00797000E+00;
    WT[1] = 2.01594000E+00;
    WT[2] = 1.59994000E+01;
    WT[3] = 3.19988000E+01;
    WT[4] = 1.70073700E+01;
    WT[5] = 1.80153400E+01;
    WT[6] = 2.80105500E+01;
    WT[7] = 4.40099500E+01;
    WT[8] = 1.50350600E+01;
    WT[9] = 1.60430300E+01;
    WT[10] = 3.30067700E+01;
    WT[11] = 3.40147400E+01;
    WT[12] = 3.00264900E+01;
    WT[13] = 3.00701200E+01;
    WT[14] = 2.80541800E+01;
    WT[15] = 2.90621500E+01;
    WT[16] = 2.50302700E+01;
    WT[17] = 2.60382400E+01;
    WT[18] = 3.20424300E+01;
    WT[19] = 4.20376400E+01;
    WT[20] = 4.10296700E+01;
    WT[21] = 4.40535800E+01;
    WT[22] = 4.00653300E+01;
    WT[23] = 4.00653300E+01;
    WT[24] = 4.20812700E+01;
    WT[25] = 5.40924200E+01;
    WT[26] = 5.51003900E+01;
    WT[27] = 5.61083600E+01;
    WT[28] = 5.71163300E+01;
    WT[29] = 5.70727000E+01;
    WT[30] = 5.80806700E+01;
    WT[31] = 6.91274800E+01;
    WT[32] = 7.01354500E+01;
    WT[33] = 7.11434200E+01;
    WT[34] = 4.70338600E+01;
    WT[35] = 4.80418300E+01;
    WT[36] = 5.50567600E+01;
    WT[37] = 5.60647300E+01;
    WT[38] = 4.10733000E+01;
    WT[39] = 3.90573600E+01;
    WT[40] = 7.21077600E+01;
    WT[41] = 7.10997900E+01;
    WT[42] = 5.60647300E+01;
    WT[43] = 8.61348500E+01;
    WT[44] = 8.51268800E+01;
    WT[45] = 1.00205570E+02;
    WT[46] = 9.91976000E+01;
    WT[47] = 1.31196400E+02;
    WT[48] = 1.63195200E+02;
    WT[49] = 1.14189030E+02;
    WT[50] = 1.46187830E+02;
    WT[51] = 2.80134000E+01;
}


/*the lennard-jones potential well depth eps/kb in K */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetEPS(amrex::Real* EPS ) {
    EPS[0] = 1.45000000E+02;
    EPS[1] = 3.80000000E+01;
    EPS[2] = 8.00000000E+01;
    EPS[3] = 1.07400000E+02;
    EPS[4] = 8.00000000E+01;
    EPS[5] = 5.72400000E+02;
    EPS[6] = 9.81000000E+01;
    EPS[7] = 2.44000000E+02;
    EPS[8] = 1.44000000E+02;
    EPS[9] = 1.41400000E+02;
    EPS[10] = 1.07400000E+02;
    EPS[11] = 1.07400000E+02;
    EPS[12] = 4.98000000E+02;
    EPS[13] = 2.47500000E+02;
    EPS[14] = 2.38400000E+02;
    EPS[15] = 2.47500000E+02;
    EPS[16] = 2.65300000E+02;
    EPS[17] = 2.65300000E+02;
    EPS[18] = 4.81800000E+02;
    EPS[19] = 4.36000000E+02;
    EPS[20] = 1.50000000E+02;
    EPS[21] = 4.36000000E+02;
    EPS[22] = 3.24800000E+02;
    EPS[23] = 3.24800000E+02;
    EPS[24] = 3.07800000E+02;
    EPS[25] = 3.57000000E+02;
    EPS[26] = 3.55000000E+02;
    EPS[27] = 3.55000000E+02;
    EPS[28] = 3.52000000E+02;
    EPS[29] = 4.35500000E+02;
    EPS[30] = 4.35200000E+02;
    EPS[31] = 3.96800000E+02;
    EPS[32] = 3.86200000E+02;
    EPS[33] = 4.40735000E+02;
    EPS[34] = 4.81800000E+02;
    EPS[35] = 4.81800000E+02;
    EPS[36] = 4.43200000E+02;
    EPS[37] = 4.28800000E+02;
    EPS[38] = 3.16000000E+02;
    EPS[39] = 3.24800000E+02;
    EPS[40] = 4.64200000E+02;
    EPS[41] = 4.54000000E+02;
    EPS[42] = 4.43200000E+02;
    EPS[43] = 4.77000000E+02;
    EPS[44] = 4.89084000E+02;
    EPS[45] = 4.59600000E+02;
    EPS[46] = 4.59600000E+02;
    EPS[47] = 5.61000000E+02;
    EPS[48] = 6.00600000E+02;
    EPS[49] = 5.11500000E+02;
    EPS[50] = 5.81300000E+02;
    EPS[51] = 9.75300000E+01;
}


/*the lennard-jones collision diameter in Angstroms */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetSIG(amrex::Real* SIG ) {
    SIG[0] = 2.05000000E+00;
    SIG[1] = 2.92000000E+00;
    SIG[2] = 2.75000000E+00;
    SIG[3] = 3.45800000E+00;
    SIG[4] = 2.75000000E+00;
    SIG[5] = 2.60500000E+00;
    SIG[6] = 3.65000000E+00;
    SIG[7] = 3.76300000E+00;
    SIG[8] = 3.80000000E+00;
    SIG[9] = 3.74600000E+00;
    SIG[10] = 3.45800000E+00;
    SIG[11] = 3.45800000E+00;
    SIG[12] = 3.59000000E+00;
    SIG[13] = 4.35000000E+00;
    SIG[14] = 3.49600000E+00;
    SIG[15] = 4.35000000E+00;
    SIG[16] = 3.72100000E+00;
    SIG[17] = 3.72100000E+00;
    SIG[18] = 3.62600000E+00;
    SIG[19] = 3.97000000E+00;
    SIG[20] = 2.50000000E+00;
    SIG[21] = 3.97000000E+00;
    SIG[22] = 4.29000000E+00;
    SIG[23] = 4.29000000E+00;
    SIG[24] = 4.14000000E+00;
    SIG[25] = 4.72000000E+00;
    SIG[26] = 4.65000000E+00;
    SIG[27] = 4.65000000E+00;
    SIG[28] = 5.24000000E+00;
    SIG[29] = 4.86000000E+00;
    SIG[30] = 4.66200000E+00;
    SIG[31] = 5.45800000E+00;
    SIG[32] = 5.48900000E+00;
    SIG[33] = 5.04100000E+00;
    SIG[34] = 3.62600000E+00;
    SIG[35] = 3.62600000E+00;
    SIG[36] = 4.12000000E+00;
    SIG[37] = 4.95800000E+00;
    SIG[38] = 4.22000000E+00;
    SIG[39] = 4.29000000E+00;
    SIG[40] = 5.00900000E+00;
    SIG[41] = 5.41300000E+00;
    SIG[42] = 4.12000000E+00;
    SIG[43] = 5.73500000E+00;
    SIG[44] = 5.34800000E+00;
    SIG[45] = 6.25300000E+00;
    SIG[46] = 6.25300000E+00;
    SIG[47] = 6.31700000E+00;
    SIG[48] = 7.22900000E+00;
    SIG[49] = 6.29700000E+00;
    SIG[50] = 6.50600000E+00;
    SIG[51] = 3.62100000E+00;
}


/*the dipole moment in Debye */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetDIP(amrex::Real* DIP ) {
    DIP[0] = 0.00000000E+00;
    DIP[1] = 0.00000000E+00;
    DIP[2] = 0.00000000E+00;
    DIP[3] = 0.00000000E+00;
    DIP[4] = 0.00000000E+00;
    DIP[5] = 1.84400000E+00;
    DIP[6] = 0.00000000E+00;
    DIP[7] = 0.00000000E+00;
    DIP[8] = 0.00000000E+00;
    DIP[9] = 0.00000000E+00;
    DIP[10] = 0.00000000E+00;
    DIP[11] = 0.00000000E+00;
    DIP[12] = 0.00000000E+00;
    DIP[13] = 0.00000000E+00;
    DIP[14] = 0.00000000E+00;
    DIP[15] = 0.00000000E+00;
    DIP[16] = 0.00000000E+00;
    DIP[17] = 0.00000000E+00;
    DIP[18] = 0.00000000E+00;
    DIP[19] = 0.00000000E+00;
    DIP[20] = 0.00000000E+00;
    DIP[21] = 0.00000000E+00;
    DIP[22] = 0.00000000E+00;
    DIP[23] = 0.00000000E+00;
    DIP[24] = 0.00000000E+00;
    DIP[25] = 0.00000000E+00;
    DIP[26] = 0.00000000E+00;
    DIP[27] = 0.00000000E+00;
    DIP[28] = 0.00000000E+00;
    DIP[29] = 0.00000000E+00;
    DIP[30] = 2.70000000E+00;
    DIP[31] = 0.00000000E+00;
    DIP[32] = 4.00000000E-01;
    DIP[33] = 0.00000000E+00;
    DIP[34] = 0.00000000E+00;
    DIP[35] = 0.00000000E+00;
    DIP[36] = 0.00000000E+00;
    DIP[37] = 2.90000000E+00;
    DIP[38] = 0.00000000E+00;
    DIP[39] = 0.00000000E+00;
    DIP[40] = 2.60000000E+00;
    DIP[41] = 3.30000000E+00;
    DIP[42] = 0.00000000E+00;
    DIP[43] = 2.50000000E+00;
    DIP[44] = 0.00000000E+00;
    DIP[45] = 0.00000000E+00;
    DIP[46] = 0.00000000E+00;
    DIP[47] = 1.70000000E+00;
    DIP[48] = 1.80000000E+00;
    DIP[49] = 0.00000000E+00;
    DIP[50] = 2.00000000E+00;
    DIP[51] = 0.00000000E+00;
}


/*the polarizability in cubic Angstroms */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetPOL(amrex::Real* POL ) {
    POL[0] = 0.00000000E+00;
    POL[1] = 7.90000000E-01;
    POL[2] = 0.00000000E+00;
    POL[3] = 1.60000000E+00;
    POL[4] = 0.00000000E+00;
    POL[5] = 0.00000000E+00;
    POL[6] = 1.95000000E+00;
    POL[7] = 2.65000000E+00;
    POL[8] = 0.00000000E+00;
    POL[9] = 2.60000000E+00;
    POL[10] = 0.00000000E+00;
    POL[11] = 0.00000000E+00;
    POL[12] = 0.00000000E+00;
    POL[13] = 0.00000000E+00;
    POL[14] = 0.00000000E+00;
    POL[15] = 0.00000000E+00;
    POL[16] = 0.00000000E+00;
    POL[17] = 0.00000000E+00;
    POL[18] = 0.00000000E+00;
    POL[19] = 0.00000000E+00;
    POL[20] = 0.00000000E+00;
    POL[21] = 0.00000000E+00;
    POL[22] = 0.00000000E+00;
    POL[23] = 0.00000000E+00;
    POL[24] = 0.00000000E+00;
    POL[25] = 0.00000000E+00;
    POL[26] = 0.00000000E+00;
    POL[27] = 0.00000000E+00;
    POL[28] = 0.00000000E+00;
    POL[29] = 0.00000000E+00;
    POL[30] = 0.00000000E+00;
    POL[31] = 0.00000000E+00;
    POL[32] = 0.00000000E+00;
    POL[33] = 0.00000000E+00;
    POL[34] = 0.00000000E+00;
    POL[35] = 0.00000000E+00;
    POL[36] = 0.00000000E+00;
    POL[37] = 0.00000000E+00;
    POL[38] = 0.00000000E+00;
    POL[39] = 0.00000000E+00;
    POL[40] = 0.00000000E+00;
    POL[41] = 0.00000000E+00;
    POL[42] = 0.00000000E+00;
    POL[43] = 0.00000000E+00;
    POL[44] = 0.00000000E+00;
    POL[45] = 0.00000000E+00;
    POL[46] = 0.00000000E+00;
    POL[47] = 0.00000000E+00;
    POL[48] = 0.00000000E+00;
    POL[49] = 0.00000000E+00;
    POL[50] = 0.00000000E+00;
    POL[51] = 1.76000000E+00;
}


/*the rotational relaxation collision number at 298 K */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetZROT(amrex::Real* ZROT ) {
    ZROT[0] = 0.00000000E+00;
    ZROT[1] = 2.80000000E+02;
    ZROT[2] = 0.00000000E+00;
    ZROT[3] = 3.80000000E+00;
    ZROT[4] = 0.00000000E+00;
    ZROT[5] = 4.00000000E+00;
    ZROT[6] = 1.80000000E+00;
    ZROT[7] = 2.10000000E+00;
    ZROT[8] = 0.00000000E+00;
    ZROT[9] = 1.30000000E+01;
    ZROT[10] = 1.00000000E+00;
    ZROT[11] = 3.80000000E+00;
    ZROT[12] = 2.00000000E+00;
    ZROT[13] = 1.50000000E+00;
    ZROT[14] = 1.50000000E+00;
    ZROT[15] = 1.50000000E+00;
    ZROT[16] = 2.50000000E+00;
    ZROT[17] = 2.50000000E+00;
    ZROT[18] = 1.00000000E+00;
    ZROT[19] = 2.00000000E+00;
    ZROT[20] = 1.00000000E+00;
    ZROT[21] = 2.00000000E+00;
    ZROT[22] = 1.00000000E+00;
    ZROT[23] = 1.00000000E+00;
    ZROT[24] = 1.00000000E+00;
    ZROT[25] = 1.00000000E+00;
    ZROT[26] = 1.00000000E+00;
    ZROT[27] = 1.00000000E+00;
    ZROT[28] = 1.00000000E+00;
    ZROT[29] = 1.00000000E+00;
    ZROT[30] = 1.00000000E+00;
    ZROT[31] = 1.00000000E+00;
    ZROT[32] = 1.00000000E+00;
    ZROT[33] = 0.00000000E+00;
    ZROT[34] = 1.00000000E+00;
    ZROT[35] = 1.00000000E+00;
    ZROT[36] = 1.00000000E+00;
    ZROT[37] = 1.00000000E+00;
    ZROT[38] = 1.00000000E+00;
    ZROT[39] = 1.00000000E+00;
    ZROT[40] = 1.00000000E+00;
    ZROT[41] = 1.00000000E+00;
    ZROT[42] = 1.00000000E+00;
    ZROT[43] = 1.00000000E+00;
    ZROT[44] = 0.00000000E+00;
    ZROT[45] = 1.00000000E+00;
    ZROT[46] = 1.00000000E+00;
    ZROT[47] = 1.00000000E+00;
    ZROT[48] = 1.00000000E+00;
    ZROT[49] = 1.00000000E+00;
    ZROT[50] = 1.00000000E+00;
    ZROT[51] = 4.00000000E+00;
}


/*0: monoatomic, 1: linear, 2: nonlinear */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetNLIN(int* NLIN) {
    NLIN[0] = 0;
    NLIN[1] = 1;
    NLIN[2] = 0;
    NLIN[3] = 1;
    NLIN[4] = 1;
    NLIN[5] = 2;
    NLIN[6] = 1;
    NLIN[7] = 1;
    NLIN[8] = 1;
    NLIN[9] = 2;
    NLIN[10] = 2;
    NLIN[11] = 2;
    NLIN[12] = 2;
    NLIN[13] = 2;
    NLIN[14] = 2;
    NLIN[15] = 2;
    NLIN[16] = 1;
    NLIN[17] = 1;
    NLIN[18] = 2;
    NLIN[19] = 2;
    NLIN[20] = 2;
    NLIN[21] = 2;
    NLIN[22] = 1;
    NLIN[23] = 1;
    NLIN[24] = 2;
    NLIN[25] = 2;
    NLIN[26] = 2;
    NLIN[27] = 2;
    NLIN[28] = 2;
    NLIN[29] = 2;
    NLIN[30] = 2;
    NLIN[31] = 2;
    NLIN[32] = 2;
    NLIN[33] = 2;
    NLIN[34] = 2;
    NLIN[35] = 2;
    NLIN[36] = 2;
    NLIN[37] = 2;
    NLIN[38] = 2;
    NLIN[39] = 1;
    NLIN[40] = 2;
    NLIN[41] = 2;
    NLIN[42] = 2;
    NLIN[43] = 2;
    NLIN[44] = 2;
    NLIN[45] = 2;
    NLIN[46] = 2;
    NLIN[47] = 2;
    NLIN[48] = 2;
    NLIN[49] = 2;
    NLIN[50] = 2;
    NLIN[51] = 1;
}


/*Poly fits for the viscosities, dim NO*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFETA(amrex::Real* COFETA) {
    COFETA[0] = -1.95281568E+01;
    COFETA[1] = 3.26212115E+00;
    COFETA[2] = -3.40404685E-01;
    COFETA[3] = 1.47532537E-02;
    COFETA[4] = -1.34911581E+01;
    COFETA[5] = 8.49217713E-01;
    COFETA[6] = -2.75872442E-02;
    COFETA[7] = 1.23722694E-03;
    COFETA[8] = -1.45327736E+01;
    COFETA[9] = 1.67628298E+00;
    COFETA[10] = -1.36586011E-01;
    COFETA[11] = 6.02009005E-03;
    COFETA[12] = -1.65122758E+01;
    COFETA[13] = 2.39046525E+00;
    COFETA[14] = -2.29715738E-01;
    COFETA[15] = 1.00664463E-02;
    COFETA[16] = -1.45022258E+01;
    COFETA[17] = 1.67628298E+00;
    COFETA[18] = -1.36586011E-01;
    COFETA[19] = 6.02009005E-03;
    COFETA[20] = -1.29340489E+01;
    COFETA[21] = -3.14082959E-01;
    COFETA[22] = 2.63797240E-01;
    COFETA[23] = -1.64132505E-02;
    COFETA[24] = -1.60225925E+01;
    COFETA[25] = 2.13662660E+00;
    COFETA[26] = -1.96738977E-01;
    COFETA[27] = 8.63898199E-03;
    COFETA[28] = -2.33572324E+01;
    COFETA[29] = 4.85724263E+00;
    COFETA[30] = -5.32102483E-01;
    COFETA[31] = 2.24362161E-02;
    COFETA[32] = -1.93615554E+01;
    COFETA[33] = 3.24398914E+00;
    COFETA[34] = -3.38151420E-01;
    COFETA[35] = 1.46598067E-02;
    COFETA[36] = -1.91726480E+01;
    COFETA[37] = 3.19752754E+00;
    COFETA[38] = -3.32398830E-01;
    COFETA[39] = 1.44222031E-02;
    COFETA[40] = -1.64967686E+01;
    COFETA[41] = 2.39046525E+00;
    COFETA[42] = -2.29715738E-01;
    COFETA[43] = 1.00664463E-02;
    COFETA[44] = -1.64817280E+01;
    COFETA[45] = 2.39046525E+00;
    COFETA[46] = -2.29715738E-01;
    COFETA[47] = 1.00664463E-02;
    COFETA[48] = -2.20006672E+01;
    COFETA[49] = 3.66530272E+00;
    COFETA[50] = -3.09756445E-01;
    COFETA[51] = 1.02090749E-02;
    COFETA[52] = -2.39150298E+01;
    COFETA[53] = 4.87957850E+00;
    COFETA[54] = -5.34166762E-01;
    COFETA[55] = 2.24923677E-02;
    COFETA[56] = -2.32853624E+01;
    COFETA[57] = 4.81068153E+00;
    COFETA[58] = -5.27301048E-01;
    COFETA[59] = 2.22776802E-02;
    COFETA[60] = -2.39320775E+01;
    COFETA[61] = 4.87957850E+00;
    COFETA[62] = -5.34166762E-01;
    COFETA[63] = 2.24923677E-02;
    COFETA[64] = -2.42271099E+01;
    COFETA[65] = 5.05175437E+00;
    COFETA[66] = -5.52836130E-01;
    COFETA[67] = 2.31559377E-02;
    COFETA[68] = -2.42073697E+01;
    COFETA[69] = 5.05175437E+00;
    COFETA[70] = -5.52836130E-01;
    COFETA[71] = 2.31559377E-02;
    COFETA[72] = -2.24921638E+01;
    COFETA[73] = 3.90337385E+00;
    COFETA[74] = -3.45176993E-01;
    COFETA[75] = 1.19122568E-02;
    COFETA[76] = -2.37089384E+01;
    COFETA[77] = 4.47860513E+00;
    COFETA[78] = -4.32776310E-01;
    COFETA[79] = 1.61943631E-02;
    COFETA[80] = -1.83221188E+01;
    COFETA[81] = 3.35360209E+00;
    COFETA[82] = -3.51788995E-01;
    COFETA[83] = 1.52259804E-02;
    COFETA[84] = -2.36855178E+01;
    COFETA[85] = 4.47860513E+00;
    COFETA[86] = -4.32776310E-01;
    COFETA[87] = 1.61943631E-02;
    COFETA[88] = -2.49619466E+01;
    COFETA[89] = 5.17942495E+00;
    COFETA[90] = -5.54420660E-01;
    COFETA[91] = 2.26381154E-02;
    COFETA[92] = -2.49619466E+01;
    COFETA[93] = 5.17942495E+00;
    COFETA[94] = -5.54420660E-01;
    COFETA[95] = 2.26381154E-02;
    COFETA[96] = -2.47298373E+01;
    COFETA[97] = 5.16628217E+00;
    COFETA[98] = -5.57016710E-01;
    COFETA[99] = 2.29176713E-02;
    COFETA[100] = -2.49972039E+01;
    COFETA[101] = 5.09673056E+00;
    COFETA[102] = -5.34946406E-01;
    COFETA[103] = 2.14562401E-02;
    COFETA[104] = -2.49657309E+01;
    COFETA[105] = 5.10453889E+00;
    COFETA[106] = -5.36480276E-01;
    COFETA[107] = 2.15426062E-02;
    COFETA[108] = -2.49566669E+01;
    COFETA[109] = 5.10453889E+00;
    COFETA[110] = -5.36480276E-01;
    COFETA[111] = 2.15426062E-02;
    COFETA[112] = -2.51973546E+01;
    COFETA[113] = 5.11598056E+00;
    COFETA[114] = -5.38749098E-01;
    COFETA[115] = 2.16708786E-02;
    COFETA[116] = -2.39706955E+01;
    COFETA[117] = 4.48381911E+00;
    COFETA[118] = -4.33592214E-01;
    COFETA[119] = 1.62349637E-02;
    COFETA[120] = -2.07232330E+01;
    COFETA[121] = 3.09576972E+00;
    COFETA[122] = -2.36384540E-01;
    COFETA[123] = 7.09310859E-03;
    COFETA[124] = -2.47825183E+01;
    COFETA[125] = 4.84621755E+00;
    COFETA[126] = -4.91412112E-01;
    COFETA[127] = 1.91489049E-02;
    COFETA[128] = -2.49010581E+01;
    COFETA[129] = 4.91845648E+00;
    COFETA[130] = -5.03809621E-01;
    COFETA[131] = 1.98018057E-02;
    COFETA[132] = -2.38227845E+01;
    COFETA[133] = 4.42705793E+00;
    COFETA[134] = -4.24748121E-01;
    COFETA[135] = 1.57959821E-02;
    COFETA[136] = -2.23002604E+01;
    COFETA[137] = 3.90337385E+00;
    COFETA[138] = -3.45176993E-01;
    COFETA[139] = 1.19122568E-02;
    COFETA[140] = -2.22896583E+01;
    COFETA[141] = 3.90337385E+00;
    COFETA[142] = -3.45176993E-01;
    COFETA[143] = 1.19122568E-02;
    COFETA[144] = -2.34917309E+01;
    COFETA[145] = 4.39885838E+00;
    COFETA[146] = -4.20378963E-01;
    COFETA[147] = 1.55798335E-02;
    COFETA[148] = -2.11256515E+01;
    COFETA[149] = 3.22130837E+00;
    COFETA[150] = -2.55162995E-01;
    COFETA[151] = 7.99590805E-03;
    COFETA[152] = -2.48372497E+01;
    COFETA[153] = 5.16888194E+00;
    COFETA[154] = -5.55239216E-01;
    COFETA[155] = 2.27583887E-02;
    COFETA[156] = -2.49746866E+01;
    COFETA[157] = 5.17942495E+00;
    COFETA[158] = -5.54420660E-01;
    COFETA[159] = 2.26381154E-02;
    COFETA[160] = -2.10282452E+01;
    COFETA[161] = 3.20441952E+00;
    COFETA[162] = -2.50942950E-01;
    COFETA[163] = 7.72807867E-03;
    COFETA[164] = -2.06838690E+01;
    COFETA[165] = 2.98076056E+00;
    COFETA[166] = -2.19099978E-01;
    COFETA[167] = 6.25476429E-03;
    COFETA[168] = -2.34826598E+01;
    COFETA[169] = 4.39885838E+00;
    COFETA[170] = -4.20378963E-01;
    COFETA[171] = 1.55798335E-02;
    COFETA[172] = -2.19266754E+01;
    COFETA[173] = 3.50222222E+00;
    COFETA[174] = -2.91485443E-01;
    COFETA[175] = 9.53907459E-03;
    COFETA[176] = -2.25608335E+01;
    COFETA[177] = 3.79896931E+00;
    COFETA[178] = -3.29599952E-01;
    COFETA[179] = 1.11618162E-02;
    COFETA[180] = -2.36172714E+01;
    COFETA[181] = 4.19624204E+00;
    COFETA[182] = -3.89379568E-01;
    COFETA[183] = 1.40593024E-02;
    COFETA[184] = -2.36223264E+01;
    COFETA[185] = 4.19624204E+00;
    COFETA[186] = -3.89379568E-01;
    COFETA[187] = 1.40593024E-02;
    COFETA[188] = -2.00604654E+01;
    COFETA[189] = 2.60194078E+00;
    COFETA[190] = -1.55261704E-01;
    COFETA[191] = 2.91379503E-03;
    COFETA[192] = -1.87654031E+01;
    COFETA[193] = 1.94195749E+00;
    COFETA[194] = -5.96682391E-02;
    COFETA[195] = -1.59051690E-03;
    COFETA[196] = -2.20476793E+01;
    COFETA[197] = 3.47183888E+00;
    COFETA[198] = -2.80951710E-01;
    COFETA[199] = 8.82358958E-03;
    COFETA[200] = -1.92165798E+01;
    COFETA[201] = 2.22203568E+00;
    COFETA[202] = -1.00720551E-01;
    COFETA[203] = 3.61477459E-04;
    COFETA[204] = -1.59733786E+01;
    COFETA[205] = 2.12416328E+00;
    COFETA[206] = -1.95139471E-01;
    COFETA[207] = 8.57058122E-03;
}


/*Poly fits for the conductivities, dim NO*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFLAM(amrex::Real* COFLAM) {
    COFLAM[0] = 2.17536347E-02;
    COFLAM[1] = 3.26212115E+00;
    COFLAM[2] = -3.40404685E-01;
    COFLAM[3] = 1.47532537E-02;
    COFLAM[4] = 1.25921627E+01;
    COFLAM[5] = -1.98118570E+00;
    COFLAM[6] = 3.41585728E-01;
    COFLAM[7] = -1.37541082E-02;
    COFLAM[8] = 2.25252403E+00;
    COFLAM[9] = 1.67628298E+00;
    COFLAM[10] = -1.36586011E-01;
    COFLAM[11] = 6.02009005E-03;
    COFLAM[12] = -2.61206080E+00;
    COFLAM[13] = 3.19791630E+00;
    COFLAM[14] = -3.17109344E-01;
    COFLAM[15] = 1.38112807E-02;
    COFLAM[16] = 1.50122657E+01;
    COFLAM[17] = -3.77386870E+00;
    COFLAM[18] = 6.33322103E-01;
    COFLAM[19] = -2.91909452E-02;
    COFLAM[20] = 2.03772403E+01;
    COFLAM[21] = -7.67826677E+00;
    COFLAM[22] = 1.34610048E+00;
    COFLAM[23] = -6.73692973E-02;
    COFLAM[24] = 8.34784014E+00;
    COFLAM[25] = -1.58645601E+00;
    COFLAM[26] = 3.71139553E-01;
    COFLAM[27] = -1.90604654E-02;
    COFLAM[28] = -1.29294568E+01;
    COFLAM[29] = 6.58853915E+00;
    COFLAM[30] = -6.74340577E-01;
    COFLAM[31] = 2.55894378E-02;
    COFLAM[32] = 1.02094592E+01;
    COFLAM[33] = -2.93698993E+00;
    COFLAM[34] = 6.53896256E-01;
    COFLAM[35] = -3.52423209E-02;
    COFLAM[36] = 5.28309555E+00;
    COFLAM[37] = -1.44253852E+00;
    COFLAM[38] = 5.26102426E-01;
    COFLAM[39] = -3.24374803E-02;
    COFLAM[40] = 2.84329786E+00;
    COFLAM[41] = 4.28805451E-01;
    COFLAM[42] = 1.41908857E-01;
    COFLAM[43] = -1.03645530E-02;
    COFLAM[44] = 1.51123736E+00;
    COFLAM[45] = 1.05082768E+00;
    COFLAM[46] = 5.89158608E-02;
    COFLAM[47] = -6.47068228E-03;
    COFLAM[48] = -3.78625901E+00;
    COFLAM[49] = 1.71406876E+00;
    COFLAM[50] = 1.31217187E-01;
    COFLAM[51] = -1.60027645E-02;
    COFLAM[52] = -1.73779253E+01;
    COFLAM[53] = 7.50260709E+00;
    COFLAM[54] = -6.58047056E-01;
    COFLAM[55] = 1.96731960E-02;
    COFLAM[56] = -1.59123020E+01;
    COFLAM[57] = 7.22151973E+00;
    COFLAM[58] = -6.48384550E-01;
    COFLAM[59] = 2.03544291E-02;
    COFLAM[60] = -1.15783951E+01;
    COFLAM[61] = 5.08669844E+00;
    COFLAM[62] = -3.32743856E-01;
    COFLAM[63] = 5.17681405E-03;
    COFLAM[64] = -8.97918157E+00;
    COFLAM[65] = 5.26969187E+00;
    COFLAM[66] = -5.30903818E-01;
    COFLAM[67] = 2.11696907E-02;
    COFLAM[68] = -8.47004097E+00;
    COFLAM[69] = 4.79841607E+00;
    COFLAM[70] = -4.18356186E-01;
    COFLAM[71] = 1.40201875E-02;
    COFLAM[72] = -7.46310021E+00;
    COFLAM[73] = 3.12048160E+00;
    COFLAM[74] = -3.53572264E-02;
    COFLAM[75] = -9.26347553E-03;
    COFLAM[76] = -1.02388434E+01;
    COFLAM[77] = 4.85175579E+00;
    COFLAM[78] = -3.55720194E-01;
    COFLAM[79] = 8.26691128E-03;
    COFLAM[80] = -8.56201588E-01;
    COFLAM[81] = 2.36542065E+00;
    COFLAM[82] = -1.28936207E-01;
    COFLAM[83] = 2.10883480E-03;
    COFLAM[84] = -1.31233827E+01;
    COFLAM[85] = 5.61069079E+00;
    COFLAM[86] = -3.96021793E-01;
    COFLAM[87] = 7.70108783E-03;
    COFLAM[88] = -1.53195929E+01;
    COFLAM[89] = 6.83501133E+00;
    COFLAM[90] = -6.00074934E-01;
    COFLAM[91] = 1.83171185E-02;
    COFLAM[92] = -1.28900959E+01;
    COFLAM[93] = 5.84742944E+00;
    COFLAM[94] = -4.67171657E-01;
    COFLAM[95] = 1.23718921E-02;
    COFLAM[96] = -1.78898978E+01;
    COFLAM[97] = 7.76779880E+00;
    COFLAM[98] = -6.99718925E-01;
    COFLAM[99] = 2.17864510E-02;
    COFLAM[100] = -2.27201162E+01;
    COFLAM[101] = 9.81386769E+00;
    COFLAM[102] = -9.98452158E-01;
    COFLAM[103] = 3.60221102E-02;
    COFLAM[104] = -2.07926038E+01;
    COFLAM[105] = 8.91865701E+00;
    COFLAM[106] = -8.59588762E-01;
    COFLAM[107] = 2.91231839E-02;
    COFLAM[108] = -2.15359921E+01;
    COFLAM[109] = 9.15859129E+00;
    COFLAM[110] = -8.81784967E-01;
    COFLAM[111] = 2.97335771E-02;
    COFLAM[112] = -2.02208492E+01;
    COFLAM[113] = 8.58160889E+00;
    COFLAM[114] = -8.07480954E-01;
    COFLAM[115] = 2.65519716E-02;
    COFLAM[116] = -1.61698448E+01;
    COFLAM[117] = 6.91606213E+00;
    COFLAM[118] = -5.91956548E-01;
    COFLAM[119] = 1.72065795E-02;
    COFLAM[120] = -1.27025613E+01;
    COFLAM[121] = 5.33085865E+00;
    COFLAM[122] = -3.52133931E-01;
    COFLAM[123] = 5.57337063E-03;
    COFLAM[124] = -2.33498943E+01;
    COFLAM[125] = 9.74613300E+00;
    COFLAM[126] = -9.56820780E-01;
    COFLAM[127] = 3.29353240E-02;
    COFLAM[128] = -2.19994180E+01;
    COFLAM[129] = 9.21723676E+00;
    COFLAM[130] = -8.85870988E-01;
    COFLAM[131] = 2.97842381E-02;
    COFLAM[132] = -2.03112685E+01;
    COFLAM[133] = 8.52078401E+00;
    COFLAM[134] = -7.82685295E-01;
    COFLAM[135] = 2.48440321E-02;
    COFLAM[136] = -4.43895254E-02;
    COFLAM[137] = 4.46027685E-01;
    COFLAM[138] = 2.78709951E-01;
    COFLAM[139] = -2.15815666E-02;
    COFLAM[140] = -9.95685769E+00;
    COFLAM[141] = 4.63788986E+00;
    COFLAM[142] = -2.96039549E-01;
    COFLAM[143] = 4.40600450E-03;
    COFLAM[144] = -1.60659105E+01;
    COFLAM[145] = 7.26067097E+00;
    COFLAM[146] = -6.78836336E-01;
    COFLAM[147] = 2.27000439E-02;
    COFLAM[148] = -1.95704643E+01;
    COFLAM[149] = 8.30295520E+00;
    COFLAM[150] = -7.86163857E-01;
    COFLAM[151] = 2.61814324E-02;
    COFLAM[152] = -2.13992571E+01;
    COFLAM[153] = 9.33379574E+00;
    COFLAM[154] = -9.32773227E-01;
    COFLAM[155] = 3.30421693E-02;
    COFLAM[156] = -1.01860015E+01;
    COFLAM[157] = 5.17249797E+00;
    COFLAM[158] = -4.31613491E-01;
    COFLAM[159] = 1.28965891E-02;
    COFLAM[160] = -1.49697244E+01;
    COFLAM[161] = 6.26466263E+00;
    COFLAM[162] = -4.78896529E-01;
    COFLAM[163] = 1.12767966E-02;
    COFLAM[164] = -1.61832427E+01;
    COFLAM[165] = 6.76038051E+00;
    COFLAM[166] = -5.54437038E-01;
    COFLAM[167] = 1.49070696E-02;
    COFLAM[168] = -1.61893803E+01;
    COFLAM[169] = 7.20974544E+00;
    COFLAM[170] = -6.54511436E-01;
    COFLAM[171] = 2.09184219E-02;
    COFLAM[172] = -1.65182314E+01;
    COFLAM[173] = 6.77409685E+00;
    COFLAM[174] = -5.39725559E-01;
    COFLAM[175] = 1.36600079E-02;
    COFLAM[176] = -1.67150787E+01;
    COFLAM[177] = 6.99048809E+00;
    COFLAM[178] = -5.80985376E-01;
    COFLAM[179] = 1.59720875E-02;
    COFLAM[180] = -2.10271413E+01;
    COFLAM[181] = 8.65052441E+00;
    COFLAM[182] = -7.92479996E-01;
    COFLAM[183] = 2.49902167E-02;
    COFLAM[184] = -2.02655631E+01;
    COFLAM[185] = 8.36202130E+00;
    COFLAM[186] = -7.57796108E-01;
    COFLAM[187] = 2.35933234E-02;
    COFLAM[188] = -1.66908746E+01;
    COFLAM[189] = 6.81710536E+00;
    COFLAM[190] = -5.44157083E-01;
    COFLAM[191] = 1.38399468E-02;
    COFLAM[192] = -1.64947810E+01;
    COFLAM[193] = 6.74239206E+00;
    COFLAM[194] = -5.48574109E-01;
    COFLAM[195] = 1.45634965E-02;
    COFLAM[196] = -2.58993779E+01;
    COFLAM[197] = 1.05612702E+01;
    COFLAM[198] = -1.04820266E+00;
    COFLAM[199] = 3.63660308E-02;
    COFLAM[200] = -1.62662396E+01;
    COFLAM[201] = 6.66893722E+00;
    COFLAM[202] = -5.31770392E-01;
    COFLAM[203] = 1.35382162E-02;
    COFLAM[204] = 1.00878735E+01;
    COFLAM[205] = -2.26868946E+00;
    COFLAM[206] = 4.60600802E-01;
    COFLAM[207] = -2.29420040E-02;
}


/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFD(amrex::Real* COFD) {
    COFD[0] = -1.41218687E+01;
    COFD[1] = 3.92946227E+00;
    COFD[2] = -2.91729340E-01;
    COFD[3] = 1.24783320E-02;
    COFD[4] = -1.09419258E+01;
    COFD[5] = 2.56365104E+00;
    COFD[6] = -1.18913262E-01;
    COFD[7] = 5.17669825E-03;
    COFD[8] = -1.29709581E+01;
    COFD[9] = 3.28440910E+00;
    COFD[10] = -2.11707534E-01;
    COFD[11] = 9.16419558E-03;
    COFD[12] = -1.41183595E+01;
    COFD[13] = 3.59720720E+00;
    COFD[14] = -2.50852412E-01;
    COFD[15] = 1.07992653E-02;
    COFD[16] = -1.29727175E+01;
    COFD[17] = 3.28440910E+00;
    COFD[18] = -2.11707534E-01;
    COFD[19] = 9.16419558E-03;
    COFD[20] = -1.91952961E+01;
    COFD[21] = 5.44274135E+00;
    COFD[22] = -4.65397833E-01;
    COFD[23] = 1.91171363E-02;
    COFD[24] = -1.38467617E+01;
    COFD[25] = 3.47073417E+00;
    COFD[26] = -2.34706158E-01;
    COFD[27] = 1.01108517E-02;
    COFD[28] = -1.69855694E+01;
    COFD[29] = 4.58049527E+00;
    COFD[30] = -3.70921885E-01;
    COFD[31] = 1.56999324E-02;
    COFD[32] = -1.51283685E+01;
    COFD[33] = 3.92295205E+00;
    COFD[34] = -2.90942506E-01;
    COFD[35] = 1.24465450E-02;
    COFD[36] = -1.50643750E+01;
    COFD[37] = 3.90641370E+00;
    COFD[38] = -2.88957677E-01;
    COFD[39] = 1.23669987E-02;
    COFD[40] = -1.41188260E+01;
    COFD[41] = 3.59720720E+00;
    COFD[42] = -2.50852412E-01;
    COFD[43] = 1.07992653E-02;
    COFD[44] = -1.41192653E+01;
    COFD[45] = 3.59720720E+00;
    COFD[46] = -2.50852412E-01;
    COFD[47] = 1.07992653E-02;
    COFD[48] = -1.93172760E+01;
    COFD[49] = 5.37588635E+00;
    COFD[50] = -4.60762834E-01;
    COFD[51] = 1.90722493E-02;
    COFD[52] = -1.72388534E+01;
    COFD[53] = 4.60447813E+00;
    COFD[54] = -3.73889368E-01;
    COFD[55] = 1.58229235E-02;
    COFD[56] = -1.68005136E+01;
    COFD[57] = 4.55072147E+00;
    COFD[58] = -3.67343400E-01;
    COFD[59] = 1.55560012E-02;
    COFD[60] = -1.72382912E+01;
    COFD[61] = 4.60447813E+00;
    COFD[62] = -3.73889368E-01;
    COFD[63] = 1.58229235E-02;
    COFD[64] = -1.73275531E+01;
    COFD[65] = 4.71167685E+00;
    COFD[66] = -3.87004255E-01;
    COFD[67] = 1.63602477E-02;
    COFD[68] = -1.73283030E+01;
    COFD[69] = 4.71167685E+00;
    COFD[70] = -3.87004255E-01;
    COFD[71] = 1.63602477E-02;
    COFD[72] = -1.92320503E+01;
    COFD[73] = 5.34516253E+00;
    COFD[74] = -4.57548788E-01;
    COFD[75] = 1.89625374E-02;
    COFD[76] = -1.89960613E+01;
    COFD[77] = 5.22783910E+00;
    COFD[78] = -4.44450099E-01;
    COFD[79] = 1.84764493E-02;
    COFD[80] = -1.47582497E+01;
    COFD[81] = 3.96253155E+00;
    COFD[82] = -2.95745106E-01;
    COFD[83] = 1.26413618E-02;
    COFD[84] = -1.89965974E+01;
    COFD[85] = 5.22783910E+00;
    COFD[86] = -4.44450099E-01;
    COFD[87] = 1.84764493E-02;
    COFD[88] = -1.82115937E+01;
    COFD[89] = 4.94453558E+00;
    COFD[90] = -4.13840584E-01;
    COFD[91] = 1.73898654E-02;
    COFD[92] = -1.82115937E+01;
    COFD[93] = 4.94453558E+00;
    COFD[94] = -4.13840584E-01;
    COFD[95] = 1.73898654E-02;
    COFD[96] = -1.79863237E+01;
    COFD[97] = 4.88559216E+00;
    COFD[98] = -4.07208104E-01;
    COFD[99] = 1.71425982E-02;
    COFD[100] = -1.86950023E+01;
    COFD[101] = 5.06329467E+00;
    COFD[102] = -4.27615122E-01;
    COFD[103] = 1.79228042E-02;
    COFD[104] = -1.86585461E+01;
    COFD[105] = 5.05830652E+00;
    COFD[106] = -4.27083712E-01;
    COFD[107] = 1.79042712E-02;
    COFD[108] = -1.86587075E+01;
    COFD[109] = 5.05830652E+00;
    COFD[110] = -4.27083712E-01;
    COFD[111] = 1.79042712E-02;
    COFD[112] = -1.88039398E+01;
    COFD[113] = 5.05088527E+00;
    COFD[114] = -4.26300667E-01;
    COFD[115] = 1.78773589E-02;
    COFD[116] = -1.92706041E+01;
    COFD[117] = 5.22641306E+00;
    COFD[118] = -4.44290196E-01;
    COFD[119] = 1.84704961E-02;
    COFD[120] = -1.92100304E+01;
    COFD[121] = 5.22556324E+00;
    COFD[122] = -4.44195247E-01;
    COFD[123] = 1.84669783E-02;
    COFD[124] = -1.92110979E+01;
    COFD[125] = 5.16046404E+00;
    COFD[126] = -4.38008342E-01;
    COFD[127] = 1.82870041E-02;
    COFD[128] = -1.91375595E+01;
    COFD[129] = 5.13436741E+00;
    COFD[130] = -4.35193261E-01;
    COFD[131] = 1.81872196E-02;
    COFD[132] = -1.93697973E+01;
    COFD[133] = 5.24163725E+00;
    COFD[134] = -4.46011560E-01;
    COFD[135] = 1.85352610E-02;
    COFD[136] = -1.92369345E+01;
    COFD[137] = 5.34516253E+00;
    COFD[138] = -4.57548788E-01;
    COFD[139] = 1.89625374E-02;
    COFD[140] = -1.92371546E+01;
    COFD[141] = 5.34516253E+00;
    COFD[142] = -4.57548788E-01;
    COFD[143] = 1.89625374E-02;
    COFD[144] = -1.91106350E+01;
    COFD[145] = 5.24865128E+00;
    COFD[146] = -4.46802764E-01;
    COFD[147] = 1.85649476E-02;
    COFD[148] = -1.92438114E+01;
    COFD[149] = 5.20849012E+00;
    COFD[150] = -4.42302201E-01;
    COFD[151] = 1.83974704E-02;
    COFD[152] = -1.80817998E+01;
    COFD[153] = 4.90738815E+00;
    COFD[154] = -4.09483894E-01;
    COFD[155] = 1.72192910E-02;
    COFD[156] = -1.82112771E+01;
    COFD[157] = 4.94453558E+00;
    COFD[158] = -4.13840584E-01;
    COFD[159] = 1.73898654E-02;
    COFD[160] = -1.95514332E+01;
    COFD[161] = 5.30461369E+00;
    COFD[162] = -4.53088496E-01;
    COFD[163] = 1.87998409E-02;
    COFD[164] = -1.95829229E+01;
    COFD[165] = 5.27835392E+00;
    COFD[166] = -4.50143769E-01;
    COFD[167] = 1.86899020E-02;
    COFD[168] = -1.91107966E+01;
    COFD[169] = 5.24865128E+00;
    COFD[170] = -4.46802764E-01;
    COFD[171] = 1.85649476E-02;
    COFD[172] = -1.98427244E+01;
    COFD[173] = 5.33546471E+00;
    COFD[174] = -4.56520851E-01;
    COFD[175] = 1.89268199E-02;
    COFD[176] = -1.98161406E+01;
    COFD[177] = 5.35897956E+00;
    COFD[178] = -4.58992472E-01;
    COFD[179] = 1.90117416E-02;
    COFD[180] = -1.98430923E+01;
    COFD[181] = 5.29318575E+00;
    COFD[182] = -4.51815332E-01;
    COFD[183] = 1.87526763E-02;
    COFD[184] = -1.98430417E+01;
    COFD[185] = 5.29318575E+00;
    COFD[186] = -4.51815332E-01;
    COFD[187] = 1.87526763E-02;
    COFD[188] = -2.03662784E+01;
    COFD[189] = 5.43917732E+00;
    COFD[190] = -4.65585169E-01;
    COFD[191] = 1.91504194E-02;
    COFD[192] = -2.06638821E+01;
    COFD[193] = 5.45347668E+00;
    COFD[194] = -4.65289819E-01;
    COFD[195] = 1.90547288E-02;
    COFD[196] = -2.01926972E+01;
    COFD[197] = 5.40010667E+00;
    COFD[198] = -4.63239430E-01;
    COFD[199] = 1.91539217E-02;
    COFD[200] = -2.04559021E+01;
    COFD[201] = 5.44558859E+00;
    COFD[202] = -4.65276039E-01;
    COFD[203] = 1.90929031E-02;
    COFD[204] = -1.38221692E+01;
    COFD[205] = 3.46560000E+00;
    COFD[206] = -2.34075114E-01;
    COFD[207] = 1.00849661E-02;
    COFD[208] = -1.09419258E+01;
    COFD[209] = 2.56365104E+00;
    COFD[210] = -1.18913262E-01;
    COFD[211] = 5.17669825E-03;
    COFD[212] = -9.78881019E+00;
    COFD[213] = 1.95527707E+00;
    COFD[214] = -4.07233312E-02;
    COFD[215] = 1.83285363E-03;
    COFD[216] = -1.06253201E+01;
    COFD[217] = 2.15819366E+00;
    COFD[218] = -6.52830580E-02;
    COFD[219] = 2.80657258E-03;
    COFD[220] = -1.12814216E+01;
    COFD[221] = 2.30017030E+00;
    COFD[222] = -8.34928497E-02;
    COFD[223] = 3.58611589E-03;
    COFD[224] = -1.06286471E+01;
    COFD[225] = 2.15819366E+00;
    COFD[226] = -6.52830580E-02;
    COFD[227] = 2.80657258E-03;
    COFD[228] = -1.65320223E+01;
    COFD[229] = 4.34135676E+00;
    COFD[230] = -3.42166619E-01;
    COFD[231] = 1.45425431E-02;
    COFD[232] = -1.11288595E+01;
    COFD[233] = 2.21964693E+00;
    COFD[234] = -7.27132181E-02;
    COFD[235] = 3.10467640E-03;
    COFD[236] = -1.33494365E+01;
    COFD[237] = 3.04769365E+00;
    COFD[238] = -1.81586731E-01;
    COFD[239] = 7.88540063E-03;
    COFD[240] = -1.20163448E+01;
    COFD[241] = 2.55772522E+00;
    COFD[242] = -1.18119456E-01;
    COFD[243] = 5.14119183E-03;
    COFD[244] = -1.19609333E+01;
    COFD[245] = 2.54091639E+00;
    COFD[246] = -1.15850475E-01;
    COFD[247] = 5.03893346E-03;
    COFD[248] = -1.12823274E+01;
    COFD[249] = 2.30017030E+00;
    COFD[250] = -8.34928497E-02;
    COFD[251] = 3.58611589E-03;
    COFD[252] = -1.12831810E+01;
    COFD[253] = 2.30017030E+00;
    COFD[254] = -8.34928497E-02;
    COFD[255] = 3.58611589E-03;
    COFD[256] = -1.54361025E+01;
    COFD[257] = 3.83529069E+00;
    COFD[258] = -2.80481778E-01;
    COFD[259] = 1.20302770E-02;
    COFD[260] = -1.35410895E+01;
    COFD[261] = 3.05984388E+00;
    COFD[262] = -1.83114865E-01;
    COFD[263] = 7.94950068E-03;
    COFD[264] = -1.32011826E+01;
    COFD[265] = 3.02784198E+00;
    COFD[266] = -1.79087955E-01;
    COFD[267] = 7.78048913E-03;
    COFD[268] = -1.35400012E+01;
    COFD[269] = 3.05984388E+00;
    COFD[270] = -1.83114865E-01;
    COFD[271] = 7.94950068E-03;
    COFD[272] = -1.35202960E+01;
    COFD[273] = 3.12061076E+00;
    COFD[274] = -1.90761802E-01;
    COFD[275] = 8.27051223E-03;
    COFD[276] = -1.35217408E+01;
    COFD[277] = 3.12061076E+00;
    COFD[278] = -1.90761802E-01;
    COFD[279] = 8.27051223E-03;
    COFD[280] = -1.53596434E+01;
    COFD[281] = 3.80403213E+00;
    COFD[282] = -2.76736813E-01;
    COFD[283] = 1.18808568E-02;
    COFD[284] = -1.51423689E+01;
    COFD[285] = 3.68565858E+00;
    COFD[286] = -2.62077329E-01;
    COFD[287] = 1.12749663E-02;
    COFD[288] = -1.17243965E+01;
    COFD[289] = 2.59606889E+00;
    COFD[290] = -1.23237035E-01;
    COFD[291] = 5.36924411E-03;
    COFD[292] = -1.51434170E+01;
    COFD[293] = 3.68565858E+00;
    COFD[294] = -2.62077329E-01;
    COFD[295] = 1.12749663E-02;
    COFD[296] = -1.43067326E+01;
    COFD[297] = 3.34609996E+00;
    COFD[298] = -2.19417637E-01;
    COFD[299] = 9.48565577E-03;
    COFD[300] = -1.43067326E+01;
    COFD[301] = 3.34609996E+00;
    COFD[302] = -2.19417637E-01;
    COFD[303] = 9.48565577E-03;
    COFD[304] = -1.41186747E+01;
    COFD[305] = 3.29251025E+00;
    COFD[306] = -2.12716771E-01;
    COFD[307] = 9.20613624E-03;
    COFD[308] = -1.46609678E+01;
    COFD[309] = 3.42877141E+00;
    COFD[310] = -2.29550711E-01;
    COFD[311] = 9.89957140E-03;
    COFD[312] = -1.46294931E+01;
    COFD[313] = 3.42403711E+00;
    COFD[314] = -2.28974532E-01;
    COFD[315] = 9.87621386E-03;
    COFD[316] = -1.46298102E+01;
    COFD[317] = 3.42403711E+00;
    COFD[318] = -2.28974532E-01;
    COFD[319] = 9.87621386E-03;
    COFD[320] = -1.47601286E+01;
    COFD[321] = 3.41694605E+00;
    COFD[322] = -2.28112628E-01;
    COFD[323] = 9.84132390E-03;
    COFD[324] = -1.53871573E+01;
    COFD[325] = 3.68408734E+00;
    COFD[326] = -2.61879015E-01;
    COFD[327] = 1.12666094E-02;
    COFD[328] = -1.57444919E+01;
    COFD[329] = 3.82742326E+00;
    COFD[330] = -2.79536566E-01;
    COFD[331] = 1.19924367E-02;
    COFD[332] = -1.51818187E+01;
    COFD[333] = 3.55111766E+00;
    COFD[334] = -2.44966429E-01;
    COFD[335] = 1.05482272E-02;
    COFD[336] = -1.50908350E+01;
    COFD[337] = 3.51404074E+00;
    COFD[338] = -2.40221542E-01;
    COFD[339] = 1.03454651E-02;
    COFD[340] = -1.54795771E+01;
    COFD[341] = 3.69992299E+00;
    COFD[342] = -2.63867460E-01;
    COFD[343] = 1.13499698E-02;
    COFD[344] = -1.53691667E+01;
    COFD[345] = 3.80403213E+00;
    COFD[346] = -2.76736813E-01;
    COFD[347] = 1.18808568E-02;
    COFD[348] = -1.53695980E+01;
    COFD[349] = 3.80403213E+00;
    COFD[350] = -2.76736813E-01;
    COFD[351] = 1.18808568E-02;
    COFD[352] = -1.52493207E+01;
    COFD[353] = 3.70713307E+00;
    COFD[354] = -2.64771252E-01;
    COFD[355] = 1.13878078E-02;
    COFD[356] = -1.57696917E+01;
    COFD[357] = 3.80988496E+00;
    COFD[358] = -2.77445338E-01;
    COFD[359] = 1.19094487E-02;
    COFD[360] = -1.42128070E+01;
    COFD[361] = 3.31879365E+00;
    COFD[362] = -2.16004062E-01;
    COFD[363] = 9.34330595E-03;
    COFD[364] = -1.43061148E+01;
    COFD[365] = 3.34609996E+00;
    COFD[366] = -2.19417637E-01;
    COFD[367] = 9.48565577E-03;
    COFD[368] = -1.59135880E+01;
    COFD[369] = 3.85500637E+00;
    COFD[370] = -2.82834219E-01;
    COFD[371] = 1.21237924E-02;
    COFD[372] = -1.60254653E+01;
    COFD[373] = 3.85862177E+00;
    COFD[374] = -2.83262563E-01;
    COFD[375] = 1.21406839E-02;
    COFD[376] = -1.52496383E+01;
    COFD[377] = 3.70713307E+00;
    COFD[378] = -2.64771252E-01;
    COFD[379] = 1.13878078E-02;
    COFD[380] = -1.60659819E+01;
    COFD[381] = 3.84737458E+00;
    COFD[382] = -2.81927023E-01;
    COFD[383] = 1.20878691E-02;
    COFD[384] = -1.58862636E+01;
    COFD[385] = 3.81826954E+00;
    COFD[386] = -2.78442735E-01;
    COFD[387] = 1.19489210E-02;
    COFD[388] = -1.59092604E+01;
    COFD[389] = 3.75176869E+00;
    COFD[390] = -2.70323254E-01;
    COFD[391] = 1.16183021E-02;
    COFD[392] = -1.59091602E+01;
    COFD[393] = 3.75176869E+00;
    COFD[394] = -2.70323254E-01;
    COFD[395] = 1.16183021E-02;
    COFD[396] = -1.65211535E+01;
    COFD[397] = 3.96045057E+00;
    COFD[398] = -2.95491083E-01;
    COFD[399] = 1.26309912E-02;
    COFD[400] = -1.69048097E+01;
    COFD[401] = 4.02910047E+00;
    COFD[402] = -3.03801059E-01;
    COFD[403] = 1.29668774E-02;
    COFD[404] = -1.62247094E+01;
    COFD[405] = 3.85939026E+00;
    COFD[406] = -2.83353647E-01;
    COFD[407] = 1.21442775E-02;
    COFD[408] = -1.66789079E+01;
    COFD[409] = 4.00162523E+00;
    COFD[410] = -3.00504266E-01;
    COFD[411] = 1.28349600E-02;
    COFD[412] = -1.11071542E+01;
    COFD[413] = 2.21467976E+00;
    COFD[414] = -7.20501534E-02;
    COFD[415] = 3.07515437E-03;
    COFD[416] = -1.29709581E+01;
    COFD[417] = 3.28440910E+00;
    COFD[418] = -2.11707534E-01;
    COFD[419] = 9.16419558E-03;
    COFD[420] = -1.06253201E+01;
    COFD[421] = 2.15819366E+00;
    COFD[422] = -6.52830580E-02;
    COFD[423] = 2.80657258E-03;
    COFD[424] = -1.27611493E+01;
    COFD[425] = 2.70797975E+00;
    COFD[426] = -1.37895013E-01;
    COFD[427] = 6.01030858E-03;
    COFD[428] = -1.38802699E+01;
    COFD[429] = 2.98298723E+00;
    COFD[430] = -1.73428215E-01;
    COFD[431] = 7.54224684E-03;
    COFD[432] = -1.27761900E+01;
    COFD[433] = 2.70797975E+00;
    COFD[434] = -1.37895013E-01;
    COFD[435] = 6.01030858E-03;
    COFD[436] = -1.88412847E+01;
    COFD[437] = 4.90686422E+00;
    COFD[438] = -4.09428407E-01;
    COFD[439] = 1.72173805E-02;
    COFD[440] = -1.37003145E+01;
    COFD[441] = 2.90243493E+00;
    COFD[442] = -1.63148841E-01;
    COFD[443] = 7.10447920E-03;
    COFD[444] = -1.64305407E+01;
    COFD[445] = 3.86315631E+00;
    COFD[446] = -2.83800980E-01;
    COFD[447] = 1.21619703E-02;
    COFD[448] = -1.46233077E+01;
    COFD[449] = 3.27714313E+00;
    COFD[450] = -2.10795461E-01;
    COFD[451] = 9.12600427E-03;
    COFD[452] = -1.45714634E+01;
    COFD[453] = 3.25814060E+00;
    COFD[454] = -2.08402129E-01;
    COFD[455] = 9.02544298E-03;
    COFD[456] = -1.38853857E+01;
    COFD[457] = 2.98298723E+00;
    COFD[458] = -1.73428215E-01;
    COFD[459] = 7.54224684E-03;
    COFD[460] = -1.38902466E+01;
    COFD[461] = 2.98298723E+00;
    COFD[462] = -1.73428215E-01;
    COFD[463] = 7.54224684E-03;
    COFD[464] = -1.88359692E+01;
    COFD[465] = 4.75778339E+00;
    COFD[466] = -3.92510587E-01;
    COFD[467] = 1.65802295E-02;
    COFD[468] = -1.65809526E+01;
    COFD[469] = 3.87564052E+00;
    COFD[470] = -2.85285750E-01;
    COFD[471] = 1.22208235E-02;
    COFD[472] = -1.62171827E+01;
    COFD[473] = 3.84275248E+00;
    COFD[474] = -2.81379017E-01;
    COFD[475] = 1.20662349E-02;
    COFD[476] = -1.65749660E+01;
    COFD[477] = 3.87564052E+00;
    COFD[478] = -2.85285750E-01;
    COFD[479] = 1.22208235E-02;
    COFD[480] = -1.65414090E+01;
    COFD[481] = 3.93824541E+00;
    COFD[482] = -2.92787590E-01;
    COFD[483] = 1.25209211E-02;
    COFD[484] = -1.65490142E+01;
    COFD[485] = 3.93824541E+00;
    COFD[486] = -2.92787590E-01;
    COFD[487] = 1.25209211E-02;
    COFD[488] = -1.87339672E+01;
    COFD[489] = 4.71430273E+00;
    COFD[490] = -3.87319126E-01;
    COFD[491] = 1.63728857E-02;
    COFD[492] = -1.84519843E+01;
    COFD[493] = 4.56059591E+00;
    COFD[494] = -3.68504274E-01;
    COFD[495] = 1.56015944E-02;
    COFD[496] = -1.44910354E+01;
    COFD[497] = 3.31811351E+00;
    COFD[498] = -2.15918664E-01;
    COFD[499] = 9.33972748E-03;
    COFD[500] = -1.84583320E+01;
    COFD[501] = 4.56059591E+00;
    COFD[502] = -3.68504274E-01;
    COFD[503] = 1.56015944E-02;
    COFD[504] = -1.74524239E+01;
    COFD[505] = 4.17306000E+00;
    COFD[506] = -3.21208070E-01;
    COFD[507] = 1.36692112E-02;
    COFD[508] = -1.74524239E+01;
    COFD[509] = 4.17306000E+00;
    COFD[510] = -3.21208070E-01;
    COFD[511] = 1.36692112E-02;
    COFD[512] = -1.72189894E+01;
    COFD[513] = 4.10247595E+00;
    COFD[514] = -3.12590710E-01;
    COFD[515] = 1.33179236E-02;
    COFD[516] = -1.80017096E+01;
    COFD[517] = 4.31633461E+00;
    COFD[518] = -3.39053638E-01;
    COFD[519] = 1.44128540E-02;
    COFD[520] = -1.79605052E+01;
    COFD[521] = 4.30733766E+00;
    COFD[522] = -3.37925627E-01;
    COFD[523] = 1.43655069E-02;
    COFD[524] = -1.79625305E+01;
    COFD[525] = 4.30733766E+00;
    COFD[526] = -3.37925627E-01;
    COFD[527] = 1.43655069E-02;
    COFD[528] = -1.80806236E+01;
    COFD[529] = 4.29360798E+00;
    COFD[530] = -3.36202347E-01;
    COFD[531] = 1.42930967E-02;
    COFD[532] = -1.87347447E+01;
    COFD[533] = 4.55932865E+00;
    COFD[534] = -3.68355198E-01;
    COFD[535] = 1.55957356E-02;
    COFD[536] = -1.86817105E+01;
    COFD[537] = 4.55856850E+00;
    COFD[538] = -3.68265791E-01;
    COFD[539] = 1.55922225E-02;
    COFD[540] = -1.86097280E+01;
    COFD[541] = 4.45757911E+00;
    COFD[542] = -3.56368733E-01;
    COFD[543] = 1.51238590E-02;
    COFD[544] = -1.85284648E+01;
    COFD[545] = 4.42603366E+00;
    COFD[546] = -3.52586895E-01;
    COFD[547] = 1.49720921E-02;
    COFD[548] = -1.88465938E+01;
    COFD[549] = 4.57429837E+00;
    COFD[550] = -3.70149872E-01;
    COFD[551] = 1.56677324E-02;
    COFD[552] = -1.87900753E+01;
    COFD[553] = 4.71430273E+00;
    COFD[554] = -3.87319126E-01;
    COFD[555] = 1.63728857E-02;
    COFD[556] = -1.87927452E+01;
    COFD[557] = 4.71430273E+00;
    COFD[558] = -3.87319126E-01;
    COFD[559] = 1.63728857E-02;
    COFD[560] = -1.85962044E+01;
    COFD[561] = 4.58430295E+00;
    COFD[562] = -3.71395187E-01;
    COFD[563] = 1.57196293E-02;
    COFD[564] = -1.87089834E+01;
    COFD[565] = 4.54237822E+00;
    COFD[566] = -3.66363712E-01;
    COFD[567] = 1.55175729E-02;
    COFD[568] = -1.73332914E+01;
    COFD[569] = 4.13619942E+00;
    COFD[570] = -3.16693330E-01;
    COFD[571] = 1.34845235E-02;
    COFD[572] = -1.74487550E+01;
    COFD[573] = 4.17306000E+00;
    COFD[574] = -3.21208070E-01;
    COFD[575] = 1.36692112E-02;
    COFD[576] = -1.90756172E+01;
    COFD[577] = 4.65976652E+00;
    COFD[578] = -3.80698349E-01;
    COFD[579] = 1.61038060E-02;
    COFD[580] = -1.90773916E+01;
    COFD[581] = 4.62426562E+00;
    COFD[582] = -3.76333900E-01;
    COFD[583] = 1.59240937E-02;
    COFD[584] = -1.85982326E+01;
    COFD[585] = 4.58430295E+00;
    COFD[586] = -3.71395187E-01;
    COFD[587] = 1.57196293E-02;
    COFD[588] = -1.93826457E+01;
    COFD[589] = 4.70007538E+00;
    COFD[590] = -3.85598099E-01;
    COFD[591] = 1.63031809E-02;
    COFD[592] = -1.93847925E+01;
    COFD[593] = 4.73396422E+00;
    COFD[594] = -3.89665033E-01;
    COFD[595] = 1.64665068E-02;
    COFD[596] = -1.93555562E+01;
    COFD[597] = 4.64408224E+00;
    COFD[598] = -3.78774559E-01;
    COFD[599] = 1.60247845E-02;
    COFD[600] = -1.93548571E+01;
    COFD[601] = 4.64408224E+00;
    COFD[602] = -3.78774559E-01;
    COFD[603] = 1.60247845E-02;
    COFD[604] = -2.01011659E+01;
    COFD[605] = 4.89020393E+00;
    COFD[606] = -4.07690029E-01;
    COFD[607] = 1.71588628E-02;
    COFD[608] = -2.05451641E+01;
    COFD[609] = 4.97155662E+00;
    COFD[610] = -4.17014903E-01;
    COFD[611] = 1.75144222E-02;
    COFD[612] = -1.97866869E+01;
    COFD[613] = 4.78902303E+00;
    COFD[614] = -3.96144974E-01;
    COFD[615] = 1.67211738E-02;
    COFD[616] = -2.02611614E+01;
    COFD[617] = 4.92728595E+00;
    COFD[618] = -4.11813748E-01;
    COFD[619] = 1.73103280E-02;
    COFD[620] = -1.36761960E+01;
    COFD[621] = 2.89678405E+00;
    COFD[622] = -1.62418864E-01;
    COFD[623] = 7.07300755E-03;
    COFD[624] = -1.41183595E+01;
    COFD[625] = 3.59720720E+00;
    COFD[626] = -2.50852412E-01;
    COFD[627] = 1.07992653E-02;
    COFD[628] = -1.12814216E+01;
    COFD[629] = 2.30017030E+00;
    COFD[630] = -8.34928497E-02;
    COFD[631] = 3.58611589E-03;
    COFD[632] = -1.38802699E+01;
    COFD[633] = 2.98298723E+00;
    COFD[634] = -1.73428215E-01;
    COFD[635] = 7.54224684E-03;
    COFD[636] = -1.50981868E+01;
    COFD[637] = 3.27853241E+00;
    COFD[638] = -2.10970532E-01;
    COFD[639] = 9.13336403E-03;
    COFD[640] = -1.39004264E+01;
    COFD[641] = 2.98298723E+00;
    COFD[642] = -1.73428215E-01;
    COFD[643] = 7.54224684E-03;
    COFD[644] = -2.08507625E+01;
    COFD[645] = 5.41592504E+00;
    COFD[646] = -4.64640986E-01;
    COFD[647] = 1.91895970E-02;
    COFD[648] = -1.48406405E+01;
    COFD[649] = 3.17494843E+00;
    COFD[650] = -1.97752410E-01;
    COFD[651] = 8.57054106E-03;
    COFD[652] = -1.77768991E+01;
    COFD[653] = 4.18471443E+00;
    COFD[654] = -3.22640353E-01;
    COFD[655] = 1.37280239E-02;
    COFD[656] = -1.58176414E+01;
    COFD[657] = 3.58720055E+00;
    COFD[658] = -2.49572308E-01;
    COFD[659] = 1.07445741E-02;
    COFD[660] = -1.57555897E+01;
    COFD[661] = 3.56140952E+00;
    COFD[662] = -2.46282892E-01;
    COFD[663] = 1.06044757E-02;
    COFD[664] = -1.51058803E+01;
    COFD[665] = 3.27853241E+00;
    COFD[666] = -2.10970532E-01;
    COFD[667] = 9.13336403E-03;
    COFD[668] = -1.51132274E+01;
    COFD[669] = 3.27853241E+00;
    COFD[670] = -2.10970532E-01;
    COFD[671] = 9.13336403E-03;
    COFD[672] = -2.02419373E+01;
    COFD[673] = 5.09185655E+00;
    COFD[674] = -4.30641631E-01;
    COFD[675] = 1.80275784E-02;
    COFD[676] = -1.79002510E+01;
    COFD[677] = 4.20508193E+00;
    COFD[678] = -3.25158699E-01;
    COFD[679] = 1.38321131E-02;
    COFD[680] = -1.75072049E+01;
    COFD[681] = 4.15347096E+00;
    COFD[682] = -3.18813734E-01;
    COFD[683] = 1.35714974E-02;
    COFD[684] = -1.78913897E+01;
    COFD[685] = 4.20508193E+00;
    COFD[686] = -3.25158699E-01;
    COFD[687] = 1.38321131E-02;
    COFD[688] = -1.79764577E+01;
    COFD[689] = 4.31259444E+00;
    COFD[690] = -3.38584827E-01;
    COFD[691] = 1.43931810E-02;
    COFD[692] = -1.79874377E+01;
    COFD[693] = 4.31259444E+00;
    COFD[694] = -3.38584827E-01;
    COFD[695] = 1.43931810E-02;
    COFD[696] = -2.01764857E+01;
    COFD[697] = 5.06295779E+00;
    COFD[698] = -4.27579258E-01;
    COFD[699] = 1.79215546E-02;
    COFD[700] = -1.99640254E+01;
    COFD[701] = 4.93673368E+00;
    COFD[702] = -4.12924135E-01;
    COFD[703] = 1.73539126E-02;
    COFD[704] = -1.58631082E+01;
    COFD[705] = 3.64628858E+00;
    COFD[706] = -2.57092933E-01;
    COFD[707] = 1.10642427E-02;
    COFD[708] = -1.99740136E+01;
    COFD[709] = 4.93673368E+00;
    COFD[710] = -4.12924135E-01;
    COFD[711] = 1.73539126E-02;
    COFD[712] = -1.89513110E+01;
    COFD[713] = 4.56070745E+00;
    COFD[714] = -3.68517396E-01;
    COFD[715] = 1.56021102E-02;
    COFD[716] = -1.89513110E+01;
    COFD[717] = 4.56070745E+00;
    COFD[718] = -3.68517396E-01;
    COFD[719] = 1.56021102E-02;
    COFD[720] = -1.87553343E+01;
    COFD[721] = 4.50343635E+00;
    COFD[722] = -3.61818471E-01;
    COFD[723] = 1.53404608E-02;
    COFD[724] = -1.95280746E+01;
    COFD[725] = 4.70710890E+00;
    COFD[726] = -3.86453789E-01;
    COFD[727] = 1.63380393E-02;
    COFD[728] = -1.94909629E+01;
    COFD[729] = 4.69879097E+00;
    COFD[730] = -3.85441745E-01;
    COFD[731] = 1.62968073E-02;
    COFD[732] = -1.94942738E+01;
    COFD[733] = 4.69879097E+00;
    COFD[734] = -3.85441745E-01;
    COFD[735] = 1.62968073E-02;
    COFD[736] = -1.96030292E+01;
    COFD[737] = 4.68635682E+00;
    COFD[738] = -3.83931845E-01;
    COFD[739] = 1.62354222E-02;
    COFD[740] = -2.02462660E+01;
    COFD[741] = 4.93517265E+00;
    COFD[742] = -4.12740828E-01;
    COFD[743] = 1.73467245E-02;
    COFD[744] = -2.05570536E+01;
    COFD[745] = 5.05447966E+00;
    COFD[746] = -4.26675417E-01;
    COFD[747] = 1.78900033E-02;
    COFD[748] = -2.01150260E+01;
    COFD[749] = 4.83483781E+00;
    COFD[750] = -4.01422871E-01;
    COFD[751] = 1.69235370E-02;
    COFD[752] = -2.00372657E+01;
    COFD[753] = 4.80531581E+00;
    COFD[754] = -3.98010400E-01;
    COFD[755] = 1.67921710E-02;
    COFD[756] = -2.03732076E+01;
    COFD[757] = 4.95138566E+00;
    COFD[758] = -4.14645044E-01;
    COFD[759] = 1.74214162E-02;
    COFD[760] = -2.02632221E+01;
    COFD[761] = 5.06295779E+00;
    COFD[762] = -4.27579258E-01;
    COFD[763] = 1.79215546E-02;
    COFD[764] = -2.02674877E+01;
    COFD[765] = 5.06295779E+00;
    COFD[766] = -4.27579258E-01;
    COFD[767] = 1.79215546E-02;
    COFD[768] = -2.01222360E+01;
    COFD[769] = 4.95889869E+00;
    COFD[770] = -4.15526934E-01;
    COFD[771] = 1.74559870E-02;
    COFD[772] = -2.05600561E+01;
    COFD[773] = 5.03374066E+00;
    COFD[774] = -4.24314609E-01;
    COFD[775] = 1.78006508E-02;
    COFD[776] = -1.88511274E+01;
    COFD[777] = 4.53065425E+00;
    COFD[778] = -3.64983473E-01;
    COFD[779] = 1.54632547E-02;
    COFD[780] = -1.89456139E+01;
    COFD[781] = 4.56070745E+00;
    COFD[782] = -3.68517396E-01;
    COFD[783] = 1.56021102E-02;
    COFD[784] = -2.07843393E+01;
    COFD[785] = 5.08900550E+00;
    COFD[786] = -4.30341928E-01;
    COFD[787] = 1.80173213E-02;
    COFD[788] = -2.08687893E+01;
    COFD[789] = 5.08643773E+00;
    COFD[790] = -4.30070504E-01;
    COFD[791] = 1.80079562E-02;
    COFD[792] = -2.01255512E+01;
    COFD[793] = 4.95889869E+00;
    COFD[794] = -4.15526934E-01;
    COFD[795] = 1.74559870E-02;
    COFD[796] = -2.09753274E+01;
    COFD[797] = 5.08993111E+00;
    COFD[798] = -4.30439721E-01;
    COFD[799] = 1.80206932E-02;
    COFD[800] = -2.08403997E+01;
    COFD[801] = 5.07617318E+00;
    COFD[802] = -4.28983517E-01;
    COFD[803] = 1.79703565E-02;
    COFD[804] = -2.08510201E+01;
    COFD[805] = 5.00799007E+00;
    COFD[806] = -4.21303846E-01;
    COFD[807] = 1.76831398E-02;
    COFD[808] = -2.08497919E+01;
    COFD[809] = 5.00799007E+00;
    COFD[810] = -4.21303846E-01;
    COFD[811] = 1.76831398E-02;
    COFD[812] = -2.15001622E+01;
    COFD[813] = 5.19791528E+00;
    COFD[814] = -4.41398852E-01;
    COFD[815] = 1.83766038E-02;
    COFD[816] = -2.19077574E+01;
    COFD[817] = 5.26494284E+00;
    COFD[818] = -4.48633355E-01;
    COFD[819] = 1.86333022E-02;
    COFD[820] = -2.12048425E+01;
    COFD[821] = 5.11479735E+00;
    COFD[822] = -4.33050578E-01;
    COFD[823] = 1.81098360E-02;
    COFD[824] = -2.16612294E+01;
    COFD[825] = 5.23259254E+00;
    COFD[826] = -4.44988191E-01;
    COFD[827] = 1.84967231E-02;
    COFD[828] = -1.48133672E+01;
    COFD[829] = 3.16773813E+00;
    COFD[830] = -1.96823168E-01;
    COFD[831] = 8.53058072E-03;
    COFD[832] = -1.29727175E+01;
    COFD[833] = 3.28440910E+00;
    COFD[834] = -2.11707534E-01;
    COFD[835] = 9.16419558E-03;
    COFD[836] = -1.06286471E+01;
    COFD[837] = 2.15819366E+00;
    COFD[838] = -6.52830580E-02;
    COFD[839] = 2.80657258E-03;
    COFD[840] = -1.27761900E+01;
    COFD[841] = 2.70797975E+00;
    COFD[842] = -1.37895013E-01;
    COFD[843] = 6.01030858E-03;
    COFD[844] = -1.39004264E+01;
    COFD[845] = 2.98298723E+00;
    COFD[846] = -1.73428215E-01;
    COFD[847] = 7.54224684E-03;
    COFD[848] = -1.27916971E+01;
    COFD[849] = 2.70797975E+00;
    COFD[850] = -1.37895013E-01;
    COFD[851] = 6.01030858E-03;
    COFD[852] = -1.88572311E+01;
    COFD[853] = 4.90686422E+00;
    COFD[854] = -4.09428407E-01;
    COFD[855] = 1.72173805E-02;
    COFD[856] = -1.37195399E+01;
    COFD[857] = 2.90243493E+00;
    COFD[858] = -1.63148841E-01;
    COFD[859] = 7.10447920E-03;
    COFD[860] = -1.64527598E+01;
    COFD[861] = 3.86315631E+00;
    COFD[862] = -2.83800980E-01;
    COFD[863] = 1.21619703E-02;
    COFD[864] = -1.46378741E+01;
    COFD[865] = 3.27714313E+00;
    COFD[866] = -2.10795461E-01;
    COFD[867] = 9.12600427E-03;
    COFD[868] = -1.45865248E+01;
    COFD[869] = 3.25814060E+00;
    COFD[870] = -2.08402129E-01;
    COFD[871] = 9.02544298E-03;
    COFD[872] = -1.39057537E+01;
    COFD[873] = 2.98298723E+00;
    COFD[874] = -1.73428215E-01;
    COFD[875] = 7.54224684E-03;
    COFD[876] = -1.39108177E+01;
    COFD[877] = 2.98298723E+00;
    COFD[878] = -1.73428215E-01;
    COFD[879] = 7.54224684E-03;
    COFD[880] = -1.88556852E+01;
    COFD[881] = 4.75778339E+00;
    COFD[882] = -3.92510587E-01;
    COFD[883] = 1.65802295E-02;
    COFD[884] = -1.66006787E+01;
    COFD[885] = 3.87564052E+00;
    COFD[886] = -2.85285750E-01;
    COFD[887] = 1.22208235E-02;
    COFD[888] = -1.62364191E+01;
    COFD[889] = 3.84275248E+00;
    COFD[890] = -2.81379017E-01;
    COFD[891] = 1.20662349E-02;
    COFD[892] = -1.65944527E+01;
    COFD[893] = 3.87564052E+00;
    COFD[894] = -2.85285750E-01;
    COFD[895] = 1.22208235E-02;
    COFD[896] = -1.65598218E+01;
    COFD[897] = 3.93824541E+00;
    COFD[898] = -2.92787590E-01;
    COFD[899] = 1.25209211E-02;
    COFD[900] = -1.65677145E+01;
    COFD[901] = 3.93824541E+00;
    COFD[902] = -2.92787590E-01;
    COFD[903] = 1.25209211E-02;
    COFD[904] = -1.87541329E+01;
    COFD[905] = 4.71430273E+00;
    COFD[906] = -3.87319126E-01;
    COFD[907] = 1.63728857E-02;
    COFD[908] = -1.84739227E+01;
    COFD[909] = 4.56059591E+00;
    COFD[910] = -3.68504274E-01;
    COFD[911] = 1.56015944E-02;
    COFD[912] = -1.45128230E+01;
    COFD[913] = 3.31811351E+00;
    COFD[914] = -2.15918664E-01;
    COFD[915] = 9.33972748E-03;
    COFD[916] = -1.84805571E+01;
    COFD[917] = 4.56059591E+00;
    COFD[918] = -3.68504274E-01;
    COFD[919] = 1.56015944E-02;
    COFD[920] = -1.74740622E+01;
    COFD[921] = 4.17306000E+00;
    COFD[922] = -3.21208070E-01;
    COFD[923] = 1.36692112E-02;
    COFD[924] = -1.74740622E+01;
    COFD[925] = 4.17306000E+00;
    COFD[926] = -3.21208070E-01;
    COFD[927] = 1.36692112E-02;
    COFD[928] = -1.72409342E+01;
    COFD[929] = 4.10247595E+00;
    COFD[930] = -3.12590710E-01;
    COFD[931] = 1.33179236E-02;
    COFD[932] = -1.80251183E+01;
    COFD[933] = 4.31633461E+00;
    COFD[934] = -3.39053638E-01;
    COFD[935] = 1.44128540E-02;
    COFD[936] = -1.79840143E+01;
    COFD[937] = 4.30733766E+00;
    COFD[938] = -3.37925627E-01;
    COFD[939] = 1.43655069E-02;
    COFD[940] = -1.79861374E+01;
    COFD[941] = 4.30733766E+00;
    COFD[942] = -3.37925627E-01;
    COFD[943] = 1.43655069E-02;
    COFD[944] = -1.81043255E+01;
    COFD[945] = 4.29360798E+00;
    COFD[946] = -3.36202347E-01;
    COFD[947] = 1.42930967E-02;
    COFD[948] = -1.87584426E+01;
    COFD[949] = 4.55932865E+00;
    COFD[950] = -3.68355198E-01;
    COFD[951] = 1.55957356E-02;
    COFD[952] = -1.87055009E+01;
    COFD[953] = 4.55856850E+00;
    COFD[954] = -3.68265791E-01;
    COFD[955] = 1.55922225E-02;
    COFD[956] = -1.86343902E+01;
    COFD[957] = 4.45757911E+00;
    COFD[958] = -3.56368733E-01;
    COFD[959] = 1.51238590E-02;
    COFD[960] = -1.85531954E+01;
    COFD[961] = 4.42603366E+00;
    COFD[962] = -3.52586895E-01;
    COFD[963] = 1.49720921E-02;
    COFD[964] = -1.88713913E+01;
    COFD[965] = 4.57429837E+00;
    COFD[966] = -3.70149872E-01;
    COFD[967] = 1.56677324E-02;
    COFD[968] = -1.88126908E+01;
    COFD[969] = 4.71430273E+00;
    COFD[970] = -3.87319126E-01;
    COFD[971] = 1.63728857E-02;
    COFD[972] = -1.88154846E+01;
    COFD[973] = 4.71430273E+00;
    COFD[974] = -3.87319126E-01;
    COFD[975] = 1.63728857E-02;
    COFD[976] = -1.86197092E+01;
    COFD[977] = 4.58430295E+00;
    COFD[978] = -3.71395187E-01;
    COFD[979] = 1.57196293E-02;
    COFD[980] = -1.87325861E+01;
    COFD[981] = 4.54237822E+00;
    COFD[982] = -3.66363712E-01;
    COFD[983] = 1.55175729E-02;
    COFD[984] = -1.73550857E+01;
    COFD[985] = 4.13619942E+00;
    COFD[986] = -3.16693330E-01;
    COFD[987] = 1.34845235E-02;
    COFD[988] = -1.74702316E+01;
    COFD[989] = 4.17306000E+00;
    COFD[990] = -3.21208070E-01;
    COFD[991] = 1.36692112E-02;
    COFD[992] = -1.91004773E+01;
    COFD[993] = 4.65976652E+00;
    COFD[994] = -3.80698349E-01;
    COFD[995] = 1.61038060E-02;
    COFD[996] = -1.91021863E+01;
    COFD[997] = 4.62426562E+00;
    COFD[998] = -3.76333900E-01;
    COFD[999] = 1.59240937E-02;
    COFD[1000] = -1.86218352E+01;
    COFD[1001] = 4.58430295E+00;
    COFD[1002] = -3.71395187E-01;
    COFD[1003] = 1.57196293E-02;
    COFD[1004] = -1.94082832E+01;
    COFD[1005] = 4.70007538E+00;
    COFD[1006] = -3.85598099E-01;
    COFD[1007] = 1.63031809E-02;
    COFD[1008] = -1.94103813E+01;
    COFD[1009] = 4.73396422E+00;
    COFD[1010] = -3.89665033E-01;
    COFD[1011] = 1.64665068E-02;
    COFD[1012] = -1.93817856E+01;
    COFD[1013] = 4.64408224E+00;
    COFD[1014] = -3.78774559E-01;
    COFD[1015] = 1.60247845E-02;
    COFD[1016] = -1.93810490E+01;
    COFD[1017] = 4.64408224E+00;
    COFD[1018] = -3.78774559E-01;
    COFD[1019] = 1.60247845E-02;
    COFD[1020] = -2.01283014E+01;
    COFD[1021] = 4.89020393E+00;
    COFD[1022] = -4.07690029E-01;
    COFD[1023] = 1.71588628E-02;
    COFD[1024] = -2.05729073E+01;
    COFD[1025] = 4.97155662E+00;
    COFD[1026] = -4.17014903E-01;
    COFD[1027] = 1.75144222E-02;
    COFD[1028] = -1.98133783E+01;
    COFD[1029] = 4.78902303E+00;
    COFD[1030] = -3.96144974E-01;
    COFD[1031] = 1.67211738E-02;
    COFD[1032] = -2.02886114E+01;
    COFD[1033] = 4.92728595E+00;
    COFD[1034] = -4.11813748E-01;
    COFD[1035] = 1.73103280E-02;
    COFD[1036] = -1.36954221E+01;
    COFD[1037] = 2.89678405E+00;
    COFD[1038] = -1.62418864E-01;
    COFD[1039] = 7.07300755E-03;
    COFD[1040] = -1.91952961E+01;
    COFD[1041] = 5.44274135E+00;
    COFD[1042] = -4.65397833E-01;
    COFD[1043] = 1.91171363E-02;
    COFD[1044] = -1.65320223E+01;
    COFD[1045] = 4.34135676E+00;
    COFD[1046] = -3.42166619E-01;
    COFD[1047] = 1.45425431E-02;
    COFD[1048] = -1.88412847E+01;
    COFD[1049] = 4.90686422E+00;
    COFD[1050] = -4.09428407E-01;
    COFD[1051] = 1.72173805E-02;
    COFD[1052] = -2.08507625E+01;
    COFD[1053] = 5.41592504E+00;
    COFD[1054] = -4.64640986E-01;
    COFD[1055] = 1.91895970E-02;
    COFD[1056] = -1.88572311E+01;
    COFD[1057] = 4.90686422E+00;
    COFD[1058] = -4.09428407E-01;
    COFD[1059] = 1.72173805E-02;
    COFD[1060] = -1.42050531E+01;
    COFD[1061] = 1.94854490E+00;
    COFD[1062] = 9.16761580E-02;
    COFD[1063] = -8.62876628E-03;
    COFD[1064] = -2.06859096E+01;
    COFD[1065] = 5.35499363E+00;
    COFD[1066] = -4.58577026E-01;
    COFD[1067] = 1.89976312E-02;
    COFD[1068] = -2.14098987E+01;
    COFD[1069] = 5.30121376E+00;
    COFD[1070] = -4.21263191E-01;
    COFD[1071] = 1.62329486E-02;
    COFD[1072] = -2.08999981E+01;
    COFD[1073] = 5.44151563E+00;
    COFD[1074] = -4.65460498E-01;
    COFD[1075] = 1.91284180E-02;
    COFD[1076] = -2.11989257E+01;
    COFD[1077] = 5.48294376E+00;
    COFD[1078] = -4.63910522E-01;
    COFD[1079] = 1.87967250E-02;
    COFD[1080] = -2.01625839E+01;
    COFD[1081] = 5.20184446E+00;
    COFD[1082] = -4.41727448E-01;
    COFD[1083] = 1.83837722E-02;
    COFD[1084] = -2.01678431E+01;
    COFD[1085] = 5.20184446E+00;
    COFD[1086] = -4.41727448E-01;
    COFD[1087] = 1.83837722E-02;
    COFD[1088] = -1.94329144E+01;
    COFD[1089] = 4.32601172E+00;
    COFD[1090] = -2.67841680E-01;
    COFD[1091] = 8.59623460E-03;
    COFD[1092] = -2.16940491E+01;
    COFD[1093] = 5.44394964E+00;
    COFD[1094] = -4.47152201E-01;
    COFD[1095] = 1.76308159E-02;
    COFD[1096] = -2.14160385E+01;
    COFD[1097] = 5.45640844E+00;
    COFD[1098] = -4.50149956E-01;
    COFD[1099] = 1.78119073E-02;
    COFD[1100] = -2.16875939E+01;
    COFD[1101] = 5.44394964E+00;
    COFD[1102] = -4.47152201E-01;
    COFD[1103] = 1.76308159E-02;
    COFD[1104] = -2.14146533E+01;
    COFD[1105] = 5.39388015E+00;
    COFD[1106] = -4.37574572E-01;
    COFD[1107] = 1.70997728E-02;
    COFD[1108] = -2.14228203E+01;
    COFD[1109] = 5.39388015E+00;
    COFD[1110] = -4.37574572E-01;
    COFD[1111] = 1.70997728E-02;
    COFD[1112] = -1.96338867E+01;
    COFD[1113] = 4.41245751E+00;
    COFD[1114] = -2.80855500E-01;
    COFD[1115] = 9.22505777E-03;
    COFD[1116] = -2.02357367E+01;
    COFD[1117] = 4.63893949E+00;
    COFD[1118] = -3.15589295E-01;
    COFD[1119] = 1.09248423E-02;
    COFD[1120] = -2.07110846E+01;
    COFD[1121] = 5.44973198E+00;
    COFD[1122] = -4.65238176E-01;
    COFD[1123] = 1.90691386E-02;
    COFD[1124] = -2.02426482E+01;
    COFD[1125] = 4.63893949E+00;
    COFD[1126] = -3.15589295E-01;
    COFD[1127] = 1.09248423E-02;
    COFD[1128] = -2.13313704E+01;
    COFD[1129] = 5.17422070E+00;
    COFD[1130] = -3.99997068E-01;
    COFD[1131] = 1.51336871E-02;
    COFD[1132] = -2.13313704E+01;
    COFD[1133] = 5.17422070E+00;
    COFD[1134] = -3.99997068E-01;
    COFD[1135] = 1.51336871E-02;
    COFD[1136] = -2.14155602E+01;
    COFD[1137] = 5.24536650E+00;
    COFD[1138] = -4.11756786E-01;
    COFD[1139] = 1.57371608E-02;
    COFD[1140] = -2.12712683E+01;
    COFD[1141] = 5.04351440E+00;
    COFD[1142] = -3.78432548E-01;
    COFD[1143] = 1.40280477E-02;
    COFD[1144] = -2.12710109E+01;
    COFD[1145] = 5.05274021E+00;
    COFD[1146] = -3.79920248E-01;
    COFD[1147] = 1.41033861E-02;
    COFD[1148] = -2.12732290E+01;
    COFD[1149] = 5.05274021E+00;
    COFD[1150] = -3.79920248E-01;
    COFD[1151] = 1.41033861E-02;
    COFD[1152] = -2.14557574E+01;
    COFD[1153] = 5.06615192E+00;
    COFD[1154] = -3.82089787E-01;
    COFD[1155] = 1.42134191E-02;
    COFD[1156] = -2.05361209E+01;
    COFD[1157] = 4.64159173E+00;
    COFD[1158] = -3.15993927E-01;
    COFD[1159] = 1.09445813E-02;
    COFD[1160] = -1.83638021E+01;
    COFD[1161] = 3.61550901E+00;
    COFD[1162] = -1.62122786E-01;
    COFD[1163] = 3.57292746E-03;
    COFD[1164] = -2.11196901E+01;
    COFD[1165] = 4.84742105E+00;
    COFD[1166] = -3.47604370E-01;
    COFD[1167] = 1.24938668E-02;
    COFD[1168] = -2.10582825E+01;
    COFD[1169] = 4.83010534E+00;
    COFD[1170] = -3.46423002E-01;
    COFD[1171] = 1.24833814E-02;
    COFD[1172] = -2.05625147E+01;
    COFD[1173] = 4.61780386E+00;
    COFD[1174] = -3.12290044E-01;
    COFD[1175] = 1.07615587E-02;
    COFD[1176] = -1.96948070E+01;
    COFD[1177] = 4.41245751E+00;
    COFD[1178] = -2.80855500E-01;
    COFD[1179] = 9.22505777E-03;
    COFD[1180] = -1.96977208E+01;
    COFD[1181] = 4.41245751E+00;
    COFD[1182] = -2.80855500E-01;
    COFD[1183] = 9.22505777E-03;
    COFD[1184] = -2.02528193E+01;
    COFD[1185] = 4.60540477E+00;
    COFD[1186] = -3.10373488E-01;
    COFD[1187] = 1.06672049E-02;
    COFD[1188] = -1.86106672E+01;
    COFD[1189] = 3.70369165E+00;
    COFD[1190] = -1.75747020E-01;
    COFD[1191] = 4.23774985E-03;
    COFD[1192] = -2.13752043E+01;
    COFD[1193] = 5.21010770E+00;
    COFD[1194] = -4.05945421E-01;
    COFD[1195] = 1.54393383E-02;
    COFD[1196] = -2.13273839E+01;
    COFD[1197] = 5.17422070E+00;
    COFD[1198] = -3.99997068E-01;
    COFD[1199] = 1.51336871E-02;
    COFD[1200] = -1.85864727E+01;
    COFD[1201] = 3.68972017E+00;
    COFD[1202] = -1.75261310E-01;
    COFD[1203] = 4.25906382E-03;
    COFD[1204] = -1.86138117E+01;
    COFD[1205] = 3.64451102E+00;
    COFD[1206] = -1.67695039E-01;
    COFD[1207] = 3.87428023E-03;
    COFD[1208] = -2.02550404E+01;
    COFD[1209] = 4.60540477E+00;
    COFD[1210] = -3.10373488E-01;
    COFD[1211] = 1.06672049E-02;
    COFD[1212] = -1.91771822E+01;
    COFD[1213] = 3.89325739E+00;
    COFD[1214] = -2.06939169E-01;
    COFD[1215] = 5.80832388E-03;
    COFD[1216] = -2.01687636E+01;
    COFD[1217] = 4.37329337E+00;
    COFD[1218] = -2.74949941E-01;
    COFD[1219] = 8.93930138E-03;
    COFD[1220] = -2.07027595E+01;
    COFD[1221] = 4.52425469E+00;
    COFD[1222] = -2.97869560E-01;
    COFD[1223] = 1.00531031E-02;
    COFD[1224] = -2.07019859E+01;
    COFD[1225] = 4.52425469E+00;
    COFD[1226] = -2.97869560E-01;
    COFD[1227] = 1.00531031E-02;
    COFD[1228] = -1.90849698E+01;
    COFD[1229] = 3.76037078E+00;
    COFD[1230] = -1.86344485E-01;
    COFD[1231] = 4.78342458E-03;
    COFD[1232] = -1.89362802E+01;
    COFD[1233] = 3.59051852E+00;
    COFD[1234] = -1.61030577E-01;
    COFD[1235] = 3.56695151E-03;
    COFD[1236] = -2.01714866E+01;
    COFD[1237] = 4.25399524E+00;
    COFD[1238] = -2.57013670E-01;
    COFD[1239] = 8.07339307E-03;
    COFD[1240] = -1.88216419E+01;
    COFD[1241] = 3.61672004E+00;
    COFD[1242] = -1.65415347E-01;
    COFD[1243] = 3.79600441E-03;
    COFD[1244] = -2.06077924E+01;
    COFD[1245] = 5.33386529E+00;
    COFD[1246] = -4.56343427E-01;
    COFD[1247] = 1.89202788E-02;
    COFD[1248] = -1.38467617E+01;
    COFD[1249] = 3.47073417E+00;
    COFD[1250] = -2.34706158E-01;
    COFD[1251] = 1.01108517E-02;
    COFD[1252] = -1.11288595E+01;
    COFD[1253] = 2.21964693E+00;
    COFD[1254] = -7.27132181E-02;
    COFD[1255] = 3.10467640E-03;
    COFD[1256] = -1.37003145E+01;
    COFD[1257] = 2.90243493E+00;
    COFD[1258] = -1.63148841E-01;
    COFD[1259] = 7.10447920E-03;
    COFD[1260] = -1.48406405E+01;
    COFD[1261] = 3.17494843E+00;
    COFD[1262] = -1.97752410E-01;
    COFD[1263] = 8.57054106E-03;
    COFD[1264] = -1.37195399E+01;
    COFD[1265] = 2.90243493E+00;
    COFD[1266] = -1.63148841E-01;
    COFD[1267] = 7.10447920E-03;
    COFD[1268] = -2.06859096E+01;
    COFD[1269] = 5.35499363E+00;
    COFD[1270] = -4.58577026E-01;
    COFD[1271] = 1.89976312E-02;
    COFD[1272] = -1.46011865E+01;
    COFD[1273] = 3.07793917E+00;
    COFD[1274] = -1.85367944E-01;
    COFD[1275] = 8.04304684E-03;
    COFD[1276] = -1.74636259E+01;
    COFD[1277] = 4.06842532E+00;
    COFD[1278] = -3.08492475E-01;
    COFD[1279] = 1.31534275E-02;
    COFD[1280] = -1.55192394E+01;
    COFD[1281] = 3.46464286E+00;
    COFD[1282] = -2.33957556E-01;
    COFD[1283] = 1.00801477E-02;
    COFD[1284] = -1.54796758E+01;
    COFD[1285] = 3.44821689E+00;
    COFD[1286] = -2.31932994E-01;
    COFD[1287] = 9.99688061E-03;
    COFD[1288] = -1.48478190E+01;
    COFD[1289] = 3.17494843E+00;
    COFD[1290] = -1.97752410E-01;
    COFD[1291] = 8.57054106E-03;
    COFD[1292] = -1.48546674E+01;
    COFD[1293] = 3.17494843E+00;
    COFD[1294] = -1.97752410E-01;
    COFD[1295] = 8.57054106E-03;
    COFD[1296] = -1.99660141E+01;
    COFD[1297] = 4.99411384E+00;
    COFD[1298] = -4.19672196E-01;
    COFD[1299] = 1.76190456E-02;
    COFD[1300] = -1.75803754E+01;
    COFD[1301] = 4.08497111E+00;
    COFD[1302] = -3.10470821E-01;
    COFD[1303] = 1.32322450E-02;
    COFD[1304] = -1.72146872E+01;
    COFD[1305] = 4.04198813E+00;
    COFD[1306] = -3.05338303E-01;
    COFD[1307] = 1.30280077E-02;
    COFD[1308] = -1.75720813E+01;
    COFD[1309] = 4.08497111E+00;
    COFD[1310] = -3.10470821E-01;
    COFD[1311] = 1.32322450E-02;
    COFD[1312] = -1.76229592E+01;
    COFD[1313] = 4.17526194E+00;
    COFD[1314] = -3.21478247E-01;
    COFD[1315] = 1.36802836E-02;
    COFD[1316] = -1.76332867E+01;
    COFD[1317] = 4.17526194E+00;
    COFD[1318] = -3.21478247E-01;
    COFD[1319] = 1.36802836E-02;
    COFD[1320] = -1.98617110E+01;
    COFD[1321] = 4.94937905E+00;
    COFD[1322] = -4.14409534E-01;
    COFD[1323] = 1.74121859E-02;
    COFD[1324] = -1.96817244E+01;
    COFD[1325] = 4.83908504E+00;
    COFD[1326] = -4.01914941E-01;
    COFD[1327] = 1.69425288E-02;
    COFD[1328] = -1.55380429E+01;
    COFD[1329] = 3.51529689E+00;
    COFD[1330] = -2.40382392E-01;
    COFD[1331] = 1.03523419E-02;
    COFD[1332] = -1.96909586E+01;
    COFD[1333] = 4.83908504E+00;
    COFD[1334] = -4.01914941E-01;
    COFD[1335] = 1.69425288E-02;
    COFD[1336] = -1.86750745E+01;
    COFD[1337] = 4.46203550E+00;
    COFD[1338] = -3.56905086E-01;
    COFD[1339] = 1.51454720E-02;
    COFD[1340] = -1.86750745E+01;
    COFD[1341] = 4.46203550E+00;
    COFD[1342] = -3.56905086E-01;
    COFD[1343] = 1.51454720E-02;
    COFD[1344] = -1.84591064E+01;
    COFD[1345] = 4.39561236E+00;
    COFD[1346] = -3.48858889E-01;
    COFD[1347] = 1.48190509E-02;
    COFD[1348] = -1.91368159E+01;
    COFD[1349] = 4.56508921E+00;
    COFD[1350] = -3.69033080E-01;
    COFD[1351] = 1.56223865E-02;
    COFD[1352] = -1.91050585E+01;
    COFD[1353] = 4.55886937E+00;
    COFD[1354] = -3.68301178E-01;
    COFD[1355] = 1.55936129E-02;
    COFD[1356] = -1.91080950E+01;
    COFD[1357] = 4.55886937E+00;
    COFD[1358] = -3.68301178E-01;
    COFD[1359] = 1.55936129E-02;
    COFD[1360] = -1.92212680E+01;
    COFD[1361] = 4.54955804E+00;
    COFD[1362] = -3.67206725E-01;
    COFD[1363] = 1.55506376E-02;
    COFD[1364] = -1.99543652E+01;
    COFD[1365] = 4.83774371E+00;
    COFD[1366] = -4.01759538E-01;
    COFD[1367] = 1.69365310E-02;
    COFD[1368] = -2.02463222E+01;
    COFD[1369] = 4.94688197E+00;
    COFD[1370] = -4.14116294E-01;
    COFD[1371] = 1.74006856E-02;
    COFD[1372] = -1.97966756E+01;
    COFD[1373] = 4.72725308E+00;
    COFD[1374] = -3.88865048E-01;
    COFD[1375] = 1.64346151E-02;
    COFD[1376] = -1.97045701E+01;
    COFD[1377] = 4.69156097E+00;
    COFD[1378] = -3.84562356E-01;
    COFD[1379] = 1.62609922E-02;
    COFD[1380] = -2.00697375E+01;
    COFD[1381] = 4.85107798E+00;
    COFD[1382] = -4.03289235E-01;
    COFD[1383] = 1.69949072E-02;
    COFD[1384] = -1.99421880E+01;
    COFD[1385] = 4.94937905E+00;
    COFD[1386] = -4.14409534E-01;
    COFD[1387] = 1.74121859E-02;
    COFD[1388] = -1.99461190E+01;
    COFD[1389] = 4.94937905E+00;
    COFD[1390] = -4.14409534E-01;
    COFD[1391] = 1.74121859E-02;
    COFD[1392] = -1.98245633E+01;
    COFD[1393] = 4.85730922E+00;
    COFD[1394] = -4.04004009E-01;
    COFD[1395] = 1.70221852E-02;
    COFD[1396] = -2.02385185E+01;
    COFD[1397] = 4.92216751E+00;
    COFD[1398] = -4.11213729E-01;
    COFD[1399] = 1.72868483E-02;
    COFD[1400] = -1.85705503E+01;
    COFD[1401] = 4.42984539E+00;
    COFD[1402] = -3.53041284E-01;
    COFD[1403] = 1.49902085E-02;
    COFD[1404] = -1.86697931E+01;
    COFD[1405] = 4.46203550E+00;
    COFD[1406] = -3.56905086E-01;
    COFD[1407] = 1.51454720E-02;
    COFD[1408] = -2.05043391E+01;
    COFD[1409] = 4.99719781E+00;
    COFD[1410] = -4.20035873E-01;
    COFD[1411] = 1.76333794E-02;
    COFD[1412] = -2.05899148E+01;
    COFD[1413] = 4.99533607E+00;
    COFD[1414] = -4.19816362E-01;
    COFD[1415] = 1.76247292E-02;
    COFD[1416] = -1.98276037E+01;
    COFD[1417] = 4.85730922E+00;
    COFD[1418] = -4.04004009E-01;
    COFD[1419] = 1.70221852E-02;
    COFD[1420] = -2.06819784E+01;
    COFD[1421] = 4.99569563E+00;
    COFD[1422] = -4.19858765E-01;
    COFD[1423] = 1.76264005E-02;
    COFD[1424] = -2.05171137E+01;
    COFD[1425] = 4.96961006E+00;
    COFD[1426] = -4.16785706E-01;
    COFD[1427] = 1.75054037E-02;
    COFD[1428] = -2.05041077E+01;
    COFD[1429] = 4.89401297E+00;
    COFD[1430] = -4.08087781E-01;
    COFD[1431] = 1.71722700E-02;
    COFD[1432] = -2.05029990E+01;
    COFD[1433] = 4.89401297E+00;
    COFD[1434] = -4.08087781E-01;
    COFD[1435] = 1.71722700E-02;
    COFD[1436] = -2.12596231E+01;
    COFD[1437] = 5.13107283E+00;
    COFD[1438] = -4.34832558E-01;
    COFD[1439] = 1.81741958E-02;
    COFD[1440] = -2.16246599E+01;
    COFD[1441] = 5.18371903E+00;
    COFD[1442] = -4.40201151E-01;
    COFD[1443] = 1.83498534E-02;
    COFD[1444] = -2.09238627E+01;
    COFD[1445] = 5.02943557E+00;
    COFD[1446] = -4.23811924E-01;
    COFD[1447] = 1.77810575E-02;
    COFD[1448] = -2.14243807E+01;
    COFD[1449] = 5.16859564E+00;
    COFD[1450] = -4.38877619E-01;
    COFD[1451] = 1.83174566E-02;
    COFD[1452] = -1.45806444E+01;
    COFD[1453] = 3.07340947E+00;
    COFD[1454] = -1.84805140E-01;
    COFD[1455] = 8.01973216E-03;
    COFD[1456] = -1.69855694E+01;
    COFD[1457] = 4.58049527E+00;
    COFD[1458] = -3.70921885E-01;
    COFD[1459] = 1.56999324E-02;
    COFD[1460] = -1.33494365E+01;
    COFD[1461] = 3.04769365E+00;
    COFD[1462] = -1.81586731E-01;
    COFD[1463] = 7.88540063E-03;
    COFD[1464] = -1.64305407E+01;
    COFD[1465] = 3.86315631E+00;
    COFD[1466] = -2.83800980E-01;
    COFD[1467] = 1.21619703E-02;
    COFD[1468] = -1.77768991E+01;
    COFD[1469] = 4.18471443E+00;
    COFD[1470] = -3.22640353E-01;
    COFD[1471] = 1.37280239E-02;
    COFD[1472] = -1.64527598E+01;
    COFD[1473] = 3.86315631E+00;
    COFD[1474] = -2.83800980E-01;
    COFD[1475] = 1.21619703E-02;
    COFD[1476] = -2.14098987E+01;
    COFD[1477] = 5.30121376E+00;
    COFD[1478] = -4.21263191E-01;
    COFD[1479] = 1.62329486E-02;
    COFD[1480] = -1.74636259E+01;
    COFD[1481] = 4.06842532E+00;
    COFD[1482] = -3.08492475E-01;
    COFD[1483] = 1.31534275E-02;
    COFD[1484] = -2.08457859E+01;
    COFD[1485] = 5.18382353E+00;
    COFD[1486] = -4.40210428E-01;
    COFD[1487] = 1.83500882E-02;
    COFD[1488] = -1.86951919E+01;
    COFD[1489] = 4.56871889E+00;
    COFD[1490] = -3.69460525E-01;
    COFD[1491] = 1.56392048E-02;
    COFD[1492] = -1.86463588E+01;
    COFD[1493] = 4.54862236E+00;
    COFD[1494] = -3.67096821E-01;
    COFD[1495] = 1.55463252E-02;
    COFD[1496] = -1.77858192E+01;
    COFD[1497] = 4.18471443E+00;
    COFD[1498] = -3.22640353E-01;
    COFD[1499] = 1.37280239E-02;
    COFD[1500] = -1.77943584E+01;
    COFD[1501] = 4.18471443E+00;
    COFD[1502] = -3.22640353E-01;
    COFD[1503] = 1.37280239E-02;
    COFD[1504] = -2.19981784E+01;
    COFD[1505] = 5.48342675E+00;
    COFD[1506] = -4.57782070E-01;
    COFD[1507] = 1.82957850E-02;
    COFD[1508] = -2.09219827E+01;
    COFD[1509] = 5.19204548E+00;
    COFD[1510] = -4.40915128E-01;
    COFD[1511] = 1.83664930E-02;
    COFD[1512] = -2.05980858E+01;
    COFD[1513] = 5.16997013E+00;
    COFD[1514] = -4.38998744E-01;
    COFD[1515] = 1.83204571E-02;
    COFD[1516] = -2.09117850E+01;
    COFD[1517] = 5.19204548E+00;
    COFD[1518] = -4.40915128E-01;
    COFD[1519] = 1.83664930E-02;
    COFD[1520] = -2.09138237E+01;
    COFD[1521] = 5.25757923E+00;
    COFD[1522] = -4.47803407E-01;
    COFD[1523] = 1.86021883E-02;
    COFD[1524] = -2.09263167E+01;
    COFD[1525] = 5.25757923E+00;
    COFD[1526] = -4.47803407E-01;
    COFD[1527] = 1.86021883E-02;
    COFD[1528] = -2.20081706E+01;
    COFD[1529] = 5.48782481E+00;
    COFD[1530] = -4.59567900E-01;
    COFD[1531] = 1.84190115E-02;
    COFD[1532] = -2.20938193E+01;
    COFD[1533] = 5.49242695E+00;
    COFD[1534] = -4.63877105E-01;
    COFD[1535] = 1.87472502E-02;
    COFD[1536] = -1.88249737E+01;
    COFD[1537] = 4.63673045E+00;
    COFD[1538] = -3.77870004E-01;
    COFD[1539] = 1.59875055E-02;
    COFD[1540] = -2.21056610E+01;
    COFD[1541] = 5.49242695E+00;
    COFD[1542] = -4.63877105E-01;
    COFD[1543] = 1.87472502E-02;
    COFD[1544] = -2.17754519E+01;
    COFD[1545] = 5.43321422E+00;
    COFD[1546] = -4.65643875E-01;
    COFD[1547] = 1.91852144E-02;
    COFD[1548] = -2.17754519E+01;
    COFD[1549] = 5.43321422E+00;
    COFD[1550] = -4.65643875E-01;
    COFD[1551] = 1.91852144E-02;
    COFD[1552] = -2.16482216E+01;
    COFD[1553] = 5.40849808E+00;
    COFD[1554] = -4.63989845E-01;
    COFD[1555] = 1.91734161E-02;
    COFD[1556] = -2.20785179E+01;
    COFD[1557] = 5.45354415E+00;
    COFD[1558] = -4.65291244E-01;
    COFD[1559] = 1.90545098E-02;
    COFD[1560] = -2.20577824E+01;
    COFD[1561] = 5.45198096E+00;
    COFD[1562] = -4.65259344E-01;
    COFD[1563] = 1.90596774E-02;
    COFD[1564] = -2.20617871E+01;
    COFD[1565] = 5.45198096E+00;
    COFD[1566] = -4.65259344E-01;
    COFD[1567] = 1.90596774E-02;
    COFD[1568] = -2.21891699E+01;
    COFD[1569] = 5.44981030E+00;
    COFD[1570] = -4.65238257E-01;
    COFD[1571] = 1.90687553E-02;
    COFD[1572] = -2.23824354E+01;
    COFD[1573] = 5.49217181E+00;
    COFD[1574] = -4.63882677E-01;
    COFD[1575] = 1.87489631E-02;
    COFD[1576] = -2.24145023E+01;
    COFD[1577] = 5.49116789E+00;
    COFD[1578] = -4.60872367E-01;
    COFD[1579] = 1.85086095E-02;
    COFD[1580] = -2.24163735E+01;
    COFD[1581] = 5.46720411E+00;
    COFD[1582] = -4.63664606E-01;
    COFD[1583] = 1.88543495E-02;
    COFD[1584] = -2.24050898E+01;
    COFD[1585] = 5.46747549E+00;
    COFD[1586] = -4.64570865E-01;
    COFD[1587] = 1.89267998E-02;
    COFD[1588] = -2.24829635E+01;
    COFD[1589] = 5.49361964E+00;
    COFD[1590] = -4.63657707E-01;
    COFD[1591] = 1.87235131E-02;
    COFD[1592] = -2.21101148E+01;
    COFD[1593] = 5.48782481E+00;
    COFD[1594] = -4.59567900E-01;
    COFD[1595] = 1.84190115E-02;
    COFD[1596] = -2.21152117E+01;
    COFD[1597] = 5.48782481E+00;
    COFD[1598] = -4.59567900E-01;
    COFD[1599] = 1.84190115E-02;
    COFD[1600] = -2.22150103E+01;
    COFD[1601] = 5.49406269E+00;
    COFD[1602] = -4.63519837E-01;
    COFD[1603] = 1.87101015E-02;
    COFD[1604] = -2.24591465E+01;
    COFD[1605] = 5.49093616E+00;
    COFD[1606] = -4.61489559E-01;
    COFD[1607] = 1.85600323E-02;
    COFD[1608] = -2.17175223E+01;
    COFD[1609] = 5.42282789E+00;
    COFD[1610] = -4.65108864E-01;
    COFD[1611] = 1.91933632E-02;
    COFD[1612] = -2.17687426E+01;
    COFD[1613] = 5.43321422E+00;
    COFD[1614] = -4.65643875E-01;
    COFD[1615] = 1.91852144E-02;
    COFD[1616] = -2.25606430E+01;
    COFD[1617] = 5.48480058E+00;
    COFD[1618] = -4.58431586E-01;
    COFD[1619] = 1.83416243E-02;
    COFD[1620] = -2.26462274E+01;
    COFD[1621] = 5.48561354E+00;
    COFD[1622] = -4.58748903E-01;
    COFD[1623] = 1.83633678E-02;
    COFD[1624] = -2.22190198E+01;
    COFD[1625] = 5.49406269E+00;
    COFD[1626] = -4.63519837E-01;
    COFD[1627] = 1.87101015E-02;
    COFD[1628] = -2.27524699E+01;
    COFD[1629] = 5.48417137E+00;
    COFD[1630] = -4.58172226E-01;
    COFD[1631] = 1.83236937E-02;
    COFD[1632] = -2.26588618E+01;
    COFD[1633] = 5.48550227E+00;
    COFD[1634] = -4.58704694E-01;
    COFD[1635] = 1.83603287E-02;
    COFD[1636] = -2.28305301E+01;
    COFD[1637] = 5.49066995E+00;
    COFD[1638] = -4.61705340E-01;
    COFD[1639] = 1.85789333E-02;
    COFD[1640] = -2.28289821E+01;
    COFD[1641] = 5.49066995E+00;
    COFD[1642] = -4.61705340E-01;
    COFD[1643] = 1.85789333E-02;
    COFD[1644] = -2.29703215E+01;
    COFD[1645] = 5.45099312E+00;
    COFD[1646] = -4.48927722E-01;
    COFD[1647] = 1.77397382E-02;
    COFD[1648] = -2.31435861E+01;
    COFD[1649] = 5.41522402E+00;
    COFD[1650] = -4.41568977E-01;
    COFD[1651] = 1.73190088E-02;
    COFD[1652] = -2.29289604E+01;
    COFD[1653] = 5.48025735E+00;
    COFD[1654] = -4.56422639E-01;
    COFD[1655] = 1.82014258E-02;
    COFD[1656] = -2.30145662E+01;
    COFD[1657] = 5.43378566E+00;
    COFD[1658] = -4.45159774E-01;
    COFD[1659] = 1.75190628E-02;
    COFD[1660] = -1.74370371E+01;
    COFD[1661] = 4.06188586E+00;
    COFD[1662] = -3.07715468E-01;
    COFD[1663] = 1.31226885E-02;
    COFD[1664] = -1.51283685E+01;
    COFD[1665] = 3.92295205E+00;
    COFD[1666] = -2.90942506E-01;
    COFD[1667] = 1.24465450E-02;
    COFD[1668] = -1.20163448E+01;
    COFD[1669] = 2.55772522E+00;
    COFD[1670] = -1.18119456E-01;
    COFD[1671] = 5.14119183E-03;
    COFD[1672] = -1.46233077E+01;
    COFD[1673] = 3.27714313E+00;
    COFD[1674] = -2.10795461E-01;
    COFD[1675] = 9.12600427E-03;
    COFD[1676] = -1.58176414E+01;
    COFD[1677] = 3.58720055E+00;
    COFD[1678] = -2.49572308E-01;
    COFD[1679] = 1.07445741E-02;
    COFD[1680] = -1.46378741E+01;
    COFD[1681] = 3.27714313E+00;
    COFD[1682] = -2.10795461E-01;
    COFD[1683] = 9.12600427E-03;
    COFD[1684] = -2.08999981E+01;
    COFD[1685] = 5.44151563E+00;
    COFD[1686] = -4.65460498E-01;
    COFD[1687] = 1.91284180E-02;
    COFD[1688] = -1.55192394E+01;
    COFD[1689] = 3.46464286E+00;
    COFD[1690] = -2.33957556E-01;
    COFD[1691] = 1.00801477E-02;
    COFD[1692] = -1.86951919E+01;
    COFD[1693] = 4.56871889E+00;
    COFD[1694] = -3.69460525E-01;
    COFD[1695] = 1.56392048E-02;
    COFD[1696] = -1.66707532E+01;
    COFD[1697] = 3.91662496E+00;
    COFD[1698] = -2.90182401E-01;
    COFD[1699] = 1.24160479E-02;
    COFD[1700] = -1.66255694E+01;
    COFD[1701] = 3.90031906E+00;
    COFD[1702] = -2.88230661E-01;
    COFD[1703] = 1.23380648E-02;
    COFD[1704] = -1.58225464E+01;
    COFD[1705] = 3.58720055E+00;
    COFD[1706] = -2.49572308E-01;
    COFD[1707] = 1.07445741E-02;
    COFD[1708] = -1.58272050E+01;
    COFD[1709] = 3.58720055E+00;
    COFD[1710] = -2.49572308E-01;
    COFD[1711] = 1.07445741E-02;
    COFD[1712] = -2.10007851E+01;
    COFD[1713] = 5.36906158E+00;
    COFD[1714] = -4.60039482E-01;
    COFD[1715] = 1.90471208E-02;
    COFD[1716] = -1.88552742E+01;
    COFD[1717] = 4.59285496E+00;
    COFD[1718] = -3.72450186E-01;
    COFD[1719] = 1.57632217E-02;
    COFD[1720] = -1.84813514E+01;
    COFD[1721] = 4.54318703E+00;
    COFD[1722] = -3.66458643E-01;
    COFD[1723] = 1.55212948E-02;
    COFD[1724] = -1.88495268E+01;
    COFD[1725] = 4.59285496E+00;
    COFD[1726] = -3.72450186E-01;
    COFD[1727] = 1.57632217E-02;
    COFD[1728] = -1.89650398E+01;
    COFD[1729] = 4.70175320E+00;
    COFD[1730] = -3.85802299E-01;
    COFD[1731] = 1.63115029E-02;
    COFD[1732] = -1.89723565E+01;
    COFD[1733] = 4.70175320E+00;
    COFD[1734] = -3.85802299E-01;
    COFD[1735] = 1.63115029E-02;
    COFD[1736] = -2.09231900E+01;
    COFD[1737] = 5.33874924E+00;
    COFD[1738] = -4.56879272E-01;
    COFD[1739] = 1.89397603E-02;
    COFD[1740] = -2.06886834E+01;
    COFD[1741] = 5.21843630E+00;
    COFD[1742] = -4.43374339E-01;
    COFD[1743] = 1.84354402E-02;
    COFD[1744] = -1.65966953E+01;
    COFD[1745] = 3.95555732E+00;
    COFD[1746] = -2.94894134E-01;
    COFD[1747] = 1.26066350E-02;
    COFD[1748] = -2.06947476E+01;
    COFD[1749] = 5.21843630E+00;
    COFD[1750] = -4.43374339E-01;
    COFD[1751] = 1.84354402E-02;
    COFD[1752] = -1.98761303E+01;
    COFD[1753] = 4.93512898E+00;
    COFD[1754] = -4.12735698E-01;
    COFD[1755] = 1.73465233E-02;
    COFD[1756] = -1.98761303E+01;
    COFD[1757] = 4.93512898E+00;
    COFD[1758] = -4.12735698E-01;
    COFD[1759] = 1.73465233E-02;
    COFD[1760] = -1.96749452E+01;
    COFD[1761] = 4.87934939E+00;
    COFD[1762] = -4.06539237E-01;
    COFD[1763] = 1.71192321E-02;
    COFD[1764] = -2.03731422E+01;
    COFD[1765] = 5.05714666E+00;
    COFD[1766] = -4.26960034E-01;
    COFD[1767] = 1.78999526E-02;
    COFD[1768] = -2.03431655E+01;
    COFD[1769] = 5.05232294E+00;
    COFD[1770] = -4.26452915E-01;
    COFD[1771] = 1.78825996E-02;
    COFD[1772] = -2.03450948E+01;
    COFD[1773] = 5.05232294E+00;
    COFD[1774] = -4.26452915E-01;
    COFD[1775] = 1.78825996E-02;
    COFD[1776] = -2.04523945E+01;
    COFD[1777] = 5.04244535E+00;
    COFD[1778] = -4.25330213E-01;
    COFD[1779] = 1.78402047E-02;
    COFD[1780] = -2.09365593E+01;
    COFD[1781] = 5.21672969E+00;
    COFD[1782] = -4.43175115E-01;
    COFD[1783] = 1.84276681E-02;
    COFD[1784] = -2.08891203E+01;
    COFD[1785] = 5.21570257E+00;
    COFD[1786] = -4.43055160E-01;
    COFD[1787] = 1.84229862E-02;
    COFD[1788] = -2.08689918E+01;
    COFD[1789] = 5.15366259E+00;
    COFD[1790] = -4.37268432E-01;
    COFD[1791] = 1.82604618E-02;
    COFD[1792] = -2.07942327E+01;
    COFD[1793] = 5.12723540E+00;
    COFD[1794] = -4.34411432E-01;
    COFD[1795] = 1.81589442E-02;
    COFD[1796] = -2.10468267E+01;
    COFD[1797] = 5.23281177E+00;
    COFD[1798] = -4.45013021E-01;
    COFD[1799] = 1.84976592E-02;
    COFD[1800] = -2.09768682E+01;
    COFD[1801] = 5.33874924E+00;
    COFD[1802] = -4.56879272E-01;
    COFD[1803] = 1.89397603E-02;
    COFD[1804] = -2.09794158E+01;
    COFD[1805] = 5.33874924E+00;
    COFD[1806] = -4.56879272E-01;
    COFD[1807] = 1.89397603E-02;
    COFD[1808] = -2.08233373E+01;
    COFD[1809] = 5.23993213E+00;
    COFD[1810] = -4.45818789E-01;
    COFD[1811] = 1.85280084E-02;
    COFD[1812] = -2.09160351E+01;
    COFD[1813] = 5.20438984E+00;
    COFD[1814] = -4.41945569E-01;
    COFD[1815] = 1.83888515E-02;
    COFD[1816] = -1.97629737E+01;
    COFD[1817] = 4.90161091E+00;
    COFD[1818] = -4.08878655E-01;
    COFD[1819] = 1.71987880E-02;
    COFD[1820] = -1.98726217E+01;
    COFD[1821] = 4.93512898E+00;
    COFD[1822] = -4.12735698E-01;
    COFD[1823] = 1.73465233E-02;
    COFD[1824] = -2.12336858E+01;
    COFD[1825] = 5.29671000E+00;
    COFD[1826] = -4.52208613E-01;
    COFD[1827] = 1.87672717E-02;
    COFD[1828] = -2.12409720E+01;
    COFD[1829] = 5.26987677E+00;
    COFD[1830] = -4.49189697E-01;
    COFD[1831] = 1.86541737E-02;
    COFD[1832] = -2.08252693E+01;
    COFD[1833] = 5.23993213E+00;
    COFD[1834] = -4.45818789E-01;
    COFD[1835] = 1.85280084E-02;
    COFD[1836] = -2.15010463E+01;
    COFD[1837] = 5.32769529E+00;
    COFD[1838] = -4.55658342E-01;
    COFD[1839] = 1.88949932E-02;
    COFD[1840] = -2.14961087E+01;
    COFD[1841] = 5.35264246E+00;
    COFD[1842] = -4.58331575E-01;
    COFD[1843] = 1.89892757E-02;
    COFD[1844] = -2.14866930E+01;
    COFD[1845] = 5.28483859E+00;
    COFD[1846] = -4.50875157E-01;
    COFD[1847] = 1.87173937E-02;
    COFD[1848] = -2.14860306E+01;
    COFD[1849] = 5.28483859E+00;
    COFD[1850] = -4.50875157E-01;
    COFD[1851] = 1.87173937E-02;
    COFD[1852] = -2.20377114E+01;
    COFD[1853] = 5.43771555E+00;
    COFD[1854] = -4.65617512E-01;
    COFD[1855] = 1.91604198E-02;
    COFD[1856] = -2.23095159E+01;
    COFD[1857] = 5.45157330E+00;
    COFD[1858] = -4.65253665E-01;
    COFD[1859] = 1.90612411E-02;
    COFD[1860] = -2.18504787E+01;
    COFD[1861] = 5.39546703E+00;
    COFD[1862] = -4.62825371E-01;
    COFD[1863] = 1.91432238E-02;
    COFD[1864] = -2.21249716E+01;
    COFD[1865] = 5.44426165E+00;
    COFD[1866] = -4.65323191E-01;
    COFD[1867] = 1.91034051E-02;
    COFD[1868] = -1.54972198E+01;
    COFD[1869] = 3.45957322E+00;
    COFD[1870] = -2.33335381E-01;
    COFD[1871] = 1.00546680E-02;
    COFD[1872] = -1.50643750E+01;
    COFD[1873] = 3.90641370E+00;
    COFD[1874] = -2.88957677E-01;
    COFD[1875] = 1.23669987E-02;
    COFD[1876] = -1.19609333E+01;
    COFD[1877] = 2.54091639E+00;
    COFD[1878] = -1.15850475E-01;
    COFD[1879] = 5.03893346E-03;
    COFD[1880] = -1.45714634E+01;
    COFD[1881] = 3.25814060E+00;
    COFD[1882] = -2.08402129E-01;
    COFD[1883] = 9.02544298E-03;
    COFD[1884] = -1.57555897E+01;
    COFD[1885] = 3.56140952E+00;
    COFD[1886] = -2.46282892E-01;
    COFD[1887] = 1.06044757E-02;
    COFD[1888] = -1.45865248E+01;
    COFD[1889] = 3.25814060E+00;
    COFD[1890] = -2.08402129E-01;
    COFD[1891] = 9.02544298E-03;
    COFD[1892] = -2.11989257E+01;
    COFD[1893] = 5.48294376E+00;
    COFD[1894] = -4.63910522E-01;
    COFD[1895] = 1.87967250E-02;
    COFD[1896] = -1.54796758E+01;
    COFD[1897] = 3.44821689E+00;
    COFD[1898] = -2.31932994E-01;
    COFD[1899] = 9.99688061E-03;
    COFD[1900] = -1.86463588E+01;
    COFD[1901] = 4.54862236E+00;
    COFD[1902] = -3.67096821E-01;
    COFD[1903] = 1.55463252E-02;
    COFD[1904] = -1.66255694E+01;
    COFD[1905] = 3.90031906E+00;
    COFD[1906] = -2.88230661E-01;
    COFD[1907] = 1.23380648E-02;
    COFD[1908] = -1.65813658E+01;
    COFD[1909] = 3.88425001E+00;
    COFD[1910] = -2.86312476E-01;
    COFD[1911] = 1.22616627E-02;
    COFD[1912] = -1.57607149E+01;
    COFD[1913] = 3.56140952E+00;
    COFD[1914] = -2.46282892E-01;
    COFD[1915] = 1.06044757E-02;
    COFD[1916] = -1.57655847E+01;
    COFD[1917] = 3.56140952E+00;
    COFD[1918] = -2.46282892E-01;
    COFD[1919] = 1.06044757E-02;
    COFD[1920] = -2.09538981E+01;
    COFD[1921] = 5.35250010E+00;
    COFD[1922] = -4.58316703E-01;
    COFD[1923] = 1.89887690E-02;
    COFD[1924] = -1.87841274E+01;
    COFD[1925] = 4.56428962E+00;
    COFD[1926] = -3.68938951E-01;
    COFD[1927] = 1.56186844E-02;
    COFD[1928] = -1.84296831E+01;
    COFD[1929] = 4.52338181E+00;
    COFD[1930] = -3.64129854E-01;
    COFD[1931] = 1.54297700E-02;
    COFD[1932] = -1.87781303E+01;
    COFD[1933] = 4.56428962E+00;
    COFD[1934] = -3.68938951E-01;
    COFD[1935] = 1.56186844E-02;
    COFD[1936] = -1.88959377E+01;
    COFD[1937] = 4.67519664E+00;
    COFD[1938] = -3.82581158E-01;
    COFD[1939] = 1.61807120E-02;
    COFD[1940] = -1.89035556E+01;
    COFD[1941] = 4.67519664E+00;
    COFD[1942] = -3.82581158E-01;
    COFD[1943] = 1.61807120E-02;
    COFD[1944] = -2.08679661E+01;
    COFD[1945] = 5.31839685E+00;
    COFD[1946] = -4.54623948E-01;
    COFD[1947] = 1.88567295E-02;
    COFD[1948] = -2.06489560E+01;
    COFD[1949] = 5.20348257E+00;
    COFD[1950] = -4.41867814E-01;
    COFD[1951] = 1.83870405E-02;
    COFD[1952] = -1.65522787E+01;
    COFD[1953] = 3.93761009E+00;
    COFD[1954] = -2.92710565E-01;
    COFD[1955] = 1.25177995E-02;
    COFD[1956] = -2.06553163E+01;
    COFD[1957] = 5.20348257E+00;
    COFD[1958] = -4.41867814E-01;
    COFD[1959] = 1.83870405E-02;
    COFD[1960] = -1.98121403E+01;
    COFD[1961] = 4.90949838E+00;
    COFD[1962] = -4.09710987E-01;
    COFD[1963] = 1.72272850E-02;
    COFD[1964] = -1.98121403E+01;
    COFD[1965] = 4.90949838E+00;
    COFD[1966] = -4.09710987E-01;
    COFD[1967] = 1.72272850E-02;
    COFD[1968] = -1.96231480E+01;
    COFD[1969] = 4.85846969E+00;
    COFD[1970] = -4.04137136E-01;
    COFD[1971] = 1.70272658E-02;
    COFD[1972] = -2.03251784E+01;
    COFD[1973] = 5.03716590E+00;
    COFD[1974] = -4.24714370E-01;
    COFD[1975] = 1.78162251E-02;
    COFD[1976] = -2.02895141E+01;
    COFD[1977] = 5.02992150E+00;
    COFD[1978] = -4.23868677E-01;
    COFD[1979] = 1.77832701E-02;
    COFD[1980] = -2.02915437E+01;
    COFD[1981] = 5.02992150E+00;
    COFD[1982] = -4.23868677E-01;
    COFD[1983] = 1.77832701E-02;
    COFD[1984] = -2.03971631E+01;
    COFD[1985] = 5.01889159E+00;
    COFD[1986] = -4.22579631E-01;
    COFD[1987] = 1.77329812E-02;
    COFD[1988] = -2.09026212E+01;
    COFD[1989] = 5.20282016E+00;
    COFD[1990] = -4.41811086E-01;
    COFD[1991] = 1.83857216E-02;
    COFD[1992] = -2.12234786E+01;
    COFD[1993] = 5.31936752E+00;
    COFD[1994] = -4.54732039E-01;
    COFD[1995] = 1.88607327E-02;
    COFD[1996] = -2.08290042E+01;
    COFD[1997] = 5.13624657E+00;
    COFD[1998] = -4.35395745E-01;
    COFD[1999] = 1.81943797E-02;
    COFD[2000] = -2.07580489E+01;
    COFD[2001] = 5.11080865E+00;
    COFD[2002] = -4.32626175E-01;
    COFD[2003] = 1.80950687E-02;
    COFD[2004] = -2.09922482E+01;
    COFD[2005] = 5.20990101E+00;
    COFD[2006] = -4.42427060E-01;
    COFD[2007] = 1.84006140E-02;
    COFD[2008] = -2.09241822E+01;
    COFD[2009] = 5.31839685E+00;
    COFD[2010] = -4.54623948E-01;
    COFD[2011] = 1.88567295E-02;
    COFD[2012] = -2.09268575E+01;
    COFD[2013] = 5.31839685E+00;
    COFD[2014] = -4.54623948E-01;
    COFD[2015] = 1.88567295E-02;
    COFD[2016] = -2.07627166E+01;
    COFD[2017] = 5.21569495E+00;
    COFD[2018] = -4.43054270E-01;
    COFD[2019] = 1.84229515E-02;
    COFD[2020] = -2.12249561E+01;
    COFD[2021] = 5.29857429E+00;
    COFD[2022] = -4.52417024E-01;
    COFD[2023] = 1.87750282E-02;
    COFD[2024] = -1.97245207E+01;
    COFD[2025] = 4.88654569E+00;
    COFD[2026] = -4.07307743E-01;
    COFD[2027] = 1.71459608E-02;
    COFD[2028] = -1.98084643E+01;
    COFD[2029] = 4.90949838E+00;
    COFD[2030] = -4.09710987E-01;
    COFD[2031] = 1.72272850E-02;
    COFD[2032] = -2.14398976E+01;
    COFD[2033] = 5.35638301E+00;
    COFD[2034] = -4.58721935E-01;
    COFD[2035] = 1.90025575E-02;
    COFD[2036] = -2.15282219E+01;
    COFD[2037] = 5.35564147E+00;
    COFD[2038] = -4.58644607E-01;
    COFD[2039] = 1.89999292E-02;
    COFD[2040] = -2.07647491E+01;
    COFD[2041] = 5.21569495E+00;
    COFD[2042] = -4.43054270E-01;
    COFD[2043] = 1.84229515E-02;
    COFD[2044] = -2.16063225E+01;
    COFD[2045] = 5.35467220E+00;
    COFD[2046] = -4.58543488E-01;
    COFD[2047] = 1.89964903E-02;
    COFD[2048] = -2.14562897E+01;
    COFD[2049] = 5.33532691E+00;
    COFD[2050] = -4.56505568E-01;
    COFD[2051] = 1.89262566E-02;
    COFD[2052] = -2.14364392E+01;
    COFD[2053] = 5.26245213E+00;
    COFD[2054] = -4.48352237E-01;
    COFD[2055] = 1.86227443E-02;
    COFD[2056] = -2.14357385E+01;
    COFD[2057] = 5.26245213E+00;
    COFD[2058] = -4.48352237E-01;
    COFD[2059] = 1.86227443E-02;
    COFD[2060] = -2.20523530E+01;
    COFD[2061] = 5.43694233E+00;
    COFD[2062] = -4.65634840E-01;
    COFD[2063] = 1.91657287E-02;
    COFD[2064] = -2.23193198E+01;
    COFD[2065] = 5.44963079E+00;
    COFD[2066] = -4.65238095E-01;
    COFD[2067] = 1.90696359E-02;
    COFD[2068] = -2.18110440E+01;
    COFD[2069] = 5.37758536E+00;
    COFD[2070] = -4.60945701E-01;
    COFD[2071] = 1.90787306E-02;
    COFD[2072] = -2.21451593E+01;
    COFD[2073] = 5.44423992E+00;
    COFD[2074] = -4.65323901E-01;
    COFD[2075] = 1.91035720E-02;
    COFD[2076] = -1.54573340E+01;
    COFD[2077] = 3.44303981E+00;
    COFD[2078] = -2.31296721E-01;
    COFD[2079] = 9.97079557E-03;
    COFD[2080] = -1.41188260E+01;
    COFD[2081] = 3.59720720E+00;
    COFD[2082] = -2.50852412E-01;
    COFD[2083] = 1.07992653E-02;
    COFD[2084] = -1.12823274E+01;
    COFD[2085] = 2.30017030E+00;
    COFD[2086] = -8.34928497E-02;
    COFD[2087] = 3.58611589E-03;
    COFD[2088] = -1.38853857E+01;
    COFD[2089] = 2.98298723E+00;
    COFD[2090] = -1.73428215E-01;
    COFD[2091] = 7.54224684E-03;
    COFD[2092] = -1.51058803E+01;
    COFD[2093] = 3.27853241E+00;
    COFD[2094] = -2.10970532E-01;
    COFD[2095] = 9.13336403E-03;
    COFD[2096] = -1.39057537E+01;
    COFD[2097] = 2.98298723E+00;
    COFD[2098] = -1.73428215E-01;
    COFD[2099] = 7.54224684E-03;
    COFD[2100] = -2.01625839E+01;
    COFD[2101] = 5.20184446E+00;
    COFD[2102] = -4.41727448E-01;
    COFD[2103] = 1.83837722E-02;
    COFD[2104] = -1.48478190E+01;
    COFD[2105] = 3.17494843E+00;
    COFD[2106] = -1.97752410E-01;
    COFD[2107] = 8.57054106E-03;
    COFD[2108] = -1.77858192E+01;
    COFD[2109] = 4.18471443E+00;
    COFD[2110] = -3.22640353E-01;
    COFD[2111] = 1.37280239E-02;
    COFD[2112] = -1.58225464E+01;
    COFD[2113] = 3.58720055E+00;
    COFD[2114] = -2.49572308E-01;
    COFD[2115] = 1.07445741E-02;
    COFD[2116] = -1.57607149E+01;
    COFD[2117] = 3.56140952E+00;
    COFD[2118] = -2.46282892E-01;
    COFD[2119] = 1.06044757E-02;
    COFD[2120] = -1.51136939E+01;
    COFD[2121] = 3.27853241E+00;
    COFD[2122] = -2.10970532E-01;
    COFD[2123] = 9.13336403E-03;
    COFD[2124] = -1.51211577E+01;
    COFD[2125] = 3.27853241E+00;
    COFD[2126] = -2.10970532E-01;
    COFD[2127] = 9.13336403E-03;
    COFD[2128] = -2.02493843E+01;
    COFD[2129] = 5.09185655E+00;
    COFD[2130] = -4.30641631E-01;
    COFD[2131] = 1.80275784E-02;
    COFD[2132] = -1.79077036E+01;
    COFD[2133] = 4.20508193E+00;
    COFD[2134] = -3.25158699E-01;
    COFD[2135] = 1.38321131E-02;
    COFD[2136] = -1.75143894E+01;
    COFD[2137] = 4.15347096E+00;
    COFD[2138] = -3.18813734E-01;
    COFD[2139] = 1.35714974E-02;
    COFD[2140] = -1.78987105E+01;
    COFD[2141] = 4.20508193E+00;
    COFD[2142] = -3.25158699E-01;
    COFD[2143] = 1.38321131E-02;
    COFD[2144] = -1.79832047E+01;
    COFD[2145] = 4.31259444E+00;
    COFD[2146] = -3.38584827E-01;
    COFD[2147] = 1.43931810E-02;
    COFD[2148] = -1.79943355E+01;
    COFD[2149] = 4.31259444E+00;
    COFD[2150] = -3.38584827E-01;
    COFD[2151] = 1.43931810E-02;
    COFD[2152] = -2.01841845E+01;
    COFD[2153] = 5.06295779E+00;
    COFD[2154] = -4.27579258E-01;
    COFD[2155] = 1.79215546E-02;
    COFD[2156] = -1.99727712E+01;
    COFD[2157] = 4.93673368E+00;
    COFD[2158] = -4.12924135E-01;
    COFD[2159] = 1.73539126E-02;
    COFD[2160] = -1.58717613E+01;
    COFD[2161] = 3.64628858E+00;
    COFD[2162] = -2.57092933E-01;
    COFD[2163] = 1.10642427E-02;
    COFD[2164] = -1.99829374E+01;
    COFD[2165] = 4.93673368E+00;
    COFD[2166] = -4.12924135E-01;
    COFD[2167] = 1.73539126E-02;
    COFD[2168] = -1.89598731E+01;
    COFD[2169] = 4.56070745E+00;
    COFD[2170] = -3.68517396E-01;
    COFD[2171] = 1.56021102E-02;
    COFD[2172] = -1.89598731E+01;
    COFD[2173] = 4.56070745E+00;
    COFD[2174] = -3.68517396E-01;
    COFD[2175] = 1.56021102E-02;
    COFD[2176] = -1.87640840E+01;
    COFD[2177] = 4.50343635E+00;
    COFD[2178] = -3.61818471E-01;
    COFD[2179] = 1.53404608E-02;
    COFD[2180] = -1.95377616E+01;
    COFD[2181] = 4.70710890E+00;
    COFD[2182] = -3.86453789E-01;
    COFD[2183] = 1.63380393E-02;
    COFD[2184] = -1.95007169E+01;
    COFD[2185] = 4.69879097E+00;
    COFD[2186] = -3.85441745E-01;
    COFD[2187] = 1.62968073E-02;
    COFD[2188] = -1.95040932E+01;
    COFD[2189] = 4.69879097E+00;
    COFD[2190] = -3.85441745E-01;
    COFD[2191] = 1.62968073E-02;
    COFD[2192] = -1.96129127E+01;
    COFD[2193] = 4.68635682E+00;
    COFD[2194] = -3.83931845E-01;
    COFD[2195] = 1.62354222E-02;
    COFD[2196] = -2.02561467E+01;
    COFD[2197] = 4.93517265E+00;
    COFD[2198] = -4.12740828E-01;
    COFD[2199] = 1.73467245E-02;
    COFD[2200] = -2.02084294E+01;
    COFD[2201] = 4.93423363E+00;
    COFD[2202] = -4.12630529E-01;
    COFD[2203] = 1.73423980E-02;
    COFD[2204] = -2.01255741E+01;
    COFD[2205] = 4.83483781E+00;
    COFD[2206] = -4.01422871E-01;
    COFD[2207] = 1.69235370E-02;
    COFD[2208] = -2.00433341E+01;
    COFD[2209] = 4.80382495E+00;
    COFD[2210] = -3.97838890E-01;
    COFD[2211] = 1.67856063E-02;
    COFD[2212] = -2.03838522E+01;
    COFD[2213] = 4.95138566E+00;
    COFD[2214] = -4.14645044E-01;
    COFD[2215] = 1.74214162E-02;
    COFD[2216] = -2.02723926E+01;
    COFD[2217] = 5.06295779E+00;
    COFD[2218] = -4.27579258E-01;
    COFD[2219] = 1.79215546E-02;
    COFD[2220] = -2.02767375E+01;
    COFD[2221] = 5.06295779E+00;
    COFD[2222] = -4.27579258E-01;
    COFD[2223] = 1.79215546E-02;
    COFD[2224] = -2.01319872E+01;
    COFD[2225] = 4.95889869E+00;
    COFD[2226] = -4.15526934E-01;
    COFD[2227] = 1.74559870E-02;
    COFD[2228] = -2.02137759E+01;
    COFD[2229] = 4.91354418E+00;
    COFD[2230] = -4.10191969E-01;
    COFD[2231] = 1.72463969E-02;
    COFD[2232] = -1.88597846E+01;
    COFD[2233] = 4.53065425E+00;
    COFD[2234] = -3.64983473E-01;
    COFD[2235] = 1.54632547E-02;
    COFD[2236] = -1.89540782E+01;
    COFD[2237] = 4.56070745E+00;
    COFD[2238] = -3.68517396E-01;
    COFD[2239] = 1.56021102E-02;
    COFD[2240] = -2.05805117E+01;
    COFD[2241] = 5.02104665E+00;
    COFD[2242] = -4.22831627E-01;
    COFD[2243] = 1.77428177E-02;
    COFD[2244] = -2.05856649E+01;
    COFD[2245] = 4.99149256E+00;
    COFD[2246] = -4.19363176E-01;
    COFD[2247] = 1.76068685E-02;
    COFD[2248] = -2.01353679E+01;
    COFD[2249] = 4.95889869E+00;
    COFD[2250] = -4.15526934E-01;
    COFD[2251] = 1.74559870E-02;
    COFD[2252] = -2.08689562E+01;
    COFD[2253] = 5.05404465E+00;
    COFD[2254] = -4.26628831E-01;
    COFD[2255] = 1.78883654E-02;
    COFD[2256] = -2.08516224E+01;
    COFD[2257] = 5.07617318E+00;
    COFD[2258] = -4.28983517E-01;
    COFD[2259] = 1.79703565E-02;
    COFD[2260] = -2.08627296E+01;
    COFD[2261] = 5.00799007E+00;
    COFD[2262] = -4.21303846E-01;
    COFD[2263] = 1.76831398E-02;
    COFD[2264] = -2.08614723E+01;
    COFD[2265] = 5.00799007E+00;
    COFD[2266] = -4.21303846E-01;
    COFD[2267] = 1.76831398E-02;
    COFD[2268] = -2.14847798E+01;
    COFD[2269] = 5.19080403E+00;
    COFD[2270] = -4.40813300E-01;
    COFD[2271] = 1.83643944E-02;
    COFD[2272] = -2.18845663E+01;
    COFD[2273] = 5.25324658E+00;
    COFD[2274] = -4.47318036E-01;
    COFD[2275] = 1.85841342E-02;
    COFD[2276] = -2.12169140E+01;
    COFD[2277] = 5.11479735E+00;
    COFD[2278] = -4.33050578E-01;
    COFD[2279] = 1.81098360E-02;
    COFD[2280] = -2.16082141E+01;
    COFD[2281] = 5.21099602E+00;
    COFD[2282] = -4.42524031E-01;
    COFD[2283] = 1.84030593E-02;
    COFD[2284] = -1.48205460E+01;
    COFD[2285] = 3.16773813E+00;
    COFD[2286] = -1.96823168E-01;
    COFD[2287] = 8.53058072E-03;
    COFD[2288] = -1.41192653E+01;
    COFD[2289] = 3.59720720E+00;
    COFD[2290] = -2.50852412E-01;
    COFD[2291] = 1.07992653E-02;
    COFD[2292] = -1.12831810E+01;
    COFD[2293] = 2.30017030E+00;
    COFD[2294] = -8.34928497E-02;
    COFD[2295] = 3.58611589E-03;
    COFD[2296] = -1.38902466E+01;
    COFD[2297] = 2.98298723E+00;
    COFD[2298] = -1.73428215E-01;
    COFD[2299] = 7.54224684E-03;
    COFD[2300] = -1.51132274E+01;
    COFD[2301] = 3.27853241E+00;
    COFD[2302] = -2.10970532E-01;
    COFD[2303] = 9.13336403E-03;
    COFD[2304] = -1.39108177E+01;
    COFD[2305] = 2.98298723E+00;
    COFD[2306] = -1.73428215E-01;
    COFD[2307] = 7.54224684E-03;
    COFD[2308] = -2.01678431E+01;
    COFD[2309] = 5.20184446E+00;
    COFD[2310] = -4.41727448E-01;
    COFD[2311] = 1.83837722E-02;
    COFD[2312] = -1.48546674E+01;
    COFD[2313] = 3.17494843E+00;
    COFD[2314] = -1.97752410E-01;
    COFD[2315] = 8.57054106E-03;
    COFD[2316] = -1.77943584E+01;
    COFD[2317] = 4.18471443E+00;
    COFD[2318] = -3.22640353E-01;
    COFD[2319] = 1.37280239E-02;
    COFD[2320] = -1.58272050E+01;
    COFD[2321] = 3.58720055E+00;
    COFD[2322] = -2.49572308E-01;
    COFD[2323] = 1.07445741E-02;
    COFD[2324] = -1.57655847E+01;
    COFD[2325] = 3.56140952E+00;
    COFD[2326] = -2.46282892E-01;
    COFD[2327] = 1.06044757E-02;
    COFD[2328] = -1.51211577E+01;
    COFD[2329] = 3.27853241E+00;
    COFD[2330] = -2.10970532E-01;
    COFD[2331] = 9.13336403E-03;
    COFD[2332] = -1.51287346E+01;
    COFD[2333] = 3.27853241E+00;
    COFD[2334] = -2.10970532E-01;
    COFD[2335] = 9.13336403E-03;
    COFD[2336] = -2.02564926E+01;
    COFD[2337] = 5.09185655E+00;
    COFD[2338] = -4.30641631E-01;
    COFD[2339] = 1.80275784E-02;
    COFD[2340] = -1.79148174E+01;
    COFD[2341] = 4.20508193E+00;
    COFD[2342] = -3.25158699E-01;
    COFD[2343] = 1.38321131E-02;
    COFD[2344] = -1.75212436E+01;
    COFD[2345] = 4.15347096E+00;
    COFD[2346] = -3.18813734E-01;
    COFD[2347] = 1.35714974E-02;
    COFD[2348] = -1.79056966E+01;
    COFD[2349] = 4.20508193E+00;
    COFD[2350] = -3.25158699E-01;
    COFD[2351] = 1.38321131E-02;
    COFD[2352] = -1.79896360E+01;
    COFD[2353] = 4.31259444E+00;
    COFD[2354] = -3.38584827E-01;
    COFD[2355] = 1.43931810E-02;
    COFD[2356] = -1.80009126E+01;
    COFD[2357] = 4.31259444E+00;
    COFD[2358] = -3.38584827E-01;
    COFD[2359] = 1.43931810E-02;
    COFD[2360] = -2.01915368E+01;
    COFD[2361] = 5.06295779E+00;
    COFD[2362] = -4.27579258E-01;
    COFD[2363] = 1.79215546E-02;
    COFD[2364] = -1.99811407E+01;
    COFD[2365] = 4.93673368E+00;
    COFD[2366] = -4.12924135E-01;
    COFD[2367] = 1.73539126E-02;
    COFD[2368] = -1.58800406E+01;
    COFD[2369] = 3.64628858E+00;
    COFD[2370] = -2.57092933E-01;
    COFD[2371] = 1.10642427E-02;
    COFD[2372] = -1.99914803E+01;
    COFD[2373] = 4.93673368E+00;
    COFD[2374] = -4.12924135E-01;
    COFD[2375] = 1.73539126E-02;
    COFD[2376] = -1.89680638E+01;
    COFD[2377] = 4.56070745E+00;
    COFD[2378] = -3.68517396E-01;
    COFD[2379] = 1.56021102E-02;
    COFD[2380] = -1.89680638E+01;
    COFD[2381] = 4.56070745E+00;
    COFD[2382] = -3.68517396E-01;
    COFD[2383] = 1.56021102E-02;
    COFD[2384] = -1.87724574E+01;
    COFD[2385] = 4.50343635E+00;
    COFD[2386] = -3.61818471E-01;
    COFD[2387] = 1.53404608E-02;
    COFD[2388] = -1.95470492E+01;
    COFD[2389] = 4.70710890E+00;
    COFD[2390] = -3.86453789E-01;
    COFD[2391] = 1.63380393E-02;
    COFD[2392] = -1.95100699E+01;
    COFD[2393] = 4.69879097E+00;
    COFD[2394] = -3.85441745E-01;
    COFD[2395] = 1.62968073E-02;
    COFD[2396] = -1.95135102E+01;
    COFD[2397] = 4.69879097E+00;
    COFD[2398] = -3.85441745E-01;
    COFD[2399] = 1.62968073E-02;
    COFD[2400] = -1.96223921E+01;
    COFD[2401] = 4.68635682E+00;
    COFD[2402] = -3.83931845E-01;
    COFD[2403] = 1.62354222E-02;
    COFD[2404] = -2.02656235E+01;
    COFD[2405] = 4.93517265E+00;
    COFD[2406] = -4.12740828E-01;
    COFD[2407] = 1.73467245E-02;
    COFD[2408] = -2.02179674E+01;
    COFD[2409] = 4.93423363E+00;
    COFD[2410] = -4.12630529E-01;
    COFD[2411] = 1.73423980E-02;
    COFD[2412] = -2.01357044E+01;
    COFD[2413] = 4.83483781E+00;
    COFD[2414] = -4.01422871E-01;
    COFD[2415] = 1.69235370E-02;
    COFD[2416] = -2.00535121E+01;
    COFD[2417] = 4.80382495E+00;
    COFD[2418] = -3.97838890E-01;
    COFD[2419] = 1.67856063E-02;
    COFD[2420] = -2.03940771E+01;
    COFD[2421] = 4.95138566E+00;
    COFD[2422] = -4.14645044E-01;
    COFD[2423] = 1.74214162E-02;
    COFD[2424] = -2.02811760E+01;
    COFD[2425] = 5.06295779E+00;
    COFD[2426] = -4.27579258E-01;
    COFD[2427] = 1.79215546E-02;
    COFD[2428] = -2.02855982E+01;
    COFD[2429] = 5.06295779E+00;
    COFD[2430] = -4.27579258E-01;
    COFD[2431] = 1.79215546E-02;
    COFD[2432] = -2.01413374E+01;
    COFD[2433] = 4.95889869E+00;
    COFD[2434] = -4.15526934E-01;
    COFD[2435] = 1.74559870E-02;
    COFD[2436] = -2.02231901E+01;
    COFD[2437] = 4.91354418E+00;
    COFD[2438] = -4.10191969E-01;
    COFD[2439] = 1.72463969E-02;
    COFD[2440] = -1.88680678E+01;
    COFD[2441] = 4.53065425E+00;
    COFD[2442] = -3.64983473E-01;
    COFD[2443] = 1.54632547E-02;
    COFD[2444] = -1.89621737E+01;
    COFD[2445] = 4.56070745E+00;
    COFD[2446] = -3.68517396E-01;
    COFD[2447] = 1.56021102E-02;
    COFD[2448] = -2.05907805E+01;
    COFD[2449] = 5.02104665E+00;
    COFD[2450] = -4.22831627E-01;
    COFD[2451] = 1.77428177E-02;
    COFD[2452] = -2.05958878E+01;
    COFD[2453] = 4.99149256E+00;
    COFD[2454] = -4.19363176E-01;
    COFD[2455] = 1.76068685E-02;
    COFD[2456] = -2.01447821E+01;
    COFD[2457] = 4.95889869E+00;
    COFD[2458] = -4.15526934E-01;
    COFD[2459] = 1.74559870E-02;
    COFD[2460] = -2.08797845E+01;
    COFD[2461] = 5.05404465E+00;
    COFD[2462] = -4.26628831E-01;
    COFD[2463] = 1.78883654E-02;
    COFD[2464] = -2.08624149E+01;
    COFD[2465] = 5.07617318E+00;
    COFD[2466] = -4.28983517E-01;
    COFD[2467] = 1.79703565E-02;
    COFD[2468] = -2.08740011E+01;
    COFD[2469] = 5.00799007E+00;
    COFD[2470] = -4.21303846E-01;
    COFD[2471] = 1.76831398E-02;
    COFD[2472] = -2.08727152E+01;
    COFD[2473] = 5.00799007E+00;
    COFD[2474] = -4.21303846E-01;
    COFD[2475] = 1.76831398E-02;
    COFD[2476] = -2.14967605E+01;
    COFD[2477] = 5.19080403E+00;
    COFD[2478] = -4.40813300E-01;
    COFD[2479] = 1.83643944E-02;
    COFD[2480] = -2.18970448E+01;
    COFD[2481] = 5.25324658E+00;
    COFD[2482] = -4.47318036E-01;
    COFD[2483] = 1.85841342E-02;
    COFD[2484] = -2.12285424E+01;
    COFD[2485] = 5.11479735E+00;
    COFD[2486] = -4.33050578E-01;
    COFD[2487] = 1.81098360E-02;
    COFD[2488] = -2.16204501E+01;
    COFD[2489] = 5.21099602E+00;
    COFD[2490] = -4.42524031E-01;
    COFD[2491] = 1.84030593E-02;
    COFD[2492] = -1.48273948E+01;
    COFD[2493] = 3.16773813E+00;
    COFD[2494] = -1.96823168E-01;
    COFD[2495] = 8.53058072E-03;
    COFD[2496] = -1.93172760E+01;
    COFD[2497] = 5.37588635E+00;
    COFD[2498] = -4.60762834E-01;
    COFD[2499] = 1.90722493E-02;
    COFD[2500] = -1.54361025E+01;
    COFD[2501] = 3.83529069E+00;
    COFD[2502] = -2.80481778E-01;
    COFD[2503] = 1.20302770E-02;
    COFD[2504] = -1.88359692E+01;
    COFD[2505] = 4.75778339E+00;
    COFD[2506] = -3.92510587E-01;
    COFD[2507] = 1.65802295E-02;
    COFD[2508] = -2.02419373E+01;
    COFD[2509] = 5.09185655E+00;
    COFD[2510] = -4.30641631E-01;
    COFD[2511] = 1.80275784E-02;
    COFD[2512] = -1.88556852E+01;
    COFD[2513] = 4.75778339E+00;
    COFD[2514] = -3.92510587E-01;
    COFD[2515] = 1.65802295E-02;
    COFD[2516] = -1.94329144E+01;
    COFD[2517] = 4.32601172E+00;
    COFD[2518] = -2.67841680E-01;
    COFD[2519] = 8.59623460E-03;
    COFD[2520] = -1.99660141E+01;
    COFD[2521] = 4.99411384E+00;
    COFD[2522] = -4.19672196E-01;
    COFD[2523] = 1.76190456E-02;
    COFD[2524] = -2.19981784E+01;
    COFD[2525] = 5.48342675E+00;
    COFD[2526] = -4.57782070E-01;
    COFD[2527] = 1.82957850E-02;
    COFD[2528] = -2.10007851E+01;
    COFD[2529] = 5.36906158E+00;
    COFD[2530] = -4.60039482E-01;
    COFD[2531] = 1.90471208E-02;
    COFD[2532] = -2.09538981E+01;
    COFD[2533] = 5.35250010E+00;
    COFD[2534] = -4.58316703E-01;
    COFD[2535] = 1.89887690E-02;
    COFD[2536] = -2.02493843E+01;
    COFD[2537] = 5.09185655E+00;
    COFD[2538] = -4.30641631E-01;
    COFD[2539] = 1.80275784E-02;
    COFD[2540] = -2.02564926E+01;
    COFD[2541] = 5.09185655E+00;
    COFD[2542] = -4.30641631E-01;
    COFD[2543] = 1.80275784E-02;
    COFD[2544] = -2.05279942E+01;
    COFD[2545] = 4.65341631E+00;
    COFD[2546] = -3.17798333E-01;
    COFD[2547] = 1.10326172E-02;
    COFD[2548] = -2.20751052E+01;
    COFD[2549] = 5.48219421E+00;
    COFD[2550] = -4.57118011E-01;
    COFD[2551] = 1.82481275E-02;
    COFD[2552] = -2.18067905E+01;
    COFD[2553] = 5.48629379E+00;
    COFD[2554] = -4.58998623E-01;
    COFD[2555] = 1.83802900E-02;
    COFD[2556] = -2.20665149E+01;
    COFD[2557] = 5.48219421E+00;
    COFD[2558] = -4.57118011E-01;
    COFD[2559] = 1.82481275E-02;
    COFD[2560] = -2.18872558E+01;
    COFD[2561] = 5.46783880E+00;
    COFD[2562] = -4.52790454E-01;
    COFD[2563] = 1.79689550E-02;
    COFD[2564] = -2.18979249E+01;
    COFD[2565] = 5.46783880E+00;
    COFD[2566] = -4.52790454E-01;
    COFD[2567] = 1.79689550E-02;
    COFD[2568] = -2.07135807E+01;
    COFD[2569] = 4.73262810E+00;
    COFD[2570] = -3.29863658E-01;
    COFD[2571] = 1.16208974E-02;
    COFD[2572] = -2.12524252E+01;
    COFD[2573] = 4.93163464E+00;
    COFD[2574] = -3.60831131E-01;
    COFD[2575] = 1.31516943E-02;
    COFD[2576] = -2.10061234E+01;
    COFD[2577] = 5.40483692E+00;
    COFD[2578] = -4.63658873E-01;
    COFD[2579] = 1.91646114E-02;
    COFD[2580] = -2.12620508E+01;
    COFD[2581] = 4.93163464E+00;
    COFD[2582] = -3.60831131E-01;
    COFD[2583] = 1.31516943E-02;
    COFD[2584] = -2.20218901E+01;
    COFD[2585] = 5.33662478E+00;
    COFD[2586] = -4.27371229E-01;
    COFD[2587] = 1.65539339E-02;
    COFD[2588] = -2.20218901E+01;
    COFD[2589] = 5.33662478E+00;
    COFD[2590] = -4.27371229E-01;
    COFD[2591] = 1.65539339E-02;
    COFD[2592] = -2.20651382E+01;
    COFD[2593] = 5.38704575E+00;
    COFD[2594] = -4.36283370E-01;
    COFD[2595] = 1.70286543E-02;
    COFD[2596] = -2.20279970E+01;
    COFD[2597] = 5.23448541E+00;
    COFD[2598] = -4.09927386E-01;
    COFD[2599] = 1.56423231E-02;
    COFD[2600] = -2.20252840E+01;
    COFD[2601] = 5.24121166E+00;
    COFD[2602] = -4.11059025E-01;
    COFD[2603] = 1.57010104E-02;
    COFD[2604] = -2.20284624E+01;
    COFD[2605] = 5.24121166E+00;
    COFD[2606] = -4.11059025E-01;
    COFD[2607] = 1.57010104E-02;
    COFD[2608] = -2.21864158E+01;
    COFD[2609] = 5.25139402E+00;
    COFD[2610] = -4.12771525E-01;
    COFD[2611] = 1.57897958E-02;
    COFD[2612] = -2.15372889E+01;
    COFD[2613] = 4.93380640E+00;
    COFD[2614] = -3.61171699E-01;
    COFD[2615] = 1.31686212E-02;
    COFD[2616] = -2.14953513E+01;
    COFD[2617] = 4.93510998E+00;
    COFD[2618] = -3.61376137E-01;
    COFD[2619] = 1.31787827E-02;
    COFD[2620] = -2.20062771E+01;
    COFD[2621] = 5.09608960E+00;
    COFD[2622] = -3.86955479E-01;
    COFD[2623] = 1.44609083E-02;
    COFD[2624] = -2.20820648E+01;
    COFD[2625] = 5.13480519E+00;
    COFD[2626] = -3.93316744E-01;
    COFD[2627] = 1.47863438E-02;
    COFD[2628] = -2.15725574E+01;
    COFD[2629] = 4.91148766E+00;
    COFD[2630] = -3.57668763E-01;
    COFD[2631] = 1.29944614E-02;
    COFD[2632] = -2.07973123E+01;
    COFD[2633] = 4.73262810E+00;
    COFD[2634] = -3.29863658E-01;
    COFD[2635] = 1.16208974E-02;
    COFD[2636] = -2.08014167E+01;
    COFD[2637] = 4.73262810E+00;
    COFD[2638] = -3.29863658E-01;
    COFD[2639] = 1.16208974E-02;
    COFD[2640] = -2.12869434E+01;
    COFD[2641] = 4.90168525E+00;
    COFD[2642] = -3.56120032E-01;
    COFD[2643] = 1.29171867E-02;
    COFD[2644] = -2.16143949E+01;
    COFD[2645] = 4.96367791E+00;
    COFD[2646] = -3.65848555E-01;
    COFD[2647] = 1.34008603E-02;
    COFD[2648] = -2.20470531E+01;
    COFD[2649] = 5.36287304E+00;
    COFD[2650] = -4.31985459E-01;
    COFD[2651] = 1.67989944E-02;
    COFD[2652] = -2.20163926E+01;
    COFD[2653] = 5.33662478E+00;
    COFD[2654] = -4.27371229E-01;
    COFD[2655] = 1.65539339E-02;
    COFD[2656] = -2.13808065E+01;
    COFD[2657] = 4.81362806E+00;
    COFD[2658] = -3.42326990E-01;
    COFD[2659] = 1.22325012E-02;
    COFD[2660] = -2.15521873E+01;
    COFD[2661] = 4.85623963E+00;
    COFD[2662] = -3.48984100E-01;
    COFD[2663] = 1.25623211E-02;
    COFD[2664] = -2.12901259E+01;
    COFD[2665] = 4.90168525E+00;
    COFD[2666] = -3.56120032E-01;
    COFD[2667] = 1.29171867E-02;
    COFD[2668] = -2.14531444E+01;
    COFD[2669] = 4.75554226E+00;
    COFD[2670] = -3.33372867E-01;
    COFD[2671] = 1.17926183E-02;
    COFD[2672] = -2.12488262E+01;
    COFD[2673] = 4.69672574E+00;
    COFD[2674] = -3.24393524E-01;
    COFD[2675] = 1.13541249E-02;
    COFD[2676] = -2.17305043E+01;
    COFD[2677] = 4.83264094E+00;
    COFD[2678] = -3.45294821E-01;
    COFD[2679] = 1.23794137E-02;
    COFD[2680] = -2.17293342E+01;
    COFD[2681] = 4.83264094E+00;
    COFD[2682] = -3.45294821E-01;
    COFD[2683] = 1.23794137E-02;
    COFD[2684] = -2.08674131E+01;
    COFD[2685] = 4.37860613E+00;
    COFD[2686] = -2.75749639E-01;
    COFD[2687] = 8.97792922E-03;
    COFD[2688] = -2.06803913E+01;
    COFD[2689] = 4.19360028E+00;
    COFD[2690] = -2.47975502E-01;
    COFD[2691] = 7.63834548E-03;
    COFD[2692] = -2.12861576E+01;
    COFD[2693] = 4.59638767E+00;
    COFD[2694] = -3.08980273E-01;
    COFD[2695] = 1.05986522E-02;
    COFD[2696] = -2.07241021E+01;
    COFD[2697] = 4.28574979E+00;
    COFD[2698] = -2.61781425E-01;
    COFD[2699] = 8.30349840E-03;
    COFD[2700] = -1.99350648E+01;
    COFD[2701] = 4.98617709E+00;
    COFD[2702] = -4.18736053E-01;
    COFD[2703] = 1.75821344E-02;
    COFD[2704] = -1.72388534E+01;
    COFD[2705] = 4.60447813E+00;
    COFD[2706] = -3.73889368E-01;
    COFD[2707] = 1.58229235E-02;
    COFD[2708] = -1.35410895E+01;
    COFD[2709] = 3.05984388E+00;
    COFD[2710] = -1.83114865E-01;
    COFD[2711] = 7.94950068E-03;
    COFD[2712] = -1.65809526E+01;
    COFD[2713] = 3.87564052E+00;
    COFD[2714] = -2.85285750E-01;
    COFD[2715] = 1.22208235E-02;
    COFD[2716] = -1.79002510E+01;
    COFD[2717] = 4.20508193E+00;
    COFD[2718] = -3.25158699E-01;
    COFD[2719] = 1.38321131E-02;
    COFD[2720] = -1.66006787E+01;
    COFD[2721] = 3.87564052E+00;
    COFD[2722] = -2.85285750E-01;
    COFD[2723] = 1.22208235E-02;
    COFD[2724] = -2.16940491E+01;
    COFD[2725] = 5.44394964E+00;
    COFD[2726] = -4.47152201E-01;
    COFD[2727] = 1.76308159E-02;
    COFD[2728] = -1.75803754E+01;
    COFD[2729] = 4.08497111E+00;
    COFD[2730] = -3.10470821E-01;
    COFD[2731] = 1.32322450E-02;
    COFD[2732] = -2.09219827E+01;
    COFD[2733] = 5.19204548E+00;
    COFD[2734] = -4.40915128E-01;
    COFD[2735] = 1.83664930E-02;
    COFD[2736] = -1.88552742E+01;
    COFD[2737] = 4.59285496E+00;
    COFD[2738] = -3.72450186E-01;
    COFD[2739] = 1.57632217E-02;
    COFD[2740] = -1.87841274E+01;
    COFD[2741] = 4.56428962E+00;
    COFD[2742] = -3.68938951E-01;
    COFD[2743] = 1.56186844E-02;
    COFD[2744] = -1.79077036E+01;
    COFD[2745] = 4.20508193E+00;
    COFD[2746] = -3.25158699E-01;
    COFD[2747] = 1.38321131E-02;
    COFD[2748] = -1.79148174E+01;
    COFD[2749] = 4.20508193E+00;
    COFD[2750] = -3.25158699E-01;
    COFD[2751] = 1.38321131E-02;
    COFD[2752] = -2.20751052E+01;
    COFD[2753] = 5.48219421E+00;
    COFD[2754] = -4.57118011E-01;
    COFD[2755] = 1.82481275E-02;
    COFD[2756] = -2.10046868E+01;
    COFD[2757] = 5.19982587E+00;
    COFD[2758] = -4.41557143E-01;
    COFD[2759] = 1.83799644E-02;
    COFD[2760] = -2.07012761E+01;
    COFD[2761] = 5.17836635E+00;
    COFD[2762] = -4.39728088E-01;
    COFD[2763] = 1.83380035E-02;
    COFD[2764] = -2.09960904E+01;
    COFD[2765] = 5.19982587E+00;
    COFD[2766] = -4.41557143E-01;
    COFD[2767] = 1.83799644E-02;
    COFD[2768] = -2.10405011E+01;
    COFD[2769] = 5.27526035E+00;
    COFD[2770] = -4.49795919E-01;
    COFD[2771] = 1.86768812E-02;
    COFD[2772] = -2.10511773E+01;
    COFD[2773] = 5.27526035E+00;
    COFD[2774] = -4.49795919E-01;
    COFD[2775] = 1.86768812E-02;
    COFD[2776] = -2.20793238E+01;
    COFD[2777] = 5.48560730E+00;
    COFD[2778] = -4.58746417E-01;
    COFD[2779] = 1.83631968E-02;
    COFD[2780] = -2.21540958E+01;
    COFD[2781] = 5.49388984E+00;
    COFD[2782] = -4.63572253E-01;
    COFD[2783] = 1.87152221E-02;
    COFD[2784] = -1.89659293E+01;
    COFD[2785] = 4.65931164E+00;
    COFD[2786] = -3.80642811E-01;
    COFD[2787] = 1.61015362E-02;
    COFD[2788] = -2.21637296E+01;
    COFD[2789] = 5.49388984E+00;
    COFD[2790] = -4.63572253E-01;
    COFD[2791] = 1.87152221E-02;
    COFD[2792] = -2.18363989E+01;
    COFD[2793] = 5.43667692E+00;
    COFD[2794] = -4.65640530E-01;
    COFD[2795] = 1.91675299E-02;
    COFD[2796] = -2.18363989E+01;
    COFD[2797] = 5.43667692E+00;
    COFD[2798] = -4.65640530E-01;
    COFD[2799] = 1.91675299E-02;
    COFD[2800] = -2.17225130E+01;
    COFD[2801] = 5.41757449E+00;
    COFD[2802] = -4.64778880E-01;
    COFD[2803] = 1.91926324E-02;
    COFD[2804] = -2.21190597E+01;
    COFD[2805] = 5.45729646E+00;
    COFD[2806] = -4.65340425E-01;
    COFD[2807] = 1.90398312E-02;
    COFD[2808] = -2.20987680E+01;
    COFD[2809] = 5.45582845E+00;
    COFD[2810] = -4.65322899E-01;
    COFD[2811] = 1.90457166E-02;
    COFD[2812] = -2.21019494E+01;
    COFD[2813] = 5.45582845E+00;
    COFD[2814] = -4.65322899E-01;
    COFD[2815] = 1.90457166E-02;
    COFD[2816] = -2.22197724E+01;
    COFD[2817] = 5.45358272E+00;
    COFD[2818] = -4.65292061E-01;
    COFD[2819] = 1.90543847E-02;
    COFD[2820] = -2.24143090E+01;
    COFD[2821] = 5.49379864E+00;
    COFD[2822] = -4.63600612E-01;
    COFD[2823] = 1.87179811E-02;
    COFD[2824] = -2.23731099E+01;
    COFD[2825] = 5.49374441E+00;
    COFD[2826] = -4.63617711E-01;
    COFD[2827] = 1.87196409E-02;
    COFD[2828] = -2.24306557E+01;
    COFD[2829] = 5.46873098E+00;
    COFD[2830] = -4.63395699E-01;
    COFD[2831] = 1.88249267E-02;
    COFD[2832] = -2.24132666E+01;
    COFD[2833] = 5.46734823E+00;
    COFD[2834] = -4.64115641E-01;
    COFD[2835] = 1.88903818E-02;
    COFD[2836] = -2.24986477E+01;
    COFD[2837] = 5.49353251E+00;
    COFD[2838] = -4.63127751E-01;
    COFD[2839] = 1.86807325E-02;
    COFD[2840] = -2.21631238E+01;
    COFD[2841] = 5.48560730E+00;
    COFD[2842] = -4.58746417E-01;
    COFD[2843] = 1.83631968E-02;
    COFD[2844] = -2.21672319E+01;
    COFD[2845] = 5.48560730E+00;
    COFD[2846] = -4.58746417E-01;
    COFD[2847] = 1.83631968E-02;
    COFD[2848] = -2.22550223E+01;
    COFD[2849] = 5.49328085E+00;
    COFD[2850] = -4.62887295E-01;
    COFD[2851] = 1.86623580E-02;
    COFD[2852] = -2.24142003E+01;
    COFD[2853] = 5.49189132E+00;
    COFD[2854] = -4.63888808E-01;
    COFD[2855] = 1.87508465E-02;
    COFD[2856] = -2.17838478E+01;
    COFD[2857] = 5.42865606E+00;
    COFD[2858] = -4.65447368E-01;
    COFD[2859] = 1.91919264E-02;
    COFD[2860] = -2.18308969E+01;
    COFD[2861] = 5.43667692E+00;
    COFD[2862] = -4.65640530E-01;
    COFD[2863] = 1.91675299E-02;
    COFD[2864] = -2.25355062E+01;
    COFD[2865] = 5.49116738E+00;
    COFD[2866] = -4.60872204E-01;
    COFD[2867] = 1.85085987E-02;
    COFD[2868] = -2.25966091E+01;
    COFD[2869] = 5.49074590E+00;
    COFD[2870] = -4.61641528E-01;
    COFD[2871] = 1.85733572E-02;
    COFD[2872] = -2.22582078E+01;
    COFD[2873] = 5.49328085E+00;
    COFD[2874] = -4.62887295E-01;
    COFD[2875] = 1.86623580E-02;
    COFD[2876] = -2.27245540E+01;
    COFD[2877] = 5.48717471E+00;
    COFD[2878] = -4.59324333E-01;
    COFD[2879] = 1.84024197E-02;
    COFD[2880] = -2.26588661E+01;
    COFD[2881] = 5.48378580E+00;
    COFD[2882] = -4.57964754E-01;
    COFD[2883] = 1.83088022E-02;
    COFD[2884] = -2.28191084E+01;
    COFD[2885] = 5.49118368E+00;
    COFD[2886] = -4.61244560E-01;
    COFD[2887] = 1.85388488E-02;
    COFD[2888] = -2.28179371E+01;
    COFD[2889] = 5.49118368E+00;
    COFD[2890] = -4.61244560E-01;
    COFD[2891] = 1.85388488E-02;
    COFD[2892] = -2.29332780E+01;
    COFD[2893] = 5.44974657E+00;
    COFD[2894] = -4.48651670E-01;
    COFD[2895] = 1.77235430E-02;
    COFD[2896] = -2.30861220E+01;
    COFD[2897] = 5.41058317E+00;
    COFD[2898] = -4.40692244E-01;
    COFD[2899] = 1.72706793E-02;
    COFD[2900] = -2.29039602E+01;
    COFD[2901] = 5.47759909E+00;
    COFD[2902] = -4.55570838E-01;
    COFD[2903] = 1.81456235E-02;
    COFD[2904] = -2.29719722E+01;
    COFD[2905] = 5.43353924E+00;
    COFD[2906] = -4.45111358E-01;
    COFD[2907] = 1.75163465E-02;
    COFD[2908] = -1.75535683E+01;
    COFD[2909] = 4.07809359E+00;
    COFD[2910] = -3.09646036E-01;
    COFD[2911] = 1.31992743E-02;
    COFD[2912] = -1.68005136E+01;
    COFD[2913] = 4.55072147E+00;
    COFD[2914] = -3.67343400E-01;
    COFD[2915] = 1.55560012E-02;
    COFD[2916] = -1.32011826E+01;
    COFD[2917] = 3.02784198E+00;
    COFD[2918] = -1.79087955E-01;
    COFD[2919] = 7.78048913E-03;
    COFD[2920] = -1.62171827E+01;
    COFD[2921] = 3.84275248E+00;
    COFD[2922] = -2.81379017E-01;
    COFD[2923] = 1.20662349E-02;
    COFD[2924] = -1.75072049E+01;
    COFD[2925] = 4.15347096E+00;
    COFD[2926] = -3.18813734E-01;
    COFD[2927] = 1.35714974E-02;
    COFD[2928] = -1.62364191E+01;
    COFD[2929] = 3.84275248E+00;
    COFD[2930] = -2.81379017E-01;
    COFD[2931] = 1.20662349E-02;
    COFD[2932] = -2.14160385E+01;
    COFD[2933] = 5.45640844E+00;
    COFD[2934] = -4.50149956E-01;
    COFD[2935] = 1.78119073E-02;
    COFD[2936] = -1.72146872E+01;
    COFD[2937] = 4.04198813E+00;
    COFD[2938] = -3.05338303E-01;
    COFD[2939] = 1.30280077E-02;
    COFD[2940] = -2.05980858E+01;
    COFD[2941] = 5.16997013E+00;
    COFD[2942] = -4.38998744E-01;
    COFD[2943] = 1.83204571E-02;
    COFD[2944] = -1.84813514E+01;
    COFD[2945] = 4.54318703E+00;
    COFD[2946] = -3.66458643E-01;
    COFD[2947] = 1.55212948E-02;
    COFD[2948] = -1.84296831E+01;
    COFD[2949] = 4.52338181E+00;
    COFD[2950] = -3.64129854E-01;
    COFD[2951] = 1.54297700E-02;
    COFD[2952] = -1.75143894E+01;
    COFD[2953] = 4.15347096E+00;
    COFD[2954] = -3.18813734E-01;
    COFD[2955] = 1.35714974E-02;
    COFD[2956] = -1.75212436E+01;
    COFD[2957] = 4.15347096E+00;
    COFD[2958] = -3.18813734E-01;
    COFD[2959] = 1.35714974E-02;
    COFD[2960] = -2.18067905E+01;
    COFD[2961] = 5.48629379E+00;
    COFD[2962] = -4.58998623E-01;
    COFD[2963] = 1.83802900E-02;
    COFD[2964] = -2.07012761E+01;
    COFD[2965] = 5.17836635E+00;
    COFD[2966] = -4.39728088E-01;
    COFD[2967] = 1.83380035E-02;
    COFD[2968] = -2.03551591E+01;
    COFD[2969] = 5.14854026E+00;
    COFD[2970] = -4.36716588E-01;
    COFD[2971] = 1.82409313E-02;
    COFD[2972] = -2.06929753E+01;
    COFD[2973] = 5.17836635E+00;
    COFD[2974] = -4.39728088E-01;
    COFD[2975] = 1.83380035E-02;
    COFD[2976] = -2.06594639E+01;
    COFD[2977] = 5.22838636E+00;
    COFD[2978] = -4.44511745E-01;
    COFD[2979] = 1.84787582E-02;
    COFD[2980] = -2.06697990E+01;
    COFD[2981] = 5.22838636E+00;
    COFD[2982] = -4.44511745E-01;
    COFD[2983] = 1.84787582E-02;
    COFD[2984] = -2.18142322E+01;
    COFD[2985] = 5.49121786E+00;
    COFD[2986] = -4.60888505E-01;
    COFD[2987] = 1.85096783E-02;
    COFD[2988] = -2.18686233E+01;
    COFD[2989] = 5.48743241E+00;
    COFD[2990] = -4.64013540E-01;
    COFD[2991] = 1.87830225E-02;
    COFD[2992] = -1.85111018E+01;
    COFD[2993] = 4.59851812E+00;
    COFD[2994] = -3.73151692E-01;
    COFD[2995] = 1.57923348E-02;
    COFD[2996] = -2.18778661E+01;
    COFD[2997] = 5.48743241E+00;
    COFD[2998] = -4.64013540E-01;
    COFD[2999] = 1.87830225E-02;
    COFD[3000] = -2.15488001E+01;
    COFD[3001] = 5.42458950E+00;
    COFD[3002] = -4.65212987E-01;
    COFD[3003] = 1.91930766E-02;
    COFD[3004] = -2.15488001E+01;
    COFD[3005] = 5.42458950E+00;
    COFD[3006] = -4.65212987E-01;
    COFD[3007] = 1.91930766E-02;
    COFD[3008] = -2.14001824E+01;
    COFD[3009] = 5.39261033E+00;
    COFD[3010] = -4.62555501E-01;
    COFD[3011] = 1.91354360E-02;
    COFD[3012] = -2.18452486E+01;
    COFD[3013] = 5.44764923E+00;
    COFD[3014] = -4.65239944E-01;
    COFD[3015] = 1.90796594E-02;
    COFD[3016] = -2.18235959E+01;
    COFD[3017] = 5.44634823E+00;
    COFD[3018] = -4.65245448E-01;
    COFD[3019] = 1.90865967E-02;
    COFD[3020] = -2.18266355E+01;
    COFD[3021] = 5.44634823E+00;
    COFD[3022] = -4.65245448E-01;
    COFD[3023] = 1.90865967E-02;
    COFD[3024] = -2.19586673E+01;
    COFD[3025] = 5.44474614E+00;
    COFD[3026] = -4.65306860E-01;
    COFD[3027] = 1.90996430E-02;
    COFD[3028] = -2.21479968E+01;
    COFD[3029] = 5.48719572E+00;
    COFD[3030] = -4.64022066E-01;
    COFD[3031] = 1.87848905E-02;
    COFD[3032] = -2.21019004E+01;
    COFD[3033] = 5.48702704E+00;
    COFD[3034] = -4.64023035E-01;
    COFD[3035] = 1.87858024E-02;
    COFD[3036] = -2.21898994E+01;
    COFD[3037] = 5.46750298E+00;
    COFD[3038] = -4.64484282E-01;
    COFD[3039] = 1.89196260E-02;
    COFD[3040] = -2.21759649E+01;
    COFD[3041] = 5.46726363E+00;
    COFD[3042] = -4.65349411E-01;
    COFD[3043] = 1.89911285E-02;
    COFD[3044] = -2.22409758E+01;
    COFD[3045] = 5.48966728E+00;
    COFD[3046] = -4.63937028E-01;
    COFD[3047] = 1.87657460E-02;
    COFD[3048] = -2.18947817E+01;
    COFD[3049] = 5.49121786E+00;
    COFD[3050] = -4.60888505E-01;
    COFD[3051] = 1.85096783E-02;
    COFD[3052] = -2.18987166E+01;
    COFD[3053] = 5.49121786E+00;
    COFD[3054] = -4.60888505E-01;
    COFD[3055] = 1.85096783E-02;
    COFD[3056] = -2.19808137E+01;
    COFD[3057] = 5.49090680E+00;
    COFD[3058] = -4.63910285E-01;
    COFD[3059] = 1.87574539E-02;
    COFD[3060] = -2.21435675E+01;
    COFD[3061] = 5.48236053E+00;
    COFD[3062] = -4.63887012E-01;
    COFD[3063] = 1.87976815E-02;
    COFD[3064] = -2.14802351E+01;
    COFD[3065] = 5.41058359E+00;
    COFD[3066] = -4.64178328E-01;
    COFD[3067] = 1.91784294E-02;
    COFD[3068] = -2.15435139E+01;
    COFD[3069] = 5.42458950E+00;
    COFD[3070] = -4.65212987E-01;
    COFD[3071] = 1.91930766E-02;
    COFD[3072] = -2.22850608E+01;
    COFD[3073] = 5.49115033E+00;
    COFD[3074] = -4.62266752E-01;
    COFD[3075] = 1.86223135E-02;
    COFD[3076] = -2.23629227E+01;
    COFD[3077] = 5.49413323E+00;
    COFD[3078] = -4.63498940E-01;
    COFD[3079] = 1.87080519E-02;
    COFD[3080] = -2.19838573E+01;
    COFD[3081] = 5.49090680E+00;
    COFD[3082] = -4.63910285E-01;
    COFD[3083] = 1.87574539E-02;
    COFD[3084] = -2.24972209E+01;
    COFD[3085] = 5.49117520E+00;
    COFD[3086] = -4.61254571E-01;
    COFD[3087] = 1.85397061E-02;
    COFD[3088] = -2.24277170E+01;
    COFD[3089] = 5.48901287E+00;
    COFD[3090] = -4.60031399E-01;
    COFD[3091] = 1.84508313E-02;
    COFD[3092] = -2.25951559E+01;
    COFD[3093] = 5.49332573E+00;
    COFD[3094] = -4.62935000E-01;
    COFD[3095] = 1.86660283E-02;
    COFD[3096] = -2.25940458E+01;
    COFD[3097] = 5.49332573E+00;
    COFD[3098] = -4.62935000E-01;
    COFD[3099] = 1.86660283E-02;
    COFD[3100] = -2.27385006E+01;
    COFD[3101] = 5.46322829E+00;
    COFD[3102] = -4.51749940E-01;
    COFD[3103] = 1.79075310E-02;
    COFD[3104] = -2.29301063E+01;
    COFD[3105] = 5.43667678E+00;
    COFD[3106] = -4.45729487E-01;
    COFD[3107] = 1.75510667E-02;
    COFD[3108] = -2.26898752E+01;
    COFD[3109] = 5.48310931E+00;
    COFD[3110] = -4.57615590E-01;
    COFD[3111] = 1.82838764E-02;
    COFD[3112] = -2.27854551E+01;
    COFD[3113] = 5.45051999E+00;
    COFD[3114] = -4.48822771E-01;
    COFD[3115] = 1.77335774E-02;
    COFD[3116] = -1.71876409E+01;
    COFD[3117] = 4.03536932E+00;
    COFD[3118] = -3.04548174E-01;
    COFD[3119] = 1.29965591E-02;
    COFD[3120] = -1.72382912E+01;
    COFD[3121] = 4.60447813E+00;
    COFD[3122] = -3.73889368E-01;
    COFD[3123] = 1.58229235E-02;
    COFD[3124] = -1.35400012E+01;
    COFD[3125] = 3.05984388E+00;
    COFD[3126] = -1.83114865E-01;
    COFD[3127] = 7.94950068E-03;
    COFD[3128] = -1.65749660E+01;
    COFD[3129] = 3.87564052E+00;
    COFD[3130] = -2.85285750E-01;
    COFD[3131] = 1.22208235E-02;
    COFD[3132] = -1.78913897E+01;
    COFD[3133] = 4.20508193E+00;
    COFD[3134] = -3.25158699E-01;
    COFD[3135] = 1.38321131E-02;
    COFD[3136] = -1.65944527E+01;
    COFD[3137] = 3.87564052E+00;
    COFD[3138] = -2.85285750E-01;
    COFD[3139] = 1.22208235E-02;
    COFD[3140] = -2.16875939E+01;
    COFD[3141] = 5.44394964E+00;
    COFD[3142] = -4.47152201E-01;
    COFD[3143] = 1.76308159E-02;
    COFD[3144] = -1.75720813E+01;
    COFD[3145] = 4.08497111E+00;
    COFD[3146] = -3.10470821E-01;
    COFD[3147] = 1.32322450E-02;
    COFD[3148] = -2.09117850E+01;
    COFD[3149] = 5.19204548E+00;
    COFD[3150] = -4.40915128E-01;
    COFD[3151] = 1.83664930E-02;
    COFD[3152] = -1.88495268E+01;
    COFD[3153] = 4.59285496E+00;
    COFD[3154] = -3.72450186E-01;
    COFD[3155] = 1.57632217E-02;
    COFD[3156] = -1.87781303E+01;
    COFD[3157] = 4.56428962E+00;
    COFD[3158] = -3.68938951E-01;
    COFD[3159] = 1.56186844E-02;
    COFD[3160] = -1.78987105E+01;
    COFD[3161] = 4.20508193E+00;
    COFD[3162] = -3.25158699E-01;
    COFD[3163] = 1.38321131E-02;
    COFD[3164] = -1.79056966E+01;
    COFD[3165] = 4.20508193E+00;
    COFD[3166] = -3.25158699E-01;
    COFD[3167] = 1.38321131E-02;
    COFD[3168] = -2.20665149E+01;
    COFD[3169] = 5.48219421E+00;
    COFD[3170] = -4.57118011E-01;
    COFD[3171] = 1.82481275E-02;
    COFD[3172] = -2.09960904E+01;
    COFD[3173] = 5.19982587E+00;
    COFD[3174] = -4.41557143E-01;
    COFD[3175] = 1.83799644E-02;
    COFD[3176] = -2.06929753E+01;
    COFD[3177] = 5.17836635E+00;
    COFD[3178] = -4.39728088E-01;
    COFD[3179] = 1.83380035E-02;
    COFD[3180] = -2.09876392E+01;
    COFD[3181] = 5.19982587E+00;
    COFD[3182] = -4.41557143E-01;
    COFD[3183] = 1.83799644E-02;
    COFD[3184] = -2.10326848E+01;
    COFD[3185] = 5.27526035E+00;
    COFD[3186] = -4.49795919E-01;
    COFD[3187] = 1.86768812E-02;
    COFD[3188] = -2.10431936E+01;
    COFD[3189] = 5.27526035E+00;
    COFD[3190] = -4.49795919E-01;
    COFD[3191] = 1.86768812E-02;
    COFD[3192] = -2.20704568E+01;
    COFD[3193] = 5.48560730E+00;
    COFD[3194] = -4.58746417E-01;
    COFD[3195] = 1.83631968E-02;
    COFD[3196] = -2.21440868E+01;
    COFD[3197] = 5.49388984E+00;
    COFD[3198] = -4.63572253E-01;
    COFD[3199] = 1.87152221E-02;
    COFD[3200] = -1.89560208E+01;
    COFD[3201] = 4.65931164E+00;
    COFD[3202] = -3.80642811E-01;
    COFD[3203] = 1.61015362E-02;
    COFD[3204] = -2.21535278E+01;
    COFD[3205] = 5.49388984E+00;
    COFD[3206] = -4.63572253E-01;
    COFD[3207] = 1.87152221E-02;
    COFD[3208] = -2.18265893E+01;
    COFD[3209] = 5.43667692E+00;
    COFD[3210] = -4.65640530E-01;
    COFD[3211] = 1.91675299E-02;
    COFD[3212] = -2.18265893E+01;
    COFD[3213] = 5.43667692E+00;
    COFD[3214] = -4.65640530E-01;
    COFD[3215] = 1.91675299E-02;
    COFD[3216] = -2.17124997E+01;
    COFD[3217] = 5.41757449E+00;
    COFD[3218] = -4.64778880E-01;
    COFD[3219] = 1.91926324E-02;
    COFD[3220] = -2.21080365E+01;
    COFD[3221] = 5.45729646E+00;
    COFD[3222] = -4.65340425E-01;
    COFD[3223] = 1.90398312E-02;
    COFD[3224] = -2.20876730E+01;
    COFD[3225] = 5.45582845E+00;
    COFD[3226] = -4.65322899E-01;
    COFD[3227] = 1.90457166E-02;
    COFD[3228] = -2.20907844E+01;
    COFD[3229] = 5.45582845E+00;
    COFD[3230] = -4.65322899E-01;
    COFD[3231] = 1.90457166E-02;
    COFD[3232] = -2.22085390E+01;
    COFD[3233] = 5.45358272E+00;
    COFD[3234] = -4.65292061E-01;
    COFD[3235] = 1.90543847E-02;
    COFD[3236] = -2.24030784E+01;
    COFD[3237] = 5.49379864E+00;
    COFD[3238] = -4.63600612E-01;
    COFD[3239] = 1.87179811E-02;
    COFD[3240] = -2.23618125E+01;
    COFD[3241] = 5.49374441E+00;
    COFD[3242] = -4.63617711E-01;
    COFD[3243] = 1.87196409E-02;
    COFD[3244] = -2.24187146E+01;
    COFD[3245] = 5.46873098E+00;
    COFD[3246] = -4.63395699E-01;
    COFD[3247] = 1.88249267E-02;
    COFD[3248] = -2.24012739E+01;
    COFD[3249] = 5.46734823E+00;
    COFD[3250] = -4.64115641E-01;
    COFD[3251] = 1.88903818E-02;
    COFD[3252] = -2.24866045E+01;
    COFD[3253] = 5.49353251E+00;
    COFD[3254] = -4.63127751E-01;
    COFD[3255] = 1.86807325E-02;
    COFD[3256] = -2.21526557E+01;
    COFD[3257] = 5.48560730E+00;
    COFD[3258] = -4.58746417E-01;
    COFD[3259] = 1.83631968E-02;
    COFD[3260] = -2.21566783E+01;
    COFD[3261] = 5.48560730E+00;
    COFD[3262] = -4.58746417E-01;
    COFD[3263] = 1.83631968E-02;
    COFD[3264] = -2.22439303E+01;
    COFD[3265] = 5.49328085E+00;
    COFD[3266] = -4.62887295E-01;
    COFD[3267] = 1.86623580E-02;
    COFD[3268] = -2.24030382E+01;
    COFD[3269] = 5.49189132E+00;
    COFD[3270] = -4.63888808E-01;
    COFD[3271] = 1.87508465E-02;
    COFD[3272] = -2.17739349E+01;
    COFD[3273] = 5.42865606E+00;
    COFD[3274] = -4.65447368E-01;
    COFD[3275] = 1.91919264E-02;
    COFD[3276] = -2.18211936E+01;
    COFD[3277] = 5.43667692E+00;
    COFD[3278] = -4.65640530E-01;
    COFD[3279] = 1.91675299E-02;
    COFD[3280] = -2.25234155E+01;
    COFD[3281] = 5.49116738E+00;
    COFD[3282] = -4.60872204E-01;
    COFD[3283] = 1.85085987E-02;
    COFD[3284] = -2.25845680E+01;
    COFD[3285] = 5.49074590E+00;
    COFD[3286] = -4.61641528E-01;
    COFD[3287] = 1.85733572E-02;
    COFD[3288] = -2.22470457E+01;
    COFD[3289] = 5.49328085E+00;
    COFD[3290] = -4.62887295E-01;
    COFD[3291] = 1.86623580E-02;
    COFD[3292] = -2.27118623E+01;
    COFD[3293] = 5.48717471E+00;
    COFD[3294] = -4.59324333E-01;
    COFD[3295] = 1.84024197E-02;
    COFD[3296] = -2.26462126E+01;
    COFD[3297] = 5.48378580E+00;
    COFD[3298] = -4.57964754E-01;
    COFD[3299] = 1.83088022E-02;
    COFD[3300] = -2.28059444E+01;
    COFD[3301] = 5.49118368E+00;
    COFD[3302] = -4.61244560E-01;
    COFD[3303] = 1.85388488E-02;
    COFD[3304] = -2.28048034E+01;
    COFD[3305] = 5.49118368E+00;
    COFD[3306] = -4.61244560E-01;
    COFD[3307] = 1.85388488E-02;
    COFD[3308] = -2.29193653E+01;
    COFD[3309] = 5.44974657E+00;
    COFD[3310] = -4.48651670E-01;
    COFD[3311] = 1.77235430E-02;
    COFD[3312] = -2.30716889E+01;
    COFD[3313] = 5.41058317E+00;
    COFD[3314] = -4.40692244E-01;
    COFD[3315] = 1.72706793E-02;
    COFD[3316] = -2.28904184E+01;
    COFD[3317] = 5.47759909E+00;
    COFD[3318] = -4.55570838E-01;
    COFD[3319] = 1.81456235E-02;
    COFD[3320] = -2.29577921E+01;
    COFD[3321] = 5.43353924E+00;
    COFD[3322] = -4.45111358E-01;
    COFD[3323] = 1.75163465E-02;
    COFD[3324] = -1.75452737E+01;
    COFD[3325] = 4.07809359E+00;
    COFD[3326] = -3.09646036E-01;
    COFD[3327] = 1.31992743E-02;
    COFD[3328] = -1.73275531E+01;
    COFD[3329] = 4.71167685E+00;
    COFD[3330] = -3.87004255E-01;
    COFD[3331] = 1.63602477E-02;
    COFD[3332] = -1.35202960E+01;
    COFD[3333] = 3.12061076E+00;
    COFD[3334] = -1.90761802E-01;
    COFD[3335] = 8.27051223E-03;
    COFD[3336] = -1.65414090E+01;
    COFD[3337] = 3.93824541E+00;
    COFD[3338] = -2.92787590E-01;
    COFD[3339] = 1.25209211E-02;
    COFD[3340] = -1.79764577E+01;
    COFD[3341] = 4.31259444E+00;
    COFD[3342] = -3.38584827E-01;
    COFD[3343] = 1.43931810E-02;
    COFD[3344] = -1.65598218E+01;
    COFD[3345] = 3.93824541E+00;
    COFD[3346] = -2.92787590E-01;
    COFD[3347] = 1.25209211E-02;
    COFD[3348] = -2.14146533E+01;
    COFD[3349] = 5.39388015E+00;
    COFD[3350] = -4.37574572E-01;
    COFD[3351] = 1.70997728E-02;
    COFD[3352] = -1.76229592E+01;
    COFD[3353] = 4.17526194E+00;
    COFD[3354] = -3.21478247E-01;
    COFD[3355] = 1.36802836E-02;
    COFD[3356] = -2.09138237E+01;
    COFD[3357] = 5.25757923E+00;
    COFD[3358] = -4.47803407E-01;
    COFD[3359] = 1.86021883E-02;
    COFD[3360] = -1.89650398E+01;
    COFD[3361] = 4.70175320E+00;
    COFD[3362] = -3.85802299E-01;
    COFD[3363] = 1.63115029E-02;
    COFD[3364] = -1.88959377E+01;
    COFD[3365] = 4.67519664E+00;
    COFD[3366] = -3.82581158E-01;
    COFD[3367] = 1.61807120E-02;
    COFD[3368] = -1.79832047E+01;
    COFD[3369] = 4.31259444E+00;
    COFD[3370] = -3.38584827E-01;
    COFD[3371] = 1.43931810E-02;
    COFD[3372] = -1.79896360E+01;
    COFD[3373] = 4.31259444E+00;
    COFD[3374] = -3.38584827E-01;
    COFD[3375] = 1.43931810E-02;
    COFD[3376] = -2.18872558E+01;
    COFD[3377] = 5.46783880E+00;
    COFD[3378] = -4.52790454E-01;
    COFD[3379] = 1.79689550E-02;
    COFD[3380] = -2.10405011E+01;
    COFD[3381] = 5.27526035E+00;
    COFD[3382] = -4.49795919E-01;
    COFD[3383] = 1.86768812E-02;
    COFD[3384] = -2.06594639E+01;
    COFD[3385] = 5.22838636E+00;
    COFD[3386] = -4.44511745E-01;
    COFD[3387] = 1.84787582E-02;
    COFD[3388] = -2.10326848E+01;
    COFD[3389] = 5.27526035E+00;
    COFD[3390] = -4.49795919E-01;
    COFD[3391] = 1.86768812E-02;
    COFD[3392] = -2.10710613E+01;
    COFD[3393] = 5.35206716E+00;
    COFD[3394] = -4.58271471E-01;
    COFD[3395] = 1.89872276E-02;
    COFD[3396] = -2.10808340E+01;
    COFD[3397] = 5.35206716E+00;
    COFD[3398] = -4.58271471E-01;
    COFD[3399] = 1.89872276E-02;
    COFD[3400] = -2.19028963E+01;
    COFD[3401] = 5.47597765E+00;
    COFD[3402] = -4.55020526E-01;
    COFD[3403] = 1.81091205E-02;
    COFD[3404] = -2.19954815E+01;
    COFD[3405] = 5.49007189E+00;
    COFD[3406] = -4.60469654E-01;
    COFD[3407] = 1.84812405E-02;
    COFD[3408] = -1.89929565E+01;
    COFD[3409] = 4.75627153E+00;
    COFD[3410] = -3.92329566E-01;
    COFD[3411] = 1.65729779E-02;
    COFD[3412] = -2.20040945E+01;
    COFD[3413] = 5.49007189E+00;
    COFD[3414] = -4.60469654E-01;
    COFD[3415] = 1.84812405E-02;
    COFD[3416] = -2.17210686E+01;
    COFD[3417] = 5.45061574E+00;
    COFD[3418] = -4.65241893E-01;
    COFD[3419] = 1.90650429E-02;
    COFD[3420] = -2.17210686E+01;
    COFD[3421] = 5.45061574E+00;
    COFD[3422] = -4.65241893E-01;
    COFD[3423] = 1.90650429E-02;
    COFD[3424] = -2.16232094E+01;
    COFD[3425] = 5.44000785E+00;
    COFD[3426] = -4.65569459E-01;
    COFD[3427] = 1.91449572E-02;
    COFD[3428] = -2.19960701E+01;
    COFD[3429] = 5.46751578E+00;
    COFD[3430] = -4.64445737E-01;
    COFD[3431] = 1.89164296E-02;
    COFD[3432] = -2.19775435E+01;
    COFD[3433] = 5.46745929E+00;
    COFD[3434] = -4.64624358E-01;
    COFD[3435] = 1.89312280E-02;
    COFD[3436] = -2.19803572E+01;
    COFD[3437] = 5.46745929E+00;
    COFD[3438] = -4.64624358E-01;
    COFD[3439] = 1.89312280E-02;
    COFD[3440] = -2.21121170E+01;
    COFD[3441] = 5.46738366E+00;
    COFD[3442] = -4.64895439E-01;
    COFD[3443] = 1.89536362E-02;
    COFD[3444] = -2.22655327E+01;
    COFD[3445] = 5.49021632E+00;
    COFD[3446] = -4.60531670E-01;
    COFD[3447] = 1.84855708E-02;
    COFD[3448] = -2.22210765E+01;
    COFD[3449] = 5.49030211E+00;
    COFD[3450] = -4.60568793E-01;
    COFD[3451] = 1.84881663E-02;
    COFD[3452] = -2.23433381E+01;
    COFD[3453] = 5.49008914E+00;
    COFD[3454] = -4.63927959E-01;
    COFD[3455] = 1.87629280E-02;
    COFD[3456] = -2.23119605E+01;
    COFD[3457] = 5.48315923E+00;
    COFD[3458] = -4.63919102E-01;
    COFD[3459] = 1.87963630E-02;
    COFD[3460] = -2.23453522E+01;
    COFD[3461] = 5.48861596E+00;
    COFD[3462] = -4.59873623E-01;
    COFD[3463] = 1.84399621E-02;
    COFD[3464] = -2.19781845E+01;
    COFD[3465] = 5.47597765E+00;
    COFD[3466] = -4.55020526E-01;
    COFD[3467] = 1.81091205E-02;
    COFD[3468] = -2.19818415E+01;
    COFD[3469] = 5.47597765E+00;
    COFD[3470] = -4.55020526E-01;
    COFD[3471] = 1.81091205E-02;
    COFD[3472] = -2.20896596E+01;
    COFD[3473] = 5.48779691E+00;
    COFD[3474] = -4.59557311E-01;
    COFD[3475] = 1.84182884E-02;
    COFD[3476] = -2.22736645E+01;
    COFD[3477] = 5.49118605E+00;
    COFD[3478] = -4.61241736E-01;
    COFD[3479] = 1.85386070E-02;
    COFD[3480] = -2.16714153E+01;
    COFD[3481] = 5.44455741E+00;
    COFD[3482] = -4.65313340E-01;
    COFD[3483] = 1.91011181E-02;
    COFD[3484] = -2.17161314E+01;
    COFD[3485] = 5.45061574E+00;
    COFD[3486] = -4.65241893E-01;
    COFD[3487] = 1.90650429E-02;
    COFD[3488] = -2.23705009E+01;
    COFD[3489] = 5.48227325E+00;
    COFD[3490] = -4.57157805E-01;
    COFD[3491] = 1.82509581E-02;
    COFD[3492] = -2.24444141E+01;
    COFD[3493] = 5.48421194E+00;
    COFD[3494] = -4.58197253E-01;
    COFD[3495] = 1.83255177E-02;
    COFD[3496] = -2.20924770E+01;
    COFD[3497] = 5.48779691E+00;
    COFD[3498] = -4.59557311E-01;
    COFD[3499] = 1.84182884E-02;
    COFD[3500] = -2.25641242E+01;
    COFD[3501] = 5.47765924E+00;
    COFD[3502] = -4.55592027E-01;
    COFD[3503] = 1.81470399E-02;
    COFD[3504] = -2.24904947E+01;
    COFD[3505] = 5.47343682E+00;
    COFD[3506] = -4.54166709E-01;
    COFD[3507] = 1.80526774E-02;
    COFD[3508] = -2.26677760E+01;
    COFD[3509] = 5.48311530E+00;
    COFD[3510] = -4.57619182E-01;
    COFD[3511] = 1.82841374E-02;
    COFD[3512] = -2.26667616E+01;
    COFD[3513] = 5.48311530E+00;
    COFD[3514] = -4.57619182E-01;
    COFD[3515] = 1.82841374E-02;
    COFD[3516] = -2.27048922E+01;
    COFD[3517] = 5.40961679E+00;
    COFD[3518] = -4.40511249E-01;
    COFD[3519] = 1.72607396E-02;
    COFD[3520] = -2.28184844E+01;
    COFD[3521] = 5.35124813E+00;
    COFD[3522] = -4.29932425E-01;
    COFD[3523] = 1.66897210E-02;
    COFD[3524] = -2.27210131E+01;
    COFD[3525] = 5.45828076E+00;
    COFD[3526] = -4.50591566E-01;
    COFD[3527] = 1.78383531E-02;
    COFD[3528] = -2.27171217E+01;
    COFD[3529] = 5.38171422E+00;
    COFD[3530] = -4.35323468E-01;
    COFD[3531] = 1.69770246E-02;
    COFD[3532] = -1.75930902E+01;
    COFD[3533] = 4.16734978E+00;
    COFD[3534] = -3.20508152E-01;
    COFD[3535] = 1.36405589E-02;
    COFD[3536] = -1.73283030E+01;
    COFD[3537] = 4.71167685E+00;
    COFD[3538] = -3.87004255E-01;
    COFD[3539] = 1.63602477E-02;
    COFD[3540] = -1.35217408E+01;
    COFD[3541] = 3.12061076E+00;
    COFD[3542] = -1.90761802E-01;
    COFD[3543] = 8.27051223E-03;
    COFD[3544] = -1.65490142E+01;
    COFD[3545] = 3.93824541E+00;
    COFD[3546] = -2.92787590E-01;
    COFD[3547] = 1.25209211E-02;
    COFD[3548] = -1.79874377E+01;
    COFD[3549] = 4.31259444E+00;
    COFD[3550] = -3.38584827E-01;
    COFD[3551] = 1.43931810E-02;
    COFD[3552] = -1.65677145E+01;
    COFD[3553] = 3.93824541E+00;
    COFD[3554] = -2.92787590E-01;
    COFD[3555] = 1.25209211E-02;
    COFD[3556] = -2.14228203E+01;
    COFD[3557] = 5.39388015E+00;
    COFD[3558] = -4.37574572E-01;
    COFD[3559] = 1.70997728E-02;
    COFD[3560] = -1.76332867E+01;
    COFD[3561] = 4.17526194E+00;
    COFD[3562] = -3.21478247E-01;
    COFD[3563] = 1.36802836E-02;
    COFD[3564] = -2.09263167E+01;
    COFD[3565] = 5.25757923E+00;
    COFD[3566] = -4.47803407E-01;
    COFD[3567] = 1.86021883E-02;
    COFD[3568] = -1.89723565E+01;
    COFD[3569] = 4.70175320E+00;
    COFD[3570] = -3.85802299E-01;
    COFD[3571] = 1.63115029E-02;
    COFD[3572] = -1.89035556E+01;
    COFD[3573] = 4.67519664E+00;
    COFD[3574] = -3.82581158E-01;
    COFD[3575] = 1.61807120E-02;
    COFD[3576] = -1.79943355E+01;
    COFD[3577] = 4.31259444E+00;
    COFD[3578] = -3.38584827E-01;
    COFD[3579] = 1.43931810E-02;
    COFD[3580] = -1.80009126E+01;
    COFD[3581] = 4.31259444E+00;
    COFD[3582] = -3.38584827E-01;
    COFD[3583] = 1.43931810E-02;
    COFD[3584] = -2.18979249E+01;
    COFD[3585] = 5.46783880E+00;
    COFD[3586] = -4.52790454E-01;
    COFD[3587] = 1.79689550E-02;
    COFD[3588] = -2.10511773E+01;
    COFD[3589] = 5.27526035E+00;
    COFD[3590] = -4.49795919E-01;
    COFD[3591] = 1.86768812E-02;
    COFD[3592] = -2.06697990E+01;
    COFD[3593] = 5.22838636E+00;
    COFD[3594] = -4.44511745E-01;
    COFD[3595] = 1.84787582E-02;
    COFD[3596] = -2.10431936E+01;
    COFD[3597] = 5.27526035E+00;
    COFD[3598] = -4.49795919E-01;
    COFD[3599] = 1.86768812E-02;
    COFD[3600] = -2.10808340E+01;
    COFD[3601] = 5.35206716E+00;
    COFD[3602] = -4.58271471E-01;
    COFD[3603] = 1.89872276E-02;
    COFD[3604] = -2.10908015E+01;
    COFD[3605] = 5.35206716E+00;
    COFD[3606] = -4.58271471E-01;
    COFD[3607] = 1.89872276E-02;
    COFD[3608] = -2.19138830E+01;
    COFD[3609] = 5.47597765E+00;
    COFD[3610] = -4.55020526E-01;
    COFD[3611] = 1.81091205E-02;
    COFD[3612] = -2.20077630E+01;
    COFD[3613] = 5.49007189E+00;
    COFD[3614] = -4.60469654E-01;
    COFD[3615] = 1.84812405E-02;
    COFD[3616] = -1.90051251E+01;
    COFD[3617] = 4.75627153E+00;
    COFD[3618] = -3.92329566E-01;
    COFD[3619] = 1.65729779E-02;
    COFD[3620] = -2.20165921E+01;
    COFD[3621] = 5.49007189E+00;
    COFD[3622] = -4.60469654E-01;
    COFD[3623] = 1.84812405E-02;
    COFD[3624] = -2.17331258E+01;
    COFD[3625] = 5.45061574E+00;
    COFD[3626] = -4.65241893E-01;
    COFD[3627] = 1.90650429E-02;
    COFD[3628] = -2.17331258E+01;
    COFD[3629] = 5.45061574E+00;
    COFD[3630] = -4.65241893E-01;
    COFD[3631] = 1.90650429E-02;
    COFD[3632] = -2.16354957E+01;
    COFD[3633] = 5.44000785E+00;
    COFD[3634] = -4.65569459E-01;
    COFD[3635] = 1.91449572E-02;
    COFD[3636] = -2.20094808E+01;
    COFD[3637] = 5.46751578E+00;
    COFD[3638] = -4.64445737E-01;
    COFD[3639] = 1.89164296E-02;
    COFD[3640] = -2.19910334E+01;
    COFD[3641] = 5.46745929E+00;
    COFD[3642] = -4.64624358E-01;
    COFD[3643] = 1.89312280E-02;
    COFD[3644] = -2.19939242E+01;
    COFD[3645] = 5.46745929E+00;
    COFD[3646] = -4.64624358E-01;
    COFD[3647] = 1.89312280E-02;
    COFD[3648] = -2.21257593E+01;
    COFD[3649] = 5.46738366E+00;
    COFD[3650] = -4.64895439E-01;
    COFD[3651] = 1.89536362E-02;
    COFD[3652] = -2.22791718E+01;
    COFD[3653] = 5.49021632E+00;
    COFD[3654] = -4.60531670E-01;
    COFD[3655] = 1.84855708E-02;
    COFD[3656] = -2.22347892E+01;
    COFD[3657] = 5.49030211E+00;
    COFD[3658] = -4.60568793E-01;
    COFD[3659] = 1.84881663E-02;
    COFD[3660] = -2.23577541E+01;
    COFD[3661] = 5.49008914E+00;
    COFD[3662] = -4.63927959E-01;
    COFD[3663] = 1.87629280E-02;
    COFD[3664] = -2.23264327E+01;
    COFD[3665] = 5.48315923E+00;
    COFD[3666] = -4.63919102E-01;
    COFD[3667] = 1.87963630E-02;
    COFD[3668] = -2.23598793E+01;
    COFD[3669] = 5.48861596E+00;
    COFD[3670] = -4.59873623E-01;
    COFD[3671] = 1.84399621E-02;
    COFD[3672] = -2.19909795E+01;
    COFD[3673] = 5.47597765E+00;
    COFD[3674] = -4.55020526E-01;
    COFD[3675] = 1.81091205E-02;
    COFD[3676] = -2.19947318E+01;
    COFD[3677] = 5.47597765E+00;
    COFD[3678] = -4.55020526E-01;
    COFD[3679] = 1.81091205E-02;
    COFD[3680] = -2.21031461E+01;
    COFD[3681] = 5.48779691E+00;
    COFD[3682] = -4.59557311E-01;
    COFD[3683] = 1.84182884E-02;
    COFD[3684] = -2.22872282E+01;
    COFD[3685] = 5.49118605E+00;
    COFD[3686] = -4.61241736E-01;
    COFD[3687] = 1.85386070E-02;
    COFD[3688] = -2.16835888E+01;
    COFD[3689] = 5.44455741E+00;
    COFD[3690] = -4.65313340E-01;
    COFD[3691] = 1.91011181E-02;
    COFD[3692] = -2.17280687E+01;
    COFD[3693] = 5.45061574E+00;
    COFD[3694] = -4.65241893E-01;
    COFD[3695] = 1.90650429E-02;
    COFD[3696] = -2.23850795E+01;
    COFD[3697] = 5.48227325E+00;
    COFD[3698] = -4.57157805E-01;
    COFD[3699] = 1.82509581E-02;
    COFD[3700] = -2.24589388E+01;
    COFD[3701] = 5.48421194E+00;
    COFD[3702] = -4.58197253E-01;
    COFD[3703] = 1.83255177E-02;
    COFD[3704] = -2.21060408E+01;
    COFD[3705] = 5.48779691E+00;
    COFD[3706] = -4.59557311E-01;
    COFD[3707] = 1.84182884E-02;
    COFD[3708] = -2.25793511E+01;
    COFD[3709] = 5.47765924E+00;
    COFD[3710] = -4.55592027E-01;
    COFD[3711] = 1.81470399E-02;
    COFD[3712] = -2.25056806E+01;
    COFD[3713] = 5.47343682E+00;
    COFD[3714] = -4.54166709E-01;
    COFD[3715] = 1.80526774E-02;
    COFD[3716] = -2.26835080E+01;
    COFD[3717] = 5.48311530E+00;
    COFD[3718] = -4.57619182E-01;
    COFD[3719] = 1.82841374E-02;
    COFD[3720] = -2.26824612E+01;
    COFD[3721] = 5.48311530E+00;
    COFD[3722] = -4.57619182E-01;
    COFD[3723] = 1.82841374E-02;
    COFD[3724] = -2.27214167E+01;
    COFD[3725] = 5.40961679E+00;
    COFD[3726] = -4.40511249E-01;
    COFD[3727] = 1.72607396E-02;
    COFD[3728] = -2.28355542E+01;
    COFD[3729] = 5.35124813E+00;
    COFD[3730] = -4.29932425E-01;
    COFD[3731] = 1.66897210E-02;
    COFD[3732] = -2.27371462E+01;
    COFD[3733] = 5.45828076E+00;
    COFD[3734] = -4.50591566E-01;
    COFD[3735] = 1.78383531E-02;
    COFD[3736] = -2.27339269E+01;
    COFD[3737] = 5.38171422E+00;
    COFD[3738] = -4.35323468E-01;
    COFD[3739] = 1.69770246E-02;
    COFD[3740] = -1.76034182E+01;
    COFD[3741] = 4.16734978E+00;
    COFD[3742] = -3.20508152E-01;
    COFD[3743] = 1.36405589E-02;
    COFD[3744] = -1.92320503E+01;
    COFD[3745] = 5.34516253E+00;
    COFD[3746] = -4.57548788E-01;
    COFD[3747] = 1.89625374E-02;
    COFD[3748] = -1.53596434E+01;
    COFD[3749] = 3.80403213E+00;
    COFD[3750] = -2.76736813E-01;
    COFD[3751] = 1.18808568E-02;
    COFD[3752] = -1.87339672E+01;
    COFD[3753] = 4.71430273E+00;
    COFD[3754] = -3.87319126E-01;
    COFD[3755] = 1.63728857E-02;
    COFD[3756] = -2.01764857E+01;
    COFD[3757] = 5.06295779E+00;
    COFD[3758] = -4.27579258E-01;
    COFD[3759] = 1.79215546E-02;
    COFD[3760] = -1.87541329E+01;
    COFD[3761] = 4.71430273E+00;
    COFD[3762] = -3.87319126E-01;
    COFD[3763] = 1.63728857E-02;
    COFD[3764] = -1.96338867E+01;
    COFD[3765] = 4.41245751E+00;
    COFD[3766] = -2.80855500E-01;
    COFD[3767] = 9.22505777E-03;
    COFD[3768] = -1.98617110E+01;
    COFD[3769] = 4.94937905E+00;
    COFD[3770] = -4.14409534E-01;
    COFD[3771] = 1.74121859E-02;
    COFD[3772] = -2.20081706E+01;
    COFD[3773] = 5.48782481E+00;
    COFD[3774] = -4.59567900E-01;
    COFD[3775] = 1.84190115E-02;
    COFD[3776] = -2.09231900E+01;
    COFD[3777] = 5.33874924E+00;
    COFD[3778] = -4.56879272E-01;
    COFD[3779] = 1.89397603E-02;
    COFD[3780] = -2.08679661E+01;
    COFD[3781] = 5.31839685E+00;
    COFD[3782] = -4.54623948E-01;
    COFD[3783] = 1.88567295E-02;
    COFD[3784] = -2.01841845E+01;
    COFD[3785] = 5.06295779E+00;
    COFD[3786] = -4.27579258E-01;
    COFD[3787] = 1.79215546E-02;
    COFD[3788] = -2.01915368E+01;
    COFD[3789] = 5.06295779E+00;
    COFD[3790] = -4.27579258E-01;
    COFD[3791] = 1.79215546E-02;
    COFD[3792] = -2.07135807E+01;
    COFD[3793] = 4.73262810E+00;
    COFD[3794] = -3.29863658E-01;
    COFD[3795] = 1.16208974E-02;
    COFD[3796] = -2.20793238E+01;
    COFD[3797] = 5.48560730E+00;
    COFD[3798] = -4.58746417E-01;
    COFD[3799] = 1.83631968E-02;
    COFD[3800] = -2.18142322E+01;
    COFD[3801] = 5.49121786E+00;
    COFD[3802] = -4.60888505E-01;
    COFD[3803] = 1.85096783E-02;
    COFD[3804] = -2.20704568E+01;
    COFD[3805] = 5.48560730E+00;
    COFD[3806] = -4.58746417E-01;
    COFD[3807] = 1.83631968E-02;
    COFD[3808] = -2.19028963E+01;
    COFD[3809] = 5.47597765E+00;
    COFD[3810] = -4.55020526E-01;
    COFD[3811] = 1.81091205E-02;
    COFD[3812] = -2.19138830E+01;
    COFD[3813] = 5.47597765E+00;
    COFD[3814] = -4.55020526E-01;
    COFD[3815] = 1.81091205E-02;
    COFD[3816] = -2.08843468E+01;
    COFD[3817] = 4.80548676E+00;
    COFD[3818] = -3.41059239E-01;
    COFD[3819] = 1.21698161E-02;
    COFD[3820] = -2.13998650E+01;
    COFD[3821] = 4.99422177E+00;
    COFD[3822] = -3.70645569E-01;
    COFD[3823] = 1.36394075E-02;
    COFD[3824] = -2.09428491E+01;
    COFD[3825] = 5.37676571E+00;
    COFD[3826] = -4.60857509E-01;
    COFD[3827] = 1.90756062E-02;
    COFD[3828] = -2.14098610E+01;
    COFD[3829] = 4.99422177E+00;
    COFD[3830] = -3.70645569E-01;
    COFD[3831] = 1.36394075E-02;
    COFD[3832] = -2.20943989E+01;
    COFD[3833] = 5.36809191E+00;
    COFD[3834] = -4.32908721E-01;
    COFD[3835] = 1.68481857E-02;
    COFD[3836] = -2.20943989E+01;
    COFD[3837] = 5.36809191E+00;
    COFD[3838] = -4.32908721E-01;
    COFD[3839] = 1.68481857E-02;
    COFD[3840] = -2.21226089E+01;
    COFD[3841] = 5.41240338E+00;
    COFD[3842] = -4.41033923E-01;
    COFD[3843] = 1.72894624E-02;
    COFD[3844] = -2.21202988E+01;
    COFD[3845] = 5.27340167E+00;
    COFD[3846] = -4.16479301E-01;
    COFD[3847] = 1.59820820E-02;
    COFD[3848] = -2.21169561E+01;
    COFD[3849] = 5.27973457E+00;
    COFD[3850] = -4.17551922E-01;
    COFD[3851] = 1.60378887E-02;
    COFD[3852] = -2.21202699E+01;
    COFD[3853] = 5.27973457E+00;
    COFD[3854] = -4.17551922E-01;
    COFD[3855] = 1.60378887E-02;
    COFD[3856] = -2.22744971E+01;
    COFD[3857] = 5.28846317E+00;
    COFD[3858] = -4.19054070E-01;
    COFD[3859] = 1.61166909E-02;
    COFD[3860] = -2.16858752E+01;
    COFD[3861] = 4.99630533E+00;
    COFD[3862] = -3.70973757E-01;
    COFD[3863] = 1.36557558E-02;
    COFD[3864] = -2.16441409E+01;
    COFD[3865] = 4.99755422E+00;
    COFD[3866] = -3.71170522E-01;
    COFD[3867] = 1.36655589E-02;
    COFD[3868] = -2.21143598E+01;
    COFD[3869] = 5.14170132E+00;
    COFD[3870] = -3.94484367E-01;
    COFD[3871] = 1.48469922E-02;
    COFD[3872] = -2.21735855E+01;
    COFD[3873] = 5.17327578E+00;
    COFD[3874] = -3.99835717E-01;
    COFD[3875] = 1.51252623E-02;
    COFD[3876] = -2.17232952E+01;
    COFD[3877] = 4.97429231E+00;
    COFD[3878] = -3.67512538E-01;
    COFD[3879] = 1.34835182E-02;
    COFD[3880] = -2.09711478E+01;
    COFD[3881] = 4.80548676E+00;
    COFD[3882] = -3.41059239E-01;
    COFD[3883] = 1.21698161E-02;
    COFD[3884] = -2.09754168E+01;
    COFD[3885] = 4.80548676E+00;
    COFD[3886] = -3.41059239E-01;
    COFD[3887] = 1.21698161E-02;
    COFD[3888] = -2.14350985E+01;
    COFD[3889] = 4.96375446E+00;
    COFD[3890] = -3.65860553E-01;
    COFD[3891] = 1.34014562E-02;
    COFD[3892] = -2.17575745E+01;
    COFD[3893] = 5.02398766E+00;
    COFD[3894] = -3.75343697E-01;
    COFD[3895] = 1.38737246E-02;
    COFD[3896] = -2.21138917E+01;
    COFD[3897] = 5.39199985E+00;
    COFD[3898] = -4.37219136E-01;
    COFD[3899] = 1.70801890E-02;
    COFD[3900] = -2.20886975E+01;
    COFD[3901] = 5.36809191E+00;
    COFD[3902] = -4.32908721E-01;
    COFD[3903] = 1.68481857E-02;
    COFD[3904] = -2.15333456E+01;
    COFD[3905] = 4.87709635E+00;
    COFD[3906] = -3.52250416E-01;
    COFD[3907] = 1.27244832E-02;
    COFD[3908] = -2.16990357E+01;
    COFD[3909] = 4.91762397E+00;
    COFD[3910] = -3.58634830E-01;
    COFD[3911] = 1.30425594E-02;
    COFD[3912] = -2.14384166E+01;
    COFD[3913] = 4.96375446E+00;
    COFD[3914] = -3.65860553E-01;
    COFD[3915] = 1.34014562E-02;
    COFD[3916] = -2.16206093E+01;
    COFD[3917] = 4.82487053E+00;
    COFD[3918] = -3.44081206E-01;
    COFD[3919] = 1.23193175E-02;
    COFD[3920] = -2.14352422E+01;
    COFD[3921] = 4.77388530E+00;
    COFD[3922] = -3.36183750E-01;
    COFD[3923] = 1.19302239E-02;
    COFD[3924] = -2.18839486E+01;
    COFD[3925] = 4.89616559E+00;
    COFD[3926] = -3.55241555E-01;
    COFD[3927] = 1.28731550E-02;
    COFD[3928] = -2.18827191E+01;
    COFD[3929] = 4.89616559E+00;
    COFD[3930] = -3.55241555E-01;
    COFD[3931] = 1.28731550E-02;
    COFD[3932] = -2.10721493E+01;
    COFD[3933] = 4.46261834E+00;
    COFD[3934] = -2.88453410E-01;
    COFD[3935] = 9.59378935E-03;
    COFD[3936] = -2.09076816E+01;
    COFD[3937] = 4.28688746E+00;
    COFD[3938] = -2.61952393E-01;
    COFD[3939] = 8.31175423E-03;
    COFD[3940] = -2.14610939E+01;
    COFD[3941] = 4.66858934E+00;
    COFD[3942] = -3.20108476E-01;
    COFD[3943] = 1.11452199E-02;
    COFD[3944] = -2.09356172E+01;
    COFD[3945] = 4.37219452E+00;
    COFD[3946] = -2.74784562E-01;
    COFD[3947] = 8.93131387E-03;
    COFD[3948] = -1.98307830E+01;
    COFD[3949] = 4.94145734E+00;
    COFD[3950] = -4.13478429E-01;
    COFD[3951] = 1.73756310E-02;
    COFD[3952] = -1.89960613E+01;
    COFD[3953] = 5.22783910E+00;
    COFD[3954] = -4.44450099E-01;
    COFD[3955] = 1.84764493E-02;
    COFD[3956] = -1.51423689E+01;
    COFD[3957] = 3.68565858E+00;
    COFD[3958] = -2.62077329E-01;
    COFD[3959] = 1.12749663E-02;
    COFD[3960] = -1.84519843E+01;
    COFD[3961] = 4.56059591E+00;
    COFD[3962] = -3.68504274E-01;
    COFD[3963] = 1.56015944E-02;
    COFD[3964] = -1.99640254E+01;
    COFD[3965] = 4.93673368E+00;
    COFD[3966] = -4.12924135E-01;
    COFD[3967] = 1.73539126E-02;
    COFD[3968] = -1.84739227E+01;
    COFD[3969] = 4.56059591E+00;
    COFD[3970] = -3.68504274E-01;
    COFD[3971] = 1.56015944E-02;
    COFD[3972] = -2.02357367E+01;
    COFD[3973] = 4.63893949E+00;
    COFD[3974] = -3.15589295E-01;
    COFD[3975] = 1.09248423E-02;
    COFD[3976] = -1.96817244E+01;
    COFD[3977] = 4.83908504E+00;
    COFD[3978] = -4.01914941E-01;
    COFD[3979] = 1.69425288E-02;
    COFD[3980] = -2.20938193E+01;
    COFD[3981] = 5.49242695E+00;
    COFD[3982] = -4.63877105E-01;
    COFD[3983] = 1.87472502E-02;
    COFD[3984] = -2.06886834E+01;
    COFD[3985] = 5.21843630E+00;
    COFD[3986] = -4.43374339E-01;
    COFD[3987] = 1.84354402E-02;
    COFD[3988] = -2.06489560E+01;
    COFD[3989] = 5.20348257E+00;
    COFD[3990] = -4.41867814E-01;
    COFD[3991] = 1.83870405E-02;
    COFD[3992] = -1.99727712E+01;
    COFD[3993] = 4.93673368E+00;
    COFD[3994] = -4.12924135E-01;
    COFD[3995] = 1.73539126E-02;
    COFD[3996] = -1.99811407E+01;
    COFD[3997] = 4.93673368E+00;
    COFD[3998] = -4.12924135E-01;
    COFD[3999] = 1.73539126E-02;
    COFD[4000] = -2.12524252E+01;
    COFD[4001] = 4.93163464E+00;
    COFD[4002] = -3.60831131E-01;
    COFD[4003] = 1.31516943E-02;
    COFD[4004] = -2.21540958E+01;
    COFD[4005] = 5.49388984E+00;
    COFD[4006] = -4.63572253E-01;
    COFD[4007] = 1.87152221E-02;
    COFD[4008] = -2.18686233E+01;
    COFD[4009] = 5.48743241E+00;
    COFD[4010] = -4.64013540E-01;
    COFD[4011] = 1.87830225E-02;
    COFD[4012] = -2.21440868E+01;
    COFD[4013] = 5.49388984E+00;
    COFD[4014] = -4.63572253E-01;
    COFD[4015] = 1.87152221E-02;
    COFD[4016] = -2.19954815E+01;
    COFD[4017] = 5.49007189E+00;
    COFD[4018] = -4.60469654E-01;
    COFD[4019] = 1.84812405E-02;
    COFD[4020] = -2.20077630E+01;
    COFD[4021] = 5.49007189E+00;
    COFD[4022] = -4.60469654E-01;
    COFD[4023] = 1.84812405E-02;
    COFD[4024] = -2.13998650E+01;
    COFD[4025] = 4.99422177E+00;
    COFD[4026] = -3.70645569E-01;
    COFD[4027] = 1.36394075E-02;
    COFD[4028] = -2.18359714E+01;
    COFD[4029] = 5.14861103E+00;
    COFD[4030] = -3.95646185E-01;
    COFD[4031] = 1.49071576E-02;
    COFD[4032] = -2.07927921E+01;
    COFD[4033] = 5.27033141E+00;
    COFD[4034] = -4.49240926E-01;
    COFD[4035] = 1.86560940E-02;
    COFD[4036] = -2.18475446E+01;
    COFD[4037] = 5.14861103E+00;
    COFD[4038] = -3.95646185E-01;
    COFD[4039] = 1.49071576E-02;
    COFD[4040] = -2.23377280E+01;
    COFD[4041] = 5.44415012E+00;
    COFD[4042] = -4.47195848E-01;
    COFD[4043] = 1.76333776E-02;
    COFD[4044] = -2.23377280E+01;
    COFD[4045] = 5.44415012E+00;
    COFD[4046] = -4.47195848E-01;
    COFD[4047] = 1.76333776E-02;
    COFD[4048] = -2.23100112E+01;
    COFD[4049] = 5.46199540E+00;
    COFD[4050] = -4.51467141E-01;
    COFD[4051] = 1.78907506E-02;
    COFD[4052] = -2.24302648E+01;
    COFD[4053] = 5.37317627E+00;
    COFD[4054] = -4.33806826E-01;
    COFD[4055] = 1.68960241E-02;
    COFD[4056] = -2.24254137E+01;
    COFD[4057] = 5.37839633E+00;
    COFD[4058] = -4.34732804E-01;
    COFD[4059] = 1.69454459E-02;
    COFD[4060] = -2.24293161E+01;
    COFD[4061] = 5.37839633E+00;
    COFD[4062] = -4.34732804E-01;
    COFD[4063] = 1.69454459E-02;
    COFD[4064] = -2.25761948E+01;
    COFD[4065] = 5.38613682E+00;
    COFD[4066] = -4.36112477E-01;
    COFD[4067] = 1.70192630E-02;
    COFD[4068] = -2.21212336E+01;
    COFD[4069] = 5.14996158E+00;
    COFD[4070] = -3.95874432E-01;
    COFD[4071] = 1.49190064E-02;
    COFD[4072] = -2.20808599E+01;
    COFD[4073] = 5.15076979E+00;
    COFD[4074] = -3.96011114E-01;
    COFD[4075] = 1.49261042E-02;
    COFD[4076] = -2.24838323E+01;
    COFD[4077] = 5.26658797E+00;
    COFD[4078] = -4.15335910E-01;
    COFD[4079] = 1.59229288E-02;
    COFD[4080] = -2.25382696E+01;
    COFD[4081] = 5.29559685E+00;
    COFD[4082] = -4.20288187E-01;
    COFD[4083] = 1.61815942E-02;
    COFD[4084] = -2.21828607E+01;
    COFD[4085] = 5.13579028E+00;
    COFD[4086] = -3.93483830E-01;
    COFD[4087] = 1.47950287E-02;
    COFD[4088] = -2.14996220E+01;
    COFD[4089] = 4.99422177E+00;
    COFD[4090] = -3.70645569E-01;
    COFD[4091] = 1.36394075E-02;
    COFD[4092] = -2.15045977E+01;
    COFD[4093] = 4.99422177E+00;
    COFD[4094] = -3.70645569E-01;
    COFD[4095] = 1.36394075E-02;
    COFD[4096] = -2.19030787E+01;
    COFD[4097] = 5.12829880E+00;
    COFD[4098] = -3.92242339E-01;
    COFD[4099] = 1.47312337E-02;
    COFD[4100] = -2.21678146E+01;
    COFD[4101] = 5.16786857E+00;
    COFD[4102] = -3.98913959E-01;
    COFD[4103] = 1.50771754E-02;
    COFD[4104] = -2.23238237E+01;
    COFD[4105] = 5.45282659E+00;
    COFD[4106] = -4.49336536E-01;
    COFD[4107] = 1.77637787E-02;
    COFD[4108] = -2.23311644E+01;
    COFD[4109] = 5.44415012E+00;
    COFD[4110] = -4.47195848E-01;
    COFD[4111] = 1.76333776E-02;
    COFD[4112] = -2.20442423E+01;
    COFD[4113] = 5.05908558E+00;
    COFD[4114] = -3.80945929E-01;
    COFD[4115] = 1.41553816E-02;
    COFD[4116] = -2.21907779E+01;
    COFD[4117] = 5.09350157E+00;
    COFD[4118] = -3.86533585E-01;
    COFD[4119] = 1.44394093E-02;
    COFD[4120] = -2.19069859E+01;
    COFD[4121] = 5.12829880E+00;
    COFD[4122] = -3.92242339E-01;
    COFD[4123] = 1.47312337E-02;
    COFD[4124] = -2.21441058E+01;
    COFD[4125] = 5.01223834E+00;
    COFD[4126] = -3.73486752E-01;
    COFD[4127] = 1.37810367E-02;
    COFD[4128] = -2.19740179E+01;
    COFD[4129] = 4.96635252E+00;
    COFD[4130] = -3.66267779E-01;
    COFD[4131] = 1.34216831E-02;
    COFD[4132] = -2.23860697E+01;
    COFD[4133] = 5.07468436E+00;
    COFD[4134] = -3.83473296E-01;
    COFD[4135] = 1.42836887E-02;
    COFD[4136] = -2.23845704E+01;
    COFD[4137] = 5.07468436E+00;
    COFD[4138] = -3.83473296E-01;
    COFD[4139] = 1.42836887E-02;
    COFD[4140] = -2.16898608E+01;
    COFD[4141] = 4.68658569E+00;
    COFD[4142] = -3.22849054E-01;
    COFD[4143] = 1.12788192E-02;
    COFD[4144] = -2.15840392E+01;
    COFD[4145] = 4.53462922E+00;
    COFD[4146] = -2.99463362E-01;
    COFD[4147] = 1.01311942E-02;
    COFD[4148] = -2.20529119E+01;
    COFD[4149] = 4.88259016E+00;
    COFD[4150] = -3.53111373E-01;
    COFD[4151] = 1.27672516E-02;
    COFD[4152] = -2.15822613E+01;
    COFD[4153] = 4.60749458E+00;
    COFD[4154] = -3.10696534E-01;
    COFD[4155] = 1.06831054E-02;
    COFD[4156] = -1.96541192E+01;
    COFD[4157] = 4.83228605E+00;
    COFD[4158] = -4.01127235E-01;
    COFD[4159] = 1.69121265E-02;
    COFD[4160] = -1.47582497E+01;
    COFD[4161] = 3.96253155E+00;
    COFD[4162] = -2.95745106E-01;
    COFD[4163] = 1.26413618E-02;
    COFD[4164] = -1.17243965E+01;
    COFD[4165] = 2.59606889E+00;
    COFD[4166] = -1.23237035E-01;
    COFD[4167] = 5.36924411E-03;
    COFD[4168] = -1.44910354E+01;
    COFD[4169] = 3.31811351E+00;
    COFD[4170] = -2.15918664E-01;
    COFD[4171] = 9.33972748E-03;
    COFD[4172] = -1.58631082E+01;
    COFD[4173] = 3.64628858E+00;
    COFD[4174] = -2.57092933E-01;
    COFD[4175] = 1.10642427E-02;
    COFD[4176] = -1.45128230E+01;
    COFD[4177] = 3.31811351E+00;
    COFD[4178] = -2.15918664E-01;
    COFD[4179] = 9.33972748E-03;
    COFD[4180] = -2.07110846E+01;
    COFD[4181] = 5.44973198E+00;
    COFD[4182] = -4.65238176E-01;
    COFD[4183] = 1.90691386E-02;
    COFD[4184] = -1.55380429E+01;
    COFD[4185] = 3.51529689E+00;
    COFD[4186] = -2.40382392E-01;
    COFD[4187] = 1.03523419E-02;
    COFD[4188] = -1.88249737E+01;
    COFD[4189] = 4.63673045E+00;
    COFD[4190] = -3.77870004E-01;
    COFD[4191] = 1.59875055E-02;
    COFD[4192] = -1.65966953E+01;
    COFD[4193] = 3.95555732E+00;
    COFD[4194] = -2.94894134E-01;
    COFD[4195] = 1.26066350E-02;
    COFD[4196] = -1.65522787E+01;
    COFD[4197] = 3.93761009E+00;
    COFD[4198] = -2.92710565E-01;
    COFD[4199] = 1.25177995E-02;
    COFD[4200] = -1.58717613E+01;
    COFD[4201] = 3.64628858E+00;
    COFD[4202] = -2.57092933E-01;
    COFD[4203] = 1.10642427E-02;
    COFD[4204] = -1.58800406E+01;
    COFD[4205] = 3.64628858E+00;
    COFD[4206] = -2.57092933E-01;
    COFD[4207] = 1.10642427E-02;
    COFD[4208] = -2.10061234E+01;
    COFD[4209] = 5.40483692E+00;
    COFD[4210] = -4.63658873E-01;
    COFD[4211] = 1.91646114E-02;
    COFD[4212] = -1.89659293E+01;
    COFD[4213] = 4.65931164E+00;
    COFD[4214] = -3.80642811E-01;
    COFD[4215] = 1.61015362E-02;
    COFD[4216] = -1.85111018E+01;
    COFD[4217] = 4.59851812E+00;
    COFD[4218] = -3.73151692E-01;
    COFD[4219] = 1.57923348E-02;
    COFD[4220] = -1.89560208E+01;
    COFD[4221] = 4.65931164E+00;
    COFD[4222] = -3.80642811E-01;
    COFD[4223] = 1.61015362E-02;
    COFD[4224] = -1.89929565E+01;
    COFD[4225] = 4.75627153E+00;
    COFD[4226] = -3.92329566E-01;
    COFD[4227] = 1.65729779E-02;
    COFD[4228] = -1.90051251E+01;
    COFD[4229] = 4.75627153E+00;
    COFD[4230] = -3.92329566E-01;
    COFD[4231] = 1.65729779E-02;
    COFD[4232] = -2.09428491E+01;
    COFD[4233] = 5.37676571E+00;
    COFD[4234] = -4.60857509E-01;
    COFD[4235] = 1.90756062E-02;
    COFD[4236] = -2.07927921E+01;
    COFD[4237] = 5.27033141E+00;
    COFD[4238] = -4.49240926E-01;
    COFD[4239] = 1.86560940E-02;
    COFD[4240] = -1.65593994E+01;
    COFD[4241] = 3.99552704E+00;
    COFD[4242] = -2.99734066E-01;
    COFD[4243] = 1.28024419E-02;
    COFD[4244] = -2.08042232E+01;
    COFD[4245] = 5.27033141E+00;
    COFD[4246] = -4.49240926E-01;
    COFD[4247] = 1.86560940E-02;
    COFD[4248] = -1.99946981E+01;
    COFD[4249] = 4.99038775E+00;
    COFD[4250] = -4.19232899E-01;
    COFD[4251] = 1.76017333E-02;
    COFD[4252] = -1.99946981E+01;
    COFD[4253] = 4.99038775E+00;
    COFD[4254] = -4.19232899E-01;
    COFD[4255] = 1.76017333E-02;
    COFD[4256] = -1.97503281E+01;
    COFD[4257] = 4.91724977E+00;
    COFD[4258] = -4.10632023E-01;
    COFD[4259] = 1.72638615E-02;
    COFD[4260] = -2.04982298E+01;
    COFD[4261] = 5.09287074E+00;
    COFD[4262] = -4.30747337E-01;
    COFD[4263] = 1.80311502E-02;
    COFD[4264] = -2.04673093E+01;
    COFD[4265] = 5.08807429E+00;
    COFD[4266] = -4.30243518E-01;
    COFD[4267] = 1.80139269E-02;
    COFD[4268] = -2.04711578E+01;
    COFD[4269] = 5.08807429E+00;
    COFD[4270] = -4.30243518E-01;
    COFD[4271] = 1.80139269E-02;
    COFD[4272] = -2.06099144E+01;
    COFD[4273] = 5.08069320E+00;
    COFD[4274] = -4.29462569E-01;
    COFD[4275] = 1.79869461E-02;
    COFD[4276] = -2.11159694E+01;
    COFD[4277] = 5.26890914E+00;
    COFD[4278] = -4.49080646E-01;
    COFD[4279] = 1.86500850E-02;
    COFD[4280] = -2.10624904E+01;
    COFD[4281] = 5.26805363E+00;
    COFD[4282] = -4.48984207E-01;
    COFD[4283] = 1.86464682E-02;
    COFD[4284] = -2.10353566E+01;
    COFD[4285] = 5.18365882E+00;
    COFD[4286] = -4.40195806E-01;
    COFD[4287] = 1.83497181E-02;
    COFD[4288] = -2.09861375E+01;
    COFD[4289] = 5.16708806E+00;
    COFD[4290] = -4.38725884E-01;
    COFD[4291] = 1.83125969E-02;
    COFD[4292] = -2.12518544E+01;
    COFD[4293] = 5.28351879E+00;
    COFD[4294] = -4.50726414E-01;
    COFD[4295] = 1.87118076E-02;
    COFD[4296] = -2.10414466E+01;
    COFD[4297] = 5.37676571E+00;
    COFD[4298] = -4.60857509E-01;
    COFD[4299] = 1.90756062E-02;
    COFD[4300] = -2.10463583E+01;
    COFD[4301] = 5.37676571E+00;
    COFD[4302] = -4.60857509E-01;
    COFD[4303] = 1.90756062E-02;
    COFD[4304] = -2.09608951E+01;
    COFD[4305] = 5.29024160E+00;
    COFD[4306] = -4.51483508E-01;
    COFD[4307] = 1.87402157E-02;
    COFD[4308] = -2.10807241E+01;
    COFD[4309] = 5.24972538E+00;
    COFD[4310] = -4.46923248E-01;
    COFD[4311] = 1.85694358E-02;
    COFD[4312] = -1.98725502E+01;
    COFD[4313] = 4.95323342E+00;
    COFD[4314] = -4.14861803E-01;
    COFD[4315] = 1.74299070E-02;
    COFD[4316] = -1.99882118E+01;
    COFD[4317] = 4.99038775E+00;
    COFD[4318] = -4.19232899E-01;
    COFD[4319] = 1.76017333E-02;
    COFD[4320] = -2.14253720E+01;
    COFD[4321] = 5.34211391E+00;
    COFD[4322] = -4.57230216E-01;
    COFD[4323] = 1.89516739E-02;
    COFD[4324] = -2.14529319E+01;
    COFD[4325] = 5.31795346E+00;
    COFD[4326] = -4.54574565E-01;
    COFD[4327] = 1.88549001E-02;
    COFD[4328] = -2.09647485E+01;
    COFD[4329] = 5.29024160E+00;
    COFD[4330] = -4.51483508E-01;
    COFD[4331] = 1.87402157E-02;
    COFD[4332] = -2.17208609E+01;
    COFD[4333] = 5.36702791E+00;
    COFD[4334] = -4.59828733E-01;
    COFD[4335] = 1.90400208E-02;
    COFD[4336] = -2.17030104E+01;
    COFD[4337] = 5.39237617E+00;
    COFD[4338] = -4.62530602E-01;
    COFD[4339] = 1.91345681E-02;
    COFD[4340] = -2.17533197E+01;
    COFD[4341] = 5.33181444E+00;
    COFD[4342] = -4.56115824E-01;
    COFD[4343] = 1.89118833E-02;
    COFD[4344] = -2.17518459E+01;
    COFD[4345] = 5.33181444E+00;
    COFD[4346] = -4.56115824E-01;
    COFD[4347] = 1.89118833E-02;
    COFD[4348] = -2.22344688E+01;
    COFD[4349] = 5.44527718E+00;
    COFD[4350] = -4.65287801E-01;
    COFD[4351] = 1.90954247E-02;
    COFD[4352] = -2.25529616E+01;
    COFD[4353] = 5.46297654E+00;
    COFD[4354] = -4.65485734E-01;
    COFD[4355] = 1.90234280E-02;
    COFD[4356] = -2.20783474E+01;
    COFD[4357] = 5.42074110E+00;
    COFD[4358] = -4.64982008E-01;
    COFD[4359] = 1.91934150E-02;
    COFD[4360] = -2.23385572E+01;
    COFD[4361] = 5.45382388E+00;
    COFD[4362] = -4.65297200E-01;
    COFD[4363] = 1.90536055E-02;
    COFD[4364] = -1.55055897E+01;
    COFD[4365] = 3.50663983E+00;
    COFD[4366] = -2.39273841E-01;
    COFD[4367] = 1.03049488E-02;
    COFD[4368] = -1.89965974E+01;
    COFD[4369] = 5.22783910E+00;
    COFD[4370] = -4.44450099E-01;
    COFD[4371] = 1.84764493E-02;
    COFD[4372] = -1.51434170E+01;
    COFD[4373] = 3.68565858E+00;
    COFD[4374] = -2.62077329E-01;
    COFD[4375] = 1.12749663E-02;
    COFD[4376] = -1.84583320E+01;
    COFD[4377] = 4.56059591E+00;
    COFD[4378] = -3.68504274E-01;
    COFD[4379] = 1.56015944E-02;
    COFD[4380] = -1.99740136E+01;
    COFD[4381] = 4.93673368E+00;
    COFD[4382] = -4.12924135E-01;
    COFD[4383] = 1.73539126E-02;
    COFD[4384] = -1.84805571E+01;
    COFD[4385] = 4.56059591E+00;
    COFD[4386] = -3.68504274E-01;
    COFD[4387] = 1.56015944E-02;
    COFD[4388] = -2.02426482E+01;
    COFD[4389] = 4.63893949E+00;
    COFD[4390] = -3.15589295E-01;
    COFD[4391] = 1.09248423E-02;
    COFD[4392] = -1.96909586E+01;
    COFD[4393] = 4.83908504E+00;
    COFD[4394] = -4.01914941E-01;
    COFD[4395] = 1.69425288E-02;
    COFD[4396] = -2.21056610E+01;
    COFD[4397] = 5.49242695E+00;
    COFD[4398] = -4.63877105E-01;
    COFD[4399] = 1.87472502E-02;
    COFD[4400] = -2.06947476E+01;
    COFD[4401] = 5.21843630E+00;
    COFD[4402] = -4.43374339E-01;
    COFD[4403] = 1.84354402E-02;
    COFD[4404] = -2.06553163E+01;
    COFD[4405] = 5.20348257E+00;
    COFD[4406] = -4.41867814E-01;
    COFD[4407] = 1.83870405E-02;
    COFD[4408] = -1.99829374E+01;
    COFD[4409] = 4.93673368E+00;
    COFD[4410] = -4.12924135E-01;
    COFD[4411] = 1.73539126E-02;
    COFD[4412] = -1.99914803E+01;
    COFD[4413] = 4.93673368E+00;
    COFD[4414] = -4.12924135E-01;
    COFD[4415] = 1.73539126E-02;
    COFD[4416] = -2.12620508E+01;
    COFD[4417] = 4.93163464E+00;
    COFD[4418] = -3.60831131E-01;
    COFD[4419] = 1.31516943E-02;
    COFD[4420] = -2.21637296E+01;
    COFD[4421] = 5.49388984E+00;
    COFD[4422] = -4.63572253E-01;
    COFD[4423] = 1.87152221E-02;
    COFD[4424] = -2.18778661E+01;
    COFD[4425] = 5.48743241E+00;
    COFD[4426] = -4.64013540E-01;
    COFD[4427] = 1.87830225E-02;
    COFD[4428] = -2.21535278E+01;
    COFD[4429] = 5.49388984E+00;
    COFD[4430] = -4.63572253E-01;
    COFD[4431] = 1.87152221E-02;
    COFD[4432] = -2.20040945E+01;
    COFD[4433] = 5.49007189E+00;
    COFD[4434] = -4.60469654E-01;
    COFD[4435] = 1.84812405E-02;
    COFD[4436] = -2.20165921E+01;
    COFD[4437] = 5.49007189E+00;
    COFD[4438] = -4.60469654E-01;
    COFD[4439] = 1.84812405E-02;
    COFD[4440] = -2.14098610E+01;
    COFD[4441] = 4.99422177E+00;
    COFD[4442] = -3.70645569E-01;
    COFD[4443] = 1.36394075E-02;
    COFD[4444] = -2.18475446E+01;
    COFD[4445] = 5.14861103E+00;
    COFD[4446] = -3.95646185E-01;
    COFD[4447] = 1.49071576E-02;
    COFD[4448] = -2.08042232E+01;
    COFD[4449] = 5.27033141E+00;
    COFD[4450] = -4.49240926E-01;
    COFD[4451] = 1.86560940E-02;
    COFD[4452] = -2.18593920E+01;
    COFD[4453] = 5.14861103E+00;
    COFD[4454] = -3.95646185E-01;
    COFD[4455] = 1.49071576E-02;
    COFD[4456] = -2.23490200E+01;
    COFD[4457] = 5.44415012E+00;
    COFD[4458] = -4.47195848E-01;
    COFD[4459] = 1.76333776E-02;
    COFD[4460] = -2.23490200E+01;
    COFD[4461] = 5.44415012E+00;
    COFD[4462] = -4.47195848E-01;
    COFD[4463] = 1.76333776E-02;
    COFD[4464] = -2.23215905E+01;
    COFD[4465] = 5.46199540E+00;
    COFD[4466] = -4.51467141E-01;
    COFD[4467] = 1.78907506E-02;
    COFD[4468] = -2.24433083E+01;
    COFD[4469] = 5.37317627E+00;
    COFD[4470] = -4.33806826E-01;
    COFD[4471] = 1.68960241E-02;
    COFD[4472] = -2.24385639E+01;
    COFD[4473] = 5.37839633E+00;
    COFD[4474] = -4.34732804E-01;
    COFD[4475] = 1.69454459E-02;
    COFD[4476] = -2.24425707E+01;
    COFD[4477] = 5.37839633E+00;
    COFD[4478] = -4.34732804E-01;
    COFD[4479] = 1.69454459E-02;
    COFD[4480] = -2.25895516E+01;
    COFD[4481] = 5.38613682E+00;
    COFD[4482] = -4.36112477E-01;
    COFD[4483] = 1.70192630E-02;
    COFD[4484] = -2.21345861E+01;
    COFD[4485] = 5.14996158E+00;
    COFD[4486] = -3.95874432E-01;
    COFD[4487] = 1.49190064E-02;
    COFD[4488] = -2.20943127E+01;
    COFD[4489] = 5.15076979E+00;
    COFD[4490] = -3.96011114E-01;
    COFD[4491] = 1.49261042E-02;
    COFD[4492] = -2.24982668E+01;
    COFD[4493] = 5.26658797E+00;
    COFD[4494] = -4.15335910E-01;
    COFD[4495] = 1.59229288E-02;
    COFD[4496] = -2.25527842E+01;
    COFD[4497] = 5.29559685E+00;
    COFD[4498] = -4.20288187E-01;
    COFD[4499] = 1.61815942E-02;
    COFD[4500] = -2.21974539E+01;
    COFD[4501] = 5.13579028E+00;
    COFD[4502] = -3.93483830E-01;
    COFD[4503] = 1.47950287E-02;
    COFD[4504] = -2.15118523E+01;
    COFD[4505] = 4.99422177E+00;
    COFD[4506] = -3.70645569E-01;
    COFD[4507] = 1.36394075E-02;
    COFD[4508] = -2.15169519E+01;
    COFD[4509] = 4.99422177E+00;
    COFD[4510] = -3.70645569E-01;
    COFD[4511] = 1.36394075E-02;
    COFD[4512] = -2.19162243E+01;
    COFD[4513] = 5.12829880E+00;
    COFD[4514] = -3.92242339E-01;
    COFD[4515] = 1.47312337E-02;
    COFD[4516] = -2.21810647E+01;
    COFD[4517] = 5.16786857E+00;
    COFD[4518] = -3.98913959E-01;
    COFD[4519] = 1.50771754E-02;
    COFD[4520] = -2.23352610E+01;
    COFD[4521] = 5.45282659E+00;
    COFD[4522] = -4.49336536E-01;
    COFD[4523] = 1.77637787E-02;
    COFD[4524] = -2.23423075E+01;
    COFD[4525] = 5.44415012E+00;
    COFD[4526] = -4.47195848E-01;
    COFD[4527] = 1.76333776E-02;
    COFD[4528] = -2.20589094E+01;
    COFD[4529] = 5.05908558E+00;
    COFD[4530] = -3.80945929E-01;
    COFD[4531] = 1.41553816E-02;
    COFD[4532] = -2.22053677E+01;
    COFD[4533] = 5.09350157E+00;
    COFD[4534] = -3.86533585E-01;
    COFD[4535] = 1.44394093E-02;
    COFD[4536] = -2.19202360E+01;
    COFD[4537] = 5.12829880E+00;
    COFD[4538] = -3.92242339E-01;
    COFD[4539] = 1.47312337E-02;
    COFD[4540] = -2.21597235E+01;
    COFD[4541] = 5.01223834E+00;
    COFD[4542] = -3.73486752E-01;
    COFD[4543] = 1.37810367E-02;
    COFD[4544] = -2.19895742E+01;
    COFD[4545] = 4.96635252E+00;
    COFD[4546] = -3.66267779E-01;
    COFD[4547] = 1.34216831E-02;
    COFD[4548] = -2.24024538E+01;
    COFD[4549] = 5.07468436E+00;
    COFD[4550] = -3.83473296E-01;
    COFD[4551] = 1.42836887E-02;
    COFD[4552] = -2.24009047E+01;
    COFD[4553] = 5.07468436E+00;
    COFD[4554] = -3.83473296E-01;
    COFD[4555] = 1.42836887E-02;
    COFD[4556] = -2.17074965E+01;
    COFD[4557] = 4.68658569E+00;
    COFD[4558] = -3.22849054E-01;
    COFD[4559] = 1.12788192E-02;
    COFD[4560] = -2.16025724E+01;
    COFD[4561] = 4.53462922E+00;
    COFD[4562] = -2.99463362E-01;
    COFD[4563] = 1.01311942E-02;
    COFD[4564] = -2.20699218E+01;
    COFD[4565] = 4.88259016E+00;
    COFD[4566] = -3.53111373E-01;
    COFD[4567] = 1.27672516E-02;
    COFD[4568] = -2.16003553E+01;
    COFD[4569] = 4.60749458E+00;
    COFD[4570] = -3.10696534E-01;
    COFD[4571] = 1.06831054E-02;
    COFD[4572] = -1.96633539E+01;
    COFD[4573] = 4.83228605E+00;
    COFD[4574] = -4.01127235E-01;
    COFD[4575] = 1.69121265E-02;
    COFD[4576] = -1.82115937E+01;
    COFD[4577] = 4.94453558E+00;
    COFD[4578] = -4.13840584E-01;
    COFD[4579] = 1.73898654E-02;
    COFD[4580] = -1.43067326E+01;
    COFD[4581] = 3.34609996E+00;
    COFD[4582] = -2.19417637E-01;
    COFD[4583] = 9.48565577E-03;
    COFD[4584] = -1.74524239E+01;
    COFD[4585] = 4.17306000E+00;
    COFD[4586] = -3.21208070E-01;
    COFD[4587] = 1.36692112E-02;
    COFD[4588] = -1.89513110E+01;
    COFD[4589] = 4.56070745E+00;
    COFD[4590] = -3.68517396E-01;
    COFD[4591] = 1.56021102E-02;
    COFD[4592] = -1.74740622E+01;
    COFD[4593] = 4.17306000E+00;
    COFD[4594] = -3.21208070E-01;
    COFD[4595] = 1.36692112E-02;
    COFD[4596] = -2.13313704E+01;
    COFD[4597] = 5.17422070E+00;
    COFD[4598] = -3.99997068E-01;
    COFD[4599] = 1.51336871E-02;
    COFD[4600] = -1.86750745E+01;
    COFD[4601] = 4.46203550E+00;
    COFD[4602] = -3.56905086E-01;
    COFD[4603] = 1.51454720E-02;
    COFD[4604] = -2.17754519E+01;
    COFD[4605] = 5.43321422E+00;
    COFD[4606] = -4.65643875E-01;
    COFD[4607] = 1.91852144E-02;
    COFD[4608] = -1.98761303E+01;
    COFD[4609] = 4.93512898E+00;
    COFD[4610] = -4.12735698E-01;
    COFD[4611] = 1.73465233E-02;
    COFD[4612] = -1.98121403E+01;
    COFD[4613] = 4.90949838E+00;
    COFD[4614] = -4.09710987E-01;
    COFD[4615] = 1.72272850E-02;
    COFD[4616] = -1.89598731E+01;
    COFD[4617] = 4.56070745E+00;
    COFD[4618] = -3.68517396E-01;
    COFD[4619] = 1.56021102E-02;
    COFD[4620] = -1.89680638E+01;
    COFD[4621] = 4.56070745E+00;
    COFD[4622] = -3.68517396E-01;
    COFD[4623] = 1.56021102E-02;
    COFD[4624] = -2.20218901E+01;
    COFD[4625] = 5.33662478E+00;
    COFD[4626] = -4.27371229E-01;
    COFD[4627] = 1.65539339E-02;
    COFD[4628] = -2.18363989E+01;
    COFD[4629] = 5.43667692E+00;
    COFD[4630] = -4.65640530E-01;
    COFD[4631] = 1.91675299E-02;
    COFD[4632] = -2.15488001E+01;
    COFD[4633] = 5.42458950E+00;
    COFD[4634] = -4.65212987E-01;
    COFD[4635] = 1.91930766E-02;
    COFD[4636] = -2.18265893E+01;
    COFD[4637] = 5.43667692E+00;
    COFD[4638] = -4.65640530E-01;
    COFD[4639] = 1.91675299E-02;
    COFD[4640] = -2.17210686E+01;
    COFD[4641] = 5.45061574E+00;
    COFD[4642] = -4.65241893E-01;
    COFD[4643] = 1.90650429E-02;
    COFD[4644] = -2.17331258E+01;
    COFD[4645] = 5.45061574E+00;
    COFD[4646] = -4.65241893E-01;
    COFD[4647] = 1.90650429E-02;
    COFD[4648] = -2.20943989E+01;
    COFD[4649] = 5.36809191E+00;
    COFD[4650] = -4.32908721E-01;
    COFD[4651] = 1.68481857E-02;
    COFD[4652] = -2.23377280E+01;
    COFD[4653] = 5.44415012E+00;
    COFD[4654] = -4.47195848E-01;
    COFD[4655] = 1.76333776E-02;
    COFD[4656] = -1.99946981E+01;
    COFD[4657] = 4.99038775E+00;
    COFD[4658] = -4.19232899E-01;
    COFD[4659] = 1.76017333E-02;
    COFD[4660] = -2.23490200E+01;
    COFD[4661] = 5.44415012E+00;
    COFD[4662] = -4.47195848E-01;
    COFD[4663] = 1.76333776E-02;
    COFD[4664] = -2.22545120E+01;
    COFD[4665] = 5.49056182E+00;
    COFD[4666] = -4.63917767E-01;
    COFD[4667] = 1.87597656E-02;
    COFD[4668] = -2.22545120E+01;
    COFD[4669] = 5.49056182E+00;
    COFD[4670] = -4.63917767E-01;
    COFD[4671] = 1.87597656E-02;
    COFD[4672] = -2.21467458E+01;
    COFD[4673] = 5.47453801E+00;
    COFD[4674] = -4.63573634E-01;
    COFD[4675] = 1.88106950E-02;
    COFD[4676] = -2.25051541E+01;
    COFD[4677] = 5.48975794E+00;
    COFD[4678] = -4.60336840E-01;
    COFD[4679] = 1.84719899E-02;
    COFD[4680] = -2.24900659E+01;
    COFD[4681] = 5.49046566E+00;
    COFD[4682] = -4.60640519E-01;
    COFD[4683] = 1.84931917E-02;
    COFD[4684] = -2.24938619E+01;
    COFD[4685] = 5.49046566E+00;
    COFD[4686] = -4.60640519E-01;
    COFD[4687] = 1.84931917E-02;
    COFD[4688] = -2.26196129E+01;
    COFD[4689] = 5.49130890E+00;
    COFD[4690] = -4.61073466E-01;
    COFD[4691] = 1.85242949E-02;
    COFD[4692] = -2.26109161E+01;
    COFD[4693] = 5.44450212E+00;
    COFD[4694] = -4.47289567E-01;
    COFD[4695] = 1.76392556E-02;
    COFD[4696] = -2.25702108E+01;
    COFD[4697] = 5.44454618E+00;
    COFD[4698] = -4.47321672E-01;
    COFD[4699] = 1.76416355E-02;
    COFD[4700] = -2.27673763E+01;
    COFD[4701] = 5.47456971E+00;
    COFD[4702] = -4.54548739E-01;
    COFD[4703] = 1.80779412E-02;
    COFD[4704] = -2.27646473E+01;
    COFD[4705] = 5.47937827E+00;
    COFD[4706] = -4.56129749E-01;
    COFD[4707] = 1.81820758E-02;
    COFD[4708] = -2.26859156E+01;
    COFD[4709] = 5.43682302E+00;
    COFD[4710] = -4.45758386E-01;
    COFD[4711] = 1.75526921E-02;
    COFD[4712] = -2.21918593E+01;
    COFD[4713] = 5.36809191E+00;
    COFD[4714] = -4.32908721E-01;
    COFD[4715] = 1.68481857E-02;
    COFD[4716] = -2.21967084E+01;
    COFD[4717] = 5.36809191E+00;
    COFD[4718] = -4.32908721E-01;
    COFD[4719] = 1.68481857E-02;
    COFD[4720] = -2.24246556E+01;
    COFD[4721] = 5.43317403E+00;
    COFD[4722] = -4.45039640E-01;
    COFD[4723] = 1.75123239E-02;
    COFD[4724] = -2.26244858E+01;
    COFD[4725] = 5.44849005E+00;
    COFD[4726] = -4.48374875E-01;
    COFD[4727] = 1.77073338E-02;
    COFD[4728] = -2.22041980E+01;
    COFD[4729] = 5.48376307E+00;
    COFD[4730] = -4.63942837E-01;
    COFD[4731] = 1.87953234E-02;
    COFD[4732] = -2.22481014E+01;
    COFD[4733] = 5.49056182E+00;
    COFD[4734] = -4.63917767E-01;
    COFD[4735] = 1.87597656E-02;
    COFD[4736] = -2.26414192E+01;
    COFD[4737] = 5.39941661E+00;
    COFD[4738] = -4.38608992E-01;
    COFD[4739] = 1.71564833E-02;
    COFD[4740] = -2.27438518E+01;
    COFD[4741] = 5.41656161E+00;
    COFD[4742] = -4.41824209E-01;
    COFD[4743] = 1.73331395E-02;
    COFD[4744] = -2.24284563E+01;
    COFD[4745] = 5.43317403E+00;
    COFD[4746] = -4.45039640E-01;
    COFD[4747] = 1.75123239E-02;
    COFD[4748] = -2.27946566E+01;
    COFD[4749] = 5.37750827E+00;
    COFD[4750] = -4.34574933E-01;
    COFD[4751] = 1.69370116E-02;
    COFD[4752] = -2.26809417E+01;
    COFD[4753] = 5.35408130E+00;
    COFD[4754] = -4.30433039E-01;
    COFD[4755] = 1.67163724E-02;
    COFD[4756] = -2.29546228E+01;
    COFD[4757] = 5.40726595E+00;
    COFD[4758] = -4.40072331E-01;
    COFD[4759] = 1.72366684E-02;
    COFD[4760] = -2.29531738E+01;
    COFD[4761] = 5.40726595E+00;
    COFD[4762] = -4.40072331E-01;
    COFD[4763] = 1.72366684E-02;
    COFD[4764] = -2.26847799E+01;
    COFD[4765] = 5.19870800E+00;
    COFD[4766] = -4.04106072E-01;
    COFD[4767] = 1.53462074E-02;
    COFD[4768] = -2.27363668E+01;
    COFD[4769] = 5.11476665E+00;
    COFD[4770] = -3.90017950E-01;
    COFD[4771] = 1.46174143E-02;
    COFD[4772] = -2.28434554E+01;
    COFD[4773] = 5.30874092E+00;
    COFD[4774] = -4.22570517E-01;
    COFD[4775] = 1.63018518E-02;
    COFD[4776] = -2.26645375E+01;
    COFD[4777] = 5.15653775E+00;
    COFD[4778] = -3.96988537E-01;
    COFD[4779] = 1.49769158E-02;
    COFD[4780] = -1.86480131E+01;
    COFD[4781] = 4.45508831E+00;
    COFD[4782] = -3.56068927E-01;
    COFD[4783] = 1.51117769E-02;
    COFD[4784] = -1.82115937E+01;
    COFD[4785] = 4.94453558E+00;
    COFD[4786] = -4.13840584E-01;
    COFD[4787] = 1.73898654E-02;
    COFD[4788] = -1.43067326E+01;
    COFD[4789] = 3.34609996E+00;
    COFD[4790] = -2.19417637E-01;
    COFD[4791] = 9.48565577E-03;
    COFD[4792] = -1.74524239E+01;
    COFD[4793] = 4.17306000E+00;
    COFD[4794] = -3.21208070E-01;
    COFD[4795] = 1.36692112E-02;
    COFD[4796] = -1.89513110E+01;
    COFD[4797] = 4.56070745E+00;
    COFD[4798] = -3.68517396E-01;
    COFD[4799] = 1.56021102E-02;
    COFD[4800] = -1.74740622E+01;
    COFD[4801] = 4.17306000E+00;
    COFD[4802] = -3.21208070E-01;
    COFD[4803] = 1.36692112E-02;
    COFD[4804] = -2.13313704E+01;
    COFD[4805] = 5.17422070E+00;
    COFD[4806] = -3.99997068E-01;
    COFD[4807] = 1.51336871E-02;
    COFD[4808] = -1.86750745E+01;
    COFD[4809] = 4.46203550E+00;
    COFD[4810] = -3.56905086E-01;
    COFD[4811] = 1.51454720E-02;
    COFD[4812] = -2.17754519E+01;
    COFD[4813] = 5.43321422E+00;
    COFD[4814] = -4.65643875E-01;
    COFD[4815] = 1.91852144E-02;
    COFD[4816] = -1.98761303E+01;
    COFD[4817] = 4.93512898E+00;
    COFD[4818] = -4.12735698E-01;
    COFD[4819] = 1.73465233E-02;
    COFD[4820] = -1.98121403E+01;
    COFD[4821] = 4.90949838E+00;
    COFD[4822] = -4.09710987E-01;
    COFD[4823] = 1.72272850E-02;
    COFD[4824] = -1.89598731E+01;
    COFD[4825] = 4.56070745E+00;
    COFD[4826] = -3.68517396E-01;
    COFD[4827] = 1.56021102E-02;
    COFD[4828] = -1.89680638E+01;
    COFD[4829] = 4.56070745E+00;
    COFD[4830] = -3.68517396E-01;
    COFD[4831] = 1.56021102E-02;
    COFD[4832] = -2.20218901E+01;
    COFD[4833] = 5.33662478E+00;
    COFD[4834] = -4.27371229E-01;
    COFD[4835] = 1.65539339E-02;
    COFD[4836] = -2.18363989E+01;
    COFD[4837] = 5.43667692E+00;
    COFD[4838] = -4.65640530E-01;
    COFD[4839] = 1.91675299E-02;
    COFD[4840] = -2.15488001E+01;
    COFD[4841] = 5.42458950E+00;
    COFD[4842] = -4.65212987E-01;
    COFD[4843] = 1.91930766E-02;
    COFD[4844] = -2.18265893E+01;
    COFD[4845] = 5.43667692E+00;
    COFD[4846] = -4.65640530E-01;
    COFD[4847] = 1.91675299E-02;
    COFD[4848] = -2.17210686E+01;
    COFD[4849] = 5.45061574E+00;
    COFD[4850] = -4.65241893E-01;
    COFD[4851] = 1.90650429E-02;
    COFD[4852] = -2.17331258E+01;
    COFD[4853] = 5.45061574E+00;
    COFD[4854] = -4.65241893E-01;
    COFD[4855] = 1.90650429E-02;
    COFD[4856] = -2.20943989E+01;
    COFD[4857] = 5.36809191E+00;
    COFD[4858] = -4.32908721E-01;
    COFD[4859] = 1.68481857E-02;
    COFD[4860] = -2.23377280E+01;
    COFD[4861] = 5.44415012E+00;
    COFD[4862] = -4.47195848E-01;
    COFD[4863] = 1.76333776E-02;
    COFD[4864] = -1.99946981E+01;
    COFD[4865] = 4.99038775E+00;
    COFD[4866] = -4.19232899E-01;
    COFD[4867] = 1.76017333E-02;
    COFD[4868] = -2.23490200E+01;
    COFD[4869] = 5.44415012E+00;
    COFD[4870] = -4.47195848E-01;
    COFD[4871] = 1.76333776E-02;
    COFD[4872] = -2.22545120E+01;
    COFD[4873] = 5.49056182E+00;
    COFD[4874] = -4.63917767E-01;
    COFD[4875] = 1.87597656E-02;
    COFD[4876] = -2.22545120E+01;
    COFD[4877] = 5.49056182E+00;
    COFD[4878] = -4.63917767E-01;
    COFD[4879] = 1.87597656E-02;
    COFD[4880] = -2.21467458E+01;
    COFD[4881] = 5.47453801E+00;
    COFD[4882] = -4.63573634E-01;
    COFD[4883] = 1.88106950E-02;
    COFD[4884] = -2.25051541E+01;
    COFD[4885] = 5.48975794E+00;
    COFD[4886] = -4.60336840E-01;
    COFD[4887] = 1.84719899E-02;
    COFD[4888] = -2.24900659E+01;
    COFD[4889] = 5.49046566E+00;
    COFD[4890] = -4.60640519E-01;
    COFD[4891] = 1.84931917E-02;
    COFD[4892] = -2.24938619E+01;
    COFD[4893] = 5.49046566E+00;
    COFD[4894] = -4.60640519E-01;
    COFD[4895] = 1.84931917E-02;
    COFD[4896] = -2.26196129E+01;
    COFD[4897] = 5.49130890E+00;
    COFD[4898] = -4.61073466E-01;
    COFD[4899] = 1.85242949E-02;
    COFD[4900] = -2.26109161E+01;
    COFD[4901] = 5.44450212E+00;
    COFD[4902] = -4.47289567E-01;
    COFD[4903] = 1.76392556E-02;
    COFD[4904] = -2.25702108E+01;
    COFD[4905] = 5.44454618E+00;
    COFD[4906] = -4.47321672E-01;
    COFD[4907] = 1.76416355E-02;
    COFD[4908] = -2.27673763E+01;
    COFD[4909] = 5.47456971E+00;
    COFD[4910] = -4.54548739E-01;
    COFD[4911] = 1.80779412E-02;
    COFD[4912] = -2.27646473E+01;
    COFD[4913] = 5.47937827E+00;
    COFD[4914] = -4.56129749E-01;
    COFD[4915] = 1.81820758E-02;
    COFD[4916] = -2.26859156E+01;
    COFD[4917] = 5.43682302E+00;
    COFD[4918] = -4.45758386E-01;
    COFD[4919] = 1.75526921E-02;
    COFD[4920] = -2.21918593E+01;
    COFD[4921] = 5.36809191E+00;
    COFD[4922] = -4.32908721E-01;
    COFD[4923] = 1.68481857E-02;
    COFD[4924] = -2.21967084E+01;
    COFD[4925] = 5.36809191E+00;
    COFD[4926] = -4.32908721E-01;
    COFD[4927] = 1.68481857E-02;
    COFD[4928] = -2.24246556E+01;
    COFD[4929] = 5.43317403E+00;
    COFD[4930] = -4.45039640E-01;
    COFD[4931] = 1.75123239E-02;
    COFD[4932] = -2.26244858E+01;
    COFD[4933] = 5.44849005E+00;
    COFD[4934] = -4.48374875E-01;
    COFD[4935] = 1.77073338E-02;
    COFD[4936] = -2.22041980E+01;
    COFD[4937] = 5.48376307E+00;
    COFD[4938] = -4.63942837E-01;
    COFD[4939] = 1.87953234E-02;
    COFD[4940] = -2.22481014E+01;
    COFD[4941] = 5.49056182E+00;
    COFD[4942] = -4.63917767E-01;
    COFD[4943] = 1.87597656E-02;
    COFD[4944] = -2.26414192E+01;
    COFD[4945] = 5.39941661E+00;
    COFD[4946] = -4.38608992E-01;
    COFD[4947] = 1.71564833E-02;
    COFD[4948] = -2.27438518E+01;
    COFD[4949] = 5.41656161E+00;
    COFD[4950] = -4.41824209E-01;
    COFD[4951] = 1.73331395E-02;
    COFD[4952] = -2.24284563E+01;
    COFD[4953] = 5.43317403E+00;
    COFD[4954] = -4.45039640E-01;
    COFD[4955] = 1.75123239E-02;
    COFD[4956] = -2.27946566E+01;
    COFD[4957] = 5.37750827E+00;
    COFD[4958] = -4.34574933E-01;
    COFD[4959] = 1.69370116E-02;
    COFD[4960] = -2.26809417E+01;
    COFD[4961] = 5.35408130E+00;
    COFD[4962] = -4.30433039E-01;
    COFD[4963] = 1.67163724E-02;
    COFD[4964] = -2.29546228E+01;
    COFD[4965] = 5.40726595E+00;
    COFD[4966] = -4.40072331E-01;
    COFD[4967] = 1.72366684E-02;
    COFD[4968] = -2.29531738E+01;
    COFD[4969] = 5.40726595E+00;
    COFD[4970] = -4.40072331E-01;
    COFD[4971] = 1.72366684E-02;
    COFD[4972] = -2.26847799E+01;
    COFD[4973] = 5.19870800E+00;
    COFD[4974] = -4.04106072E-01;
    COFD[4975] = 1.53462074E-02;
    COFD[4976] = -2.27363668E+01;
    COFD[4977] = 5.11476665E+00;
    COFD[4978] = -3.90017950E-01;
    COFD[4979] = 1.46174143E-02;
    COFD[4980] = -2.28434554E+01;
    COFD[4981] = 5.30874092E+00;
    COFD[4982] = -4.22570517E-01;
    COFD[4983] = 1.63018518E-02;
    COFD[4984] = -2.26645375E+01;
    COFD[4985] = 5.15653775E+00;
    COFD[4986] = -3.96988537E-01;
    COFD[4987] = 1.49769158E-02;
    COFD[4988] = -1.86480131E+01;
    COFD[4989] = 4.45508831E+00;
    COFD[4990] = -3.56068927E-01;
    COFD[4991] = 1.51117769E-02;
    COFD[4992] = -1.79863237E+01;
    COFD[4993] = 4.88559216E+00;
    COFD[4994] = -4.07208104E-01;
    COFD[4995] = 1.71425982E-02;
    COFD[4996] = -1.41186747E+01;
    COFD[4997] = 3.29251025E+00;
    COFD[4998] = -2.12716771E-01;
    COFD[4999] = 9.20613624E-03;
    COFD[5000] = -1.72189894E+01;
    COFD[5001] = 4.10247595E+00;
    COFD[5002] = -3.12590710E-01;
    COFD[5003] = 1.33179236E-02;
    COFD[5004] = -1.87553343E+01;
    COFD[5005] = 4.50343635E+00;
    COFD[5006] = -3.61818471E-01;
    COFD[5007] = 1.53404608E-02;
    COFD[5008] = -1.72409342E+01;
    COFD[5009] = 4.10247595E+00;
    COFD[5010] = -3.12590710E-01;
    COFD[5011] = 1.33179236E-02;
    COFD[5012] = -2.14155602E+01;
    COFD[5013] = 5.24536650E+00;
    COFD[5014] = -4.11756786E-01;
    COFD[5015] = 1.57371608E-02;
    COFD[5016] = -1.84591064E+01;
    COFD[5017] = 4.39561236E+00;
    COFD[5018] = -3.48858889E-01;
    COFD[5019] = 1.48190509E-02;
    COFD[5020] = -2.16482216E+01;
    COFD[5021] = 5.40849808E+00;
    COFD[5022] = -4.63989845E-01;
    COFD[5023] = 1.91734161E-02;
    COFD[5024] = -1.96749452E+01;
    COFD[5025] = 4.87934939E+00;
    COFD[5026] = -4.06539237E-01;
    COFD[5027] = 1.71192321E-02;
    COFD[5028] = -1.96231480E+01;
    COFD[5029] = 4.85846969E+00;
    COFD[5030] = -4.04137136E-01;
    COFD[5031] = 1.70272658E-02;
    COFD[5032] = -1.87640840E+01;
    COFD[5033] = 4.50343635E+00;
    COFD[5034] = -3.61818471E-01;
    COFD[5035] = 1.53404608E-02;
    COFD[5036] = -1.87724574E+01;
    COFD[5037] = 4.50343635E+00;
    COFD[5038] = -3.61818471E-01;
    COFD[5039] = 1.53404608E-02;
    COFD[5040] = -2.20651382E+01;
    COFD[5041] = 5.38704575E+00;
    COFD[5042] = -4.36283370E-01;
    COFD[5043] = 1.70286543E-02;
    COFD[5044] = -2.17225130E+01;
    COFD[5045] = 5.41757449E+00;
    COFD[5046] = -4.64778880E-01;
    COFD[5047] = 1.91926324E-02;
    COFD[5048] = -2.14001824E+01;
    COFD[5049] = 5.39261033E+00;
    COFD[5050] = -4.62555501E-01;
    COFD[5051] = 1.91354360E-02;
    COFD[5052] = -2.17124997E+01;
    COFD[5053] = 5.41757449E+00;
    COFD[5054] = -4.64778880E-01;
    COFD[5055] = 1.91926324E-02;
    COFD[5056] = -2.16232094E+01;
    COFD[5057] = 5.44000785E+00;
    COFD[5058] = -4.65569459E-01;
    COFD[5059] = 1.91449572E-02;
    COFD[5060] = -2.16354957E+01;
    COFD[5061] = 5.44000785E+00;
    COFD[5062] = -4.65569459E-01;
    COFD[5063] = 1.91449572E-02;
    COFD[5064] = -2.21226089E+01;
    COFD[5065] = 5.41240338E+00;
    COFD[5066] = -4.41033923E-01;
    COFD[5067] = 1.72894624E-02;
    COFD[5068] = -2.23100112E+01;
    COFD[5069] = 5.46199540E+00;
    COFD[5070] = -4.51467141E-01;
    COFD[5071] = 1.78907506E-02;
    COFD[5072] = -1.97503281E+01;
    COFD[5073] = 4.91724977E+00;
    COFD[5074] = -4.10632023E-01;
    COFD[5075] = 1.72638615E-02;
    COFD[5076] = -2.23215905E+01;
    COFD[5077] = 5.46199540E+00;
    COFD[5078] = -4.51467141E-01;
    COFD[5079] = 1.78907506E-02;
    COFD[5080] = -2.21467458E+01;
    COFD[5081] = 5.47453801E+00;
    COFD[5082] = -4.63573634E-01;
    COFD[5083] = 1.88106950E-02;
    COFD[5084] = -2.21467458E+01;
    COFD[5085] = 5.47453801E+00;
    COFD[5086] = -4.63573634E-01;
    COFD[5087] = 1.88106950E-02;
    COFD[5088] = -2.20611504E+01;
    COFD[5089] = 5.46751904E+00;
    COFD[5090] = -4.64436117E-01;
    COFD[5091] = 1.89156315E-02;
    COFD[5092] = -2.24444255E+01;
    COFD[5093] = 5.49310818E+00;
    COFD[5094] = -4.62798521E-01;
    COFD[5095] = 1.86559695E-02;
    COFD[5096] = -2.24283088E+01;
    COFD[5097] = 5.49344748E+00;
    COFD[5098] = -4.63052816E-01;
    COFD[5099] = 1.86750389E-02;
    COFD[5100] = -2.24322136E+01;
    COFD[5101] = 5.49344748E+00;
    COFD[5102] = -4.63052816E-01;
    COFD[5103] = 1.86750389E-02;
    COFD[5104] = -2.25598937E+01;
    COFD[5105] = 5.49413184E+00;
    COFD[5106] = -4.63461183E-01;
    COFD[5107] = 1.87049804E-02;
    COFD[5108] = -2.25889692E+01;
    COFD[5109] = 5.46241873E+00;
    COFD[5110] = -4.51563227E-01;
    COFD[5111] = 1.78964292E-02;
    COFD[5112] = -2.25481812E+01;
    COFD[5113] = 5.46266513E+00;
    COFD[5114] = -4.51619963E-01;
    COFD[5115] = 1.78998009E-02;
    COFD[5116] = -2.27255842E+01;
    COFD[5117] = 5.48298766E+00;
    COFD[5118] = -4.57543838E-01;
    COFD[5119] = 1.82786713E-02;
    COFD[5120] = -2.27182958E+01;
    COFD[5121] = 5.48606693E+00;
    COFD[5122] = -4.58921993E-01;
    COFD[5123] = 1.83751829E-02;
    COFD[5124] = -2.26730180E+01;
    COFD[5125] = 5.45838393E+00;
    COFD[5126] = -4.50616109E-01;
    COFD[5127] = 1.78398266E-02;
    COFD[5128] = -2.22224154E+01;
    COFD[5129] = 5.41240338E+00;
    COFD[5130] = -4.41033923E-01;
    COFD[5131] = 1.72894624E-02;
    COFD[5132] = -2.22273939E+01;
    COFD[5133] = 5.41240338E+00;
    COFD[5134] = -4.41033923E-01;
    COFD[5135] = 1.72894624E-02;
    COFD[5136] = -2.24102868E+01;
    COFD[5137] = 5.45651702E+00;
    COFD[5138] = -4.50175370E-01;
    COFD[5139] = 1.78134258E-02;
    COFD[5140] = -2.26076410E+01;
    COFD[5141] = 5.46824702E+00;
    COFD[5142] = -4.52880728E-01;
    COFD[5143] = 1.79742471E-02;
    COFD[5144] = -2.20957576E+01;
    COFD[5145] = 5.46736313E+00;
    COFD[5146] = -4.63533201E-01;
    COFD[5147] = 1.88428681E-02;
    COFD[5148] = -2.21401790E+01;
    COFD[5149] = 5.47453801E+00;
    COFD[5150] = -4.63573634E-01;
    COFD[5151] = 1.88106950E-02;
    COFD[5152] = -2.26656749E+01;
    COFD[5153] = 5.43804144E+00;
    COFD[5154] = -4.45999485E-01;
    COFD[5155] = 1.75662609E-02;
    COFD[5156] = -2.27506175E+01;
    COFD[5157] = 5.44709568E+00;
    COFD[5158] = -4.48069354E-01;
    COFD[5159] = 1.76894758E-02;
    COFD[5160] = -2.24141963E+01;
    COFD[5161] = 5.45651702E+00;
    COFD[5162] = -4.50175370E-01;
    COFD[5163] = 1.78134258E-02;
    COFD[5164] = -2.28313929E+01;
    COFD[5165] = 5.41966149E+00;
    COFD[5166] = -4.42417688E-01;
    COFD[5167] = 1.73660474E-02;
    COFD[5168] = -2.27276855E+01;
    COFD[5169] = 5.40059960E+00;
    COFD[5170] = -4.38830427E-01;
    COFD[5171] = 1.71686343E-02;
    COFD[5172] = -2.29808472E+01;
    COFD[5173] = 5.44447552E+00;
    COFD[5174] = -4.47281890E-01;
    COFD[5175] = 1.76387634E-02;
    COFD[5176] = -2.29793469E+01;
    COFD[5177] = 5.44447552E+00;
    COFD[5178] = -4.47281890E-01;
    COFD[5179] = 1.76387634E-02;
    COFD[5180] = -2.27935663E+01;
    COFD[5181] = 5.26880815E+00;
    COFD[5182] = -4.15708807E-01;
    COFD[5183] = 1.59422297E-02;
    COFD[5184] = -2.28367149E+01;
    COFD[5185] = 5.18065235E+00;
    COFD[5186] = -4.01095736E-01;
    COFD[5187] = 1.51910791E-02;
    COFD[5188] = -2.29126248E+01;
    COFD[5189] = 5.36245215E+00;
    COFD[5190] = -4.31910941E-01;
    COFD[5191] = 1.67950232E-02;
    COFD[5192] = -2.27726207E+01;
    COFD[5193] = 5.22626014E+00;
    COFD[5194] = -4.08564584E-01;
    COFD[5195] = 1.55722943E-02;
    COFD[5196] = -1.84298463E+01;
    COFD[5197] = 4.38777139E+00;
    COFD[5198] = -3.47895358E-01;
    COFD[5199] = 1.47793796E-02;
    COFD[5200] = -1.86950023E+01;
    COFD[5201] = 5.06329467E+00;
    COFD[5202] = -4.27615122E-01;
    COFD[5203] = 1.79228042E-02;
    COFD[5204] = -1.46609678E+01;
    COFD[5205] = 3.42877141E+00;
    COFD[5206] = -2.29550711E-01;
    COFD[5207] = 9.89957140E-03;
    COFD[5208] = -1.80017096E+01;
    COFD[5209] = 4.31633461E+00;
    COFD[5210] = -3.39053638E-01;
    COFD[5211] = 1.44128540E-02;
    COFD[5212] = -1.95280746E+01;
    COFD[5213] = 4.70710890E+00;
    COFD[5214] = -3.86453789E-01;
    COFD[5215] = 1.63380393E-02;
    COFD[5216] = -1.80251183E+01;
    COFD[5217] = 4.31633461E+00;
    COFD[5218] = -3.39053638E-01;
    COFD[5219] = 1.44128540E-02;
    COFD[5220] = -2.12712683E+01;
    COFD[5221] = 5.04351440E+00;
    COFD[5222] = -3.78432548E-01;
    COFD[5223] = 1.40280477E-02;
    COFD[5224] = -1.91368159E+01;
    COFD[5225] = 4.56508921E+00;
    COFD[5226] = -3.69033080E-01;
    COFD[5227] = 1.56223865E-02;
    COFD[5228] = -2.20785179E+01;
    COFD[5229] = 5.45354415E+00;
    COFD[5230] = -4.65291244E-01;
    COFD[5231] = 1.90545098E-02;
    COFD[5232] = -2.03731422E+01;
    COFD[5233] = 5.05714666E+00;
    COFD[5234] = -4.26960034E-01;
    COFD[5235] = 1.78999526E-02;
    COFD[5236] = -2.03251784E+01;
    COFD[5237] = 5.03716590E+00;
    COFD[5238] = -4.24714370E-01;
    COFD[5239] = 1.78162251E-02;
    COFD[5240] = -1.95377616E+01;
    COFD[5241] = 4.70710890E+00;
    COFD[5242] = -3.86453789E-01;
    COFD[5243] = 1.63380393E-02;
    COFD[5244] = -1.95470492E+01;
    COFD[5245] = 4.70710890E+00;
    COFD[5246] = -3.86453789E-01;
    COFD[5247] = 1.63380393E-02;
    COFD[5248] = -2.20279970E+01;
    COFD[5249] = 5.23448541E+00;
    COFD[5250] = -4.09927386E-01;
    COFD[5251] = 1.56423231E-02;
    COFD[5252] = -2.21190597E+01;
    COFD[5253] = 5.45729646E+00;
    COFD[5254] = -4.65340425E-01;
    COFD[5255] = 1.90398312E-02;
    COFD[5256] = -2.18452486E+01;
    COFD[5257] = 5.44764923E+00;
    COFD[5258] = -4.65239944E-01;
    COFD[5259] = 1.90796594E-02;
    COFD[5260] = -2.21080365E+01;
    COFD[5261] = 5.45729646E+00;
    COFD[5262] = -4.65340425E-01;
    COFD[5263] = 1.90398312E-02;
    COFD[5264] = -2.19960701E+01;
    COFD[5265] = 5.46751578E+00;
    COFD[5266] = -4.64445737E-01;
    COFD[5267] = 1.89164296E-02;
    COFD[5268] = -2.20094808E+01;
    COFD[5269] = 5.46751578E+00;
    COFD[5270] = -4.64445737E-01;
    COFD[5271] = 1.89164296E-02;
    COFD[5272] = -2.21202988E+01;
    COFD[5273] = 5.27340167E+00;
    COFD[5274] = -4.16479301E-01;
    COFD[5275] = 1.59820820E-02;
    COFD[5276] = -2.24302648E+01;
    COFD[5277] = 5.37317627E+00;
    COFD[5278] = -4.33806826E-01;
    COFD[5279] = 1.68960241E-02;
    COFD[5280] = -2.04982298E+01;
    COFD[5281] = 5.09287074E+00;
    COFD[5282] = -4.30747337E-01;
    COFD[5283] = 1.80311502E-02;
    COFD[5284] = -2.24433083E+01;
    COFD[5285] = 5.37317627E+00;
    COFD[5286] = -4.33806826E-01;
    COFD[5287] = 1.68960241E-02;
    COFD[5288] = -2.25051541E+01;
    COFD[5289] = 5.48975794E+00;
    COFD[5290] = -4.60336840E-01;
    COFD[5291] = 1.84719899E-02;
    COFD[5292] = -2.25051541E+01;
    COFD[5293] = 5.48975794E+00;
    COFD[5294] = -4.60336840E-01;
    COFD[5295] = 1.84719899E-02;
    COFD[5296] = -2.24444255E+01;
    COFD[5297] = 5.49310818E+00;
    COFD[5298] = -4.62798521E-01;
    COFD[5299] = 1.86559695E-02;
    COFD[5300] = -2.27303053E+01;
    COFD[5301] = 5.47648899E+00;
    COFD[5302] = -4.55190091E-01;
    COFD[5303] = 1.81203119E-02;
    COFD[5304] = -2.27174313E+01;
    COFD[5305] = 5.47742424E+00;
    COFD[5306] = -4.55509593E-01;
    COFD[5307] = 1.81415342E-02;
    COFD[5308] = -2.27219009E+01;
    COFD[5309] = 5.47742424E+00;
    COFD[5310] = -4.55509593E-01;
    COFD[5311] = 1.81415342E-02;
    COFD[5312] = -2.28449688E+01;
    COFD[5313] = 5.47905347E+00;
    COFD[5314] = -4.56023683E-01;
    COFD[5315] = 1.81750993E-02;
    COFD[5316] = -2.27070219E+01;
    COFD[5317] = 5.37425628E+00;
    COFD[5318] = -4.33997939E-01;
    COFD[5319] = 1.69062130E-02;
    COFD[5320] = -2.26704052E+01;
    COFD[5321] = 5.37490349E+00;
    COFD[5322] = -4.34112523E-01;
    COFD[5323] = 1.69123235E-02;
    COFD[5324] = -2.29510526E+01;
    COFD[5325] = 5.44399037E+00;
    COFD[5326] = -4.47161064E-01;
    COFD[5327] = 1.76313359E-02;
    COFD[5328] = -2.29571084E+01;
    COFD[5329] = 5.45255106E+00;
    COFD[5330] = -4.49274880E-01;
    COFD[5331] = 1.77601485E-02;
    COFD[5332] = -2.27791634E+01;
    COFD[5333] = 5.36302327E+00;
    COFD[5334] = -4.32012062E-01;
    COFD[5335] = 1.68004123E-02;
    COFD[5336] = -2.22319743E+01;
    COFD[5337] = 5.27340167E+00;
    COFD[5338] = -4.16479301E-01;
    COFD[5339] = 1.59820820E-02;
    COFD[5340] = -2.22376174E+01;
    COFD[5341] = 5.27340167E+00;
    COFD[5342] = -4.16479301E-01;
    COFD[5343] = 1.59820820E-02;
    COFD[5344] = -2.25142044E+01;
    COFD[5345] = 5.35783672E+00;
    COFD[5346] = -4.31095103E-01;
    COFD[5347] = 1.67515802E-02;
    COFD[5348] = -2.27416253E+01;
    COFD[5349] = 5.38802436E+00;
    COFD[5350] = -4.36467632E-01;
    COFD[5351] = 1.70387872E-02;
    COFD[5352] = -2.24743159E+01;
    COFD[5353] = 5.49090532E+00;
    COFD[5354] = -4.61515055E-01;
    COFD[5355] = 1.85622624E-02;
    COFD[5356] = -2.24977955E+01;
    COFD[5357] = 5.48975794E+00;
    COFD[5358] = -4.60336840E-01;
    COFD[5359] = 1.84719899E-02;
    COFD[5360] = -2.26989068E+01;
    COFD[5361] = 5.31130174E+00;
    COFD[5362] = -4.23012509E-01;
    COFD[5363] = 1.63250839E-02;
    COFD[5364] = -2.28127150E+01;
    COFD[5365] = 5.33452791E+00;
    COFD[5366] = -4.27009745E-01;
    COFD[5367] = 1.65349406E-02;
    COFD[5368] = -2.25186794E+01;
    COFD[5369] = 5.35783672E+00;
    COFD[5370] = -4.31095103E-01;
    COFD[5371] = 1.67515802E-02;
    COFD[5372] = -2.28386812E+01;
    COFD[5373] = 5.28428478E+00;
    COFD[5374] = -4.18334137E-01;
    COFD[5375] = 1.60789015E-02;
    COFD[5376] = -2.27156370E+01;
    COFD[5377] = 5.25593075E+00;
    COFD[5378] = -4.13540252E-01;
    COFD[5379] = 1.58298021E-02;
    COFD[5380] = -2.30198383E+01;
    COFD[5381] = 5.32186706E+00;
    COFD[5382] = -4.24826256E-01;
    COFD[5383] = 1.64201917E-02;
    COFD[5384] = -2.30180603E+01;
    COFD[5385] = 5.32186706E+00;
    COFD[5386] = -4.24826256E-01;
    COFD[5387] = 1.64201917E-02;
    COFD[5388] = -2.26694016E+01;
    COFD[5389] = 5.07552523E+00;
    COFD[5390] = -3.83609779E-01;
    COFD[5391] = 1.42906252E-02;
    COFD[5392] = -2.26349181E+01;
    COFD[5393] = 4.95583598E+00;
    COFD[5394] = -3.64620833E-01;
    COFD[5395] = 1.33399033E-02;
    COFD[5396] = -2.28440138E+01;
    COFD[5397] = 5.19534599E+00;
    COFD[5398] = -4.03564468E-01;
    COFD[5399] = 1.53188143E-02;
    COFD[5400] = -2.26098548E+01;
    COFD[5401] = 5.01607052E+00;
    COFD[5402] = -3.74092059E-01;
    COFD[5403] = 1.38112395E-02;
    COFD[5404] = -1.91111450E+01;
    COFD[5405] = 4.55863849E+00;
    COFD[5406] = -3.68274023E-01;
    COFD[5407] = 1.55925459E-02;
    COFD[5408] = -1.86585461E+01;
    COFD[5409] = 5.05830652E+00;
    COFD[5410] = -4.27083712E-01;
    COFD[5411] = 1.79042712E-02;
    COFD[5412] = -1.46294931E+01;
    COFD[5413] = 3.42403711E+00;
    COFD[5414] = -2.28974532E-01;
    COFD[5415] = 9.87621386E-03;
    COFD[5416] = -1.79605052E+01;
    COFD[5417] = 4.30733766E+00;
    COFD[5418] = -3.37925627E-01;
    COFD[5419] = 1.43655069E-02;
    COFD[5420] = -1.94909629E+01;
    COFD[5421] = 4.69879097E+00;
    COFD[5422] = -3.85441745E-01;
    COFD[5423] = 1.62968073E-02;
    COFD[5424] = -1.79840143E+01;
    COFD[5425] = 4.30733766E+00;
    COFD[5426] = -3.37925627E-01;
    COFD[5427] = 1.43655069E-02;
    COFD[5428] = -2.12710109E+01;
    COFD[5429] = 5.05274021E+00;
    COFD[5430] = -3.79920248E-01;
    COFD[5431] = 1.41033861E-02;
    COFD[5432] = -1.91050585E+01;
    COFD[5433] = 4.55886937E+00;
    COFD[5434] = -3.68301178E-01;
    COFD[5435] = 1.55936129E-02;
    COFD[5436] = -2.20577824E+01;
    COFD[5437] = 5.45198096E+00;
    COFD[5438] = -4.65259344E-01;
    COFD[5439] = 1.90596774E-02;
    COFD[5440] = -2.03431655E+01;
    COFD[5441] = 5.05232294E+00;
    COFD[5442] = -4.26452915E-01;
    COFD[5443] = 1.78825996E-02;
    COFD[5444] = -2.02895141E+01;
    COFD[5445] = 5.02992150E+00;
    COFD[5446] = -4.23868677E-01;
    COFD[5447] = 1.77832701E-02;
    COFD[5448] = -1.95007169E+01;
    COFD[5449] = 4.69879097E+00;
    COFD[5450] = -3.85441745E-01;
    COFD[5451] = 1.62968073E-02;
    COFD[5452] = -1.95100699E+01;
    COFD[5453] = 4.69879097E+00;
    COFD[5454] = -3.85441745E-01;
    COFD[5455] = 1.62968073E-02;
    COFD[5456] = -2.20252840E+01;
    COFD[5457] = 5.24121166E+00;
    COFD[5458] = -4.11059025E-01;
    COFD[5459] = 1.57010104E-02;
    COFD[5460] = -2.20987680E+01;
    COFD[5461] = 5.45582845E+00;
    COFD[5462] = -4.65322899E-01;
    COFD[5463] = 1.90457166E-02;
    COFD[5464] = -2.18235959E+01;
    COFD[5465] = 5.44634823E+00;
    COFD[5466] = -4.65245448E-01;
    COFD[5467] = 1.90865967E-02;
    COFD[5468] = -2.20876730E+01;
    COFD[5469] = 5.45582845E+00;
    COFD[5470] = -4.65322899E-01;
    COFD[5471] = 1.90457166E-02;
    COFD[5472] = -2.19775435E+01;
    COFD[5473] = 5.46745929E+00;
    COFD[5474] = -4.64624358E-01;
    COFD[5475] = 1.89312280E-02;
    COFD[5476] = -2.19910334E+01;
    COFD[5477] = 5.46745929E+00;
    COFD[5478] = -4.64624358E-01;
    COFD[5479] = 1.89312280E-02;
    COFD[5480] = -2.21169561E+01;
    COFD[5481] = 5.27973457E+00;
    COFD[5482] = -4.17551922E-01;
    COFD[5483] = 1.60378887E-02;
    COFD[5484] = -2.24254137E+01;
    COFD[5485] = 5.37839633E+00;
    COFD[5486] = -4.34732804E-01;
    COFD[5487] = 1.69454459E-02;
    COFD[5488] = -2.04673093E+01;
    COFD[5489] = 5.08807429E+00;
    COFD[5490] = -4.30243518E-01;
    COFD[5491] = 1.80139269E-02;
    COFD[5492] = -2.24385639E+01;
    COFD[5493] = 5.37839633E+00;
    COFD[5494] = -4.34732804E-01;
    COFD[5495] = 1.69454459E-02;
    COFD[5496] = -2.24900659E+01;
    COFD[5497] = 5.49046566E+00;
    COFD[5498] = -4.60640519E-01;
    COFD[5499] = 1.84931917E-02;
    COFD[5500] = -2.24900659E+01;
    COFD[5501] = 5.49046566E+00;
    COFD[5502] = -4.60640519E-01;
    COFD[5503] = 1.84931917E-02;
    COFD[5504] = -2.24283088E+01;
    COFD[5505] = 5.49344748E+00;
    COFD[5506] = -4.63052816E-01;
    COFD[5507] = 1.86750389E-02;
    COFD[5508] = -2.27174313E+01;
    COFD[5509] = 5.47742424E+00;
    COFD[5510] = -4.55509593E-01;
    COFD[5511] = 1.81415342E-02;
    COFD[5512] = -2.27048016E+01;
    COFD[5513] = 5.47849754E+00;
    COFD[5514] = -4.55849160E-01;
    COFD[5515] = 1.81637187E-02;
    COFD[5516] = -2.27093131E+01;
    COFD[5517] = 5.47849754E+00;
    COFD[5518] = -4.55849160E-01;
    COFD[5519] = 1.81637187E-02;
    COFD[5520] = -2.28331297E+01;
    COFD[5521] = 5.48004688E+00;
    COFD[5522] = -4.56351694E-01;
    COFD[5523] = 1.81967269E-02;
    COFD[5524] = -2.27043153E+01;
    COFD[5525] = 5.37944812E+00;
    COFD[5526] = -4.34919903E-01;
    COFD[5527] = 1.69554449E-02;
    COFD[5528] = -2.26673882E+01;
    COFD[5529] = 5.38007821E+00;
    COFD[5530] = -4.35032051E-01;
    COFD[5531] = 1.69614401E-02;
    COFD[5532] = -2.29382506E+01;
    COFD[5533] = 5.44442303E+00;
    COFD[5534] = -4.47427940E-01;
    COFD[5535] = 1.76507657E-02;
    COFD[5536] = -2.29489166E+01;
    COFD[5537] = 5.45471156E+00;
    COFD[5538] = -4.49760531E-01;
    COFD[5539] = 1.77887865E-02;
    COFD[5540] = -2.27773775E+01;
    COFD[5541] = 5.36826093E+00;
    COFD[5542] = -4.32938535E-01;
    COFD[5543] = 1.68497727E-02;
    COFD[5544] = -2.22294897E+01;
    COFD[5545] = 5.27973457E+00;
    COFD[5546] = -4.17551922E-01;
    COFD[5547] = 1.60378887E-02;
    COFD[5548] = -2.22351815E+01;
    COFD[5549] = 5.27973457E+00;
    COFD[5550] = -4.17551922E-01;
    COFD[5551] = 1.60378887E-02;
    COFD[5552] = -2.25103101E+01;
    COFD[5553] = 5.36306078E+00;
    COFD[5554] = -4.32018705E-01;
    COFD[5555] = 1.68007664E-02;
    COFD[5556] = -2.27363406E+01;
    COFD[5557] = 5.39211717E+00;
    COFD[5558] = -4.37241388E-01;
    COFD[5559] = 1.70814168E-02;
    COFD[5560] = -2.24568662E+01;
    COFD[5561] = 5.49068313E+00;
    COFD[5562] = -4.61694593E-01;
    COFD[5563] = 1.85779923E-02;
    COFD[5564] = -2.24826500E+01;
    COFD[5565] = 5.49046566E+00;
    COFD[5566] = -4.60640519E-01;
    COFD[5567] = 1.84931917E-02;
    COFD[5568] = -2.26989595E+01;
    COFD[5569] = 5.31729293E+00;
    COFD[5570] = -4.24040190E-01;
    COFD[5571] = 1.63789533E-02;
    COFD[5572] = -2.28126636E+01;
    COFD[5573] = 5.34025709E+00;
    COFD[5574] = -4.28002604E-01;
    COFD[5575] = 1.65872781E-02;
    COFD[5576] = -2.25148268E+01;
    COFD[5577] = 5.36306078E+00;
    COFD[5578] = -4.32018705E-01;
    COFD[5579] = 1.68007664E-02;
    COFD[5580] = -2.28395602E+01;
    COFD[5581] = 5.29000716E+00;
    COFD[5582] = -4.19320503E-01;
    COFD[5583] = 1.61306865E-02;
    COFD[5584] = -2.27176356E+01;
    COFD[5585] = 5.26237767E+00;
    COFD[5586] = -4.14629620E-01;
    COFD[5587] = 1.58863936E-02;
    COFD[5588] = -2.30219288E+01;
    COFD[5589] = 5.32771617E+00;
    COFD[5590] = -4.25833456E-01;
    COFD[5591] = 1.64730816E-02;
    COFD[5592] = -2.30201295E+01;
    COFD[5593] = 5.32771617E+00;
    COFD[5594] = -4.25833456E-01;
    COFD[5595] = 1.64730816E-02;
    COFD[5596] = -2.26785584E+01;
    COFD[5597] = 5.08419101E+00;
    COFD[5598] = -3.85017801E-01;
    COFD[5599] = 1.43622302E-02;
    COFD[5600] = -2.26506465E+01;
    COFD[5601] = 4.96659745E+00;
    COFD[5602] = -3.66306172E-01;
    COFD[5603] = 1.34235902E-02;
    COFD[5604] = -2.28532517E+01;
    COFD[5605] = 5.20408596E+00;
    COFD[5606] = -4.04973723E-01;
    COFD[5607] = 1.53901363E-02;
    COFD[5608] = -2.26226017E+01;
    COFD[5609] = 5.02598120E+00;
    COFD[5610] = -3.75659093E-01;
    COFD[5611] = 1.38894769E-02;
    COFD[5612] = -1.90794464E+01;
    COFD[5613] = 4.55246544E+00;
    COFD[5614] = -3.67548311E-01;
    COFD[5615] = 1.55640443E-02;
    COFD[5616] = -1.86587075E+01;
    COFD[5617] = 5.05830652E+00;
    COFD[5618] = -4.27083712E-01;
    COFD[5619] = 1.79042712E-02;
    COFD[5620] = -1.46298102E+01;
    COFD[5621] = 3.42403711E+00;
    COFD[5622] = -2.28974532E-01;
    COFD[5623] = 9.87621386E-03;
    COFD[5624] = -1.79625305E+01;
    COFD[5625] = 4.30733766E+00;
    COFD[5626] = -3.37925627E-01;
    COFD[5627] = 1.43655069E-02;
    COFD[5628] = -1.94942738E+01;
    COFD[5629] = 4.69879097E+00;
    COFD[5630] = -3.85441745E-01;
    COFD[5631] = 1.62968073E-02;
    COFD[5632] = -1.79861374E+01;
    COFD[5633] = 4.30733766E+00;
    COFD[5634] = -3.37925627E-01;
    COFD[5635] = 1.43655069E-02;
    COFD[5636] = -2.12732290E+01;
    COFD[5637] = 5.05274021E+00;
    COFD[5638] = -3.79920248E-01;
    COFD[5639] = 1.41033861E-02;
    COFD[5640] = -1.91080950E+01;
    COFD[5641] = 4.55886937E+00;
    COFD[5642] = -3.68301178E-01;
    COFD[5643] = 1.55936129E-02;
    COFD[5644] = -2.20617871E+01;
    COFD[5645] = 5.45198096E+00;
    COFD[5646] = -4.65259344E-01;
    COFD[5647] = 1.90596774E-02;
    COFD[5648] = -2.03450948E+01;
    COFD[5649] = 5.05232294E+00;
    COFD[5650] = -4.26452915E-01;
    COFD[5651] = 1.78825996E-02;
    COFD[5652] = -2.02915437E+01;
    COFD[5653] = 5.02992150E+00;
    COFD[5654] = -4.23868677E-01;
    COFD[5655] = 1.77832701E-02;
    COFD[5656] = -1.95040932E+01;
    COFD[5657] = 4.69879097E+00;
    COFD[5658] = -3.85441745E-01;
    COFD[5659] = 1.62968073E-02;
    COFD[5660] = -1.95135102E+01;
    COFD[5661] = 4.69879097E+00;
    COFD[5662] = -3.85441745E-01;
    COFD[5663] = 1.62968073E-02;
    COFD[5664] = -2.20284624E+01;
    COFD[5665] = 5.24121166E+00;
    COFD[5666] = -4.11059025E-01;
    COFD[5667] = 1.57010104E-02;
    COFD[5668] = -2.21019494E+01;
    COFD[5669] = 5.45582845E+00;
    COFD[5670] = -4.65322899E-01;
    COFD[5671] = 1.90457166E-02;
    COFD[5672] = -2.18266355E+01;
    COFD[5673] = 5.44634823E+00;
    COFD[5674] = -4.65245448E-01;
    COFD[5675] = 1.90865967E-02;
    COFD[5676] = -2.20907844E+01;
    COFD[5677] = 5.45582845E+00;
    COFD[5678] = -4.65322899E-01;
    COFD[5679] = 1.90457166E-02;
    COFD[5680] = -2.19803572E+01;
    COFD[5681] = 5.46745929E+00;
    COFD[5682] = -4.64624358E-01;
    COFD[5683] = 1.89312280E-02;
    COFD[5684] = -2.19939242E+01;
    COFD[5685] = 5.46745929E+00;
    COFD[5686] = -4.64624358E-01;
    COFD[5687] = 1.89312280E-02;
    COFD[5688] = -2.21202699E+01;
    COFD[5689] = 5.27973457E+00;
    COFD[5690] = -4.17551922E-01;
    COFD[5691] = 1.60378887E-02;
    COFD[5692] = -2.24293161E+01;
    COFD[5693] = 5.37839633E+00;
    COFD[5694] = -4.34732804E-01;
    COFD[5695] = 1.69454459E-02;
    COFD[5696] = -2.04711578E+01;
    COFD[5697] = 5.08807429E+00;
    COFD[5698] = -4.30243518E-01;
    COFD[5699] = 1.80139269E-02;
    COFD[5700] = -2.24425707E+01;
    COFD[5701] = 5.37839633E+00;
    COFD[5702] = -4.34732804E-01;
    COFD[5703] = 1.69454459E-02;
    COFD[5704] = -2.24938619E+01;
    COFD[5705] = 5.49046566E+00;
    COFD[5706] = -4.60640519E-01;
    COFD[5707] = 1.84931917E-02;
    COFD[5708] = -2.24938619E+01;
    COFD[5709] = 5.49046566E+00;
    COFD[5710] = -4.60640519E-01;
    COFD[5711] = 1.84931917E-02;
    COFD[5712] = -2.24322136E+01;
    COFD[5713] = 5.49344748E+00;
    COFD[5714] = -4.63052816E-01;
    COFD[5715] = 1.86750389E-02;
    COFD[5716] = -2.27219009E+01;
    COFD[5717] = 5.47742424E+00;
    COFD[5718] = -4.55509593E-01;
    COFD[5719] = 1.81415342E-02;
    COFD[5720] = -2.27093131E+01;
    COFD[5721] = 5.47849754E+00;
    COFD[5722] = -4.55849160E-01;
    COFD[5723] = 1.81637187E-02;
    COFD[5724] = -2.27138656E+01;
    COFD[5725] = 5.47849754E+00;
    COFD[5726] = -4.55849160E-01;
    COFD[5727] = 1.81637187E-02;
    COFD[5728] = -2.28377226E+01;
    COFD[5729] = 5.48004688E+00;
    COFD[5730] = -4.56351694E-01;
    COFD[5731] = 1.81967269E-02;
    COFD[5732] = -2.27089064E+01;
    COFD[5733] = 5.37944812E+00;
    COFD[5734] = -4.34919903E-01;
    COFD[5735] = 1.69554449E-02;
    COFD[5736] = -2.26720190E+01;
    COFD[5737] = 5.38007821E+00;
    COFD[5738] = -4.35032051E-01;
    COFD[5739] = 1.69614401E-02;
    COFD[5740] = -2.29432741E+01;
    COFD[5741] = 5.44442303E+00;
    COFD[5742] = -4.47427940E-01;
    COFD[5743] = 1.76507657E-02;
    COFD[5744] = -2.29539724E+01;
    COFD[5745] = 5.45471156E+00;
    COFD[5746] = -4.49760531E-01;
    COFD[5747] = 1.77887865E-02;
    COFD[5748] = -2.27824652E+01;
    COFD[5749] = 5.36826093E+00;
    COFD[5750] = -4.32938535E-01;
    COFD[5751] = 1.68497727E-02;
    COFD[5752] = -2.22336433E+01;
    COFD[5753] = 5.27973457E+00;
    COFD[5754] = -4.17551922E-01;
    COFD[5755] = 1.60378887E-02;
    COFD[5756] = -2.22393829E+01;
    COFD[5757] = 5.27973457E+00;
    COFD[5758] = -4.17551922E-01;
    COFD[5759] = 1.60378887E-02;
    COFD[5760] = -2.25148197E+01;
    COFD[5761] = 5.36306078E+00;
    COFD[5762] = -4.32018705E-01;
    COFD[5763] = 1.68007664E-02;
    COFD[5764] = -2.27408914E+01;
    COFD[5765] = 5.39211717E+00;
    COFD[5766] = -4.37241388E-01;
    COFD[5767] = 1.70814168E-02;
    COFD[5768] = -2.24607172E+01;
    COFD[5769] = 5.49068313E+00;
    COFD[5770] = -4.61694593E-01;
    COFD[5771] = 1.85779923E-02;
    COFD[5772] = -2.24863899E+01;
    COFD[5773] = 5.49046566E+00;
    COFD[5774] = -4.60640519E-01;
    COFD[5775] = 1.84931917E-02;
    COFD[5776] = -2.27040773E+01;
    COFD[5777] = 5.31729293E+00;
    COFD[5778] = -4.24040190E-01;
    COFD[5779] = 1.63789533E-02;
    COFD[5780] = -2.28177500E+01;
    COFD[5781] = 5.34025709E+00;
    COFD[5782] = -4.28002604E-01;
    COFD[5783] = 1.65872781E-02;
    COFD[5784] = -2.25193776E+01;
    COFD[5785] = 5.36306078E+00;
    COFD[5786] = -4.32018705E-01;
    COFD[5787] = 1.68007664E-02;
    COFD[5788] = -2.28450684E+01;
    COFD[5789] = 5.29000716E+00;
    COFD[5790] = -4.19320503E-01;
    COFD[5791] = 1.61306865E-02;
    COFD[5792] = -2.27231185E+01;
    COFD[5793] = 5.26237767E+00;
    COFD[5794] = -4.14629620E-01;
    COFD[5795] = 1.58863936E-02;
    COFD[5796] = -2.30277582E+01;
    COFD[5797] = 5.32771617E+00;
    COFD[5798] = -4.25833456E-01;
    COFD[5799] = 1.64730816E-02;
    COFD[5800] = -2.30259379E+01;
    COFD[5801] = 5.32771617E+00;
    COFD[5802] = -4.25833456E-01;
    COFD[5803] = 1.64730816E-02;
    COFD[5804] = -2.26849244E+01;
    COFD[5805] = 5.08419101E+00;
    COFD[5806] = -3.85017801E-01;
    COFD[5807] = 1.43622302E-02;
    COFD[5808] = -2.26574071E+01;
    COFD[5809] = 4.96659745E+00;
    COFD[5810] = -3.66306172E-01;
    COFD[5811] = 1.34235902E-02;
    COFD[5812] = -2.28593475E+01;
    COFD[5813] = 5.20408596E+00;
    COFD[5814] = -4.04973723E-01;
    COFD[5815] = 1.53901363E-02;
    COFD[5816] = -2.26291682E+01;
    COFD[5817] = 5.02598120E+00;
    COFD[5818] = -3.75659093E-01;
    COFD[5819] = 1.38894769E-02;
    COFD[5820] = -1.90824831E+01;
    COFD[5821] = 4.55246544E+00;
    COFD[5822] = -3.67548311E-01;
    COFD[5823] = 1.55640443E-02;
    COFD[5824] = -1.88039398E+01;
    COFD[5825] = 5.05088527E+00;
    COFD[5826] = -4.26300667E-01;
    COFD[5827] = 1.78773589E-02;
    COFD[5828] = -1.47601286E+01;
    COFD[5829] = 3.41694605E+00;
    COFD[5830] = -2.28112628E-01;
    COFD[5831] = 9.84132390E-03;
    COFD[5832] = -1.80806236E+01;
    COFD[5833] = 4.29360798E+00;
    COFD[5834] = -3.36202347E-01;
    COFD[5835] = 1.42930967E-02;
    COFD[5836] = -1.96030292E+01;
    COFD[5837] = 4.68635682E+00;
    COFD[5838] = -3.83931845E-01;
    COFD[5839] = 1.62354222E-02;
    COFD[5840] = -1.81043255E+01;
    COFD[5841] = 4.29360798E+00;
    COFD[5842] = -3.36202347E-01;
    COFD[5843] = 1.42930967E-02;
    COFD[5844] = -2.14557574E+01;
    COFD[5845] = 5.06615192E+00;
    COFD[5846] = -3.82089787E-01;
    COFD[5847] = 1.42134191E-02;
    COFD[5848] = -1.92212680E+01;
    COFD[5849] = 4.54955804E+00;
    COFD[5850] = -3.67206725E-01;
    COFD[5851] = 1.55506376E-02;
    COFD[5852] = -2.21891699E+01;
    COFD[5853] = 5.44981030E+00;
    COFD[5854] = -4.65238257E-01;
    COFD[5855] = 1.90687553E-02;
    COFD[5856] = -2.04523945E+01;
    COFD[5857] = 5.04244535E+00;
    COFD[5858] = -4.25330213E-01;
    COFD[5859] = 1.78402047E-02;
    COFD[5860] = -2.03971631E+01;
    COFD[5861] = 5.01889159E+00;
    COFD[5862] = -4.22579631E-01;
    COFD[5863] = 1.77329812E-02;
    COFD[5864] = -1.96129127E+01;
    COFD[5865] = 4.68635682E+00;
    COFD[5866] = -3.83931845E-01;
    COFD[5867] = 1.62354222E-02;
    COFD[5868] = -1.96223921E+01;
    COFD[5869] = 4.68635682E+00;
    COFD[5870] = -3.83931845E-01;
    COFD[5871] = 1.62354222E-02;
    COFD[5872] = -2.21864158E+01;
    COFD[5873] = 5.25139402E+00;
    COFD[5874] = -4.12771525E-01;
    COFD[5875] = 1.57897958E-02;
    COFD[5876] = -2.22197724E+01;
    COFD[5877] = 5.45358272E+00;
    COFD[5878] = -4.65292061E-01;
    COFD[5879] = 1.90543847E-02;
    COFD[5880] = -2.19586673E+01;
    COFD[5881] = 5.44474614E+00;
    COFD[5882] = -4.65306860E-01;
    COFD[5883] = 1.90996430E-02;
    COFD[5884] = -2.22085390E+01;
    COFD[5885] = 5.45358272E+00;
    COFD[5886] = -4.65292061E-01;
    COFD[5887] = 1.90543847E-02;
    COFD[5888] = -2.21121170E+01;
    COFD[5889] = 5.46738366E+00;
    COFD[5890] = -4.64895439E-01;
    COFD[5891] = 1.89536362E-02;
    COFD[5892] = -2.21257593E+01;
    COFD[5893] = 5.46738366E+00;
    COFD[5894] = -4.64895439E-01;
    COFD[5895] = 1.89536362E-02;
    COFD[5896] = -2.22744971E+01;
    COFD[5897] = 5.28846317E+00;
    COFD[5898] = -4.19054070E-01;
    COFD[5899] = 1.61166909E-02;
    COFD[5900] = -2.25761948E+01;
    COFD[5901] = 5.38613682E+00;
    COFD[5902] = -4.36112477E-01;
    COFD[5903] = 1.70192630E-02;
    COFD[5904] = -2.06099144E+01;
    COFD[5905] = 5.08069320E+00;
    COFD[5906] = -4.29462569E-01;
    COFD[5907] = 1.79869461E-02;
    COFD[5908] = -2.25895516E+01;
    COFD[5909] = 5.38613682E+00;
    COFD[5910] = -4.36112477E-01;
    COFD[5911] = 1.70192630E-02;
    COFD[5912] = -2.26196129E+01;
    COFD[5913] = 5.49130890E+00;
    COFD[5914] = -4.61073466E-01;
    COFD[5915] = 1.85242949E-02;
    COFD[5916] = -2.26196129E+01;
    COFD[5917] = 5.49130890E+00;
    COFD[5918] = -4.61073466E-01;
    COFD[5919] = 1.85242949E-02;
    COFD[5920] = -2.25598937E+01;
    COFD[5921] = 5.49413184E+00;
    COFD[5922] = -4.63461183E-01;
    COFD[5923] = 1.87049804E-02;
    COFD[5924] = -2.28449688E+01;
    COFD[5925] = 5.47905347E+00;
    COFD[5926] = -4.56023683E-01;
    COFD[5927] = 1.81750993E-02;
    COFD[5928] = -2.28331297E+01;
    COFD[5929] = 5.48004688E+00;
    COFD[5930] = -4.56351694E-01;
    COFD[5931] = 1.81967269E-02;
    COFD[5932] = -2.28377226E+01;
    COFD[5933] = 5.48004688E+00;
    COFD[5934] = -4.56351694E-01;
    COFD[5935] = 1.81967269E-02;
    COFD[5936] = -2.29541620E+01;
    COFD[5937] = 5.48145821E+00;
    COFD[5938] = -4.56834882E-01;
    COFD[5939] = 1.82288359E-02;
    COFD[5940] = -2.28439474E+01;
    COFD[5941] = 5.38698965E+00;
    COFD[5942] = -4.36272815E-01;
    COFD[5943] = 1.70280741E-02;
    COFD[5944] = -2.28092883E+01;
    COFD[5945] = 5.38749986E+00;
    COFD[5946] = -4.36368839E-01;
    COFD[5947] = 1.70333535E-02;
    COFD[5948] = -2.30611542E+01;
    COFD[5949] = 5.44729879E+00;
    COFD[5950] = -4.48113754E-01;
    COFD[5951] = 1.76920689E-02;
    COFD[5952] = -2.30719702E+01;
    COFD[5953] = 5.45769390E+00;
    COFD[5954] = -4.50452338E-01;
    COFD[5955] = 1.78300009E-02;
    COFD[5956] = -2.29169419E+01;
    COFD[5957] = 5.37626442E+00;
    COFD[5958] = -4.34353973E-01;
    COFD[5959] = 1.69252109E-02;
    COFD[5960] = -2.23886928E+01;
    COFD[5961] = 5.28846317E+00;
    COFD[5962] = -4.19054070E-01;
    COFD[5963] = 1.61166909E-02;
    COFD[5964] = -2.23944792E+01;
    COFD[5965] = 5.28846317E+00;
    COFD[5966] = -4.19054070E-01;
    COFD[5967] = 1.61166909E-02;
    COFD[5968] = -2.26606931E+01;
    COFD[5969] = 5.37102704E+00;
    COFD[5970] = -4.33426866E-01;
    COFD[5971] = 1.68757765E-02;
    COFD[5972] = -2.28720058E+01;
    COFD[5973] = 5.39857883E+00;
    COFD[5974] = -4.38452012E-01;
    COFD[5975] = 1.71478655E-02;
    COFD[5976] = -2.25843967E+01;
    COFD[5977] = 5.49024820E+00;
    COFD[5978] = -4.61953415E-01;
    COFD[5979] = 1.86012452E-02;
    COFD[5980] = -2.26120860E+01;
    COFD[5981] = 5.49130890E+00;
    COFD[5982] = -4.61073466E-01;
    COFD[5983] = 1.85242949E-02;
    COFD[5984] = -2.28412555E+01;
    COFD[5985] = 5.32619107E+00;
    COFD[5986] = -4.25570588E-01;
    COFD[5987] = 1.64592712E-02;
    COFD[5988] = -2.29487149E+01;
    COFD[5989] = 5.34854008E+00;
    COFD[5990] = -4.29453170E-01;
    COFD[5991] = 1.66641874E-02;
    COFD[5992] = -2.26652914E+01;
    COFD[5993] = 5.37102704E+00;
    COFD[5994] = -4.33426866E-01;
    COFD[5995] = 1.68757765E-02;
    COFD[5996] = -2.29733635E+01;
    COFD[5997] = 5.29835260E+00;
    COFD[5998] = -4.20767237E-01;
    COFD[5999] = 1.62068436E-02;
    COFD[6000] = -2.28592728E+01;
    COFD[6001] = 5.27229756E+00;
    COFD[6002] = -4.16295172E-01;
    COFD[6003] = 1.59725857E-02;
    COFD[6004] = -2.31520360E+01;
    COFD[6005] = 5.33643523E+00;
    COFD[6006] = -4.27338696E-01;
    COFD[6007] = 1.65522282E-02;
    COFD[6008] = -2.31501948E+01;
    COFD[6009] = 5.33643523E+00;
    COFD[6010] = -4.27338696E-01;
    COFD[6011] = 1.65522282E-02;
    COFD[6012] = -2.28187224E+01;
    COFD[6013] = 5.09710954E+00;
    COFD[6014] = -3.87121977E-01;
    COFD[6015] = 1.44693985E-02;
    COFD[6016] = -2.27912671E+01;
    COFD[6017] = 4.98252735E+00;
    COFD[6018] = -3.68805402E-01;
    COFD[6019] = 1.35477963E-02;
    COFD[6020] = -2.29938969E+01;
    COFD[6021] = 5.21722096E+00;
    COFD[6022] = -4.07093564E-01;
    COFD[6023] = 1.54974812E-02;
    COFD[6024] = -2.27668090E+01;
    COFD[6025] = 5.04103112E+00;
    COFD[6026] = -3.78039122E-01;
    COFD[6027] = 1.40083345E-02;
    COFD[6028] = -1.91962804E+01;
    COFD[6029] = 4.54321736E+00;
    COFD[6030] = -3.66462203E-01;
    COFD[6031] = 1.55214344E-02;
    COFD[6032] = -1.92706041E+01;
    COFD[6033] = 5.22641306E+00;
    COFD[6034] = -4.44290196E-01;
    COFD[6035] = 1.84704961E-02;
    COFD[6036] = -1.53871573E+01;
    COFD[6037] = 3.68408734E+00;
    COFD[6038] = -2.61879015E-01;
    COFD[6039] = 1.12666094E-02;
    COFD[6040] = -1.87347447E+01;
    COFD[6041] = 4.55932865E+00;
    COFD[6042] = -3.68355198E-01;
    COFD[6043] = 1.55957356E-02;
    COFD[6044] = -2.02462660E+01;
    COFD[6045] = 4.93517265E+00;
    COFD[6046] = -4.12740828E-01;
    COFD[6047] = 1.73467245E-02;
    COFD[6048] = -1.87584426E+01;
    COFD[6049] = 4.55932865E+00;
    COFD[6050] = -3.68355198E-01;
    COFD[6051] = 1.55957356E-02;
    COFD[6052] = -2.05361209E+01;
    COFD[6053] = 4.64159173E+00;
    COFD[6054] = -3.15993927E-01;
    COFD[6055] = 1.09445813E-02;
    COFD[6056] = -1.99543652E+01;
    COFD[6057] = 4.83774371E+00;
    COFD[6058] = -4.01759538E-01;
    COFD[6059] = 1.69365310E-02;
    COFD[6060] = -2.23824354E+01;
    COFD[6061] = 5.49217181E+00;
    COFD[6062] = -4.63882677E-01;
    COFD[6063] = 1.87489631E-02;
    COFD[6064] = -2.09365593E+01;
    COFD[6065] = 5.21672969E+00;
    COFD[6066] = -4.43175115E-01;
    COFD[6067] = 1.84276681E-02;
    COFD[6068] = -2.09026212E+01;
    COFD[6069] = 5.20282016E+00;
    COFD[6070] = -4.41811086E-01;
    COFD[6071] = 1.83857216E-02;
    COFD[6072] = -2.02561467E+01;
    COFD[6073] = 4.93517265E+00;
    COFD[6074] = -4.12740828E-01;
    COFD[6075] = 1.73467245E-02;
    COFD[6076] = -2.02656235E+01;
    COFD[6077] = 4.93517265E+00;
    COFD[6078] = -4.12740828E-01;
    COFD[6079] = 1.73467245E-02;
    COFD[6080] = -2.15372889E+01;
    COFD[6081] = 4.93380640E+00;
    COFD[6082] = -3.61171699E-01;
    COFD[6083] = 1.31686212E-02;
    COFD[6084] = -2.24143090E+01;
    COFD[6085] = 5.49379864E+00;
    COFD[6086] = -4.63600612E-01;
    COFD[6087] = 1.87179811E-02;
    COFD[6088] = -2.21479968E+01;
    COFD[6089] = 5.48719572E+00;
    COFD[6090] = -4.64022066E-01;
    COFD[6091] = 1.87848905E-02;
    COFD[6092] = -2.24030784E+01;
    COFD[6093] = 5.49379864E+00;
    COFD[6094] = -4.63600612E-01;
    COFD[6095] = 1.87179811E-02;
    COFD[6096] = -2.22655327E+01;
    COFD[6097] = 5.49021632E+00;
    COFD[6098] = -4.60531670E-01;
    COFD[6099] = 1.84855708E-02;
    COFD[6100] = -2.22791718E+01;
    COFD[6101] = 5.49021632E+00;
    COFD[6102] = -4.60531670E-01;
    COFD[6103] = 1.84855708E-02;
    COFD[6104] = -2.16858752E+01;
    COFD[6105] = 4.99630533E+00;
    COFD[6106] = -3.70973757E-01;
    COFD[6107] = 1.36557558E-02;
    COFD[6108] = -2.21212336E+01;
    COFD[6109] = 5.14996158E+00;
    COFD[6110] = -3.95874432E-01;
    COFD[6111] = 1.49190064E-02;
    COFD[6112] = -2.11159694E+01;
    COFD[6113] = 5.26890914E+00;
    COFD[6114] = -4.49080646E-01;
    COFD[6115] = 1.86500850E-02;
    COFD[6116] = -2.21345861E+01;
    COFD[6117] = 5.14996158E+00;
    COFD[6118] = -3.95874432E-01;
    COFD[6119] = 1.49190064E-02;
    COFD[6120] = -2.26109161E+01;
    COFD[6121] = 5.44450212E+00;
    COFD[6122] = -4.47289567E-01;
    COFD[6123] = 1.76392556E-02;
    COFD[6124] = -2.26109161E+01;
    COFD[6125] = 5.44450212E+00;
    COFD[6126] = -4.47289567E-01;
    COFD[6127] = 1.76392556E-02;
    COFD[6128] = -2.25889692E+01;
    COFD[6129] = 5.46241873E+00;
    COFD[6130] = -4.51563227E-01;
    COFD[6131] = 1.78964292E-02;
    COFD[6132] = -2.27070219E+01;
    COFD[6133] = 5.37425628E+00;
    COFD[6134] = -4.33997939E-01;
    COFD[6135] = 1.69062130E-02;
    COFD[6136] = -2.27043153E+01;
    COFD[6137] = 5.37944812E+00;
    COFD[6138] = -4.34919903E-01;
    COFD[6139] = 1.69554449E-02;
    COFD[6140] = -2.27089064E+01;
    COFD[6141] = 5.37944812E+00;
    COFD[6142] = -4.34919903E-01;
    COFD[6143] = 1.69554449E-02;
    COFD[6144] = -2.28439474E+01;
    COFD[6145] = 5.38698965E+00;
    COFD[6146] = -4.36272815E-01;
    COFD[6147] = 1.70280741E-02;
    COFD[6148] = -2.23977002E+01;
    COFD[6149] = 5.15130619E+00;
    COFD[6150] = -3.96101863E-01;
    COFD[6151] = 1.49308178E-02;
    COFD[6152] = -2.23621775E+01;
    COFD[6153] = 5.15211086E+00;
    COFD[6154] = -3.96238056E-01;
    COFD[6155] = 1.49378934E-02;
    COFD[6156] = -2.27558767E+01;
    COFD[6157] = 5.26792797E+00;
    COFD[6158] = -4.15560936E-01;
    COFD[6159] = 1.59345751E-02;
    COFD[6160] = -2.28097485E+01;
    COFD[6161] = 5.29670997E+00;
    COFD[6162] = -4.20481579E-01;
    COFD[6163] = 1.61917846E-02;
    COFD[6164] = -2.24637960E+01;
    COFD[6165] = 5.13708543E+00;
    COFD[6166] = -3.93703669E-01;
    COFD[6167] = 1.48064600E-02;
    COFD[6168] = -2.18000357E+01;
    COFD[6169] = 4.99630533E+00;
    COFD[6170] = -3.70973757E-01;
    COFD[6171] = 1.36557558E-02;
    COFD[6172] = -2.18058201E+01;
    COFD[6173] = 4.99630533E+00;
    COFD[6174] = -3.70973757E-01;
    COFD[6175] = 1.36557558E-02;
    COFD[6176] = -2.21956028E+01;
    COFD[6177] = 5.12996650E+00;
    COFD[6178] = -3.92516016E-01;
    COFD[6179] = 1.47452270E-02;
    COFD[6180] = -2.24415132E+01;
    COFD[6181] = 5.16918377E+00;
    COFD[6182] = -3.99137871E-01;
    COFD[6183] = 1.50888488E-02;
    COFD[6184] = -2.25999700E+01;
    COFD[6185] = 5.45327166E+00;
    COFD[6186] = -4.49436298E-01;
    COFD[6187] = 1.77696558E-02;
    COFD[6188] = -2.26033915E+01;
    COFD[6189] = 5.44450212E+00;
    COFD[6190] = -4.47289567E-01;
    COFD[6191] = 1.76392556E-02;
    COFD[6192] = -2.23275170E+01;
    COFD[6193] = 5.06089774E+00;
    COFD[6194] = -3.81239105E-01;
    COFD[6195] = 1.41702517E-02;
    COFD[6196] = -2.24655743E+01;
    COFD[6197] = 5.09524210E+00;
    COFD[6198] = -3.86817295E-01;
    COFD[6199] = 1.44538658E-02;
    COFD[6200] = -2.22001994E+01;
    COFD[6201] = 5.12996650E+00;
    COFD[6202] = -3.92516016E-01;
    COFD[6203] = 1.47452270E-02;
    COFD[6204] = -2.24208370E+01;
    COFD[6205] = 5.01427991E+00;
    COFD[6206] = -3.73809183E-01;
    COFD[6207] = 1.37971237E-02;
    COFD[6208] = -2.22576038E+01;
    COFD[6209] = 4.96852134E+00;
    COFD[6210] = -3.66607756E-01;
    COFD[6211] = 1.34385708E-02;
    COFD[6212] = -2.26588219E+01;
    COFD[6213] = 5.07646321E+00;
    COFD[6214] = -3.83762054E-01;
    COFD[6215] = 1.42983651E-02;
    COFD[6216] = -2.26569816E+01;
    COFD[6217] = 5.07646321E+00;
    COFD[6218] = -3.83762054E-01;
    COFD[6219] = 1.42983651E-02;
    COFD[6220] = -2.19726300E+01;
    COFD[6221] = 4.68933904E+00;
    COFD[6222] = -3.23268395E-01;
    COFD[6223] = 1.12992644E-02;
    COFD[6224] = -2.18595544E+01;
    COFD[6225] = 4.53718327E+00;
    COFD[6226] = -2.99856359E-01;
    COFD[6227] = 1.01504725E-02;
    COFD[6228] = -2.23303179E+01;
    COFD[6229] = 4.88479479E+00;
    COFD[6230] = -3.53456980E-01;
    COFD[6231] = 1.27844233E-02;
    COFD[6232] = -2.18647351E+01;
    COFD[6233] = 4.61001619E+00;
    COFD[6234] = -3.11086406E-01;
    COFD[6235] = 1.07022976E-02;
    COFD[6236] = -1.99275800E+01;
    COFD[6237] = 4.83095096E+00;
    COFD[6238] = -4.00972558E-01;
    COFD[6239] = 1.69061563E-02;
    COFD[6240] = -1.92100304E+01;
    COFD[6241] = 5.22556324E+00;
    COFD[6242] = -4.44195247E-01;
    COFD[6243] = 1.84669783E-02;
    COFD[6244] = -1.57444919E+01;
    COFD[6245] = 3.82742326E+00;
    COFD[6246] = -2.79536566E-01;
    COFD[6247] = 1.19924367E-02;
    COFD[6248] = -1.86817105E+01;
    COFD[6249] = 4.55856850E+00;
    COFD[6250] = -3.68265791E-01;
    COFD[6251] = 1.55922225E-02;
    COFD[6252] = -2.05570536E+01;
    COFD[6253] = 5.05447966E+00;
    COFD[6254] = -4.26675417E-01;
    COFD[6255] = 1.78900033E-02;
    COFD[6256] = -1.87055009E+01;
    COFD[6257] = 4.55856850E+00;
    COFD[6258] = -3.68265791E-01;
    COFD[6259] = 1.55922225E-02;
    COFD[6260] = -1.83638021E+01;
    COFD[6261] = 3.61550901E+00;
    COFD[6262] = -1.62122786E-01;
    COFD[6263] = 3.57292746E-03;
    COFD[6264] = -2.02463222E+01;
    COFD[6265] = 4.94688197E+00;
    COFD[6266] = -4.14116294E-01;
    COFD[6267] = 1.74006856E-02;
    COFD[6268] = -2.24145023E+01;
    COFD[6269] = 5.49116789E+00;
    COFD[6270] = -4.60872367E-01;
    COFD[6271] = 1.85086095E-02;
    COFD[6272] = -2.08891203E+01;
    COFD[6273] = 5.21570257E+00;
    COFD[6274] = -4.43055160E-01;
    COFD[6275] = 1.84229862E-02;
    COFD[6276] = -2.12234786E+01;
    COFD[6277] = 5.31936752E+00;
    COFD[6278] = -4.54732039E-01;
    COFD[6279] = 1.88607327E-02;
    COFD[6280] = -2.02084294E+01;
    COFD[6281] = 4.93423363E+00;
    COFD[6282] = -4.12630529E-01;
    COFD[6283] = 1.73423980E-02;
    COFD[6284] = -2.02179674E+01;
    COFD[6285] = 4.93423363E+00;
    COFD[6286] = -4.12630529E-01;
    COFD[6287] = 1.73423980E-02;
    COFD[6288] = -2.14953513E+01;
    COFD[6289] = 4.93510998E+00;
    COFD[6290] = -3.61376137E-01;
    COFD[6291] = 1.31787827E-02;
    COFD[6292] = -2.23731099E+01;
    COFD[6293] = 5.49374441E+00;
    COFD[6294] = -4.63617711E-01;
    COFD[6295] = 1.87196409E-02;
    COFD[6296] = -2.21019004E+01;
    COFD[6297] = 5.48702704E+00;
    COFD[6298] = -4.64023035E-01;
    COFD[6299] = 1.87858024E-02;
    COFD[6300] = -2.23618125E+01;
    COFD[6301] = 5.49374441E+00;
    COFD[6302] = -4.63617711E-01;
    COFD[6303] = 1.87196409E-02;
    COFD[6304] = -2.22210765E+01;
    COFD[6305] = 5.49030211E+00;
    COFD[6306] = -4.60568793E-01;
    COFD[6307] = 1.84881663E-02;
    COFD[6308] = -2.22347892E+01;
    COFD[6309] = 5.49030211E+00;
    COFD[6310] = -4.60568793E-01;
    COFD[6311] = 1.84881663E-02;
    COFD[6312] = -2.16441409E+01;
    COFD[6313] = 4.99755422E+00;
    COFD[6314] = -3.71170522E-01;
    COFD[6315] = 1.36655589E-02;
    COFD[6316] = -2.20808599E+01;
    COFD[6317] = 5.15076979E+00;
    COFD[6318] = -3.96011114E-01;
    COFD[6319] = 1.49261042E-02;
    COFD[6320] = -2.10624904E+01;
    COFD[6321] = 5.26805363E+00;
    COFD[6322] = -4.48984207E-01;
    COFD[6323] = 1.86464682E-02;
    COFD[6324] = -2.20943127E+01;
    COFD[6325] = 5.15076979E+00;
    COFD[6326] = -3.96011114E-01;
    COFD[6327] = 1.49261042E-02;
    COFD[6328] = -2.25702108E+01;
    COFD[6329] = 5.44454618E+00;
    COFD[6330] = -4.47321672E-01;
    COFD[6331] = 1.76416355E-02;
    COFD[6332] = -2.25702108E+01;
    COFD[6333] = 5.44454618E+00;
    COFD[6334] = -4.47321672E-01;
    COFD[6335] = 1.76416355E-02;
    COFD[6336] = -2.25481812E+01;
    COFD[6337] = 5.46266513E+00;
    COFD[6338] = -4.51619963E-01;
    COFD[6339] = 1.78998009E-02;
    COFD[6340] = -2.26704052E+01;
    COFD[6341] = 5.37490349E+00;
    COFD[6342] = -4.34112523E-01;
    COFD[6343] = 1.69123235E-02;
    COFD[6344] = -2.26673882E+01;
    COFD[6345] = 5.38007821E+00;
    COFD[6346] = -4.35032051E-01;
    COFD[6347] = 1.69614401E-02;
    COFD[6348] = -2.26720190E+01;
    COFD[6349] = 5.38007821E+00;
    COFD[6350] = -4.35032051E-01;
    COFD[6351] = 1.69614401E-02;
    COFD[6352] = -2.28092883E+01;
    COFD[6353] = 5.38749986E+00;
    COFD[6354] = -4.36368839E-01;
    COFD[6355] = 1.70333535E-02;
    COFD[6356] = -2.23621775E+01;
    COFD[6357] = 5.15211086E+00;
    COFD[6358] = -3.96238056E-01;
    COFD[6359] = 1.49378934E-02;
    COFD[6360] = -2.09249833E+01;
    COFD[6361] = 4.47619102E+00;
    COFD[6362] = -2.94481770E-01;
    COFD[6363] = 1.00314822E-02;
    COFD[6364] = -2.27232029E+01;
    COFD[6365] = 5.26873030E+00;
    COFD[6366] = -4.15695726E-01;
    COFD[6367] = 1.59415526E-02;
    COFD[6368] = -2.26522041E+01;
    COFD[6369] = 5.24607014E+00;
    COFD[6370] = -4.13634549E-01;
    COFD[6371] = 1.58853580E-02;
    COFD[6372] = -2.24294514E+01;
    COFD[6373] = 5.13786054E+00;
    COFD[6374] = -3.93835323E-01;
    COFD[6375] = 1.48133082E-02;
    COFD[6376] = -2.17591073E+01;
    COFD[6377] = 4.99755422E+00;
    COFD[6378] = -3.71170522E-01;
    COFD[6379] = 1.36655589E-02;
    COFD[6380] = -2.17649377E+01;
    COFD[6381] = 4.99755422E+00;
    COFD[6382] = -3.71170522E-01;
    COFD[6383] = 1.36655589E-02;
    COFD[6384] = -2.21570195E+01;
    COFD[6385] = 5.13096622E+00;
    COFD[6386] = -3.92680132E-01;
    COFD[6387] = 1.47536199E-02;
    COFD[6388] = -2.10421374E+01;
    COFD[6389] = 4.51136123E+00;
    COFD[6390] = -3.00118189E-01;
    COFD[6391] = 1.03138650E-02;
    COFD[6392] = -2.25595681E+01;
    COFD[6393] = 5.45353757E+00;
    COFD[6394] = -4.49496004E-01;
    COFD[6395] = 1.77731752E-02;
    COFD[6396] = -2.25626324E+01;
    COFD[6397] = 5.44454618E+00;
    COFD[6398] = -4.47321672E-01;
    COFD[6399] = 1.76416355E-02;
    COFD[6400] = -2.10466968E+01;
    COFD[6401] = 4.48381723E+00;
    COFD[6402] = -2.96509945E-01;
    COFD[6403] = 1.01535197E-02;
    COFD[6404] = -2.10487877E+01;
    COFD[6405] = 4.43820990E+00;
    COFD[6406] = -2.88877613E-01;
    COFD[6407] = 9.76344368E-03;
    COFD[6408] = -2.21616557E+01;
    COFD[6409] = 5.13096622E+00;
    COFD[6410] = -3.92680132E-01;
    COFD[6411] = 1.47536199E-02;
    COFD[6412] = -2.14543740E+01;
    COFD[6413] = 4.59523034E+00;
    COFD[6414] = -3.13534126E-01;
    COFD[6415] = 1.09798752E-02;
    COFD[6416] = -2.22261321E+01;
    COFD[6417] = 4.96982327E+00;
    COFD[6418] = -3.66811857E-01;
    COFD[6419] = 1.34487095E-02;
    COFD[6420] = -2.26303193E+01;
    COFD[6421] = 5.07752972E+00;
    COFD[6422] = -3.83935232E-01;
    COFD[6423] = 1.43071687E-02;
    COFD[6424] = -2.26284586E+01;
    COFD[6425] = 5.07752972E+00;
    COFD[6426] = -3.83935232E-01;
    COFD[6427] = 1.43071687E-02;
    COFD[6428] = -2.14841962E+01;
    COFD[6429] = 4.49540264E+00;
    COFD[6430] = -2.96514053E-01;
    COFD[6431] = 1.00871155E-02;
    COFD[6432] = -2.14584112E+01;
    COFD[6433] = 4.38117808E+00;
    COFD[6434] = -2.78570982E-01;
    COFD[6435] = 9.19430194E-03;
    COFD[6436] = -2.23028472E+01;
    COFD[6437] = 4.88611782E+00;
    COFD[6438] = -3.53664413E-01;
    COFD[6439] = 1.27947308E-02;
    COFD[6440] = -2.13113937E+01;
    COFD[6441] = 4.38760834E+00;
    COFD[6442] = -2.80341923E-01;
    COFD[6443] = 9.30748203E-03;
    COFD[6444] = -2.01873367E+01;
    COFD[6445] = 4.92914806E+00;
    COFD[6446] = -4.12032782E-01;
    COFD[6447] = 1.73189335E-02;
    COFD[6448] = -1.92110979E+01;
    COFD[6449] = 5.16046404E+00;
    COFD[6450] = -4.38008342E-01;
    COFD[6451] = 1.82870041E-02;
    COFD[6452] = -1.51818187E+01;
    COFD[6453] = 3.55111766E+00;
    COFD[6454] = -2.44966429E-01;
    COFD[6455] = 1.05482272E-02;
    COFD[6456] = -1.86097280E+01;
    COFD[6457] = 4.45757911E+00;
    COFD[6458] = -3.56368733E-01;
    COFD[6459] = 1.51238590E-02;
    COFD[6460] = -2.01150260E+01;
    COFD[6461] = 4.83483781E+00;
    COFD[6462] = -4.01422871E-01;
    COFD[6463] = 1.69235370E-02;
    COFD[6464] = -1.86343902E+01;
    COFD[6465] = 4.45757911E+00;
    COFD[6466] = -3.56368733E-01;
    COFD[6467] = 1.51238590E-02;
    COFD[6468] = -2.11196901E+01;
    COFD[6469] = 4.84742105E+00;
    COFD[6470] = -3.47604370E-01;
    COFD[6471] = 1.24938668E-02;
    COFD[6472] = -1.97966756E+01;
    COFD[6473] = 4.72725308E+00;
    COFD[6474] = -3.88865048E-01;
    COFD[6475] = 1.64346151E-02;
    COFD[6476] = -2.24163735E+01;
    COFD[6477] = 5.46720411E+00;
    COFD[6478] = -4.63664606E-01;
    COFD[6479] = 1.88543495E-02;
    COFD[6480] = -2.08689918E+01;
    COFD[6481] = 5.15366259E+00;
    COFD[6482] = -4.37268432E-01;
    COFD[6483] = 1.82604618E-02;
    COFD[6484] = -2.08290042E+01;
    COFD[6485] = 5.13624657E+00;
    COFD[6486] = -4.35395745E-01;
    COFD[6487] = 1.81943797E-02;
    COFD[6488] = -2.01255741E+01;
    COFD[6489] = 4.83483781E+00;
    COFD[6490] = -4.01422871E-01;
    COFD[6491] = 1.69235370E-02;
    COFD[6492] = -2.01357044E+01;
    COFD[6493] = 4.83483781E+00;
    COFD[6494] = -4.01422871E-01;
    COFD[6495] = 1.69235370E-02;
    COFD[6496] = -2.20062771E+01;
    COFD[6497] = 5.09608960E+00;
    COFD[6498] = -3.86955479E-01;
    COFD[6499] = 1.44609083E-02;
    COFD[6500] = -2.24306557E+01;
    COFD[6501] = 5.46873098E+00;
    COFD[6502] = -4.63395699E-01;
    COFD[6503] = 1.88249267E-02;
    COFD[6504] = -2.21898994E+01;
    COFD[6505] = 5.46750298E+00;
    COFD[6506] = -4.64484282E-01;
    COFD[6507] = 1.89196260E-02;
    COFD[6508] = -2.24187146E+01;
    COFD[6509] = 5.46873098E+00;
    COFD[6510] = -4.63395699E-01;
    COFD[6511] = 1.88249267E-02;
    COFD[6512] = -2.23433381E+01;
    COFD[6513] = 5.49008914E+00;
    COFD[6514] = -4.63927959E-01;
    COFD[6515] = 1.87629280E-02;
    COFD[6516] = -2.23577541E+01;
    COFD[6517] = 5.49008914E+00;
    COFD[6518] = -4.63927959E-01;
    COFD[6519] = 1.87629280E-02;
    COFD[6520] = -2.21143598E+01;
    COFD[6521] = 5.14170132E+00;
    COFD[6522] = -3.94484367E-01;
    COFD[6523] = 1.48469922E-02;
    COFD[6524] = -2.24838323E+01;
    COFD[6525] = 5.26658797E+00;
    COFD[6526] = -4.15335910E-01;
    COFD[6527] = 1.59229288E-02;
    COFD[6528] = -2.10353566E+01;
    COFD[6529] = 5.18365882E+00;
    COFD[6530] = -4.40195806E-01;
    COFD[6531] = 1.83497181E-02;
    COFD[6532] = -2.24982668E+01;
    COFD[6533] = 5.26658797E+00;
    COFD[6534] = -4.15335910E-01;
    COFD[6535] = 1.59229288E-02;
    COFD[6536] = -2.27673763E+01;
    COFD[6537] = 5.47456971E+00;
    COFD[6538] = -4.54548739E-01;
    COFD[6539] = 1.80779412E-02;
    COFD[6540] = -2.27673763E+01;
    COFD[6541] = 5.47456971E+00;
    COFD[6542] = -4.54548739E-01;
    COFD[6543] = 1.80779412E-02;
    COFD[6544] = -2.27255842E+01;
    COFD[6545] = 5.48298766E+00;
    COFD[6546] = -4.57543838E-01;
    COFD[6547] = 1.82786713E-02;
    COFD[6548] = -2.29510526E+01;
    COFD[6549] = 5.44399037E+00;
    COFD[6550] = -4.47161064E-01;
    COFD[6551] = 1.76313359E-02;
    COFD[6552] = -2.29382506E+01;
    COFD[6553] = 5.44442303E+00;
    COFD[6554] = -4.47427940E-01;
    COFD[6555] = 1.76507657E-02;
    COFD[6556] = -2.29432741E+01;
    COFD[6557] = 5.44442303E+00;
    COFD[6558] = -4.47427940E-01;
    COFD[6559] = 1.76507657E-02;
    COFD[6560] = -2.30611542E+01;
    COFD[6561] = 5.44729879E+00;
    COFD[6562] = -4.48113754E-01;
    COFD[6563] = 1.76920689E-02;
    COFD[6564] = -2.27558767E+01;
    COFD[6565] = 5.26792797E+00;
    COFD[6566] = -4.15560936E-01;
    COFD[6567] = 1.59345751E-02;
    COFD[6568] = -2.27232029E+01;
    COFD[6569] = 5.26873030E+00;
    COFD[6570] = -4.15695726E-01;
    COFD[6571] = 1.59415526E-02;
    COFD[6572] = -2.30607499E+01;
    COFD[6573] = 5.36239026E+00;
    COFD[6574] = -4.31899984E-01;
    COFD[6575] = 1.67944393E-02;
    COFD[6576] = -2.31045202E+01;
    COFD[6577] = 5.38724026E+00;
    COFD[6578] = -4.36319971E-01;
    COFD[6579] = 1.70306665E-02;
    COFD[6580] = -2.28257638E+01;
    COFD[6581] = 5.25412374E+00;
    COFD[6582] = -4.13233123E-01;
    COFD[6583] = 1.58137928E-02;
    COFD[6584] = -2.22371738E+01;
    COFD[6585] = 5.14170132E+00;
    COFD[6586] = -3.94484367E-01;
    COFD[6587] = 1.48469922E-02;
    COFD[6588] = -2.22434560E+01;
    COFD[6589] = 5.14170132E+00;
    COFD[6590] = -3.94484367E-01;
    COFD[6591] = 1.48469922E-02;
    COFD[6592] = -2.25643461E+01;
    COFD[6593] = 5.24758220E+00;
    COFD[6594] = -4.12129116E-01;
    COFD[6595] = 1.57564573E-02;
    COFD[6596] = -2.27971335E+01;
    COFD[6597] = 5.28513609E+00;
    COFD[6598] = -4.18480689E-01;
    COFD[6599] = 1.60865907E-02;
    COFD[6600] = -2.27488997E+01;
    COFD[6601] = 5.47944787E+00;
    COFD[6602] = -4.56152622E-01;
    COFD[6603] = 1.81835823E-02;
    COFD[6604] = -2.27592733E+01;
    COFD[6605] = 5.47456971E+00;
    COFD[6606] = -4.54548739E-01;
    COFD[6607] = 1.80779412E-02;
    COFD[6608] = -2.27023831E+01;
    COFD[6609] = 5.18383798E+00;
    COFD[6610] = -4.01666785E-01;
    COFD[6611] = 1.52216204E-02;
    COFD[6612] = -2.28320950E+01;
    COFD[6613] = 5.21637493E+00;
    COFD[6614] = -4.06956958E-01;
    COFD[6615] = 1.54905633E-02;
    COFD[6616] = -2.25693752E+01;
    COFD[6617] = 5.24758220E+00;
    COFD[6618] = -4.12129116E-01;
    COFD[6619] = 1.57564573E-02;
    COFD[6620] = -2.28314712E+01;
    COFD[6621] = 5.15368301E+00;
    COFD[6622] = -3.96504343E-01;
    COFD[6623] = 1.49517329E-02;
    COFD[6624] = -2.27064705E+01;
    COFD[6625] = 5.12211189E+00;
    COFD[6626] = -3.91225928E-01;
    COFD[6627] = 1.46792378E-02;
    COFD[6628] = -2.30240568E+01;
    COFD[6629] = 5.19737465E+00;
    COFD[6630] = -4.03890980E-01;
    COFD[6631] = 1.53353180E-02;
    COFD[6632] = -2.30219870E+01;
    COFD[6633] = 5.19737465E+00;
    COFD[6634] = -4.03890980E-01;
    COFD[6635] = 1.53353180E-02;
    COFD[6636] = -2.25275251E+01;
    COFD[6637] = 4.88611782E+00;
    COFD[6638] = -3.53664413E-01;
    COFD[6639] = 1.27947308E-02;
    COFD[6640] = -2.24440124E+01;
    COFD[6641] = 4.74813432E+00;
    COFD[6642] = -3.32238042E-01;
    COFD[6643] = 1.17370757E-02;
    COFD[6644] = -2.28166735E+01;
    COFD[6645] = 5.05461151E+00;
    COFD[6646] = -3.80222586E-01;
    COFD[6647] = 1.41187082E-02;
    COFD[6648] = -2.24468874E+01;
    COFD[6649] = 4.81788504E+00;
    COFD[6650] = -3.42990993E-01;
    COFD[6651] = 1.22653559E-02;
    COFD[6652] = -1.97684750E+01;
    COFD[6653] = 4.71961665E+00;
    COFD[6654] = -3.87953877E-01;
    COFD[6655] = 1.63982506E-02;
    COFD[6656] = -1.91375595E+01;
    COFD[6657] = 5.13436741E+00;
    COFD[6658] = -4.35193261E-01;
    COFD[6659] = 1.81872196E-02;
    COFD[6660] = -1.50908350E+01;
    COFD[6661] = 3.51404074E+00;
    COFD[6662] = -2.40221542E-01;
    COFD[6663] = 1.03454651E-02;
    COFD[6664] = -1.85284648E+01;
    COFD[6665] = 4.42603366E+00;
    COFD[6666] = -3.52586895E-01;
    COFD[6667] = 1.49720921E-02;
    COFD[6668] = -2.00372657E+01;
    COFD[6669] = 4.80531581E+00;
    COFD[6670] = -3.98010400E-01;
    COFD[6671] = 1.67921710E-02;
    COFD[6672] = -1.85531954E+01;
    COFD[6673] = 4.42603366E+00;
    COFD[6674] = -3.52586895E-01;
    COFD[6675] = 1.49720921E-02;
    COFD[6676] = -2.10582825E+01;
    COFD[6677] = 4.83010534E+00;
    COFD[6678] = -3.46423002E-01;
    COFD[6679] = 1.24833814E-02;
    COFD[6680] = -1.97045701E+01;
    COFD[6681] = 4.69156097E+00;
    COFD[6682] = -3.84562356E-01;
    COFD[6683] = 1.62609922E-02;
    COFD[6684] = -2.24050898E+01;
    COFD[6685] = 5.46747549E+00;
    COFD[6686] = -4.64570865E-01;
    COFD[6687] = 1.89267998E-02;
    COFD[6688] = -2.07942327E+01;
    COFD[6689] = 5.12723540E+00;
    COFD[6690] = -4.34411432E-01;
    COFD[6691] = 1.81589442E-02;
    COFD[6692] = -2.07580489E+01;
    COFD[6693] = 5.11080865E+00;
    COFD[6694] = -4.32626175E-01;
    COFD[6695] = 1.80950687E-02;
    COFD[6696] = -2.00433341E+01;
    COFD[6697] = 4.80382495E+00;
    COFD[6698] = -3.97838890E-01;
    COFD[6699] = 1.67856063E-02;
    COFD[6700] = -2.00535121E+01;
    COFD[6701] = 4.80382495E+00;
    COFD[6702] = -3.97838890E-01;
    COFD[6703] = 1.67856063E-02;
    COFD[6704] = -2.20820648E+01;
    COFD[6705] = 5.13480519E+00;
    COFD[6706] = -3.93316744E-01;
    COFD[6707] = 1.47863438E-02;
    COFD[6708] = -2.24132666E+01;
    COFD[6709] = 5.46734823E+00;
    COFD[6710] = -4.64115641E-01;
    COFD[6711] = 1.88903818E-02;
    COFD[6712] = -2.21759649E+01;
    COFD[6713] = 5.46726363E+00;
    COFD[6714] = -4.65349411E-01;
    COFD[6715] = 1.89911285E-02;
    COFD[6716] = -2.24012739E+01;
    COFD[6717] = 5.46734823E+00;
    COFD[6718] = -4.64115641E-01;
    COFD[6719] = 1.88903818E-02;
    COFD[6720] = -2.23119605E+01;
    COFD[6721] = 5.48315923E+00;
    COFD[6722] = -4.63919102E-01;
    COFD[6723] = 1.87963630E-02;
    COFD[6724] = -2.23264327E+01;
    COFD[6725] = 5.48315923E+00;
    COFD[6726] = -4.63919102E-01;
    COFD[6727] = 1.87963630E-02;
    COFD[6728] = -2.21735855E+01;
    COFD[6729] = 5.17327578E+00;
    COFD[6730] = -3.99835717E-01;
    COFD[6731] = 1.51252623E-02;
    COFD[6732] = -2.25382696E+01;
    COFD[6733] = 5.29559685E+00;
    COFD[6734] = -4.20288187E-01;
    COFD[6735] = 1.61815942E-02;
    COFD[6736] = -2.09861375E+01;
    COFD[6737] = 5.16708806E+00;
    COFD[6738] = -4.38725884E-01;
    COFD[6739] = 1.83125969E-02;
    COFD[6740] = -2.25527842E+01;
    COFD[6741] = 5.29559685E+00;
    COFD[6742] = -4.20288187E-01;
    COFD[6743] = 1.61815942E-02;
    COFD[6744] = -2.27646473E+01;
    COFD[6745] = 5.47937827E+00;
    COFD[6746] = -4.56129749E-01;
    COFD[6747] = 1.81820758E-02;
    COFD[6748] = -2.27646473E+01;
    COFD[6749] = 5.47937827E+00;
    COFD[6750] = -4.56129749E-01;
    COFD[6751] = 1.81820758E-02;
    COFD[6752] = -2.27182958E+01;
    COFD[6753] = 5.48606693E+00;
    COFD[6754] = -4.58921993E-01;
    COFD[6755] = 1.83751829E-02;
    COFD[6756] = -2.29571084E+01;
    COFD[6757] = 5.45255106E+00;
    COFD[6758] = -4.49274880E-01;
    COFD[6759] = 1.77601485E-02;
    COFD[6760] = -2.29489166E+01;
    COFD[6761] = 5.45471156E+00;
    COFD[6762] = -4.49760531E-01;
    COFD[6763] = 1.77887865E-02;
    COFD[6764] = -2.29539724E+01;
    COFD[6765] = 5.45471156E+00;
    COFD[6766] = -4.49760531E-01;
    COFD[6767] = 1.77887865E-02;
    COFD[6768] = -2.30719702E+01;
    COFD[6769] = 5.45769390E+00;
    COFD[6770] = -4.50452338E-01;
    COFD[6771] = 1.78300009E-02;
    COFD[6772] = -2.28097485E+01;
    COFD[6773] = 5.29670997E+00;
    COFD[6774] = -4.20481579E-01;
    COFD[6775] = 1.61917846E-02;
    COFD[6776] = -2.26522041E+01;
    COFD[6777] = 5.24607014E+00;
    COFD[6778] = -4.13634549E-01;
    COFD[6779] = 1.58853580E-02;
    COFD[6780] = -2.31045202E+01;
    COFD[6781] = 5.38724026E+00;
    COFD[6782] = -4.36319971E-01;
    COFD[6783] = 1.70306665E-02;
    COFD[6784] = -2.31235339E+01;
    COFD[6785] = 5.40224339E+00;
    COFD[6786] = -4.39442969E-01;
    COFD[6787] = 1.72102199E-02;
    COFD[6788] = -2.28842083E+01;
    COFD[6789] = 5.28475021E+00;
    COFD[6790] = -4.18414253E-01;
    COFD[6791] = 1.60831047E-02;
    COFD[6792] = -2.22970365E+01;
    COFD[6793] = 5.17327578E+00;
    COFD[6794] = -3.99835717E-01;
    COFD[6795] = 1.51252623E-02;
    COFD[6796] = -2.23033557E+01;
    COFD[6797] = 5.17327578E+00;
    COFD[6798] = -3.99835717E-01;
    COFD[6799] = 1.51252623E-02;
    COFD[6800] = -2.26246057E+01;
    COFD[6801] = 5.27897272E+00;
    COFD[6802] = -4.17421132E-01;
    COFD[6803] = 1.60310357E-02;
    COFD[6804] = -2.27232104E+01;
    COFD[6805] = 5.26126980E+00;
    COFD[6806] = -4.16242057E-01;
    COFD[6807] = 1.60218084E-02;
    COFD[6808] = -2.27431590E+01;
    COFD[6809] = 5.48304850E+00;
    COFD[6810] = -4.57579455E-01;
    COFD[6811] = 1.82812528E-02;
    COFD[6812] = -2.27565017E+01;
    COFD[6813] = 5.47937827E+00;
    COFD[6814] = -4.56129749E-01;
    COFD[6815] = 1.81820758E-02;
    COFD[6816] = -2.26781027E+01;
    COFD[6817] = 5.18129009E+00;
    COFD[6818] = -4.02449417E-01;
    COFD[6819] = 1.52976330E-02;
    COFD[6820] = -2.27790113E+01;
    COFD[6821] = 5.20168946E+00;
    COFD[6822] = -4.06049775E-01;
    COFD[6823] = 1.54891495E-02;
    COFD[6824] = -2.26296673E+01;
    COFD[6825] = 5.27897272E+00;
    COFD[6826] = -4.17421132E-01;
    COFD[6827] = 1.60310357E-02;
    COFD[6828] = -2.28089000E+01;
    COFD[6829] = 5.15221867E+00;
    COFD[6830] = -3.97461937E-01;
    COFD[6831] = 1.50362563E-02;
    COFD[6832] = -2.27719473E+01;
    COFD[6833] = 5.15606268E+00;
    COFD[6834] = -3.96907901E-01;
    COFD[6835] = 1.49727202E-02;
    COFD[6836] = -2.31022549E+01;
    COFD[6837] = 5.23643800E+00;
    COFD[6838] = -4.10255644E-01;
    COFD[6839] = 1.56593397E-02;
    COFD[6840] = -2.31001674E+01;
    COFD[6841] = 5.23643800E+00;
    COFD[6842] = -4.10255644E-01;
    COFD[6843] = 1.56593397E-02;
    COFD[6844] = -2.25844328E+01;
    COFD[6845] = 4.91727638E+00;
    COFD[6846] = -3.59056497E-01;
    COFD[6847] = 1.30783726E-02;
    COFD[6848] = -2.25319271E+01;
    COFD[6849] = 4.79201228E+00;
    COFD[6850] = -3.39345465E-01;
    COFD[6851] = 1.20970888E-02;
    COFD[6852] = -2.29010051E+01;
    COFD[6853] = 5.09658829E+00;
    COFD[6854] = -3.87036802E-01;
    COFD[6855] = 1.44650532E-02;
    COFD[6856] = -2.25039575E+01;
    COFD[6857] = 4.84915044E+00;
    COFD[6858] = -3.48402726E-01;
    COFD[6859] = 1.25499430E-02;
    COFD[6860] = -1.96740282E+01;
    COFD[6861] = 4.68296253E+00;
    COFD[6862] = -3.83520383E-01;
    COFD[6863] = 1.62187256E-02;
    COFD[6864] = -1.93697973E+01;
    COFD[6865] = 5.24163725E+00;
    COFD[6866] = -4.46011560E-01;
    COFD[6867] = 1.85352610E-02;
    COFD[6868] = -1.54795771E+01;
    COFD[6869] = 3.69992299E+00;
    COFD[6870] = -2.63867460E-01;
    COFD[6871] = 1.13499698E-02;
    COFD[6872] = -1.88465938E+01;
    COFD[6873] = 4.57429837E+00;
    COFD[6874] = -3.70149872E-01;
    COFD[6875] = 1.56677324E-02;
    COFD[6876] = -2.03732076E+01;
    COFD[6877] = 4.95138566E+00;
    COFD[6878] = -4.14645044E-01;
    COFD[6879] = 1.74214162E-02;
    COFD[6880] = -1.88713913E+01;
    COFD[6881] = 4.57429837E+00;
    COFD[6882] = -3.70149872E-01;
    COFD[6883] = 1.56677324E-02;
    COFD[6884] = -2.05625147E+01;
    COFD[6885] = 4.61780386E+00;
    COFD[6886] = -3.12290044E-01;
    COFD[6887] = 1.07615587E-02;
    COFD[6888] = -2.00697375E+01;
    COFD[6889] = 4.85107798E+00;
    COFD[6890] = -4.03289235E-01;
    COFD[6891] = 1.69949072E-02;
    COFD[6892] = -2.24829635E+01;
    COFD[6893] = 5.49361964E+00;
    COFD[6894] = -4.63657707E-01;
    COFD[6895] = 1.87235131E-02;
    COFD[6896] = -2.10468267E+01;
    COFD[6897] = 5.23281177E+00;
    COFD[6898] = -4.45013021E-01;
    COFD[6899] = 1.84976592E-02;
    COFD[6900] = -2.09922482E+01;
    COFD[6901] = 5.20990101E+00;
    COFD[6902] = -4.42427060E-01;
    COFD[6903] = 1.84006140E-02;
    COFD[6904] = -2.03838522E+01;
    COFD[6905] = 4.95138566E+00;
    COFD[6906] = -4.14645044E-01;
    COFD[6907] = 1.74214162E-02;
    COFD[6908] = -2.03940771E+01;
    COFD[6909] = 4.95138566E+00;
    COFD[6910] = -4.14645044E-01;
    COFD[6911] = 1.74214162E-02;
    COFD[6912] = -2.15725574E+01;
    COFD[6913] = 4.91148766E+00;
    COFD[6914] = -3.57668763E-01;
    COFD[6915] = 1.29944614E-02;
    COFD[6916] = -2.24986477E+01;
    COFD[6917] = 5.49353251E+00;
    COFD[6918] = -4.63127751E-01;
    COFD[6919] = 1.86807325E-02;
    COFD[6920] = -2.22409758E+01;
    COFD[6921] = 5.48966728E+00;
    COFD[6922] = -4.63937028E-01;
    COFD[6923] = 1.87657460E-02;
    COFD[6924] = -2.24866045E+01;
    COFD[6925] = 5.49353251E+00;
    COFD[6926] = -4.63127751E-01;
    COFD[6927] = 1.86807325E-02;
    COFD[6928] = -2.23453522E+01;
    COFD[6929] = 5.48861596E+00;
    COFD[6930] = -4.59873623E-01;
    COFD[6931] = 1.84399621E-02;
    COFD[6932] = -2.23598793E+01;
    COFD[6933] = 5.48861596E+00;
    COFD[6934] = -4.59873623E-01;
    COFD[6935] = 1.84399621E-02;
    COFD[6936] = -2.17232952E+01;
    COFD[6937] = 4.97429231E+00;
    COFD[6938] = -3.67512538E-01;
    COFD[6939] = 1.34835182E-02;
    COFD[6940] = -2.21828607E+01;
    COFD[6941] = 5.13579028E+00;
    COFD[6942] = -3.93483830E-01;
    COFD[6943] = 1.47950287E-02;
    COFD[6944] = -2.12518544E+01;
    COFD[6945] = 5.28351879E+00;
    COFD[6946] = -4.50726414E-01;
    COFD[6947] = 1.87118076E-02;
    COFD[6948] = -2.21974539E+01;
    COFD[6949] = 5.13579028E+00;
    COFD[6950] = -3.93483830E-01;
    COFD[6951] = 1.47950287E-02;
    COFD[6952] = -2.26859156E+01;
    COFD[6953] = 5.43682302E+00;
    COFD[6954] = -4.45758386E-01;
    COFD[6955] = 1.75526921E-02;
    COFD[6956] = -2.26859156E+01;
    COFD[6957] = 5.43682302E+00;
    COFD[6958] = -4.45758386E-01;
    COFD[6959] = 1.75526921E-02;
    COFD[6960] = -2.26730180E+01;
    COFD[6961] = 5.45838393E+00;
    COFD[6962] = -4.50616109E-01;
    COFD[6963] = 1.78398266E-02;
    COFD[6964] = -2.27791634E+01;
    COFD[6965] = 5.36302327E+00;
    COFD[6966] = -4.32012062E-01;
    COFD[6967] = 1.68004123E-02;
    COFD[6968] = -2.27773775E+01;
    COFD[6969] = 5.36826093E+00;
    COFD[6970] = -4.32938535E-01;
    COFD[6971] = 1.68497727E-02;
    COFD[6972] = -2.27824652E+01;
    COFD[6973] = 5.36826093E+00;
    COFD[6974] = -4.32938535E-01;
    COFD[6975] = 1.68497727E-02;
    COFD[6976] = -2.29169419E+01;
    COFD[6977] = 5.37626442E+00;
    COFD[6978] = -4.34353973E-01;
    COFD[6979] = 1.69252109E-02;
    COFD[6980] = -2.24637960E+01;
    COFD[6981] = 5.13708543E+00;
    COFD[6982] = -3.93703669E-01;
    COFD[6983] = 1.48064600E-02;
    COFD[6984] = -2.24294514E+01;
    COFD[6985] = 5.13786054E+00;
    COFD[6986] = -3.93835323E-01;
    COFD[6987] = 1.48133082E-02;
    COFD[6988] = -2.28257638E+01;
    COFD[6989] = 5.25412374E+00;
    COFD[6990] = -4.13233123E-01;
    COFD[6991] = 1.58137928E-02;
    COFD[6992] = -2.28842083E+01;
    COFD[6993] = 5.28475021E+00;
    COFD[6994] = -4.18414253E-01;
    COFD[6995] = 1.60831047E-02;
    COFD[6996] = -2.25305663E+01;
    COFD[6997] = 5.12081054E+00;
    COFD[6998] = -3.91011719E-01;
    COFD[6999] = 1.46682690E-02;
    COFD[7000] = -2.18473715E+01;
    COFD[7001] = 4.97429231E+00;
    COFD[7002] = -3.67512538E-01;
    COFD[7003] = 1.34835182E-02;
    COFD[7004] = -2.18537271E+01;
    COFD[7005] = 4.97429231E+00;
    COFD[7006] = -3.67512538E-01;
    COFD[7007] = 1.34835182E-02;
    COFD[7008] = -2.22571178E+01;
    COFD[7009] = 5.11286874E+00;
    COFD[7010] = -3.89706255E-01;
    COFD[7011] = 1.46014745E-02;
    COFD[7012] = -2.25079550E+01;
    COFD[7013] = 5.15543496E+00;
    COFD[7014] = -3.96801390E-01;
    COFD[7015] = 1.49671794E-02;
    COFD[7016] = -2.26814510E+01;
    COFD[7017] = 5.44849592E+00;
    COFD[7018] = -4.48376165E-01;
    COFD[7019] = 1.77074093E-02;
    COFD[7020] = -2.26777281E+01;
    COFD[7021] = 5.43682302E+00;
    COFD[7022] = -4.45758386E-01;
    COFD[7023] = 1.75526921E-02;
    COFD[7024] = -2.23869377E+01;
    COFD[7025] = 5.04130371E+00;
    COFD[7026] = -3.78082012E-01;
    COFD[7027] = 1.40104739E-02;
    COFD[7028] = -2.25262686E+01;
    COFD[7029] = 5.07694371E+00;
    COFD[7030] = -3.83840072E-01;
    COFD[7031] = 1.43023310E-02;
    COFD[7032] = -2.22622113E+01;
    COFD[7033] = 5.11286874E+00;
    COFD[7034] = -3.89706255E-01;
    COFD[7035] = 1.46014745E-02;
    COFD[7036] = -2.24776133E+01;
    COFD[7037] = 4.99278627E+00;
    COFD[7038] = -3.70419518E-01;
    COFD[7039] = 1.36281487E-02;
    COFD[7040] = -2.23118869E+01;
    COFD[7041] = 4.94556893E+00;
    COFD[7042] = -3.63013720E-01;
    COFD[7043] = 1.32601068E-02;
    COFD[7044] = -2.27251033E+01;
    COFD[7045] = 5.05775235E+00;
    COFD[7046] = -3.80730307E-01;
    COFD[7047] = 1.41444474E-02;
    COFD[7048] = -2.27229983E+01;
    COFD[7049] = 5.05775235E+00;
    COFD[7050] = -3.80730307E-01;
    COFD[7051] = 1.41444474E-02;
    COFD[7052] = -2.20211813E+01;
    COFD[7053] = 4.66069984E+00;
    COFD[7054] = -3.18907247E-01;
    COFD[7055] = 1.10866672E-02;
    COFD[7056] = -2.19159785E+01;
    COFD[7057] = 4.51065350E+00;
    COFD[7058] = -2.95784060E-01;
    COFD[7059] = 9.95114964E-03;
    COFD[7060] = -2.23893523E+01;
    COFD[7061] = 4.86181814E+00;
    COFD[7062] = -3.49857477E-01;
    COFD[7063] = 1.26056712E-02;
    COFD[7064] = -2.19212431E+01;
    COFD[7065] = 4.58361215E+00;
    COFD[7066] = -3.07008194E-01;
    COFD[7067] = 1.05016776E-02;
    COFD[7068] = -2.00440183E+01;
    COFD[7069] = 4.84467881E+00;
    COFD[7070] = -4.02557342E-01;
    COFD[7071] = 1.69670693E-02;
    COFD[7072] = -1.92369345E+01;
    COFD[7073] = 5.34516253E+00;
    COFD[7074] = -4.57548788E-01;
    COFD[7075] = 1.89625374E-02;
    COFD[7076] = -1.53691667E+01;
    COFD[7077] = 3.80403213E+00;
    COFD[7078] = -2.76736813E-01;
    COFD[7079] = 1.18808568E-02;
    COFD[7080] = -1.87900753E+01;
    COFD[7081] = 4.71430273E+00;
    COFD[7082] = -3.87319126E-01;
    COFD[7083] = 1.63728857E-02;
    COFD[7084] = -2.02632221E+01;
    COFD[7085] = 5.06295779E+00;
    COFD[7086] = -4.27579258E-01;
    COFD[7087] = 1.79215546E-02;
    COFD[7088] = -1.88126908E+01;
    COFD[7089] = 4.71430273E+00;
    COFD[7090] = -3.87319126E-01;
    COFD[7091] = 1.63728857E-02;
    COFD[7092] = -1.96948070E+01;
    COFD[7093] = 4.41245751E+00;
    COFD[7094] = -2.80855500E-01;
    COFD[7095] = 9.22505777E-03;
    COFD[7096] = -1.99421880E+01;
    COFD[7097] = 4.94937905E+00;
    COFD[7098] = -4.14409534E-01;
    COFD[7099] = 1.74121859E-02;
    COFD[7100] = -2.21101148E+01;
    COFD[7101] = 5.48782481E+00;
    COFD[7102] = -4.59567900E-01;
    COFD[7103] = 1.84190115E-02;
    COFD[7104] = -2.09768682E+01;
    COFD[7105] = 5.33874924E+00;
    COFD[7106] = -4.56879272E-01;
    COFD[7107] = 1.89397603E-02;
    COFD[7108] = -2.09241822E+01;
    COFD[7109] = 5.31839685E+00;
    COFD[7110] = -4.54623948E-01;
    COFD[7111] = 1.88567295E-02;
    COFD[7112] = -2.02723926E+01;
    COFD[7113] = 5.06295779E+00;
    COFD[7114] = -4.27579258E-01;
    COFD[7115] = 1.79215546E-02;
    COFD[7116] = -2.02811760E+01;
    COFD[7117] = 5.06295779E+00;
    COFD[7118] = -4.27579258E-01;
    COFD[7119] = 1.79215546E-02;
    COFD[7120] = -2.07973123E+01;
    COFD[7121] = 4.73262810E+00;
    COFD[7122] = -3.29863658E-01;
    COFD[7123] = 1.16208974E-02;
    COFD[7124] = -2.21631238E+01;
    COFD[7125] = 5.48560730E+00;
    COFD[7126] = -4.58746417E-01;
    COFD[7127] = 1.83631968E-02;
    COFD[7128] = -2.18947817E+01;
    COFD[7129] = 5.49121786E+00;
    COFD[7130] = -4.60888505E-01;
    COFD[7131] = 1.85096783E-02;
    COFD[7132] = -2.21526557E+01;
    COFD[7133] = 5.48560730E+00;
    COFD[7134] = -4.58746417E-01;
    COFD[7135] = 1.83631968E-02;
    COFD[7136] = -2.19781845E+01;
    COFD[7137] = 5.47597765E+00;
    COFD[7138] = -4.55020526E-01;
    COFD[7139] = 1.81091205E-02;
    COFD[7140] = -2.19909795E+01;
    COFD[7141] = 5.47597765E+00;
    COFD[7142] = -4.55020526E-01;
    COFD[7143] = 1.81091205E-02;
    COFD[7144] = -2.09711478E+01;
    COFD[7145] = 4.80548676E+00;
    COFD[7146] = -3.41059239E-01;
    COFD[7147] = 1.21698161E-02;
    COFD[7148] = -2.14996220E+01;
    COFD[7149] = 4.99422177E+00;
    COFD[7150] = -3.70645569E-01;
    COFD[7151] = 1.36394075E-02;
    COFD[7152] = -2.10414466E+01;
    COFD[7153] = 5.37676571E+00;
    COFD[7154] = -4.60857509E-01;
    COFD[7155] = 1.90756062E-02;
    COFD[7156] = -2.15118523E+01;
    COFD[7157] = 4.99422177E+00;
    COFD[7158] = -3.70645569E-01;
    COFD[7159] = 1.36394075E-02;
    COFD[7160] = -2.21918593E+01;
    COFD[7161] = 5.36809191E+00;
    COFD[7162] = -4.32908721E-01;
    COFD[7163] = 1.68481857E-02;
    COFD[7164] = -2.21918593E+01;
    COFD[7165] = 5.36809191E+00;
    COFD[7166] = -4.32908721E-01;
    COFD[7167] = 1.68481857E-02;
    COFD[7168] = -2.22224154E+01;
    COFD[7169] = 5.41240338E+00;
    COFD[7170] = -4.41033923E-01;
    COFD[7171] = 1.72894624E-02;
    COFD[7172] = -2.22319743E+01;
    COFD[7173] = 5.27340167E+00;
    COFD[7174] = -4.16479301E-01;
    COFD[7175] = 1.59820820E-02;
    COFD[7176] = -2.22294897E+01;
    COFD[7177] = 5.27973457E+00;
    COFD[7178] = -4.17551922E-01;
    COFD[7179] = 1.60378887E-02;
    COFD[7180] = -2.22336433E+01;
    COFD[7181] = 5.27973457E+00;
    COFD[7182] = -4.17551922E-01;
    COFD[7183] = 1.60378887E-02;
    COFD[7184] = -2.23886928E+01;
    COFD[7185] = 5.28846317E+00;
    COFD[7186] = -4.19054070E-01;
    COFD[7187] = 1.61166909E-02;
    COFD[7188] = -2.18000357E+01;
    COFD[7189] = 4.99630533E+00;
    COFD[7190] = -3.70973757E-01;
    COFD[7191] = 1.36557558E-02;
    COFD[7192] = -2.17591073E+01;
    COFD[7193] = 4.99755422E+00;
    COFD[7194] = -3.71170522E-01;
    COFD[7195] = 1.36655589E-02;
    COFD[7196] = -2.22371738E+01;
    COFD[7197] = 5.14170132E+00;
    COFD[7198] = -3.94484367E-01;
    COFD[7199] = 1.48469922E-02;
    COFD[7200] = -2.22970365E+01;
    COFD[7201] = 5.17327578E+00;
    COFD[7202] = -3.99835717E-01;
    COFD[7203] = 1.51252623E-02;
    COFD[7204] = -2.18473715E+01;
    COFD[7205] = 4.97429231E+00;
    COFD[7206] = -3.67512538E-01;
    COFD[7207] = 1.34835182E-02;
    COFD[7208] = -2.10762502E+01;
    COFD[7209] = 4.80548676E+00;
    COFD[7210] = -3.41059239E-01;
    COFD[7211] = 1.21698161E-02;
    COFD[7212] = -2.10815232E+01;
    COFD[7213] = 4.80548676E+00;
    COFD[7214] = -3.41059239E-01;
    COFD[7215] = 1.21698161E-02;
    COFD[7216] = -2.15475953E+01;
    COFD[7217] = 4.96375446E+00;
    COFD[7218] = -3.65860553E-01;
    COFD[7219] = 1.34014562E-02;
    COFD[7220] = -2.18709120E+01;
    COFD[7221] = 5.02398766E+00;
    COFD[7222] = -3.75343697E-01;
    COFD[7223] = 1.38737246E-02;
    COFD[7224] = -2.22125399E+01;
    COFD[7225] = 5.39199985E+00;
    COFD[7226] = -4.37219136E-01;
    COFD[7227] = 1.70801890E-02;
    COFD[7228] = -2.21849394E+01;
    COFD[7229] = 5.36809191E+00;
    COFD[7230] = -4.32908721E-01;
    COFD[7231] = 1.68481857E-02;
    COFD[7232] = -2.16580095E+01;
    COFD[7233] = 4.87709635E+00;
    COFD[7234] = -3.52250416E-01;
    COFD[7235] = 1.27244832E-02;
    COFD[7236] = -2.18230852E+01;
    COFD[7237] = 4.91762397E+00;
    COFD[7238] = -3.58634830E-01;
    COFD[7239] = 1.30425594E-02;
    COFD[7240] = -2.15517540E+01;
    COFD[7241] = 4.96375446E+00;
    COFD[7242] = -3.65860553E-01;
    COFD[7243] = 1.34014562E-02;
    COFD[7244] = -2.17527972E+01;
    COFD[7245] = 4.82487053E+00;
    COFD[7246] = -3.44081206E-01;
    COFD[7247] = 1.23193175E-02;
    COFD[7248] = -2.15669461E+01;
    COFD[7249] = 4.77388530E+00;
    COFD[7250] = -3.36183750E-01;
    COFD[7251] = 1.19302239E-02;
    COFD[7252] = -2.20221614E+01;
    COFD[7253] = 4.89616559E+00;
    COFD[7254] = -3.55241555E-01;
    COFD[7255] = 1.28731550E-02;
    COFD[7256] = -2.20205411E+01;
    COFD[7257] = 4.89616559E+00;
    COFD[7258] = -3.55241555E-01;
    COFD[7259] = 1.28731550E-02;
    COFD[7260] = -2.12201218E+01;
    COFD[7261] = 4.46261834E+00;
    COFD[7262] = -2.88453410E-01;
    COFD[7263] = 9.59378935E-03;
    COFD[7264] = -2.10625949E+01;
    COFD[7265] = 4.28688746E+00;
    COFD[7266] = -2.61952393E-01;
    COFD[7267] = 8.31175423E-03;
    COFD[7268] = -2.16041988E+01;
    COFD[7269] = 4.66858934E+00;
    COFD[7270] = -3.20108476E-01;
    COFD[7271] = 1.11452199E-02;
    COFD[7272] = -2.10871397E+01;
    COFD[7273] = 4.37219452E+00;
    COFD[7274] = -2.74784562E-01;
    COFD[7275] = 8.93131387E-03;
    COFD[7276] = -1.99112647E+01;
    COFD[7277] = 4.94145734E+00;
    COFD[7278] = -4.13478429E-01;
    COFD[7279] = 1.73756310E-02;
    COFD[7280] = -1.92371546E+01;
    COFD[7281] = 5.34516253E+00;
    COFD[7282] = -4.57548788E-01;
    COFD[7283] = 1.89625374E-02;
    COFD[7284] = -1.53695980E+01;
    COFD[7285] = 3.80403213E+00;
    COFD[7286] = -2.76736813E-01;
    COFD[7287] = 1.18808568E-02;
    COFD[7288] = -1.87927452E+01;
    COFD[7289] = 4.71430273E+00;
    COFD[7290] = -3.87319126E-01;
    COFD[7291] = 1.63728857E-02;
    COFD[7292] = -2.02674877E+01;
    COFD[7293] = 5.06295779E+00;
    COFD[7294] = -4.27579258E-01;
    COFD[7295] = 1.79215546E-02;
    COFD[7296] = -1.88154846E+01;
    COFD[7297] = 4.71430273E+00;
    COFD[7298] = -3.87319126E-01;
    COFD[7299] = 1.63728857E-02;
    COFD[7300] = -1.96977208E+01;
    COFD[7301] = 4.41245751E+00;
    COFD[7302] = -2.80855500E-01;
    COFD[7303] = 9.22505777E-03;
    COFD[7304] = -1.99461190E+01;
    COFD[7305] = 4.94937905E+00;
    COFD[7306] = -4.14409534E-01;
    COFD[7307] = 1.74121859E-02;
    COFD[7308] = -2.21152117E+01;
    COFD[7309] = 5.48782481E+00;
    COFD[7310] = -4.59567900E-01;
    COFD[7311] = 1.84190115E-02;
    COFD[7312] = -2.09794158E+01;
    COFD[7313] = 5.33874924E+00;
    COFD[7314] = -4.56879272E-01;
    COFD[7315] = 1.89397603E-02;
    COFD[7316] = -2.09268575E+01;
    COFD[7317] = 5.31839685E+00;
    COFD[7318] = -4.54623948E-01;
    COFD[7319] = 1.88567295E-02;
    COFD[7320] = -2.02767375E+01;
    COFD[7321] = 5.06295779E+00;
    COFD[7322] = -4.27579258E-01;
    COFD[7323] = 1.79215546E-02;
    COFD[7324] = -2.02855982E+01;
    COFD[7325] = 5.06295779E+00;
    COFD[7326] = -4.27579258E-01;
    COFD[7327] = 1.79215546E-02;
    COFD[7328] = -2.08014167E+01;
    COFD[7329] = 4.73262810E+00;
    COFD[7330] = -3.29863658E-01;
    COFD[7331] = 1.16208974E-02;
    COFD[7332] = -2.21672319E+01;
    COFD[7333] = 5.48560730E+00;
    COFD[7334] = -4.58746417E-01;
    COFD[7335] = 1.83631968E-02;
    COFD[7336] = -2.18987166E+01;
    COFD[7337] = 5.49121786E+00;
    COFD[7338] = -4.60888505E-01;
    COFD[7339] = 1.85096783E-02;
    COFD[7340] = -2.21566783E+01;
    COFD[7341] = 5.48560730E+00;
    COFD[7342] = -4.58746417E-01;
    COFD[7343] = 1.83631968E-02;
    COFD[7344] = -2.19818415E+01;
    COFD[7345] = 5.47597765E+00;
    COFD[7346] = -4.55020526E-01;
    COFD[7347] = 1.81091205E-02;
    COFD[7348] = -2.19947318E+01;
    COFD[7349] = 5.47597765E+00;
    COFD[7350] = -4.55020526E-01;
    COFD[7351] = 1.81091205E-02;
    COFD[7352] = -2.09754168E+01;
    COFD[7353] = 4.80548676E+00;
    COFD[7354] = -3.41059239E-01;
    COFD[7355] = 1.21698161E-02;
    COFD[7356] = -2.15045977E+01;
    COFD[7357] = 4.99422177E+00;
    COFD[7358] = -3.70645569E-01;
    COFD[7359] = 1.36394075E-02;
    COFD[7360] = -2.10463583E+01;
    COFD[7361] = 5.37676571E+00;
    COFD[7362] = -4.60857509E-01;
    COFD[7363] = 1.90756062E-02;
    COFD[7364] = -2.15169519E+01;
    COFD[7365] = 4.99422177E+00;
    COFD[7366] = -3.70645569E-01;
    COFD[7367] = 1.36394075E-02;
    COFD[7368] = -2.21967084E+01;
    COFD[7369] = 5.36809191E+00;
    COFD[7370] = -4.32908721E-01;
    COFD[7371] = 1.68481857E-02;
    COFD[7372] = -2.21967084E+01;
    COFD[7373] = 5.36809191E+00;
    COFD[7374] = -4.32908721E-01;
    COFD[7375] = 1.68481857E-02;
    COFD[7376] = -2.22273939E+01;
    COFD[7377] = 5.41240338E+00;
    COFD[7378] = -4.41033923E-01;
    COFD[7379] = 1.72894624E-02;
    COFD[7380] = -2.22376174E+01;
    COFD[7381] = 5.27340167E+00;
    COFD[7382] = -4.16479301E-01;
    COFD[7383] = 1.59820820E-02;
    COFD[7384] = -2.22351815E+01;
    COFD[7385] = 5.27973457E+00;
    COFD[7386] = -4.17551922E-01;
    COFD[7387] = 1.60378887E-02;
    COFD[7388] = -2.22393829E+01;
    COFD[7389] = 5.27973457E+00;
    COFD[7390] = -4.17551922E-01;
    COFD[7391] = 1.60378887E-02;
    COFD[7392] = -2.23944792E+01;
    COFD[7393] = 5.28846317E+00;
    COFD[7394] = -4.19054070E-01;
    COFD[7395] = 1.61166909E-02;
    COFD[7396] = -2.18058201E+01;
    COFD[7397] = 4.99630533E+00;
    COFD[7398] = -3.70973757E-01;
    COFD[7399] = 1.36557558E-02;
    COFD[7400] = -2.17649377E+01;
    COFD[7401] = 4.99755422E+00;
    COFD[7402] = -3.71170522E-01;
    COFD[7403] = 1.36655589E-02;
    COFD[7404] = -2.22434560E+01;
    COFD[7405] = 5.14170132E+00;
    COFD[7406] = -3.94484367E-01;
    COFD[7407] = 1.48469922E-02;
    COFD[7408] = -2.23033557E+01;
    COFD[7409] = 5.17327578E+00;
    COFD[7410] = -3.99835717E-01;
    COFD[7411] = 1.51252623E-02;
    COFD[7412] = -2.18537271E+01;
    COFD[7413] = 4.97429231E+00;
    COFD[7414] = -3.67512538E-01;
    COFD[7415] = 1.34835182E-02;
    COFD[7416] = -2.10815232E+01;
    COFD[7417] = 4.80548676E+00;
    COFD[7418] = -3.41059239E-01;
    COFD[7419] = 1.21698161E-02;
    COFD[7420] = -2.10868524E+01;
    COFD[7421] = 4.80548676E+00;
    COFD[7422] = -3.41059239E-01;
    COFD[7423] = 1.21698161E-02;
    COFD[7424] = -2.15532850E+01;
    COFD[7425] = 4.96375446E+00;
    COFD[7426] = -3.65860553E-01;
    COFD[7427] = 1.34014562E-02;
    COFD[7428] = -2.18766495E+01;
    COFD[7429] = 5.02398766E+00;
    COFD[7430] = -3.75343697E-01;
    COFD[7431] = 1.38737246E-02;
    COFD[7432] = -2.22174544E+01;
    COFD[7433] = 5.39199985E+00;
    COFD[7434] = -4.37219136E-01;
    COFD[7435] = 1.70801890E-02;
    COFD[7436] = -2.21897214E+01;
    COFD[7437] = 5.36809191E+00;
    COFD[7438] = -4.32908721E-01;
    COFD[7439] = 1.68481857E-02;
    COFD[7440] = -2.16643993E+01;
    COFD[7441] = 4.87709635E+00;
    COFD[7442] = -3.52250416E-01;
    COFD[7443] = 1.27244832E-02;
    COFD[7444] = -2.18294392E+01;
    COFD[7445] = 4.91762397E+00;
    COFD[7446] = -3.58634830E-01;
    COFD[7447] = 1.30425594E-02;
    COFD[7448] = -2.15574915E+01;
    COFD[7449] = 4.96375446E+00;
    COFD[7450] = -3.65860553E-01;
    COFD[7451] = 1.34014562E-02;
    COFD[7452] = -2.17596291E+01;
    COFD[7453] = 4.82487053E+00;
    COFD[7454] = -3.44081206E-01;
    COFD[7455] = 1.23193175E-02;
    COFD[7456] = -2.15737493E+01;
    COFD[7457] = 4.77388530E+00;
    COFD[7458] = -3.36183750E-01;
    COFD[7459] = 1.19302239E-02;
    COFD[7460] = -2.20293523E+01;
    COFD[7461] = 4.89616559E+00;
    COFD[7462] = -3.55241555E-01;
    COFD[7463] = 1.28731550E-02;
    COFD[7464] = -2.20277086E+01;
    COFD[7465] = 4.89616559E+00;
    COFD[7466] = -3.55241555E-01;
    COFD[7467] = 1.28731550E-02;
    COFD[7468] = -2.12279042E+01;
    COFD[7469] = 4.46261834E+00;
    COFD[7470] = -2.88453410E-01;
    COFD[7471] = 9.59378935E-03;
    COFD[7472] = -2.10708054E+01;
    COFD[7473] = 4.28688746E+00;
    COFD[7474] = -2.61952393E-01;
    COFD[7475] = 8.31175423E-03;
    COFD[7476] = -2.16116846E+01;
    COFD[7477] = 4.66858934E+00;
    COFD[7478] = -3.20108476E-01;
    COFD[7479] = 1.11452199E-02;
    COFD[7480] = -2.10951403E+01;
    COFD[7481] = 4.37219452E+00;
    COFD[7482] = -2.74784562E-01;
    COFD[7483] = 8.93131387E-03;
    COFD[7484] = -1.99151960E+01;
    COFD[7485] = 4.94145734E+00;
    COFD[7486] = -4.13478429E-01;
    COFD[7487] = 1.73756310E-02;
    COFD[7488] = -1.91106350E+01;
    COFD[7489] = 5.24865128E+00;
    COFD[7490] = -4.46802764E-01;
    COFD[7491] = 1.85649476E-02;
    COFD[7492] = -1.52493207E+01;
    COFD[7493] = 3.70713307E+00;
    COFD[7494] = -2.64771252E-01;
    COFD[7495] = 1.13878078E-02;
    COFD[7496] = -1.85962044E+01;
    COFD[7497] = 4.58430295E+00;
    COFD[7498] = -3.71395187E-01;
    COFD[7499] = 1.57196293E-02;
    COFD[7500] = -2.01222360E+01;
    COFD[7501] = 4.95889869E+00;
    COFD[7502] = -4.15526934E-01;
    COFD[7503] = 1.74559870E-02;
    COFD[7504] = -1.86197092E+01;
    COFD[7505] = 4.58430295E+00;
    COFD[7506] = -3.71395187E-01;
    COFD[7507] = 1.57196293E-02;
    COFD[7508] = -2.02528193E+01;
    COFD[7509] = 4.60540477E+00;
    COFD[7510] = -3.10373488E-01;
    COFD[7511] = 1.06672049E-02;
    COFD[7512] = -1.98245633E+01;
    COFD[7513] = 4.85730922E+00;
    COFD[7514] = -4.04004009E-01;
    COFD[7515] = 1.70221852E-02;
    COFD[7516] = -2.22150103E+01;
    COFD[7517] = 5.49406269E+00;
    COFD[7518] = -4.63519837E-01;
    COFD[7519] = 1.87101015E-02;
    COFD[7520] = -2.08233373E+01;
    COFD[7521] = 5.23993213E+00;
    COFD[7522] = -4.45818789E-01;
    COFD[7523] = 1.85280084E-02;
    COFD[7524] = -2.07627166E+01;
    COFD[7525] = 5.21569495E+00;
    COFD[7526] = -4.43054270E-01;
    COFD[7527] = 1.84229515E-02;
    COFD[7528] = -2.01319872E+01;
    COFD[7529] = 4.95889869E+00;
    COFD[7530] = -4.15526934E-01;
    COFD[7531] = 1.74559870E-02;
    COFD[7532] = -2.01413374E+01;
    COFD[7533] = 4.95889869E+00;
    COFD[7534] = -4.15526934E-01;
    COFD[7535] = 1.74559870E-02;
    COFD[7536] = -2.12869434E+01;
    COFD[7537] = 4.90168525E+00;
    COFD[7538] = -3.56120032E-01;
    COFD[7539] = 1.29171867E-02;
    COFD[7540] = -2.22550223E+01;
    COFD[7541] = 5.49328085E+00;
    COFD[7542] = -4.62887295E-01;
    COFD[7543] = 1.86623580E-02;
    COFD[7544] = -2.19808137E+01;
    COFD[7545] = 5.49090680E+00;
    COFD[7546] = -4.63910285E-01;
    COFD[7547] = 1.87574539E-02;
    COFD[7548] = -2.22439303E+01;
    COFD[7549] = 5.49328085E+00;
    COFD[7550] = -4.62887295E-01;
    COFD[7551] = 1.86623580E-02;
    COFD[7552] = -2.20896596E+01;
    COFD[7553] = 5.48779691E+00;
    COFD[7554] = -4.59557311E-01;
    COFD[7555] = 1.84182884E-02;
    COFD[7556] = -2.21031461E+01;
    COFD[7557] = 5.48779691E+00;
    COFD[7558] = -4.59557311E-01;
    COFD[7559] = 1.84182884E-02;
    COFD[7560] = -2.14350985E+01;
    COFD[7561] = 4.96375446E+00;
    COFD[7562] = -3.65860553E-01;
    COFD[7563] = 1.34014562E-02;
    COFD[7564] = -2.19030787E+01;
    COFD[7565] = 5.12829880E+00;
    COFD[7566] = -3.92242339E-01;
    COFD[7567] = 1.47312337E-02;
    COFD[7568] = -2.09608951E+01;
    COFD[7569] = 5.29024160E+00;
    COFD[7570] = -4.51483508E-01;
    COFD[7571] = 1.87402157E-02;
    COFD[7572] = -2.19162243E+01;
    COFD[7573] = 5.12829880E+00;
    COFD[7574] = -3.92242339E-01;
    COFD[7575] = 1.47312337E-02;
    COFD[7576] = -2.24246556E+01;
    COFD[7577] = 5.43317403E+00;
    COFD[7578] = -4.45039640E-01;
    COFD[7579] = 1.75123239E-02;
    COFD[7580] = -2.24246556E+01;
    COFD[7581] = 5.43317403E+00;
    COFD[7582] = -4.45039640E-01;
    COFD[7583] = 1.75123239E-02;
    COFD[7584] = -2.24102868E+01;
    COFD[7585] = 5.45651702E+00;
    COFD[7586] = -4.50175370E-01;
    COFD[7587] = 1.78134258E-02;
    COFD[7588] = -2.25142044E+01;
    COFD[7589] = 5.35783672E+00;
    COFD[7590] = -4.31095103E-01;
    COFD[7591] = 1.67515802E-02;
    COFD[7592] = -2.25103101E+01;
    COFD[7593] = 5.36306078E+00;
    COFD[7594] = -4.32018705E-01;
    COFD[7595] = 1.68007664E-02;
    COFD[7596] = -2.25148197E+01;
    COFD[7597] = 5.36306078E+00;
    COFD[7598] = -4.32018705E-01;
    COFD[7599] = 1.68007664E-02;
    COFD[7600] = -2.26606931E+01;
    COFD[7601] = 5.37102704E+00;
    COFD[7602] = -4.33426866E-01;
    COFD[7603] = 1.68757765E-02;
    COFD[7604] = -2.21956028E+01;
    COFD[7605] = 5.12996650E+00;
    COFD[7606] = -3.92516016E-01;
    COFD[7607] = 1.47452270E-02;
    COFD[7608] = -2.21570195E+01;
    COFD[7609] = 5.13096622E+00;
    COFD[7610] = -3.92680132E-01;
    COFD[7611] = 1.47536199E-02;
    COFD[7612] = -2.25643461E+01;
    COFD[7613] = 5.24758220E+00;
    COFD[7614] = -4.12129116E-01;
    COFD[7615] = 1.57564573E-02;
    COFD[7616] = -2.26246057E+01;
    COFD[7617] = 5.27897272E+00;
    COFD[7618] = -4.17421132E-01;
    COFD[7619] = 1.60310357E-02;
    COFD[7620] = -2.22571178E+01;
    COFD[7621] = 5.11286874E+00;
    COFD[7622] = -3.89706255E-01;
    COFD[7623] = 1.46014745E-02;
    COFD[7624] = -2.15475953E+01;
    COFD[7625] = 4.96375446E+00;
    COFD[7626] = -3.65860553E-01;
    COFD[7627] = 1.34014562E-02;
    COFD[7628] = -2.15532850E+01;
    COFD[7629] = 4.96375446E+00;
    COFD[7630] = -3.65860553E-01;
    COFD[7631] = 1.34014562E-02;
    COFD[7632] = -2.19714779E+01;
    COFD[7633] = 5.10490214E+00;
    COFD[7634] = -3.88397529E-01;
    COFD[7635] = 1.45345306E-02;
    COFD[7636] = -2.22437818E+01;
    COFD[7637] = 5.14893259E+00;
    COFD[7638] = -3.95700512E-01;
    COFD[7639] = 1.49099774E-02;
    COFD[7640] = -2.24203640E+01;
    COFD[7641] = 5.44615881E+00;
    COFD[7642] = -4.47865005E-01;
    COFD[7643] = 1.76775502E-02;
    COFD[7644] = -2.24172422E+01;
    COFD[7645] = 5.43317403E+00;
    COFD[7646] = -4.45039640E-01;
    COFD[7647] = 1.75123239E-02;
    COFD[7648] = -2.21073752E+01;
    COFD[7649] = 5.03131423E+00;
    COFD[7650] = -3.76503291E-01;
    COFD[7651] = 1.39316539E-02;
    COFD[7652] = -2.22582230E+01;
    COFD[7653] = 5.06826700E+00;
    COFD[7654] = -3.82432503E-01;
    COFD[7655] = 1.42308185E-02;
    COFD[7656] = -2.19759929E+01;
    COFD[7657] = 5.10490214E+00;
    COFD[7658] = -3.88397529E-01;
    COFD[7659] = 1.45345306E-02;
    COFD[7660] = -2.22053962E+01;
    COFD[7661] = 4.98256919E+00;
    COFD[7662] = -3.68811981E-01;
    COFD[7663] = 1.35481236E-02;
    COFD[7664] = -2.20318339E+01;
    COFD[7665] = 4.93482700E+00;
    COFD[7666] = -3.61331756E-01;
    COFD[7667] = 1.31765767E-02;
    COFD[7668] = -2.24600223E+01;
    COFD[7669] = 5.04880598E+00;
    COFD[7670] = -3.79285322E-01;
    COFD[7671] = 1.40712210E-02;
    COFD[7672] = -2.24582240E+01;
    COFD[7673] = 5.04880598E+00;
    COFD[7674] = -3.79285322E-01;
    COFD[7675] = 1.40712210E-02;
    COFD[7676] = -2.17388998E+01;
    COFD[7677] = 4.64756828E+00;
    COFD[7678] = -3.16905849E-01;
    COFD[7679] = 1.09890709E-02;
    COFD[7680] = -2.16429075E+01;
    COFD[7681] = 4.49834811E+00;
    COFD[7682] = -2.93898563E-01;
    COFD[7683] = 9.85901966E-03;
    COFD[7684] = -2.21166468E+01;
    COFD[7685] = 4.85109413E+00;
    COFD[7686] = -3.48178913E-01;
    COFD[7687] = 1.25223681E-02;
    COFD[7688] = -2.16405632E+01;
    COFD[7689] = 4.57106634E+00;
    COFD[7690] = -3.05074221E-01;
    COFD[7691] = 1.04066420E-02;
    COFD[7692] = -1.97977739E+01;
    COFD[7693] = 4.85079945E+00;
    COFD[7694] = -4.03257284E-01;
    COFD[7695] = 1.69936877E-02;
    COFD[7696] = -1.92438114E+01;
    COFD[7697] = 5.20849012E+00;
    COFD[7698] = -4.42302201E-01;
    COFD[7699] = 1.83974704E-02;
    COFD[7700] = -1.57696917E+01;
    COFD[7701] = 3.80988496E+00;
    COFD[7702] = -2.77445338E-01;
    COFD[7703] = 1.19094487E-02;
    COFD[7704] = -1.87089834E+01;
    COFD[7705] = 4.54237822E+00;
    COFD[7706] = -3.66363712E-01;
    COFD[7707] = 1.55175729E-02;
    COFD[7708] = -2.05600561E+01;
    COFD[7709] = 5.03374066E+00;
    COFD[7710] = -4.24314609E-01;
    COFD[7711] = 1.78006508E-02;
    COFD[7712] = -1.87325861E+01;
    COFD[7713] = 4.54237822E+00;
    COFD[7714] = -3.66363712E-01;
    COFD[7715] = 1.55175729E-02;
    COFD[7716] = -1.86106672E+01;
    COFD[7717] = 3.70369165E+00;
    COFD[7718] = -1.75747020E-01;
    COFD[7719] = 4.23774985E-03;
    COFD[7720] = -2.02385185E+01;
    COFD[7721] = 4.92216751E+00;
    COFD[7722] = -4.11213729E-01;
    COFD[7723] = 1.72868483E-02;
    COFD[7724] = -2.24591465E+01;
    COFD[7725] = 5.49093616E+00;
    COFD[7726] = -4.61489559E-01;
    COFD[7727] = 1.85600323E-02;
    COFD[7728] = -2.09160351E+01;
    COFD[7729] = 5.20438984E+00;
    COFD[7730] = -4.41945569E-01;
    COFD[7731] = 1.83888515E-02;
    COFD[7732] = -2.12249561E+01;
    COFD[7733] = 5.29857429E+00;
    COFD[7734] = -4.52417024E-01;
    COFD[7735] = 1.87750282E-02;
    COFD[7736] = -2.02137759E+01;
    COFD[7737] = 4.91354418E+00;
    COFD[7738] = -4.10191969E-01;
    COFD[7739] = 1.72463969E-02;
    COFD[7740] = -2.02231901E+01;
    COFD[7741] = 4.91354418E+00;
    COFD[7742] = -4.10191969E-01;
    COFD[7743] = 1.72463969E-02;
    COFD[7744] = -2.16143949E+01;
    COFD[7745] = 4.96367791E+00;
    COFD[7746] = -3.65848555E-01;
    COFD[7747] = 1.34008603E-02;
    COFD[7748] = -2.24142003E+01;
    COFD[7749] = 5.49189132E+00;
    COFD[7750] = -4.63888808E-01;
    COFD[7751] = 1.87508465E-02;
    COFD[7752] = -2.21435675E+01;
    COFD[7753] = 5.48236053E+00;
    COFD[7754] = -4.63887012E-01;
    COFD[7755] = 1.87976815E-02;
    COFD[7756] = -2.24030382E+01;
    COFD[7757] = 5.49189132E+00;
    COFD[7758] = -4.63888808E-01;
    COFD[7759] = 1.87508465E-02;
    COFD[7760] = -2.22736645E+01;
    COFD[7761] = 5.49118605E+00;
    COFD[7762] = -4.61241736E-01;
    COFD[7763] = 1.85386070E-02;
    COFD[7764] = -2.22872282E+01;
    COFD[7765] = 5.49118605E+00;
    COFD[7766] = -4.61241736E-01;
    COFD[7767] = 1.85386070E-02;
    COFD[7768] = -2.17575745E+01;
    COFD[7769] = 5.02398766E+00;
    COFD[7770] = -3.75343697E-01;
    COFD[7771] = 1.38737246E-02;
    COFD[7772] = -2.21678146E+01;
    COFD[7773] = 5.16786857E+00;
    COFD[7774] = -3.98913959E-01;
    COFD[7775] = 1.50771754E-02;
    COFD[7776] = -2.10807241E+01;
    COFD[7777] = 5.24972538E+00;
    COFD[7778] = -4.46923248E-01;
    COFD[7779] = 1.85694358E-02;
    COFD[7780] = -2.21810647E+01;
    COFD[7781] = 5.16786857E+00;
    COFD[7782] = -3.98913959E-01;
    COFD[7783] = 1.50771754E-02;
    COFD[7784] = -2.26244858E+01;
    COFD[7785] = 5.44849005E+00;
    COFD[7786] = -4.48374875E-01;
    COFD[7787] = 1.77073338E-02;
    COFD[7788] = -2.26244858E+01;
    COFD[7789] = 5.44849005E+00;
    COFD[7790] = -4.48374875E-01;
    COFD[7791] = 1.77073338E-02;
    COFD[7792] = -2.26076410E+01;
    COFD[7793] = 5.46824702E+00;
    COFD[7794] = -4.52880728E-01;
    COFD[7795] = 1.79742471E-02;
    COFD[7796] = -2.27416253E+01;
    COFD[7797] = 5.38802436E+00;
    COFD[7798] = -4.36467632E-01;
    COFD[7799] = 1.70387872E-02;
    COFD[7800] = -2.27363406E+01;
    COFD[7801] = 5.39211717E+00;
    COFD[7802] = -4.37241388E-01;
    COFD[7803] = 1.70814168E-02;
    COFD[7804] = -2.27408914E+01;
    COFD[7805] = 5.39211717E+00;
    COFD[7806] = -4.37241388E-01;
    COFD[7807] = 1.70814168E-02;
    COFD[7808] = -2.28720058E+01;
    COFD[7809] = 5.39857883E+00;
    COFD[7810] = -4.38452012E-01;
    COFD[7811] = 1.71478655E-02;
    COFD[7812] = -2.24415132E+01;
    COFD[7813] = 5.16918377E+00;
    COFD[7814] = -3.99137871E-01;
    COFD[7815] = 1.50888488E-02;
    COFD[7816] = -2.10421374E+01;
    COFD[7817] = 4.51136123E+00;
    COFD[7818] = -3.00118189E-01;
    COFD[7819] = 1.03138650E-02;
    COFD[7820] = -2.27971335E+01;
    COFD[7821] = 5.28513609E+00;
    COFD[7822] = -4.18480689E-01;
    COFD[7823] = 1.60865907E-02;
    COFD[7824] = -2.27232104E+01;
    COFD[7825] = 5.26126980E+00;
    COFD[7826] = -4.16242057E-01;
    COFD[7827] = 1.60218084E-02;
    COFD[7828] = -2.25079550E+01;
    COFD[7829] = 5.15543496E+00;
    COFD[7830] = -3.96801390E-01;
    COFD[7831] = 1.49671794E-02;
    COFD[7832] = -2.18709120E+01;
    COFD[7833] = 5.02398766E+00;
    COFD[7834] = -3.75343697E-01;
    COFD[7835] = 1.38737246E-02;
    COFD[7836] = -2.18766495E+01;
    COFD[7837] = 5.02398766E+00;
    COFD[7838] = -3.75343697E-01;
    COFD[7839] = 1.38737246E-02;
    COFD[7840] = -2.22437818E+01;
    COFD[7841] = 5.14893259E+00;
    COFD[7842] = -3.95700512E-01;
    COFD[7843] = 1.49099774E-02;
    COFD[7844] = -2.11492935E+01;
    COFD[7845] = 4.54244419E+00;
    COFD[7846] = -3.05131457E-01;
    COFD[7847] = 1.05659771E-02;
    COFD[7848] = -2.26176957E+01;
    COFD[7849] = 5.45876240E+00;
    COFD[7850] = -4.50706326E-01;
    COFD[7851] = 1.78452467E-02;
    COFD[7852] = -2.26170162E+01;
    COFD[7853] = 5.44849005E+00;
    COFD[7854] = -4.48374875E-01;
    COFD[7855] = 1.77073338E-02;
    COFD[7856] = -2.11645228E+01;
    COFD[7857] = 4.51897995E+00;
    COFD[7858] = -3.01998503E-01;
    COFD[7859] = 1.04247720E-02;
    COFD[7860] = -2.11597907E+01;
    COFD[7861] = 4.47157782E+00;
    COFD[7862] = -2.94123466E-01;
    COFD[7863] = 1.00236418E-02;
    COFD[7864] = -2.22483379E+01;
    COFD[7865] = 5.14893259E+00;
    COFD[7866] = -3.95700512E-01;
    COFD[7867] = 1.49099774E-02;
    COFD[7868] = -2.15608572E+01;
    COFD[7869] = 4.62655831E+00;
    COFD[7870] = -3.18380357E-01;
    COFD[7871] = 1.12180276E-02;
    COFD[7872] = -2.23257344E+01;
    COFD[7873] = 4.99720991E+00;
    COFD[7874] = -3.71116272E-01;
    COFD[7875] = 1.36628560E-02;
    COFD[7876] = -2.27123285E+01;
    COFD[7877] = 5.10007015E+00;
    COFD[7878] = -3.87606161E-01;
    COFD[7879] = 1.44941102E-02;
    COFD[7880] = -2.27105091E+01;
    COFD[7881] = 5.10007015E+00;
    COFD[7882] = -3.87606161E-01;
    COFD[7883] = 1.44941102E-02;
    COFD[7884] = -2.16032600E+01;
    COFD[7885] = 4.53358880E+00;
    COFD[7886] = -3.02280026E-01;
    COFD[7887] = 1.03666188E-02;
    COFD[7888] = -2.15585553E+01;
    COFD[7889] = 4.41358685E+00;
    COFD[7890] = -2.83552193E-01;
    COFD[7891] = 9.43843173E-03;
    COFD[7892] = -2.23943879E+01;
    COFD[7893] = 4.91274983E+00;
    COFD[7894] = -3.57867526E-01;
    COFD[7895] = 1.30043584E-02;
    COFD[7896] = -2.14143742E+01;
    COFD[7897] = 4.41973109E+00;
    COFD[7898] = -2.85282560E-01;
    COFD[7899] = 9.54967542E-03;
    COFD[7900] = -2.01834870E+01;
    COFD[7901] = 4.90590035E+00;
    COFD[7902] = -4.09326358E-01;
    COFD[7903] = 1.72138685E-02;
    COFD[7904] = -1.80817998E+01;
    COFD[7905] = 4.90738815E+00;
    COFD[7906] = -4.09483894E-01;
    COFD[7907] = 1.72192910E-02;
    COFD[7908] = -1.42128070E+01;
    COFD[7909] = 3.31879365E+00;
    COFD[7910] = -2.16004062E-01;
    COFD[7911] = 9.34330595E-03;
    COFD[7912] = -1.73332914E+01;
    COFD[7913] = 4.13619942E+00;
    COFD[7914] = -3.16693330E-01;
    COFD[7915] = 1.34845235E-02;
    COFD[7916] = -1.88511274E+01;
    COFD[7917] = 4.53065425E+00;
    COFD[7918] = -3.64983473E-01;
    COFD[7919] = 1.54632547E-02;
    COFD[7920] = -1.73550857E+01;
    COFD[7921] = 4.13619942E+00;
    COFD[7922] = -3.16693330E-01;
    COFD[7923] = 1.34845235E-02;
    COFD[7924] = -2.13752043E+01;
    COFD[7925] = 5.21010770E+00;
    COFD[7926] = -4.05945421E-01;
    COFD[7927] = 1.54393383E-02;
    COFD[7928] = -1.85705503E+01;
    COFD[7929] = 4.42984539E+00;
    COFD[7930] = -3.53041284E-01;
    COFD[7931] = 1.49902085E-02;
    COFD[7932] = -2.17175223E+01;
    COFD[7933] = 5.42282789E+00;
    COFD[7934] = -4.65108864E-01;
    COFD[7935] = 1.91933632E-02;
    COFD[7936] = -1.97629737E+01;
    COFD[7937] = 4.90161091E+00;
    COFD[7938] = -4.08878655E-01;
    COFD[7939] = 1.71987880E-02;
    COFD[7940] = -1.97245207E+01;
    COFD[7941] = 4.88654569E+00;
    COFD[7942] = -4.07307743E-01;
    COFD[7943] = 1.71459608E-02;
    COFD[7944] = -1.88597846E+01;
    COFD[7945] = 4.53065425E+00;
    COFD[7946] = -3.64983473E-01;
    COFD[7947] = 1.54632547E-02;
    COFD[7948] = -1.88680678E+01;
    COFD[7949] = 4.53065425E+00;
    COFD[7950] = -3.64983473E-01;
    COFD[7951] = 1.54632547E-02;
    COFD[7952] = -2.20470531E+01;
    COFD[7953] = 5.36287304E+00;
    COFD[7954] = -4.31985459E-01;
    COFD[7955] = 1.67989944E-02;
    COFD[7956] = -2.17838478E+01;
    COFD[7957] = 5.42865606E+00;
    COFD[7958] = -4.65447368E-01;
    COFD[7959] = 1.91919264E-02;
    COFD[7960] = -2.14802351E+01;
    COFD[7961] = 5.41058359E+00;
    COFD[7962] = -4.64178328E-01;
    COFD[7963] = 1.91784294E-02;
    COFD[7964] = -2.17739349E+01;
    COFD[7965] = 5.42865606E+00;
    COFD[7966] = -4.65447368E-01;
    COFD[7967] = 1.91919264E-02;
    COFD[7968] = -2.16714153E+01;
    COFD[7969] = 5.44455741E+00;
    COFD[7970] = -4.65313340E-01;
    COFD[7971] = 1.91011181E-02;
    COFD[7972] = -2.16835888E+01;
    COFD[7973] = 5.44455741E+00;
    COFD[7974] = -4.65313340E-01;
    COFD[7975] = 1.91011181E-02;
    COFD[7976] = -2.21138917E+01;
    COFD[7977] = 5.39199985E+00;
    COFD[7978] = -4.37219136E-01;
    COFD[7979] = 1.70801890E-02;
    COFD[7980] = -2.23238237E+01;
    COFD[7981] = 5.45282659E+00;
    COFD[7982] = -4.49336536E-01;
    COFD[7983] = 1.77637787E-02;
    COFD[7984] = -1.98725502E+01;
    COFD[7985] = 4.95323342E+00;
    COFD[7986] = -4.14861803E-01;
    COFD[7987] = 1.74299070E-02;
    COFD[7988] = -2.23352610E+01;
    COFD[7989] = 5.45282659E+00;
    COFD[7990] = -4.49336536E-01;
    COFD[7991] = 1.77637787E-02;
    COFD[7992] = -2.22041980E+01;
    COFD[7993] = 5.48376307E+00;
    COFD[7994] = -4.63942837E-01;
    COFD[7995] = 1.87953234E-02;
    COFD[7996] = -2.22041980E+01;
    COFD[7997] = 5.48376307E+00;
    COFD[7998] = -4.63942837E-01;
    COFD[7999] = 1.87953234E-02;
    COFD[8000] = -2.20957576E+01;
    COFD[8001] = 5.46736313E+00;
    COFD[8002] = -4.63533201E-01;
    COFD[8003] = 1.88428681E-02;
    COFD[8004] = -2.24743159E+01;
    COFD[8005] = 5.49090532E+00;
    COFD[8006] = -4.61515055E-01;
    COFD[8007] = 1.85622624E-02;
    COFD[8008] = -2.24568662E+01;
    COFD[8009] = 5.49068313E+00;
    COFD[8010] = -4.61694593E-01;
    COFD[8011] = 1.85779923E-02;
    COFD[8012] = -2.24607172E+01;
    COFD[8013] = 5.49068313E+00;
    COFD[8014] = -4.61694593E-01;
    COFD[8015] = 1.85779923E-02;
    COFD[8016] = -2.25843967E+01;
    COFD[8017] = 5.49024820E+00;
    COFD[8018] = -4.61953415E-01;
    COFD[8019] = 1.86012452E-02;
    COFD[8020] = -2.25999700E+01;
    COFD[8021] = 5.45327166E+00;
    COFD[8022] = -4.49436298E-01;
    COFD[8023] = 1.77696558E-02;
    COFD[8024] = -2.25595681E+01;
    COFD[8025] = 5.45353757E+00;
    COFD[8026] = -4.49496004E-01;
    COFD[8027] = 1.77731752E-02;
    COFD[8028] = -2.27488997E+01;
    COFD[8029] = 5.47944787E+00;
    COFD[8030] = -4.56152622E-01;
    COFD[8031] = 1.81835823E-02;
    COFD[8032] = -2.27431590E+01;
    COFD[8033] = 5.48304850E+00;
    COFD[8034] = -4.57579455E-01;
    COFD[8035] = 1.82812528E-02;
    COFD[8036] = -2.26814510E+01;
    COFD[8037] = 5.44849592E+00;
    COFD[8038] = -4.48376165E-01;
    COFD[8039] = 1.77074093E-02;
    COFD[8040] = -2.22125399E+01;
    COFD[8041] = 5.39199985E+00;
    COFD[8042] = -4.37219136E-01;
    COFD[8043] = 1.70801890E-02;
    COFD[8044] = -2.22174544E+01;
    COFD[8045] = 5.39199985E+00;
    COFD[8046] = -4.37219136E-01;
    COFD[8047] = 1.70801890E-02;
    COFD[8048] = -2.24203640E+01;
    COFD[8049] = 5.44615881E+00;
    COFD[8050] = -4.47865005E-01;
    COFD[8051] = 1.76775502E-02;
    COFD[8052] = -2.26176957E+01;
    COFD[8053] = 5.45876240E+00;
    COFD[8054] = -4.50706326E-01;
    COFD[8055] = 1.78452467E-02;
    COFD[8056] = -2.21475192E+01;
    COFD[8057] = 5.47416855E+00;
    COFD[8058] = -4.63563799E-01;
    COFD[8059] = 1.88117222E-02;
    COFD[8060] = -2.21977082E+01;
    COFD[8061] = 5.48376307E+00;
    COFD[8062] = -4.63942837E-01;
    COFD[8063] = 1.87953234E-02;
    COFD[8064] = -2.26581212E+01;
    COFD[8065] = 5.42030479E+00;
    COFD[8066] = -4.42541204E-01;
    COFD[8067] = 1.73729052E-02;
    COFD[8068] = -2.27561195E+01;
    COFD[8069] = 5.43540157E+00;
    COFD[8070] = -4.45477820E-01;
    COFD[8071] = 1.75369197E-02;
    COFD[8072] = -2.24242197E+01;
    COFD[8073] = 5.44615881E+00;
    COFD[8074] = -4.47865005E-01;
    COFD[8075] = 1.76775502E-02;
    COFD[8076] = -2.28163992E+01;
    COFD[8077] = 5.39961965E+00;
    COFD[8078] = -4.38647062E-01;
    COFD[8079] = 1.71585740E-02;
    COFD[8080] = -2.27109978E+01;
    COFD[8081] = 5.37975560E+00;
    COFD[8082] = -4.34974624E-01;
    COFD[8083] = 1.69583701E-02;
    COFD[8084] = -2.29717221E+01;
    COFD[8085] = 5.42721943E+00;
    COFD[8086] = -4.43876872E-01;
    COFD[8087] = 1.74472689E-02;
    COFD[8088] = -2.29702473E+01;
    COFD[8089] = 5.42721943E+00;
    COFD[8090] = -4.43876872E-01;
    COFD[8091] = 1.74472689E-02;
    COFD[8092] = -2.27497015E+01;
    COFD[8093] = 5.23793490E+00;
    COFD[8094] = -4.10507434E-01;
    COFD[8095] = 1.56723961E-02;
    COFD[8096] = -2.27937067E+01;
    COFD[8097] = 5.15050798E+00;
    COFD[8098] = -3.95966829E-01;
    COFD[8099] = 1.49238043E-02;
    COFD[8100] = -2.28830519E+01;
    COFD[8101] = 5.33738936E+00;
    COFD[8102] = -4.27502526E-01;
    COFD[8103] = 1.65608197E-02;
    COFD[8104] = -2.27063862E+01;
    COFD[8105] = 5.18619255E+00;
    COFD[8106] = -4.02090863E-01;
    COFD[8107] = 1.52443533E-02;
    COFD[8108] = -1.85439247E+01;
    COFD[8109] = 4.42312522E+00;
    COFD[8110] = -3.52239033E-01;
    COFD[8111] = 1.49581799E-02;
    COFD[8112] = -1.82112771E+01;
    COFD[8113] = 4.94453558E+00;
    COFD[8114] = -4.13840584E-01;
    COFD[8115] = 1.73898654E-02;
    COFD[8116] = -1.43061148E+01;
    COFD[8117] = 3.34609996E+00;
    COFD[8118] = -2.19417637E-01;
    COFD[8119] = 9.48565577E-03;
    COFD[8120] = -1.74487550E+01;
    COFD[8121] = 4.17306000E+00;
    COFD[8122] = -3.21208070E-01;
    COFD[8123] = 1.36692112E-02;
    COFD[8124] = -1.89456139E+01;
    COFD[8125] = 4.56070745E+00;
    COFD[8126] = -3.68517396E-01;
    COFD[8127] = 1.56021102E-02;
    COFD[8128] = -1.74702316E+01;
    COFD[8129] = 4.17306000E+00;
    COFD[8130] = -3.21208070E-01;
    COFD[8131] = 1.36692112E-02;
    COFD[8132] = -2.13273839E+01;
    COFD[8133] = 5.17422070E+00;
    COFD[8134] = -3.99997068E-01;
    COFD[8135] = 1.51336871E-02;
    COFD[8136] = -1.86697931E+01;
    COFD[8137] = 4.46203550E+00;
    COFD[8138] = -3.56905086E-01;
    COFD[8139] = 1.51454720E-02;
    COFD[8140] = -2.17687426E+01;
    COFD[8141] = 5.43321422E+00;
    COFD[8142] = -4.65643875E-01;
    COFD[8143] = 1.91852144E-02;
    COFD[8144] = -1.98726217E+01;
    COFD[8145] = 4.93512898E+00;
    COFD[8146] = -4.12735698E-01;
    COFD[8147] = 1.73465233E-02;
    COFD[8148] = -1.98084643E+01;
    COFD[8149] = 4.90949838E+00;
    COFD[8150] = -4.09710987E-01;
    COFD[8151] = 1.72272850E-02;
    COFD[8152] = -1.89540782E+01;
    COFD[8153] = 4.56070745E+00;
    COFD[8154] = -3.68517396E-01;
    COFD[8155] = 1.56021102E-02;
    COFD[8156] = -1.89621737E+01;
    COFD[8157] = 4.56070745E+00;
    COFD[8158] = -3.68517396E-01;
    COFD[8159] = 1.56021102E-02;
    COFD[8160] = -2.20163926E+01;
    COFD[8161] = 5.33662478E+00;
    COFD[8162] = -4.27371229E-01;
    COFD[8163] = 1.65539339E-02;
    COFD[8164] = -2.18308969E+01;
    COFD[8165] = 5.43667692E+00;
    COFD[8166] = -4.65640530E-01;
    COFD[8167] = 1.91675299E-02;
    COFD[8168] = -2.15435139E+01;
    COFD[8169] = 5.42458950E+00;
    COFD[8170] = -4.65212987E-01;
    COFD[8171] = 1.91930766E-02;
    COFD[8172] = -2.18211936E+01;
    COFD[8173] = 5.43667692E+00;
    COFD[8174] = -4.65640530E-01;
    COFD[8175] = 1.91675299E-02;
    COFD[8176] = -2.17161314E+01;
    COFD[8177] = 5.45061574E+00;
    COFD[8178] = -4.65241893E-01;
    COFD[8179] = 1.90650429E-02;
    COFD[8180] = -2.17280687E+01;
    COFD[8181] = 5.45061574E+00;
    COFD[8182] = -4.65241893E-01;
    COFD[8183] = 1.90650429E-02;
    COFD[8184] = -2.20886975E+01;
    COFD[8185] = 5.36809191E+00;
    COFD[8186] = -4.32908721E-01;
    COFD[8187] = 1.68481857E-02;
    COFD[8188] = -2.23311644E+01;
    COFD[8189] = 5.44415012E+00;
    COFD[8190] = -4.47195848E-01;
    COFD[8191] = 1.76333776E-02;
    COFD[8192] = -1.99882118E+01;
    COFD[8193] = 4.99038775E+00;
    COFD[8194] = -4.19232899E-01;
    COFD[8195] = 1.76017333E-02;
    COFD[8196] = -2.23423075E+01;
    COFD[8197] = 5.44415012E+00;
    COFD[8198] = -4.47195848E-01;
    COFD[8199] = 1.76333776E-02;
    COFD[8200] = -2.22481014E+01;
    COFD[8201] = 5.49056182E+00;
    COFD[8202] = -4.63917767E-01;
    COFD[8203] = 1.87597656E-02;
    COFD[8204] = -2.22481014E+01;
    COFD[8205] = 5.49056182E+00;
    COFD[8206] = -4.63917767E-01;
    COFD[8207] = 1.87597656E-02;
    COFD[8208] = -2.21401790E+01;
    COFD[8209] = 5.47453801E+00;
    COFD[8210] = -4.63573634E-01;
    COFD[8211] = 1.88106950E-02;
    COFD[8212] = -2.24977955E+01;
    COFD[8213] = 5.48975794E+00;
    COFD[8214] = -4.60336840E-01;
    COFD[8215] = 1.84719899E-02;
    COFD[8216] = -2.24826500E+01;
    COFD[8217] = 5.49046566E+00;
    COFD[8218] = -4.60640519E-01;
    COFD[8219] = 1.84931917E-02;
    COFD[8220] = -2.24863899E+01;
    COFD[8221] = 5.49046566E+00;
    COFD[8222] = -4.60640519E-01;
    COFD[8223] = 1.84931917E-02;
    COFD[8224] = -2.26120860E+01;
    COFD[8225] = 5.49130890E+00;
    COFD[8226] = -4.61073466E-01;
    COFD[8227] = 1.85242949E-02;
    COFD[8228] = -2.26033915E+01;
    COFD[8229] = 5.44450212E+00;
    COFD[8230] = -4.47289567E-01;
    COFD[8231] = 1.76392556E-02;
    COFD[8232] = -2.25626324E+01;
    COFD[8233] = 5.44454618E+00;
    COFD[8234] = -4.47321672E-01;
    COFD[8235] = 1.76416355E-02;
    COFD[8236] = -2.27592733E+01;
    COFD[8237] = 5.47456971E+00;
    COFD[8238] = -4.54548739E-01;
    COFD[8239] = 1.80779412E-02;
    COFD[8240] = -2.27565017E+01;
    COFD[8241] = 5.47937827E+00;
    COFD[8242] = -4.56129749E-01;
    COFD[8243] = 1.81820758E-02;
    COFD[8244] = -2.26777281E+01;
    COFD[8245] = 5.43682302E+00;
    COFD[8246] = -4.45758386E-01;
    COFD[8247] = 1.75526921E-02;
    COFD[8248] = -2.21849394E+01;
    COFD[8249] = 5.36809191E+00;
    COFD[8250] = -4.32908721E-01;
    COFD[8251] = 1.68481857E-02;
    COFD[8252] = -2.21897214E+01;
    COFD[8253] = 5.36809191E+00;
    COFD[8254] = -4.32908721E-01;
    COFD[8255] = 1.68481857E-02;
    COFD[8256] = -2.24172422E+01;
    COFD[8257] = 5.43317403E+00;
    COFD[8258] = -4.45039640E-01;
    COFD[8259] = 1.75123239E-02;
    COFD[8260] = -2.26170162E+01;
    COFD[8261] = 5.44849005E+00;
    COFD[8262] = -4.48374875E-01;
    COFD[8263] = 1.77073338E-02;
    COFD[8264] = -2.21977082E+01;
    COFD[8265] = 5.48376307E+00;
    COFD[8266] = -4.63942837E-01;
    COFD[8267] = 1.87953234E-02;
    COFD[8268] = -2.22417720E+01;
    COFD[8269] = 5.49056182E+00;
    COFD[8270] = -4.63917767E-01;
    COFD[8271] = 1.87597656E-02;
    COFD[8272] = -2.26331925E+01;
    COFD[8273] = 5.39941661E+00;
    COFD[8274] = -4.38608992E-01;
    COFD[8275] = 1.71564833E-02;
    COFD[8276] = -2.27356661E+01;
    COFD[8277] = 5.41656161E+00;
    COFD[8278] = -4.41824209E-01;
    COFD[8279] = 1.73331395E-02;
    COFD[8280] = -2.24209867E+01;
    COFD[8281] = 5.43317403E+00;
    COFD[8282] = -4.45039640E-01;
    COFD[8283] = 1.75123239E-02;
    COFD[8284] = -2.27859261E+01;
    COFD[8285] = 5.37750827E+00;
    COFD[8286] = -4.34574933E-01;
    COFD[8287] = 1.69370116E-02;
    COFD[8288] = -2.26722436E+01;
    COFD[8289] = 5.35408130E+00;
    COFD[8290] = -4.30433039E-01;
    COFD[8291] = 1.67163724E-02;
    COFD[8292] = -2.29454887E+01;
    COFD[8293] = 5.40726595E+00;
    COFD[8294] = -4.40072331E-01;
    COFD[8295] = 1.72366684E-02;
    COFD[8296] = -2.29440659E+01;
    COFD[8297] = 5.40726595E+00;
    COFD[8298] = -4.40072331E-01;
    COFD[8299] = 1.72366684E-02;
    COFD[8300] = -2.26749914E+01;
    COFD[8301] = 5.19870800E+00;
    COFD[8302] = -4.04106072E-01;
    COFD[8303] = 1.53462074E-02;
    COFD[8304] = -2.27261124E+01;
    COFD[8305] = 5.11476665E+00;
    COFD[8306] = -3.90017950E-01;
    COFD[8307] = 1.46174143E-02;
    COFD[8308] = -2.28339934E+01;
    COFD[8309] = 5.30874092E+00;
    COFD[8310] = -4.22570517E-01;
    COFD[8311] = 1.63018518E-02;
    COFD[8312] = -2.26545108E+01;
    COFD[8313] = 5.15653775E+00;
    COFD[8314] = -3.96988537E-01;
    COFD[8315] = 1.49769158E-02;
    COFD[8316] = -1.86427314E+01;
    COFD[8317] = 4.45508831E+00;
    COFD[8318] = -3.56068927E-01;
    COFD[8319] = 1.51117769E-02;
    COFD[8320] = -1.95514332E+01;
    COFD[8321] = 5.30461369E+00;
    COFD[8322] = -4.53088496E-01;
    COFD[8323] = 1.87998409E-02;
    COFD[8324] = -1.59135880E+01;
    COFD[8325] = 3.85500637E+00;
    COFD[8326] = -2.82834219E-01;
    COFD[8327] = 1.21237924E-02;
    COFD[8328] = -1.90756172E+01;
    COFD[8329] = 4.65976652E+00;
    COFD[8330] = -3.80698349E-01;
    COFD[8331] = 1.61038060E-02;
    COFD[8332] = -2.07843393E+01;
    COFD[8333] = 5.08900550E+00;
    COFD[8334] = -4.30341928E-01;
    COFD[8335] = 1.80173213E-02;
    COFD[8336] = -1.91004773E+01;
    COFD[8337] = 4.65976652E+00;
    COFD[8338] = -3.80698349E-01;
    COFD[8339] = 1.61038060E-02;
    COFD[8340] = -1.85864727E+01;
    COFD[8341] = 3.68972017E+00;
    COFD[8342] = -1.75261310E-01;
    COFD[8343] = 4.25906382E-03;
    COFD[8344] = -2.05043391E+01;
    COFD[8345] = 4.99719781E+00;
    COFD[8346] = -4.20035873E-01;
    COFD[8347] = 1.76333794E-02;
    COFD[8348] = -2.25606430E+01;
    COFD[8349] = 5.48480058E+00;
    COFD[8350] = -4.58431586E-01;
    COFD[8351] = 1.83416243E-02;
    COFD[8352] = -2.12336858E+01;
    COFD[8353] = 5.29671000E+00;
    COFD[8354] = -4.52208613E-01;
    COFD[8355] = 1.87672717E-02;
    COFD[8356] = -2.14398976E+01;
    COFD[8357] = 5.35638301E+00;
    COFD[8358] = -4.58721935E-01;
    COFD[8359] = 1.90025575E-02;
    COFD[8360] = -2.05805117E+01;
    COFD[8361] = 5.02104665E+00;
    COFD[8362] = -4.22831627E-01;
    COFD[8363] = 1.77428177E-02;
    COFD[8364] = -2.05907805E+01;
    COFD[8365] = 5.02104665E+00;
    COFD[8366] = -4.22831627E-01;
    COFD[8367] = 1.77428177E-02;
    COFD[8368] = -2.13808065E+01;
    COFD[8369] = 4.81362806E+00;
    COFD[8370] = -3.42326990E-01;
    COFD[8371] = 1.22325012E-02;
    COFD[8372] = -2.25355062E+01;
    COFD[8373] = 5.49116738E+00;
    COFD[8374] = -4.60872204E-01;
    COFD[8375] = 1.85085987E-02;
    COFD[8376] = -2.22850608E+01;
    COFD[8377] = 5.49115033E+00;
    COFD[8378] = -4.62266752E-01;
    COFD[8379] = 1.86223135E-02;
    COFD[8380] = -2.25234155E+01;
    COFD[8381] = 5.49116738E+00;
    COFD[8382] = -4.60872204E-01;
    COFD[8383] = 1.85085987E-02;
    COFD[8384] = -2.23705009E+01;
    COFD[8385] = 5.48227325E+00;
    COFD[8386] = -4.57157805E-01;
    COFD[8387] = 1.82509581E-02;
    COFD[8388] = -2.23850795E+01;
    COFD[8389] = 5.48227325E+00;
    COFD[8390] = -4.57157805E-01;
    COFD[8391] = 1.82509581E-02;
    COFD[8392] = -2.15333456E+01;
    COFD[8393] = 4.87709635E+00;
    COFD[8394] = -3.52250416E-01;
    COFD[8395] = 1.27244832E-02;
    COFD[8396] = -2.20442423E+01;
    COFD[8397] = 5.05908558E+00;
    COFD[8398] = -3.80945929E-01;
    COFD[8399] = 1.41553816E-02;
    COFD[8400] = -2.14253720E+01;
    COFD[8401] = 5.34211391E+00;
    COFD[8402] = -4.57230216E-01;
    COFD[8403] = 1.89516739E-02;
    COFD[8404] = -2.20589094E+01;
    COFD[8405] = 5.05908558E+00;
    COFD[8406] = -3.80945929E-01;
    COFD[8407] = 1.41553816E-02;
    COFD[8408] = -2.26414192E+01;
    COFD[8409] = 5.39941661E+00;
    COFD[8410] = -4.38608992E-01;
    COFD[8411] = 1.71564833E-02;
    COFD[8412] = -2.26414192E+01;
    COFD[8413] = 5.39941661E+00;
    COFD[8414] = -4.38608992E-01;
    COFD[8415] = 1.71564833E-02;
    COFD[8416] = -2.26656749E+01;
    COFD[8417] = 5.43804144E+00;
    COFD[8418] = -4.45999485E-01;
    COFD[8419] = 1.75662609E-02;
    COFD[8420] = -2.26989068E+01;
    COFD[8421] = 5.31130174E+00;
    COFD[8422] = -4.23012509E-01;
    COFD[8423] = 1.63250839E-02;
    COFD[8424] = -2.26989595E+01;
    COFD[8425] = 5.31729293E+00;
    COFD[8426] = -4.24040190E-01;
    COFD[8427] = 1.63789533E-02;
    COFD[8428] = -2.27040773E+01;
    COFD[8429] = 5.31729293E+00;
    COFD[8430] = -4.24040190E-01;
    COFD[8431] = 1.63789533E-02;
    COFD[8432] = -2.28412555E+01;
    COFD[8433] = 5.32619107E+00;
    COFD[8434] = -4.25570588E-01;
    COFD[8435] = 1.64592712E-02;
    COFD[8436] = -2.23275170E+01;
    COFD[8437] = 5.06089774E+00;
    COFD[8438] = -3.81239105E-01;
    COFD[8439] = 1.41702517E-02;
    COFD[8440] = -2.10466968E+01;
    COFD[8441] = 4.48381723E+00;
    COFD[8442] = -2.96509945E-01;
    COFD[8443] = 1.01535197E-02;
    COFD[8444] = -2.27023831E+01;
    COFD[8445] = 5.18383798E+00;
    COFD[8446] = -4.01666785E-01;
    COFD[8447] = 1.52216204E-02;
    COFD[8448] = -2.26781027E+01;
    COFD[8449] = 5.18129009E+00;
    COFD[8450] = -4.02449417E-01;
    COFD[8451] = 1.52976330E-02;
    COFD[8452] = -2.23869377E+01;
    COFD[8453] = 5.04130371E+00;
    COFD[8454] = -3.78082012E-01;
    COFD[8455] = 1.40104739E-02;
    COFD[8456] = -2.16580095E+01;
    COFD[8457] = 4.87709635E+00;
    COFD[8458] = -3.52250416E-01;
    COFD[8459] = 1.27244832E-02;
    COFD[8460] = -2.16643993E+01;
    COFD[8461] = 4.87709635E+00;
    COFD[8462] = -3.52250416E-01;
    COFD[8463] = 1.27244832E-02;
    COFD[8464] = -2.21073752E+01;
    COFD[8465] = 5.03131423E+00;
    COFD[8466] = -3.76503291E-01;
    COFD[8467] = 1.39316539E-02;
    COFD[8468] = -2.11645228E+01;
    COFD[8469] = 4.51897995E+00;
    COFD[8470] = -3.01998503E-01;
    COFD[8471] = 1.04247720E-02;
    COFD[8472] = -2.26581212E+01;
    COFD[8473] = 5.42030479E+00;
    COFD[8474] = -4.42541204E-01;
    COFD[8475] = 1.73729052E-02;
    COFD[8476] = -2.26331925E+01;
    COFD[8477] = 5.39941661E+00;
    COFD[8478] = -4.38608992E-01;
    COFD[8479] = 1.71564833E-02;
    COFD[8480] = -2.11409832E+01;
    COFD[8481] = 4.46557987E+00;
    COFD[8482] = -2.93692029E-01;
    COFD[8483] = 1.00124437E-02;
    COFD[8484] = -2.11452590E+01;
    COFD[8485] = 4.43017828E+00;
    COFD[8486] = -2.88253942E-01;
    COFD[8487] = 9.74852180E-03;
    COFD[8488] = -2.21124987E+01;
    COFD[8489] = 5.03131423E+00;
    COFD[8490] = -3.76503291E-01;
    COFD[8491] = 1.39316539E-02;
    COFD[8492] = -2.15053166E+01;
    COFD[8493] = 4.55066138E+00;
    COFD[8494] = -3.06272562E-01;
    COFD[8495] = 1.06092478E-02;
    COFD[8496] = -2.21244697E+01;
    COFD[8497] = 4.84823431E+00;
    COFD[8498] = -3.47731563E-01;
    COFD[8499] = 1.25001760E-02;
    COFD[8500] = -2.25506266E+01;
    COFD[8501] = 4.96541553E+00;
    COFD[8502] = -3.66120907E-01;
    COFD[8503] = 1.34143879E-02;
    COFD[8504] = -2.25485050E+01;
    COFD[8505] = 4.96541553E+00;
    COFD[8506] = -3.66120907E-01;
    COFD[8507] = 1.34143879E-02;
    COFD[8508] = -2.13962129E+01;
    COFD[8509] = 4.37956579E+00;
    COFD[8510] = -2.78446618E-01;
    COFD[8511] = 9.19232241E-03;
    COFD[8512] = -2.13456084E+01;
    COFD[8513] = 4.25479441E+00;
    COFD[8514] = -2.59093497E-01;
    COFD[8515] = 8.23869170E-03;
    COFD[8516] = -2.21833064E+01;
    COFD[8517] = 4.75657673E+00;
    COFD[8518] = -3.33531354E-01;
    COFD[8519] = 1.18003758E-02;
    COFD[8520] = -2.12434067E+01;
    COFD[8521] = 4.28015384E+00;
    COFD[8522] = -2.63514048E-01;
    COFD[8523] = 8.47277117E-03;
    COFD[8524] = -2.04535647E+01;
    COFD[8525] = 4.98209169E+00;
    COFD[8526] = -4.18253583E-01;
    COFD[8527] = 1.75630850E-02;
    COFD[8528] = -1.95829229E+01;
    COFD[8529] = 5.27835392E+00;
    COFD[8530] = -4.50143769E-01;
    COFD[8531] = 1.86899020E-02;
    COFD[8532] = -1.60254653E+01;
    COFD[8533] = 3.85862177E+00;
    COFD[8534] = -2.83262563E-01;
    COFD[8535] = 1.21406839E-02;
    COFD[8536] = -1.90773916E+01;
    COFD[8537] = 4.62426562E+00;
    COFD[8538] = -3.76333900E-01;
    COFD[8539] = 1.59240937E-02;
    COFD[8540] = -2.08687893E+01;
    COFD[8541] = 5.08643773E+00;
    COFD[8542] = -4.30070504E-01;
    COFD[8543] = 1.80079562E-02;
    COFD[8544] = -1.91021863E+01;
    COFD[8545] = 4.62426562E+00;
    COFD[8546] = -3.76333900E-01;
    COFD[8547] = 1.59240937E-02;
    COFD[8548] = -1.86138117E+01;
    COFD[8549] = 3.64451102E+00;
    COFD[8550] = -1.67695039E-01;
    COFD[8551] = 3.87428023E-03;
    COFD[8552] = -2.05899148E+01;
    COFD[8553] = 4.99533607E+00;
    COFD[8554] = -4.19816362E-01;
    COFD[8555] = 1.76247292E-02;
    COFD[8556] = -2.26462274E+01;
    COFD[8557] = 5.48561354E+00;
    COFD[8558] = -4.58748903E-01;
    COFD[8559] = 1.83633678E-02;
    COFD[8560] = -2.12409720E+01;
    COFD[8561] = 5.26987677E+00;
    COFD[8562] = -4.49189697E-01;
    COFD[8563] = 1.86541737E-02;
    COFD[8564] = -2.15282219E+01;
    COFD[8565] = 5.35564147E+00;
    COFD[8566] = -4.58644607E-01;
    COFD[8567] = 1.89999292E-02;
    COFD[8568] = -2.05856649E+01;
    COFD[8569] = 4.99149256E+00;
    COFD[8570] = -4.19363176E-01;
    COFD[8571] = 1.76068685E-02;
    COFD[8572] = -2.05958878E+01;
    COFD[8573] = 4.99149256E+00;
    COFD[8574] = -4.19363176E-01;
    COFD[8575] = 1.76068685E-02;
    COFD[8576] = -2.15521873E+01;
    COFD[8577] = 4.85623963E+00;
    COFD[8578] = -3.48984100E-01;
    COFD[8579] = 1.25623211E-02;
    COFD[8580] = -2.25966091E+01;
    COFD[8581] = 5.49074590E+00;
    COFD[8582] = -4.61641528E-01;
    COFD[8583] = 1.85733572E-02;
    COFD[8584] = -2.23629227E+01;
    COFD[8585] = 5.49413323E+00;
    COFD[8586] = -4.63498940E-01;
    COFD[8587] = 1.87080519E-02;
    COFD[8588] = -2.25845680E+01;
    COFD[8589] = 5.49074590E+00;
    COFD[8590] = -4.61641528E-01;
    COFD[8591] = 1.85733572E-02;
    COFD[8592] = -2.24444141E+01;
    COFD[8593] = 5.48421194E+00;
    COFD[8594] = -4.58197253E-01;
    COFD[8595] = 1.83255177E-02;
    COFD[8596] = -2.24589388E+01;
    COFD[8597] = 5.48421194E+00;
    COFD[8598] = -4.58197253E-01;
    COFD[8599] = 1.83255177E-02;
    COFD[8600] = -2.16990357E+01;
    COFD[8601] = 4.91762397E+00;
    COFD[8602] = -3.58634830E-01;
    COFD[8603] = 1.30425594E-02;
    COFD[8604] = -2.21907779E+01;
    COFD[8605] = 5.09350157E+00;
    COFD[8606] = -3.86533585E-01;
    COFD[8607] = 1.44394093E-02;
    COFD[8608] = -2.14529319E+01;
    COFD[8609] = 5.31795346E+00;
    COFD[8610] = -4.54574565E-01;
    COFD[8611] = 1.88549001E-02;
    COFD[8612] = -2.22053677E+01;
    COFD[8613] = 5.09350157E+00;
    COFD[8614] = -3.86533585E-01;
    COFD[8615] = 1.44394093E-02;
    COFD[8616] = -2.27438518E+01;
    COFD[8617] = 5.41656161E+00;
    COFD[8618] = -4.41824209E-01;
    COFD[8619] = 1.73331395E-02;
    COFD[8620] = -2.27438518E+01;
    COFD[8621] = 5.41656161E+00;
    COFD[8622] = -4.41824209E-01;
    COFD[8623] = 1.73331395E-02;
    COFD[8624] = -2.27506175E+01;
    COFD[8625] = 5.44709568E+00;
    COFD[8626] = -4.48069354E-01;
    COFD[8627] = 1.76894758E-02;
    COFD[8628] = -2.28127150E+01;
    COFD[8629] = 5.33452791E+00;
    COFD[8630] = -4.27009745E-01;
    COFD[8631] = 1.65349406E-02;
    COFD[8632] = -2.28126636E+01;
    COFD[8633] = 5.34025709E+00;
    COFD[8634] = -4.28002604E-01;
    COFD[8635] = 1.65872781E-02;
    COFD[8636] = -2.28177500E+01;
    COFD[8637] = 5.34025709E+00;
    COFD[8638] = -4.28002604E-01;
    COFD[8639] = 1.65872781E-02;
    COFD[8640] = -2.29487149E+01;
    COFD[8641] = 5.34854008E+00;
    COFD[8642] = -4.29453170E-01;
    COFD[8643] = 1.66641874E-02;
    COFD[8644] = -2.24655743E+01;
    COFD[8645] = 5.09524210E+00;
    COFD[8646] = -3.86817295E-01;
    COFD[8647] = 1.44538658E-02;
    COFD[8648] = -2.10487877E+01;
    COFD[8649] = 4.43820990E+00;
    COFD[8650] = -2.88877613E-01;
    COFD[8651] = 9.76344368E-03;
    COFD[8652] = -2.28320950E+01;
    COFD[8653] = 5.21637493E+00;
    COFD[8654] = -4.06956958E-01;
    COFD[8655] = 1.54905633E-02;
    COFD[8656] = -2.27790113E+01;
    COFD[8657] = 5.20168946E+00;
    COFD[8658] = -4.06049775E-01;
    COFD[8659] = 1.54891495E-02;
    COFD[8660] = -2.25262686E+01;
    COFD[8661] = 5.07694371E+00;
    COFD[8662] = -3.83840072E-01;
    COFD[8663] = 1.43023310E-02;
    COFD[8664] = -2.18230852E+01;
    COFD[8665] = 4.91762397E+00;
    COFD[8666] = -3.58634830E-01;
    COFD[8667] = 1.30425594E-02;
    COFD[8668] = -2.18294392E+01;
    COFD[8669] = 4.91762397E+00;
    COFD[8670] = -3.58634830E-01;
    COFD[8671] = 1.30425594E-02;
    COFD[8672] = -2.22582230E+01;
    COFD[8673] = 5.06826700E+00;
    COFD[8674] = -3.82432503E-01;
    COFD[8675] = 1.42308185E-02;
    COFD[8676] = -2.11597907E+01;
    COFD[8677] = 4.47157782E+00;
    COFD[8678] = -2.94123466E-01;
    COFD[8679] = 1.00236418E-02;
    COFD[8680] = -2.27561195E+01;
    COFD[8681] = 5.43540157E+00;
    COFD[8682] = -4.45477820E-01;
    COFD[8683] = 1.75369197E-02;
    COFD[8684] = -2.27356661E+01;
    COFD[8685] = 5.41656161E+00;
    COFD[8686] = -4.41824209E-01;
    COFD[8687] = 1.73331395E-02;
    COFD[8688] = -2.11452590E+01;
    COFD[8689] = 4.43017828E+00;
    COFD[8690] = -2.88253942E-01;
    COFD[8691] = 9.74852180E-03;
    COFD[8692] = -2.11277860E+01;
    COFD[8693] = 4.38003380E+00;
    COFD[8694] = -2.80097327E-01;
    COFD[8695] = 9.33767565E-03;
    COFD[8696] = -2.22633150E+01;
    COFD[8697] = 5.06826700E+00;
    COFD[8698] = -3.82432503E-01;
    COFD[8699] = 1.42308185E-02;
    COFD[8700] = -2.15153764E+01;
    COFD[8701] = 4.52500980E+00;
    COFD[8702] = -3.02670241E-01;
    COFD[8703] = 1.04452726E-02;
    COFD[8704] = -2.22786912E+01;
    COFD[8705] = 4.89093124E+00;
    COFD[8706] = -3.54419517E-01;
    COFD[8707] = 1.28322610E-02;
    COFD[8708] = -2.26945043E+01;
    COFD[8709] = 5.00634977E+00;
    COFD[8710] = -3.72557302E-01;
    COFD[8711] = 1.37346802E-02;
    COFD[8712] = -2.26924001E+01;
    COFD[8713] = 5.00634977E+00;
    COFD[8714] = -3.72557302E-01;
    COFD[8715] = 1.37346802E-02;
    COFD[8716] = -2.14928955E+01;
    COFD[8717] = 4.39966514E+00;
    COFD[8718] = -2.81886234E-01;
    COFD[8719] = 9.37271319E-03;
    COFD[8720] = -2.14606335E+01;
    COFD[8721] = 4.28507327E+00;
    COFD[8722] = -2.63941348E-01;
    COFD[8723] = 8.48270771E-03;
    COFD[8724] = -2.23442480E+01;
    COFD[8725] = 4.80471822E+00;
    COFD[8726] = -3.40939669E-01;
    COFD[8727] = 1.21639062E-02;
    COFD[8728] = -2.13309187E+01;
    COFD[8729] = 4.29619647E+00;
    COFD[8730] = -2.66307880E-01;
    COFD[8731] = 8.62039924E-03;
    COFD[8732] = -2.05335806E+01;
    COFD[8733] = 4.97826004E+00;
    COFD[8734] = -4.17803418E-01;
    COFD[8735] = 1.75454151E-02;
    COFD[8736] = -1.91107966E+01;
    COFD[8737] = 5.24865128E+00;
    COFD[8738] = -4.46802764E-01;
    COFD[8739] = 1.85649476E-02;
    COFD[8740] = -1.52496383E+01;
    COFD[8741] = 3.70713307E+00;
    COFD[8742] = -2.64771252E-01;
    COFD[8743] = 1.13878078E-02;
    COFD[8744] = -1.85982326E+01;
    COFD[8745] = 4.58430295E+00;
    COFD[8746] = -3.71395187E-01;
    COFD[8747] = 1.57196293E-02;
    COFD[8748] = -2.01255512E+01;
    COFD[8749] = 4.95889869E+00;
    COFD[8750] = -4.15526934E-01;
    COFD[8751] = 1.74559870E-02;
    COFD[8752] = -1.86218352E+01;
    COFD[8753] = 4.58430295E+00;
    COFD[8754] = -3.71395187E-01;
    COFD[8755] = 1.57196293E-02;
    COFD[8756] = -2.02550404E+01;
    COFD[8757] = 4.60540477E+00;
    COFD[8758] = -3.10373488E-01;
    COFD[8759] = 1.06672049E-02;
    COFD[8760] = -1.98276037E+01;
    COFD[8761] = 4.85730922E+00;
    COFD[8762] = -4.04004009E-01;
    COFD[8763] = 1.70221852E-02;
    COFD[8764] = -2.22190198E+01;
    COFD[8765] = 5.49406269E+00;
    COFD[8766] = -4.63519837E-01;
    COFD[8767] = 1.87101015E-02;
    COFD[8768] = -2.08252693E+01;
    COFD[8769] = 5.23993213E+00;
    COFD[8770] = -4.45818789E-01;
    COFD[8771] = 1.85280084E-02;
    COFD[8772] = -2.07647491E+01;
    COFD[8773] = 5.21569495E+00;
    COFD[8774] = -4.43054270E-01;
    COFD[8775] = 1.84229515E-02;
    COFD[8776] = -2.01353679E+01;
    COFD[8777] = 4.95889869E+00;
    COFD[8778] = -4.15526934E-01;
    COFD[8779] = 1.74559870E-02;
    COFD[8780] = -2.01447821E+01;
    COFD[8781] = 4.95889869E+00;
    COFD[8782] = -4.15526934E-01;
    COFD[8783] = 1.74559870E-02;
    COFD[8784] = -2.12901259E+01;
    COFD[8785] = 4.90168525E+00;
    COFD[8786] = -3.56120032E-01;
    COFD[8787] = 1.29171867E-02;
    COFD[8788] = -2.22582078E+01;
    COFD[8789] = 5.49328085E+00;
    COFD[8790] = -4.62887295E-01;
    COFD[8791] = 1.86623580E-02;
    COFD[8792] = -2.19838573E+01;
    COFD[8793] = 5.49090680E+00;
    COFD[8794] = -4.63910285E-01;
    COFD[8795] = 1.87574539E-02;
    COFD[8796] = -2.22470457E+01;
    COFD[8797] = 5.49328085E+00;
    COFD[8798] = -4.62887295E-01;
    COFD[8799] = 1.86623580E-02;
    COFD[8800] = -2.20924770E+01;
    COFD[8801] = 5.48779691E+00;
    COFD[8802] = -4.59557311E-01;
    COFD[8803] = 1.84182884E-02;
    COFD[8804] = -2.21060408E+01;
    COFD[8805] = 5.48779691E+00;
    COFD[8806] = -4.59557311E-01;
    COFD[8807] = 1.84182884E-02;
    COFD[8808] = -2.14384166E+01;
    COFD[8809] = 4.96375446E+00;
    COFD[8810] = -3.65860553E-01;
    COFD[8811] = 1.34014562E-02;
    COFD[8812] = -2.19069859E+01;
    COFD[8813] = 5.12829880E+00;
    COFD[8814] = -3.92242339E-01;
    COFD[8815] = 1.47312337E-02;
    COFD[8816] = -2.09647485E+01;
    COFD[8817] = 5.29024160E+00;
    COFD[8818] = -4.51483508E-01;
    COFD[8819] = 1.87402157E-02;
    COFD[8820] = -2.19202360E+01;
    COFD[8821] = 5.12829880E+00;
    COFD[8822] = -3.92242339E-01;
    COFD[8823] = 1.47312337E-02;
    COFD[8824] = -2.24284563E+01;
    COFD[8825] = 5.43317403E+00;
    COFD[8826] = -4.45039640E-01;
    COFD[8827] = 1.75123239E-02;
    COFD[8828] = -2.24284563E+01;
    COFD[8829] = 5.43317403E+00;
    COFD[8830] = -4.45039640E-01;
    COFD[8831] = 1.75123239E-02;
    COFD[8832] = -2.24141963E+01;
    COFD[8833] = 5.45651702E+00;
    COFD[8834] = -4.50175370E-01;
    COFD[8835] = 1.78134258E-02;
    COFD[8836] = -2.25186794E+01;
    COFD[8837] = 5.35783672E+00;
    COFD[8838] = -4.31095103E-01;
    COFD[8839] = 1.67515802E-02;
    COFD[8840] = -2.25148268E+01;
    COFD[8841] = 5.36306078E+00;
    COFD[8842] = -4.32018705E-01;
    COFD[8843] = 1.68007664E-02;
    COFD[8844] = -2.25193776E+01;
    COFD[8845] = 5.36306078E+00;
    COFD[8846] = -4.32018705E-01;
    COFD[8847] = 1.68007664E-02;
    COFD[8848] = -2.26652914E+01;
    COFD[8849] = 5.37102704E+00;
    COFD[8850] = -4.33426866E-01;
    COFD[8851] = 1.68757765E-02;
    COFD[8852] = -2.22001994E+01;
    COFD[8853] = 5.12996650E+00;
    COFD[8854] = -3.92516016E-01;
    COFD[8855] = 1.47452270E-02;
    COFD[8856] = -2.21616557E+01;
    COFD[8857] = 5.13096622E+00;
    COFD[8858] = -3.92680132E-01;
    COFD[8859] = 1.47536199E-02;
    COFD[8860] = -2.25693752E+01;
    COFD[8861] = 5.24758220E+00;
    COFD[8862] = -4.12129116E-01;
    COFD[8863] = 1.57564573E-02;
    COFD[8864] = -2.26296673E+01;
    COFD[8865] = 5.27897272E+00;
    COFD[8866] = -4.17421132E-01;
    COFD[8867] = 1.60310357E-02;
    COFD[8868] = -2.22622113E+01;
    COFD[8869] = 5.11286874E+00;
    COFD[8870] = -3.89706255E-01;
    COFD[8871] = 1.46014745E-02;
    COFD[8872] = -2.15517540E+01;
    COFD[8873] = 4.96375446E+00;
    COFD[8874] = -3.65860553E-01;
    COFD[8875] = 1.34014562E-02;
    COFD[8876] = -2.15574915E+01;
    COFD[8877] = 4.96375446E+00;
    COFD[8878] = -3.65860553E-01;
    COFD[8879] = 1.34014562E-02;
    COFD[8880] = -2.19759929E+01;
    COFD[8881] = 5.10490214E+00;
    COFD[8882] = -3.88397529E-01;
    COFD[8883] = 1.45345306E-02;
    COFD[8884] = -2.22483379E+01;
    COFD[8885] = 5.14893259E+00;
    COFD[8886] = -3.95700512E-01;
    COFD[8887] = 1.49099774E-02;
    COFD[8888] = -2.24242197E+01;
    COFD[8889] = 5.44615881E+00;
    COFD[8890] = -4.47865005E-01;
    COFD[8891] = 1.76775502E-02;
    COFD[8892] = -2.24209867E+01;
    COFD[8893] = 5.43317403E+00;
    COFD[8894] = -4.45039640E-01;
    COFD[8895] = 1.75123239E-02;
    COFD[8896] = -2.21124987E+01;
    COFD[8897] = 5.03131423E+00;
    COFD[8898] = -3.76503291E-01;
    COFD[8899] = 1.39316539E-02;
    COFD[8900] = -2.22633150E+01;
    COFD[8901] = 5.06826700E+00;
    COFD[8902] = -3.82432503E-01;
    COFD[8903] = 1.42308185E-02;
    COFD[8904] = -2.19805490E+01;
    COFD[8905] = 5.10490214E+00;
    COFD[8906] = -3.88397529E-01;
    COFD[8907] = 1.45345306E-02;
    COFD[8908] = -2.22109105E+01;
    COFD[8909] = 4.98256919E+00;
    COFD[8910] = -3.68811981E-01;
    COFD[8911] = 1.35481236E-02;
    COFD[8912] = -2.20373227E+01;
    COFD[8913] = 4.93482700E+00;
    COFD[8914] = -3.61331756E-01;
    COFD[8915] = 1.31765767E-02;
    COFD[8916] = -2.24658579E+01;
    COFD[8917] = 5.04880598E+00;
    COFD[8918] = -3.79285322E-01;
    COFD[8919] = 1.40712210E-02;
    COFD[8920] = -2.24640385E+01;
    COFD[8921] = 5.04880598E+00;
    COFD[8922] = -3.79285322E-01;
    COFD[8923] = 1.40712210E-02;
    COFD[8924] = -2.17452723E+01;
    COFD[8925] = 4.64756828E+00;
    COFD[8926] = -3.16905849E-01;
    COFD[8927] = 1.09890709E-02;
    COFD[8928] = -2.16496748E+01;
    COFD[8929] = 4.49834811E+00;
    COFD[8930] = -2.93898563E-01;
    COFD[8931] = 9.85901966E-03;
    COFD[8932] = -2.21227489E+01;
    COFD[8933] = 4.85109413E+00;
    COFD[8934] = -3.48178913E-01;
    COFD[8935] = 1.25223681E-02;
    COFD[8936] = -2.16471363E+01;
    COFD[8937] = 4.57106634E+00;
    COFD[8938] = -3.05074221E-01;
    COFD[8939] = 1.04066420E-02;
    COFD[8940] = -1.98008145E+01;
    COFD[8941] = 4.85079945E+00;
    COFD[8942] = -4.03257284E-01;
    COFD[8943] = 1.69936877E-02;
    COFD[8944] = -1.98427244E+01;
    COFD[8945] = 5.33546471E+00;
    COFD[8946] = -4.56520851E-01;
    COFD[8947] = 1.89268199E-02;
    COFD[8948] = -1.60659819E+01;
    COFD[8949] = 3.84737458E+00;
    COFD[8950] = -2.81927023E-01;
    COFD[8951] = 1.20878691E-02;
    COFD[8952] = -1.93826457E+01;
    COFD[8953] = 4.70007538E+00;
    COFD[8954] = -3.85598099E-01;
    COFD[8955] = 1.63031809E-02;
    COFD[8956] = -2.09753274E+01;
    COFD[8957] = 5.08993111E+00;
    COFD[8958] = -4.30439721E-01;
    COFD[8959] = 1.80206932E-02;
    COFD[8960] = -1.94082832E+01;
    COFD[8961] = 4.70007538E+00;
    COFD[8962] = -3.85598099E-01;
    COFD[8963] = 1.63031809E-02;
    COFD[8964] = -1.91771822E+01;
    COFD[8965] = 3.89325739E+00;
    COFD[8966] = -2.06939169E-01;
    COFD[8967] = 5.80832388E-03;
    COFD[8968] = -2.06819784E+01;
    COFD[8969] = 4.99569563E+00;
    COFD[8970] = -4.19858765E-01;
    COFD[8971] = 1.76264005E-02;
    COFD[8972] = -2.27524699E+01;
    COFD[8973] = 5.48417137E+00;
    COFD[8974] = -4.58172226E-01;
    COFD[8975] = 1.83236937E-02;
    COFD[8976] = -2.15010463E+01;
    COFD[8977] = 5.32769529E+00;
    COFD[8978] = -4.55658342E-01;
    COFD[8979] = 1.88949932E-02;
    COFD[8980] = -2.16063225E+01;
    COFD[8981] = 5.35467220E+00;
    COFD[8982] = -4.58543488E-01;
    COFD[8983] = 1.89964903E-02;
    COFD[8984] = -2.08689562E+01;
    COFD[8985] = 5.05404465E+00;
    COFD[8986] = -4.26628831E-01;
    COFD[8987] = 1.78883654E-02;
    COFD[8988] = -2.08797845E+01;
    COFD[8989] = 5.05404465E+00;
    COFD[8990] = -4.26628831E-01;
    COFD[8991] = 1.78883654E-02;
    COFD[8992] = -2.14531444E+01;
    COFD[8993] = 4.75554226E+00;
    COFD[8994] = -3.33372867E-01;
    COFD[8995] = 1.17926183E-02;
    COFD[8996] = -2.27245540E+01;
    COFD[8997] = 5.48717471E+00;
    COFD[8998] = -4.59324333E-01;
    COFD[8999] = 1.84024197E-02;
    COFD[9000] = -2.24972209E+01;
    COFD[9001] = 5.49117520E+00;
    COFD[9002] = -4.61254571E-01;
    COFD[9003] = 1.85397061E-02;
    COFD[9004] = -2.27118623E+01;
    COFD[9005] = 5.48717471E+00;
    COFD[9006] = -4.59324333E-01;
    COFD[9007] = 1.84024197E-02;
    COFD[9008] = -2.25641242E+01;
    COFD[9009] = 5.47765924E+00;
    COFD[9010] = -4.55592027E-01;
    COFD[9011] = 1.81470399E-02;
    COFD[9012] = -2.25793511E+01;
    COFD[9013] = 5.47765924E+00;
    COFD[9014] = -4.55592027E-01;
    COFD[9015] = 1.81470399E-02;
    COFD[9016] = -2.16206093E+01;
    COFD[9017] = 4.82487053E+00;
    COFD[9018] = -3.44081206E-01;
    COFD[9019] = 1.23193175E-02;
    COFD[9020] = -2.21441058E+01;
    COFD[9021] = 5.01223834E+00;
    COFD[9022] = -3.73486752E-01;
    COFD[9023] = 1.37810367E-02;
    COFD[9024] = -2.17208609E+01;
    COFD[9025] = 5.36702791E+00;
    COFD[9026] = -4.59828733E-01;
    COFD[9027] = 1.90400208E-02;
    COFD[9028] = -2.21597235E+01;
    COFD[9029] = 5.01223834E+00;
    COFD[9030] = -3.73486752E-01;
    COFD[9031] = 1.37810367E-02;
    COFD[9032] = -2.27946566E+01;
    COFD[9033] = 5.37750827E+00;
    COFD[9034] = -4.34574933E-01;
    COFD[9035] = 1.69370116E-02;
    COFD[9036] = -2.27946566E+01;
    COFD[9037] = 5.37750827E+00;
    COFD[9038] = -4.34574933E-01;
    COFD[9039] = 1.69370116E-02;
    COFD[9040] = -2.28313929E+01;
    COFD[9041] = 5.41966149E+00;
    COFD[9042] = -4.42417688E-01;
    COFD[9043] = 1.73660474E-02;
    COFD[9044] = -2.28386812E+01;
    COFD[9045] = 5.28428478E+00;
    COFD[9046] = -4.18334137E-01;
    COFD[9047] = 1.60789015E-02;
    COFD[9048] = -2.28395602E+01;
    COFD[9049] = 5.29000716E+00;
    COFD[9050] = -4.19320503E-01;
    COFD[9051] = 1.61306865E-02;
    COFD[9052] = -2.28450684E+01;
    COFD[9053] = 5.29000716E+00;
    COFD[9054] = -4.19320503E-01;
    COFD[9055] = 1.61306865E-02;
    COFD[9056] = -2.29733635E+01;
    COFD[9057] = 5.29835260E+00;
    COFD[9058] = -4.20767237E-01;
    COFD[9059] = 1.62068436E-02;
    COFD[9060] = -2.24208370E+01;
    COFD[9061] = 5.01427991E+00;
    COFD[9062] = -3.73809183E-01;
    COFD[9063] = 1.37971237E-02;
    COFD[9064] = -2.14543740E+01;
    COFD[9065] = 4.59523034E+00;
    COFD[9066] = -3.13534126E-01;
    COFD[9067] = 1.09798752E-02;
    COFD[9068] = -2.28314712E+01;
    COFD[9069] = 5.15368301E+00;
    COFD[9070] = -3.96504343E-01;
    COFD[9071] = 1.49517329E-02;
    COFD[9072] = -2.28089000E+01;
    COFD[9073] = 5.15221867E+00;
    COFD[9074] = -3.97461937E-01;
    COFD[9075] = 1.50362563E-02;
    COFD[9076] = -2.24776133E+01;
    COFD[9077] = 4.99278627E+00;
    COFD[9078] = -3.70419518E-01;
    COFD[9079] = 1.36281487E-02;
    COFD[9080] = -2.17527972E+01;
    COFD[9081] = 4.82487053E+00;
    COFD[9082] = -3.44081206E-01;
    COFD[9083] = 1.23193175E-02;
    COFD[9084] = -2.17596291E+01;
    COFD[9085] = 4.82487053E+00;
    COFD[9086] = -3.44081206E-01;
    COFD[9087] = 1.23193175E-02;
    COFD[9088] = -2.22053962E+01;
    COFD[9089] = 4.98256919E+00;
    COFD[9090] = -3.68811981E-01;
    COFD[9091] = 1.35481236E-02;
    COFD[9092] = -2.15608572E+01;
    COFD[9093] = 4.62655831E+00;
    COFD[9094] = -3.18380357E-01;
    COFD[9095] = 1.12180276E-02;
    COFD[9096] = -2.28163992E+01;
    COFD[9097] = 5.39961965E+00;
    COFD[9098] = -4.38647062E-01;
    COFD[9099] = 1.71585740E-02;
    COFD[9100] = -2.27859261E+01;
    COFD[9101] = 5.37750827E+00;
    COFD[9102] = -4.34574933E-01;
    COFD[9103] = 1.69370116E-02;
    COFD[9104] = -2.15053166E+01;
    COFD[9105] = 4.55066138E+00;
    COFD[9106] = -3.06272562E-01;
    COFD[9107] = 1.06092478E-02;
    COFD[9108] = -2.15153764E+01;
    COFD[9109] = 4.52500980E+00;
    COFD[9110] = -3.02670241E-01;
    COFD[9111] = 1.04452726E-02;
    COFD[9112] = -2.22109105E+01;
    COFD[9113] = 4.98256919E+00;
    COFD[9114] = -3.68811981E-01;
    COFD[9115] = 1.35481236E-02;
    COFD[9116] = -2.18081829E+01;
    COFD[9117] = 4.60265738E+00;
    COFD[9118] = -3.13370476E-01;
    COFD[9119] = 1.09264302E-02;
    COFD[9120] = -2.22044849E+01;
    COFD[9121] = 4.79539928E+00;
    COFD[9122] = -3.39494765E-01;
    COFD[9123] = 1.20926677E-02;
    COFD[9124] = -2.26263872E+01;
    COFD[9125] = 4.91358137E+00;
    COFD[9126] = -3.57998484E-01;
    COFD[9127] = 1.30108794E-02;
    COFD[9128] = -2.26240442E+01;
    COFD[9129] = 4.91358137E+00;
    COFD[9130] = -3.57998484E-01;
    COFD[9131] = 1.30108794E-02;
    COFD[9132] = -2.15478797E+01;
    COFD[9133] = 4.35915869E+00;
    COFD[9134] = -2.74825098E-01;
    COFD[9135] = 8.99908615E-03;
    COFD[9136] = -2.14478435E+01;
    COFD[9137] = 4.21514791E+00;
    COFD[9138] = -2.52791771E-01;
    COFD[9139] = 7.92322551E-03;
    COFD[9140] = -2.22319277E+01;
    COFD[9141] = 4.69258044E+00;
    COFD[9142] = -3.23762098E-01;
    COFD[9143] = 1.13233362E-02;
    COFD[9144] = -2.13912488E+01;
    COFD[9145] = 4.25856149E+00;
    COFD[9146] = -2.59776530E-01;
    COFD[9147] = 8.27569932E-03;
    COFD[9148] = -2.06399268E+01;
    COFD[9149] = 4.98336829E+00;
    COFD[9150] = -4.18404333E-01;
    COFD[9151] = 1.75690366E-02;
    COFD[9152] = -1.98161406E+01;
    COFD[9153] = 5.35897956E+00;
    COFD[9154] = -4.58992472E-01;
    COFD[9155] = 1.90117416E-02;
    COFD[9156] = -1.58862636E+01;
    COFD[9157] = 3.81826954E+00;
    COFD[9158] = -2.78442735E-01;
    COFD[9159] = 1.19489210E-02;
    COFD[9160] = -1.93847925E+01;
    COFD[9161] = 4.73396422E+00;
    COFD[9162] = -3.89665033E-01;
    COFD[9163] = 1.64665068E-02;
    COFD[9164] = -2.08403997E+01;
    COFD[9165] = 5.07617318E+00;
    COFD[9166] = -4.28983517E-01;
    COFD[9167] = 1.79703565E-02;
    COFD[9168] = -1.94103813E+01;
    COFD[9169] = 4.73396422E+00;
    COFD[9170] = -3.89665033E-01;
    COFD[9171] = 1.64665068E-02;
    COFD[9172] = -2.01687636E+01;
    COFD[9173] = 4.37329337E+00;
    COFD[9174] = -2.74949941E-01;
    COFD[9175] = 8.93930138E-03;
    COFD[9176] = -2.05171137E+01;
    COFD[9177] = 4.96961006E+00;
    COFD[9178] = -4.16785706E-01;
    COFD[9179] = 1.75054037E-02;
    COFD[9180] = -2.26588618E+01;
    COFD[9181] = 5.48550227E+00;
    COFD[9182] = -4.58704694E-01;
    COFD[9183] = 1.83603287E-02;
    COFD[9184] = -2.14961087E+01;
    COFD[9185] = 5.35264246E+00;
    COFD[9186] = -4.58331575E-01;
    COFD[9187] = 1.89892757E-02;
    COFD[9188] = -2.14562897E+01;
    COFD[9189] = 5.33532691E+00;
    COFD[9190] = -4.56505568E-01;
    COFD[9191] = 1.89262566E-02;
    COFD[9192] = -2.08516224E+01;
    COFD[9193] = 5.07617318E+00;
    COFD[9194] = -4.28983517E-01;
    COFD[9195] = 1.79703565E-02;
    COFD[9196] = -2.08624149E+01;
    COFD[9197] = 5.07617318E+00;
    COFD[9198] = -4.28983517E-01;
    COFD[9199] = 1.79703565E-02;
    COFD[9200] = -2.12488262E+01;
    COFD[9201] = 4.69672574E+00;
    COFD[9202] = -3.24393524E-01;
    COFD[9203] = 1.13541249E-02;
    COFD[9204] = -2.26588661E+01;
    COFD[9205] = 5.48378580E+00;
    COFD[9206] = -4.57964754E-01;
    COFD[9207] = 1.83088022E-02;
    COFD[9208] = -2.24277170E+01;
    COFD[9209] = 5.48901287E+00;
    COFD[9210] = -4.60031399E-01;
    COFD[9211] = 1.84508313E-02;
    COFD[9212] = -2.26462126E+01;
    COFD[9213] = 5.48378580E+00;
    COFD[9214] = -4.57964754E-01;
    COFD[9215] = 1.83088022E-02;
    COFD[9216] = -2.24904947E+01;
    COFD[9217] = 5.47343682E+00;
    COFD[9218] = -4.54166709E-01;
    COFD[9219] = 1.80526774E-02;
    COFD[9220] = -2.25056806E+01;
    COFD[9221] = 5.47343682E+00;
    COFD[9222] = -4.54166709E-01;
    COFD[9223] = 1.80526774E-02;
    COFD[9224] = -2.14352422E+01;
    COFD[9225] = 4.77388530E+00;
    COFD[9226] = -3.36183750E-01;
    COFD[9227] = 1.19302239E-02;
    COFD[9228] = -2.19740179E+01;
    COFD[9229] = 4.96635252E+00;
    COFD[9230] = -3.66267779E-01;
    COFD[9231] = 1.34216831E-02;
    COFD[9232] = -2.17030104E+01;
    COFD[9233] = 5.39237617E+00;
    COFD[9234] = -4.62530602E-01;
    COFD[9235] = 1.91345681E-02;
    COFD[9236] = -2.19895742E+01;
    COFD[9237] = 4.96635252E+00;
    COFD[9238] = -3.66267779E-01;
    COFD[9239] = 1.34216831E-02;
    COFD[9240] = -2.26809417E+01;
    COFD[9241] = 5.35408130E+00;
    COFD[9242] = -4.30433039E-01;
    COFD[9243] = 1.67163724E-02;
    COFD[9244] = -2.26809417E+01;
    COFD[9245] = 5.35408130E+00;
    COFD[9246] = -4.30433039E-01;
    COFD[9247] = 1.67163724E-02;
    COFD[9248] = -2.27276855E+01;
    COFD[9249] = 5.40059960E+00;
    COFD[9250] = -4.38830427E-01;
    COFD[9251] = 1.71686343E-02;
    COFD[9252] = -2.27156370E+01;
    COFD[9253] = 5.25593075E+00;
    COFD[9254] = -4.13540252E-01;
    COFD[9255] = 1.58298021E-02;
    COFD[9256] = -2.27176356E+01;
    COFD[9257] = 5.26237767E+00;
    COFD[9258] = -4.14629620E-01;
    COFD[9259] = 1.58863936E-02;
    COFD[9260] = -2.27231185E+01;
    COFD[9261] = 5.26237767E+00;
    COFD[9262] = -4.14629620E-01;
    COFD[9263] = 1.58863936E-02;
    COFD[9264] = -2.28592728E+01;
    COFD[9265] = 5.27229756E+00;
    COFD[9266] = -4.16295172E-01;
    COFD[9267] = 1.59725857E-02;
    COFD[9268] = -2.22576038E+01;
    COFD[9269] = 4.96852134E+00;
    COFD[9270] = -3.66607756E-01;
    COFD[9271] = 1.34385708E-02;
    COFD[9272] = -2.22261321E+01;
    COFD[9273] = 4.96982327E+00;
    COFD[9274] = -3.66811857E-01;
    COFD[9275] = 1.34487095E-02;
    COFD[9276] = -2.27064705E+01;
    COFD[9277] = 5.12211189E+00;
    COFD[9278] = -3.91225928E-01;
    COFD[9279] = 1.46792378E-02;
    COFD[9280] = -2.27719473E+01;
    COFD[9281] = 5.15606268E+00;
    COFD[9282] = -3.96907901E-01;
    COFD[9283] = 1.49727202E-02;
    COFD[9284] = -2.23118869E+01;
    COFD[9285] = 4.94556893E+00;
    COFD[9286] = -3.63013720E-01;
    COFD[9287] = 1.32601068E-02;
    COFD[9288] = -2.15669461E+01;
    COFD[9289] = 4.77388530E+00;
    COFD[9290] = -3.36183750E-01;
    COFD[9291] = 1.19302239E-02;
    COFD[9292] = -2.15737493E+01;
    COFD[9293] = 4.77388530E+00;
    COFD[9294] = -3.36183750E-01;
    COFD[9295] = 1.19302239E-02;
    COFD[9296] = -2.20318339E+01;
    COFD[9297] = 4.93482700E+00;
    COFD[9298] = -3.61331756E-01;
    COFD[9299] = 1.31765767E-02;
    COFD[9300] = -2.23257344E+01;
    COFD[9301] = 4.99720991E+00;
    COFD[9302] = -3.71116272E-01;
    COFD[9303] = 1.36628560E-02;
    COFD[9304] = -2.27109978E+01;
    COFD[9305] = 5.37975560E+00;
    COFD[9306] = -4.34974624E-01;
    COFD[9307] = 1.69583701E-02;
    COFD[9308] = -2.26722436E+01;
    COFD[9309] = 5.35408130E+00;
    COFD[9310] = -4.30433039E-01;
    COFD[9311] = 1.67163724E-02;
    COFD[9312] = -2.21244697E+01;
    COFD[9313] = 4.84823431E+00;
    COFD[9314] = -3.47731563E-01;
    COFD[9315] = 1.25001760E-02;
    COFD[9316] = -2.22786912E+01;
    COFD[9317] = 4.89093124E+00;
    COFD[9318] = -3.54419517E-01;
    COFD[9319] = 1.28322610E-02;
    COFD[9320] = -2.20373227E+01;
    COFD[9321] = 4.93482700E+00;
    COFD[9322] = -3.61331756E-01;
    COFD[9323] = 1.31765767E-02;
    COFD[9324] = -2.22044849E+01;
    COFD[9325] = 4.79539928E+00;
    COFD[9326] = -3.39494765E-01;
    COFD[9327] = 1.20926677E-02;
    COFD[9328] = -2.20198043E+01;
    COFD[9329] = 4.73965868E+00;
    COFD[9330] = -3.30939918E-01;
    COFD[9331] = 1.16735493E-02;
    COFD[9332] = -2.24702000E+01;
    COFD[9333] = 4.86736132E+00;
    COFD[9334] = -3.50725707E-01;
    COFD[9335] = 1.26487745E-02;
    COFD[9336] = -2.24678718E+01;
    COFD[9337] = 4.86736132E+00;
    COFD[9338] = -3.50725707E-01;
    COFD[9339] = 1.26487745E-02;
    COFD[9340] = -2.16636303E+01;
    COFD[9341] = 4.42530942E+00;
    COFD[9342] = -2.82800211E-01;
    COFD[9343] = 9.31937239E-03;
    COFD[9344] = -2.14832299E+01;
    COFD[9345] = 4.24469152E+00;
    COFD[9346] = -2.55617920E-01;
    COFD[9347] = 8.00605387E-03;
    COFD[9348] = -2.20449981E+01;
    COFD[9349] = 4.63376397E+00;
    COFD[9350] = -3.14793416E-01;
    COFD[9351] = 1.08858273E-02;
    COFD[9352] = -2.15271080E+01;
    COFD[9353] = 4.33297139E+00;
    COFD[9354] = -2.68886993E-01;
    COFD[9355] = 8.64665258E-03;
    COFD[9356] = -2.04878867E+01;
    COFD[9357] = 4.96174500E+00;
    COFD[9358] = -4.15861207E-01;
    COFD[9359] = 1.74690994E-02;
    COFD[9360] = -1.98430923E+01;
    COFD[9361] = 5.29318575E+00;
    COFD[9362] = -4.51815332E-01;
    COFD[9363] = 1.87526763E-02;
    COFD[9364] = -1.59092604E+01;
    COFD[9365] = 3.75176869E+00;
    COFD[9366] = -2.70323254E-01;
    COFD[9367] = 1.16183021E-02;
    COFD[9368] = -1.93555562E+01;
    COFD[9369] = 4.64408224E+00;
    COFD[9370] = -3.78774559E-01;
    COFD[9371] = 1.60247845E-02;
    COFD[9372] = -2.08510201E+01;
    COFD[9373] = 5.00799007E+00;
    COFD[9374] = -4.21303846E-01;
    COFD[9375] = 1.76831398E-02;
    COFD[9376] = -1.93817856E+01;
    COFD[9377] = 4.64408224E+00;
    COFD[9378] = -3.78774559E-01;
    COFD[9379] = 1.60247845E-02;
    COFD[9380] = -2.07027595E+01;
    COFD[9381] = 4.52425469E+00;
    COFD[9382] = -2.97869560E-01;
    COFD[9383] = 1.00531031E-02;
    COFD[9384] = -2.05041077E+01;
    COFD[9385] = 4.89401297E+00;
    COFD[9386] = -4.08087781E-01;
    COFD[9387] = 1.71722700E-02;
    COFD[9388] = -2.28305301E+01;
    COFD[9389] = 5.49066995E+00;
    COFD[9390] = -4.61705340E-01;
    COFD[9391] = 1.85789333E-02;
    COFD[9392] = -2.14866930E+01;
    COFD[9393] = 5.28483859E+00;
    COFD[9394] = -4.50875157E-01;
    COFD[9395] = 1.87173937E-02;
    COFD[9396] = -2.14364392E+01;
    COFD[9397] = 5.26245213E+00;
    COFD[9398] = -4.48352237E-01;
    COFD[9399] = 1.86227443E-02;
    COFD[9400] = -2.08627296E+01;
    COFD[9401] = 5.00799007E+00;
    COFD[9402] = -4.21303846E-01;
    COFD[9403] = 1.76831398E-02;
    COFD[9404] = -2.08740011E+01;
    COFD[9405] = 5.00799007E+00;
    COFD[9406] = -4.21303846E-01;
    COFD[9407] = 1.76831398E-02;
    COFD[9408] = -2.17305043E+01;
    COFD[9409] = 4.83264094E+00;
    COFD[9410] = -3.45294821E-01;
    COFD[9411] = 1.23794137E-02;
    COFD[9412] = -2.28191084E+01;
    COFD[9413] = 5.49118368E+00;
    COFD[9414] = -4.61244560E-01;
    COFD[9415] = 1.85388488E-02;
    COFD[9416] = -2.25951559E+01;
    COFD[9417] = 5.49332573E+00;
    COFD[9418] = -4.62935000E-01;
    COFD[9419] = 1.86660283E-02;
    COFD[9420] = -2.28059444E+01;
    COFD[9421] = 5.49118368E+00;
    COFD[9422] = -4.61244560E-01;
    COFD[9423] = 1.85388488E-02;
    COFD[9424] = -2.26677760E+01;
    COFD[9425] = 5.48311530E+00;
    COFD[9426] = -4.57619182E-01;
    COFD[9427] = 1.82841374E-02;
    COFD[9428] = -2.26835080E+01;
    COFD[9429] = 5.48311530E+00;
    COFD[9430] = -4.57619182E-01;
    COFD[9431] = 1.82841374E-02;
    COFD[9432] = -2.18839486E+01;
    COFD[9433] = 4.89616559E+00;
    COFD[9434] = -3.55241555E-01;
    COFD[9435] = 1.28731550E-02;
    COFD[9436] = -2.23860697E+01;
    COFD[9437] = 5.07468436E+00;
    COFD[9438] = -3.83473296E-01;
    COFD[9439] = 1.42836887E-02;
    COFD[9440] = -2.17533197E+01;
    COFD[9441] = 5.33181444E+00;
    COFD[9442] = -4.56115824E-01;
    COFD[9443] = 1.89118833E-02;
    COFD[9444] = -2.24024538E+01;
    COFD[9445] = 5.07468436E+00;
    COFD[9446] = -3.83473296E-01;
    COFD[9447] = 1.42836887E-02;
    COFD[9448] = -2.29546228E+01;
    COFD[9449] = 5.40726595E+00;
    COFD[9450] = -4.40072331E-01;
    COFD[9451] = 1.72366684E-02;
    COFD[9452] = -2.29546228E+01;
    COFD[9453] = 5.40726595E+00;
    COFD[9454] = -4.40072331E-01;
    COFD[9455] = 1.72366684E-02;
    COFD[9456] = -2.29808472E+01;
    COFD[9457] = 5.44447552E+00;
    COFD[9458] = -4.47281890E-01;
    COFD[9459] = 1.76387634E-02;
    COFD[9460] = -2.30198383E+01;
    COFD[9461] = 5.32186706E+00;
    COFD[9462] = -4.24826256E-01;
    COFD[9463] = 1.64201917E-02;
    COFD[9464] = -2.30219288E+01;
    COFD[9465] = 5.32771617E+00;
    COFD[9466] = -4.25833456E-01;
    COFD[9467] = 1.64730816E-02;
    COFD[9468] = -2.30277582E+01;
    COFD[9469] = 5.32771617E+00;
    COFD[9470] = -4.25833456E-01;
    COFD[9471] = 1.64730816E-02;
    COFD[9472] = -2.31520360E+01;
    COFD[9473] = 5.33643523E+00;
    COFD[9474] = -4.27338696E-01;
    COFD[9475] = 1.65522282E-02;
    COFD[9476] = -2.26588219E+01;
    COFD[9477] = 5.07646321E+00;
    COFD[9478] = -3.83762054E-01;
    COFD[9479] = 1.42983651E-02;
    COFD[9480] = -2.26303193E+01;
    COFD[9481] = 5.07752972E+00;
    COFD[9482] = -3.83935232E-01;
    COFD[9483] = 1.43071687E-02;
    COFD[9484] = -2.30240568E+01;
    COFD[9485] = 5.19737465E+00;
    COFD[9486] = -4.03890980E-01;
    COFD[9487] = 1.53353180E-02;
    COFD[9488] = -2.31022549E+01;
    COFD[9489] = 5.23643800E+00;
    COFD[9490] = -4.10255644E-01;
    COFD[9491] = 1.56593397E-02;
    COFD[9492] = -2.27251033E+01;
    COFD[9493] = 5.05775235E+00;
    COFD[9494] = -3.80730307E-01;
    COFD[9495] = 1.41444474E-02;
    COFD[9496] = -2.20221614E+01;
    COFD[9497] = 4.89616559E+00;
    COFD[9498] = -3.55241555E-01;
    COFD[9499] = 1.28731550E-02;
    COFD[9500] = -2.20293523E+01;
    COFD[9501] = 4.89616559E+00;
    COFD[9502] = -3.55241555E-01;
    COFD[9503] = 1.28731550E-02;
    COFD[9504] = -2.24600223E+01;
    COFD[9505] = 5.04880598E+00;
    COFD[9506] = -3.79285322E-01;
    COFD[9507] = 1.40712210E-02;
    COFD[9508] = -2.27123285E+01;
    COFD[9509] = 5.10007015E+00;
    COFD[9510] = -3.87606161E-01;
    COFD[9511] = 1.44941102E-02;
    COFD[9512] = -2.29717221E+01;
    COFD[9513] = 5.42721943E+00;
    COFD[9514] = -4.43876872E-01;
    COFD[9515] = 1.74472689E-02;
    COFD[9516] = -2.29454887E+01;
    COFD[9517] = 5.40726595E+00;
    COFD[9518] = -4.40072331E-01;
    COFD[9519] = 1.72366684E-02;
    COFD[9520] = -2.25506266E+01;
    COFD[9521] = 4.96541553E+00;
    COFD[9522] = -3.66120907E-01;
    COFD[9523] = 1.34143879E-02;
    COFD[9524] = -2.26945043E+01;
    COFD[9525] = 5.00634977E+00;
    COFD[9526] = -3.72557302E-01;
    COFD[9527] = 1.37346802E-02;
    COFD[9528] = -2.24658579E+01;
    COFD[9529] = 5.04880598E+00;
    COFD[9530] = -3.79285322E-01;
    COFD[9531] = 1.40712210E-02;
    COFD[9532] = -2.26263872E+01;
    COFD[9533] = 4.91358137E+00;
    COFD[9534] = -3.57998484E-01;
    COFD[9535] = 1.30108794E-02;
    COFD[9536] = -2.24702000E+01;
    COFD[9537] = 4.86736132E+00;
    COFD[9538] = -3.50725707E-01;
    COFD[9539] = 1.26487745E-02;
    COFD[9540] = -2.28846873E+01;
    COFD[9541] = 4.98408738E+00;
    COFD[9542] = -3.69050697E-01;
    COFD[9543] = 1.35600028E-02;
    COFD[9544] = -2.28821534E+01;
    COFD[9545] = 4.98408738E+00;
    COFD[9546] = -3.69050697E-01;
    COFD[9547] = 1.35600028E-02;
    COFD[9548] = -2.21485602E+01;
    COFD[9549] = 4.56927551E+00;
    COFD[9550] = -3.04798371E-01;
    COFD[9551] = 1.03930919E-02;
    COFD[9552] = -2.20144983E+01;
    COFD[9553] = 4.41014015E+00;
    COFD[9554] = -2.80504988E-01;
    COFD[9555] = 9.20806328E-03;
    COFD[9556] = -2.25315042E+01;
    COFD[9557] = 4.77923223E+00;
    COFD[9558] = -3.37003339E-01;
    COFD[9559] = 1.19703492E-02;
    COFD[9560] = -2.20447991E+01;
    COFD[9561] = 4.49024164E+00;
    COFD[9562] = -2.92656701E-01;
    COFD[9563] = 9.79834494E-03;
    COFD[9564] = -2.04827540E+01;
    COFD[9565] = 4.88918294E+00;
    COFD[9566] = -4.07583331E-01;
    COFD[9567] = 1.71552617E-02;
    COFD[9568] = -1.98430417E+01;
    COFD[9569] = 5.29318575E+00;
    COFD[9570] = -4.51815332E-01;
    COFD[9571] = 1.87526763E-02;
    COFD[9572] = -1.59091602E+01;
    COFD[9573] = 3.75176869E+00;
    COFD[9574] = -2.70323254E-01;
    COFD[9575] = 1.16183021E-02;
    COFD[9576] = -1.93548571E+01;
    COFD[9577] = 4.64408224E+00;
    COFD[9578] = -3.78774559E-01;
    COFD[9579] = 1.60247845E-02;
    COFD[9580] = -2.08497919E+01;
    COFD[9581] = 5.00799007E+00;
    COFD[9582] = -4.21303846E-01;
    COFD[9583] = 1.76831398E-02;
    COFD[9584] = -1.93810490E+01;
    COFD[9585] = 4.64408224E+00;
    COFD[9586] = -3.78774559E-01;
    COFD[9587] = 1.60247845E-02;
    COFD[9588] = -2.07019859E+01;
    COFD[9589] = 4.52425469E+00;
    COFD[9590] = -2.97869560E-01;
    COFD[9591] = 1.00531031E-02;
    COFD[9592] = -2.05029990E+01;
    COFD[9593] = 4.89401297E+00;
    COFD[9594] = -4.08087781E-01;
    COFD[9595] = 1.71722700E-02;
    COFD[9596] = -2.28289821E+01;
    COFD[9597] = 5.49066995E+00;
    COFD[9598] = -4.61705340E-01;
    COFD[9599] = 1.85789333E-02;
    COFD[9600] = -2.14860306E+01;
    COFD[9601] = 5.28483859E+00;
    COFD[9602] = -4.50875157E-01;
    COFD[9603] = 1.87173937E-02;
    COFD[9604] = -2.14357385E+01;
    COFD[9605] = 5.26245213E+00;
    COFD[9606] = -4.48352237E-01;
    COFD[9607] = 1.86227443E-02;
    COFD[9608] = -2.08614723E+01;
    COFD[9609] = 5.00799007E+00;
    COFD[9610] = -4.21303846E-01;
    COFD[9611] = 1.76831398E-02;
    COFD[9612] = -2.08727152E+01;
    COFD[9613] = 5.00799007E+00;
    COFD[9614] = -4.21303846E-01;
    COFD[9615] = 1.76831398E-02;
    COFD[9616] = -2.17293342E+01;
    COFD[9617] = 4.83264094E+00;
    COFD[9618] = -3.45294821E-01;
    COFD[9619] = 1.23794137E-02;
    COFD[9620] = -2.28179371E+01;
    COFD[9621] = 5.49118368E+00;
    COFD[9622] = -4.61244560E-01;
    COFD[9623] = 1.85388488E-02;
    COFD[9624] = -2.25940458E+01;
    COFD[9625] = 5.49332573E+00;
    COFD[9626] = -4.62935000E-01;
    COFD[9627] = 1.86660283E-02;
    COFD[9628] = -2.28048034E+01;
    COFD[9629] = 5.49118368E+00;
    COFD[9630] = -4.61244560E-01;
    COFD[9631] = 1.85388488E-02;
    COFD[9632] = -2.26667616E+01;
    COFD[9633] = 5.48311530E+00;
    COFD[9634] = -4.57619182E-01;
    COFD[9635] = 1.82841374E-02;
    COFD[9636] = -2.26824612E+01;
    COFD[9637] = 5.48311530E+00;
    COFD[9638] = -4.57619182E-01;
    COFD[9639] = 1.82841374E-02;
    COFD[9640] = -2.18827191E+01;
    COFD[9641] = 4.89616559E+00;
    COFD[9642] = -3.55241555E-01;
    COFD[9643] = 1.28731550E-02;
    COFD[9644] = -2.23845704E+01;
    COFD[9645] = 5.07468436E+00;
    COFD[9646] = -3.83473296E-01;
    COFD[9647] = 1.42836887E-02;
    COFD[9648] = -2.17518459E+01;
    COFD[9649] = 5.33181444E+00;
    COFD[9650] = -4.56115824E-01;
    COFD[9651] = 1.89118833E-02;
    COFD[9652] = -2.24009047E+01;
    COFD[9653] = 5.07468436E+00;
    COFD[9654] = -3.83473296E-01;
    COFD[9655] = 1.42836887E-02;
    COFD[9656] = -2.29531738E+01;
    COFD[9657] = 5.40726595E+00;
    COFD[9658] = -4.40072331E-01;
    COFD[9659] = 1.72366684E-02;
    COFD[9660] = -2.29531738E+01;
    COFD[9661] = 5.40726595E+00;
    COFD[9662] = -4.40072331E-01;
    COFD[9663] = 1.72366684E-02;
    COFD[9664] = -2.29793469E+01;
    COFD[9665] = 5.44447552E+00;
    COFD[9666] = -4.47281890E-01;
    COFD[9667] = 1.76387634E-02;
    COFD[9668] = -2.30180603E+01;
    COFD[9669] = 5.32186706E+00;
    COFD[9670] = -4.24826256E-01;
    COFD[9671] = 1.64201917E-02;
    COFD[9672] = -2.30201295E+01;
    COFD[9673] = 5.32771617E+00;
    COFD[9674] = -4.25833456E-01;
    COFD[9675] = 1.64730816E-02;
    COFD[9676] = -2.30259379E+01;
    COFD[9677] = 5.32771617E+00;
    COFD[9678] = -4.25833456E-01;
    COFD[9679] = 1.64730816E-02;
    COFD[9680] = -2.31501948E+01;
    COFD[9681] = 5.33643523E+00;
    COFD[9682] = -4.27338696E-01;
    COFD[9683] = 1.65522282E-02;
    COFD[9684] = -2.26569816E+01;
    COFD[9685] = 5.07646321E+00;
    COFD[9686] = -3.83762054E-01;
    COFD[9687] = 1.42983651E-02;
    COFD[9688] = -2.26284586E+01;
    COFD[9689] = 5.07752972E+00;
    COFD[9690] = -3.83935232E-01;
    COFD[9691] = 1.43071687E-02;
    COFD[9692] = -2.30219870E+01;
    COFD[9693] = 5.19737465E+00;
    COFD[9694] = -4.03890980E-01;
    COFD[9695] = 1.53353180E-02;
    COFD[9696] = -2.31001674E+01;
    COFD[9697] = 5.23643800E+00;
    COFD[9698] = -4.10255644E-01;
    COFD[9699] = 1.56593397E-02;
    COFD[9700] = -2.27229983E+01;
    COFD[9701] = 5.05775235E+00;
    COFD[9702] = -3.80730307E-01;
    COFD[9703] = 1.41444474E-02;
    COFD[9704] = -2.20205411E+01;
    COFD[9705] = 4.89616559E+00;
    COFD[9706] = -3.55241555E-01;
    COFD[9707] = 1.28731550E-02;
    COFD[9708] = -2.20277086E+01;
    COFD[9709] = 4.89616559E+00;
    COFD[9710] = -3.55241555E-01;
    COFD[9711] = 1.28731550E-02;
    COFD[9712] = -2.24582240E+01;
    COFD[9713] = 5.04880598E+00;
    COFD[9714] = -3.79285322E-01;
    COFD[9715] = 1.40712210E-02;
    COFD[9716] = -2.27105091E+01;
    COFD[9717] = 5.10007015E+00;
    COFD[9718] = -3.87606161E-01;
    COFD[9719] = 1.44941102E-02;
    COFD[9720] = -2.29702473E+01;
    COFD[9721] = 5.42721943E+00;
    COFD[9722] = -4.43876872E-01;
    COFD[9723] = 1.74472689E-02;
    COFD[9724] = -2.29440659E+01;
    COFD[9725] = 5.40726595E+00;
    COFD[9726] = -4.40072331E-01;
    COFD[9727] = 1.72366684E-02;
    COFD[9728] = -2.25485050E+01;
    COFD[9729] = 4.96541553E+00;
    COFD[9730] = -3.66120907E-01;
    COFD[9731] = 1.34143879E-02;
    COFD[9732] = -2.26924001E+01;
    COFD[9733] = 5.00634977E+00;
    COFD[9734] = -3.72557302E-01;
    COFD[9735] = 1.37346802E-02;
    COFD[9736] = -2.24640385E+01;
    COFD[9737] = 5.04880598E+00;
    COFD[9738] = -3.79285322E-01;
    COFD[9739] = 1.40712210E-02;
    COFD[9740] = -2.26240442E+01;
    COFD[9741] = 4.91358137E+00;
    COFD[9742] = -3.57998484E-01;
    COFD[9743] = 1.30108794E-02;
    COFD[9744] = -2.24678718E+01;
    COFD[9745] = 4.86736132E+00;
    COFD[9746] = -3.50725707E-01;
    COFD[9747] = 1.26487745E-02;
    COFD[9748] = -2.28821534E+01;
    COFD[9749] = 4.98408738E+00;
    COFD[9750] = -3.69050697E-01;
    COFD[9751] = 1.35600028E-02;
    COFD[9752] = -2.28796323E+01;
    COFD[9753] = 4.98408738E+00;
    COFD[9754] = -3.69050697E-01;
    COFD[9755] = 1.35600028E-02;
    COFD[9756] = -2.21456879E+01;
    COFD[9757] = 4.56927551E+00;
    COFD[9758] = -3.04798371E-01;
    COFD[9759] = 1.03930919E-02;
    COFD[9760] = -2.20113604E+01;
    COFD[9761] = 4.41014015E+00;
    COFD[9762] = -2.80504988E-01;
    COFD[9763] = 9.20806328E-03;
    COFD[9764] = -2.25288055E+01;
    COFD[9765] = 4.77923223E+00;
    COFD[9766] = -3.37003339E-01;
    COFD[9767] = 1.19703492E-02;
    COFD[9768] = -2.20417938E+01;
    COFD[9769] = 4.49024164E+00;
    COFD[9770] = -2.92656701E-01;
    COFD[9771] = 9.79834494E-03;
    COFD[9772] = -2.04816453E+01;
    COFD[9773] = 4.88918294E+00;
    COFD[9774] = -4.07583331E-01;
    COFD[9775] = 1.71552617E-02;
    COFD[9776] = -2.03662784E+01;
    COFD[9777] = 5.43917732E+00;
    COFD[9778] = -4.65585169E-01;
    COFD[9779] = 1.91504194E-02;
    COFD[9780] = -1.65211535E+01;
    COFD[9781] = 3.96045057E+00;
    COFD[9782] = -2.95491083E-01;
    COFD[9783] = 1.26309912E-02;
    COFD[9784] = -2.01011659E+01;
    COFD[9785] = 4.89020393E+00;
    COFD[9786] = -4.07690029E-01;
    COFD[9787] = 1.71588628E-02;
    COFD[9788] = -2.15001622E+01;
    COFD[9789] = 5.19791528E+00;
    COFD[9790] = -4.41398852E-01;
    COFD[9791] = 1.83766038E-02;
    COFD[9792] = -2.01283014E+01;
    COFD[9793] = 4.89020393E+00;
    COFD[9794] = -4.07690029E-01;
    COFD[9795] = 1.71588628E-02;
    COFD[9796] = -1.90849698E+01;
    COFD[9797] = 3.76037078E+00;
    COFD[9798] = -1.86344485E-01;
    COFD[9799] = 4.78342458E-03;
    COFD[9800] = -2.12596231E+01;
    COFD[9801] = 5.13107283E+00;
    COFD[9802] = -4.34832558E-01;
    COFD[9803] = 1.81741958E-02;
    COFD[9804] = -2.29703215E+01;
    COFD[9805] = 5.45099312E+00;
    COFD[9806] = -4.48927722E-01;
    COFD[9807] = 1.77397382E-02;
    COFD[9808] = -2.20377114E+01;
    COFD[9809] = 5.43771555E+00;
    COFD[9810] = -4.65617512E-01;
    COFD[9811] = 1.91604198E-02;
    COFD[9812] = -2.20523530E+01;
    COFD[9813] = 5.43694233E+00;
    COFD[9814] = -4.65634840E-01;
    COFD[9815] = 1.91657287E-02;
    COFD[9816] = -2.14847798E+01;
    COFD[9817] = 5.19080403E+00;
    COFD[9818] = -4.40813300E-01;
    COFD[9819] = 1.83643944E-02;
    COFD[9820] = -2.14967605E+01;
    COFD[9821] = 5.19080403E+00;
    COFD[9822] = -4.40813300E-01;
    COFD[9823] = 1.83643944E-02;
    COFD[9824] = -2.08674131E+01;
    COFD[9825] = 4.37860613E+00;
    COFD[9826] = -2.75749639E-01;
    COFD[9827] = 8.97792922E-03;
    COFD[9828] = -2.29332780E+01;
    COFD[9829] = 5.44974657E+00;
    COFD[9830] = -4.48651670E-01;
    COFD[9831] = 1.77235430E-02;
    COFD[9832] = -2.27385006E+01;
    COFD[9833] = 5.46322829E+00;
    COFD[9834] = -4.51749940E-01;
    COFD[9835] = 1.79075310E-02;
    COFD[9836] = -2.29193653E+01;
    COFD[9837] = 5.44974657E+00;
    COFD[9838] = -4.48651670E-01;
    COFD[9839] = 1.77235430E-02;
    COFD[9840] = -2.27048922E+01;
    COFD[9841] = 5.40961679E+00;
    COFD[9842] = -4.40511249E-01;
    COFD[9843] = 1.72607396E-02;
    COFD[9844] = -2.27214167E+01;
    COFD[9845] = 5.40961679E+00;
    COFD[9846] = -4.40511249E-01;
    COFD[9847] = 1.72607396E-02;
    COFD[9848] = -2.10721493E+01;
    COFD[9849] = 4.46261834E+00;
    COFD[9850] = -2.88453410E-01;
    COFD[9851] = 9.59378935E-03;
    COFD[9852] = -2.16898608E+01;
    COFD[9853] = 4.68658569E+00;
    COFD[9854] = -3.22849054E-01;
    COFD[9855] = 1.12788192E-02;
    COFD[9856] = -2.22344688E+01;
    COFD[9857] = 5.44527718E+00;
    COFD[9858] = -4.65287801E-01;
    COFD[9859] = 1.90954247E-02;
    COFD[9860] = -2.17074965E+01;
    COFD[9861] = 4.68658569E+00;
    COFD[9862] = -3.22849054E-01;
    COFD[9863] = 1.12788192E-02;
    COFD[9864] = -2.26847799E+01;
    COFD[9865] = 5.19870800E+00;
    COFD[9866] = -4.04106072E-01;
    COFD[9867] = 1.53462074E-02;
    COFD[9868] = -2.26847799E+01;
    COFD[9869] = 5.19870800E+00;
    COFD[9870] = -4.04106072E-01;
    COFD[9871] = 1.53462074E-02;
    COFD[9872] = -2.27935663E+01;
    COFD[9873] = 5.26880815E+00;
    COFD[9874] = -4.15708807E-01;
    COFD[9875] = 1.59422297E-02;
    COFD[9876] = -2.26694016E+01;
    COFD[9877] = 5.07552523E+00;
    COFD[9878] = -3.83609779E-01;
    COFD[9879] = 1.42906252E-02;
    COFD[9880] = -2.26785584E+01;
    COFD[9881] = 5.08419101E+00;
    COFD[9882] = -3.85017801E-01;
    COFD[9883] = 1.43622302E-02;
    COFD[9884] = -2.26849244E+01;
    COFD[9885] = 5.08419101E+00;
    COFD[9886] = -3.85017801E-01;
    COFD[9887] = 1.43622302E-02;
    COFD[9888] = -2.28187224E+01;
    COFD[9889] = 5.09710954E+00;
    COFD[9890] = -3.87121977E-01;
    COFD[9891] = 1.44693985E-02;
    COFD[9892] = -2.19726300E+01;
    COFD[9893] = 4.68933904E+00;
    COFD[9894] = -3.23268395E-01;
    COFD[9895] = 1.12992644E-02;
    COFD[9896] = -2.14841962E+01;
    COFD[9897] = 4.49540264E+00;
    COFD[9898] = -2.96514053E-01;
    COFD[9899] = 1.00871155E-02;
    COFD[9900] = -2.25275251E+01;
    COFD[9901] = 4.88611782E+00;
    COFD[9902] = -3.53664413E-01;
    COFD[9903] = 1.27947308E-02;
    COFD[9904] = -2.25844328E+01;
    COFD[9905] = 4.91727638E+00;
    COFD[9906] = -3.59056497E-01;
    COFD[9907] = 1.30783726E-02;
    COFD[9908] = -2.20211813E+01;
    COFD[9909] = 4.66069984E+00;
    COFD[9910] = -3.18907247E-01;
    COFD[9911] = 1.10866672E-02;
    COFD[9912] = -2.12201218E+01;
    COFD[9913] = 4.46261834E+00;
    COFD[9914] = -2.88453410E-01;
    COFD[9915] = 9.59378935E-03;
    COFD[9916] = -2.12279042E+01;
    COFD[9917] = 4.46261834E+00;
    COFD[9918] = -2.88453410E-01;
    COFD[9919] = 9.59378935E-03;
    COFD[9920] = -2.17388998E+01;
    COFD[9921] = 4.64756828E+00;
    COFD[9922] = -3.16905849E-01;
    COFD[9923] = 1.09890709E-02;
    COFD[9924] = -2.16032600E+01;
    COFD[9925] = 4.53358880E+00;
    COFD[9926] = -3.02280026E-01;
    COFD[9927] = 1.03666188E-02;
    COFD[9928] = -2.27497015E+01;
    COFD[9929] = 5.23793490E+00;
    COFD[9930] = -4.10507434E-01;
    COFD[9931] = 1.56723961E-02;
    COFD[9932] = -2.26749914E+01;
    COFD[9933] = 5.19870800E+00;
    COFD[9934] = -4.04106072E-01;
    COFD[9935] = 1.53462074E-02;
    COFD[9936] = -2.13962129E+01;
    COFD[9937] = 4.37956579E+00;
    COFD[9938] = -2.78446618E-01;
    COFD[9939] = 9.19232241E-03;
    COFD[9940] = -2.14928955E+01;
    COFD[9941] = 4.39966514E+00;
    COFD[9942] = -2.81886234E-01;
    COFD[9943] = 9.37271319E-03;
    COFD[9944] = -2.17452723E+01;
    COFD[9945] = 4.64756828E+00;
    COFD[9946] = -3.16905849E-01;
    COFD[9947] = 1.09890709E-02;
    COFD[9948] = -2.15478797E+01;
    COFD[9949] = 4.35915869E+00;
    COFD[9950] = -2.74825098E-01;
    COFD[9951] = 8.99908615E-03;
    COFD[9952] = -2.16636303E+01;
    COFD[9953] = 4.42530942E+00;
    COFD[9954] = -2.82800211E-01;
    COFD[9955] = 9.31937239E-03;
    COFD[9956] = -2.21485602E+01;
    COFD[9957] = 4.56927551E+00;
    COFD[9958] = -3.04798371E-01;
    COFD[9959] = 1.03930919E-02;
    COFD[9960] = -2.21456879E+01;
    COFD[9961] = 4.56927551E+00;
    COFD[9962] = -3.04798371E-01;
    COFD[9963] = 1.03930919E-02;
    COFD[9964] = -2.10080955E+01;
    COFD[9965] = 3.98704561E+00;
    COFD[9966] = -2.18285980E-01;
    COFD[9967] = 6.25051699E-03;
    COFD[9968] = -2.07325259E+01;
    COFD[9969] = 3.76716506E+00;
    COFD[9970] = -1.85775467E-01;
    COFD[9971] = 4.69798135E-03;
    COFD[9972] = -2.16640062E+01;
    COFD[9973] = 4.30841644E+00;
    COFD[9974] = -2.65200479E-01;
    COFD[9975] = 8.46889236E-03;
    COFD[9976] = -2.07979763E+01;
    COFD[9977] = 3.86302522E+00;
    COFD[9978] = -2.00078966E-01;
    COFD[9979] = 5.38570099E-03;
    COFD[9980] = -2.12317255E+01;
    COFD[9981] = 5.12393233E+00;
    COFD[9982] = -4.34047487E-01;
    COFD[9983] = 1.81456921E-02;
    COFD[9984] = -2.06638821E+01;
    COFD[9985] = 5.45347668E+00;
    COFD[9986] = -4.65289819E-01;
    COFD[9987] = 1.90547288E-02;
    COFD[9988] = -1.69048097E+01;
    COFD[9989] = 4.02910047E+00;
    COFD[9990] = -3.03801059E-01;
    COFD[9991] = 1.29668774E-02;
    COFD[9992] = -2.05451641E+01;
    COFD[9993] = 4.97155662E+00;
    COFD[9994] = -4.17014903E-01;
    COFD[9995] = 1.75144222E-02;
    COFD[9996] = -2.19077574E+01;
    COFD[9997] = 5.26494284E+00;
    COFD[9998] = -4.48633355E-01;
    COFD[9999] = 1.86333022E-02;
    COFD[10000] = -2.05729073E+01;
    COFD[10001] = 4.97155662E+00;
    COFD[10002] = -4.17014903E-01;
    COFD[10003] = 1.75144222E-02;
    COFD[10004] = -1.89362802E+01;
    COFD[10005] = 3.59051852E+00;
    COFD[10006] = -1.61030577E-01;
    COFD[10007] = 3.56695151E-03;
    COFD[10008] = -2.16246599E+01;
    COFD[10009] = 5.18371903E+00;
    COFD[10010] = -4.40201151E-01;
    COFD[10011] = 1.83498534E-02;
    COFD[10012] = -2.31435861E+01;
    COFD[10013] = 5.41522402E+00;
    COFD[10014] = -4.41568977E-01;
    COFD[10015] = 1.73190088E-02;
    COFD[10016] = -2.23095159E+01;
    COFD[10017] = 5.45157330E+00;
    COFD[10018] = -4.65253665E-01;
    COFD[10019] = 1.90612411E-02;
    COFD[10020] = -2.23193198E+01;
    COFD[10021] = 5.44963079E+00;
    COFD[10022] = -4.65238095E-01;
    COFD[10023] = 1.90696359E-02;
    COFD[10024] = -2.18845663E+01;
    COFD[10025] = 5.25324658E+00;
    COFD[10026] = -4.47318036E-01;
    COFD[10027] = 1.85841342E-02;
    COFD[10028] = -2.18970448E+01;
    COFD[10029] = 5.25324658E+00;
    COFD[10030] = -4.47318036E-01;
    COFD[10031] = 1.85841342E-02;
    COFD[10032] = -2.06803913E+01;
    COFD[10033] = 4.19360028E+00;
    COFD[10034] = -2.47975502E-01;
    COFD[10035] = 7.63834548E-03;
    COFD[10036] = -2.30861220E+01;
    COFD[10037] = 5.41058317E+00;
    COFD[10038] = -4.40692244E-01;
    COFD[10039] = 1.72706793E-02;
    COFD[10040] = -2.29301063E+01;
    COFD[10041] = 5.43667678E+00;
    COFD[10042] = -4.45729487E-01;
    COFD[10043] = 1.75510667E-02;
    COFD[10044] = -2.30716889E+01;
    COFD[10045] = 5.41058317E+00;
    COFD[10046] = -4.40692244E-01;
    COFD[10047] = 1.72706793E-02;
    COFD[10048] = -2.28184844E+01;
    COFD[10049] = 5.35124813E+00;
    COFD[10050] = -4.29932425E-01;
    COFD[10051] = 1.66897210E-02;
    COFD[10052] = -2.28355542E+01;
    COFD[10053] = 5.35124813E+00;
    COFD[10054] = -4.29932425E-01;
    COFD[10055] = 1.66897210E-02;
    COFD[10056] = -2.09076816E+01;
    COFD[10057] = 4.28688746E+00;
    COFD[10058] = -2.61952393E-01;
    COFD[10059] = 8.31175423E-03;
    COFD[10060] = -2.15840392E+01;
    COFD[10061] = 4.53462922E+00;
    COFD[10062] = -2.99463362E-01;
    COFD[10063] = 1.01311942E-02;
    COFD[10064] = -2.25529616E+01;
    COFD[10065] = 5.46297654E+00;
    COFD[10066] = -4.65485734E-01;
    COFD[10067] = 1.90234280E-02;
    COFD[10068] = -2.16025724E+01;
    COFD[10069] = 4.53462922E+00;
    COFD[10070] = -2.99463362E-01;
    COFD[10071] = 1.01311942E-02;
    COFD[10072] = -2.27363668E+01;
    COFD[10073] = 5.11476665E+00;
    COFD[10074] = -3.90017950E-01;
    COFD[10075] = 1.46174143E-02;
    COFD[10076] = -2.27363668E+01;
    COFD[10077] = 5.11476665E+00;
    COFD[10078] = -3.90017950E-01;
    COFD[10079] = 1.46174143E-02;
    COFD[10080] = -2.28367149E+01;
    COFD[10081] = 5.18065235E+00;
    COFD[10082] = -4.01095736E-01;
    COFD[10083] = 1.51910791E-02;
    COFD[10084] = -2.26349181E+01;
    COFD[10085] = 4.95583598E+00;
    COFD[10086] = -3.64620833E-01;
    COFD[10087] = 1.33399033E-02;
    COFD[10088] = -2.26506465E+01;
    COFD[10089] = 4.96659745E+00;
    COFD[10090] = -3.66306172E-01;
    COFD[10091] = 1.34235902E-02;
    COFD[10092] = -2.26574071E+01;
    COFD[10093] = 4.96659745E+00;
    COFD[10094] = -3.66306172E-01;
    COFD[10095] = 1.34235902E-02;
    COFD[10096] = -2.27912671E+01;
    COFD[10097] = 4.98252735E+00;
    COFD[10098] = -3.68805402E-01;
    COFD[10099] = 1.35477963E-02;
    COFD[10100] = -2.18595544E+01;
    COFD[10101] = 4.53718327E+00;
    COFD[10102] = -2.99856359E-01;
    COFD[10103] = 1.01504725E-02;
    COFD[10104] = -2.14584112E+01;
    COFD[10105] = 4.38117808E+00;
    COFD[10106] = -2.78570982E-01;
    COFD[10107] = 9.19430194E-03;
    COFD[10108] = -2.24440124E+01;
    COFD[10109] = 4.74813432E+00;
    COFD[10110] = -3.32238042E-01;
    COFD[10111] = 1.17370757E-02;
    COFD[10112] = -2.25319271E+01;
    COFD[10113] = 4.79201228E+00;
    COFD[10114] = -3.39345465E-01;
    COFD[10115] = 1.20970888E-02;
    COFD[10116] = -2.19159785E+01;
    COFD[10117] = 4.51065350E+00;
    COFD[10118] = -2.95784060E-01;
    COFD[10119] = 9.95114964E-03;
    COFD[10120] = -2.10625949E+01;
    COFD[10121] = 4.28688746E+00;
    COFD[10122] = -2.61952393E-01;
    COFD[10123] = 8.31175423E-03;
    COFD[10124] = -2.10708054E+01;
    COFD[10125] = 4.28688746E+00;
    COFD[10126] = -2.61952393E-01;
    COFD[10127] = 8.31175423E-03;
    COFD[10128] = -2.16429075E+01;
    COFD[10129] = 4.49834811E+00;
    COFD[10130] = -2.93898563E-01;
    COFD[10131] = 9.85901966E-03;
    COFD[10132] = -2.15585553E+01;
    COFD[10133] = 4.41358685E+00;
    COFD[10134] = -2.83552193E-01;
    COFD[10135] = 9.43843173E-03;
    COFD[10136] = -2.27937067E+01;
    COFD[10137] = 5.15050798E+00;
    COFD[10138] = -3.95966829E-01;
    COFD[10139] = 1.49238043E-02;
    COFD[10140] = -2.27261124E+01;
    COFD[10141] = 5.11476665E+00;
    COFD[10142] = -3.90017950E-01;
    COFD[10143] = 1.46174143E-02;
    COFD[10144] = -2.13456084E+01;
    COFD[10145] = 4.25479441E+00;
    COFD[10146] = -2.59093497E-01;
    COFD[10147] = 8.23869170E-03;
    COFD[10148] = -2.14606335E+01;
    COFD[10149] = 4.28507327E+00;
    COFD[10150] = -2.63941348E-01;
    COFD[10151] = 8.48270771E-03;
    COFD[10152] = -2.16496748E+01;
    COFD[10153] = 4.49834811E+00;
    COFD[10154] = -2.93898563E-01;
    COFD[10155] = 9.85901966E-03;
    COFD[10156] = -2.14478435E+01;
    COFD[10157] = 4.21514791E+00;
    COFD[10158] = -2.52791771E-01;
    COFD[10159] = 7.92322551E-03;
    COFD[10160] = -2.14832299E+01;
    COFD[10161] = 4.24469152E+00;
    COFD[10162] = -2.55617920E-01;
    COFD[10163] = 8.00605387E-03;
    COFD[10164] = -2.20144983E+01;
    COFD[10165] = 4.41014015E+00;
    COFD[10166] = -2.80504988E-01;
    COFD[10167] = 9.20806328E-03;
    COFD[10168] = -2.20113604E+01;
    COFD[10169] = 4.41014015E+00;
    COFD[10170] = -2.80504988E-01;
    COFD[10171] = 9.20806328E-03;
    COFD[10172] = -2.07325259E+01;
    COFD[10173] = 3.76716506E+00;
    COFD[10174] = -1.85775467E-01;
    COFD[10175] = 4.69798135E-03;
    COFD[10176] = -2.04162071E+01;
    COFD[10177] = 3.53236186E+00;
    COFD[10178] = -1.51329616E-01;
    COFD[10179] = 3.06303488E-03;
    COFD[10180] = -2.14573505E+01;
    COFD[10181] = 4.11830569E+00;
    COFD[10182] = -2.36735638E-01;
    COFD[10183] = 7.09853565E-03;
    COFD[10184] = -2.05116161E+01;
    COFD[10185] = 3.63890356E+00;
    COFD[10186] = -1.67076072E-01;
    COFD[10187] = 3.81485236E-03;
    COFD[10188] = -2.16048217E+01;
    COFD[10189] = 5.17969892E+00;
    COFD[10190] = -4.39844545E-01;
    COFD[10191] = 1.83408439E-02;
    COFD[10192] = -2.01926972E+01;
    COFD[10193] = 5.40010667E+00;
    COFD[10194] = -4.63239430E-01;
    COFD[10195] = 1.91539217E-02;
    COFD[10196] = -1.62247094E+01;
    COFD[10197] = 3.85939026E+00;
    COFD[10198] = -2.83353647E-01;
    COFD[10199] = 1.21442775E-02;
    COFD[10200] = -1.97866869E+01;
    COFD[10201] = 4.78902303E+00;
    COFD[10202] = -3.96144974E-01;
    COFD[10203] = 1.67211738E-02;
    COFD[10204] = -2.12048425E+01;
    COFD[10205] = 5.11479735E+00;
    COFD[10206] = -4.33050578E-01;
    COFD[10207] = 1.81098360E-02;
    COFD[10208] = -1.98133783E+01;
    COFD[10209] = 4.78902303E+00;
    COFD[10210] = -3.96144974E-01;
    COFD[10211] = 1.67211738E-02;
    COFD[10212] = -2.01714866E+01;
    COFD[10213] = 4.25399524E+00;
    COFD[10214] = -2.57013670E-01;
    COFD[10215] = 8.07339307E-03;
    COFD[10216] = -2.09238627E+01;
    COFD[10217] = 5.02943557E+00;
    COFD[10218] = -4.23811924E-01;
    COFD[10219] = 1.77810575E-02;
    COFD[10220] = -2.29289604E+01;
    COFD[10221] = 5.48025735E+00;
    COFD[10222] = -4.56422639E-01;
    COFD[10223] = 1.82014258E-02;
    COFD[10224] = -2.18504787E+01;
    COFD[10225] = 5.39546703E+00;
    COFD[10226] = -4.62825371E-01;
    COFD[10227] = 1.91432238E-02;
    COFD[10228] = -2.18110440E+01;
    COFD[10229] = 5.37758536E+00;
    COFD[10230] = -4.60945701E-01;
    COFD[10231] = 1.90787306E-02;
    COFD[10232] = -2.12169140E+01;
    COFD[10233] = 5.11479735E+00;
    COFD[10234] = -4.33050578E-01;
    COFD[10235] = 1.81098360E-02;
    COFD[10236] = -2.12285424E+01;
    COFD[10237] = 5.11479735E+00;
    COFD[10238] = -4.33050578E-01;
    COFD[10239] = 1.81098360E-02;
    COFD[10240] = -2.12861576E+01;
    COFD[10241] = 4.59638767E+00;
    COFD[10242] = -3.08980273E-01;
    COFD[10243] = 1.05986522E-02;
    COFD[10244] = -2.29039602E+01;
    COFD[10245] = 5.47759909E+00;
    COFD[10246] = -4.55570838E-01;
    COFD[10247] = 1.81456235E-02;
    COFD[10248] = -2.26898752E+01;
    COFD[10249] = 5.48310931E+00;
    COFD[10250] = -4.57615590E-01;
    COFD[10251] = 1.82838764E-02;
    COFD[10252] = -2.28904184E+01;
    COFD[10253] = 5.47759909E+00;
    COFD[10254] = -4.55570838E-01;
    COFD[10255] = 1.81456235E-02;
    COFD[10256] = -2.27210131E+01;
    COFD[10257] = 5.45828076E+00;
    COFD[10258] = -4.50591566E-01;
    COFD[10259] = 1.78383531E-02;
    COFD[10260] = -2.27371462E+01;
    COFD[10261] = 5.45828076E+00;
    COFD[10262] = -4.50591566E-01;
    COFD[10263] = 1.78383531E-02;
    COFD[10264] = -2.14610939E+01;
    COFD[10265] = 4.66858934E+00;
    COFD[10266] = -3.20108476E-01;
    COFD[10267] = 1.11452199E-02;
    COFD[10268] = -2.20529119E+01;
    COFD[10269] = 4.88259016E+00;
    COFD[10270] = -3.53111373E-01;
    COFD[10271] = 1.27672516E-02;
    COFD[10272] = -2.20783474E+01;
    COFD[10273] = 5.42074110E+00;
    COFD[10274] = -4.64982008E-01;
    COFD[10275] = 1.91934150E-02;
    COFD[10276] = -2.20699218E+01;
    COFD[10277] = 4.88259016E+00;
    COFD[10278] = -3.53111373E-01;
    COFD[10279] = 1.27672516E-02;
    COFD[10280] = -2.28434554E+01;
    COFD[10281] = 5.30874092E+00;
    COFD[10282] = -4.22570517E-01;
    COFD[10283] = 1.63018518E-02;
    COFD[10284] = -2.28434554E+01;
    COFD[10285] = 5.30874092E+00;
    COFD[10286] = -4.22570517E-01;
    COFD[10287] = 1.63018518E-02;
    COFD[10288] = -2.29126248E+01;
    COFD[10289] = 5.36245215E+00;
    COFD[10290] = -4.31910941E-01;
    COFD[10291] = 1.67950232E-02;
    COFD[10292] = -2.28440138E+01;
    COFD[10293] = 5.19534599E+00;
    COFD[10294] = -4.03564468E-01;
    COFD[10295] = 1.53188143E-02;
    COFD[10296] = -2.28532517E+01;
    COFD[10297] = 5.20408596E+00;
    COFD[10298] = -4.04973723E-01;
    COFD[10299] = 1.53901363E-02;
    COFD[10300] = -2.28593475E+01;
    COFD[10301] = 5.20408596E+00;
    COFD[10302] = -4.04973723E-01;
    COFD[10303] = 1.53901363E-02;
    COFD[10304] = -2.29938969E+01;
    COFD[10305] = 5.21722096E+00;
    COFD[10306] = -4.07093564E-01;
    COFD[10307] = 1.54974812E-02;
    COFD[10308] = -2.23303179E+01;
    COFD[10309] = 4.88479479E+00;
    COFD[10310] = -3.53456980E-01;
    COFD[10311] = 1.27844233E-02;
    COFD[10312] = -2.23028472E+01;
    COFD[10313] = 4.88611782E+00;
    COFD[10314] = -3.53664413E-01;
    COFD[10315] = 1.27947308E-02;
    COFD[10316] = -2.28166735E+01;
    COFD[10317] = 5.05461151E+00;
    COFD[10318] = -3.80222586E-01;
    COFD[10319] = 1.41187082E-02;
    COFD[10320] = -2.29010051E+01;
    COFD[10321] = 5.09658829E+00;
    COFD[10322] = -3.87036802E-01;
    COFD[10323] = 1.44650532E-02;
    COFD[10324] = -2.23893523E+01;
    COFD[10325] = 4.86181814E+00;
    COFD[10326] = -3.49857477E-01;
    COFD[10327] = 1.26056712E-02;
    COFD[10328] = -2.16041988E+01;
    COFD[10329] = 4.66858934E+00;
    COFD[10330] = -3.20108476E-01;
    COFD[10331] = 1.11452199E-02;
    COFD[10332] = -2.16116846E+01;
    COFD[10333] = 4.66858934E+00;
    COFD[10334] = -3.20108476E-01;
    COFD[10335] = 1.11452199E-02;
    COFD[10336] = -2.21166468E+01;
    COFD[10337] = 4.85109413E+00;
    COFD[10338] = -3.48178913E-01;
    COFD[10339] = 1.25223681E-02;
    COFD[10340] = -2.23943879E+01;
    COFD[10341] = 4.91274983E+00;
    COFD[10342] = -3.57867526E-01;
    COFD[10343] = 1.30043584E-02;
    COFD[10344] = -2.28830519E+01;
    COFD[10345] = 5.33738936E+00;
    COFD[10346] = -4.27502526E-01;
    COFD[10347] = 1.65608197E-02;
    COFD[10348] = -2.28339934E+01;
    COFD[10349] = 5.30874092E+00;
    COFD[10350] = -4.22570517E-01;
    COFD[10351] = 1.63018518E-02;
    COFD[10352] = -2.21833064E+01;
    COFD[10353] = 4.75657673E+00;
    COFD[10354] = -3.33531354E-01;
    COFD[10355] = 1.18003758E-02;
    COFD[10356] = -2.23442480E+01;
    COFD[10357] = 4.80471822E+00;
    COFD[10358] = -3.40939669E-01;
    COFD[10359] = 1.21639062E-02;
    COFD[10360] = -2.21227489E+01;
    COFD[10361] = 4.85109413E+00;
    COFD[10362] = -3.48178913E-01;
    COFD[10363] = 1.25223681E-02;
    COFD[10364] = -2.22319277E+01;
    COFD[10365] = 4.69258044E+00;
    COFD[10366] = -3.23762098E-01;
    COFD[10367] = 1.13233362E-02;
    COFD[10368] = -2.20449981E+01;
    COFD[10369] = 4.63376397E+00;
    COFD[10370] = -3.14793416E-01;
    COFD[10371] = 1.08858273E-02;
    COFD[10372] = -2.25315042E+01;
    COFD[10373] = 4.77923223E+00;
    COFD[10374] = -3.37003339E-01;
    COFD[10375] = 1.19703492E-02;
    COFD[10376] = -2.25288055E+01;
    COFD[10377] = 4.77923223E+00;
    COFD[10378] = -3.37003339E-01;
    COFD[10379] = 1.19703492E-02;
    COFD[10380] = -2.16640062E+01;
    COFD[10381] = 4.30841644E+00;
    COFD[10382] = -2.65200479E-01;
    COFD[10383] = 8.46889236E-03;
    COFD[10384] = -2.14573505E+01;
    COFD[10385] = 4.11830569E+00;
    COFD[10386] = -2.36735638E-01;
    COFD[10387] = 7.09853565E-03;
    COFD[10388] = -2.20903525E+01;
    COFD[10389] = 4.53657991E+00;
    COFD[10390] = -2.99763510E-01;
    COFD[10391] = 1.01459176E-02;
    COFD[10392] = -2.15157070E+01;
    COFD[10393] = 4.21035526E+00;
    COFD[10394] = -2.50479118E-01;
    COFD[10395] = 7.75869538E-03;
    COFD[10396] = -2.08959939E+01;
    COFD[10397] = 5.02187007E+00;
    COFD[10398] = -4.22927894E-01;
    COFD[10399] = 1.77465747E-02;
    COFD[10400] = -2.04559021E+01;
    COFD[10401] = 5.44558859E+00;
    COFD[10402] = -4.65276039E-01;
    COFD[10403] = 1.90929031E-02;
    COFD[10404] = -1.66789079E+01;
    COFD[10405] = 4.00162523E+00;
    COFD[10406] = -3.00504266E-01;
    COFD[10407] = 1.28349600E-02;
    COFD[10408] = -2.02611614E+01;
    COFD[10409] = 4.92728595E+00;
    COFD[10410] = -4.11813748E-01;
    COFD[10411] = 1.73103280E-02;
    COFD[10412] = -2.16612294E+01;
    COFD[10413] = 5.23259254E+00;
    COFD[10414] = -4.44988191E-01;
    COFD[10415] = 1.84967231E-02;
    COFD[10416] = -2.02886114E+01;
    COFD[10417] = 4.92728595E+00;
    COFD[10418] = -4.11813748E-01;
    COFD[10419] = 1.73103280E-02;
    COFD[10420] = -1.88216419E+01;
    COFD[10421] = 3.61672004E+00;
    COFD[10422] = -1.65415347E-01;
    COFD[10423] = 3.79600441E-03;
    COFD[10424] = -2.14243807E+01;
    COFD[10425] = 5.16859564E+00;
    COFD[10426] = -4.38877619E-01;
    COFD[10427] = 1.83174566E-02;
    COFD[10428] = -2.30145662E+01;
    COFD[10429] = 5.43378566E+00;
    COFD[10430] = -4.45159774E-01;
    COFD[10431] = 1.75190628E-02;
    COFD[10432] = -2.21249716E+01;
    COFD[10433] = 5.44426165E+00;
    COFD[10434] = -4.65323191E-01;
    COFD[10435] = 1.91034051E-02;
    COFD[10436] = -2.21451593E+01;
    COFD[10437] = 5.44423992E+00;
    COFD[10438] = -4.65323901E-01;
    COFD[10439] = 1.91035720E-02;
    COFD[10440] = -2.16082141E+01;
    COFD[10441] = 5.21099602E+00;
    COFD[10442] = -4.42524031E-01;
    COFD[10443] = 1.84030593E-02;
    COFD[10444] = -2.16204501E+01;
    COFD[10445] = 5.21099602E+00;
    COFD[10446] = -4.42524031E-01;
    COFD[10447] = 1.84030593E-02;
    COFD[10448] = -2.07241021E+01;
    COFD[10449] = 4.28574979E+00;
    COFD[10450] = -2.61781425E-01;
    COFD[10451] = 8.30349840E-03;
    COFD[10452] = -2.29719722E+01;
    COFD[10453] = 5.43353924E+00;
    COFD[10454] = -4.45111358E-01;
    COFD[10455] = 1.75163465E-02;
    COFD[10456] = -2.27854551E+01;
    COFD[10457] = 5.45051999E+00;
    COFD[10458] = -4.48822771E-01;
    COFD[10459] = 1.77335774E-02;
    COFD[10460] = -2.29577921E+01;
    COFD[10461] = 5.43353924E+00;
    COFD[10462] = -4.45111358E-01;
    COFD[10463] = 1.75163465E-02;
    COFD[10464] = -2.27171217E+01;
    COFD[10465] = 5.38171422E+00;
    COFD[10466] = -4.35323468E-01;
    COFD[10467] = 1.69770246E-02;
    COFD[10468] = -2.27339269E+01;
    COFD[10469] = 5.38171422E+00;
    COFD[10470] = -4.35323468E-01;
    COFD[10471] = 1.69770246E-02;
    COFD[10472] = -2.09356172E+01;
    COFD[10473] = 4.37219452E+00;
    COFD[10474] = -2.74784562E-01;
    COFD[10475] = 8.93131387E-03;
    COFD[10476] = -2.15822613E+01;
    COFD[10477] = 4.60749458E+00;
    COFD[10478] = -3.10696534E-01;
    COFD[10479] = 1.06831054E-02;
    COFD[10480] = -2.23385572E+01;
    COFD[10481] = 5.45382388E+00;
    COFD[10482] = -4.65297200E-01;
    COFD[10483] = 1.90536055E-02;
    COFD[10484] = -2.16003553E+01;
    COFD[10485] = 4.60749458E+00;
    COFD[10486] = -3.10696534E-01;
    COFD[10487] = 1.06831054E-02;
    COFD[10488] = -2.26645375E+01;
    COFD[10489] = 5.15653775E+00;
    COFD[10490] = -3.96988537E-01;
    COFD[10491] = 1.49769158E-02;
    COFD[10492] = -2.26645375E+01;
    COFD[10493] = 5.15653775E+00;
    COFD[10494] = -3.96988537E-01;
    COFD[10495] = 1.49769158E-02;
    COFD[10496] = -2.27726207E+01;
    COFD[10497] = 5.22626014E+00;
    COFD[10498] = -4.08564584E-01;
    COFD[10499] = 1.55722943E-02;
    COFD[10500] = -2.26098548E+01;
    COFD[10501] = 5.01607052E+00;
    COFD[10502] = -3.74092059E-01;
    COFD[10503] = 1.38112395E-02;
    COFD[10504] = -2.26226017E+01;
    COFD[10505] = 5.02598120E+00;
    COFD[10506] = -3.75659093E-01;
    COFD[10507] = 1.38894769E-02;
    COFD[10508] = -2.26291682E+01;
    COFD[10509] = 5.02598120E+00;
    COFD[10510] = -3.75659093E-01;
    COFD[10511] = 1.38894769E-02;
    COFD[10512] = -2.27668090E+01;
    COFD[10513] = 5.04103112E+00;
    COFD[10514] = -3.78039122E-01;
    COFD[10515] = 1.40083345E-02;
    COFD[10516] = -2.18647351E+01;
    COFD[10517] = 4.61001619E+00;
    COFD[10518] = -3.11086406E-01;
    COFD[10519] = 1.07022976E-02;
    COFD[10520] = -2.13113937E+01;
    COFD[10521] = 4.38760834E+00;
    COFD[10522] = -2.80341923E-01;
    COFD[10523] = 9.30748203E-03;
    COFD[10524] = -2.24468874E+01;
    COFD[10525] = 4.81788504E+00;
    COFD[10526] = -3.42990993E-01;
    COFD[10527] = 1.22653559E-02;
    COFD[10528] = -2.25039575E+01;
    COFD[10529] = 4.84915044E+00;
    COFD[10530] = -3.48402726E-01;
    COFD[10531] = 1.25499430E-02;
    COFD[10532] = -2.19212431E+01;
    COFD[10533] = 4.58361215E+00;
    COFD[10534] = -3.07008194E-01;
    COFD[10535] = 1.05016776E-02;
    COFD[10536] = -2.10871397E+01;
    COFD[10537] = 4.37219452E+00;
    COFD[10538] = -2.74784562E-01;
    COFD[10539] = 8.93131387E-03;
    COFD[10540] = -2.10951403E+01;
    COFD[10541] = 4.37219452E+00;
    COFD[10542] = -2.74784562E-01;
    COFD[10543] = 8.93131387E-03;
    COFD[10544] = -2.16405632E+01;
    COFD[10545] = 4.57106634E+00;
    COFD[10546] = -3.05074221E-01;
    COFD[10547] = 1.04066420E-02;
    COFD[10548] = -2.14143742E+01;
    COFD[10549] = 4.41973109E+00;
    COFD[10550] = -2.85282560E-01;
    COFD[10551] = 9.54967542E-03;
    COFD[10552] = -2.27063862E+01;
    COFD[10553] = 5.18619255E+00;
    COFD[10554] = -4.02090863E-01;
    COFD[10555] = 1.52443533E-02;
    COFD[10556] = -2.26545108E+01;
    COFD[10557] = 5.15653775E+00;
    COFD[10558] = -3.96988537E-01;
    COFD[10559] = 1.49769158E-02;
    COFD[10560] = -2.12434067E+01;
    COFD[10561] = 4.28015384E+00;
    COFD[10562] = -2.63514048E-01;
    COFD[10563] = 8.47277117E-03;
    COFD[10564] = -2.13309187E+01;
    COFD[10565] = 4.29619647E+00;
    COFD[10566] = -2.66307880E-01;
    COFD[10567] = 8.62039924E-03;
    COFD[10568] = -2.16471363E+01;
    COFD[10569] = 4.57106634E+00;
    COFD[10570] = -3.05074221E-01;
    COFD[10571] = 1.04066420E-02;
    COFD[10572] = -2.13912488E+01;
    COFD[10573] = 4.25856149E+00;
    COFD[10574] = -2.59776530E-01;
    COFD[10575] = 8.27569932E-03;
    COFD[10576] = -2.15271080E+01;
    COFD[10577] = 4.33297139E+00;
    COFD[10578] = -2.68886993E-01;
    COFD[10579] = 8.64665258E-03;
    COFD[10580] = -2.20447991E+01;
    COFD[10581] = 4.49024164E+00;
    COFD[10582] = -2.92656701E-01;
    COFD[10583] = 9.79834494E-03;
    COFD[10584] = -2.20417938E+01;
    COFD[10585] = 4.49024164E+00;
    COFD[10586] = -2.92656701E-01;
    COFD[10587] = 9.79834494E-03;
    COFD[10588] = -2.07979763E+01;
    COFD[10589] = 3.86302522E+00;
    COFD[10590] = -2.00078966E-01;
    COFD[10591] = 5.38570099E-03;
    COFD[10592] = -2.05116161E+01;
    COFD[10593] = 3.63890356E+00;
    COFD[10594] = -1.67076072E-01;
    COFD[10595] = 3.81485236E-03;
    COFD[10596] = -2.15157070E+01;
    COFD[10597] = 4.21035526E+00;
    COFD[10598] = -2.50479118E-01;
    COFD[10599] = 7.75869538E-03;
    COFD[10600] = -2.05695271E+01;
    COFD[10601] = 3.73118862E+00;
    COFD[10602] = -1.80835826E-01;
    COFD[10603] = 4.47560913E-03;
    COFD[10604] = -2.13977914E+01;
    COFD[10605] = 5.16203659E+00;
    COFD[10606] = -4.38179184E-01;
    COFD[10607] = 1.82931217E-02;
    COFD[10608] = -1.38221692E+01;
    COFD[10609] = 3.46560000E+00;
    COFD[10610] = -2.34075114E-01;
    COFD[10611] = 1.00849661E-02;
    COFD[10612] = -1.11071542E+01;
    COFD[10613] = 2.21467976E+00;
    COFD[10614] = -7.20501534E-02;
    COFD[10615] = 3.07515437E-03;
    COFD[10616] = -1.36761960E+01;
    COFD[10617] = 2.89678405E+00;
    COFD[10618] = -1.62418864E-01;
    COFD[10619] = 7.07300755E-03;
    COFD[10620] = -1.48133672E+01;
    COFD[10621] = 3.16773813E+00;
    COFD[10622] = -1.96823168E-01;
    COFD[10623] = 8.53058072E-03;
    COFD[10624] = -1.36954221E+01;
    COFD[10625] = 2.89678405E+00;
    COFD[10626] = -1.62418864E-01;
    COFD[10627] = 7.07300755E-03;
    COFD[10628] = -2.06077924E+01;
    COFD[10629] = 5.33386529E+00;
    COFD[10630] = -4.56343427E-01;
    COFD[10631] = 1.89202788E-02;
    COFD[10632] = -1.45806444E+01;
    COFD[10633] = 3.07340947E+00;
    COFD[10634] = -1.84805140E-01;
    COFD[10635] = 8.01973216E-03;
    COFD[10636] = -1.74370371E+01;
    COFD[10637] = 4.06188586E+00;
    COFD[10638] = -3.07715468E-01;
    COFD[10639] = 1.31226885E-02;
    COFD[10640] = -1.54972198E+01;
    COFD[10641] = 3.45957322E+00;
    COFD[10642] = -2.33335381E-01;
    COFD[10643] = 1.00546680E-02;
    COFD[10644] = -1.54573340E+01;
    COFD[10645] = 3.44303981E+00;
    COFD[10646] = -2.31296721E-01;
    COFD[10647] = 9.97079557E-03;
    COFD[10648] = -1.48205460E+01;
    COFD[10649] = 3.16773813E+00;
    COFD[10650] = -1.96823168E-01;
    COFD[10651] = 8.53058072E-03;
    COFD[10652] = -1.48273948E+01;
    COFD[10653] = 3.16773813E+00;
    COFD[10654] = -1.96823168E-01;
    COFD[10655] = 8.53058072E-03;
    COFD[10656] = -1.99350648E+01;
    COFD[10657] = 4.98617709E+00;
    COFD[10658] = -4.18736053E-01;
    COFD[10659] = 1.75821344E-02;
    COFD[10660] = -1.75535683E+01;
    COFD[10661] = 4.07809359E+00;
    COFD[10662] = -3.09646036E-01;
    COFD[10663] = 1.31992743E-02;
    COFD[10664] = -1.71876409E+01;
    COFD[10665] = 4.03536932E+00;
    COFD[10666] = -3.04548174E-01;
    COFD[10667] = 1.29965591E-02;
    COFD[10668] = -1.75452737E+01;
    COFD[10669] = 4.07809359E+00;
    COFD[10670] = -3.09646036E-01;
    COFD[10671] = 1.31992743E-02;
    COFD[10672] = -1.75930902E+01;
    COFD[10673] = 4.16734978E+00;
    COFD[10674] = -3.20508152E-01;
    COFD[10675] = 1.36405589E-02;
    COFD[10676] = -1.76034182E+01;
    COFD[10677] = 4.16734978E+00;
    COFD[10678] = -3.20508152E-01;
    COFD[10679] = 1.36405589E-02;
    COFD[10680] = -1.98307830E+01;
    COFD[10681] = 4.94145734E+00;
    COFD[10682] = -4.13478429E-01;
    COFD[10683] = 1.73756310E-02;
    COFD[10684] = -1.96541192E+01;
    COFD[10685] = 4.83228605E+00;
    COFD[10686] = -4.01127235E-01;
    COFD[10687] = 1.69121265E-02;
    COFD[10688] = -1.55055897E+01;
    COFD[10689] = 3.50663983E+00;
    COFD[10690] = -2.39273841E-01;
    COFD[10691] = 1.03049488E-02;
    COFD[10692] = -1.96633539E+01;
    COFD[10693] = 4.83228605E+00;
    COFD[10694] = -4.01127235E-01;
    COFD[10695] = 1.69121265E-02;
    COFD[10696] = -1.86480131E+01;
    COFD[10697] = 4.45508831E+00;
    COFD[10698] = -3.56068927E-01;
    COFD[10699] = 1.51117769E-02;
    COFD[10700] = -1.86480131E+01;
    COFD[10701] = 4.45508831E+00;
    COFD[10702] = -3.56068927E-01;
    COFD[10703] = 1.51117769E-02;
    COFD[10704] = -1.84298463E+01;
    COFD[10705] = 4.38777139E+00;
    COFD[10706] = -3.47895358E-01;
    COFD[10707] = 1.47793796E-02;
    COFD[10708] = -1.91111450E+01;
    COFD[10709] = 4.55863849E+00;
    COFD[10710] = -3.68274023E-01;
    COFD[10711] = 1.55925459E-02;
    COFD[10712] = -1.90794464E+01;
    COFD[10713] = 4.55246544E+00;
    COFD[10714] = -3.67548311E-01;
    COFD[10715] = 1.55640443E-02;
    COFD[10716] = -1.90824831E+01;
    COFD[10717] = 4.55246544E+00;
    COFD[10718] = -3.67548311E-01;
    COFD[10719] = 1.55640443E-02;
    COFD[10720] = -1.91962804E+01;
    COFD[10721] = 4.54321736E+00;
    COFD[10722] = -3.66462203E-01;
    COFD[10723] = 1.55214344E-02;
    COFD[10724] = -1.99275800E+01;
    COFD[10725] = 4.83095096E+00;
    COFD[10726] = -4.00972558E-01;
    COFD[10727] = 1.69061563E-02;
    COFD[10728] = -2.01873367E+01;
    COFD[10729] = 4.92914806E+00;
    COFD[10730] = -4.12032782E-01;
    COFD[10731] = 1.73189335E-02;
    COFD[10732] = -1.97684750E+01;
    COFD[10733] = 4.71961665E+00;
    COFD[10734] = -3.87953877E-01;
    COFD[10735] = 1.63982506E-02;
    COFD[10736] = -1.96740282E+01;
    COFD[10737] = 4.68296253E+00;
    COFD[10738] = -3.83520383E-01;
    COFD[10739] = 1.62187256E-02;
    COFD[10740] = -2.00440183E+01;
    COFD[10741] = 4.84467881E+00;
    COFD[10742] = -4.02557342E-01;
    COFD[10743] = 1.69670693E-02;
    COFD[10744] = -1.99112647E+01;
    COFD[10745] = 4.94145734E+00;
    COFD[10746] = -4.13478429E-01;
    COFD[10747] = 1.73756310E-02;
    COFD[10748] = -1.99151960E+01;
    COFD[10749] = 4.94145734E+00;
    COFD[10750] = -4.13478429E-01;
    COFD[10751] = 1.73756310E-02;
    COFD[10752] = -1.97977739E+01;
    COFD[10753] = 4.85079945E+00;
    COFD[10754] = -4.03257284E-01;
    COFD[10755] = 1.69936877E-02;
    COFD[10756] = -2.01834870E+01;
    COFD[10757] = 4.90590035E+00;
    COFD[10758] = -4.09326358E-01;
    COFD[10759] = 1.72138685E-02;
    COFD[10760] = -1.85439247E+01;
    COFD[10761] = 4.42312522E+00;
    COFD[10762] = -3.52239033E-01;
    COFD[10763] = 1.49581799E-02;
    COFD[10764] = -1.86427314E+01;
    COFD[10765] = 4.45508831E+00;
    COFD[10766] = -3.56068927E-01;
    COFD[10767] = 1.51117769E-02;
    COFD[10768] = -2.04535647E+01;
    COFD[10769] = 4.98209169E+00;
    COFD[10770] = -4.18253583E-01;
    COFD[10771] = 1.75630850E-02;
    COFD[10772] = -2.05335806E+01;
    COFD[10773] = 4.97826004E+00;
    COFD[10774] = -4.17803418E-01;
    COFD[10775] = 1.75454151E-02;
    COFD[10776] = -1.98008145E+01;
    COFD[10777] = 4.85079945E+00;
    COFD[10778] = -4.03257284E-01;
    COFD[10779] = 1.69936877E-02;
    COFD[10780] = -2.06399268E+01;
    COFD[10781] = 4.98336829E+00;
    COFD[10782] = -4.18404333E-01;
    COFD[10783] = 1.75690366E-02;
    COFD[10784] = -2.04878867E+01;
    COFD[10785] = 4.96174500E+00;
    COFD[10786] = -4.15861207E-01;
    COFD[10787] = 1.74690994E-02;
    COFD[10788] = -2.04827540E+01;
    COFD[10789] = 4.88918294E+00;
    COFD[10790] = -4.07583331E-01;
    COFD[10791] = 1.71552617E-02;
    COFD[10792] = -2.04816453E+01;
    COFD[10793] = 4.88918294E+00;
    COFD[10794] = -4.07583331E-01;
    COFD[10795] = 1.71552617E-02;
    COFD[10796] = -2.12317255E+01;
    COFD[10797] = 5.12393233E+00;
    COFD[10798] = -4.34047487E-01;
    COFD[10799] = 1.81456921E-02;
    COFD[10800] = -2.16048217E+01;
    COFD[10801] = 5.17969892E+00;
    COFD[10802] = -4.39844545E-01;
    COFD[10803] = 1.83408439E-02;
    COFD[10804] = -2.08959939E+01;
    COFD[10805] = 5.02187007E+00;
    COFD[10806] = -4.22927894E-01;
    COFD[10807] = 1.77465747E-02;
    COFD[10808] = -2.13977914E+01;
    COFD[10809] = 5.16203659E+00;
    COFD[10810] = -4.38179184E-01;
    COFD[10811] = 1.82931217E-02;
    COFD[10812] = -1.45599006E+01;
    COFD[10813] = 3.06880733E+00;
    COFD[10814] = -1.84232163E-01;
    COFD[10815] = 7.99594639E-03;
}


/*List of specs with small weight, dim NLITE */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetKTDIF(int* KTDIF) {
    KTDIF[0] = 0;
    KTDIF[1] = 1;
}


/*Poly fits for thermal diff ratios, dim NO*NLITE*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFTD(amrex::Real* COFTD) {
    COFTD[0] = 0.00000000E+00;
    COFTD[1] = 0.00000000E+00;
    COFTD[2] = 0.00000000E+00;
    COFTD[3] = 0.00000000E+00;
    COFTD[4] = 1.58593068E-01;
    COFTD[5] = 3.77534860E-05;
    COFTD[6] = -1.76801591E-08;
    COFTD[7] = 2.63675347E-12;
    COFTD[8] = 2.97686725E-01;
    COFTD[9] = 2.85458827E-04;
    COFTD[10] = -1.28554387E-07;
    COFTD[11] = 1.75751301E-11;
    COFTD[12] = 2.53688960E-01;
    COFTD[13] = 3.89161371E-04;
    COFTD[14] = -1.70616430E-07;
    COFTD[15] = 2.27713392E-11;
    COFTD[16] = 2.99926473E-01;
    COFTD[17] = 2.87606574E-04;
    COFTD[18] = -1.29521610E-07;
    COFTD[19] = 1.77073627E-11;
    COFTD[20] = -1.11292235E-01;
    COFTD[21] = 6.84139751E-04;
    COFTD[22] = -2.50608599E-07;
    COFTD[23] = 2.96932609E-11;
    COFTD[24] = 2.71278225E-01;
    COFTD[25] = 3.59778238E-04;
    COFTD[26] = -1.59065631E-07;
    COFTD[27] = 2.13804877E-11;
    COFTD[28] = 6.30512606E-02;
    COFTD[29] = 6.13216460E-04;
    COFTD[30] = -2.48082221E-07;
    COFTD[31] = 3.12631540E-11;
    COFTD[32] = 1.73585553E-01;
    COFTD[33] = 4.39398024E-04;
    COFTD[34] = -1.87350240E-07;
    COFTD[35] = 2.44721912E-11;
    COFTD[36] = 1.79062640E-01;
    COFTD[37] = 4.38432550E-04;
    COFTD[38] = -1.87267083E-07;
    COFTD[39] = 2.44928705E-11;
    COFTD[40] = 2.54177980E-01;
    COFTD[41] = 3.89911532E-04;
    COFTD[42] = -1.70945317E-07;
    COFTD[43] = 2.28152340E-11;
    COFTD[44] = 2.54638851E-01;
    COFTD[45] = 3.90618514E-04;
    COFTD[46] = -1.71255272E-07;
    COFTD[47] = 2.28566022E-11;
    COFTD[48] = -9.12115192E-02;
    COFTD[49] = 7.07340708E-04;
    COFTD[50] = -2.63898211E-07;
    COFTD[51] = 3.16206874E-11;
    COFTD[52] = 5.84180917E-02;
    COFTD[53] = 6.03362935E-04;
    COFTD[54] = -2.43735599E-07;
    COFTD[55] = 3.06860991E-11;
    COFTD[56] = 6.68057417E-02;
    COFTD[57] = 5.92435093E-04;
    COFTD[58] = -2.40250790E-07;
    COFTD[59] = 3.03233580E-11;
    COFTD[60] = 5.82822574E-02;
    COFTD[61] = 6.01959990E-04;
    COFTD[62] = -2.43168862E-07;
    COFTD[63] = 3.06147475E-11;
    COFTD[64] = 4.18090302E-02;
    COFTD[65] = 6.09450110E-04;
    COFTD[66] = -2.44413437E-07;
    COFTD[67] = 3.06283377E-11;
    COFTD[68] = 4.19397901E-02;
    COFTD[69] = 6.11356195E-04;
    COFTD[70] = -2.45177852E-07;
    COFTD[71] = 3.07241293E-11;
    COFTD[72] = -8.53031345E-02;
    COFTD[73] = 7.07704735E-04;
    COFTD[74] = -2.65142408E-07;
    COFTD[75] = 3.18513094E-11;
    COFTD[76] = -6.66682584E-02;
    COFTD[77] = 7.08752133E-04;
    COFTD[78] = -2.68826880E-07;
    COFTD[79] = 3.25385022E-11;
    COFTD[80] = 1.79302797E-01;
    COFTD[81] = 4.89725732E-04;
    COFTD[82] = -2.07986653E-07;
    COFTD[83] = 2.70899527E-11;
    COFTD[84] = -6.68148032E-02;
    COFTD[85] = 7.10310054E-04;
    COFTD[86] = -2.69417794E-07;
    COFTD[87] = 3.26100257E-11;
    COFTD[88] = -3.32042076E-03;
    COFTD[89] = 6.66106286E-04;
    COFTD[90] = -2.61376693E-07;
    COFTD[91] = 3.23069179E-11;
    COFTD[92] = -3.32042076E-03;
    COFTD[93] = 6.66106286E-04;
    COFTD[94] = -2.61376693E-07;
    COFTD[95] = 3.23069179E-11;
    COFTD[96] = 8.84153069E-03;
    COFTD[97] = 6.58331966E-04;
    COFTD[98] = -2.59852601E-07;
    COFTD[99] = 3.22371458E-11;
    COFTD[100] = -2.45819473E-02;
    COFTD[101] = 6.90165396E-04;
    COFTD[102] = -2.67970561E-07;
    COFTD[103] = 3.29051275E-11;
    COFTD[104] = -2.33526384E-02;
    COFTD[105] = 6.89776192E-04;
    COFTD[106] = -2.67989830E-07;
    COFTD[107] = 3.29204573E-11;
    COFTD[108] = -2.33679974E-02;
    COFTD[109] = 6.90229858E-04;
    COFTD[110] = -2.68166087E-07;
    COFTD[111] = 3.29421091E-11;
    COFTD[112] = -2.14940840E-02;
    COFTD[113] = 6.89354633E-04;
    COFTD[114] = -2.68082952E-07;
    COFTD[115] = 3.29514659E-11;
    COFTD[116] = -6.72792134E-02;
    COFTD[117] = 7.17641289E-04;
    COFTD[118] = -2.72236296E-07;
    COFTD[119] = 3.29539785E-11;
    COFTD[120] = -6.71780088E-02;
    COFTD[121] = 7.18005938E-04;
    COFTD[122] = -2.72397329E-07;
    COFTD[123] = 3.29751562E-11;
    COFTD[124] = -4.79844727E-02;
    COFTD[125] = 7.10843780E-04;
    COFTD[126] = -2.72660435E-07;
    COFTD[127] = 3.32305029E-11;
    COFTD[128] = -4.21365511E-02;
    COFTD[129] = 7.07510653E-04;
    COFTD[130] = -2.72239988E-07;
    COFTD[131] = 3.32438058E-11;
    COFTD[132] = -7.02306173E-02;
    COFTD[133] = 7.23974404E-04;
    COFTD[134] = -2.74241180E-07;
    COFTD[135] = 3.31672060E-11;
    COFTD[136] = -8.70322417E-02;
    COFTD[137] = 7.22050015E-04;
    COFTD[138] = -2.70516884E-07;
    COFTD[139] = 3.24969401E-11;
    COFTD[140] = -8.71105785E-02;
    COFTD[141] = 7.22699924E-04;
    COFTD[142] = -2.70760373E-07;
    COFTD[143] = 3.25261903E-11;
    COFTD[144] = -7.07949543E-02;
    COFTD[145] = 7.18594456E-04;
    COFTD[146] = -2.72018661E-07;
    COFTD[147] = 3.28847379E-11;
    COFTD[148] = -6.40250531E-02;
    COFTD[149] = 7.15461906E-04;
    COFTD[150] = -2.71917633E-07;
    COFTD[151] = 3.29531717E-11;
    COFTD[152] = 2.86681699E-03;
    COFTD[153] = 6.62209231E-04;
    COFTD[154] = -2.60633543E-07;
    COFTD[155] = 3.22759160E-11;
    COFTD[156] = -3.31610907E-03;
    COFTD[157] = 6.65241324E-04;
    COFTD[158] = -2.61037286E-07;
    COFTD[159] = 3.22649662E-11;
    COFTD[160] = -8.08763122E-02;
    COFTD[161] = 7.29512044E-04;
    COFTD[162] = -2.74585491E-07;
    COFTD[163] = 3.30796405E-11;
    COFTD[164] = -7.63277334E-02;
    COFTD[165] = 7.27046669E-04;
    COFTD[166] = -2.74408897E-07;
    COFTD[167] = 3.31139101E-11;
    COFTD[168] = -7.08415895E-02;
    COFTD[169] = 7.19067819E-04;
    COFTD[170] = -2.72197849E-07;
    COFTD[171] = 3.29064003E-11;
    COFTD[172] = -8.67386596E-02;
    COFTD[173] = 7.35358769E-04;
    COFTD[174] = -2.75849568E-07;
    COFTD[175] = 3.31630521E-11;
    COFTD[176] = -9.17021209E-02;
    COFTD[177] = 7.37312160E-04;
    COFTD[178] = -2.75712291E-07;
    COFTD[179] = 3.30827046E-11;
    COFTD[180] = -7.94774360E-02;
    COFTD[181] = 7.34285207E-04;
    COFTD[182] = -2.76722783E-07;
    COFTD[183] = 3.33622747E-11;
    COFTD[184] = -7.94611889E-02;
    COFTD[185] = 7.34135101E-04;
    COFTD[186] = -2.76666214E-07;
    COFTD[187] = 3.33554547E-11;
    COFTD[188] = -1.18886225E-01;
    COFTD[189] = 7.52625081E-04;
    COFTD[190] = -2.76446668E-07;
    COFTD[191] = 3.28093656E-11;
    COFTD[192] = -1.31634388E-01;
    COFTD[193] = 7.57661003E-04;
    COFTD[194] = -2.75710681E-07;
    COFTD[195] = 3.25350740E-11;
    COFTD[196] = -1.01086570E-01;
    COFTD[197] = 7.45273061E-04;
    COFTD[198] = -2.77097798E-07;
    COFTD[199] = 3.31326430E-11;
    COFTD[200] = -1.25589299E-01;
    COFTD[201] = 7.55402817E-04;
    COFTD[202] = -2.76131324E-07;
    COFTD[203] = 3.26750851E-11;
    COFTD[204] = 2.72542957E-01;
    COFTD[205] = 3.58110418E-04;
    COFTD[206] = -1.58413280E-07;
    COFTD[207] = 2.13026843E-11;
    COFTD[208] = -1.58593068E-01;
    COFTD[209] = -3.77534860E-05;
    COFTD[210] = 1.76801591E-08;
    COFTD[211] = -2.63675347E-12;
    COFTD[212] = 0.00000000E+00;
    COFTD[213] = 0.00000000E+00;
    COFTD[214] = 0.00000000E+00;
    COFTD[215] = 0.00000000E+00;
    COFTD[216] = 4.20472602E-01;
    COFTD[217] = -3.28752268E-06;
    COFTD[218] = 6.10316432E-09;
    COFTD[219] = -8.09360145E-13;
    COFTD[220] = 4.53479915E-01;
    COFTD[221] = 4.08541226E-05;
    COFTD[222] = -1.74069280E-08;
    COFTD[223] = 2.81098139E-12;
    COFTD[224] = 4.26896405E-01;
    COFTD[225] = -3.33774806E-06;
    COFTD[226] = 6.19640588E-09;
    COFTD[227] = -8.21725207E-13;
    COFTD[228] = 9.22541582E-02;
    COFTD[229] = 4.74231405E-04;
    COFTD[230] = -1.95752173E-07;
    COFTD[231] = 2.49972004E-11;
    COFTD[232] = 4.53627122E-01;
    COFTD[233] = 2.51672948E-05;
    COFTD[234] = -9.21178306E-09;
    COFTD[235] = 1.58883013E-12;
    COFTD[236] = 3.51117608E-01;
    COFTD[237] = 2.33656884E-04;
    COFTD[238] = -1.07225357E-07;
    COFTD[239] = 1.49681441E-11;
    COFTD[240] = 3.64048983E-01;
    COFTD[241] = 8.51853930E-05;
    COFTD[242] = -3.98782972E-08;
    COFTD[243] = 5.95442132E-12;
    COFTD[244] = 3.72390318E-01;
    COFTD[245] = 8.32178961E-05;
    COFTD[246] = -3.89131131E-08;
    COFTD[247] = 5.82905441E-12;
    COFTD[248] = 4.55234969E-01;
    COFTD[249] = 4.10122358E-05;
    COFTD[250] = -1.74742962E-08;
    COFTD[251] = 2.82186043E-12;
    COFTD[252] = 4.56891826E-01;
    COFTD[253] = 4.11615025E-05;
    COFTD[254] = -1.75378950E-08;
    COFTD[255] = 2.83213077E-12;
    COFTD[256] = 1.94872921E-01;
    COFTD[257] = 4.13931000E-04;
    COFTD[258] = -1.78165405E-07;
    COFTD[259] = 2.34364688E-11;
    COFTD[260] = 3.33963995E-01;
    COFTD[261] = 2.27614782E-04;
    COFTD[262] = -1.04334948E-07;
    COFTD[263] = 1.45440845E-11;
    COFTD[264] = 3.37244670E-01;
    COFTD[265] = 2.15795621E-04;
    COFTD[266] = -9.92082642E-08;
    COFTD[267] = 1.38813543E-11;
    COFTD[268] = 3.32407303E-01;
    COFTD[269] = 2.26553810E-04;
    COFTD[270] = -1.03848615E-07;
    COFTD[271] = 1.44762907E-11;
    COFTD[272] = 3.12900731E-01;
    COFTD[273] = 2.39234977E-04;
    COFTD[274] = -1.09042966E-07;
    COFTD[275] = 1.50987426E-11;
    COFTD[276] = 3.14870276E-01;
    COFTD[277] = 2.40740835E-04;
    COFTD[278] = -1.09729333E-07;
    COFTD[279] = 1.51937812E-11;
    COFTD[280] = 2.03701481E-01;
    COFTD[281] = 4.08732725E-04;
    COFTD[282] = -1.76483350E-07;
    COFTD[283] = 2.32710218E-11;
    COFTD[284] = 2.31981193E-01;
    COFTD[285] = 3.93734933E-04;
    COFTD[286] = -1.71622180E-07;
    COFTD[287] = 2.27980120E-11;
    COFTD[288] = 4.26650826E-01;
    COFTD[289] = 1.10258123E-04;
    COFTD[290] = -5.17017276E-08;
    COFTD[291] = 7.66706104E-12;
    COFTD[292] = 2.33003847E-01;
    COFTD[293] = 3.95470654E-04;
    COFTD[294] = -1.72378750E-07;
    COFTD[295] = 2.28985136E-11;
    COFTD[296] = 2.92930561E-01;
    COFTD[297] = 3.10048132E-04;
    COFTD[298] = -1.38860343E-07;
    COFTD[299] = 1.88831074E-11;
    COFTD[300] = 2.92930561E-01;
    COFTD[301] = 3.10048132E-04;
    COFTD[302] = -1.38860343E-07;
    COFTD[303] = 1.88831074E-11;
    COFTD[304] = 3.05188022E-01;
    COFTD[305] = 2.96541355E-04;
    COFTD[306] = -1.33447888E-07;
    COFTD[307] = 1.82309847E-11;
    COFTD[308] = 2.80771189E-01;
    COFTD[309] = 3.45296283E-04;
    COFTD[310] = -1.53328644E-07;
    COFTD[311] = 2.06879584E-11;
    COFTD[312] = 2.82355558E-01;
    COFTD[313] = 3.44158251E-04;
    COFTD[314] = -1.52901778E-07;
    COFTD[315] = 2.06397774E-11;
    COFTD[316] = 2.82727459E-01;
    COFTD[317] = 3.44611553E-04;
    COFTD[318] = -1.53103170E-07;
    COFTD[319] = 2.06669627E-11;
    COFTD[320] = 2.84901692E-01;
    COFTD[321] = 3.42611250E-04;
    COFTD[322] = -1.52332416E-07;
    COFTD[323] = 2.05771404E-11;
    COFTD[324] = 2.38185491E-01;
    COFTD[325] = 4.03502216E-04;
    COFTD[326] = -1.75898559E-07;
    COFTD[327] = 2.33680960E-11;
    COFTD[328] = 2.09860393E-01;
    COFTD[329] = 4.39450116E-04;
    COFTD[330] = -1.89298200E-07;
    COFTD[331] = 2.49157649E-11;
    COFTD[332] = 2.62101185E-01;
    COFTD[333] = 3.81645143E-04;
    COFTD[334] = -1.67825155E-07;
    COFTD[335] = 2.24547350E-11;
    COFTD[336] = 2.67937228E-01;
    COFTD[337] = 3.74639675E-04;
    COFTD[338] = -1.65130982E-07;
    COFTD[339] = 2.21378638E-11;
    COFTD[340] = 2.38815693E-01;
    COFTD[341] = 4.12622038E-04;
    COFTD[342] = -1.79671905E-07;
    COFTD[343] = 2.38479790E-11;
    COFTD[344] = 2.12061419E-01;
    COFTD[345] = 4.25507174E-04;
    COFTD[346] = -1.83726252E-07;
    COFTD[347] = 2.42260679E-11;
    COFTD[348] = 2.12443856E-01;
    COFTD[349] = 4.26274545E-04;
    COFTD[350] = -1.84057589E-07;
    COFTD[351] = 2.42697578E-11;
    COFTD[352] = 2.33636760E-01;
    COFTD[353] = 4.07414700E-04;
    COFTD[354] = -1.77311195E-07;
    COFTD[355] = 2.35248147E-11;
    COFTD[356] = 2.13688259E-01;
    COFTD[357] = 4.33054284E-04;
    COFTD[358] = -1.86881616E-07;
    COFTD[359] = 2.46317025E-11;
    COFTD[360] = 2.99205069E-01;
    COFTD[361] = 3.03160248E-04;
    COFTD[362] = -1.36108905E-07;
    COFTD[363] = 1.85523351E-11;
    COFTD[364] = 2.92168810E-01;
    COFTD[365] = 3.09241868E-04;
    COFTD[366] = -1.38499244E-07;
    COFTD[367] = 1.88340029E-11;
    COFTD[368] = 2.05727071E-01;
    COFTD[369] = 4.53872190E-04;
    COFTD[370] = -1.94951739E-07;
    COFTD[371] = 2.56044686E-11;
    COFTD[372] = 2.04596332E-01;
    COFTD[373] = 4.54674220E-04;
    COFTD[374] = -1.95218824E-07;
    COFTD[375] = 2.56319180E-11;
    COFTD[376] = 2.33944976E-01;
    COFTD[377] = 4.07952165E-04;
    COFTD[378] = -1.77545105E-07;
    COFTD[379] = 2.35558488E-11;
    COFTD[380] = 2.09647812E-01;
    COFTD[381] = 4.55574424E-04;
    COFTD[382] = -1.95844758E-07;
    COFTD[383] = 2.57377604E-11;
    COFTD[384] = 2.16849866E-01;
    COFTD[385] = 4.46455653E-04;
    COFTD[386] = -1.92496189E-07;
    COFTD[387] = 2.53547586E-11;
    COFTD[388] = 2.33012530E-01;
    COFTD[389] = 4.31676214E-04;
    COFTD[390] = -1.87226262E-07;
    COFTD[391] = 2.47732585E-11;
    COFTD[392] = 2.32917243E-01;
    COFTD[393] = 4.31499688E-04;
    COFTD[394] = -1.87149699E-07;
    COFTD[395] = 2.47631279E-11;
    COFTD[396] = 1.83133284E-01;
    COFTD[397] = 4.98249243E-04;
    COFTD[398] = -2.11648739E-07;
    COFTD[399] = 2.75709016E-11;
    COFTD[400] = 1.68767197E-01;
    COFTD[401] = 5.18907171E-04;
    COFTD[402] = -2.19073764E-07;
    COFTD[403] = 2.84128692E-11;
    COFTD[404] = 2.08814161E-01;
    COFTD[405] = 4.64765526E-04;
    COFTD[406] = -1.99534834E-07;
    COFTD[407] = 2.61969502E-11;
    COFTD[408] = 1.73974000E-01;
    COFTD[409] = 5.10975705E-04;
    COFTD[410] = -2.16215072E-07;
    COFTD[411] = 2.80874182E-11;
    COFTD[412] = 4.54134354E-01;
    COFTD[413] = 2.42490079E-05;
    COFTD[414] = -8.71977610E-09;
    COFTD[415] = 1.51469519E-12;
}


/*compute the critical parameters for each species */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void GET_CRITPARAMS(amrex::Real *  Tci, amrex::Real *  ai, amrex::Real *  bi, amrex::Real *  acentric_i)
{

    amrex::Real   EPS[52];
    amrex::Real   SIG[52];
    amrex::Real    wt[52];
    amrex::Real avogadro = 6.02214199e23;
    amrex::Real boltzmann = 1.3806503e-16; //we work in CGS
    amrex::Real Rcst = 83.144598; //in bar [CGS] !

    egtransetEPS(EPS);
    egtransetSIG(SIG);
    get_mw(wt);

    /*species 0: H */
    Tci[0] = 1.316 * EPS[0] ; 
    ai[0] = (5.55 * pow(avogadro,2.0) * EPS[0]*boltzmann * pow(1e-8*SIG[0],3.0) ) / (pow(wt[0],2.0)); 
    bi[0] = 0.855 * avogadro * pow(1e-8*SIG[0],3.0) / (wt[0]); 
    acentric_i[0] = 0.0 ;

    /*species 1: H2 */
    /*Imported from NIST */
    Tci[1] = 33.145000 ; 
    ai[1] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[1],2.0) / (pow(2.015880,2.0) * 12.964000); 
    bi[1] = 0.08664 * Rcst * Tci[1] / (2.015880 * 12.964000); 
    acentric_i[1] = -0.219000 ;

    /*species 2: O */
    Tci[2] = 1.316 * EPS[2] ; 
    ai[2] = (5.55 * pow(avogadro,2.0) * EPS[2]*boltzmann * pow(1e-8*SIG[2],3.0) ) / (pow(wt[2],2.0)); 
    bi[2] = 0.855 * avogadro * pow(1e-8*SIG[2],3.0) / (wt[2]); 
    acentric_i[2] = 0.0 ;

    /*species 3: O2 */
    /*Imported from NIST */
    Tci[3] = 154.581000 ; 
    ai[3] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[3],2.0) / (pow(31.998800,2.0) * 50.430466); 
    bi[3] = 0.08664 * Rcst * Tci[3] / (31.998800 * 50.430466); 
    acentric_i[3] = 0.022200 ;

    /*species 4: OH */
    Tci[4] = 1.316 * EPS[4] ; 
    ai[4] = (5.55 * pow(avogadro,2.0) * EPS[4]*boltzmann * pow(1e-8*SIG[4],3.0) ) / (pow(wt[4],2.0)); 
    bi[4] = 0.855 * avogadro * pow(1e-8*SIG[4],3.0) / (wt[4]); 
    acentric_i[4] = 0.0 ;

    /*species 5: H2O */
    /*Imported from NIST */
    Tci[5] = 647.096000 ; 
    ai[5] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[5],2.0) / (pow(18.015340,2.0) * 220.640000); 
    bi[5] = 0.08664 * Rcst * Tci[5] / (18.015340 * 220.640000); 
    acentric_i[5] = 0.344300 ;

    /*species 6: CO */
    /*Imported from NIST */
    Tci[6] = 132.850000 ; 
    ai[6] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[6],2.0) / (pow(28.010000,2.0) * 34.940000); 
    bi[6] = 0.08664 * Rcst * Tci[6] / (28.010000 * 34.940000); 
    acentric_i[6] = 0.045000 ;

    /*species 7: CO2 */
    /*Imported from NIST */
    Tci[7] = 304.120000 ; 
    ai[7] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[7],2.0) / (pow(44.009950,2.0) * 73.740000); 
    bi[7] = 0.08664 * Rcst * Tci[7] / (44.009950 * 73.740000); 
    acentric_i[7] = 0.225000 ;

    /*species 8: CH3 */
    Tci[8] = 1.316 * EPS[8] ; 
    ai[8] = (5.55 * pow(avogadro,2.0) * EPS[8]*boltzmann * pow(1e-8*SIG[8],3.0) ) / (pow(wt[8],2.0)); 
    bi[8] = 0.855 * avogadro * pow(1e-8*SIG[8],3.0) / (wt[8]); 
    acentric_i[8] = 0.0 ;

    /*species 9: CH4 */
    /*Imported from NIST */
    Tci[9] = 190.560000 ; 
    ai[9] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[9],2.0) / (pow(16.043030,2.0) * 45.990000); 
    bi[9] = 0.08664 * Rcst * Tci[9] / (16.043030 * 45.990000); 
    acentric_i[9] = 0.011000 ;

    /*species 10: HO2 */
    Tci[10] = 1.316 * EPS[10] ; 
    ai[10] = (5.55 * pow(avogadro,2.0) * EPS[10]*boltzmann * pow(1e-8*SIG[10],3.0) ) / (pow(wt[10],2.0)); 
    bi[10] = 0.855 * avogadro * pow(1e-8*SIG[10],3.0) / (wt[10]); 
    acentric_i[10] = 0.0 ;

    /*species 11: H2O2 */
    Tci[11] = 1.316 * EPS[11] ; 
    ai[11] = (5.55 * pow(avogadro,2.0) * EPS[11]*boltzmann * pow(1e-8*SIG[11],3.0) ) / (pow(wt[11],2.0)); 
    bi[11] = 0.855 * avogadro * pow(1e-8*SIG[11],3.0) / (wt[11]); 
    acentric_i[11] = 0.0 ;

    /*species 12: CH2O */
    Tci[12] = 1.316 * EPS[12] ; 
    ai[12] = (5.55 * pow(avogadro,2.0) * EPS[12]*boltzmann * pow(1e-8*SIG[12],3.0) ) / (pow(wt[12],2.0)); 
    bi[12] = 0.855 * avogadro * pow(1e-8*SIG[12],3.0) / (wt[12]); 
    acentric_i[12] = 0.0 ;

    /*species 13: C2H6 */
    /*Imported from NIST */
    Tci[13] = 305.320000 ; 
    ai[13] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[13],2.0) / (pow(30.070120,2.0) * 48.720000); 
    bi[13] = 0.08664 * Rcst * Tci[13] / (30.070120 * 48.720000); 
    acentric_i[13] = 0.099000 ;

    /*species 14: C2H4 */
    /*Imported from NIST */
    Tci[14] = 282.340000 ; 
    ai[14] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[14],2.0) / (pow(28.054000,2.0) * 50.410000); 
    bi[14] = 0.08664 * Rcst * Tci[14] / (28.054000 * 50.410000); 
    acentric_i[14] = 0.087000 ;

    /*species 15: C2H5 */
    Tci[15] = 1.316 * EPS[15] ; 
    ai[15] = (5.55 * pow(avogadro,2.0) * EPS[15]*boltzmann * pow(1e-8*SIG[15],3.0) ) / (pow(wt[15],2.0)); 
    bi[15] = 0.855 * avogadro * pow(1e-8*SIG[15],3.0) / (wt[15]); 
    acentric_i[15] = 0.0 ;

    /*species 16: C2H */
    Tci[16] = 1.316 * EPS[16] ; 
    ai[16] = (5.55 * pow(avogadro,2.0) * EPS[16]*boltzmann * pow(1e-8*SIG[16],3.0) ) / (pow(wt[16],2.0)); 
    bi[16] = 0.855 * avogadro * pow(1e-8*SIG[16],3.0) / (wt[16]); 
    acentric_i[16] = 0.0 ;

    /*species 17: C2H2 */
    /*Imported from NIST */
    Tci[17] = 308.300000 ; 
    ai[17] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[17],2.0) / (pow(26.038000,2.0) * 61.140000); 
    bi[17] = 0.08664 * Rcst * Tci[17] / (26.038000 * 61.140000); 
    acentric_i[17] = 0.189000 ;

    /*species 18: CH3OH */
    /*Imported from NIST */
    Tci[18] = 512.640000 ; 
    ai[18] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[18],2.0) / (pow(32.042000,2.0) * 80.970000); 
    bi[18] = 0.08664 * Rcst * Tci[18] / (32.042000 * 80.970000); 
    acentric_i[18] = 0.565000 ;

    /*species 19: CH2CO */
    Tci[19] = 1.316 * EPS[19] ; 
    ai[19] = (5.55 * pow(avogadro,2.0) * EPS[19]*boltzmann * pow(1e-8*SIG[19],3.0) ) / (pow(wt[19],2.0)); 
    bi[19] = 0.855 * avogadro * pow(1e-8*SIG[19],3.0) / (wt[19]); 
    acentric_i[19] = 0.0 ;

    /*species 20: HCCO */
    Tci[20] = 1.316 * EPS[20] ; 
    ai[20] = (5.55 * pow(avogadro,2.0) * EPS[20]*boltzmann * pow(1e-8*SIG[20],3.0) ) / (pow(wt[20],2.0)); 
    bi[20] = 0.855 * avogadro * pow(1e-8*SIG[20],3.0) / (wt[20]); 
    acentric_i[20] = 0.0 ;

    /*species 21: CH3CHO */
    Tci[21] = 1.316 * EPS[21] ; 
    ai[21] = (5.55 * pow(avogadro,2.0) * EPS[21]*boltzmann * pow(1e-8*SIG[21],3.0) ) / (pow(wt[21],2.0)); 
    bi[21] = 0.855 * avogadro * pow(1e-8*SIG[21],3.0) / (wt[21]); 
    acentric_i[21] = 0.0 ;

    /*species 22: C3H4-A */
    Tci[22] = 1.316 * EPS[22] ; 
    ai[22] = (5.55 * pow(avogadro,2.0) * EPS[22]*boltzmann * pow(1e-8*SIG[22],3.0) ) / (pow(wt[22],2.0)); 
    bi[22] = 0.855 * avogadro * pow(1e-8*SIG[22],3.0) / (wt[22]); 
    acentric_i[22] = 0.0 ;

    /*species 23: C3H4-P */
    Tci[23] = 1.316 * EPS[23] ; 
    ai[23] = (5.55 * pow(avogadro,2.0) * EPS[23]*boltzmann * pow(1e-8*SIG[23],3.0) ) / (pow(wt[23],2.0)); 
    bi[23] = 0.855 * avogadro * pow(1e-8*SIG[23],3.0) / (wt[23]); 
    acentric_i[23] = 0.0 ;

    /*species 24: C3H6 */
    Tci[24] = 1.316 * EPS[24] ; 
    ai[24] = (5.55 * pow(avogadro,2.0) * EPS[24]*boltzmann * pow(1e-8*SIG[24],3.0) ) / (pow(wt[24],2.0)); 
    bi[24] = 0.855 * avogadro * pow(1e-8*SIG[24],3.0) / (wt[24]); 
    acentric_i[24] = 0.0 ;

    /*species 25: C4H6 */
    Tci[25] = 1.316 * EPS[25] ; 
    ai[25] = (5.55 * pow(avogadro,2.0) * EPS[25]*boltzmann * pow(1e-8*SIG[25],3.0) ) / (pow(wt[25],2.0)); 
    bi[25] = 0.855 * avogadro * pow(1e-8*SIG[25],3.0) / (wt[25]); 
    acentric_i[25] = 0.0 ;

    /*species 26: C4H7 */
    Tci[26] = 1.316 * EPS[26] ; 
    ai[26] = (5.55 * pow(avogadro,2.0) * EPS[26]*boltzmann * pow(1e-8*SIG[26],3.0) ) / (pow(wt[26],2.0)); 
    bi[26] = 0.855 * avogadro * pow(1e-8*SIG[26],3.0) / (wt[26]); 
    acentric_i[26] = 0.0 ;

    /*species 27: C4H8-1 */
    Tci[27] = 1.316 * EPS[27] ; 
    ai[27] = (5.55 * pow(avogadro,2.0) * EPS[27]*boltzmann * pow(1e-8*SIG[27],3.0) ) / (pow(wt[27],2.0)); 
    bi[27] = 0.855 * avogadro * pow(1e-8*SIG[27],3.0) / (wt[27]); 
    acentric_i[27] = 0.0 ;

    /*species 28: PC4H9 */
    Tci[28] = 1.316 * EPS[28] ; 
    ai[28] = (5.55 * pow(avogadro,2.0) * EPS[28]*boltzmann * pow(1e-8*SIG[28],3.0) ) / (pow(wt[28],2.0)); 
    bi[28] = 0.855 * avogadro * pow(1e-8*SIG[28],3.0) / (wt[28]); 
    acentric_i[28] = 0.0 ;

    /*species 29: CH3COCH2 */
    Tci[29] = 1.316 * EPS[29] ; 
    ai[29] = (5.55 * pow(avogadro,2.0) * EPS[29]*boltzmann * pow(1e-8*SIG[29],3.0) ) / (pow(wt[29],2.0)); 
    bi[29] = 0.855 * avogadro * pow(1e-8*SIG[29],3.0) / (wt[29]); 
    acentric_i[29] = 0.0 ;

    /*species 30: C2H5CHO */
    Tci[30] = 1.316 * EPS[30] ; 
    ai[30] = (5.55 * pow(avogadro,2.0) * EPS[30]*boltzmann * pow(1e-8*SIG[30],3.0) ) / (pow(wt[30],2.0)); 
    bi[30] = 0.855 * avogadro * pow(1e-8*SIG[30],3.0) / (wt[30]); 
    acentric_i[30] = 0.0 ;

    /*species 31: C5H9 */
    Tci[31] = 1.316 * EPS[31] ; 
    ai[31] = (5.55 * pow(avogadro,2.0) * EPS[31]*boltzmann * pow(1e-8*SIG[31],3.0) ) / (pow(wt[31],2.0)); 
    bi[31] = 0.855 * avogadro * pow(1e-8*SIG[31],3.0) / (wt[31]); 
    acentric_i[31] = 0.0 ;

    /*species 32: C5H10-1 */
    Tci[32] = 1.316 * EPS[32] ; 
    ai[32] = (5.55 * pow(avogadro,2.0) * EPS[32]*boltzmann * pow(1e-8*SIG[32],3.0) ) / (pow(wt[32],2.0)); 
    bi[32] = 0.855 * avogadro * pow(1e-8*SIG[32],3.0) / (wt[32]); 
    acentric_i[32] = 0.0 ;

    /*species 33: C5H11-1 */
    Tci[33] = 1.316 * EPS[33] ; 
    ai[33] = (5.55 * pow(avogadro,2.0) * EPS[33]*boltzmann * pow(1e-8*SIG[33],3.0) ) / (pow(wt[33],2.0)); 
    bi[33] = 0.855 * avogadro * pow(1e-8*SIG[33],3.0) / (wt[33]); 
    acentric_i[33] = 0.0 ;

    /*species 34: CH3O2 */
    Tci[34] = 1.316 * EPS[34] ; 
    ai[34] = (5.55 * pow(avogadro,2.0) * EPS[34]*boltzmann * pow(1e-8*SIG[34],3.0) ) / (pow(wt[34],2.0)); 
    bi[34] = 0.855 * avogadro * pow(1e-8*SIG[34],3.0) / (wt[34]); 
    acentric_i[34] = 0.0 ;

    /*species 35: CH3O2H */
    Tci[35] = 1.316 * EPS[35] ; 
    ai[35] = (5.55 * pow(avogadro,2.0) * EPS[35]*boltzmann * pow(1e-8*SIG[35],3.0) ) / (pow(wt[35],2.0)); 
    bi[35] = 0.855 * avogadro * pow(1e-8*SIG[35],3.0) / (wt[35]); 
    acentric_i[35] = 0.0 ;

    /*species 36: C2H3CO */
    Tci[36] = 1.316 * EPS[36] ; 
    ai[36] = (5.55 * pow(avogadro,2.0) * EPS[36]*boltzmann * pow(1e-8*SIG[36],3.0) ) / (pow(wt[36],2.0)); 
    bi[36] = 0.855 * avogadro * pow(1e-8*SIG[36],3.0) / (wt[36]); 
    acentric_i[36] = 0.0 ;

    /*species 37: C2H3CHO */
    Tci[37] = 1.316 * EPS[37] ; 
    ai[37] = (5.55 * pow(avogadro,2.0) * EPS[37]*boltzmann * pow(1e-8*SIG[37],3.0) ) / (pow(wt[37],2.0)); 
    bi[37] = 0.855 * avogadro * pow(1e-8*SIG[37],3.0) / (wt[37]); 
    acentric_i[37] = 0.0 ;

    /*species 38: C3H5-A */
    Tci[38] = 1.316 * EPS[38] ; 
    ai[38] = (5.55 * pow(avogadro,2.0) * EPS[38]*boltzmann * pow(1e-8*SIG[38],3.0) ) / (pow(wt[38],2.0)); 
    bi[38] = 0.855 * avogadro * pow(1e-8*SIG[38],3.0) / (wt[38]); 
    acentric_i[38] = 0.0 ;

    /*species 39: C3H3 */
    Tci[39] = 1.316 * EPS[39] ; 
    ai[39] = (5.55 * pow(avogadro,2.0) * EPS[39]*boltzmann * pow(1e-8*SIG[39],3.0) ) / (pow(wt[39],2.0)); 
    bi[39] = 0.855 * avogadro * pow(1e-8*SIG[39],3.0) / (wt[39]); 
    acentric_i[39] = 0.0 ;

    /*species 40: NC3H7CHO */
    Tci[40] = 1.316 * EPS[40] ; 
    ai[40] = (5.55 * pow(avogadro,2.0) * EPS[40]*boltzmann * pow(1e-8*SIG[40],3.0) ) / (pow(wt[40],2.0)); 
    bi[40] = 0.855 * avogadro * pow(1e-8*SIG[40],3.0) / (wt[40]); 
    acentric_i[40] = 0.0 ;

    /*species 41: C2H5COCH2 */
    Tci[41] = 1.316 * EPS[41] ; 
    ai[41] = (5.55 * pow(avogadro,2.0) * EPS[41]*boltzmann * pow(1e-8*SIG[41],3.0) ) / (pow(wt[41],2.0)); 
    bi[41] = 0.855 * avogadro * pow(1e-8*SIG[41],3.0) / (wt[41]); 
    acentric_i[41] = 0.0 ;

    /*species 42: CH3CHCO */
    Tci[42] = 1.316 * EPS[42] ; 
    ai[42] = (5.55 * pow(avogadro,2.0) * EPS[42]*boltzmann * pow(1e-8*SIG[42],3.0) ) / (pow(wt[42],2.0)); 
    bi[42] = 0.855 * avogadro * pow(1e-8*SIG[42],3.0) / (wt[42]); 
    acentric_i[42] = 0.0 ;

    /*species 43: NC3H7COCH3 */
    Tci[43] = 1.316 * EPS[43] ; 
    ai[43] = (5.55 * pow(avogadro,2.0) * EPS[43]*boltzmann * pow(1e-8*SIG[43],3.0) ) / (pow(wt[43],2.0)); 
    bi[43] = 0.855 * avogadro * pow(1e-8*SIG[43],3.0) / (wt[43]); 
    acentric_i[43] = 0.0 ;

    /*species 44: NC3H7COCH2 */
    Tci[44] = 1.316 * EPS[44] ; 
    ai[44] = (5.55 * pow(avogadro,2.0) * EPS[44]*boltzmann * pow(1e-8*SIG[44],3.0) ) / (pow(wt[44],2.0)); 
    bi[44] = 0.855 * avogadro * pow(1e-8*SIG[44],3.0) / (wt[44]); 
    acentric_i[44] = 0.0 ;

    /*species 45: NC7H16 */
    Tci[45] = 1.316 * EPS[45] ; 
    ai[45] = (5.55 * pow(avogadro,2.0) * EPS[45]*boltzmann * pow(1e-8*SIG[45],3.0) ) / (pow(wt[45],2.0)); 
    bi[45] = 0.855 * avogadro * pow(1e-8*SIG[45],3.0) / (wt[45]); 
    acentric_i[45] = 0.0 ;

    /*species 46: C7H15-1 */
    Tci[46] = 1.316 * EPS[46] ; 
    ai[46] = (5.55 * pow(avogadro,2.0) * EPS[46]*boltzmann * pow(1e-8*SIG[46],3.0) ) / (pow(wt[46],2.0)); 
    bi[46] = 0.855 * avogadro * pow(1e-8*SIG[46],3.0) / (wt[46]); 
    acentric_i[46] = 0.0 ;

    /*species 47: C7H15O2 */
    Tci[47] = 1.316 * EPS[47] ; 
    ai[47] = (5.55 * pow(avogadro,2.0) * EPS[47]*boltzmann * pow(1e-8*SIG[47],3.0) ) / (pow(wt[47],2.0)); 
    bi[47] = 0.855 * avogadro * pow(1e-8*SIG[47],3.0) / (wt[47]); 
    acentric_i[47] = 0.0 ;

    /*species 48: C7H14OOHO2 */
    Tci[48] = 1.316 * EPS[48] ; 
    ai[48] = (5.55 * pow(avogadro,2.0) * EPS[48]*boltzmann * pow(1e-8*SIG[48],3.0) ) / (pow(wt[48],2.0)); 
    bi[48] = 0.855 * avogadro * pow(1e-8*SIG[48],3.0) / (wt[48]); 
    acentric_i[48] = 0.0 ;

    /*species 49: C7H14O */
    Tci[49] = 1.316 * EPS[49] ; 
    ai[49] = (5.55 * pow(avogadro,2.0) * EPS[49]*boltzmann * pow(1e-8*SIG[49],3.0) ) / (pow(wt[49],2.0)); 
    bi[49] = 0.855 * avogadro * pow(1e-8*SIG[49],3.0) / (wt[49]); 
    acentric_i[49] = 0.0 ;

    /*species 50: NC7KET */
    Tci[50] = 1.316 * EPS[50] ; 
    ai[50] = (5.55 * pow(avogadro,2.0) * EPS[50]*boltzmann * pow(1e-8*SIG[50],3.0) ) / (pow(wt[50],2.0)); 
    bi[50] = 0.855 * avogadro * pow(1e-8*SIG[50],3.0) / (wt[50]); 
    acentric_i[50] = 0.0 ;

    /*species 51: N2 */
    /*Imported from NIST */
    Tci[51] = 126.192000 ; 
    ai[51] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[51],2.0) / (pow(28.013400,2.0) * 33.958000); 
    bi[51] = 0.08664 * Rcst * Tci[51] / (28.013400 * 33.958000); 
    acentric_i[51] = 0.037200 ;

    return;
}

/* gauss-jordan solver external routine */
/* Replace this routine with the one generated by the Gauss Jordan solver of DW */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void sgjsolve(amrex::Real* /*A*/, amrex::Real* /*x*/, amrex::Real* /*b*/) {
    amrex::Abort("sgjsolve not implemented, choose a different solver ");
}

/* Replace this routine with the one generated by the Gauss Jordan solver of DW */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void sgjsolve_simplified(amrex::Real* /*A*/, amrex::Real* /*x*/, amrex::Real* /*b*/) {
    amrex::Abort("sgjsolve_simplified not implemented, choose a different solver ");
}

#endif

/* End of file  */
