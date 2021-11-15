#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "mechanism.H"



/*save atomic weights into array */
void atomicWeight(amrex::Real *  awt)
{
    awt[0] = 1.007970; /*H */
    awt[1] = 12.011150; /*C */
    awt[2] = 15.999400; /*O */
    awt[3] = 14.006700; /*N */
}



/*get atomic weight for all elements */
void CKAWT( amrex::Real *  awt)
{
    atomicWeight(awt);
}



/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * ncf)
{
    int id; /*loop counter */
    int kd = 4; 
    /*Zero ncf */
    for (id = 0; id < kd * 52; ++ id) {
         ncf[id] = 0; 
    }

    /*H */
    ncf[ 0 * kd + 0 ] = 1; /*H */

    /*H2 */
    ncf[ 1 * kd + 0 ] = 2; /*H */

    /*O */
    ncf[ 2 * kd + 2 ] = 1; /*O */

    /*O2 */
    ncf[ 3 * kd + 2 ] = 2; /*O */

    /*OH */
    ncf[ 4 * kd + 0 ] = 1; /*H */
    ncf[ 4 * kd + 2 ] = 1; /*O */

    /*H2O */
    ncf[ 5 * kd + 0 ] = 2; /*H */
    ncf[ 5 * kd + 2 ] = 1; /*O */

    /*CO */
    ncf[ 6 * kd + 1 ] = 1; /*C */
    ncf[ 6 * kd + 2 ] = 1; /*O */

    /*CO2 */
    ncf[ 7 * kd + 1 ] = 1; /*C */
    ncf[ 7 * kd + 2 ] = 2; /*O */

    /*CH3 */
    ncf[ 8 * kd + 1 ] = 1; /*C */
    ncf[ 8 * kd + 0 ] = 3; /*H */

    /*CH4 */
    ncf[ 9 * kd + 1 ] = 1; /*C */
    ncf[ 9 * kd + 0 ] = 4; /*H */

    /*HO2 */
    ncf[ 10 * kd + 0 ] = 1; /*H */
    ncf[ 10 * kd + 2 ] = 2; /*O */

    /*H2O2 */
    ncf[ 11 * kd + 0 ] = 2; /*H */
    ncf[ 11 * kd + 2 ] = 2; /*O */

    /*CH2O */
    ncf[ 12 * kd + 1 ] = 1; /*C */
    ncf[ 12 * kd + 0 ] = 2; /*H */
    ncf[ 12 * kd + 2 ] = 1; /*O */

    /*C2H6 */
    ncf[ 13 * kd + 1 ] = 2; /*C */
    ncf[ 13 * kd + 0 ] = 6; /*H */

    /*C2H4 */
    ncf[ 14 * kd + 1 ] = 2; /*C */
    ncf[ 14 * kd + 0 ] = 4; /*H */

    /*C2H5 */
    ncf[ 15 * kd + 1 ] = 2; /*C */
    ncf[ 15 * kd + 0 ] = 5; /*H */

    /*C2H */
    ncf[ 16 * kd + 1 ] = 2; /*C */
    ncf[ 16 * kd + 0 ] = 1; /*H */

    /*C2H2 */
    ncf[ 17 * kd + 1 ] = 2; /*C */
    ncf[ 17 * kd + 0 ] = 2; /*H */

    /*CH3OH */
    ncf[ 18 * kd + 1 ] = 1; /*C */
    ncf[ 18 * kd + 0 ] = 4; /*H */
    ncf[ 18 * kd + 2 ] = 1; /*O */

    /*CH2CO */
    ncf[ 19 * kd + 1 ] = 2; /*C */
    ncf[ 19 * kd + 0 ] = 2; /*H */
    ncf[ 19 * kd + 2 ] = 1; /*O */

    /*HCCO */
    ncf[ 20 * kd + 0 ] = 1; /*H */
    ncf[ 20 * kd + 1 ] = 2; /*C */
    ncf[ 20 * kd + 2 ] = 1; /*O */

    /*CH3CHO */
    ncf[ 21 * kd + 1 ] = 2; /*C */
    ncf[ 21 * kd + 2 ] = 1; /*O */
    ncf[ 21 * kd + 0 ] = 4; /*H */

    /*C3H4-A */
    ncf[ 22 * kd + 0 ] = 4; /*H */
    ncf[ 22 * kd + 1 ] = 3; /*C */

    /*C3H4-P */
    ncf[ 23 * kd + 0 ] = 4; /*H */
    ncf[ 23 * kd + 1 ] = 3; /*C */

    /*C3H6 */
    ncf[ 24 * kd + 1 ] = 3; /*C */
    ncf[ 24 * kd + 0 ] = 6; /*H */

    /*C4H6 */
    ncf[ 25 * kd + 1 ] = 4; /*C */
    ncf[ 25 * kd + 0 ] = 6; /*H */

    /*C4H7 */
    ncf[ 26 * kd + 1 ] = 4; /*C */
    ncf[ 26 * kd + 0 ] = 7; /*H */

    /*C4H8-1 */
    ncf[ 27 * kd + 1 ] = 4; /*C */
    ncf[ 27 * kd + 0 ] = 8; /*H */

    /*PC4H9 */
    ncf[ 28 * kd + 1 ] = 4; /*C */
    ncf[ 28 * kd + 0 ] = 9; /*H */

    /*CH3COCH2 */
    ncf[ 29 * kd + 1 ] = 3; /*C */
    ncf[ 29 * kd + 0 ] = 5; /*H */
    ncf[ 29 * kd + 2 ] = 1; /*O */

    /*C2H5CHO */
    ncf[ 30 * kd + 1 ] = 3; /*C */
    ncf[ 30 * kd + 0 ] = 6; /*H */
    ncf[ 30 * kd + 2 ] = 1; /*O */

    /*C5H9 */
    ncf[ 31 * kd + 1 ] = 5; /*C */
    ncf[ 31 * kd + 0 ] = 9; /*H */

    /*C5H10-1 */
    ncf[ 32 * kd + 1 ] = 5; /*C */
    ncf[ 32 * kd + 0 ] = 10; /*H */

    /*C5H11-1 */
    ncf[ 33 * kd + 1 ] = 5; /*C */
    ncf[ 33 * kd + 0 ] = 11; /*H */

    /*CH3O2 */
    ncf[ 34 * kd + 1 ] = 1; /*C */
    ncf[ 34 * kd + 0 ] = 3; /*H */
    ncf[ 34 * kd + 2 ] = 2; /*O */

    /*CH3O2H */
    ncf[ 35 * kd + 1 ] = 1; /*C */
    ncf[ 35 * kd + 0 ] = 4; /*H */
    ncf[ 35 * kd + 2 ] = 2; /*O */

    /*C2H3CO */
    ncf[ 36 * kd + 1 ] = 3; /*C */
    ncf[ 36 * kd + 0 ] = 3; /*H */
    ncf[ 36 * kd + 2 ] = 1; /*O */

    /*C2H3CHO */
    ncf[ 37 * kd + 1 ] = 3; /*C */
    ncf[ 37 * kd + 0 ] = 4; /*H */
    ncf[ 37 * kd + 2 ] = 1; /*O */

    /*C3H5-A */
    ncf[ 38 * kd + 1 ] = 3; /*C */
    ncf[ 38 * kd + 0 ] = 5; /*H */

    /*C3H3 */
    ncf[ 39 * kd + 1 ] = 3; /*C */
    ncf[ 39 * kd + 0 ] = 3; /*H */

    /*NC3H7CHO */
    ncf[ 40 * kd + 1 ] = 4; /*C */
    ncf[ 40 * kd + 0 ] = 8; /*H */
    ncf[ 40 * kd + 2 ] = 1; /*O */

    /*C2H5COCH2 */
    ncf[ 41 * kd + 1 ] = 4; /*C */
    ncf[ 41 * kd + 0 ] = 7; /*H */
    ncf[ 41 * kd + 2 ] = 1; /*O */

    /*CH3CHCO */
    ncf[ 42 * kd + 1 ] = 3; /*C */
    ncf[ 42 * kd + 0 ] = 4; /*H */
    ncf[ 42 * kd + 2 ] = 1; /*O */

    /*NC3H7COCH3 */
    ncf[ 43 * kd + 1 ] = 5; /*C */
    ncf[ 43 * kd + 0 ] = 10; /*H */
    ncf[ 43 * kd + 2 ] = 1; /*O */

    /*NC3H7COCH2 */
    ncf[ 44 * kd + 1 ] = 5; /*C */
    ncf[ 44 * kd + 0 ] = 9; /*H */
    ncf[ 44 * kd + 2 ] = 1; /*O */

    /*NC7H16 */
    ncf[ 45 * kd + 1 ] = 7; /*C */
    ncf[ 45 * kd + 0 ] = 16; /*H */

    /*C7H15-1 */
    ncf[ 46 * kd + 1 ] = 7; /*C */
    ncf[ 46 * kd + 0 ] = 15; /*H */

    /*C7H15O2 */
    ncf[ 47 * kd + 1 ] = 7; /*C */
    ncf[ 47 * kd + 0 ] = 15; /*H */
    ncf[ 47 * kd + 2 ] = 2; /*O */

    /*C7H14OOHO2 */
    ncf[ 48 * kd + 1 ] = 7; /*C */
    ncf[ 48 * kd + 0 ] = 15; /*H */
    ncf[ 48 * kd + 2 ] = 4; /*O */

    /*C7H14O */
    ncf[ 49 * kd + 1 ] = 7; /*C */
    ncf[ 49 * kd + 0 ] = 14; /*H */
    ncf[ 49 * kd + 2 ] = 1; /*O */

    /*NC7KET */
    ncf[ 50 * kd + 1 ] = 7; /*C */
    ncf[ 50 * kd + 0 ] = 14; /*H */
    ncf[ 50 * kd + 2 ] = 3; /*O */

    /*N2 */
    ncf[ 51 * kd + 3 ] = 2; /*N */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(4);
    ename[0] = "H";
    ename[1] = "C";
    ename[2] = "O";
    ename[3] = "N";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(52);
    kname[0] = "H";
    kname[1] = "H2";
    kname[2] = "O";
    kname[3] = "O2";
    kname[4] = "OH";
    kname[5] = "H2O";
    kname[6] = "CO";
    kname[7] = "CO2";
    kname[8] = "CH3";
    kname[9] = "CH4";
    kname[10] = "HO2";
    kname[11] = "H2O2";
    kname[12] = "CH2O";
    kname[13] = "C2H6";
    kname[14] = "C2H4";
    kname[15] = "C2H5";
    kname[16] = "C2H";
    kname[17] = "C2H2";
    kname[18] = "CH3OH";
    kname[19] = "CH2CO";
    kname[20] = "HCCO";
    kname[21] = "CH3CHO";
    kname[22] = "C3H4-A";
    kname[23] = "C3H4-P";
    kname[24] = "C3H6";
    kname[25] = "C4H6";
    kname[26] = "C4H7";
    kname[27] = "C4H8-1";
    kname[28] = "PC4H9";
    kname[29] = "CH3COCH2";
    kname[30] = "C2H5CHO";
    kname[31] = "C5H9";
    kname[32] = "C5H10-1";
    kname[33] = "C5H11-1";
    kname[34] = "CH3O2";
    kname[35] = "CH3O2H";
    kname[36] = "C2H3CO";
    kname[37] = "C2H3CHO";
    kname[38] = "C3H5-A";
    kname[39] = "C3H3";
    kname[40] = "NC3H7CHO";
    kname[41] = "C2H5COCH2";
    kname[42] = "CH3CHCO";
    kname[43] = "NC3H7COCH3";
    kname[44] = "NC3H7COCH2";
    kname[45] = "NC7H16";
    kname[46] = "C7H15-1";
    kname[47] = "C7H15O2";
    kname[48] = "C7H14OOHO2";
    kname[49] = "C7H14O";
    kname[50] = "NC7KET";
    kname[51] = "N2";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)
{
    amrex::GpuArray<amrex::Real,2809> Jac = {0.0};
    amrex::GpuArray<amrex::Real,52> conc = {0.0};
    for (int n=0; n<52; n++) {
        conc[n] = 1.0/ 52.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<53; k++) {
        for (int l=0; l<53; l++) {
            if(Jac[ 53 * k + l] != 0.0){
                nJdata_tmp = nJdata_tmp + 1;
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;
}



/*compute the sparsity pattern of the system Jacobian */
void SPARSITY_INFO_SYST( int * nJdata, const int * consP, int NCELLS)
{
    amrex::GpuArray<amrex::Real,2809> Jac = {0.0};
    amrex::GpuArray<amrex::Real,52> conc = {0.0};
    for (int n=0; n<52; n++) {
        conc[n] = 1.0/ 52.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<53; k++) {
        for (int l=0; l<53; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(Jac[ 53 * k + l] != 0.0){
                    nJdata_tmp = nJdata_tmp + 1;
                }
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;
}



/*compute the sparsity pattern of the simplified (for preconditioning) system Jacobian */
void SPARSITY_INFO_SYST_SIMPLIFIED( int * nJdata, const int * consP)
{
    amrex::GpuArray<amrex::Real,2809> Jac = {0.0};
    amrex::GpuArray<amrex::Real,52> conc = {0.0};
    for (int n=0; n<52; n++) {
        conc[n] = 1.0/ 52.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<53; k++) {
        for (int l=0; l<53; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(Jac[ 53 * k + l] != 0.0){
                    nJdata_tmp = nJdata_tmp + 1;
                }
            }
        }
    }

    nJdata[0] = nJdata_tmp;
}


/*compute the sparsity pattern of the chemistry Jacobian in CSC format -- base 0 */
void SPARSITY_PREPROC_CSC(int *  rowVals, int *  colPtrs, const int * consP, int NCELLS)
{
    amrex::GpuArray<amrex::Real,2809> Jac = {0.0};
    amrex::GpuArray<amrex::Real,52> conc = {0.0};
    for (int n=0; n<52; n++) {
        conc[n] = 1.0/ 52.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc=0; nc<NCELLS; nc++) {
        int offset_row = nc * 53;
        int offset_col = nc * 53;
        for (int k=0; k<53; k++) {
            for (int l=0; l<53; l++) {
                if(Jac[53*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l + offset_row; 
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
            colPtrs[offset_col + (k + 1)] = nJdata_tmp;
        }
    }
}

/*compute the sparsity pattern of the chemistry Jacobian in CSR format -- base 0 */
void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, const int * consP, int NCELLS, int base)
{
    amrex::GpuArray<amrex::Real,2809> Jac = {0.0};
    amrex::GpuArray<amrex::Real,52> conc = {0.0};
    for (int n=0; n<52; n++) {
        conc[n] = 1.0/ 52.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtrs[0] = 1;
        int nJdata_tmp = 1;
        for (int nc=0; nc<NCELLS; nc++) {
            int offset = nc * 53;
            for (int l=0; l<53; l++) {
                for (int k=0; k<53; k++) {
                    if(Jac[53*k + l] != 0.0) {
                        colVals[nJdata_tmp-1] = k+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    }
                }
                rowPtrs[offset + (l + 1)] = nJdata_tmp;
            }
        }
    } else {
        rowPtrs[0] = 0;
        int nJdata_tmp = 0;
        for (int nc=0; nc<NCELLS; nc++) {
            int offset = nc * 53;
            for (int l=0; l<53; l++) {
                for (int k=0; k<53; k++) {
                    if(Jac[53*k + l] != 0.0) {
                        colVals[nJdata_tmp] = k + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    }
                }
                rowPtrs[offset + (l + 1)] = nJdata_tmp;
            }
        }
    }
}

/*compute the sparsity pattern of the system Jacobian */
/*CSR format BASE is user choice */
void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtr, const int * consP, int NCELLS, int base)
{
    amrex::GpuArray<amrex::Real,2809> Jac = {0.0};
    amrex::GpuArray<amrex::Real,52> conc = {0.0};
    for (int n=0; n<52; n++) {
        conc[n] = 1.0/ 52.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int nc=0; nc<NCELLS; nc++) {
            int offset = nc * 53;
            for (int l=0; l<53; l++) {
                for (int k=0; k<53; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(Jac[53*k + l] != 0.0) {
                            colVals[nJdata_tmp-1] = k+1 + offset; 
                            nJdata_tmp = nJdata_tmp + 1; 
                        }
                    }
                }
                rowPtr[offset + (l + 1)] = nJdata_tmp;
            }
        }
    } else {
        rowPtr[0] = 0;
        int nJdata_tmp = 0;
        for (int nc=0; nc<NCELLS; nc++) {
            int offset = nc * 53;
            for (int l=0; l<53; l++) {
                for (int k=0; k<53; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(Jac[53*k + l] != 0.0) {
                            colVals[nJdata_tmp] = k + offset; 
                            nJdata_tmp = nJdata_tmp + 1; 
                        }
                    }
                }
                rowPtr[offset + (l + 1)] = nJdata_tmp;
            }
        }
    }
}

/*compute the sparsity pattern of the simplified (for precond) system Jacobian on CPU */
/*BASE 0 */
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, const int * consP)
{
    amrex::GpuArray<amrex::Real,2809> Jac = {0.0};
    amrex::GpuArray<amrex::Real,52> conc = {0.0};
    for (int n=0; n<52; n++) {
        conc[n] = 1.0/ 52.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int k=0; k<53; k++) {
        for (int l=0; l<53; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 53*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(Jac[53*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 53*k + l;
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
        }
        colPtrs[k+1] = nJdata_tmp;
    }
}

/*compute the sparsity pattern of the simplified (for precond) system Jacobian */
/*CSR format BASE is under choice */
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, const int * consP, int base)
{
    amrex::GpuArray<amrex::Real,2809> Jac = {0.0};
    amrex::GpuArray<amrex::Real,52> conc = {0.0};
    for (int n=0; n<52; n++) {
        conc[n] = 1.0/ 52.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int l=0; l<53; l++) {
            for (int k=0; k<53; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(Jac[53*k + l] != 0.0) {
                        colVals[nJdata_tmp-1] = k+1; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    }
                }
            }
            rowPtr[l+1] = nJdata_tmp;
        }
    } else {
        rowPtr[0] = 0;
        int nJdata_tmp = 0;
        for (int l=0; l<53; l++) {
            for (int k=0; k<53; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(Jac[53*k + l] != 0.0) {
                        colVals[nJdata_tmp] = k; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    }
                }
            }
            rowPtr[l+1] = nJdata_tmp;
        }
    }
}

#endif
