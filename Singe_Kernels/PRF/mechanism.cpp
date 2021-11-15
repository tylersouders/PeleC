#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "mechanism.H"



/*save atomic weights into array */
void atomicWeight(amrex::Real *  awt)
{
    awt[0] = 12.011150; /*C */
    awt[1] = 1.007970; /*H */
    awt[2] = 14.006700; /*N */
    awt[3] = 15.999400; /*O */
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
    for (id = 0; id < kd * 116; ++ id) {
         ncf[id] = 0; 
    }

    /*H */
    ncf[ 0 * kd + 1 ] = 1; /*H */

    /*H2 */
    ncf[ 1 * kd + 1 ] = 2; /*H */

    /*O */
    ncf[ 2 * kd + 3 ] = 1; /*O */

    /*O2 */
    ncf[ 3 * kd + 3 ] = 2; /*O */

    /*OH */
    ncf[ 4 * kd + 1 ] = 1; /*H */
    ncf[ 4 * kd + 3 ] = 1; /*O */

    /*H2O */
    ncf[ 5 * kd + 1 ] = 2; /*H */
    ncf[ 5 * kd + 3 ] = 1; /*O */

    /*HO2 */
    ncf[ 6 * kd + 1 ] = 1; /*H */
    ncf[ 6 * kd + 3 ] = 2; /*O */

    /*H2O2 */
    ncf[ 7 * kd + 1 ] = 2; /*H */
    ncf[ 7 * kd + 3 ] = 2; /*O */

    /*CO */
    ncf[ 8 * kd + 0 ] = 1; /*C */
    ncf[ 8 * kd + 3 ] = 1; /*O */

    /*CO2 */
    ncf[ 9 * kd + 0 ] = 1; /*C */
    ncf[ 9 * kd + 3 ] = 2; /*O */

    /*CH2O */
    ncf[ 10 * kd + 0 ] = 1; /*C */
    ncf[ 10 * kd + 1 ] = 2; /*H */
    ncf[ 10 * kd + 3 ] = 1; /*O */

    /*HCO */
    ncf[ 11 * kd + 1 ] = 1; /*H */
    ncf[ 11 * kd + 0 ] = 1; /*C */
    ncf[ 11 * kd + 3 ] = 1; /*O */

    /*HOCHO */
    ncf[ 12 * kd + 0 ] = 1; /*C */
    ncf[ 12 * kd + 1 ] = 2; /*H */
    ncf[ 12 * kd + 3 ] = 2; /*O */

    /*CH3OH */
    ncf[ 13 * kd + 0 ] = 1; /*C */
    ncf[ 13 * kd + 1 ] = 4; /*H */
    ncf[ 13 * kd + 3 ] = 1; /*O */

    /*CH3O2H */
    ncf[ 14 * kd + 0 ] = 1; /*C */
    ncf[ 14 * kd + 1 ] = 4; /*H */
    ncf[ 14 * kd + 3 ] = 2; /*O */

    /*CH3O2 */
    ncf[ 15 * kd + 0 ] = 1; /*C */
    ncf[ 15 * kd + 1 ] = 3; /*H */
    ncf[ 15 * kd + 3 ] = 2; /*O */

    /*CH4 */
    ncf[ 16 * kd + 1 ] = 4; /*H */
    ncf[ 16 * kd + 0 ] = 1; /*C */

    /*CH3 */
    ncf[ 17 * kd + 0 ] = 1; /*C */
    ncf[ 17 * kd + 1 ] = 3; /*H */

    /*C2H6 */
    ncf[ 18 * kd + 0 ] = 2; /*C */
    ncf[ 18 * kd + 1 ] = 6; /*H */

    /*C2H5 */
    ncf[ 19 * kd + 0 ] = 2; /*C */
    ncf[ 19 * kd + 1 ] = 5; /*H */

    /*C2H4 */
    ncf[ 20 * kd + 0 ] = 2; /*C */
    ncf[ 20 * kd + 1 ] = 4; /*H */

    /*C2H3 */
    ncf[ 21 * kd + 0 ] = 2; /*C */
    ncf[ 21 * kd + 1 ] = 3; /*H */

    /*C2H2 */
    ncf[ 22 * kd + 0 ] = 2; /*C */
    ncf[ 22 * kd + 1 ] = 2; /*H */

    /*CH3CHO */
    ncf[ 23 * kd + 0 ] = 2; /*C */
    ncf[ 23 * kd + 1 ] = 4; /*H */
    ncf[ 23 * kd + 3 ] = 1; /*O */

    /*CH2CHO */
    ncf[ 24 * kd + 0 ] = 2; /*C */
    ncf[ 24 * kd + 1 ] = 3; /*H */
    ncf[ 24 * kd + 3 ] = 1; /*O */

    /*CH2CO */
    ncf[ 25 * kd + 0 ] = 2; /*C */
    ncf[ 25 * kd + 1 ] = 2; /*H */
    ncf[ 25 * kd + 3 ] = 1; /*O */

    /*HCCO */
    ncf[ 26 * kd + 1 ] = 1; /*H */
    ncf[ 26 * kd + 0 ] = 2; /*C */
    ncf[ 26 * kd + 3 ] = 1; /*O */

    /*CH3CO3 */
    ncf[ 27 * kd + 0 ] = 2; /*C */
    ncf[ 27 * kd + 1 ] = 3; /*H */
    ncf[ 27 * kd + 3 ] = 3; /*O */

    /*C2H5O2H */
    ncf[ 28 * kd + 0 ] = 2; /*C */
    ncf[ 28 * kd + 1 ] = 6; /*H */
    ncf[ 28 * kd + 3 ] = 2; /*O */

    /*C2H5O2 */
    ncf[ 29 * kd + 0 ] = 2; /*C */
    ncf[ 29 * kd + 1 ] = 5; /*H */
    ncf[ 29 * kd + 3 ] = 2; /*O */

    /*CH3COCH3 */
    ncf[ 30 * kd + 0 ] = 3; /*C */
    ncf[ 30 * kd + 1 ] = 6; /*H */
    ncf[ 30 * kd + 3 ] = 1; /*O */

    /*CH3COCH2O2 */
    ncf[ 31 * kd + 0 ] = 3; /*C */
    ncf[ 31 * kd + 1 ] = 5; /*H */
    ncf[ 31 * kd + 3 ] = 3; /*O */

    /*CH3COCH2O2H */
    ncf[ 32 * kd + 0 ] = 3; /*C */
    ncf[ 32 * kd + 1 ] = 6; /*H */
    ncf[ 32 * kd + 3 ] = 3; /*O */

    /*C2H3CHO */
    ncf[ 33 * kd + 0 ] = 3; /*C */
    ncf[ 33 * kd + 1 ] = 4; /*H */
    ncf[ 33 * kd + 3 ] = 1; /*O */

    /*C2H3CO */
    ncf[ 34 * kd + 0 ] = 3; /*C */
    ncf[ 34 * kd + 1 ] = 3; /*H */
    ncf[ 34 * kd + 3 ] = 1; /*O */

    /*C2H5CHO */
    ncf[ 35 * kd + 0 ] = 3; /*C */
    ncf[ 35 * kd + 1 ] = 6; /*H */
    ncf[ 35 * kd + 3 ] = 1; /*O */

    /*C3H6 */
    ncf[ 36 * kd + 0 ] = 3; /*C */
    ncf[ 36 * kd + 1 ] = 6; /*H */

    /*C3H5-A */
    ncf[ 37 * kd + 0 ] = 3; /*C */
    ncf[ 37 * kd + 1 ] = 5; /*H */

    /*C3H4-P */
    ncf[ 38 * kd + 1 ] = 4; /*H */
    ncf[ 38 * kd + 0 ] = 3; /*C */

    /*C3H4-A */
    ncf[ 39 * kd + 1 ] = 4; /*H */
    ncf[ 39 * kd + 0 ] = 3; /*C */

    /*C3H3 */
    ncf[ 40 * kd + 0 ] = 3; /*C */
    ncf[ 40 * kd + 1 ] = 3; /*H */

    /*NC3H7O2 */
    ncf[ 41 * kd + 0 ] = 3; /*C */
    ncf[ 41 * kd + 1 ] = 7; /*H */
    ncf[ 41 * kd + 3 ] = 2; /*O */

    /*IC3H7O2 */
    ncf[ 42 * kd + 0 ] = 3; /*C */
    ncf[ 42 * kd + 1 ] = 7; /*H */
    ncf[ 42 * kd + 3 ] = 2; /*O */

    /*CH3CHCO */
    ncf[ 43 * kd + 0 ] = 3; /*C */
    ncf[ 43 * kd + 1 ] = 4; /*H */
    ncf[ 43 * kd + 3 ] = 1; /*O */

    /*C4H8-1 */
    ncf[ 44 * kd + 0 ] = 4; /*C */
    ncf[ 44 * kd + 1 ] = 8; /*H */

    /*SC4H9 */
    ncf[ 45 * kd + 0 ] = 4; /*C */
    ncf[ 45 * kd + 1 ] = 9; /*H */

    /*C4H71-3 */
    ncf[ 46 * kd + 0 ] = 4; /*C */
    ncf[ 46 * kd + 1 ] = 7; /*H */

    /*C4H6 */
    ncf[ 47 * kd + 0 ] = 4; /*C */
    ncf[ 47 * kd + 1 ] = 6; /*H */

    /*PC4H9O2 */
    ncf[ 48 * kd + 0 ] = 4; /*C */
    ncf[ 48 * kd + 1 ] = 9; /*H */
    ncf[ 48 * kd + 3 ] = 2; /*O */

    /*C2H5COCH2 */
    ncf[ 49 * kd + 0 ] = 4; /*C */
    ncf[ 49 * kd + 1 ] = 7; /*H */
    ncf[ 49 * kd + 3 ] = 1; /*O */

    /*NC3H7CHO */
    ncf[ 50 * kd + 0 ] = 4; /*C */
    ncf[ 50 * kd + 1 ] = 8; /*H */
    ncf[ 50 * kd + 3 ] = 1; /*O */

    /*IC4H8 */
    ncf[ 51 * kd + 0 ] = 4; /*C */
    ncf[ 51 * kd + 1 ] = 8; /*H */

    /*IC4H7 */
    ncf[ 52 * kd + 0 ] = 4; /*C */
    ncf[ 52 * kd + 1 ] = 7; /*H */

    /*TC4H9O2 */
    ncf[ 53 * kd + 0 ] = 4; /*C */
    ncf[ 53 * kd + 1 ] = 9; /*H */
    ncf[ 53 * kd + 3 ] = 2; /*O */

    /*IC4H9O2 */
    ncf[ 54 * kd + 0 ] = 4; /*C */
    ncf[ 54 * kd + 1 ] = 9; /*H */
    ncf[ 54 * kd + 3 ] = 2; /*O */

    /*IC3H7CHO */
    ncf[ 55 * kd + 0 ] = 4; /*C */
    ncf[ 55 * kd + 1 ] = 8; /*H */
    ncf[ 55 * kd + 3 ] = 1; /*O */

    /*TC3H6CHO */
    ncf[ 56 * kd + 0 ] = 4; /*C */
    ncf[ 56 * kd + 1 ] = 7; /*H */
    ncf[ 56 * kd + 3 ] = 1; /*O */

    /*IC4H8OOH-IO2 */
    ncf[ 57 * kd + 0 ] = 4; /*C */
    ncf[ 57 * kd + 1 ] = 9; /*H */
    ncf[ 57 * kd + 3 ] = 4; /*O */

    /*IC4KETII */
    ncf[ 58 * kd + 0 ] = 4; /*C */
    ncf[ 58 * kd + 1 ] = 8; /*H */
    ncf[ 58 * kd + 3 ] = 3; /*O */

    /*IC4H6OH */
    ncf[ 59 * kd + 0 ] = 4; /*C */
    ncf[ 59 * kd + 1 ] = 7; /*H */
    ncf[ 59 * kd + 3 ] = 1; /*O */

    /*IC3H5CHO */
    ncf[ 60 * kd + 0 ] = 4; /*C */
    ncf[ 60 * kd + 1 ] = 6; /*H */
    ncf[ 60 * kd + 3 ] = 1; /*O */

    /*IC3H5CO */
    ncf[ 61 * kd + 0 ] = 4; /*C */
    ncf[ 61 * kd + 1 ] = 5; /*H */
    ncf[ 61 * kd + 3 ] = 1; /*O */

    /*IC4H7OOH */
    ncf[ 62 * kd + 0 ] = 4; /*C */
    ncf[ 62 * kd + 1 ] = 8; /*H */
    ncf[ 62 * kd + 3 ] = 2; /*O */

    /*TC3H6O2CHO */
    ncf[ 63 * kd + 0 ] = 4; /*C */
    ncf[ 63 * kd + 1 ] = 7; /*H */
    ncf[ 63 * kd + 3 ] = 3; /*O */

    /*C5H10-1 */
    ncf[ 64 * kd + 0 ] = 5; /*C */
    ncf[ 64 * kd + 1 ] = 10; /*H */

    /*C5H91-3 */
    ncf[ 65 * kd + 0 ] = 5; /*C */
    ncf[ 65 * kd + 1 ] = 9; /*H */

    /*C5H91-4 */
    ncf[ 66 * kd + 0 ] = 5; /*C */
    ncf[ 66 * kd + 1 ] = 9; /*H */

    /*C6H13O2-1 */
    ncf[ 67 * kd + 0 ] = 6; /*C */
    ncf[ 67 * kd + 1 ] = 13; /*H */
    ncf[ 67 * kd + 3 ] = 2; /*O */

    /*NC4H9CHO */
    ncf[ 68 * kd + 0 ] = 5; /*C */
    ncf[ 68 * kd + 1 ] = 10; /*H */
    ncf[ 68 * kd + 3 ] = 1; /*O */

    /*NC7H16 */
    ncf[ 69 * kd + 0 ] = 7; /*C */
    ncf[ 69 * kd + 1 ] = 16; /*H */

    /*C7H14-1 */
    ncf[ 70 * kd + 0 ] = 7; /*C */
    ncf[ 70 * kd + 1 ] = 14; /*H */

    /*C7H14-2 */
    ncf[ 71 * kd + 0 ] = 7; /*C */
    ncf[ 71 * kd + 1 ] = 14; /*H */

    /*C7H14-3 */
    ncf[ 72 * kd + 0 ] = 7; /*C */
    ncf[ 72 * kd + 1 ] = 14; /*H */

    /*C7H132-4 */
    ncf[ 73 * kd + 0 ] = 7; /*C */
    ncf[ 73 * kd + 1 ] = 13; /*H */

    /*C7H15O2-1 */
    ncf[ 74 * kd + 0 ] = 7; /*C */
    ncf[ 74 * kd + 1 ] = 15; /*H */
    ncf[ 74 * kd + 3 ] = 2; /*O */

    /*C7H15O2-2 */
    ncf[ 75 * kd + 0 ] = 7; /*C */
    ncf[ 75 * kd + 1 ] = 15; /*H */
    ncf[ 75 * kd + 3 ] = 2; /*O */

    /*C7H15O2-3 */
    ncf[ 76 * kd + 0 ] = 7; /*C */
    ncf[ 76 * kd + 1 ] = 15; /*H */
    ncf[ 76 * kd + 3 ] = 2; /*O */

    /*C7H15O2-4 */
    ncf[ 77 * kd + 0 ] = 7; /*C */
    ncf[ 77 * kd + 1 ] = 15; /*H */
    ncf[ 77 * kd + 3 ] = 2; /*O */

    /*C7H14OOH1-3O2 */
    ncf[ 78 * kd + 0 ] = 7; /*C */
    ncf[ 78 * kd + 1 ] = 15; /*H */
    ncf[ 78 * kd + 3 ] = 4; /*O */

    /*C7H14OOH2-3O2 */
    ncf[ 79 * kd + 0 ] = 7; /*C */
    ncf[ 79 * kd + 1 ] = 15; /*H */
    ncf[ 79 * kd + 3 ] = 4; /*O */

    /*C7H14OOH2-4O2 */
    ncf[ 80 * kd + 0 ] = 7; /*C */
    ncf[ 80 * kd + 1 ] = 15; /*H */
    ncf[ 80 * kd + 3 ] = 4; /*O */

    /*C7H14OOH4-2O2 */
    ncf[ 81 * kd + 0 ] = 7; /*C */
    ncf[ 81 * kd + 1 ] = 15; /*H */
    ncf[ 81 * kd + 3 ] = 4; /*O */

    /*C7H14OOH4-3O2 */
    ncf[ 82 * kd + 0 ] = 7; /*C */
    ncf[ 82 * kd + 1 ] = 15; /*H */
    ncf[ 82 * kd + 3 ] = 4; /*O */

    /*C7H14O1-3 */
    ncf[ 83 * kd + 0 ] = 7; /*C */
    ncf[ 83 * kd + 1 ] = 14; /*H */
    ncf[ 83 * kd + 3 ] = 1; /*O */

    /*C7H14O2-4 */
    ncf[ 84 * kd + 0 ] = 7; /*C */
    ncf[ 84 * kd + 1 ] = 14; /*H */
    ncf[ 84 * kd + 3 ] = 1; /*O */

    /*C7H14O3-5 */
    ncf[ 85 * kd + 0 ] = 7; /*C */
    ncf[ 85 * kd + 1 ] = 14; /*H */
    ncf[ 85 * kd + 3 ] = 1; /*O */

    /*NC7KET13 */
    ncf[ 86 * kd + 0 ] = 7; /*C */
    ncf[ 86 * kd + 1 ] = 14; /*H */
    ncf[ 86 * kd + 3 ] = 3; /*O */

    /*C4H7OOH1-4 */
    ncf[ 87 * kd + 0 ] = 4; /*C */
    ncf[ 87 * kd + 1 ] = 8; /*H */
    ncf[ 87 * kd + 3 ] = 2; /*O */

    /*CH3CHCHO */
    ncf[ 88 * kd + 0 ] = 3; /*C */
    ncf[ 88 * kd + 1 ] = 5; /*H */
    ncf[ 88 * kd + 3 ] = 1; /*O */

    /*IC4H7-I1 */
    ncf[ 89 * kd + 0 ] = 4; /*C */
    ncf[ 89 * kd + 1 ] = 7; /*H */

    /*XC7H14 */
    ncf[ 90 * kd + 0 ] = 7; /*C */
    ncf[ 90 * kd + 1 ] = 14; /*H */

    /*YC7H14 */
    ncf[ 91 * kd + 0 ] = 7; /*C */
    ncf[ 91 * kd + 1 ] = 14; /*H */

    /*XC7H13-Z */
    ncf[ 92 * kd + 0 ] = 7; /*C */
    ncf[ 92 * kd + 1 ] = 13; /*H */

    /*YC7H13-Y2 */
    ncf[ 93 * kd + 0 ] = 7; /*C */
    ncf[ 93 * kd + 1 ] = 13; /*H */

    /*YC7H15O2 */
    ncf[ 94 * kd + 0 ] = 7; /*C */
    ncf[ 94 * kd + 1 ] = 15; /*H */
    ncf[ 94 * kd + 3 ] = 2; /*O */

    /*ACC6H10 */
    ncf[ 95 * kd + 0 ] = 6; /*C */
    ncf[ 95 * kd + 1 ] = 10; /*H */

    /*ACC6H9-A */
    ncf[ 96 * kd + 0 ] = 6; /*C */
    ncf[ 96 * kd + 1 ] = 9; /*H */

    /*ACC6H9-D */
    ncf[ 97 * kd + 0 ] = 6; /*C */
    ncf[ 97 * kd + 1 ] = 9; /*H */

    /*NEOC5H11 */
    ncf[ 98 * kd + 0 ] = 5; /*C */
    ncf[ 98 * kd + 1 ] = 11; /*H */

    /*TC4H9CHO */
    ncf[ 99 * kd + 0 ] = 5; /*C */
    ncf[ 99 * kd + 1 ] = 10; /*H */
    ncf[ 99 * kd + 3 ] = 1; /*O */

    /*IC8H18 */
    ncf[ 100 * kd + 0 ] = 8; /*C */
    ncf[ 100 * kd + 1 ] = 18; /*H */

    /*IC8H16 */
    ncf[ 101 * kd + 0 ] = 8; /*C */
    ncf[ 101 * kd + 1 ] = 16; /*H */

    /*JC8H16 */
    ncf[ 102 * kd + 0 ] = 8; /*C */
    ncf[ 102 * kd + 1 ] = 16; /*H */

    /*BC8H17O2 */
    ncf[ 103 * kd + 0 ] = 8; /*C */
    ncf[ 103 * kd + 1 ] = 17; /*H */
    ncf[ 103 * kd + 3 ] = 2; /*O */

    /*CC8H17O2 */
    ncf[ 104 * kd + 0 ] = 8; /*C */
    ncf[ 104 * kd + 1 ] = 17; /*H */
    ncf[ 104 * kd + 3 ] = 2; /*O */

    /*IC8ETERAB */
    ncf[ 105 * kd + 0 ] = 8; /*C */
    ncf[ 105 * kd + 1 ] = 16; /*H */
    ncf[ 105 * kd + 3 ] = 1; /*O */

    /*IC8ETERBD */
    ncf[ 106 * kd + 0 ] = 8; /*C */
    ncf[ 106 * kd + 1 ] = 16; /*H */
    ncf[ 106 * kd + 3 ] = 1; /*O */

    /*IC8KETDB */
    ncf[ 107 * kd + 0 ] = 8; /*C */
    ncf[ 107 * kd + 1 ] = 16; /*H */
    ncf[ 107 * kd + 3 ] = 3; /*O */

    /*iso002 */
    ncf[ 108 * kd + 0 ] = 8; /*C */
    ncf[ 108 * kd + 1 ] = 17; /*H */
    ncf[ 108 * kd + 3 ] = 4; /*O */

    /*iso003 */
    ncf[ 109 * kd + 0 ] = 8; /*C */
    ncf[ 109 * kd + 1 ] = 17; /*H */
    ncf[ 109 * kd + 3 ] = 4; /*O */

    /*iso004 */
    ncf[ 110 * kd + 0 ] = 8; /*C */
    ncf[ 110 * kd + 1 ] = 17; /*H */
    ncf[ 110 * kd + 3 ] = 2; /*O */

    /*iso008 */
    ncf[ 111 * kd + 0 ] = 8; /*C */
    ncf[ 111 * kd + 1 ] = 16; /*H */
    ncf[ 111 * kd + 3 ] = 3; /*O */

    /*iso010 */
    ncf[ 112 * kd + 0 ] = 7; /*C */
    ncf[ 112 * kd + 1 ] = 15; /*H */
    ncf[ 112 * kd + 3 ] = 4; /*O */

    /*iso011 */
    ncf[ 113 * kd + 0 ] = 7; /*C */
    ncf[ 113 * kd + 1 ] = 15; /*H */
    ncf[ 113 * kd + 3 ] = 4; /*O */

    /*iso014 */
    ncf[ 114 * kd + 0 ] = 7; /*C */
    ncf[ 114 * kd + 1 ] = 14; /*H */
    ncf[ 114 * kd + 3 ] = 3; /*O */

    /*N2 */
    ncf[ 115 * kd + 2 ] = 2; /*N */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(4);
    ename[0] = "C";
    ename[1] = "H";
    ename[2] = "N";
    ename[3] = "O";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(116);
    kname[0] = "H";
    kname[1] = "H2";
    kname[2] = "O";
    kname[3] = "O2";
    kname[4] = "OH";
    kname[5] = "H2O";
    kname[6] = "HO2";
    kname[7] = "H2O2";
    kname[8] = "CO";
    kname[9] = "CO2";
    kname[10] = "CH2O";
    kname[11] = "HCO";
    kname[12] = "HOCHO";
    kname[13] = "CH3OH";
    kname[14] = "CH3O2H";
    kname[15] = "CH3O2";
    kname[16] = "CH4";
    kname[17] = "CH3";
    kname[18] = "C2H6";
    kname[19] = "C2H5";
    kname[20] = "C2H4";
    kname[21] = "C2H3";
    kname[22] = "C2H2";
    kname[23] = "CH3CHO";
    kname[24] = "CH2CHO";
    kname[25] = "CH2CO";
    kname[26] = "HCCO";
    kname[27] = "CH3CO3";
    kname[28] = "C2H5O2H";
    kname[29] = "C2H5O2";
    kname[30] = "CH3COCH3";
    kname[31] = "CH3COCH2O2";
    kname[32] = "CH3COCH2O2H";
    kname[33] = "C2H3CHO";
    kname[34] = "C2H3CO";
    kname[35] = "C2H5CHO";
    kname[36] = "C3H6";
    kname[37] = "C3H5-A";
    kname[38] = "C3H4-P";
    kname[39] = "C3H4-A";
    kname[40] = "C3H3";
    kname[41] = "NC3H7O2";
    kname[42] = "IC3H7O2";
    kname[43] = "CH3CHCO";
    kname[44] = "C4H8-1";
    kname[45] = "SC4H9";
    kname[46] = "C4H71-3";
    kname[47] = "C4H6";
    kname[48] = "PC4H9O2";
    kname[49] = "C2H5COCH2";
    kname[50] = "NC3H7CHO";
    kname[51] = "IC4H8";
    kname[52] = "IC4H7";
    kname[53] = "TC4H9O2";
    kname[54] = "IC4H9O2";
    kname[55] = "IC3H7CHO";
    kname[56] = "TC3H6CHO";
    kname[57] = "IC4H8OOH-IO2";
    kname[58] = "IC4KETII";
    kname[59] = "IC4H6OH";
    kname[60] = "IC3H5CHO";
    kname[61] = "IC3H5CO";
    kname[62] = "IC4H7OOH";
    kname[63] = "TC3H6O2CHO";
    kname[64] = "C5H10-1";
    kname[65] = "C5H91-3";
    kname[66] = "C5H91-4";
    kname[67] = "C6H13O2-1";
    kname[68] = "NC4H9CHO";
    kname[69] = "NC7H16";
    kname[70] = "C7H14-1";
    kname[71] = "C7H14-2";
    kname[72] = "C7H14-3";
    kname[73] = "C7H132-4";
    kname[74] = "C7H15O2-1";
    kname[75] = "C7H15O2-2";
    kname[76] = "C7H15O2-3";
    kname[77] = "C7H15O2-4";
    kname[78] = "C7H14OOH1-3O2";
    kname[79] = "C7H14OOH2-3O2";
    kname[80] = "C7H14OOH2-4O2";
    kname[81] = "C7H14OOH4-2O2";
    kname[82] = "C7H14OOH4-3O2";
    kname[83] = "C7H14O1-3";
    kname[84] = "C7H14O2-4";
    kname[85] = "C7H14O3-5";
    kname[86] = "NC7KET13";
    kname[87] = "C4H7OOH1-4";
    kname[88] = "CH3CHCHO";
    kname[89] = "IC4H7-I1";
    kname[90] = "XC7H14";
    kname[91] = "YC7H14";
    kname[92] = "XC7H13-Z";
    kname[93] = "YC7H13-Y2";
    kname[94] = "YC7H15O2";
    kname[95] = "ACC6H10";
    kname[96] = "ACC6H9-A";
    kname[97] = "ACC6H9-D";
    kname[98] = "NEOC5H11";
    kname[99] = "TC4H9CHO";
    kname[100] = "IC8H18";
    kname[101] = "IC8H16";
    kname[102] = "JC8H16";
    kname[103] = "BC8H17O2";
    kname[104] = "CC8H17O2";
    kname[105] = "IC8ETERAB";
    kname[106] = "IC8ETERBD";
    kname[107] = "IC8KETDB";
    kname[108] = "iso002";
    kname[109] = "iso003";
    kname[110] = "iso004";
    kname[111] = "iso008";
    kname[112] = "iso010";
    kname[113] = "iso011";
    kname[114] = "iso014";
    kname[115] = "N2";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)
{
    amrex::GpuArray<amrex::Real,13689> Jac = {0.0};
    amrex::GpuArray<amrex::Real,116> conc = {0.0};
    for (int n=0; n<116; n++) {
        conc[n] = 1.0/ 116.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<117; k++) {
        for (int l=0; l<117; l++) {
            if(Jac[ 117 * k + l] != 0.0){
                nJdata_tmp = nJdata_tmp + 1;
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;
}



/*compute the sparsity pattern of the system Jacobian */
void SPARSITY_INFO_SYST( int * nJdata, const int * consP, int NCELLS)
{
    amrex::GpuArray<amrex::Real,13689> Jac = {0.0};
    amrex::GpuArray<amrex::Real,116> conc = {0.0};
    for (int n=0; n<116; n++) {
        conc[n] = 1.0/ 116.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<117; k++) {
        for (int l=0; l<117; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(Jac[ 117 * k + l] != 0.0){
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
    amrex::GpuArray<amrex::Real,13689> Jac = {0.0};
    amrex::GpuArray<amrex::Real,116> conc = {0.0};
    for (int n=0; n<116; n++) {
        conc[n] = 1.0/ 116.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<117; k++) {
        for (int l=0; l<117; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(Jac[ 117 * k + l] != 0.0){
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
    amrex::GpuArray<amrex::Real,13689> Jac = {0.0};
    amrex::GpuArray<amrex::Real,116> conc = {0.0};
    for (int n=0; n<116; n++) {
        conc[n] = 1.0/ 116.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc=0; nc<NCELLS; nc++) {
        int offset_row = nc * 117;
        int offset_col = nc * 117;
        for (int k=0; k<117; k++) {
            for (int l=0; l<117; l++) {
                if(Jac[117*k + l] != 0.0) {
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
    amrex::GpuArray<amrex::Real,13689> Jac = {0.0};
    amrex::GpuArray<amrex::Real,116> conc = {0.0};
    for (int n=0; n<116; n++) {
        conc[n] = 1.0/ 116.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtrs[0] = 1;
        int nJdata_tmp = 1;
        for (int nc=0; nc<NCELLS; nc++) {
            int offset = nc * 117;
            for (int l=0; l<117; l++) {
                for (int k=0; k<117; k++) {
                    if(Jac[117*k + l] != 0.0) {
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
            int offset = nc * 117;
            for (int l=0; l<117; l++) {
                for (int k=0; k<117; k++) {
                    if(Jac[117*k + l] != 0.0) {
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
    amrex::GpuArray<amrex::Real,13689> Jac = {0.0};
    amrex::GpuArray<amrex::Real,116> conc = {0.0};
    for (int n=0; n<116; n++) {
        conc[n] = 1.0/ 116.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int nc=0; nc<NCELLS; nc++) {
            int offset = nc * 117;
            for (int l=0; l<117; l++) {
                for (int k=0; k<117; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(Jac[117*k + l] != 0.0) {
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
            int offset = nc * 117;
            for (int l=0; l<117; l++) {
                for (int k=0; k<117; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(Jac[117*k + l] != 0.0) {
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
    amrex::GpuArray<amrex::Real,13689> Jac = {0.0};
    amrex::GpuArray<amrex::Real,116> conc = {0.0};
    for (int n=0; n<116; n++) {
        conc[n] = 1.0/ 116.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int k=0; k<117; k++) {
        for (int l=0; l<117; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 117*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(Jac[117*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 117*k + l;
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
    amrex::GpuArray<amrex::Real,13689> Jac = {0.0};
    amrex::GpuArray<amrex::Real,116> conc = {0.0};
    for (int n=0; n<116; n++) {
        conc[n] = 1.0/ 116.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int l=0; l<117; l++) {
            for (int k=0; k<117; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(Jac[117*k + l] != 0.0) {
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
        for (int l=0; l<117; l++) {
            for (int k=0; k<117; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(Jac[117*k + l] != 0.0) {
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
