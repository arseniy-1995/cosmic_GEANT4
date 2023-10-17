/* yn6.F -- translated by f2c (version 20230428).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Common Block Declarations */

struct data1_1_ {
    doublereal cmix[12]	/* was [6][2] */, pmasm[18], pmasb[18], pgamm[18], 
	    pgamb[18];
    integer meso[18]	/* was [9][2] */;
};

#define data1_1 (*(struct data1_1_ *) &data1_)

struct norms_1_ {
    doublereal bitr, bitel, bitc, alpha, bitrm, bitelm, bitcm;
};

#define norms_1 (*(struct norms_1_ *) &norms_)

struct cutoff_1_ {
    doublereal ba, bb, binf, pow;
    integer isl, islfor;
};

#define cutoff_1 (*(struct cutoff_1_ *) &cutoff_)

struct const_1_ {
    doublereal pi, large;
};

#define const_1 (*(struct const_1_ *) &const_)

struct comcut_1_ {
    doublereal parm, parb, swmax;
};

#define comcut_1 (*(struct comcut_1_ *) &comcut_)

struct comabm_1_ {
    doublereal alfam, betam;
};

#define comabm_1 (*(struct comabm_1_ *) &comabm_)

struct comabb_1_ {
    doublereal alfab, betab;
};

#define comabb_1 (*(struct comabb_1_ *) &comabb_)

struct data7_1_ {
    doublereal cbr[35];
    integer kdp[105]	/* was [35][3] */, idc0[18]	/* was [9][2] */;
};

#define data7_1 (*(struct data7_1_ *) &data7_)

struct data2_1_ {
    doublereal pud, ps1, sigma, cx2;
};

#define data2_1 (*(struct data2_1_ *) &data2_)

struct data4_1_ {
    doublereal qmas[9];
};

#define data4_1 (*(struct data4_1_ *) &data4_)

struct data6_1_ {
    doublereal pod810[18]	/* was [3][6] */, podsa[24]	/* was [8][3] 
	    */;
    integer kbar[36]	/* was [18][2] */;
};

#define data6_1 (*(struct data6_1_ *) &data6_)

struct data5_1_ {
    doublereal dqq[9]	/* was [3][3] */;
    integer iflm[36]	/* was [18][2] */, iflb[54]	/* was [18][3] */;
};

#define data5_1 (*(struct data5_1_ *) &data5_)

struct comcha_1_ {
    integer icham[18], ichab[18];
};

#define comcha_1 (*(struct comcha_1_ *) &comcha_)

struct comstr_1_ {
    integer istr[36];
};

#define comstr_1 (*(struct comstr_1_ *) &comstr_)

struct data3_1_ {
    doublereal popb[10];
};

#define data3_1 (*(struct data3_1_ *) &data3_)

/* Initialized data */

struct {
    doublereal e_1[84];
    integer e_2[18];
    } data1_ = { .5, .5, 1., .5, .5, 1., .25, .25, .5, 0., 0., 1., .14, .14, 
	    .492, .492, .492, .492, .14, .549, .958, .77, .77, .896, .896, 
	    .896, .896, .77, .783, 1.02, .94, .94, 1.189, 1.189, 1.189, 1.315,
	     1.315, 1.116, 1.232, 1.232, 1.232, 1.232, 1.385, 1.385, 1.385, 
	    1.53, 1.53, 1.672, 0., 0., 0., 0., 0., 0., 0., 0., 0., .152, .152,
	     .049, .049, .049, .049, .152, .01, .004, 0., 0., 0., 0., 0., 0., 
	    0., 0., .122, .122, .122, .122, .035, .042, .035, .0091, .01, 0., 
	    7, 1, 3, 2, 8, 5, 4, 6, 9, 7, 2, 4, 1, 8, 6, 3, 5, 9 };

struct {
    doublereal e_1[35];
    integer e_2[123];
    } data7_ = { 1., 1., .667, 1., .667, 1., .667, 1., .667, 1., 1., .899, 1.,
	     .485, .837, 1., 1., .333, 1., 1., .667, 1., .88, .94, 1., .88, 
	    .94, 1., .88, .94, 1., .333, 1., .333, 1., 1, 2, 5, 3, 6, 4, 3, 5,
	     4, 6, 1, 1, 1, 3, 5, 1, 37, 38, 37, 38, 38, 37, 44, 41, 39, 44, 
	    40, 41, 44, 39, 40, 42, 43, 43, 42, 7, 7, 1, 7, 2, 7, 2, 7, 1, 7, 
	    2, 2, 2, 4, 6, 2, 1, 1, 7, 2, 7, 2, 1, 1, 7, 2, 7, 2, 7, 2, 1, 7, 
	    2, 7, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 7, 0, 0, 0, 
	    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 4, 6, 8, 
	    10, 11, 13, 16, 17, 19, 20, 22, 25, 28, 31, 33 };

struct {
    doublereal e_1[4];
    } data2_ = { .43, .75, .51, .77 };

struct {
    doublereal e_1[9];
    } data4_ = { 0., 0., 0., 0., 0., 0., 0., 0., 0. };

struct {
    doublereal e_1[42];
    integer e_2[36];
    } data6_ = { 1., .667, .667, .667, 1., .667, .89, .89, .667, .89, .667, 
	    .89, .667, .89, .89, .667, .667, 1., .25, 1., .25, 0., .25, 0., 
	    1., .75, 1., .25, 0., .25, .25, 1., 0., .75, 0., 0., 1., 1., 1., 
	    .25, .25, 1., 0, 37, 39, 38, 0, 40, 37, 38, 41, 39, 44, 43, 44, 
	    40, 42, 43, 42, 0, 45, 46, 49, 48, 47, 50, 46, 48, 51, 49, 51, 53,
	     51, 50, 52, 53, 52, 54 };

struct {
    doublereal e_1[9];
    integer e_2[90];
    } data5_ = { 4., 6., 7., 6., 5., 8., 7., 8., 9., 1, 2, 1, 3, 2, 3, 1, 1, 
	    1, 1, 2, 1, 3, 2, 3, 1, 1, 3, -2, -1, -3, -1, -3, -2, -1, -1, -1, 
	    -2, -1, -3, -1, -3, -2, -1, -1, -3, 1, 1, 1, 2, 1, 2, 1, 1, 1, 1, 
	    2, 1, 1, 2, 1, 2, 1, 3, 1, 2, 1, 2, 2, 3, 3, 2, 1, 1, 2, 2, 1, 2, 
	    2, 3, 3, 3, 2, 2, 3, 3, 3, 3, 3, 3, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3 }
	    ;

struct {
    integer e_1[36];
    } comcha_ = { 1, -1, 1, -1, 0, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 0, 0, 0, 1,
	     0, 1, -1, 0, -1, 0, 0, 2, 1, -1, 0, 1, -1, 0, -1, 0, -1 };

struct {
    integer e_1[36];
    } comstr_ = { 0, 0, 1, -1, 1, -1, 0, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 0, 0,
	     0, -1, -1, -1, -2, -2, -1, 0, 0, 0, 0, -1, -1, -1, -2, -2, -3 };

struct {
    doublereal e_1[10];
    } data3_ = { 1., .89, 1., .89, .89, .89, .67, .67, .67, 1. };

struct {
    doublereal e_1[7];
    } norms_ = { 44.404, 45.075, 61.5, .5, .11404, .0045075, 61.5 };

struct {
    doublereal e_1[4];
    integer e_2[2];
    } cutoff_ = { 7.26, 7.26, 4., 2., 2, 1 };

struct {
    doublereal e_1[2];
    } const_ = { 3.14159, 1e4 };

struct {
    doublereal e_1[3];
    } comcut_ = { 1.4, 2., .35 };

struct {
    doublereal e_1[2];
    } comabm_ = { .5, 2.5 };

struct {
    doublereal e_1[2];
    } comabb_ = { .5, 3. };


/*<       real FUNCTION AMAS(IK) >*/
doublereal amas_(integer *ik)
{
    /* System generated locals */
    real ret_val;

/*  DETERMINATION OF PARTICLE MASS */
/*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
/*      REAL   AMAS */
/*<        >*/
/*<       IF(IK.EQ.55)GO TO 3 >*/
    if (*ik == 55) {
	goto L3;
    }
/*<       IF(IK-36 < 0) THEN >*/
    if (*ik - 36 < 0) {
/*<       GOTO 1 >*/
	goto L1;
/*<       ELSE IF (IK-36  == 0) THEN >*/
    } else if (*ik - 36 == 0) {
/*<       GOTO 1 >*/
	goto L1;
/*<       ELSE IF (IK-36  > 0) THEN >*/
    } else if (*ik - 36 > 0) {
/*<       GOTO 2 >*/
	goto L2;
/*<       END IF >*/
    }
/*<  1    AMAS=PMASM(IK) >*/
L1:
    ret_val = data1_1.pmasm[*ik - 1];
/*<       RETURN >*/
    return ret_val;
/*<  2    AMAS=PMASB(IK-36) >*/
L2:
    ret_val = data1_1.pmasb[*ik - 37];
/*<       RETURN >*/
    return ret_val;
/*<  3    AMAS=0. >*/
L3:
    ret_val = 0.f;
/*<       RETURN >*/
    return ret_val;
/*<       END >*/
} /* amas_ */

/*<       SUBROUTINE TC6			! ?????????????? >*/
/* Subroutine */ int tc6_(void)
{
/*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
/*<         REAL*8 LARGE >*/
/*<       COMMON /NORMS/ BITR,BITEL,BITC,ALPHA,BITRM,BITELM,BITCM >*/
/*<       COMMON /CUTOFF/BA,BB,BINF,POW,ISL,ISLFOR >*/
/*<         COMMON /CONST/ PI,LARGE >*/
/*<       COMMON /COMCUT/ PARM,PARB,SWMAX >*/

/*<       DATA PI/3.14159/, LARGE /10000./ >*/
/*<       DATA PARM/1.4/,PARB/2.0/,SWMAX/0.35/ >*/

/*   PARAMETERS FIXING THE SIZE OF ELASTIC AND RESONANT EXCITATIONS OF A */
/*   PROTON ARE BITEL,BITR,BITC */
/*<       DATA BITEL/45.075/ >*/
/*<       DATA BITR/44.404/ >*/
/*<       DATA BITC/61.50/ >*/
/*<       DATA BITELM/0.0045075/ >*/
/*<       DATA BITRM/0.11404/ >*/
/*<       DATA BITCM/61.50/ >*/
/*   ALPHA=TRAJECTORY INTERCEPT FOR HIGH MASS EXCITATIONS */
/*<       DATA ALPHA/0.5/ >*/
/*   BA,BB FIX THE ELASTIC DIFFRACTIVE SLOPE AT INFINITE ENERGY */
/*<       DATA BA/7.26/ >*/
/*<       DATA BB/7.26/ >*/
/*   BINF FIXES THE DIFFRACTIVE SLOPE FOR LARGE MASS CLUSTERS */
/*<       DATA BINF/4.0/ >*/
/*<       DATA POW/2.0/ >*/
/*   ISL=1,2,3 FIXES THE TYPE OF MASS DEPENDENCE OF THE DIFFRACTIVE SLOPE */
/*<        DATA ISL/2/ >*/
/*   ISLFOR=1 OR 2 DETERMINES THE PRESENCE OR ABSENCE OF A T(FORWARD) */
/*   CUTOFF ON THE MASS SPECTRUM */
/*<       DATA ISLFOR/1/ >*/

/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* tc6_ */

/*<       SUBROUTINE TC7 >*/
/* Subroutine */ int tc7_(void)
{
/*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
/*<       COMMON /COMABM/ ALFAM,BETAM >*/
/*<       COMMON /COMABB/ ALFAB,BETAB >*/
/*<       COMMON /DATA7/ CBR(35),KDP(35,3),IDC0(9,2) >*/
/*<       COMMON /DATA2/PUD,PS1,SIGMA,CX2 >*/
/*<        >*/
/*<       COMMON /DATA4/QMAS(9) >*/
/*<       COMMON /DATA6/POD810(3,6),PODSA(8,3),KBAR(18,2) >*/
/*<       COMMON /DATA5/DQQ(3,3),IFLM(18,2),IFLB(18,3) >*/
/*<       COMMON /COMCHA/ ICHAM(18),ICHAB(18) >*/
/*<       COMMON /COMSTR/ ISTR(36) >*/
/*<       COMMON /DATA3/ POPB(10) >*/
/*<       DATA POPB/1.0,0.89,1.0,0.89,0.89,0.89,0.67,0.67,0.67,1.0/ >*/
/*<        >*/
/*<        >*/
/*<        >*/
/*<        >*/
/*<       DATA ICHAM/1,-1,1,-1,0,0,0,0,0,1,-1,1,-1,0,0,0,0,0/ >*/
/*<       DATA ICHAB/1,0,1,-1,0,-1,0,0,2,1,-1,0,1,-1,0,-1,0,-1/ >*/
/*<       DATA ALFAM/0.5/,BETAM/2.5/,ALFAB/0.5/,BETAB/3.0/ >*/
/*<       DATA PUD/0.43/,PS1/0.75/,SIGMA/0.51/,CX2/0.77/ >*/
/*<       DATA MESO/7,1,3,2,8,5,4,6,9,7,2,4,1,8,6,3,5,9/ >*/
/*<       DATA CMIX/2*0.5,1.,2*0.5,1.,2*0.25,0.5,2*0.,1./ >*/
/*<        >*/
/*<        >*/
/*<       DATA PGAMM/9*0.,2*0.152,4*0.049,0.152,0.01,0.004/ >*/
/*<       DATA PGAMB/8*0.,4*0.122,0.035,0.042,0.035,0.0091,0.01,0./ >*/
/*<        >*/
/*<       DATA QMAS/9*0./ >*/
/*<        >*/
/*<        >*/
/*<        >*/
/*<        >*/
/*<        >*/
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* tc7_ */

