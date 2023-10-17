/* gamsu2.F -- translated by f2c (version 20230428).
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

union {
    struct {
	real pl[21];
    } _1;
    struct {
	real plagerr[21];
    } _2;
} lagerr_;

#define lagerr_1 (lagerr_._1)
#define lagerr_2 (lagerr_._2)

struct {
    char re_name__[45];
} re_name__;

#define re_name__1 re_name__

struct {
    logical eg, cmsf;
} logeg_;

#define logeg_1 logeg_

/* Table of constant values */

static integer c__6 = 6;
static integer c__9 = 9;
static integer c__14 = 14;
static integer c__19 = 19;
static integer c__26 = 26;
static integer c__33 = 33;
static integer c__42 = 42;
static integer c__1 = 1;

/*< 	REAL*4 FUNCTION TOP_GAMN(KIN_ENERGY,TOP_NUMBER) >*/
doublereal top_gamn__(real *kin_energy__, integer *top_number__)
{
    /* Format strings */
    static char fmt_1[] = "(\002 A WRONG TOPOLOGY NUMBER: \002,i4)";

    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern doublereal gamn_(real *, integer *);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 15, 0, fmt_1, 0 };


/* + */

/* FUNCTIONAL DESCRIPTION: */

/*    The function presents a calculation of gamma- */
/*    nucleon TOPOLOGICAL cross sections by means of approximation */
/*    of cross sections of certain channel and SU(2) statistical model */

/* CREATION DATE: */

/*     18.03.94       by Igor Pshenichnov */

/* FORMAL PARAMETERS: */

/*    KIN_ENERGY  - projectile kinetic energy (GeV) */

/*    TOP_NUMBER  - topology (number of particles in the final */
/*                  state (n=3,...,9)). */

/* COMMON BLOCKS: */

/*    none */

/* FUNCTION VALUE: */

/*    cross-section value (mb) */

/* - */
/*< 	REAL*4 KIN_ENERGY >*/
/*<         INTEGER*4 TOP_NUMBER >*/
/*<                 IF(TOP_NUMBER.EQ.3) THEN >*/
    if (*top_number__ == 3) {
/*<                       TOP_GAMN=GAMN(KIN_ENERGY,6)/0.7 >*/
	ret_val = gamn_(kin_energy__, &c__6) / .7f;
/*<                 ELSEIF(TOP_NUMBER.EQ.4) THEN >*/
    } else if (*top_number__ == 4) {
/*<                       TOP_GAMN=GAMN(KIN_ENERGY,9)/0.75 >*/
	ret_val = gamn_(kin_energy__, &c__9) / .75f;
/*<                 ELSEIF(TOP_NUMBER.EQ.5) THEN >*/
    } else if (*top_number__ == 5) {
/*<                       TOP_GAMN=GAMN(KIN_ENERGY,14)/0.2667 >*/
	ret_val = gamn_(kin_energy__, &c__14) / .2667f;
/*<                 ELSEIF(TOP_NUMBER.EQ.6) THEN >*/
    } else if (*top_number__ == 6) {
/*<                       TOP_GAMN=GAMN(KIN_ENERGY,19)/0.4381 >*/
	ret_val = gamn_(kin_energy__, &c__19) / .4381f;
/*<                 ELSEIF(TOP_NUMBER.EQ.7) THEN >*/
    } else if (*top_number__ == 7) {
/*<                       TOP_GAMN=GAMN(KIN_ENERGY,26)/0.125 >*/
	ret_val = gamn_(kin_energy__, &c__26) / .125f;
/*<                 ELSEIF(TOP_NUMBER.EQ.8) THEN >*/
    } else if (*top_number__ == 8) {
/*<                       TOP_GAMN=GAMN(KIN_ENERGY,33)/0.2755 >*/
	ret_val = gamn_(kin_energy__, &c__33) / .2755f;
/*<                 ELSEIF(TOP_NUMBER.EQ.9) THEN >*/
    } else if (*top_number__ == 9) {
/*<                       TOP_GAMN=GAMN(KIN_ENERGY,42)/0.05614 >*/
	ret_val = gamn_(kin_energy__, &c__42) / .05614f;
/*<                 ELSE >*/
    } else {
/*<                    WRITE(15,1) TOP_NUMBER >*/
	s_wsfe(&io___1);
	do_fio(&c__1, (char *)&(*top_number__), (ftnlen)sizeof(integer));
	e_wsfe();
/*<  1                 FORMAT(' A WRONG TOPOLOGY NUMBER: ',I4) >*/
/*<                    TOP_GAMN=0.0 >*/
	ret_val = 0.f;
/*<                    RETURN >*/
	return ret_val;
/*<                 ENDIF >*/
    }
/*<         RETURN >*/
    return ret_val;
/*< 	END >*/
} /* top_gamn__ */

/*< 	REAL*4 FUNCTION GAMN(KIN_ENERGY,TAB_NUMBER) >*/
doublereal gamn_(real *kin_energy__, integer *tab_number__)
{
    /* Initialized data */

    static char st_names__[45*7] = " 6: g p ==> p pi+ pi-                   "
	    "     " " 9: g p ==> p pi+ pi- pi0                    " "14: g p "
	    "==> p pi+ pi+ pi- pi-                " "19: g p ==> p pi+ pi+ pi"
	    "- pi- pi0            " "26: g p ==> p pi+ pi+ pi+ pi- pi- pi-   "
	    "     " "33: g p ==> p pi+ pi+ pi+ pi- pi- pi- pi0    " "42: g p "
	    "==> p pi+ pi+ pi+ pi+ pi- pi- pi- pi-";
    static real tk[7] = { 321.f,506.f,727.f,952.f,1215.f,1481.f,1788.f };
    static real a[35]	/* was [5][7] */ = { .33179f,-.08218f,.09976f,1.3e-4f,
	    0.f,.05379f,.13504f,-.03173f,0.f,0.f,-9.8e-4f,.098325f,-.0475448f,
	    0.f,0.f,.20568f,-.06299f,0.f,0.f,0.f,.0619f,-.01921f,0.f,0.f,0.f,
	    .11137f,-.04094f,0.f,0.f,0.f,.033678f,-.013025f,0.f,0.f,0.f };
    static integer max_power__[7] = { 3,2,2,1,1,1,1 };

    /* Format strings */
    static char fmt_100[] = "(\002   THIS IS NOT A BASE REACTION: \002,i4)";

    /* System generated locals */
    integer i__1;
    real ret_val, r__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double log(doublereal), exp(doublereal), pow_dd(doublereal *, doublereal *
	    );

    /* Local variables */
    static real c__, f;
    static integer n;
    static real x;
    static integer re_number__;
    extern doublereal al_(integer *, real *, real *);
    static real alv, alfa;

    /* Fortran I/O blocks */
    static cilist io___9 = { 0, 15, 0, fmt_100, 0 };


/* + */

/* FUNCTIONAL DESCRIPTION: */

/*    The function presents a calculation of approximation gamma- */
/*    nucleon partial cross-sections from experiments: */
/*    S.I.Alekhin et al., CERN-HERA 87-01 */

/* CREATION DATE: */

/*     16.03.94       by Igor Pshenichnov */

/* FORMAL PARAMETERS: */

/*    KIN_ENERGY  - projectile kinetic energy (GeV) */

/*    TAB_NUMBER  - reaction number (see table) */

/* COMMON BLOCKS: */

/*    /LAGERR/PL(0:20) - Lagerr polinom values */

/*    /RE_NAME/RE_NAME - reaction notation in output */

/* FUNCTION VALUE: */

/*    cross-section value (mb) */

/* - */
/*< 	REAL*4 PL(0:20),KIN_ENERGY,TK(7),A(0:4,7) >*/
/*<         REAL*4 ALFA,C_UNI,C_DPI >*/
/*<         INTEGER*4 TAB_NUMBER,RE_NUMBER,MAX_POWER(7) >*/
/*<         CHARACTER*45 RE_NAME,ST_NAMES(7) >*/
/*<         COMMON /LAGERR/PL >*/
/*< 	COMMON /RE_NAME/RE_NAME >*/
/* ............................................... */
/*< 	DATA  >*/
/*< 		DATA >*/
/*               .....................   6  ........................... */
/*< 	DATA  >*/
/*               .....................   9  ........................... */
/*               .....................  14  ........................... */
/*               .....................  19  ........................... */
/*               .....................  26  ........................... */
/*               .....................  33  ........................... */
/*               .....................  42  ........................... */
/*< 		DATA MAX_POWER/  3,  2,  2,  1,  1,  1,  1/     >*/
/*<                       ALFA=2.0 >*/
    alfa = 2.f;
/*<                       C=1.0 >*/
    c__ = 1.f;
/*<                 IF(TAB_NUMBER.EQ.6) THEN >*/
    if (*tab_number__ == 6) {
/*<                       RE_NUMBER=1 >*/
	re_number__ = 1;
/*<                 ELSEIF(TAB_NUMBER.EQ.9) THEN >*/
    } else if (*tab_number__ == 9) {
/*<                       RE_NUMBER=2 >*/
	re_number__ = 2;
/*<                 ELSEIF(TAB_NUMBER.EQ.14) THEN >*/
    } else if (*tab_number__ == 14) {
/*<                       RE_NUMBER=3 >*/
	re_number__ = 3;
/*<                 ELSEIF(TAB_NUMBER.EQ.19) THEN >*/
    } else if (*tab_number__ == 19) {
/*<                       RE_NUMBER=4 >*/
	re_number__ = 4;
/*<                 ELSEIF(TAB_NUMBER.EQ.26) THEN >*/
    } else if (*tab_number__ == 26) {
/*<                       RE_NUMBER=5 >*/
	re_number__ = 5;
/*<                 ELSEIF(TAB_NUMBER.EQ.33) THEN >*/
    } else if (*tab_number__ == 33) {
/*<                       RE_NUMBER=6 >*/
	re_number__ = 6;
/*<                 ELSEIF(TAB_NUMBER.EQ.42) THEN >*/
    } else if (*tab_number__ == 42) {
/*<                       RE_NUMBER=7 >*/
	re_number__ = 7;
/*<                 ELSE >*/
    } else {
/*<                    WRITE(15,100) TAB_NUMBER >*/
	s_wsfe(&io___9);
	do_fio(&c__1, (char *)&(*tab_number__), (ftnlen)sizeof(integer));
	e_wsfe();
/*<  100               FORMAT('   THIS IS NOT A BASE REACTION: ',I4) >*/
/*<                    GAMN=0.0 >*/
	ret_val = 0.f;
/*<                    RETURN >*/
	return ret_val;
/*<                 ENDIF >*/
    }
/*<                 RE_NAME=ST_NAMES(RE_NUMBER) >*/
    s_copy(re_name__1.re_name__, st_names__ + (re_number__ - 1) * 45, (ftnlen)
	    45, (ftnlen)45);
/*<                 IF(KIN_ENERGY*1000.0.LE.TK(RE_NUMBER)) THEN >*/
    if (*kin_energy__ * 1e3f <= tk[re_number__ - 1]) {
/*<                      GAMN=0.0 >*/
	ret_val = 0.f;
/*<                      RETURN >*/
	return ret_val;
/*<                 ENDIF >*/
    }
/*<         X=C*ALOG(KIN_ENERGY*1000./TK(RE_NUMBER))    >*/
    x = c__ * log(*kin_energy__ * 1e3f / tk[re_number__ - 1]);
/*< 	ALV=AL(MAX_POWER(RE_NUMBER),ALFA,X) >*/
    alv = al_(&max_power__[re_number__ - 1], &alfa, &x);
/*<                       F=0.0 >*/
    f = 0.f;
/*<         DO 1 N=0,MAX_POWER(RE_NUMBER) >*/
    i__1 = max_power__[re_number__ - 1];
    for (n = 0; n <= i__1; ++n) {
/*<         F=F+A(N,RE_NUMBER)*PL(N) >*/
	f += a[n + re_number__ * 5 - 5] * lagerr_1.pl[n];
/*<  1      CONTINUE >*/
/* L1: */
    }
/*<         F=F*(X**(ALFA/2.0))/EXP(X/2.0) >*/
    d__1 = (doublereal) x;
    d__2 = (doublereal) (alfa / 2.f);
    f = f * pow_dd(&d__1, &d__2) / exp(x / 2.f);
/*<         GAMN=F**2  >*/
/* Computing 2nd power */
    r__1 = f;
    ret_val = r__1 * r__1;
/*< 	RETURN  >*/
    return ret_val;
/*< 	END >*/
} /* gamn_ */

/*< 	REAL*4 FUNCTION AL(N,ALFA,X) >*/
doublereal al_(integer *n, real *alfa, real *x)
{
    /* System generated locals */
    integer i__1;
    real ret_val;

    /* Local variables */
    static integer m, index;

/* + */

/* FUNCTIONAL DESCRIPTION: */

/* The recurrent calculation of LAGERR's polinoms up to order=20 */
/* (See for details: Handbook of mathematical functions. Ed. by */
/* M.Abramowitz and I.A.Stegun) */

/* CREATION DATE: */

/*     08.02.93       by Igor Pshenichnov */

/* FORMAL PARAMETERS: */

/*    N    -  polinom's order */
/*    ALFA -  parameter */
/*    X    -  point to be calculated */

/* COMMON BLOCKS: */

/*    /LAGERR/PLAGERR(0:20) array contais polinom's values from */
/*                      zero oder to N */

/* FUNCTION VALUE: */

/*    AL   -  N-oder value of LAGERR's polinom */

/* - */
/*< 	COMMON /LAGERR/PLAGERR(0:20) >*/
/*<         DO INDEX=0,20 >*/
    for (index = 0; index <= 20; ++index) {
/*<           PLAGERR(INDEX)=0.0 >*/
	lagerr_2.plagerr[index] = 0.f;
/*<         END DO >*/
    }
/*< 	PLAGERR(0) = 1. >*/
    lagerr_2.plagerr[0] = 1.f;
/*<         AL = PLAGERR(0) >*/
    ret_val = lagerr_2.plagerr[0];
/*< 	IF (N.EQ.0) RETURN >*/
    if (*n == 0) {
	return ret_val;
    }

/*< 	PLAGERR(1) = ALFA+1-X >*/
    lagerr_2.plagerr[1] = *alfa + 1 - *x;
/*<         AL = PLAGERR(1) >*/
    ret_val = lagerr_2.plagerr[1];
/*<         IF (N.EQ.1) RETURN >*/
    if (*n == 1) {
	return ret_val;
    }

/*<         DO 1 M=1,N-1  >*/
    i__1 = *n - 1;
    for (m = 1; m <= i__1; ++m) {
/*<        >*/
	lagerr_2.plagerr[m + 1] = (((m << 1) + *alfa + 1 - *x) * 
		lagerr_2.plagerr[m] - (m + *alfa) * lagerr_2.plagerr[m - 1]) /
		 (m + 1);
/*<   1     CONTINUE >*/
/* L1: */
    }
/*<         AL = PLAGERR(N) >*/
    ret_val = lagerr_2.plagerr[*n];
/*<         RETURN >*/
    return ret_val;
/*< 	END >*/
} /* al_ */

/*<         SUBROUTINE GAMN_CHARG(TOP_NUMBER,CHARGE,PART_NUMB) >*/
/* Subroutine */ int gamn_charg__(integer *top_number__, integer *charge, 
	integer *part_numb__)
{
    /* Initialized data */

    static shortint charg_dat__[168]	/* was [4][42] */ = { 1,0,0,2,1,1,1,0,
	    0,1,0,1,1,0,0,3,1,1,1,1,0,1,0,2,0,2,1,0,1,0,0,4,1,1,1,2,1,2,2,0,0,
	    1,0,3,0,2,1,1,1,0,0,5,1,1,1,3,1,2,2,1,0,1,0,4,0,2,1,2,0,3,2,0,1,0,
	    0,6,1,1,1,4,1,2,2,2,1,3,3,0,0,1,0,5,0,2,1,3,0,3,2,1,1,0,0,7,1,1,1,
	    5,1,2,2,3,1,3,3,1,0,1,0,6,0,2,1,4,0,3,2,2,0,4,3,0,1,0,0,8,1,1,1,6,
	    1,2,2,4,1,3,3,2,1,4,4,0,0,1,0,7,0,2,1,5,0,3,2,3,0,4,3,1 };
    static real su2[42] = { .1f,.7f,.2f,.01f,.75f,.04f,.2f,.03333f,.3667f,
	    .2667f,.06667f,.2667f,.004762f,.2238f,.4381f,.02857f,.2095f,
	    .09524f,.00476f,.1155f,.4214f,.125f,.0119f,.1429f,.1786f,6.1e-4f,
	    .05632f,.3342f,.2755f,.004884f,.08425f,.2015f,.04273f,6.105e-4f,
	    .02624f,.226f,.3577f,.05614f,.003831f,.09195f,.3563f,.2146f };

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k, iposition, nf, ni;
    extern integer idihotomia_(real *, real *, integer *);
    static real sum, vran[1], temp1;
    static integer numpi[3];
    static real probab[9];
    static integer iabsol;
    extern /* Subroutine */ int ranmar_(real *, integer *);
    static integer numnuc[2];

/* + */

/* FUNCTIONAL DESCRIPTION: */

/*    This subroutine is a generator of charges of outgoing particles. */
/*    It contains SU(2) coefficients table for N+mpi */

/* CREATION DATE: */

/*     14.03.94       by Igor Pshenichnov */

/* FORMAL PARAMETERS: */


/*    TOP_NUMBER  - topology (number of particles in the final */
/*                  state (m=3,...,9)). */

/*    CHARGE      - charge of the taget nucleon */

/*    PART_NUMB   - outgoing particle numbers, nucleon comes first */


/* COMMON BLOCKS: */

/*    /LOGEG/EG,CMSF   - governs the mode of operation: */

/*   EG: */
/*   PARTICLE NUMBERS ARE AS IN INC (.FALSE.) OR AS IN GEANT (.TRUE.) */

/*   CMSF: */
/*   All the momenta of particles are in CMS (.TRUE.) */
/*   or in LAB system (.FALSE.) */

/* - */
/*<         INTEGER*4 TOP_NUMBER, PART_NUMB(10), NI, NF, NUMPI(3) >*/
/*<         INTEGER*4 NUMNUC(2),CHARGE,IRENUM(7) >*/
/*<         INTEGER*2 CHARG_DAT(4,42) >*/
/*<         REAL*4 SU2(42), PROBAB(9) >*/
/*<         LOGICAL EG,CMSF >*/
/*<         COMMON/LOGEG/EG,CMSF ! Insure compartibility with the INC model. >*/
/*<         real vran(1) >*/
/*< 	DATA  >*/
    /* Parameter adjustments */
    --part_numb__;

    /* Function Body */
/* --------- notation is : N(p=1, n=0) kpi+ lpi- mpi0 -------- */
/* ________________________ 3 PARTICLES ______________________ */
/* ________________________ 4 PARTICLES ______________________ */
/* ________________________ 5 PARTICLES ______________________ */
/* ________________________ 6 PARTICLES ______________________ */
/* ________________________ 7 PARTICLES ______________________ */
/* ________________________ 8 PARTICLES ______________________ */
/* ________________________ 9 PARTICLES ______________________ */
/*< 	DATA  >*/
/* ________________________ 3 PARTICLES ______________________ */
/* ________________________ 4 PARTICLES ______________________ */
/* ________________________ 5 PARTICLES ______________________ */
/* ________________________ 6 PARTICLES ______________________ */
/* ________________________ 7 PARTICLES ______________________ */
/* ________________________ 8 PARTICLES ______________________ */
/* ________________________ 9 PARTICLES ______________________ */
/*<         IF (EG) THEN >*/
    if (logeg_1.eg) {
/*<           NUMNUC(1)=14  ! proton number >*/
	numnuc[0] = 14;
/*<           NUMNUC(2)=13  ! neutron number >*/
	numnuc[1] = 13;
/*<           NUMPI(1)=8 >*/
	numpi[0] = 8;
/*<           NUMPI(2)=9 >*/
	numpi[1] = 9;
/*<           NUMPI(3)=7 >*/
	numpi[2] = 7;
/*<           SU2(4)=0.090902 >*/
	su2[3] = .090902f;
/*<           SU2(5)=0.0  ! since it is simulating via different branch >*/
	su2[4] = 0.f;
/*<           SU2(6)=0.363552 >*/
	su2[5] = .363552f;
/*<           SU2(7)=0.545464 >*/
	su2[6] = .545464f;
/*<         ELSE   ! "as was" >*/
    } else {
/*<           NUMNUC(1)=37  ! proton number >*/
	numnuc[0] = 37;
/*<           NUMNUC(2)=38  ! neutron number >*/
	numnuc[1] = 38;
/*<           NUMPI(1)=1 >*/
	numpi[0] = 1;
/*<           NUMPI(2)=2 >*/
	numpi[1] = 2;
/*<           NUMPI(3)=7 >*/
	numpi[2] = 7;
/*<           SU2(4)=0.01 >*/
	su2[3] = .01f;
/*<           SU2(5)=0.75 >*/
	su2[4] = .75f;
/*<           SU2(6)=0.04 >*/
	su2[5] = .04f;
/*<           SU2(7)=0.2      >*/
	su2[6] = .2f;
/*<         ENDIF >*/
    }
/*<         DO I=1,10 >*/
    for (i__ = 1; i__ <= 10; ++i__) {
/*<           PART_NUMB(I)=0.0 >*/
	part_numb__[i__] = 0.f;
/*<         ENDDO >*/
    }
/* ____ DETERMINE POSITIONS OF CERTAIN MULTIPLICITY IN THE TABLE ___ */
/*<         NF=TOP_NUMBER*(TOP_NUMBER+1)/2+1 >*/
    nf = *top_number__ * (*top_number__ + 1) / 2 + 1;
/*<         NI=NF-TOP_NUMBER+1 >*/
    ni = nf - *top_number__ + 1;
/*<         SUM=0. >*/
    sum = 0.f;
/*<             DO I=1,TOP_NUMBER >*/
    i__1 = *top_number__;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<                   SUM=SUM+SU2(NI+I-5) >*/
	sum += su2[ni + i__ - 6];
/*<                   PROBAB(I)=SUM >*/
	probab[i__ - 1] = sum;
/*<             ENDDO >*/
    }
/*<         call ranmar(vran,1) >*/
    ranmar_(vran, &c__1);
/*<         TEMP1=vran(1) >*/
    temp1 = vran[0];
/*<         IABSOL=IDIHOTOMIA(TEMP1,PROBAB,TOP_NUMBER)+NI-1 >*/
    iabsol = idihotomia_(&temp1, probab, top_number__) + ni - 1;
/*<         IF(CHARG_DAT(1,IABSOL-4).EQ.1) THEN >*/
    if (charg_dat__[(iabsol - 4 << 2) - 4] == 1) {
/*<             PART_NUMB(1)=NUMNUC(1) >*/
	part_numb__[1] = numnuc[0];
/*<         ELSEIF(CHARG_DAT(1,IABSOL-4).EQ.0) THEN >*/
    } else if (charg_dat__[(iabsol - 4 << 2) - 4] == 0) {
/*<             PART_NUMB(1)=NUMNUC(2) >*/
	part_numb__[1] = numnuc[1];
/*<         ENDIF >*/
    }
/*<         IPOSITION=1 >*/
    iposition = 1;
/*<         DO K=1,3 >*/
    for (k = 1; k <= 3; ++k) {
/*<         IF (CHARG_DAT(K+1,IABSOL-4).GT.0) THEN >*/
	if (charg_dat__[k + 1 + (iabsol - 4 << 2) - 5] > 0) {
/*<            DO I=1,CHARG_DAT(K+1,IABSOL-4) >*/
	    i__1 = charg_dat__[k + 1 + (iabsol - 4 << 2) - 5];
	    for (i__ = 1; i__ <= i__1; ++i__) {
/*<                IPOSITION=IPOSITION+1 >*/
		++iposition;
/*<                PART_NUMB(IPOSITION)=NUMPI(K) >*/
		part_numb__[iposition] = numpi[k - 1];
/*<            ENDDO >*/
	    }
/*<         ENDIF >*/
	}
/*<         ENDDO >*/
    }
/*<       IF(CHARGE.EQ.0) THEN   ! CASE OF GAMMA + NEUTRON >*/
    if (*charge == 0) {
/*<       I=1 >*/
	i__ = 1;
/*<       DO WHILE (PART_NUMB(I).GT.0) >*/
	while(part_numb__[i__] > 0) {
/*<         IF(PART_NUMB(I).EQ.NUMNUC(1)) THEN >*/
	    if (part_numb__[i__] == numnuc[0]) {
/*<               PART_NUMB(I)=NUMNUC(2) >*/
		part_numb__[i__] = numnuc[1];
/*<         ELSEIF(PART_NUMB(I).EQ.NUMNUC(2)) THEN >*/
	    } else if (part_numb__[i__] == numnuc[1]) {
/*<               PART_NUMB(I)=NUMNUC(1) >*/
		part_numb__[i__] = numnuc[0];
/*<         ELSEIF(PART_NUMB(I).EQ.NUMPI(2)) THEN >*/
	    } else if (part_numb__[i__] == numpi[1]) {
/*<               PART_NUMB(I)=NUMPI(1) >*/
		part_numb__[i__] = numpi[0];
/*<         ELSEIF(PART_NUMB(I).EQ.NUMPI(1)) THEN >*/
	    } else if (part_numb__[i__] == numpi[0]) {
/*<               PART_NUMB(I)=NUMPI(2) >*/
		part_numb__[i__] = numpi[1];
/*<         ENDIF >*/
	    }
/*<          I=I+1 >*/
	    ++i__;
/*<       ENDDO >*/
	}
/*<       ENDIF >*/
    }
/*<         RETURN >*/
    return 0;
/*< 	END >*/
} /* gamn_charg__ */

