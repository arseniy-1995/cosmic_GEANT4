/* prim_kine.F -- translated by f2c (version 20230428).
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

struct {
    integer itarg, iidh;
    char target[12], targfermi[12];
} target_;

#define target_1 target_

struct {
    char partic[61440]	/* was [40][4][4][8] */;
    integer jdec[40], nch[40], nplev[120]	/* was [40][3] */;
    real br[120]	/* was [40][3] */;
} channels_;

#define channels_1 channels_

struct {
    real pn[15]	/* was [5][3] */;
} mom_nucleon__;

#define mom_nucleon__1 mom_nucleon__

struct {
    integer j_channel__;
} reacch_;

#define reacch_1 reacch_

struct {
    integer i_1__[4], npart_1__;
    real pprim[90]	/* was [5][18] */, p_1__[44]	/* was [4][11] */;
} kine_1__;

#define kine_1__1 kine_1__

struct {
    real p_fermi__, rm[2];
} fermi_motion__;

#define fermi_motion__1 fermi_motion__

struct {
    real beta_cm__[4];
} betacm_;

#define betacm_1 betacm_

/* Table of constant values */

static doublereal c_b2 = 2.;

/* *********************************************************************** */
/*<       SUBROUTINE prim_kine >*/
/* Subroutine */ int prim_kine__(void)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer j, jj, jch;
    extern integer code_(char *, ftnlen);
    extern doublereal rmass_(char *, ftnlen);

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* FUNCTIONAL DESCRIPTION: */

/* Extracts kinematics of outgoing particles after first level */
/* and before decays. Channels 5->10 (delta), 11->14 and 35 (rho) */
/* are involved. */

/* MODIFICATION DATE: */

/*     05.05.99       by marco mirazita */


/* FORMAL PARAMETERS: */

/*    Npart_1    = NUMBER OF OUTGOING PARTICLES (2-18) */
/*    RMASS      = MASS OF I-TH spectator PARTICLE (REAL) */
/*    PN(1,I)    = PX OF I-TH spectator PARTICLE */
/*    PN(2,I)    = PY OF I-TH spectator PARTICLE */
/*    PN(3,I)    = PZ OF I-TH spectator PARTICLE */
/*    PN(4,I)    = ENERGY OF I-TH spectator PARTICLE */
/*    PN(5,I)    = |P| OF I-TH spectator PARTICLE */

/*    PPRIM(1,I)   = PX OF I-TH primary partICLE */
/*    PPRIM(2,I)   = PY OF I-TH primary partICLE */
/*    PPRIM(3,I)   = PZ OF I-TH primary partICLE */
/*    PPRIM(4,I)   = ENERGY OF I-TH primary partICLE */

/*    P_1(1,I)  =  PX OF I-TH primary partICLE */
/*    P_1(2,I)  =  PY OF I-TH primary partICLE */
/*    P_1(3,I)  =  PZ OF I-TH primary partICLE */
/*    P_1(4,I)  =  MASS OF I-TH primary partICLE */

/* COMMON BLOCKS: */
/*     /TARGET/ */
/*     /CHANNELS/ */
/*     /MOM_NUCLEON/ */
/* ----------------------------------------------------------------------- */
/*<       IMPLICIT NONE >*/
/*<       INTEGER ITARG,IIDH >*/
/*<       CHARACTER*12 TARGET,TARGFERMI >*/
/*<       COMMON /TARGET/ ITARG,IIDH,TARGET,TARGFERMI >*/
/*<       CHARACTER*12 PARTIC(40,0:3,0:3,8) >*/
/*<       INTEGER JDEC(40),NCH(40),NPLEV(40,3) >*/
/*<       REAL BR(40,3) >*/
/*<       COMMON /CHANNELS/ PARTIC,JDEC,NCH,NPLEV,BR >*/
/*<       COMMON /MOM_NUCLEON/PN >*/
/*<       real pn(5,3) >*/
/*<       real rmass,xm >*/
/*<       integer j,jj,jch,code >*/
/* *** ntuple variables */
/*     #include "nt_kine.inc" */
/*<       INCLUDE 'nt_kine.inc' >*/
/*<       jch=j_channel >*/
/* *** canale di reazione */
/*< 	integer j_channel >*/
/* *** common per le variabili della cinematica di primo livello */
/*<         common/kine_1/I_1,npart_1,pprim,P_1 >*/
/*<         real P_1(4,11),pprim(5,18) >*/
/*<         integer I_1(4),npart_1 >*/
/* *** Impulso di fermi */
/*< 	common/fermi_motion/p_fermi,rm >*/
/*< 	real p_fermi,rm(2) >*/
/* *** Beta(CM) of (gamma+target) system */
/*< 	common/betacm/beta_cm >*/
/*< 	real beta_cm(4) >*/
    jch = reacch_1.j_channel__;
/* *** initialization */
/*<       npart_1=0 >*/
    kine_1__1.npart_1__ = 0;
/*<       do j=1,4 >*/
    for (j = 1; j <= 4; ++j) {
/*<          do jj=1,4 >*/
	for (jj = 1; jj <= 4; ++jj) {
/*<             P_1(jj,j)=0. >*/
	    kine_1__1.p_1__[jj + (j << 2) - 5] = 0.f;
/*<          enddo >*/
	}
/*<          I_1(j)=0 >*/
	kine_1__1.i_1__[j - 1] = 0;
/*<       enddo >*/
    }
/* ***    Canali 5,6,7,8,9,10 */
/*<       IF (JCH.GE.5.AND.JCH.LE.10) THEN >*/
    if (jch >= 5 && jch <= 10) {
/*<          p_1(1,1) = pprim(1,1) >*/
	kine_1__1.p_1__[0] = kine_1__1.pprim[0];
/*<          p_1(2,1) = pprim(2,1) >*/
	kine_1__1.p_1__[1] = kine_1__1.pprim[1];
/*<          p_1(3,1) = pprim(3,1) >*/
	kine_1__1.p_1__[2] = kine_1__1.pprim[2];
/*         p_1(4,1) = pprim(4,1) */
/*<        >*/
	d__1 = (doublereal) kine_1__1.pprim[3];
	d__2 = (doublereal) kine_1__1.pprim[0];
	d__3 = (doublereal) kine_1__1.pprim[1];
	d__4 = (doublereal) kine_1__1.pprim[2];
	kine_1__1.p_1__[3] = sqrt(pow_dd(&d__1, &c_b2) - (pow_dd(&d__2, &c_b2)
		 + pow_dd(&d__3, &c_b2) + pow_dd(&d__4, &c_b2)));
/*         write (*,'(a)') PARTIC(JCH,1,0,1) */
/*<          I_1(1) = CODE(PARTIC(JCH,1,0,1)) >*/
	kine_1__1.i_1__[0] = code_(channels_1.partic + (jch + 39) * 12, (
		ftnlen)12);
/*<          p_1(1,2) = pprim(1,2) >*/
	kine_1__1.p_1__[4] = kine_1__1.pprim[5];
/*<          p_1(2,2) = pprim(2,2) >*/
	kine_1__1.p_1__[5] = kine_1__1.pprim[6];
/*<          p_1(3,2) = pprim(3,2) >*/
	kine_1__1.p_1__[6] = kine_1__1.pprim[7];
/*         p_1(4,2) = pprim(4,2) */
/*<        >*/
	d__1 = (doublereal) kine_1__1.pprim[8];
	d__2 = (doublereal) kine_1__1.pprim[5];
	d__3 = (doublereal) kine_1__1.pprim[6];
	d__4 = (doublereal) kine_1__1.pprim[7];
	kine_1__1.p_1__[7] = sqrt(pow_dd(&d__1, &c_b2) - (pow_dd(&d__2, &c_b2)
		 + pow_dd(&d__3, &c_b2) + pow_dd(&d__4, &c_b2)));
/*        write (*,'(a)') PARTIC(JCH,1,0,2) */
/*<          I_1(2) = CODE(PARTIC(JCH,1,0,2)) >*/
	kine_1__1.i_1__[1] = code_(channels_1.partic + (jch + 679) * 12, (
		ftnlen)12);
/*<          NPART_1 = 2 >*/
	kine_1__1.npart_1__ = 2;
/*<          IF (TARGET.EQ.'deuteron'.OR.TARGET.EQ.'He3') THEN >*/
	if (s_cmp(target_1.target, "deuteron", (ftnlen)12, (ftnlen)8) == 0 || 
		s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<             p_1(1,3) = PN(1,2) >*/
	    kine_1__1.p_1__[8] = mom_nucleon__1.pn[5];
/*<             p_1(2,3) = PN(2,2) >*/
	    kine_1__1.p_1__[9] = mom_nucleon__1.pn[6];
/*<             p_1(3,3) = PN(3,2) >*/
	    kine_1__1.p_1__[10] = mom_nucleon__1.pn[7];
/*<             p_1(4,3) = RMASS(PARTIC(JCH,1,0,3)) >*/
	    kine_1__1.p_1__[11] = rmass_(channels_1.partic + (jch + 1319) * 
		    12, (ftnlen)12);
/*<             I_1(3) = CODE(PARTIC(JCH,1,0,3)) >*/
	    kine_1__1.i_1__[2] = code_(channels_1.partic + (jch + 1319) * 12, 
		    (ftnlen)12);
/*<             NPART_1 = NPART_1+1 >*/
	    ++kine_1__1.npart_1__;
/*<             IF (TARGET.EQ.'He3') THEN >*/
	    if (s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<                p_1(1,4) = PN(1,3) >*/
		kine_1__1.p_1__[12] = mom_nucleon__1.pn[10];
/*<                p_1(2,4) = PN(2,3) >*/
		kine_1__1.p_1__[13] = mom_nucleon__1.pn[11];
/*<                p_1(3,4) = PN(3,3) >*/
		kine_1__1.p_1__[14] = mom_nucleon__1.pn[12];
/*<                p_1(4,4) = RMASS(PARTIC(JCH,1,0,4)) >*/
		kine_1__1.p_1__[15] = rmass_(channels_1.partic + (jch + 1959) 
			* 12, (ftnlen)12);
/*<                I_1(4) = CODE(PARTIC(JCH,1,0,4)) >*/
		kine_1__1.i_1__[3] = code_(channels_1.partic + (jch + 1959) * 
			12, (ftnlen)12);
/*<                NPART_1 = NPART_1+1 >*/
		++kine_1__1.npart_1__;
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<          GOTO 99 >*/
	goto L99;
/*<       END IF >*/
    }
/* ***    Canali 11,12,13,14 */
/*<       IF (JCH.GE.11.AND.JCH.LE.14) THEN >*/
    if (jch >= 11 && jch <= 14) {
/*<          p_1(1,1) = pprim(1,1) >*/
	kine_1__1.p_1__[0] = kine_1__1.pprim[0];
/*<          p_1(2,1) = pprim(2,1) >*/
	kine_1__1.p_1__[1] = kine_1__1.pprim[1];
/*<          p_1(3,1) = pprim(3,1) >*/
	kine_1__1.p_1__[2] = kine_1__1.pprim[2];
/*         p_1(4,1) = pprim(4,1) */
/*<        >*/
	d__1 = (doublereal) kine_1__1.pprim[3];
	d__2 = (doublereal) kine_1__1.pprim[0];
	d__3 = (doublereal) kine_1__1.pprim[1];
	d__4 = (doublereal) kine_1__1.pprim[2];
	kine_1__1.p_1__[3] = sqrt(pow_dd(&d__1, &c_b2) - (pow_dd(&d__2, &c_b2)
		 + pow_dd(&d__3, &c_b2) + pow_dd(&d__4, &c_b2)));
/*        write (*,'(a)') PARTIC(JCH,1,0,1) */
/*<          I_1(1) = CODE(PARTIC(JCH,1,0,1)) >*/
	kine_1__1.i_1__[0] = code_(channels_1.partic + (jch + 39) * 12, (
		ftnlen)12);
/*<          p_1(1,2) = pprim(1,2) >*/
	kine_1__1.p_1__[4] = kine_1__1.pprim[5];
/*<          p_1(2,2) = pprim(2,2) >*/
	kine_1__1.p_1__[5] = kine_1__1.pprim[6];
/*<          p_1(3,2) = pprim(3,2) >*/
	kine_1__1.p_1__[6] = kine_1__1.pprim[7];
/*         p_1(4,2) = pprim(4,2) */
/*<        >*/
	d__1 = (doublereal) kine_1__1.pprim[8];
	d__2 = (doublereal) kine_1__1.pprim[5];
	d__3 = (doublereal) kine_1__1.pprim[6];
	d__4 = (doublereal) kine_1__1.pprim[7];
	kine_1__1.p_1__[7] = sqrt(pow_dd(&d__1, &c_b2) - (pow_dd(&d__2, &c_b2)
		 + pow_dd(&d__3, &c_b2) + pow_dd(&d__4, &c_b2)));
/*        write (*,'(a)') PARTIC(JCH,1,0,2) */
/*<          I_1(2) = CODE(PARTIC(JCH,1,0,2)) >*/
	kine_1__1.i_1__[1] = code_(channels_1.partic + (jch + 679) * 12, (
		ftnlen)12);
/*<          NPART_1 = 2 >*/
	kine_1__1.npart_1__ = 2;
/*<          IF (TARGET.EQ.'deuteron'.OR.TARGET.EQ.'He3') THEN >*/
	if (s_cmp(target_1.target, "deuteron", (ftnlen)12, (ftnlen)8) == 0 || 
		s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<             p_1(1,3) = PN(1,2) >*/
	    kine_1__1.p_1__[8] = mom_nucleon__1.pn[5];
/*<             p_1(2,3) = PN(2,2) >*/
	    kine_1__1.p_1__[9] = mom_nucleon__1.pn[6];
/*<             p_1(3,3) = PN(3,2) >*/
	    kine_1__1.p_1__[10] = mom_nucleon__1.pn[7];
/*<             p_1(4,3) = RMASS(PARTIC(JCH,1,0,3)) >*/
	    kine_1__1.p_1__[11] = rmass_(channels_1.partic + (jch + 1319) * 
		    12, (ftnlen)12);
/*<             I_1(3) = CODE(PARTIC(JCH,1,0,3)) >*/
	    kine_1__1.i_1__[2] = code_(channels_1.partic + (jch + 1319) * 12, 
		    (ftnlen)12);
/*<             NPART_1 = NPART_1+1 >*/
	    ++kine_1__1.npart_1__;
/*<             IF (TARGET.EQ.'He3') THEN >*/
	    if (s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<                p_1(1,4) = PN(1,3) >*/
		kine_1__1.p_1__[12] = mom_nucleon__1.pn[10];
/*<                p_1(2,4) = PN(2,3) >*/
		kine_1__1.p_1__[13] = mom_nucleon__1.pn[11];
/*<                p_1(3,4) = PN(3,3) >*/
		kine_1__1.p_1__[14] = mom_nucleon__1.pn[12];
/*<                p_1(4,4) = RMASS(PARTIC(JCH,1,0,4)) >*/
		kine_1__1.p_1__[15] = rmass_(channels_1.partic + (jch + 1959) 
			* 12, (ftnlen)12);
/*<                I_1(4) = CODE(PARTIC(JCH,1,0,4)) >*/
		kine_1__1.i_1__[3] = code_(channels_1.partic + (jch + 1959) * 
			12, (ftnlen)12);
/*<                NPART_1 = NPART_1+1 >*/
		++kine_1__1.npart_1__;
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<          GOTO 99 >*/
	goto L99;
/*<       END IF >*/
    }
/* ***    Canale 35 */
/*<       IF (JCH.EQ.35) THEN >*/
    if (jch == 35) {
/*<          p_1(1,1) = pprim(1,1) >*/
	kine_1__1.p_1__[0] = kine_1__1.pprim[0];
/*<          p_1(2,1) = pprim(2,1) >*/
	kine_1__1.p_1__[1] = kine_1__1.pprim[1];
/*<          p_1(3,1) = pprim(3,1) >*/
	kine_1__1.p_1__[2] = kine_1__1.pprim[2];
/*         p_1(4,1) = pprim(4,1) */
/*<        >*/
	d__1 = (doublereal) kine_1__1.pprim[3];
	d__2 = (doublereal) kine_1__1.pprim[0];
	d__3 = (doublereal) kine_1__1.pprim[1];
	d__4 = (doublereal) kine_1__1.pprim[2];
	kine_1__1.p_1__[3] = sqrt(pow_dd(&d__1, &c_b2) - (pow_dd(&d__2, &c_b2)
		 + pow_dd(&d__3, &c_b2) + pow_dd(&d__4, &c_b2)));
/*     write (*,'(a)') PARTIC(JCH,1,0,1) */
/*<          I_1(1) = CODE(PARTIC(JCH,1,0,1)) >*/
	kine_1__1.i_1__[0] = code_(channels_1.partic + (jch + 39) * 12, (
		ftnlen)12);
/*<          p_1(1,2) = pprim(1,2) >*/
	kine_1__1.p_1__[4] = kine_1__1.pprim[5];
/*<          p_1(2,2) = pprim(2,2) >*/
	kine_1__1.p_1__[5] = kine_1__1.pprim[6];
/*<          p_1(3,2) = pprim(3,2) >*/
	kine_1__1.p_1__[6] = kine_1__1.pprim[7];
/*         p_1(4,2) = pprim(4,2) */
/*<        >*/
	d__1 = (doublereal) kine_1__1.pprim[8];
	d__2 = (doublereal) kine_1__1.pprim[5];
	d__3 = (doublereal) kine_1__1.pprim[6];
	d__4 = (doublereal) kine_1__1.pprim[7];
	kine_1__1.p_1__[7] = sqrt(pow_dd(&d__1, &c_b2) - (pow_dd(&d__2, &c_b2)
		 + pow_dd(&d__3, &c_b2) + pow_dd(&d__4, &c_b2)));
/*     write (*,'(a)') PARTIC(JCH,1,0,2) */
/*<          I_1(2) = CODE(PARTIC(JCH,1,0,2)) >*/
	kine_1__1.i_1__[1] = code_(channels_1.partic + (jch + 679) * 12, (
		ftnlen)12);
/*<          NPART_1 = 2 >*/
	kine_1__1.npart_1__ = 2;
/*<          IF (TARGET.EQ.'He3') THEN >*/
	if (s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<             p_1(1,4) = PN(1,3) >*/
	    kine_1__1.p_1__[12] = mom_nucleon__1.pn[10];
/*<             p_1(2,4) = PN(2,3) >*/
	    kine_1__1.p_1__[13] = mom_nucleon__1.pn[11];
/*<             p_1(3,4) = PN(3,3) >*/
	    kine_1__1.p_1__[14] = mom_nucleon__1.pn[12];
/*<             p_1(4,4) = RMASS(PARTIC(JCH,1,0,4)) >*/
	    kine_1__1.p_1__[15] = rmass_(channels_1.partic + (jch + 1959) * 
		    12, (ftnlen)12);
/*<             I_1(4) = CODE(PARTIC(JCH,1,0,4)) >*/
	    kine_1__1.i_1__[3] = code_(channels_1.partic + (jch + 1959) * 12, 
		    (ftnlen)12);
/*<             NPART_1 = NPART_1+1 >*/
	    ++kine_1__1.npart_1__;
/*<          END IF >*/
	}
/*<          GOTO 99 >*/
	goto L99;
/*<       END IF >*/
    }
/*<  99   continue >*/
L99:
/*<       RETURN >*/
    return 0;
/*<       END                       !*** END prim_kine >*/
} /* prim_kine__ */

