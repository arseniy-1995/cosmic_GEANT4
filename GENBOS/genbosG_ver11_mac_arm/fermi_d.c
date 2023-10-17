/* fermi_d.F -- translated by f2c (version 20230428).
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
    integer nprop[20], npron[20], nprod[4], nprohe3[4];
} flag_proc__;

#define flag_proc__1 flag_proc__

struct {
    integer nprop_tot__, npron_tot__, nprod_tot__, nprohe3_tot__;
} proc_tot__;

#define proc_tot__1 proc_tot__

struct {
    integer nchain_p__, nchain_n__, nchain_d__, nchain_he3__;
} proc_sel__;

#define proc_sel__1 proc_sel__

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

struct {
    integer ielet, ibeam, icible;
    real etag, fwhm;
    integer nchain, ichain[40], ipol, ilam;
} flag_comp__;

#define flag_comp__1 flag_comp__

struct {
    real egam;
    integer i__[11], npart;
    real p[44]	/* was [4][11] */;
} montecarl_;

#define montecarl_1 montecarl_

struct {
    integer itarg, iidh;
    char target[12], targfermi[12];
} target_;

#define target_1 target_

struct {
    real pn_n_d__[15]	/* was [5][3] */, pn_d__[15]	/* was [5][3] */, 
	    beta_n_d__[4], beta_d__[4];
} h2_;

#define h2_1 h2_

struct {
    real rmb1, rmb2, rmme, deltamass, dibmass, deumass, rhomass;
} masses_;

#define masses_1 masses_

/* Table of constant values */

static integer c__5 = 5;
static doublereal c_b7 = 2.;

/* DECK  ID>, FERMI_D. */
/* CMZ :          13/04/94  18.02.52  by  Thomas Russew */
/* -- Author : */
/* *********************************************************************** */
/*<       SUBROUTINE FERMI_D(W_N,W_DEU) >*/
/* Subroutine */ int fermi_d__(real *w_n__, real *w_deu__)
{
    /* Initialized data */

    static real distr_p__[101] = { 0.f,.0067953602f,.046238571f,.1236217f,
	    .2234159f,.3280089f,.4262296f,.5131546f,.5876556f,.6504297f,
	    .70287f,.746512f,.7827951f,.8129796f,.8381307f,.8591338f,
	    .8767187f,.891484f,.9039201f,.9144295f,.9233416f,.9309281f,
	    .937412f,.9429782f,.9477784f,.9519395f,.9555655f,.958744f,
	    .9615465f,.9640333f,.9662549f,.9682527f,.9700612f,.9717091f,
	    .9732212f,.9746166f,.9759121f,.9771215f,.9782565f,.9793263f,
	    .9803385f,.9812997f,.9822153f,.9830891f,.9839256f,.9847271f,
	    .9854963f,.9862349f,.9869449f,.9876275f,.9882837f,.9889148f,
	    .9895214f,.9901043f,.9906641f,.9912014f,.9917168f,.9922106f,
	    .9926836f,.9931361f,.9935685f,.9939811f,.9943746f,.9947497f,
	    .9951068f,.9954458f,.9957677f,.9960729f,.9963619f,.996635f,
	    .996893f,.9971365f,.9973655f,.9975811f,.9977834f,.9979731f,
	    .9981509f,.9983172f,.9984723f,.9986169f,.9987513f,.9988762f,
	    .9989921f,.9990993f,.9991982f,.9992898f,.9993736f,.999451f,
	    .9995219f,.9995865f,.9996458f,.9996991f,.9997482f,.9997926f,
	    .9998326f,.9998688f,.999901f,.99993f,.9999564f,.9999794f,1.f };
    static real distr_b__[101] = { 0.f,.0070403279f,.047886062f,.1279474f,
	    .2310553f,.33893f,.4400139f,.5292494f,.6055161f,.6695801f,
	    .7229201f,.7671533f,.8037899f,.8341466f,.8593346f,.8802744f,
	    .8977222f,.9122964f,.924503f,.934754f,.9433876f,.9506804f,
	    .9568596f,.9621121f,.9665918f,.9704258f,.9737193f,.9765596f,
	    .9790186f,.9811566f,.9830235f,.9846612f,.9861036f,.9873803f,
	    .9885156f,.9895294f,.990438f,.9912565f,.9919961f,.9926671f,
	    .9932781f,.9938356f,.9943457f,.994814f,.995244f,.9956399f,
	    .9960047f,.9963413f,.9966521f,.9969389f,.9972037f,.9974484f,
	    .9976741f,.9978822f,.9980739f,.9982503f,.9984124f,.9985616f,
	    .998698f,.9988228f,.9989371f,.9990408f,.9991356f,.9992215f,
	    .9992995f,.9993696f,.9994332f,.99949f,.999541f,.9995866f,
	    .9996274f,.9996639f,.9996958f,.9997242f,.9997495f,.9997714f,
	    .999791f,.9998078f,.9998228f,.9998358f,.9998475f,.9998578f,
	    .9998671f,.999875f,.9998825f,.9998896f,.999896f,.9999021f,
	    .9999081f,.9999138f,.99992f,.9999259f,.9999323f,.9999395f,
	    .9999465f,.9999539f,.9999623f,.9999706f,.9999802f,.9999896f,1.f };
    static integer idp = 1;

    /* System generated locals */
    real r__1, r__2, r__3;
    doublereal d__1;

    /* Builtin functions */
    double acos(doublereal), sin(doublereal), cos(doublereal), sqrt(
	    doublereal), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static real totmass1, totmass2, e_sp_fin__;
    extern doublereal mass_n_d__(real *);
    static real a, b;
    static integer j;
    static real q, s, p1;
    static integer jj;
    static real pp[101], ecm;
    static integer icn, icp;
    static real rnd[5];
    static integer ideu;
    static real etot;
    extern doublereal rmass_(char *, ftnlen);
    static real rm_fin__;
    extern /* Subroutine */ int ranmar_(real *, integer *);
    static real delta_e__, totmass;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* FUNCTIONAL DESCRIPTION: */

/*    This subroutine computes invariant mass taking into account Fermi */
/*    motion. It defines also tetramomenta and beta of nucleons in CM */

/*    The momentum distribution is calculated by Vistor, wich uses */
/*    Paris or Bonn potential */

/*     MODIFICHE */

/*      ******* by marco mirazita */
/*     29.04.99 */
/* --> Aggiunto il file di include PROCESS.INC con le variabili per il */
/*     conteggio dei canali di reazione selezionati. */

/*     05.05.99 */
/* --> Eliminato il conteggio dei canali di reazione selezionati, che ora e' */
/*     fatto una volta per tutte nella routine CARDS (che sta in IN.F). */
/* --> Eliminata la parte in cui, nel caso siano selezionati sia canali su p */
/*     che su n, veniva scelto il nucleone bersaglio perche' era sbagliata. */
/*     Adesso la massa del nucleone bersaglio e' sempre "RMASS('proton   ')" */
/*     e quella dello spettatore "RMASS('neutron   ')", tanto i loro valori */
/*     sono definiti uguali. */
/* --> Spostata l'estrazione dei numeri random subito dopo gli azzeramenti */
/*     iniziali. */

/* COMMON BLOCKS: */

/*     /MONTECARL/ */
/*     /FLAG_COMP/ */
/*     /TARGET/ */
/*     /H2/ */
/*     /FLAG_PROC/ */
/*     /MASSES/ */

/* ----------------------------------------------------------------------- */
/*<       IMPLICIT NONE >*/
/*   #include "process.inc" */
/*    #include "nt_kine.inc" */
/*     include "process.inc" */
/*<       INCLUDE 'process.inc' >*/
/*<       INCLUDE 'nt_kine.inc' >*/
/*<       COMMON /FLAG_PROC/ NPROP,NPRON,NPROD,NPROHE3 >*/
/*< 	integer nprop_tot,npron_tot,nprod_tot,nprohe3_tot >*/
/*< 	common/proc_tot/nprop_tot,npron_tot,nprod_tot,nprohe3_tot >*/
/*< 	integer nchain_p,nchain_n,nchain_d,nchain_he3 >*/
/*< 	common/proc_sel/nchain_p,nchain_n,nchain_d,nchain_he3 >*/
/*<       INTEGER IBEAM,IELET,NCHAIN,ICHAIN(40),ICIBLE,IPOL,ILAM >*/
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
/*<       REAL ETAG,FWHM >*/
/*<        >*/
/*<       REAL EGAM, P(4,11) >*/
/*<       INTEGER I(11),NPART >*/
/*<       COMMON/MONTECARL/ EGAM,I,NPART,P >*/
/*<       DOUBLE PRECISION PI,TWOPI,PIBY2,DEGRAD,RADDEG,CLIGHT,BIG,EMASS >*/
/*<       DOUBLE PRECISION EMMU,PMASS,AVO >*/

/*<       PARAMETER (PI=3.14159265358979324D0) >*/
/*<       PARAMETER (TWOPI=6.28318530717958648D0) >*/
/*<       PARAMETER (PIBY2=1.57079632679489662D0) >*/
/*<       PARAMETER (DEGRAD=0.0174532925199432958D0) >*/
/*<       PARAMETER (RADDEG=57.2957795130823209D0) >*/
/*<       PARAMETER (CLIGHT=29979245800.D0) >*/
/*<       PARAMETER (BIG=10000000000.D0) >*/
/*<       PARAMETER (EMASS=0.0005109990615D0) >*/
/*<       PARAMETER (EMMU=0.105658387D0) >*/
/*<       PARAMETER (PMASS=0.9382723128D0) >*/
/*<       PARAMETER (AVO=0.60221367D0) >*/

/*<       INTEGER ITARG,IIDH >*/
/*<       CHARACTER*12 TARGET,TARGFERMI >*/
/*<       COMMON /TARGET/ ITARG,IIDH,TARGET,TARGFERMI >*/
/*<       REAL PN_N_D(5,3),PN_D(5,3),BETA_N_D(4),BETA_D(4) >*/
/*<       COMMON /H2/ PN_N_D,PN_D,BETA_N_D,BETA_D >*/
/*<       REAL RMB1,RMB2,RMME,DELTAMASS,DIBMASS,DEUMASS,RHOMASS >*/
/*<       COMMON /MASSES/ RMB1,RMB2,RMME,DELTAMASS,DIBMASS,DEUMASS,RHOMASS >*/
/*<       INTEGER ICP,ICN,J,JJ,K,KK,IDEU,IDP >*/
/*<        >*/
/*      REAL RM(2) */
/*<       real mass_n_d,rm_fin,e_sp_fin,delta_e,totmass,totmass1,totmass2 >*/
/*<       integer kv >*/
/*<        >*/
/*<        >*/
/*<       DATA IDP /1/         !Paris potential   (IDP=0 : Bonn potential) >*/
/* ***   Inizializzazioni */
/*<       DO J=0,99 >*/
    for (j = 0; j <= 99; ++j) {
/*<           PP(J) = 0.01*J >*/
	pp[j] = j * .01f;
/*<       END DO >*/
    }
/*<       ICP = 0 >*/
    icp = 0;
/*<       ICN = 0 >*/
    icn = 0;
/*<       IDEU = 0 >*/
    ideu = 0;
/*<       W_N = 0. >*/
    *w_n__ = 0.f;
/*<       W_DEU = 0. >*/
    *w_deu__ = 0.f;
/*<       ECM = 0. >*/
    ecm = 0.f;
/*<       DO J=1,5 >*/
    for (j = 1; j <= 5; ++j) {
/*<         DO JJ=1,3 >*/
	for (jj = 1; jj <= 3; ++jj) {
/*<           PN_N_D(J,JJ) = 0. >*/
	    h2_1.pn_n_d__[j + jj * 5 - 6] = 0.f;
/*<           PN_D(J,JJ) = 0. >*/
	    h2_1.pn_d__[j + jj * 5 - 6] = 0.f;
/*<         END DO >*/
	}
/*<       END DO >*/
    }
/* ************************************* */
/* ***   Target selection */
/* ***   NCHAIN_P = number of selected channels on proton */
/* ***   NCHAIN_N = number of selected channels on neutron */
/* ************************************* */
/* *** Estrazione di 5 numeri random */
/*<  3    continue >*/
L3:
/*<       call ranmar(RND,5) >*/
    ranmar_(rnd, &c__5);
/*<  1    continue >*/
/* L1: */
/* ***   Computation of W_N_D and BETA_D for gamma-nucleon interactions */
/* ***   one nucleon is interacting, the other one is spectator */
/*<       A = ACOS(1.-2.*RND(2)) >*/
    a = acos(1.f - rnd[1] * 2.f);
/*<       B = 2.*PI*RND(3) >*/
    b = rnd[2] * 6.2831853071795862;
/*<       Q = RND(4) >*/
    q = rnd[3];
/*<       DO J=0,99 >*/
    for (j = 0; j <= 99; ++j) {
/*<          IF (IDP.EQ.1) THEN >*/
	if (idp == 1) {
/*<             IF (Q.GE.DISTR_P(J).AND.Q.LT.DISTR_P(J+1)) GOTO 2 >*/
	    if (q >= distr_p__[j] && q < distr_p__[j + 1]) {
		goto L2;
	    }
/*<          ELSE >*/
	} else {
/*<             IF (Q.GE.DISTR_B(J).AND.Q.LT.DISTR_B(J+1)) GOTO 2 >*/
	    if (q >= distr_b__[j] && q < distr_b__[j + 1]) {
		goto L2;
	    }
/*<          END IF >*/
	}
/*<       END DO >*/
    }
/*<  2    P1 = PP(J)+RND(5)*(PP(J+1)-PP(J)) >*/
L2:
    p1 = pp[j] + rnd[4] * (pp[j + 1] - pp[j]);
/* *** Check on nucleon mass inside deuteron */
/*<       totmass=deumass/2. >*/
    totmass = masses_1.deumass / 2.f;
/*<       if (abs(p1).gt.totmass) go to 3 >*/
    if (dabs(p1) > totmass) {
	goto L3;
    }
/* ************************************************** */
/* *** Calculation of nucleon masses inside the deuterons */
/* *** (effective masses) */
/* ************************************************** */
/*<       rm(1)=mass_n_d(p1) >*/
    fermi_motion__1.rm[0] = mass_n_d__(&p1);
/*<       rm(2)=mass_n_d(p1) >*/
    fermi_motion__1.rm[1] = mass_n_d__(&p1);
/* *** Target nucleon */
/*<       PN_N_D(1,2) = ABS(P1)*SIN(A)*COS(B) >*/
    h2_1.pn_n_d__[5] = dabs(p1) * sin(a) * cos(b);
/*<       PN_N_D(2,2) = ABS(P1)*SIN(A)*SIN(B) >*/
    h2_1.pn_n_d__[6] = dabs(p1) * sin(a) * sin(b);
/*<       PN_N_D(3,2) = ABS(P1)*COS(A) >*/
    h2_1.pn_n_d__[7] = dabs(p1) * cos(a);
/*<       PN_N_D(5,2) = SQRT(PN_N_D(1,2)**2+PN_N_D(2,2)**2+PN_N_D(3,2)**2) >*/
/* Computing 2nd power */
    r__1 = h2_1.pn_n_d__[5];
/* Computing 2nd power */
    r__2 = h2_1.pn_n_d__[6];
/* Computing 2nd power */
    r__3 = h2_1.pn_n_d__[7];
    h2_1.pn_n_d__[9] = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
/*<       PN_N_D(4,2) = SQRT(PN_N_D(5,2)**2+RM(2)**2) >*/
/* Computing 2nd power */
    r__1 = h2_1.pn_n_d__[9];
/* Computing 2nd power */
    r__2 = fermi_motion__1.rm[1];
    h2_1.pn_n_d__[8] = sqrt(r__1 * r__1 + r__2 * r__2);
/* *** Spectator nucleon */
/*<       PN_N_D(1,1) = -PN_N_D(1,2) >*/
    h2_1.pn_n_d__[0] = -h2_1.pn_n_d__[5];
/*<       PN_N_D(2,1) = -PN_N_D(2,2) >*/
    h2_1.pn_n_d__[1] = -h2_1.pn_n_d__[6];
/*<       PN_N_D(3,1) = -PN_N_D(3,2) >*/
    h2_1.pn_n_d__[2] = -h2_1.pn_n_d__[7];
/*<       PN_N_D(5,1) = SQRT(PN_N_D(1,1)**2+PN_N_D(2,1)**2+PN_N_D(3,1)**2) >*/
/* Computing 2nd power */
    r__1 = h2_1.pn_n_d__[0];
/* Computing 2nd power */
    r__2 = h2_1.pn_n_d__[1];
/* Computing 2nd power */
    r__3 = h2_1.pn_n_d__[2];
    h2_1.pn_n_d__[4] = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
/*<       PN_N_D(4,1) = SQRT(PN_N_D(5,1)**2+RM(1)**2) >*/
/* Computing 2nd power */
    r__1 = h2_1.pn_n_d__[4];
/* Computing 2nd power */
    r__2 = fermi_motion__1.rm[0];
    h2_1.pn_n_d__[3] = sqrt(r__1 * r__1 + r__2 * r__2);
/* *********************************************************** */
/* *** Spectator nucleon has reduced mass. I should subtract */
/* *** from total energy for the reaction the mass variation of */
/* *** the spectator */
/* *********************************************************** */
/*<       ETOT = EGAM+PN_N_D(4,1) >*/
    etot = montecarl_1.egam + h2_1.pn_n_d__[3];
/* *** Protone e neutrone hanno masse uguali */
/*<       rm_fin= RMASS('proton      ') >*/
    rm_fin__ = rmass_("proton      ", (ftnlen)12);
/*<       E_sp_fin=sqrt(PN_N_D(5,1)**2+rm_fin**2.) >*/
/* Computing 2nd power */
    r__1 = h2_1.pn_n_d__[4];
    d__1 = (doublereal) rm_fin__;
    e_sp_fin__ = sqrt(r__1 * r__1 + pow_dd(&d__1, &c_b7));
/*<       delta_e=E_sp_fin-PN_N_D(4,1) >*/
    delta_e__ = e_sp_fin__ - h2_1.pn_n_d__[3];
/*<       etot=etot-delta_e >*/
    etot -= delta_e__;
/* *** CM  beta of fotone+interacting nucleon */
/*<       BETA_N_D(1)=PN_N_D(1,1)/ETOT >*/
    h2_1.beta_n_d__[0] = h2_1.pn_n_d__[0] / etot;
/*<       BETA_N_D(2)=PN_N_D(2,1)/ETOT >*/
    h2_1.beta_n_d__[1] = h2_1.pn_n_d__[1] / etot;
/*<       BETA_N_D(3)=(PN_N_D(3,1)+EGAM)/ETOT >*/
    h2_1.beta_n_d__[2] = (h2_1.pn_n_d__[2] + montecarl_1.egam) / etot;
/*<       BETA_N_D(4)=0. >*/
    h2_1.beta_n_d__[3] = 0.f;
/* *** Variabile per impulso di FErmi nella ntupla */
/* *** Fermi motion is computed in the frame in wich deuteron is at */
/* *** rest, i.e. in LAB */
/*<       p_fermi=abs(P1) >*/
    fermi_motion__1.p_fermi__ = dabs(p1);
/* ************************************************************* */
/* ***   S is Lorentz invariant. Scalar product has only third */
/* ***   component if photon is along z */
/* ************************************************************* */
/* *** */
/* *** Square of (gamma+target nucleon) */
/* *** */
/*<       S = RM(1)**2+2.*PN_N_D(4,1)*EGAM-2.*PN_N_D(3,1)*EGAM >*/
/* Computing 2nd power */
    r__1 = fermi_motion__1.rm[0];
    s = r__1 * r__1 + h2_1.pn_n_d__[3] * 2.f * montecarl_1.egam - 
	    h2_1.pn_n_d__[2] * 2.f * montecarl_1.egam;
/* *** */
/* *** subtracting spectator mass variation */
/* *** */
/*<       S = S +delta_e**2.-2.*delta_e*(egam+PN_N_D(4,1)) >*/
    d__1 = (doublereal) delta_e__;
    s = s + pow_dd(&d__1, &c_b7) - delta_e__ * 2.f * (montecarl_1.egam + 
	    h2_1.pn_n_d__[3]);
/* ********************************************* */
/* *** Check on total energy for the reaction */
/*<       totmass1=rmass('proton      ') >*/
    totmass1 = rmass_("proton      ", (ftnlen)12);
/*<       totmass2=rmass('pi-         ') >*/
    totmass2 = rmass_("pi-         ", (ftnlen)12);
/*<       totmass=(totmass1+totmass2)**2. >*/
    d__1 = (doublereal) (totmass1 + totmass2);
    totmass = pow_dd(&d__1, &c_b7);
/*<       if (s.lt.totmass+0.000001) go to 3 >*/
    if (s < totmass + 1e-6f) {
	goto L3;
    }
/* ********************************************* */
/*<       W_N = SQRT(S) >*/
    *w_n__ = sqrt(s);
/*<       RETURN >*/
    return 0;
/*<       END                                               !*** END FERMI_D >*/
} /* fermi_d__ */

/*<       real function mass_n_d(pf) >*/
doublereal mass_n_d__(real *pf)
{
    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real x, rmn;
    extern doublereal rmass_(char *, ftnlen);

/* ****************************** */
/* *** This function computes reduced mass */
/* *** of a nucleon inside the deuteron */
/* *** from fermi motion so that energy */
/* *** is conserved */
/* ****************************** */
/*<       rmn=rmass('deuteron    ')/2. >*/
    rmn = rmass_("deuteron    ", (ftnlen)12) / 2.f;
/*<       x=sqrt(rmn*rmn-pf*pf) >*/
    x = sqrt(rmn * rmn - *pf * *pf);
/*<       mass_n_d=x >*/
    ret_val = x;
/*<       return >*/
    return ret_val;
/*<       end >*/
} /* mass_n_d__ */

