/* eg.F -- translated by f2c (version 20230428).
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
    integer nerr, loopv, ncaseg, loopc, nproc;
} flag_gen__;

#define flag_gen__1 flag_gen__

struct {
    real pn[15]	/* was [5][3] */;
} mom_nucleon__;

#define mom_nucleon__1 mom_nucleon__

struct {
    real xsect[38000]	/* was [40][50][19] */;
    integer i__;
} sigma_;

#define sigma_1 sigma_

struct {
    real p4_cdm__[3], p3_cdm__[3], phii;
} impul_;

#define impul_1 impul_

struct {
    real rmb1, rmb2, rmme, deltamass, dibmass, deumass, rhomass;
} masses_;

#define masses_1 masses_

struct {
    integer np;
    real ecm, amass[18];
    integer kgenev;
} genin_;

#define genin_1 genin_

struct {
    real pcm[90]	/* was [5][18] */, wt;
} genout_;

#define genout_1 genout_

struct {
    logical eg, cmsf;
} logeg_;

#define logeg_1 logeg_

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;
static integer c__19 = 19;
static integer c__37 = 37;
static integer c__38 = 38;
static integer c__17 = 17;
static integer c__7 = 7;
static integer c__8 = 8;
static integer c__10 = 10;
static integer c__16 = 16;
static integer c__11 = 11;

/* ********************************************************************** */
/*<       SUBROUTINE MONTECARLO >*/
/* Subroutine */ int montecarlo_(void)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    extern /* Subroutine */ int check_cl__(void);
    static real emin, emax, vran[1];
    extern /* Subroutine */ int ev_nuc__(void), ranmar_(real *, integer *);

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* FUNCTIONAL DESCRIPTION: */

/* Selects a proper energy  of gamma to make an hadron interaction, */
/* if electromagnetic event is ordered, do nothing. */

/* MODIFICATION DATE: */

/*     12.05.95       by Igor Pshenichnov */

/*  ********** BY MARCO MIRAZITA */
/*     23.04.99 */
/* --> Aggiunto il caso (IBEAM=3) di uno spettro di energia dei fotoni */
/*     uniforme tra ETAG-FWHM e ETAG+FWHM. Per il caso dello spettro di */
/*     bremstrahhlung con IBEAM=2, l'energia e' estratta tra ETAG-FWHM */
/*     e ETAG anziche' tra 0.151 GeV e ETAG. */


/* COMMON BLOCKS: */

/*     /FLAG_GEN/ */
/*     /FLAG_COMP/ */
/*     /MONTECARL/ */
/* ----------------------------------------------------------------------- */
/*<       INTEGER IBEAM,IELET,NCHAIN,ICHAIN(40),ICIBLE,IPOL,ILAM >*/
/*<       REAL ETAG,FWHM >*/
/*<        >*/
/*<       REAL EGAM,P(4,11) >*/
/*<       INTEGER I(11),NPART >*/
/*<       COMMON/MONTECARL/ EGAM,I,NPART,P >*/
/*<       INTEGER ITARG,IIDH >*/
/*<       CHARACTER*12 TARGET,TARGFERMI >*/
/*<       COMMON /TARGET/ ITARG,IIDH,TARGET,TARGFERMI >*/
/*<       real vran(1) >*/
/*<       IF (IBEAM.LE.1) THEN >*/
    if (flag_comp__1.ibeam <= 1) {
/* ***    Select the energy of tagged gamma quanta in accordance with */
/* ***    average energy and width */
/*<          call ranmar(vran,1) >*/
	ranmar_(vran, &c__1);
/*<          EGAM = (ETAG + .425*FWHM*vran(1)) >*/
	montecarl_1.egam = flag_comp__1.etag + flag_comp__1.fwhm * .425f * 
		vran[0];
/*<       ELSEIF (IBEAM.EQ.2) THEN >*/
    } else if (flag_comp__1.ibeam == 2) {
/* ***    Select the energy of gamma quanta in accordance with bremsstrahlung */
/* ***    spectra. In this case ETAG is the maximum energy and (ETAG-FWHM) is */
/* ***    the minimum energy of such a spectrum. */
/* PURE BREMSSTRAHLUNG  .... */
/*<          Emin=ETAG-FWHM >*/
	emin = flag_comp__1.etag - flag_comp__1.fwhm;
/*<          Emax=ETAG+FWHM >*/
	emax = flag_comp__1.etag + flag_comp__1.fwhm;
/*<          call ranmar(vran,1) >*/
	ranmar_(vran, &c__1);
/*<          EGAM = Emin*(Emax/Emin)**vran(1) >*/
	d__1 = (doublereal) (emax / emin);
	d__2 = (doublereal) vran[0];
	montecarl_1.egam = emin * pow_dd(&d__1, &d__2);
/*<       ELSE >*/
    } else {
/* ***   Select the energy of gamma quanta in accordance with uniform spectrum */
/* ***   between ETAG-FWHM and ETAG+FWHM */
/*<          call ranmar(vran,1) >*/
	ranmar_(vran, &c__1);
/*<          egam=etag+fwhm*(2.*vran(1)-1.) >*/
	montecarl_1.egam = flag_comp__1.etag + flag_comp__1.fwhm * (vran[0] * 
		2.f - 1.f);
/*<       ENDIF >*/
    }
/*<       IF (IELET.EQ.0) THEN >*/
    if (flag_comp__1.ielet == 0) {
/*<           CALL EV_NUC >*/
	ev_nuc__();
/*<           CALL CHECK_CL >*/
	check_cl__();
/*<       ELSE >*/
    } else {
/*<         NPART  = 1 >*/
	montecarl_1.npart = 1;
/*<         P(1,1) = 0. >*/
	montecarl_1.p[0] = 0.f;
/*<         P(2,1) = 0. >*/
	montecarl_1.p[1] = 0.f;
/*<         P(3,1) = EGAM >*/
	montecarl_1.p[2] = montecarl_1.egam;
/*<         I(1)   = 1 >*/
	montecarl_1.i__[0] = 1;
/*<       ENDIF >*/
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* montecarlo_ */

/* **************************************************************** */
/*< 	SUBROUTINE CHECK_CL >*/
/* Subroutine */ int check_cl__(void)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, il;
    static real tol, psqr, delta[4];



/* FUNCTIONAL DESCRIPTION: */

/*     Check energy-momentum conservation for an event presented */
/*     in P array. Nucleon(s) momentum (momenta) presented in PN array. */
/*     If the difference is out of tolerance, NERR=NERR+1. */

/* CREATION DATE: */

/*     22.05.95       by Igor Pshenichnov */

/* COMMON BLOCKS: */

/*     /FLAG_GEN/ */
/*     /MONTECARL/ */
/*     /LOGEG/ */

/* *** E' solo per lo stato iniziale, compreso il moto di Fermi ? */

/*<       INTEGER NERR,LOOPV,NCASEG,LOOPC,NPROC >*/
/*<       COMMON /FLAG_GEN/ NERR,LOOPV,NCASEG,LOOPC,NPROC >*/
/*<       REAL EGAM,P(4,11) >*/
/*<       INTEGER I(11),NPART >*/
/*<       COMMON/MONTECARL/ EGAM,I,NPART,P >*/
/*<       REAL PN(5,3),DELTA(4) >*/
/*<       COMMON /MOM_NUCLEON/PN >*/
/*<       TOL=1.E-04 ! 0.1 MeV accuracy >*/
    tol = 1e-4f;
/*<       IF(NPART.EQ.1) RETURN     ! If there is no hadron interaction, >*/
    if (montecarl_1.npart == 1) {
	return 0;
    }
/* nothing to check there. */
/*< 	DELTA(1)=PN(1,2)+PN(1,3) >*/
    delta[0] = mom_nucleon__1.pn[5] + mom_nucleon__1.pn[10];
/*< 	DELTA(2)=PN(2,2)+PN(2,3) >*/
    delta[1] = mom_nucleon__1.pn[6] + mom_nucleon__1.pn[11];
/*< 	DELTA(3)=-EGAM-PN(3,1) >*/
    delta[2] = -montecarl_1.egam - mom_nucleon__1.pn[2];
/*< 	DELTA(4)=-EGAM-PN(4,1) >*/
    delta[3] = -montecarl_1.egam - mom_nucleon__1.pn[3];
/*<       DO J=1,NPART >*/
    i__1 = montecarl_1.npart;
    for (j = 1; j <= i__1; ++j) {
/*<          PSQR=0.0 >*/
	psqr = 0.f;
/*<          DO IL=1,3 >*/
	for (il = 1; il <= 3; ++il) {
/*<             DELTA(IL)=DELTA(IL)+P(IL,J) >*/
	    delta[il - 1] += montecarl_1.p[il + (j << 2) - 5];
/*<             PSQR=PSQR+P(IL,J)**2  >*/
/* Computing 2nd power */
	    r__1 = montecarl_1.p[il + (j << 2) - 5];
	    psqr += r__1 * r__1;
/*<          ENDDO >*/
	}
/*<          DELTA(4)=DELTA(4)+SQRT(PSQR+P(4,J)**2) >*/
/* Computing 2nd power */
	r__1 = montecarl_1.p[(j << 2) - 1];
	delta[3] += sqrt(psqr + r__1 * r__1);
/*<       ENDDO >*/
    }
/*<        >*/
    if (dabs(delta[0]) > tol || dabs(delta[1]) > tol || dabs(delta[2]) > tol 
	    || dabs(delta[3]) > tol) {
/*<          NERR = NERR+1 >*/
	++flag_gen__1.nerr;
/*<       ENDIF >*/
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* check_cl__ */

/* *********************************************************************** */
/* *********************************************************************** */
/*<       SUBROUTINE GEN_EVT(W_S,LEV,JCH) >*/
/* Subroutine */ int gen_evt__(real *w_s__, integer *lev, integer *jch)
{
    /* System generated locals */
    real r__1, r__2, r__3;

    /* Builtin functions */
    double sqrt(doublereal), sin(doublereal), cos(doublereal);

    /* Local variables */
    static real sigmodel[19], argument;
    static integer j;
    static real w, u1, u2, u3, sigriftot, dw, sr[2], wr[2];
    extern doublereal ali_(real *, real *, real *, integer *);
    static real w_f__, p_m__, w_i__, tet[19], teta, prob[19], vran[1], 
	    phi_m__;
    extern doublereal cosdd_(real *);
    static real dtheta;
    extern /* Subroutine */ int ranmar_(real *, integer *);
    static real sigrif[19], etot_m__, theta_m__;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* FUNCTIONAL DESCRIPTION: */

/*    This subroutine extracts c.m. direction of final state particles, */
/*    accounting of differential cross sections. */

/* MODIFICATION DATE: */

/*     12.05.95       by Igor Pshenichnov */



/* COMMON BLOCKS: */

/*     /SIGMA/ */
/*     /FLAG_COMP/ */
/*     /MASSES/ */
/*     /GENIN/ */
/*     /GENOUT/ */

/* ----------------------------------------------------------------------- */
/*<       IMPLICIT NONE >*/
/*<       REAL XSECT(40,50,0:18) >*/
/*<       COMMON /SIGMA/ XSECT,I >*/
/*<       INTEGER IBEAM,IELET,NCHAIN,ICHAIN(40),ICIBLE,IPOL,ILAM >*/
/*<       REAL ETAG,FWHM >*/
/*<        >*/
/*<       REAL RNDM,ARGUMENT >*/
/*<       DOUBLE PRECISION PI,TWOPI,PIBY2,DEGRAD,RADDEG >*/

/*<       PARAMETER (PI=3.14159265358979324D0) >*/
/*<       PARAMETER (TWOPI=6.28318530717958648D0) >*/
/*<       PARAMETER (PIBY2=1.57079632679489662D0) >*/
/*<       PARAMETER (DEGRAD=0.0174532925199432958D0) >*/
/*<       PARAMETER (RADDEG=57.2957795130823209D0) >*/

/*<       COMMON/IMPUL/P4_CDM(3),P3_CDM(3),PHII >*/
/*<       REAL RMB1,RMB2,RMME,DELTAMASS,DIBMASS,DEUMASS,RHOMASS >*/
/*<       COMMON/MASSES/RMB1,RMB2,RMME,DELTAMASS,DIBMASS,DEUMASS,RHOMASS >*/
/*<       INTEGER NP,KGENEV >*/
/*<       REAL ECM,AMASS,PCM,WT >*/
/*<       COMMON/GENIN/NP,ECM,AMASS(18),KGENEV >*/
/*<       COMMON/GENOUT/PCM(5,18),WT >*/
/*<        >*/
/*<       INTEGER J,LEV,I,JCH,K >*/
/*<       REAL ALI >*/
/*<       real cosdd >*/
/*<       integer jm,jk >*/
/*<       real vran(1) >*/
/*<       DTHETA = 10. >*/
    dtheta = 10.f;
/*<       SIGRIF(0) = 0. >*/
    sigrif[0] = 0.f;
/*<       SIGMODEL(0) = 0. >*/
    sigmodel[0] = 0.f;
/*<       PROB(0) = 0. >*/
    prob[0] = 0.f;
/*<       DO J=0,18 >*/
    for (j = 0; j <= 18; ++j) {
/*<          TET(J) = J*DTHETA >*/
	tet[j] = j * dtheta;
/*<       END DO >*/
    }
/*<       IF (LEV.GT.1) GO TO 10 >*/
    if (*lev > 1) {
	goto L10;
    }
/*<       W_I = 1. >*/
    w_i__ = 1.f;
/*<       W_F = 3.5 >*/
    w_f__ = 3.5f;
/*<       DW = 0.050 >*/
    dw = .05f;
/*<       W = W_S >*/
    w = *w_s__;
/*<       SIGRIFTOT = 0. >*/
    sigriftot = 0.f;
/* *********************************************** */
/* *** Select CM angle */
/* *********************************************** */
/*<       DO I=1,50 >*/
    for (sigma_1.i__ = 1; sigma_1.i__ <= 50; ++sigma_1.i__) {
/*<          IF (W.GE.(W_I+DW*I).AND.W.LT.(W_I+DW*(I+1))) GOTO 1 >*/
	if (w >= w_i__ + dw * sigma_1.i__ && w < w_i__ + dw * (sigma_1.i__ + 
		1)) {
	    goto L1;
	}
/*<       END DO >*/
    }
/*<  1    CONTINUE >*/
L1:
/*      write (15,*)'W=',w */
/*      write (15,*) 'bin in W estratto',i */
/* *************************************** */
/* *** Correzione: I non puo' superare 50 */
/*<       if (I.gt.50) I=50 >*/
    if (sigma_1.i__ > 50) {
	sigma_1.i__ = 50;
    }
/*      write (15,*) 'bin in W corretto',i */
/* *************************************** */
/*<       IF (ILAM.EQ.0) THEN >*/
    if (flag_comp__1.ilam == 0) {
/*<          WR(1) = W_I+DW*I >*/
	wr[0] = w_i__ + dw * sigma_1.i__;
/*<          WR(2) = W_I+DW*(I+1) >*/
	wr[1] = w_i__ + dw * (sigma_1.i__ + 1);
/*         write (15,*) '*** WR=',wr(1),wr(2) */
/*<          SR(1) = XSECT(JCH,I,0) >*/
	sr[0] = sigma_1.xsect[*jch + sigma_1.i__ * 40 - 41];
/* *************************************** */
/* *** Correzione: I non puo' superare 50 */
/*<          if (I.lt.50) then >*/
	if (sigma_1.i__ < 50) {
/*<             SR(2) = XSECT(JCH,I+1,0) >*/
	    sr[1] = sigma_1.xsect[*jch + (sigma_1.i__ + 1) * 40 - 41];
/*<          else >*/
	} else {
/*<             SR(2) = XSECT(JCH,50,0) >*/
	    sr[1] = sigma_1.xsect[*jch + 1959];
/*<          endif >*/
	}
/* *************************************** */
/*         write (15,*) 'J=0 - sigma(teta)=',sr(1),sr(2) */
/*<          SIGMODEL(0) = ALI(W,WR,SR,2)  >*/
	sigmodel[0] = ali_(&w, wr, sr, &c__2);
/*<          DO J=1,18 >*/
	for (j = 1; j <= 18; ++j) {
/*<             SR(1) = XSECT(JCH,I,J) >*/
	    sr[0] = sigma_1.xsect[*jch + (sigma_1.i__ + j * 50) * 40 - 41];
/* *************************************** */
/* *** Correzione: I non puo' superare 50 */
/*<             if (I.lt.50) then >*/
	    if (sigma_1.i__ < 50) {
/*<                SR(2) = XSECT(JCH,I+1,J) >*/
		sr[1] = sigma_1.xsect[*jch + (sigma_1.i__ + 1 + j * 50) * 40 
			- 41];
/*<             else >*/
	    } else {
/*<                SR(2) = XSECT(JCH,50,J) >*/
		sr[1] = sigma_1.xsect[*jch + (j * 50 + 50) * 40 - 41];
/*<             endif >*/
	    }
/* *************************************** */
/*            write (15,*) 'J=',j,' - sigma(teta)=',sr(1),sr(2) */
/*<             SIGMODEL(J) = ALI(W,WR,SR,2)  >*/
	    sigmodel[j] = ali_(&w, wr, sr, &c__2);
/*<        >*/
	    sigrif[j] = (sigmodel[j] + sigmodel[j - 1]) / 2.f * (cosdd_(&tet[
		    j - 1]) - cosdd_(&tet[j]));
/*<             SIGRIFTOT = SIGRIFTOT + SIGRIF(J) >*/
	    sigriftot += sigrif[j];
/*<          END DO >*/
	}
/*         write (15,*) 'sigma mod.=',sigmodel */
/*         write (15,*) 'sigma rif.=',sigrif */
/*         write (15,*) 'SIGRIFTOT=',SIGRIFTOT */
/*<          PROB(0) = 0. >*/
	prob[0] = 0.f;
/*<          DO J=1,18 >*/
	for (j = 1; j <= 18; ++j) {
/*<             SIGRIF(J) = SIGRIF(J)/SIGRIFTOT >*/
	    sigrif[j] /= sigriftot;
/*<             PROB(J) = PROB(J-1)+SIGRIF(J) >*/
	    prob[j] = prob[j - 1] + sigrif[j];
/*<          END DO >*/
	}
/*         write (15,*) 'prob=',prob */
/*<          call ranmar(vran,1) >*/
	ranmar_(vran, &c__1);
/*<          TETA=ALI(vran(1),PROB,TET,19) >*/
	teta = ali_(vran, prob, tet, &c__19);
/*         write (15,*) 'teta=',teta */
/*<       ENDIF >*/
    }
/*<       THETA_M = TETA*DEGRAD >*/
    theta_m__ = teta * .0174532925199432958;
/* *** Estrazione uniforme di fi */
/*<       call ranmar(vran,1) >*/
    ranmar_(vran, &c__1);
/*<       PHII = TWOPI*vran(1) >*/
    impul_1.phii = vran[0] * 6.28318530717958648;
/*<       PHI_M = PHII >*/
    phi_m__ = impul_1.phii;
/* ********************************************* */
/* *** Final state kinematics */
/* ********************************************* */
/*<       ETOT_M = (W_S**2+RMME**2-RMB2**2)/(2*W_S) >*/
/* Computing 2nd power */
    r__1 = *w_s__;
/* Computing 2nd power */
    r__2 = masses_1.rmme;
/* Computing 2nd power */
    r__3 = masses_1.rmb2;
    etot_m__ = (r__1 * r__1 + r__2 * r__2 - r__3 * r__3) / (*w_s__ * 2);
/*<       ARGUMENT = ETOT_M**2-RMME**2 >*/
/* Computing 2nd power */
    r__1 = etot_m__;
/* Computing 2nd power */
    r__2 = masses_1.rmme;
    argument = r__1 * r__1 - r__2 * r__2;
/*<       IF (ARGUMENT.GT.0.0) THEN >*/
    if (argument > 0.f) {
/*<          P_M=SQRT(ARGUMENT) >*/
	p_m__ = sqrt(argument);
/*<       ELSE >*/
    } else {
/*<          P_M=0.0 >*/
	p_m__ = 0.f;
/*<       ENDIF >*/
    }
/*<       U1=SIN(THETA_M)*COS(PHI_M) >*/
    u1 = sin(theta_m__) * cos(phi_m__);
/*<       U2=SIN(THETA_M)*SIN(PHI_M) >*/
    u2 = sin(theta_m__) * sin(phi_m__);
/*<       U3=COS(THETA_M) >*/
    u3 = cos(theta_m__);
/* *** Mesone prodotto */
/*<       PCM(1,1)=P_M*U1 >*/
    genout_1.pcm[0] = p_m__ * u1;
/*<       PCM(2,1)=P_M*U2 >*/
    genout_1.pcm[1] = p_m__ * u2;
/*<       PCM(3,1)=P_M*U3 >*/
    genout_1.pcm[2] = p_m__ * u3;
/*<       PCM(4,1)=ETOT_M >*/
    genout_1.pcm[3] = etot_m__;
/*<       PCM(5,1)=P_M >*/
    genout_1.pcm[4] = p_m__;
/*      write (*,15) (pcm(k,1),k=1,5) */
/*<  15   format (1x,'Mes.',5(1x,f8.5)) >*/
/* L15: */
/* *** Barione diffuso */
/*<       PCM(1,2)=-PCM(1,1) >*/
    genout_1.pcm[5] = -genout_1.pcm[0];
/*<       PCM(2,2)=-PCM(2,1) >*/
    genout_1.pcm[6] = -genout_1.pcm[1];
/*<       PCM(3,2)=-PCM(3,1) >*/
    genout_1.pcm[7] = -genout_1.pcm[2];
/*<       PCM(4,2)=W_S-ETOT_M >*/
    genout_1.pcm[8] = *w_s__ - etot_m__;
/*<       PCM(5,2)=P_M >*/
    genout_1.pcm[9] = p_m__;
/*      write (*,16) (pcm(k,2),k=1,5) */
/*<  16   format (1x,'Bar.',5(1x,f8.5)) >*/
/* L16: */
/* *** */
/*<       P4_CDM(1) = PCM(1,1) >*/
    impul_1.p4_cdm__[0] = genout_1.pcm[0];
/*<       P4_CDM(2) = PCM(2,1) >*/
    impul_1.p4_cdm__[1] = genout_1.pcm[1];
/*<       P4_CDM(3) = PCM(3,1) >*/
    impul_1.p4_cdm__[2] = genout_1.pcm[2];
/*<       P3_CDM(1) = PCM(1,2) >*/
    impul_1.p3_cdm__[0] = genout_1.pcm[5];
/*<       P3_CDM(2) = PCM(2,2) >*/
    impul_1.p3_cdm__[1] = genout_1.pcm[6];
/*<       P3_CDM(3) = PCM(3,2) >*/
    impul_1.p3_cdm__[2] = genout_1.pcm[7];
/*< 10      RETURN >*/
L10:
    return 0;
/*<       END                                               !*** END GEN_EVT >*/
} /* gen_evt__ */

/* *********************************************************************** */
/* *********************************************************************** */
/*<       REAL FUNCTION RMASS(P) >*/
doublereal rmass_(char *p, ftnlen p_len)
{
    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static real q;
    extern doublereal amas_(integer *), amrez_(integer *);
    extern /* Subroutine */ int norran_(real *);

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* FUNCTIONAL DESCRIPTION: */

/*    Put a mass of particle as it was adopted by INC code, including */
/*    some resonances */

/* MODIFICATION DATE */

/*     20.05.95       by Igor Pshenichnov */


/* FUNCTION VALUE: */

/*    RMASS   (GeV) */

/* ----------------------------------------------------------------------- */
/*<       IMPLICIT NONE >*/
/*<       CHARACTER*12 P >*/
/*<       REAL Q,AMAS,AMREZ >*/
/*<       CALL NORRAN(Q) >*/
    norran_(&q);
/*<       RMASS = -1. >*/
    ret_val = -1.f;
/*<       IF(P.EQ.'gamma')  THEN       >*/
    if (s_cmp(p, "gamma", (ftnlen)12, (ftnlen)5) == 0) {
/*<                      RMASS = 0.0 >*/
	ret_val = 0.f;
/*<       ELSEIF(P.EQ.'neutrino')  THEN    >*/
    } else if (s_cmp(p, "neutrino", (ftnlen)12, (ftnlen)8) == 0) {
/*<                      RMASS = 0.0 >*/
	ret_val = 0.f;
/*<       ELSEIF(P.EQ.'nulla')  THEN       >*/
    } else if (s_cmp(p, "nulla", (ftnlen)12, (ftnlen)5) == 0) {
/*<                      RMASS = 0.0 >*/
	ret_val = 0.f;
/*<       ELSEIF(P.EQ.'electron')  THEN    >*/
    } else if (s_cmp(p, "electron", (ftnlen)12, (ftnlen)8) == 0) {
/*<                      RMASS = 0.511E-03 >*/
	ret_val = 5.11e-4f;
/*<       ELSEIF(P.EQ.'positron')  THEN    >*/
    } else if (s_cmp(p, "positron", (ftnlen)12, (ftnlen)8) == 0) {
/*<                      RMASS = 0.511E-03 >*/
	ret_val = 5.11e-4f;
/*<       ELSEIF(P.EQ.'mu+')  THEN         >*/
    } else if (s_cmp(p, "mu+", (ftnlen)12, (ftnlen)3) == 0) {
/*<                      RMASS = 0.10566 >*/
	ret_val = .10566f;
/*<       ELSEIF(P.EQ.'mu-')  THEN         >*/
    } else if (s_cmp(p, "mu-", (ftnlen)12, (ftnlen)3) == 0) {
/*<                      RMASS = 0.10566 >*/
	ret_val = .10566f;
/*<       ELSEIF(P.EQ.'proton')  THEN      >*/
    } else if (s_cmp(p, "proton", (ftnlen)12, (ftnlen)6) == 0) {
/*<                      RMASS = AMAS(37) >*/
	ret_val = amas_(&c__37);
/*<       ELSEIF(P.EQ.'neutron')  THEN       >*/
    } else if (s_cmp(p, "neutron", (ftnlen)12, (ftnlen)7) == 0) {
/*<                      RMASS = AMAS(38) >*/
	ret_val = amas_(&c__38);
/*<       ELSEIF(P.EQ.'lambda')  THEN      >*/
    } else if (s_cmp(p, "lambda", (ftnlen)12, (ftnlen)6) == 0) {
/*<                      RMASS = 1.116 >*/
	ret_val = 1.116f;
/*<       ELSEIF(P.EQ.'sigma+')  THEN      >*/
    } else if (s_cmp(p, "sigma+", (ftnlen)12, (ftnlen)6) == 0) {
/*<                      RMASS = 1.189 >*/
	ret_val = 1.189f;
/*<       ELSEIF(P.EQ.'sigma0')  THEN      >*/
    } else if (s_cmp(p, "sigma0", (ftnlen)12, (ftnlen)6) == 0) {
/*<                      RMASS = 1.192 >*/
	ret_val = 1.192f;
/*<       ELSEIF(P.EQ.'sigma-')  THEN      >*/
    } else if (s_cmp(p, "sigma-", (ftnlen)12, (ftnlen)6) == 0) {
/*<                      RMASS = 1.197 >*/
	ret_val = 1.197f;
/*<       ELSEIF(P.EQ.'csi0')  THEN        >*/
    } else if (s_cmp(p, "csi0", (ftnlen)12, (ftnlen)4) == 0) {
/*<                      RMASS = 1.3149 >*/
	ret_val = 1.3149f;
/*<       ELSEIF(P.EQ.'csi-')  THEN        >*/
    } else if (s_cmp(p, "csi-", (ftnlen)12, (ftnlen)4) == 0) {
/*<                      RMASS = 1.3213 >*/
	ret_val = 1.3213f;
/*<       ELSEIF(P.EQ.'dibaryon')  THEN    >*/
    } else if (s_cmp(p, "dibaryon", (ftnlen)12, (ftnlen)8) == 0) {
/*<        >*/
	ret_val = q * .0165f / 2.355f + 2.23f;
/*<       ELSEIF(P.EQ.'delta++')  THEN     >*/
    } else if (s_cmp(p, "delta++", (ftnlen)12, (ftnlen)7) == 0) {
/*<        >*/
/* Computing MAX */
	r__1 = 1.08f, r__2 = q * .115f / 2.355f + 1.235f;
	ret_val = dmax(r__1,r__2);
/*<       ELSEIF(P.EQ.'delta+')  THEN      >*/
    } else if (s_cmp(p, "delta+", (ftnlen)12, (ftnlen)6) == 0) {
/*<        >*/
/* Computing MAX */
	r__1 = 1.08f, r__2 = q * .115f / 2.355f + 1.235f;
	ret_val = dmax(r__1,r__2);
/*<       ELSEIF(P.EQ.'delta0')  THEN      >*/
    } else if (s_cmp(p, "delta0", (ftnlen)12, (ftnlen)6) == 0) {
/*<        >*/
/* Computing MAX */
	r__1 = 1.08f, r__2 = q * .115f / 2.355f + 1.235f;
	ret_val = dmax(r__1,r__2);
/*<       ELSEIF(P.EQ.'delta-')  THEN      >*/
    } else if (s_cmp(p, "delta-", (ftnlen)12, (ftnlen)6) == 0) {
/*<        >*/
/* Computing MAX */
	r__1 = 1.08f, r__2 = q * .115f / 2.355f + 1.235f;
	ret_val = dmax(r__1,r__2);
/*<       ELSEIF(P.EQ.'sigma+3/2')  THEN   >*/
    } else if (s_cmp(p, "sigma+3/2", (ftnlen)12, (ftnlen)9) == 0) {
/*<                      RMASS = 1.382 >*/
	ret_val = 1.382f;
/*<       ELSEIF(P.EQ.'sigma03/2')  THEN   >*/
    } else if (s_cmp(p, "sigma03/2", (ftnlen)12, (ftnlen)9) == 0) {
/*<                      RMASS = 1.387 >*/
	ret_val = 1.387f;
/*<       ELSEIF(P.EQ.'sigma-3/2')  THEN   >*/
    } else if (s_cmp(p, "sigma-3/2", (ftnlen)12, (ftnlen)9) == 0) {
/*<                      RMASS = 1.382 >*/
	ret_val = 1.382f;
/*<       ELSEIF(P.EQ.'csi03/2')  THEN     >*/
    } else if (s_cmp(p, "csi03/2", (ftnlen)12, (ftnlen)7) == 0) {
/*<                      RMASS = 1.532 >*/
	ret_val = 1.532f;
/*<       ELSEIF(P.EQ.'csi-3/2')  THEN     >*/
    } else if (s_cmp(p, "csi-3/2", (ftnlen)12, (ftnlen)7) == 0) {
/*<                      RMASS = 1.535 >*/
	ret_val = 1.535f;
/*<       ELSEIF(P.EQ.'omega-')  THEN      >*/
    } else if (s_cmp(p, "omega-", (ftnlen)12, (ftnlen)6) == 0) {
/*<                      RMASS = 1.67245 >*/
	ret_val = 1.67245f;
/*<       ELSEIF(P.EQ.'omega')  THEN       >*/
    } else if (s_cmp(p, "omega", (ftnlen)12, (ftnlen)5) == 0) {
/*<                      RMASS = AMREZ(17) >*/
	ret_val = amrez_(&c__17);
/*<       ELSEIF(P.EQ.'pi+')  THEN         >*/
    } else if (s_cmp(p, "pi+", (ftnlen)12, (ftnlen)3) == 0) {
/*<                      RMASS = AMAS(1) >*/
	ret_val = amas_(&c__1);
/*<       ELSEIF(P.EQ.'pi0')  THEN         >*/
    } else if (s_cmp(p, "pi0", (ftnlen)12, (ftnlen)3) == 0) {
/*<                      RMASS = AMAS(7) >*/
	ret_val = amas_(&c__7);
/*<       ELSEIF(P.EQ.'pi-')  THEN         >*/
    } else if (s_cmp(p, "pi-", (ftnlen)12, (ftnlen)3) == 0) {
/*<                      RMASS = AMAS(2) >*/
	ret_val = amas_(&c__2);
/*<       ELSEIF(P.EQ.'kappa+')  THEN      >*/
    } else if (s_cmp(p, "kappa+", (ftnlen)12, (ftnlen)6) == 0) {
/*<                      RMASS = 0.49367 >*/
	ret_val = .49367f;
/*<       ELSEIF(P.EQ.'kappa0')  THEN      >*/
    } else if (s_cmp(p, "kappa0", (ftnlen)12, (ftnlen)6) == 0) {
/*<                      RMASS = 0.49772 >*/
	ret_val = .49772f;
/*<       ELSEIF(P.EQ.'kappa-')  THEN      >*/
    } else if (s_cmp(p, "kappa-", (ftnlen)12, (ftnlen)6) == 0) {
/*<                      RMASS = 0.49367 >*/
	ret_val = .49367f;
/*<       ELSEIF(P.EQ.'kappas')  THEN      >*/
    } else if (s_cmp(p, "kappas", (ftnlen)12, (ftnlen)6) == 0) {
/*<                      RMASS = 0.49772 >*/
	ret_val = .49772f;
/*<       ELSEIF(P.EQ.'kappal')  THEN      >*/
    } else if (s_cmp(p, "kappal", (ftnlen)12, (ftnlen)6) == 0) {
/*<                      RMASS = 0.49772 >*/
	ret_val = .49772f;
/*<       ELSEIF(P.EQ.'akappa0')  THEN     >*/
    } else if (s_cmp(p, "akappa0", (ftnlen)12, (ftnlen)7) == 0) {
/*<                      RMASS = 0.49772 >*/
	ret_val = .49772f;
/*<       ELSEIF(P.EQ.'eta')  THEN         >*/
    } else if (s_cmp(p, "eta", (ftnlen)12, (ftnlen)3) == 0) {
/*<                      RMASS = AMAS(8) >*/
	ret_val = amas_(&c__8);
/*<       ELSEIF(P.EQ.'etap')  THEN        >*/
    } else if (s_cmp(p, "etap", (ftnlen)12, (ftnlen)4) == 0) {
/*<                      RMASS = 0.95747 >*/
	ret_val = .95747f;
/*<       ELSEIF(P.EQ.'rho+')  THEN        >*/
    } else if (s_cmp(p, "rho+", (ftnlen)12, (ftnlen)4) == 0) {
/*<                      RMASS = AMREZ(10) >*/
	ret_val = amrez_(&c__10);
/*<       ELSEIF(P.EQ.'rho0')  THEN        >*/
    } else if (s_cmp(p, "rho0", (ftnlen)12, (ftnlen)4) == 0) {
/*<                      RMASS = AMREZ(16) >*/
	ret_val = amrez_(&c__16);
/*<       ELSEIF(P.EQ.'rho-')  THEN        >*/
    } else if (s_cmp(p, "rho-", (ftnlen)12, (ftnlen)4) == 0) {
/*<                      RMASS = AMREZ(11) >*/
	ret_val = amrez_(&c__11);
/*<       ELSEIF(P.EQ.'phi')  THEN         >*/
    } else if (s_cmp(p, "phi", (ftnlen)12, (ftnlen)3) == 0) {
/*<                      RMASS = 1.020 >*/
	ret_val = 1.02f;
/*<       ELSEIF(P.EQ.'deuteron')  THEN    >*/
    } else if (s_cmp(p, "deuteron", (ftnlen)12, (ftnlen)8) == 0) {
/*<                      RMASS = 1.8755 >*/
	ret_val = 1.8755f;
/*<       ELSEIF(P.EQ.'deutheavy')  THEN   >*/
    } else if (s_cmp(p, "deutheavy", (ftnlen)12, (ftnlen)9) == 0) {
/*<                      RMASS = 1.878 >*/
	ret_val = 1.878f;
/*<       ELSEIF(P.EQ.'H3')  THEN          >*/
    } else if (s_cmp(p, "H3", (ftnlen)12, (ftnlen)2) == 0) {
/*<                      RMASS = 2.80875 >*/
	ret_val = 2.80875f;
/*<       ELSEIF(P.EQ.'He3')  THEN         >*/
    } else if (s_cmp(p, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<                      RMASS = 2.80822 >*/
	ret_val = 2.80822f;
/*<       ELSEIF(P.EQ.'He4')  THEN         >*/
    } else if (s_cmp(p, "He4", (ftnlen)12, (ftnlen)3) == 0) {
/*<                      RMASS = 3.72715 >*/
	ret_val = 3.72715f;
/*<       ENDIF >*/
    }
/*<       RETURN >*/
    return ret_val;
/*<       END                                                 !*** END RMASS >*/
} /* rmass_ */

/* *********************************************************************** */
/*<       INTEGER FUNCTION CODE(P) >*/
integer code_(char *p, ftnlen p_len)
{
    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* FUNCTIONAL DESCRIPTION: */

/*    Particle numbering according to old table (used by INC) or to GEANT */

/* MODIFICATION DATE */

/*     21.05.95       by Igor Pshenichnov */

/* FUNCTION VALUE: */

/*    CODE */

/* COMMON BLOCK: */

/*    LOGEG   -  select a proper mood for numbering */

/* ----------------------------------------------------------------------- */
/*<       IMPLICIT NONE >*/
/*<       CHARACTER*12 P >*/
/*<       LOGICAL EG,CMSF >*/
/*<       COMMON/LOGEG/EG,CMSF >*/
/*<       CODE = 0. >*/
    ret_val = 0.f;
/*<                    IF(EG) THEN >*/
    if (logeg_1.eg) {
/*<       IF(P.EQ.'gamma') THEN >*/
	if (s_cmp(p, "gamma", (ftnlen)12, (ftnlen)5) == 0) {
/*<            CODE = 1 >*/
	    ret_val = 1;
/*<       ELSEIF(P.EQ.'neutrino') THEN >*/
	} else if (s_cmp(p, "neutrino", (ftnlen)12, (ftnlen)8) == 0) {
/*<            CODE = 4 >*/
	    ret_val = 4;
/*<       ELSEIF(P.EQ.'nulla') THEN       >*/
	} else if (s_cmp(p, "nulla", (ftnlen)12, (ftnlen)5) == 0) {
/*<            CODE = 0 >*/
	    ret_val = 0;
/*<       ELSEIF(P.EQ.'electron') THEN    >*/
	} else if (s_cmp(p, "electron", (ftnlen)12, (ftnlen)8) == 0) {
/*<            CODE = 3 >*/
	    ret_val = 3;
/*<       ELSEIF(P.EQ.'positron') THEN    >*/
	} else if (s_cmp(p, "positron", (ftnlen)12, (ftnlen)8) == 0) {
/*<            CODE = 2 >*/
	    ret_val = 2;
/*<       ELSEIF(P.EQ.'mu+') THEN         >*/
	} else if (s_cmp(p, "mu+", (ftnlen)12, (ftnlen)3) == 0) {
/*<            CODE = 5 >*/
	    ret_val = 5;
/*<       ELSEIF(P.EQ.'mu-') THEN         >*/
	} else if (s_cmp(p, "mu-", (ftnlen)12, (ftnlen)3) == 0) {
/*<            CODE = 6 >*/
	    ret_val = 6;
/*<       ELSEIF(P.EQ.'proton') THEN      >*/
	} else if (s_cmp(p, "proton", (ftnlen)12, (ftnlen)6) == 0) {
/*<            CODE = 14 >*/
	    ret_val = 14;
/*<       ELSEIF(P.EQ.'neutron') THEN     >*/
	} else if (s_cmp(p, "neutron", (ftnlen)12, (ftnlen)7) == 0) {
/*<            CODE = 13 >*/
	    ret_val = 13;
/*<       ELSEIF(P.EQ.'lambda') THEN      >*/
	} else if (s_cmp(p, "lambda", (ftnlen)12, (ftnlen)6) == 0) {
/*<            CODE = 18 >*/
	    ret_val = 18;
/*<       ELSEIF(P.EQ.'sigma+') THEN      >*/
	} else if (s_cmp(p, "sigma+", (ftnlen)12, (ftnlen)6) == 0) {
/*<            CODE = 19 >*/
	    ret_val = 19;
/*<       ELSEIF(P.EQ.'sigma0') THEN      >*/
	} else if (s_cmp(p, "sigma0", (ftnlen)12, (ftnlen)6) == 0) {
/*<            CODE = 20 >*/
	    ret_val = 20;
/*<       ELSEIF(P.EQ.'sigma-') THEN      >*/
	} else if (s_cmp(p, "sigma-", (ftnlen)12, (ftnlen)6) == 0) {
/*<            CODE = 21 >*/
	    ret_val = 21;
/*<       ELSEIF(P.EQ.'csi0') THEN        >*/
	} else if (s_cmp(p, "csi0", (ftnlen)12, (ftnlen)4) == 0) {
/*<            CODE = 22 >*/
	    ret_val = 22;
/*<       ELSEIF(P.EQ.'csi-') THEN        >*/
	} else if (s_cmp(p, "csi-", (ftnlen)12, (ftnlen)4) == 0) {
/*<            CODE = 23 >*/
	    ret_val = 23;
/*<       ELSEIF(P.EQ.'dibaryon') THEN    >*/
	} else if (s_cmp(p, "dibaryon", (ftnlen)12, (ftnlen)8) == 0) {
/*<            CODE = 109 >*/
	    ret_val = 109;
/*<       ELSEIF(P.EQ.'delta') THEN       >*/
	} else if (s_cmp(p, "delta", (ftnlen)12, (ftnlen)5) == 0) {
/*<            CODE = 101 >*/
	    ret_val = 101;
/*<       ELSEIF(P.EQ.'delta++') THEN     >*/
	} else if (s_cmp(p, "delta++", (ftnlen)12, (ftnlen)7) == 0) {
/*<            CODE = 101 >*/
	    ret_val = 101;
/*<       ELSEIF(P.EQ.'delta+') THEN      >*/
	} else if (s_cmp(p, "delta+", (ftnlen)12, (ftnlen)6) == 0) {
/*<            CODE = 102 >*/
	    ret_val = 102;
/*<       ELSEIF(P.EQ.'delta0') THEN      >*/
	} else if (s_cmp(p, "delta0", (ftnlen)12, (ftnlen)6) == 0) {
/*<            CODE = 103 >*/
	    ret_val = 103;
/*<       ELSEIF(P.EQ.'delta-') THEN      >*/
	} else if (s_cmp(p, "delta-", (ftnlen)12, (ftnlen)6) == 0) {
/*<            CODE = 104 >*/
	    ret_val = 104;
/*<       ELSEIF(P.EQ.'sigma+3/2') THEN   >*/
	} else if (s_cmp(p, "sigma+3/2", (ftnlen)12, (ftnlen)9) == 0) {
/*<            CODE = 110 >*/
	    ret_val = 110;
/*<       ELSEIF(P.EQ.'sigma03/2') THEN   >*/
	} else if (s_cmp(p, "sigma03/2", (ftnlen)12, (ftnlen)9) == 0) {
/*<            CODE = 111 >*/
	    ret_val = 111;
/*<       ELSEIF(P.EQ.'sigma-3/2') THEN   >*/
	} else if (s_cmp(p, "sigma-3/2", (ftnlen)12, (ftnlen)9) == 0) {
/*<            CODE = 112 >*/
	    ret_val = 112;
/*<       ELSEIF(P.EQ.'csi03/2') THEN     >*/
	} else if (s_cmp(p, "csi03/2", (ftnlen)12, (ftnlen)7) == 0) {
/*<            CODE = 114 >*/
	    ret_val = 114;
/*<       ELSEIF(P.EQ.'csi-3/2') THEN     >*/
	} else if (s_cmp(p, "csi-3/2", (ftnlen)12, (ftnlen)7) == 0) {
/*<            CODE = 115 >*/
	    ret_val = 115;
/*<       ELSEIF(P.EQ.'omega-') THEN     !Baryon >*/
	} else if (s_cmp(p, "omega-", (ftnlen)12, (ftnlen)6) == 0) {
/*<            CODE = 24 >*/
	    ret_val = 24;
/*<       ELSEIF(P.EQ.'omega') THEN      !Meson GSIM >*/
	} else if (s_cmp(p, "omega", (ftnlen)12, (ftnlen)5) == 0) {
/*<            CODE = 60 >*/
	    ret_val = 60;
/*<       ELSEIF(P.EQ.'pi+') THEN         >*/
	} else if (s_cmp(p, "pi+", (ftnlen)12, (ftnlen)3) == 0) {
/*<            CODE = 8 >*/
	    ret_val = 8;
/*<       ELSEIF(P.EQ.'pi0') THEN         >*/
	} else if (s_cmp(p, "pi0", (ftnlen)12, (ftnlen)3) == 0) {
/*<            CODE = 7 >*/
	    ret_val = 7;
/*<       ELSEIF(P.EQ.'pi-') THEN         >*/
	} else if (s_cmp(p, "pi-", (ftnlen)12, (ftnlen)3) == 0) {
/*<            CODE = 9 >*/
	    ret_val = 9;
/*<       ELSEIF(P.EQ.'kappa+') THEN      >*/
	} else if (s_cmp(p, "kappa+", (ftnlen)12, (ftnlen)6) == 0) {
/*<            CODE = 11 >*/
	    ret_val = 11;
/*<       ELSEIF(P.EQ.'kappa0') THEN      >*/
	} else if (s_cmp(p, "kappa0", (ftnlen)12, (ftnlen)6) == 0) {
/*<            CODE = 211 >*/
	    ret_val = 211;
/*<       ELSEIF(P.EQ.'kappa-') THEN      >*/
	} else if (s_cmp(p, "kappa-", (ftnlen)12, (ftnlen)6) == 0) {
/*<            CODE = 12 >*/
	    ret_val = 12;
/*<       ELSEIF(P.EQ.'kappas') THEN      >*/
	} else if (s_cmp(p, "kappas", (ftnlen)12, (ftnlen)6) == 0) {
/*<            CODE = 16 >*/
	    ret_val = 16;
/*<       ELSEIF(P.EQ.'kappal') THEN      >*/
	} else if (s_cmp(p, "kappal", (ftnlen)12, (ftnlen)6) == 0) {
/*<            CODE = 10 >*/
	    ret_val = 10;
/*<       ELSEIF(P.EQ.'akappa0') THEN     >*/
	} else if (s_cmp(p, "akappa0", (ftnlen)12, (ftnlen)7) == 0) {
/*<            CODE = 210 >*/
	    ret_val = 210;
/*<       ELSEIF(P.EQ.'eta') THEN         >*/
	} else if (s_cmp(p, "eta", (ftnlen)12, (ftnlen)3) == 0) {
/*<            CODE = 17 >*/
	    ret_val = 17;
/*<       ELSEIF(P.EQ.'etap') THEN       ! GSIM >*/
	} else if (s_cmp(p, "etap", (ftnlen)12, (ftnlen)4) == 0) {
/*<            CODE = 61 >*/
	    ret_val = 61;
/*<       ELSEIF(P.EQ.'rho+') THEN       ! GSIM >*/
	} else if (s_cmp(p, "rho+", (ftnlen)12, (ftnlen)4) == 0) {
/*<            CODE = 58 >*/
	    ret_val = 58;
/*<       ELSEIF(P.EQ.'rho0') THEN       ! GSIM >*/
	} else if (s_cmp(p, "rho0", (ftnlen)12, (ftnlen)4) == 0) {
/*<            CODE = 57 >*/
	    ret_val = 57;
/*<       ELSEIF(P.EQ.'rho-') THEN       ! GSIM >*/
	} else if (s_cmp(p, "rho-", (ftnlen)12, (ftnlen)4) == 0) {
/*<            CODE = 59 >*/
	    ret_val = 59;
/*<       ELSEIF(P.EQ.'phi') THEN         >*/
	} else if (s_cmp(p, "phi", (ftnlen)12, (ftnlen)3) == 0) {
/*<            CODE = 117 >*/
	    ret_val = 117;
/*<       ELSEIF(P.EQ.'deuteron') THEN    >*/
	} else if (s_cmp(p, "deuteron", (ftnlen)12, (ftnlen)8) == 0) {
/*<            CODE = 45 >*/
	    ret_val = 45;
/*<       ELSEIF(P.EQ.'deutheavy') THEN   >*/
	} else if (s_cmp(p, "deutheavy", (ftnlen)12, (ftnlen)9) == 0) {
/*<            CODE = 245 >*/
	    ret_val = 245;
/*<       ELSEIF(P.EQ.'H3') THEN          >*/
	} else if (s_cmp(p, "H3", (ftnlen)12, (ftnlen)2) == 0) {
/*<            CODE = 46 >*/
	    ret_val = 46;
/*<       ELSEIF(P.EQ.'He3') THEN         >*/
	} else if (s_cmp(p, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<            CODE = 49 >*/
	    ret_val = 49;
/*<       ELSEIF(P.EQ.'He4') THEN         >*/
	} else if (s_cmp(p, "He4", (ftnlen)12, (ftnlen)3) == 0) {
/*<            CODE = 47 >*/
	    ret_val = 47;
/*<       ENDIF >*/
	}
/*<                    ELSE >*/
    } else {
/*<       IF(P.EQ.'gamma') THEN       >*/
	if (s_cmp(p, "gamma", (ftnlen)12, (ftnlen)5) == 0) {
/*<            CODE = 55 >*/
	    ret_val = 55;
/*<       ELSEIF(P.EQ.'proton') THEN      >*/
	} else if (s_cmp(p, "proton", (ftnlen)12, (ftnlen)6) == 0) {
/*<            CODE = 37 >*/
	    ret_val = 37;
/*<       ELSEIF(P.EQ.'neutron') THEN     >*/
	} else if (s_cmp(p, "neutron", (ftnlen)12, (ftnlen)7) == 0) {
/*<            CODE = 38 >*/
	    ret_val = 38;
/*<       ELSEIF(P.EQ.'omega') THEN       >*/
	} else if (s_cmp(p, "omega", (ftnlen)12, (ftnlen)5) == 0) {
/*<            CODE = 17 >*/
	    ret_val = 17;
/*<       ELSEIF(P.EQ.'pi+') THEN         >*/
	} else if (s_cmp(p, "pi+", (ftnlen)12, (ftnlen)3) == 0) {
/*<            CODE = 1 >*/
	    ret_val = 1;
/*<       ELSEIF(P.EQ.'pi0') THEN         >*/
	} else if (s_cmp(p, "pi0", (ftnlen)12, (ftnlen)3) == 0) {
/*<            CODE = 7 >*/
	    ret_val = 7;
/*<       ELSEIF(P.EQ.'pi-') THEN         >*/
	} else if (s_cmp(p, "pi-", (ftnlen)12, (ftnlen)3) == 0) {
/*<            CODE = 2 >*/
	    ret_val = 2;
/*<       ELSEIF(P.EQ.'eta') THEN         >*/
	} else if (s_cmp(p, "eta", (ftnlen)12, (ftnlen)3) == 0) {
/*<            CODE = 8 >*/
	    ret_val = 8;
/*<       ELSEIF(P.EQ.'rho+') THEN        >*/
	} else if (s_cmp(p, "rho+", (ftnlen)12, (ftnlen)4) == 0) {
/*<            CODE = 10 >*/
	    ret_val = 10;
/*<       ELSEIF(P.EQ.'rho0') THEN        >*/
	} else if (s_cmp(p, "rho0", (ftnlen)12, (ftnlen)4) == 0) {
/*<            CODE = 16 >*/
	    ret_val = 16;
/*<       ELSEIF(P.EQ.'rho-') THEN        >*/
	} else if (s_cmp(p, "rho-", (ftnlen)12, (ftnlen)4) == 0) {
/*<            CODE = 11 >*/
	    ret_val = 11;
/*<       ELSEIF(P.EQ.'kappa+') THEN      >*/
	} else if (s_cmp(p, "kappa+", (ftnlen)12, (ftnlen)6) == 0) {
/*<            CODE = 3 >*/
	    ret_val = 3;
/*<       ELSEIF(P.EQ.'kappa0') THEN      >*/
	} else if (s_cmp(p, "kappa0", (ftnlen)12, (ftnlen)6) == 0) {
/*<            CODE = 5 >*/
	    ret_val = 5;
/*<       ELSEIF(P.EQ.'kappa-') THEN      >*/
	} else if (s_cmp(p, "kappa-", (ftnlen)12, (ftnlen)6) == 0) {
/*<            CODE = 4 >*/
	    ret_val = 4;
/*<       ELSEIF(P.EQ.'akappa0') THEN     >*/
	} else if (s_cmp(p, "akappa0", (ftnlen)12, (ftnlen)7) == 0) {
/*<            CODE = 6 >*/
	    ret_val = 6;
/*<       ELSEIF(P.EQ.'phi') THEN         >*/
	} else if (s_cmp(p, "phi", (ftnlen)12, (ftnlen)3) == 0) {
/*<            CODE = 18 >*/
	    ret_val = 18;
/*<       ENDIF >*/
	}
/*<                    ENDIF >*/
    }
/*<       RETURN >*/
    return ret_val;
/*<       END                                                  !*** END CODE >*/
} /* code_ */

/* ********************************************* */
/* ********************************************* */
/* *** Trigonometric functions for angle in degrees */
/* ********************************************* */
/*<       real function acosdd(x) >*/
doublereal acosdd_(real *x)
{
    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    double acos(doublereal);

    /* Local variables */
    static real ac, pi;

/*<       real x,ac >*/
/*<       pi=acos(-1.) >*/
    pi = acos(-1.f);
/*<       ac=acos(x) >*/
    ac = acos(*x);
/*<       acosdd=ac*180./pi >*/
    ret_val = ac * 180.f / pi;
/*<       return >*/
    return ret_val;
/*<       end >*/
} /* acosdd_ */

/* ********************************************* */
/*<       real function cosdd(x) >*/
doublereal cosdd_(real *x)
{
    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    double acos(doublereal), cos(doublereal);

    /* Local variables */
    static real y, pi;

/*<       real x,y >*/
/*<       pi=acos(-1.) >*/
    pi = acos(-1.f);
/*<       y=x*pi/180. >*/
    y = *x * pi / 180.f;
/*<       cosdd = cos(y) >*/
    ret_val = cos(y);
/*<       return >*/
    return ret_val;
/*<       end >*/
} /* cosdd_ */

