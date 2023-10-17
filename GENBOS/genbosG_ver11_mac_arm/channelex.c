/* channelex.F -- translated by f2c (version 20230428).
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
    real wr[263], sigr[10520]	/* was [40][263] */;
} sigto_;

#define sigto_1 sigto_

struct {
    integer iprl;
} iprl_;

#define iprl_1 iprl_

/* Table of constant values */

static integer c__263 = 263;
static integer c__2 = 2;
static integer c__1 = 1;

/*<       SUBROUTINE CHANNELEX(E_N,E_DELTA,E_DEU,E_DIB,JCH) >*/
/* Subroutine */ int channelex_(real *e_n__, real *e_delta__, real *e_deu__, 
	real *e_dib__, integer *jch)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer ichannel, j, m, n;
    static real sigtotrif[2], xx[2];
    extern integer idihotomia_(real *, real *, integer *);
    static real ecm;
    extern doublereal ali_(real *, real *, real *, integer *);
    static integer nhe3;
    static real rnda;
    static integer ndeu;
    static real prob[41], vran[1];
    static integer ipos, npro[40], nsum;
    extern /* Subroutine */ int ranmar_(real *, integer *);
    static real sigrel[40];
    static integer nstart;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* FUNCTIONAL DESCRIPTION: */

/*    Select a channel to simulate according to the tabulated cross sections. */

/*    PROB(i)= probability for i-th selected channel */

/*    NCHAIN_P = number of selected reaction channel on proton */
/*    NCHAIN_N = number of selected reaction channel on neutron */
/*    NCHAIN_D = number of selected reaction channel on deuteron */

/* MODIFICATION DATE: */

/*     12.05.95       by Igor Pshenichnov */

/*  ************   by marco mirazita */

/*     23.04.99 */
/* --> Aggiunto l'azzeramento del vettore PROB(0:40), altrimenti non viene */
/*     calcolata correttamente la probabilita' dei vari canali di reazione. */
/*     Nella versione che girava su VMS questo probabilmente non era un */
/*     problema perche' il vettore era comunque azzerato dal compilatore. */
/*     29.04.99 */
/* --> Aggiunto il file di include PROCESS.INC con le variabili per il */
/*     conteggio dei canali di reazione selezionati */
/* --> Aggiunto il caso del canale 35 gamma d --> d rho0 (ed anche il canale */
/*     36 gamma d --> d eta, che non e' implementato completamente). */
/* --> Modificata la parte che interpola le sezioni d'urto totali per i */
/*     canali su d (era stata in precedenza modificata anche da federico). */
/*     05.05.99 */
/* --> Modificata tutta la parte di estrazione del canale di reazione, che */
/*     aveva dei giri strani. Adesso e' possibile selezionare un numero */
/*     qualsiasi di canali di reazione. */
/* --> Modificato il calcolo di NSUM (numero di canali totali) e NSTART */
/*     (indice del canale di partenza del conteggio), che ora sono calcolati */
/*     in base al numero di canali selezionati per ciascun tipo di bersaglio */
/*     (p, n, d, He3). */
/* --> Eliminata la variabile TARGFERMI. Adesso la selezione dei canali di */
/*     reazione e' fatta solo in base a TARGET. I canali su p (n) vengono */
/*     considerati solo se TARGET e' p (n), d oppure He3. */

/* COMMON BLOCKS: */

/*     /MONTECARL/ */
/*     /FLAG_COMP/ */
/*     /FLAG_PROC/ */
/*     /TARGET/ */

/* ----------------------------------------------------------------------- */
/*<       IMPLICIT NONE >*/
/*     include "process.inc" */
/*<       INCLUDE 'process.inc' >*/
/*<       INTEGER IBEAM,IELET,NCHAIN,ICHAIN(40),ICIBLE,IPOL,ILAM,IARM >*/
/*<       COMMON /FLAG_PROC/ NPROP,NPRON,NPROD,NPROHE3 >*/
/*< 	integer nprop_tot,npron_tot,nprod_tot,nprohe3_tot >*/
/*< 	common/proc_tot/nprop_tot,npron_tot,nprod_tot,nprohe3_tot >*/
/*< 	integer nchain_p,nchain_n,nchain_d,nchain_he3 >*/
/*< 	common/proc_sel/nchain_p,nchain_n,nchain_d,nchain_he3 >*/
/*<       REAL ETAG,FWHM >*/
/*<        >*/
/*<       REAL EGAM, P(4,11) >*/
/*<       INTEGER I(11),NPART >*/
/*<       COMMON/MONTECARL/ EGAM,I,NPART,P >*/
/*<       INTEGER ITARG,IIDH >*/
/*<       CHARACTER*12 TARGET,TARGFERMI >*/
/*<       COMMON /TARGET/ ITARG,IIDH,TARGET,TARGFERMI >*/
/*<       REAL WR(263),SIGR(40,263) >*/
/*<       COMMON/SIGTO/WR,SIGR >*/
/*<       INTEGER IPRL >*/
/*<       COMMON/IPRL/IPRL >*/
/*<       INTEGER NPRO(40),NDEU,NHE3,NSUM,NSTART,M,N,ICHANNEL,K,JCH,J >*/
/*<       INTEGER IPOS,IDIHOTOMIA >*/
/*<        >*/
/*<       REAL RNDM,ALI >*/
/*<       real vran(1) >*/
/* ***  Initialization */
/*<       DO J=1,40 >*/
    for (j = 1; j <= 40; ++j) {
/*<          NPRO(J) = 0 >*/
	npro[j - 1] = 0;
/*<       END DO >*/
    }
/*<       DO J=1,40 >*/
    for (j = 1; j <= 40; ++j) {
/*<          SIGREL(J) = 0. >*/
	sigrel[j - 1] = 0.f;
/*<       END DO >*/
    }
/*<       NDEU = 0 >*/
    ndeu = 0;
/*<       NHE3 = 0 >*/
    nhe3 = 0;
/*<       NSUM = 0 >*/
    nsum = 0;
/*<       NSTART = 0 >*/
    nstart = 0;
/* ***************************************** */
/* *** Initialization of the probabilities of reaction channels. */
/*<       do j=0,40 >*/
    for (j = 0; j <= 40; ++j) {
/*<          prob(j)=0. >*/
	prob[j] = 0.f;
/*<       enddo >*/
    }
/* ***************************************** */
/* ***   Different channels separation */
/* ************************************ */
/* *** Se non ci sono canali selezionati, li prendo tutti */
/*<       IF (NCHAIN.EQ.0) THEN >*/
    if (flag_comp__1.nchain == 0) {
/*<          DO J=1,40 >*/
	for (j = 1; j <= 40; ++j) {
/*<             ICHAIN(J) = NPRO(J) >*/
	    flag_comp__1.ichain[j - 1] = npro[j - 1];
/*<          END DO >*/
	}
/*<       END IF >*/
    }
/* ************************************ */
/* ************************************ */
/* *** Find energy bin for proton or neutron */
/* ************************************ */
/*<       IPOS = IDIHOTOMIA(E_N,WR,263) >*/
    ipos = idihotomia_(e_n__, sigto_1.wr, &c__263);
/* ************************************ */
/* *** Free or bound proton */
/* ******************************************** */
/*<        >*/
    if (s_cmp(target_1.target, "proton      ", (ftnlen)12, (ftnlen)12) == 0 ||
	     s_cmp(target_1.target, "deuteron    ", (ftnlen)12, (ftnlen)12) ==
	     0 || s_cmp(target_1.target, "He3      ", (ftnlen)12, (ftnlen)9) 
	    == 0) {
/*<          IF (NCHAIN.EQ.0) then >*/
	if (flag_comp__1.nchain == 0) {
/*<             NSUM = nsum + nprop_tot >*/
	    nsum += proc_tot__1.nprop_tot__;
/*<          else >*/
	} else {
/*<             NSUM = nsum + NCHAIN_p >*/
	    nsum += proc_sel__1.nchain_p__;
/*<          endif >*/
	}
/*<          DO J=nstart+1,nsum >*/
	i__1 = nsum;
	for (j = nstart + 1; j <= i__1; ++j) {
/*<             NPRO(J) = NPROP(J-nstart) >*/
	    npro[j - 1] = flag_proc__1.nprop[j - nstart - 1];
/*<          END DO >*/
	}
/*<          DO N=nstart+1,nsum >*/
	i__1 = nsum;
	for (n = nstart + 1; n <= i__1; ++n) {
/*<             ICHANNEL = NPRO(N) >*/
	    ichannel = npro[n - 1];
/* *** Linear interpolation in the energy bin */
/*<             SIGTOTRIF(1) = SIGR(ICHANNEL,IPOS-1) >*/
	    sigtotrif[0] = sigto_1.sigr[ichannel + (ipos - 1) * 40 - 41];
/*<             XX(1) = WR(IPOS-1) >*/
	    xx[0] = sigto_1.wr[ipos - 2];
/* *** IPOS deve essere < 264 */
/*<             if (ipos.le.263) then >*/
	    if (ipos <= 263) {
/*<                SIGTOTRIF(2) = SIGR(ICHANNEL,IPOS) >*/
		sigtotrif[1] = sigto_1.sigr[ichannel + ipos * 40 - 41];
/*<                XX(2) = WR(IPOS) >*/
		xx[1] = sigto_1.wr[ipos - 1];
/*<             else >*/
	    } else {
/*<                SIGTOTRIF(2) = SIGR(ICHANNEL,263) >*/
		sigtotrif[1] = sigto_1.sigr[ichannel + 10479];
/*<                XX(2) = WR(263) >*/
		xx[1] = sigto_1.wr[262];
/*<             endif >*/
	    }
/*            SIGTOTRIF(2) = SIGR(ICHANNEL,IPOS) */
/*            XX(2) = WR(IPOS) */
/*<             SIGREL(n) = ALI(E_N,XX,SIGTOTRIF,2) >*/
	    sigrel[n - 1] = ali_(e_n__, xx, sigtotrif, &c__2);
/*            write (15,*) 'N=',n,' xx=',xx */
/*            write (15,*) 'N=',n,' sigmarif=',sigtotrif */
/*            write (15,*) '  W=',e_n,' sigma=',sigrel(n) */
/*<          END DO >*/
	}
/*<          nstart=nstart+nchain_p >*/
	nstart += proc_sel__1.nchain_p__;
/*<       END IF >*/
    }
/* ************************************ */
/* *** Free or bound neutron */
/* ******************************************** */
/*<        >*/
    if (s_cmp(target_1.target, "neutron      ", (ftnlen)12, (ftnlen)13) == 0 
	    || s_cmp(target_1.target, "deuteron    ", (ftnlen)12, (ftnlen)12) 
	    == 0 || s_cmp(target_1.target, "He3      ", (ftnlen)12, (ftnlen)9)
	     == 0) {
/*<          IF (NCHAIN.EQ.0) then >*/
	if (flag_comp__1.nchain == 0) {
/*<             NSUM = nsum + npron_tot >*/
	    nsum += proc_tot__1.npron_tot__;
/*<          else >*/
	} else {
/*<             NSUM = nsum + NCHAIN_n >*/
	    nsum += proc_sel__1.nchain_n__;
/*<          endif >*/
	}
/*<          DO J=nstart+1,nsum >*/
	i__1 = nsum;
	for (j = nstart + 1; j <= i__1; ++j) {
/*<             NPRO(J) = NPRON(J-nstart) >*/
	    npro[j - 1] = flag_proc__1.npron[j - nstart - 1];
/*<          END DO >*/
	}
/*         write (15,*) 'NEUTRONE' */
/*<          DO N=nstart+1,nsum >*/
	i__1 = nsum;
	for (n = nstart + 1; n <= i__1; ++n) {
/*<             ICHANNEL = NPRO(N) >*/
	    ichannel = npro[n - 1];
/* *** Linear interpolation in the energy bin */
/*<             SIGTOTRIF(1) = SIGR(ICHANNEL,IPOS-1) >*/
	    sigtotrif[0] = sigto_1.sigr[ichannel + (ipos - 1) * 40 - 41];
/*<             XX(1) = WR(IPOS-1) >*/
	    xx[0] = sigto_1.wr[ipos - 2];
/* *** IPOS deve essere < 264 */
/*<             if (ipos.le.263) then >*/
	    if (ipos <= 263) {
/*<                SIGTOTRIF(2) = SIGR(ICHANNEL,IPOS) >*/
		sigtotrif[1] = sigto_1.sigr[ichannel + ipos * 40 - 41];
/*<                XX(2) = WR(IPOS) >*/
		xx[1] = sigto_1.wr[ipos - 1];
/*<             else >*/
	    } else {
/*<                SIGTOTRIF(2) = SIGR(ICHANNEL,263) >*/
		sigtotrif[1] = sigto_1.sigr[ichannel + 10479];
/*<                XX(2) = WR(263) >*/
		xx[1] = sigto_1.wr[262];
/*<             endif >*/
	    }
/*            SIGTOTRIF(2) = SIGR(ICHANNEL,IPOS) */
/*            XX(2) = WR(IPOS) */
/*<             SIGREL(n) = ALI(E_N,XX,SIGTOTRIF,2) >*/
	    sigrel[n - 1] = ali_(e_n__, xx, sigtotrif, &c__2);
/*            write (15,*) 'N=',n,' xx=',xx */
/*            write (15,*) 'N=',n,' sigmarif=',sigtotrif */
/*            write (15,*) '  W=',e_n,' sigma=',sigrel(n) */
/*<          END DO >*/
	}
/*<          nstart=nstart+nchain_n >*/
	nstart += proc_sel__1.nchain_n__;
/*<       END IF >*/
    }
/* ******************************************** */
/* *** Speciali deuteron channels */
/* ******************************************** */
/*<       IF (TARGET.EQ.'deuteron'.AND.NCHAIN_D.GT.0) THEN >*/
    if (s_cmp(target_1.target, "deuteron", (ftnlen)12, (ftnlen)8) == 0 && 
	    proc_sel__1.nchain_d__ > 0) {
/*<           IF (NCHAIN.EQ.0) then >*/
	if (flag_comp__1.nchain == 0) {
/*<              NSUM = nsum + nprod_tot >*/
	    nsum += proc_tot__1.nprod_tot__;
/*<           else >*/
	} else {
/*<              NSUM = nsum + NCHAIN_d >*/
	    nsum += proc_sel__1.nchain_d__;
/*<           endif >*/
	}
/* *** Aggiungo a NPRO i canali su deutone */
/*<          DO J=nstart+1,nsum >*/
	i__1 = nsum;
	for (j = nstart + 1; j <= i__1; ++j) {
/*<             NPRO(J) = NPROD(J-nstart) >*/
	    npro[j - 1] = flag_proc__1.nprod[j - nstart - 1];
/*<          END DO >*/
	}
/* *** Select energy bin */
/*<          IPOS = IDIHOTOMIA(e_deu,WR,263) ! FIND AN ENERGY INTERVAL  >*/
	ipos = idihotomia_(e_deu__, sigto_1.wr, &c__263);
/* TO WORK ON */
/* *** Select channel cross sections */
/*<          DO M=NSTART+1,Nsum >*/
	i__1 = nsum;
	for (m = nstart + 1; m <= i__1; ++m) {
/*<             ICHANNEL = NPRO(m) >*/
	    ichannel = npro[m - 1];
/* *** Linear interpolation in the energy bin */
/*<             SIGTOTRIF(1) = SIGR(ICHANNEL,IPOS-1) >*/
	    sigtotrif[0] = sigto_1.sigr[ichannel + (ipos - 1) * 40 - 41];
/*<             XX(1) = WR(IPOS-1) >*/
	    xx[0] = sigto_1.wr[ipos - 2];
/* *** IPOS deve essere < 264 */
/*<             if (ipos.le.263) then >*/
	    if (ipos <= 263) {
/*<                SIGTOTRIF(2) = SIGR(ICHANNEL,IPOS) >*/
		sigtotrif[1] = sigto_1.sigr[ichannel + ipos * 40 - 41];
/*<                XX(2) = WR(IPOS) >*/
		xx[1] = sigto_1.wr[ipos - 1];
/*<             else >*/
	    } else {
/*<                SIGTOTRIF(2) = SIGR(ICHANNEL,263) >*/
		sigtotrif[1] = sigto_1.sigr[ichannel + 10479];
/*<                XX(2) = WR(263) >*/
		xx[1] = sigto_1.wr[262];
/*<             endif >*/
	    }
/*<             SIGREL(M) = ALI(E_Deu,XX,SIGTOTRIF,2) >*/
	    sigrel[m - 1] = ali_(e_deu__, xx, sigtotrif, &c__2);
/*<          END DO >*/
	}
/*<          nstart=nstart+nchain_n >*/
	nstart += proc_sel__1.nchain_n__;
/*<       END IF >*/
    }
/* ******************************************** */
/* *** Special He3 channels - not implemented */
/* ******************************************** */
/*<       IF (TARGET.EQ.'He3'.AND.NCHAIN_HE3.GT.0) THEN >*/
    if (s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0 && 
	    proc_sel__1.nchain_he3__ > 0) {
/*<           IF (NCHAIN.EQ.0) then >*/
	if (flag_comp__1.nchain == 0) {
/*<              NSUM = nsum + nprohe3_tot >*/
	    nsum += proc_tot__1.nprohe3_tot__;
/*<           else >*/
	} else {
/*<              NSUM = nsum + NCHAIN_he3 >*/
	    nsum += proc_sel__1.nchain_he3__;
/*<           endif >*/
	}
/* *** Aggiungo a NPRO i canali du He3 */
/*<          DO J=nstart+1,nsum >*/
	i__1 = nsum;
	for (j = nstart + 1; j <= i__1; ++j) {
/*<             NPRO(J) = NPROHE3(J-nstart) >*/
	    npro[j - 1] = flag_proc__1.nprohe3[j - nstart - 1];
/*<          END DO >*/
	}
/*<          DO M=NSTART+1,NSum >*/
	i__1 = nsum;
	for (m = nstart + 1; m <= i__1; ++m) {
/*<  4          ICHANNEL = NPRO(m) >*/
/* L4: */
	    ichannel = npro[m - 1];
/* *** Select target energy */
/*<             IF (ICHANNEL.EQ.37) ECM = E_DELTA >*/
	    if (ichannel == 37) {
		ecm = *e_delta__;
	    }
/*<             IF (ICHANNEL.EQ.38.OR.ICHANNEL.EQ.39) ECM = E_DEU >*/
	    if (ichannel == 38 || ichannel == 39) {
		ecm = *e_deu__;
	    }
/*<             IF (ICHANNEL.EQ.40) ECM = E_DIB >*/
	    if (ichannel == 40) {
		ecm = *e_dib__;
	    }
/*<             IPOS = IDIHOTOMIA(Ecm,WR,263) ! FIND AN ENERGY INTERVAL  >*/
	    ipos = idihotomia_(&ecm, sigto_1.wr, &c__263);
/* TO WORK ON */
/*<             SIGTOTRIF(1) = SIGR(ICHANNEL,IPOS-1) >*/
	    sigtotrif[0] = sigto_1.sigr[ichannel + (ipos - 1) * 40 - 41];
/*<             XX(1) = WR(IPOS-1) >*/
	    xx[0] = sigto_1.wr[ipos - 2];
/*<             SIGTOTRIF(2) = SIGR(ICHANNEL,IPOS) >*/
	    sigtotrif[1] = sigto_1.sigr[ichannel + ipos * 40 - 41];
/*<             XX(2) = WR(IPOS) >*/
	    xx[1] = sigto_1.wr[ipos - 1];
/*<             SIGREL(M) = ALI(Ecm,XX,SIGTOTRIF,2) ! MAKE A LINEAR INTERPOL >*/
	    sigrel[m - 1] = ali_(&ecm, xx, sigtotrif, &c__2);
/* IN THIS INTERVAL ONLY */
/*<          END DO >*/
	}
/*<          nstart=nstart+nchain_he3 >*/
	nstart += proc_sel__1.nchain_he3__;
/*<       END IF >*/
    }
/*<       DO N=1,nsum >*/
    i__1 = nsum;
    for (n = 1; n <= i__1; ++n) {
/*<          PROB(N) = PROB(N-1)+SIGREL(N) >*/
	prob[n] = prob[n - 1] + sigrel[n - 1];
/*<       END DO >*/
    }
/*<       DO N=1,Nsum >*/
    i__1 = nsum;
    for (n = 1; n <= i__1; ++n) {
/*<          PROB(N) = PROB(N)/PROB(nsum) >*/
	prob[n] /= prob[nsum];
/*<       END DO >*/
    }
/*<       call ranmar(vran,1) >*/
    ranmar_(vran, &c__1);
/*<       rnda=vran(1) >*/
    rnda = vran[0];
/*<       DO J=1,NSum >*/
    i__1 = nsum;
    for (j = 1; j <= i__1; ++j) {
/*<          IF (RNDA.GE.PROB(J-1).AND.RNDA.LT.PROB(J)) THEN >*/
	if (rnda >= prob[j - 1] && rnda < prob[j]) {
/*<             JCH = npro(j) >*/
	    *jch = npro[j - 1];
/*<             ICHANNEL = JCH >*/
	    ichannel = *jch;
/*<             RETURN >*/
	    return 0;
/*<          END IF >*/
	}
/*<       END DO >*/
    }
/*<       RETURN >*/
    return 0;
/*<       END                       !*** END CHANNELEX >*/
} /* channelex_ */

/* RIA : moved from graph1.F */
/*<             INTEGER FUNCTION IDIHOTOMIA(TRY,ARRAY,NA) >*/
integer idihotomia_(real *try__, real *array, integer *na)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer nend, nbegin, nrange;


/*                FIND THE NUMBER ID SATISFIED: */
/*                ARRAY(ID-1).LT.TRY.LE.ARRAY(ID) */


/*<          REAL*4 ARRAY(NA) >*/
/*<          IF(TRY.LE.ARRAY(1)) THEN >*/
    /* Parameter adjustments */
    --array;

    /* Function Body */
    if (*try__ <= array[1]) {
/*<            IDIHOTOMIA=1 >*/
	ret_val = 1;
/*<          ELSE IF(TRY.GT.ARRAY(NA)) THEN >*/
    } else if (*try__ > array[*na]) {
/*<            IDIHOTOMIA=NA+1 >*/
	ret_val = *na + 1;
/*<          ELSE >*/
    } else {
/*<                              NBEGIN=1 >*/
	nbegin = 1;
/*<                              NEND=NA >*/
	nend = *na;
/*<                       NRANGE=(NEND-NBEGIN)/2+NBEGIN >*/
	nrange = (nend - nbegin) / 2 + nbegin;
/*<                  DO WHILE ((NEND-NBEGIN).GT.1) >*/
	while(nend - nbegin > 1) {
/*<                     IF(TRY.GT.ARRAY(NRANGE)) THEN >*/
	    if (*try__ > array[nrange]) {
/*<                         NBEGIN=NRANGE >*/
		nbegin = nrange;
/*<                     ELSE >*/
	    } else {
/*<                         NEND=NRANGE >*/
		nend = nrange;
/*<                     END IF >*/
	    }
/*<                     NRANGE=(NEND-NBEGIN)/2+NBEGIN >*/
	    nrange = (nend - nbegin) / 2 + nbegin;
/*<                  END DO >*/
	}
/*<          IDIHOTOMIA=NEND >*/
	ret_val = nend;
/*<          END IF >*/
    }
/*<          RETURN  >*/
    return ret_val;
/*<          END          >*/
} /* idihotomia_ */

/*<                  REAL*4 FUNCTION ALI(X,T,B,N) >*/
doublereal ali_(real *x, real *t, real *b, integer *n)
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
    static integer ik;
    extern integer idihotomia_(real *, real *, integer *);


/*          MAKE A LINEAR INTERPOLATION OF THE B(N) ARRAY */
/*                     AT X ABSCISSA VALUE */

/*<                  REAL*4 B(N),T(N) >*/
/*<         IK=IDIHOTOMIA(X,T,N) >*/
    /* Parameter adjustments */
    --b;
    --t;

    /* Function Body */
    ik = idihotomia_(x, &t[1], n);
/*<         IF (IK.EQ.1) THEN >*/
    if (ik == 1) {
/*<            ALI=B(1) >*/
	ret_val = b[1];
/*<            RETURN >*/
	return ret_val;
/*<         ELSE IF (IK.GT.N) THEN >*/
    } else if (ik > *n) {
/*<            ALI=B(N) >*/
	ret_val = b[*n];
/*<            RETURN >*/
	return ret_val;
/*<         ELSE  >*/
    } else {
/*<         ALI=((B(IK)-B(IK-1))/(T(IK)-T(IK-1)))*(X-T(IK-1))+B(IK-1) >*/
	ret_val = (b[ik] - b[ik - 1]) / (t[ik] - t[ik - 1]) * (*x - t[ik - 1])
		 + b[ik - 1];
/*<         END IF >*/
    }
/*<         RETURN >*/
    return ret_val;
/*<         END >*/
} /* ali_ */

