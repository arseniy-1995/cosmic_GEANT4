/* ev_nuc.F -- translated by f2c (version 20230428).
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
    char filen[20], filep[20];
} fname_;

#define fname_1 fname_

struct {
    real pn_n__[15]	/* was [5][3] */, beta_n__[4];
} nucleon_;

#define nucleon_1 nucleon_

struct {
    real pn_n_d__[15]	/* was [5][3] */, pn_d__[15]	/* was [5][3] */, 
	    beta_n_d__[4], beta_d__[4];
} h2_;

#define h2_1 h2_

struct {
    real pn_n_he3__[15]	/* was [5][3] */, pn_delta__[15]	/* was [5][3] 
	    */, pn_deu__[15]	/* was [5][3] */, pn_dib__[15]	/* was [5][3] 
	    */, beta_n_he3__[4], beta_delta__[4], beta_deu__[4], beta_dib__[4]
	    ;
} he3_;

#define he3_1 he3_

struct {
    real rmb1, rmb2, rmme, deltamass, dibmass, deumass, rhomass;
} masses_;

#define masses_1 masses_

struct {
    char partic[61440]	/* was [40][4][4][8] */;
    integer jdec[40], nch[40], nplev[120]	/* was [40][3] */;
    real br[120]	/* was [40][3] */;
} channels_;

#define channels_1 channels_

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
    integer nerr, loopv, ncaseg, loopc, nproc;
} flag_gen__;

#define flag_gen__1 flag_gen__

struct {
    real pn[15]	/* was [5][3] */;
} mom_nucleon__;

#define mom_nucleon__1 mom_nucleon__

struct {
    logical eg, cmsf;
} logeg_;

#define logeg_1 logeg_

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

static integer c__9 = 9;
static integer c__1 = 1;
static doublereal c_b25 = 2.;
static integer c__5 = 5;
static integer c__4 = 4;
static integer c__0 = 0;

/* *********************************************************************** */
/*<       SUBROUTINE EV_NUC >*/
/* Subroutine */ int ev_nuc__(void)
{
    /* Format strings */
    static char fmt_112[] = "(5x,\002JCH=\002,i5)";
    static char fmt_100[] = "(5x,\002CHANNEL NUMBER=\002,i4/5x,\002AMASS="
	    " \002,18f10.6)";
    static char fmt_101[] = "(5x,\002W= \002,f10.6,\002  TOTMASS= \002,f10"
	    ".6,\002  EGAM= \002,f10.6)";
    static char fmt_111[] = "(1x,\002 JD=\002,i6,\002 JCH=\002,i6,\002 JBR"
	    "=\002,i6,\002 NP=\002,i6)";

    /* System generated locals */
    integer i__1;
    real r__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    double pow_dd(doublereal *, doublereal *);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double exp(doublereal);

    /* Local variables */
    extern doublereal top_gamn__(real *, integer *);
    static real masstarg, argument, prob_top__[5];
    static integer j, k;
    static real q, s, t, u, w;
    extern /* Subroutine */ int channelex_(real *, real *, real *, real *, 
	    integer *), prim_kine__(void);
    static real p1[4], p2[4];
    static integer part_numb__[10], ia, ib, jd, jj, kk, jm, il;
    static real pl[4];
    extern /* Subroutine */ int gamn_charg__(integer *, integer *, integer *);
    static real ww;
    extern integer idihotomia_(real *, real *, integer *);
    static integer top_number__;
#define jch ((integer *)&reacch_1)
    static integer jbr;
    static real s5_9__, w_n__;
    static integer ith, lev;
    static real rmn, pcm1[90]	/* was [5][18] */, pcm2[90]	/* was [5][18]
	     */;
    extern integer code_(char *, ftnlen);
    static real beta[4];
    static integer ivar, jvar;
    static real pn_t__[4], vran[1];
#define plab1 ((real *)&kine_1__1 + 5)
    static real plab2[90]	/* was [5][18] */, w_dib__, delta[4], w_deu__;
    extern doublereal rmass_(char *, ftnlen);
    static real yrndm, amass1[18], amass2[18];
    extern /* Subroutine */ int genbod_(void), ranmar_(real *, integer *), 
	    gloren_(real *, real *, real *);
    static integer numpar;
    static real promul[9], t_coeff__;
    extern /* Subroutine */ int fermi_d__(real *, real *);
    static real w_delta__;
    extern /* Subroutine */ int gen_evt__(real *, integer *, integer *);
    static real totmass;
    static integer numplus;

    /* Fortran I/O blocks */
    static cilist io___29 = { 0, 15, 0, 0, 0 };
    static cilist io___30 = { 0, 15, 0, 0, 0 };
    static cilist io___31 = { 0, 15, 0, 0, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_112, 0 };
    static cilist io___40 = { 0, 6, 0, fmt_100, 0 };
    static cilist io___41 = { 0, 6, 0, fmt_101, 0 };
    static cilist io___54 = { 0, 15, 0, fmt_111, 0 };


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* FUNCTIONAL DESCRIPTION: */

/*    This subroutine generates the physical nuclear event, defines the */
/*    kinematics of selected process doing the proper LAB-CM transformation, */
/*    considering the decay channels and momentum distribution for processes */
/*    with bound nucleons. */

/* MODIFICATION DATE: */

/*     02.05.95       by Igor Pshenichnov */

/*   ************     by Marco mirazita */

/*     23.04.99 */
/* --> Aggiunto il file di include NT_KINE.INC con i dati delle particelle */
/*     dopo la prima fase di reazione (prima dei decadimenti). */
/* --> Aggiunta la chiamata alla routine PRIM_KINE che registra la */
/*     cinematica delle particelle prima dei decadimenti. */

/*     29.04.99 */
/* --> Aggiunto il file di include PROCESS.INC con le variabili per il */
/*     conteggio dei canali di reazione selezionati */

/*     05.05.99 */
/* --> Eliminato He3. Se TARGET = HE3 ==> errore e il programma si ferma. */
/* --> Azzerata la variabile TARGFERMI all'inizio e ricalcolata dopo */
/*     l'estrazione del canale di reazione, in base al canale estratto. */
/* --> Aggiunto il canale di reazione 35 gamma d --> rho0 d. Sono state */
/*     modificate tutti i punti in cui vengono calcolate le varie grandezze */
/*     cinematiche per il nuovo canale. (Aggiunte anche alcune chiamate per */
/*     il canale 36 gamma d --> eta d, ma e' eventualmente da finire). */

/*   ************ */

/* FORMAL PARAMETERS: */

/*    NP          = NUMBER OF OUTGOING PARTICLES (2-18) */
/*    ECM         = TOTAL CM ENERGY (REAL) */
/*    AMASS       = MASS OF I-TH PARTICLE (REAL) */
/*    KGENEV      = CONSTANT CROSS SECTION (=1), FERMI DISTRIBUTION (=2) */
/*    PCM(1,I)    = PX OF I-TH PARTICLE             (CM-SYSTEM) */
/*    PCM(2,I)    = PY OF I-TH PARTICLE             (CM-SYSTEM) */
/*    PCM(3,I)    = PZ OF I-TH PARTICLE             (CM-SYSTEM) */
/*    PCM(4,I)    = ENERGY OF I-TH PARTICLE         (CM-SYSTEM) */
/*    PCM(5,I)    = |P| OF I-TH PARTICLE            (CM-SYSTEM) */
/*    WT          = WEIGHT OF THE EVENT */
/*    PLAB(1,I)   = PX OF I-TH PARTICLE             (LAB-SYSTEM) */
/*    PLAB(2,I)   = PY OF I-TH PARTICLE             (LAB-SYSTEM) */
/*    PLAB(3,I)   = PZ OF I-TH PARTICLE             (LAB-SYSTEM) */
/*    PLAB(4,I)   = ENERGY OF I-TH PARTICLE         (LAB-SYSTEM) */
/*    PLAB(5,I)   = |P| OF I-TH PARTICLE            (LAB-SYSTEM) */
/*    KINETIC ENERGY = PLAB(4,I) - AMASS(I) */

/* COMMON BLOCKS: */

/*     /MONTECARL/ */
/*     /FLAG_COMP/ */
/*     /TARGET/ */
/*     /MASSES/ */
/*     /NUCLEUS/ */
/*     /NUCLEON/ */
/*     /H2/ */
/*     /HE3/ */
/*     /FLAG_PROC/ */
/*     /CHANNELS/ */
/*     /GENIN/ */
/*     /GENOUT/ */
/*     /MOM_NUCLEON/ */
/*     /LOGEG/ */
/* ----------------------------------------------------------------------- */
/*<       IMPLICIT NONE >*/
/*     include "process.inc" */
/*<       INCLUDE 'process.inc' >*/
/*<       INTEGER IBEAM,IELET,NCHAIN,ICHAIN(40),ICIBLE,IPOL,ILAM >*/
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
/*<       CHARACTER*20 FILEN,FILEP >*/
/*<       COMMON/FNAME/FILEN,FILEP >*/
/*<       REAL PN_N(5,3),BETA_N(4) >*/
/*<       COMMON /NUCLEON/ PN_N,BETA_N >*/
/*<       REAL PN_N_D(5,3),PN_D(5,3),BETA_N_D(4),BETA_D(4) >*/
/*<       COMMON /H2/ PN_N_D,PN_D,BETA_N_D,BETA_D >*/
/*<        >*/
/*<        >*/
/*<       REAL RMB1,RMB2,RMME,DELTAMASS,DIBMASS,DEUMASS,RHOMASS >*/
/*<        >*/
/*<       CHARACTER*12 PARTIC(40,0:3,0:3,8) >*/
/*<       INTEGER JDEC(40),NCH(40),NPLEV(40,3) >*/
/*<       REAL BR(40,3) >*/
/*<       COMMON /CHANNELS/ PARTIC,JDEC,NCH,NPLEV,BR >*/
/*<       INTEGER NP,KGENEV >*/
/*<       REAL ECM,AMASS,PCM,WT >*/
/*<       COMMON/GENIN/NP,ECM,AMASS(18),KGENEV >*/
/*<       COMMON/GENOUT/PCM(5,18),WT >*/
/*<       INTEGER NERR,LOOPV,NCASEG,LOOPC,NPROC,JCH >*/
/*<       COMMON /FLAG_GEN/ NERR,LOOPV,NCASEG,LOOPC,NPROC >*/
/*<       REAL RNDM,PROB_TOP(5),S5_9,TOP_GAMN >*/
/*<       INTEGER IPRINT,IVAR,JVAR,IDIHOTOMIA >*/
/*<       INTEGER*4 TOP_NUMBER, PART_NUMB(10) >*/
/*<        >*/
/*<        >*/
/*<       COMMON /MOM_NUCLEON/PN >*/
/*<       LOGICAL EG,CMSF >*/
/*<       COMMON/LOGEG/EG,CMSF >*/
/*<       REAL AMASS1(18),AMASS2(18), ARGUMENT >*/
/*<       real acosdd >*/
/* *** common for ntuple variables */
/*     #include "nt_kine.inc" */
/*<       INCLUDE 'nt_kine.inc' >*/
/*<       equivalence (plab1,pprim) >*/
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
/*<       equivalence (jch,j_channel) >*/
/*<       real vran(1),amx,mass_n_d,s >*/
/* ***   Inizialization */
/*<       KGENEV = 1 >*/
    genin_1.kgenev = 1;
/*<       DO J=1,11 >*/
    for (j = 1; j <= 11; ++j) {
/*<          DO JJ=1,4 >*/
	for (jj = 1; jj <= 4; ++jj) {
/*<             P(JJ,J) = 0. >*/
	    montecarl_1.p[jj + (j << 2) - 5] = 0.f;
/*<          END DO >*/
	}
/*<          I(J) = 0 >*/
	montecarl_1.i__[j - 1] = 0;
/*<       END DO >*/
    }
/*<       MASSTARG = 0. >*/
    masstarg = 0.f;
/*<       W_N = 0. >*/
    w_n__ = 0.f;
/*<       W_DEU = 0. >*/
    w_deu__ = 0.f;
/*<       W_DELTA = 0. >*/
    w_delta__ = 0.f;
/*<       W_DIB = 0. >*/
    w_dib__ = 0.f;
/*<       W = 0. >*/
    w = 0.f;
/*<       RMN = 0. >*/
    rmn = 0.f;
/*<       JCH = 0 >*/
    *jch = 0;
/*<       DO J=1,4 >*/
    for (j = 1; j <= 4; ++j) {
/*<           BETA(J) = 0. >*/
	beta[j - 1] = 0.f;
/*<           P1(J) = 0. >*/
	p1[j - 1] = 0.f;
/*<           P2(J) = 0. >*/
	p2[j - 1] = 0.f;
/*<           DELTA(J) = 0. >*/
	delta[j - 1] = 0.f;
/*<       END DO >*/
    }
/*<       DO J=1,5 >*/
    for (j = 1; j <= 5; ++j) {
/*<          DO JJ=1,3 >*/
	for (jj = 1; jj <= 3; ++jj) {
/*<             PN(J,JJ) = 0. >*/
	    mom_nucleon__1.pn[j + jj * 5 - 6] = 0.f;
/*<          END DO >*/
	}
/*<          DO JJ=1,18 >*/
	for (jj = 1; jj <= 18; ++jj) {
/*<             PCM1(J,JJ) = 0. >*/
	    pcm1[j + jj * 5 - 6] = 0.f;
/*<             PCM2(J,JJ) = 0. >*/
	    pcm2[j + jj * 5 - 6] = 0.f;
/*<             PLAB1(J,JJ) = 0. >*/
	    plab1[j + jj * 5 - 6] = 0.f;
/*<             PLAB2(J,JJ) = 0. >*/
	    plab2[j + jj * 5 - 6] = 0.f;
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       LEV = 0 >*/
    lev = 0;
/*<       JBR = 0 >*/
    jbr = 0;
/*<       TOTMASS = 0. >*/
    totmass = 0.f;
/*<       RMB1 = 0. >*/
    masses_1.rmb1 = 0.f;
/*<       RMB2 = 0. >*/
    masses_1.rmb2 = 0.f;
/*<       RMME = 0. >*/
    masses_1.rmme = 0.f;
/*<       WW = 0. >*/
    ww = 0.f;
/*<       NUMPLUS = 0 >*/
    numplus = 0;
/*<       NUMPAR = 0 >*/
    numpar = 0;
/*<       DO J=0,8 >*/
    for (j = 0; j <= 8; ++j) {
/*<          PROMUL(J) = 0. >*/
	promul[j] = 0.f;
/*<       END DO >*/
    }
/*<       ITH = 0 >*/
    ith = 0;
/*<       T_COEFF = 3. >*/
    t_coeff__ = 3.f;
/* ********************************** */
/* *** Azzero anche TARGFERMI */
/* ********************************** */
/*<       targfermi=' ' >*/
    s_copy(target_1.targfermi, " ", (ftnlen)12, (ftnlen)1);
/* *** Impulso di fermi */
/*<       p_fermi=0. >*/
    fermi_motion__1.p_fermi__ = 0.f;
/* ***  If the energy of gamma is below the pion threshold than keep */
/* ***  gamma in output array and nothing more (Target nucleon at rest) */
/*<       IF(EGAM.LE.0.151) THEN >*/
    if (montecarl_1.egam <= .151f) {
/*<          NPART  = 1 >*/
	montecarl_1.npart = 1;
/*<          P(1,1) = 0. >*/
	montecarl_1.p[0] = 0.f;
/*<          P(2,1) = 0. >*/
	montecarl_1.p[1] = 0.f;
/*<          P(3,1) = EGAM >*/
	montecarl_1.p[2] = montecarl_1.egam;
/*<          I(1)   = 1 >*/
	montecarl_1.i__[0] = 1;
/*<          RETURN >*/
	return 0;
/*<       ENDIF >*/
    }
/* ***   Definition of masses of decaying */
/* ***   particles (Delta, Rho, Omega) */
/*<  11   DELTAMASS = RMASS('delta++     ') >*/
L11:
    masses_1.deltamass = rmass_("delta++     ", (ftnlen)12);
/*<       RHOMASS =   RMASS('rho0        ') >*/
    masses_1.rhomass = rmass_("rho0        ", (ftnlen)12);
/*<       DIBMASS =   RMASS('dibaryon    ') >*/
    masses_1.dibmass = rmass_("dibaryon    ", (ftnlen)12);
/*<       DEUMASS =   RMASS('deuteron    ') >*/
    masses_1.deumass = rmass_("deuteron    ", (ftnlen)12);
/* ************************************************************* */
/* ***   Computing CM energy and beta, considering, if necessary, */
/* ***   the Fermi motion */
/* ************************************************************* */
/*<       IF (TARGET.EQ.'deuteron    ') THEN >*/
    if (s_cmp(target_1.target, "deuteron    ", (ftnlen)12, (ftnlen)12) == 0) {
/* *** Calcolo moto di fermi per i nucleoni */
/*<          if (nchain_p+nchain_n.ne.0) then >*/
	if (proc_sel__1.nchain_p__ + proc_sel__1.nchain_n__ != 0) {
/*<             CALL FERMI_D(W_N,W_DEU) >*/
	    fermi_d__(&w_n__, &w_deu__);
/*<          endif >*/
	}
/* *** Kinematical quantities for deuteron */
/*<          S = DEUMASS**2+2.*DEUMASS*EGAM >*/
/* Computing 2nd power */
	r__1 = masses_1.deumass;
	s = r__1 * r__1 + masses_1.deumass * 2.f * montecarl_1.egam;
/*<          W_DEU = SQRT(S) >*/
	w_deu__ = sqrt(s);
/*<          PN_D(1,1) = 0. >*/
	h2_1.pn_d__[0] = 0.f;
/*<          PN_D(2,1) = 0. >*/
	h2_1.pn_d__[1] = 0.f;
/*<          PN_D(3,1) = 0. >*/
	h2_1.pn_d__[2] = 0.f;
/*<          PN_D(4,1) = DEUMASS >*/
	h2_1.pn_d__[3] = masses_1.deumass;
/* *** Beta(CM) for (gamma+deuteron) system */
/*<          BETA_D(1) = 0. >*/
	h2_1.beta_d__[0] = 0.f;
/*<          BETA_D(2) = 0. >*/
	h2_1.beta_d__[1] = 0.f;
/*<          BETA_D(3) = EGAM/(EGAM+DEUMASS) >*/
	h2_1.beta_d__[2] = montecarl_1.egam / (montecarl_1.egam + 
		masses_1.deumass);
/*<          BETA_D(4) = 0. >*/
	h2_1.beta_d__[3] = 0.f;
/*<       ELSE IF (TARGET.EQ.'He3         ') THEN >*/
    } else if (s_cmp(target_1.target, "He3         ", (ftnlen)12, (ftnlen)12) 
	    == 0) {
/*<          WRITE(15,*) '***********************************' >*/
	s_wsle(&io___29);
	do_lio(&c__9, &c__1, "***********************************", (ftnlen)
		35);
	e_wsle();
/*<          WRITE(15,*) ' ERROR: Target He3 not implemented' >*/
	s_wsle(&io___30);
	do_lio(&c__9, &c__1, " ERROR: Target He3 not implemented", (ftnlen)34)
		;
	e_wsle();
/*<          WRITE(15,*) '***********************************' >*/
	s_wsle(&io___31);
	do_lio(&c__9, &c__1, "***********************************", (ftnlen)
		35);
	e_wsle();
/*<          stop >*/
	s_stop("", (ftnlen)0);
/*         CALL FERMI_HE3(W_N,W_DELTA,W_DEU,W_DIB) */
/*<       ELSE IF(TARGET.EQ.'proton      '.OR.TARGET.EQ.'neutron     ')THEN >*/
    } else if (s_cmp(target_1.target, "proton      ", (ftnlen)12, (ftnlen)12) 
	    == 0 || s_cmp(target_1.target, "neutron     ", (ftnlen)12, (
	    ftnlen)12) == 0) {
/* **** */
/*         RMN = MASSTARG */
/*<          RMN = rmass(target) >*/
	rmn = rmass_(target_1.target, (ftnlen)12);
/* **** */
/*<          PN_N(4,1) = RMN >*/
	nucleon_1.pn_n__[3] = rmn;
/*<          W_N = SQRT(RMN*RMN + 2.*RMN*EGAM) >*/
	w_n__ = sqrt(rmn * rmn + rmn * 2.f * montecarl_1.egam);
/*<          BETA_N(1) = 0. >*/
	nucleon_1.beta_n__[0] = 0.f;
/*<          BETA_N(2) = 0. >*/
	nucleon_1.beta_n__[1] = 0.f;
/*<          BETA_N(3) = EGAM/(EGAM+RMN) >*/
	nucleon_1.beta_n__[2] = montecarl_1.egam / (montecarl_1.egam + rmn);
/*<          BETA_N(4) = 0. >*/
	nucleon_1.beta_n__[3] = 0.f;
/*<       END IF >*/
    }
/* *** Computing beta(CM) of (gamma+target) system */
/*<       if (TARGET.EQ.'proton      '.OR.TARGET.EQ.'neutron     ') THEN >*/
    if (s_cmp(target_1.target, "proton      ", (ftnlen)12, (ftnlen)12) == 0 ||
	     s_cmp(target_1.target, "neutron     ", (ftnlen)12, (ftnlen)12) ==
	     0) {
/*<          do j=1,3 >*/
	for (j = 1; j <= 3; ++j) {
/*<             beta_cm(j)=beta_n(j) >*/
	    betacm_1.beta_cm__[j - 1] = nucleon_1.beta_n__[j - 1];
/*<          enddo >*/
	}
/*<       else if (TARGET.EQ.'deuteron    ') THEN >*/
    } else if (s_cmp(target_1.target, "deuteron    ", (ftnlen)12, (ftnlen)12) 
	    == 0) {
/*<          do j=1,3 >*/
	for (j = 1; j <= 3; ++j) {
/*<             beta_cm(j)=beta_d(j) >*/
	    betacm_1.beta_cm__[j - 1] = h2_1.beta_d__[j - 1];
/*<          enddo >*/
	}
/*<       endif >*/
    }
/*<        >*/
    d__1 = (doublereal) betacm_1.beta_cm__[1];
    d__2 = (doublereal) betacm_1.beta_cm__[2];
    betacm_1.beta_cm__[3] = 1.f / sqrt(1.f - (betacm_1.beta_cm__[0] * 2.f + 
	    pow_dd(&d__1, &c_b25) + pow_dd(&d__2, &c_b25)));
/* ************************************************************* */
/* ***   Selecting reaction channel */
/* ************************************************************* */
/*<       IF (NCHAIN.EQ.1) THEN >*/
    if (flag_comp__1.nchain == 1) {
/*<          JCH = ICHAIN(1) >*/
	*jch = flag_comp__1.ichain[0];
/*<       ELSE >*/
    } else {
/*<          CALL CHANNELEX(W_N,W_DELTA,W_DEU,W_DIB,JCH) >*/
	channelex_(&w_n__, &w_delta__, &w_deu__, &w_dib__, jch);
/*<       END IF >*/
    }
/* ********************************************************************** */
/* *** L'impulso di Fermi e' calcolato in ogni caso. */
/* *** Se la reazione e' su d l'impulso di Fermi deve essere nullo */
/* ********************************************************************** */
/*<       if (jch.ge.33) p_fermi=0. >*/
    if (*jch >= 33) {
	fermi_motion__1.p_fermi__ = 0.f;
    }
/* ***************************************************** */
/* *** Calcolo TARGFERMI in base al canale estratto */
/* *** TARGFERMI is the mass of the interacting particle */
/* ***************************************************** */
/* *** Protone */
/*<        >*/
    if (*jch == 1 || *jch == 2 || *jch >= 5 && *jch <= 7 || *jch == 11 || *
	    jch == 12 || *jch == 15 || *jch == 16 || *jch == 19 || *jch == 21 
	    || *jch == 23 || *jch == 24 || *jch == 27 || *jch >= 29 && *jch <=
	     32) {
/* *** Massa del bersaglio (incluso eventualmente il moto di Fermi) */
/*<          if (target.eq.'deuteron') then >*/
	if (s_cmp(target_1.target, "deuteron", (ftnlen)12, (ftnlen)8) == 0) {
/*<             masstarg=rm(1) >*/
	    masstarg = fermi_motion__1.rm[0];
/*<          else >*/
	} else {
/*<             MASSTARG = RMASS('proton      ') >*/
	    masstarg = rmass_("proton      ", (ftnlen)12);
/*<          endif >*/
	}
/*<          TARGFERMI = 'proton      ' >*/
	s_copy(target_1.targfermi, "proton      ", (ftnlen)12, (ftnlen)12);
/* *** Neutrone */
/*<        >*/
    } else if (*jch == 3 || *jch == 4 || *jch >= 8 && *jch <= 10 || *jch == 
	    13 || *jch == 14 || *jch == 17 || *jch == 18 || *jch == 20 || *
	    jch == 22 || *jch == 25 || *jch == 26 || *jch == 28) {
/* *** Massa del bersaglio (incluso eventualmente il moto di Fermi) */
/*<          if (target.eq.'deuteron') then >*/
	if (s_cmp(target_1.target, "deuteron", (ftnlen)12, (ftnlen)8) == 0) {
/*<             masstarg=rm(1) >*/
	    masstarg = fermi_motion__1.rm[0];
/*<          else >*/
	} else {
/*<             MASSTARG = RMASS('neutron     ') >*/
	    masstarg = rmass_("neutron     ", (ftnlen)12);
/*<          endif >*/
	}
/*<          TARGFERMI = 'neutron     ' >*/
	s_copy(target_1.targfermi, "neutron     ", (ftnlen)12, (ftnlen)12);
/* *** Deuterio */
/*<        >*/
    } else if (*jch >= 33 && *jch <= 36) {
/*<          MASSTARG = RMASS('deuteron    ') >*/
	masstarg = rmass_("deuteron    ", (ftnlen)12);
/*<          TARGFERMI = 'deuteron    ' >*/
	s_copy(target_1.targfermi, "deuteron    ", (ftnlen)12, (ftnlen)12);
/* *** He3 */
/*<        >*/
    } else if (*jch >= 37 && *jch <= 40) {
/*<          TARGFERMI = 'He3 ' >*/
	s_copy(target_1.targfermi, "He3 ", (ftnlen)12, (ftnlen)4);
/*<       endif >*/
    }
/* *********************************************** */
/* *** fine modifica */
/* *********************************************** */
/* *********************************************** */
/* ***    BETA calculation */
/* *** For bound nucleons (TARGET=deuteron) */
/* *** photon+interacting nucleon must be considered */
/* *********************************************** */
/*<  12   IF (JCH.LE.32) THEN >*/
/* L12: */
    if (*jch <= 32) {
/*<          W = W_N >*/
	w = w_n__;
/*<          DO J=1,4 >*/
	for (j = 1; j <= 4; ++j) {
/*<             IF (TARGET.EQ.'deuteron') BETA(J) = -BETA_N_D(J) >*/
	    if (s_cmp(target_1.target, "deuteron", (ftnlen)12, (ftnlen)8) == 
		    0) {
		beta[j - 1] = -h2_1.beta_n_d__[j - 1];
	    }
/*            IF (TARGET.EQ.'He3') BETA(J) = -BETA_N_HE3(J) */
/*<        >*/
	    if (s_cmp(target_1.target, "proton", (ftnlen)12, (ftnlen)6) == 0 
		    || s_cmp(target_1.target, "neutron", (ftnlen)12, (ftnlen)
		    7) == 0) {
		beta[j - 1] = -nucleon_1.beta_n__[j - 1];
	    }
/*<          END DO >*/
	}
/*<          DO J=1,5 >*/
	for (j = 1; j <= 5; ++j) {
/*<             DO JJ=1,3 >*/
	    for (jj = 1; jj <= 3; ++jj) {
/*<                IF (TARGET.EQ.'deuteron') PN(J,JJ)=PN_N_D(J,JJ) >*/
		if (s_cmp(target_1.target, "deuteron", (ftnlen)12, (ftnlen)8) 
			== 0) {
		    mom_nucleon__1.pn[j + jj * 5 - 6] = h2_1.pn_n_d__[j + jj *
			     5 - 6];
		}
/*<                IF (TARGET.EQ.'He3') PN(J,JJ)=PN_N_HE3(J,JJ) >*/
		if (s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0)
			 {
		    mom_nucleon__1.pn[j + jj * 5 - 6] = he3_1.pn_n_he3__[j + 
			    jj * 5 - 6];
		}
/*<        >*/
		if (s_cmp(target_1.target, "proton", (ftnlen)12, (ftnlen)6) ==
			 0 || s_cmp(target_1.target, "neutron", (ftnlen)12, (
			ftnlen)7) == 0) {
		    mom_nucleon__1.pn[j + jj * 5 - 6] = nucleon_1.pn_n__[j + 
			    jj * 5 - 6];
		}
/*<             END DO >*/
	    }
/*<          END DO >*/
	}
/*<       END IF >*/
    }
/* *************************************************************** */
/* ***    BETA calculation */
/* *** For diffractive reactions on deuteron, */
/* *** photon+deuteron must be considered */
/* *************************************************************** */
/*<       IF ( (JCH.gE.33).and.(jch.le.36) ) THEN >*/
    if (*jch >= 33 && *jch <= 36) {
/*<          W = W_DEU >*/
	w = w_deu__;
/*<          DO J=1,4 >*/
	for (j = 1; j <= 4; ++j) {
/*<             BETA(J) = -BETA_D(J) >*/
	    beta[j - 1] = -h2_1.beta_d__[j - 1];
/*<          END DO >*/
	}
/*<          DO J=1,5 >*/
	for (j = 1; j <= 5; ++j) {
/*<             DO JJ=1,3 >*/
	    for (jj = 1; jj <= 3; ++jj) {
/*<                PN(J,JJ)=PN_D(J,JJ) >*/
		mom_nucleon__1.pn[j + jj * 5 - 6] = h2_1.pn_d__[j + jj * 5 - 
			6];
/*<             END DO >*/
	    }
/*<          END DO >*/
	}
/*<       END IF >*/
    }
/* *************************************************************** */
/* ***     CALCOLO DI BETA */
/* *** Le reazioni su He3 sono da implementare */
/* *************************************************************** */
/*<       IF (JCH.EQ.37) THEN >*/
    if (*jch == 37) {
/*<          W = W_DELTA >*/
	w = w_delta__;
/*<          DO J=1,4 >*/
	for (j = 1; j <= 4; ++j) {
/*<             BETA(J) = -BETA_DELTA(J) >*/
	    beta[j - 1] = -he3_1.beta_delta__[j - 1];
/*<          END DO >*/
	}
/*<          DO J=1,5 >*/
	for (j = 1; j <= 5; ++j) {
/*<             DO JJ=1,3 >*/
	    for (jj = 1; jj <= 3; ++jj) {
/*<                PN(J,JJ)=PN_DELTA(J,JJ) >*/
		mom_nucleon__1.pn[j + jj * 5 - 6] = he3_1.pn_delta__[j + jj * 
			5 - 6];
/*<             END DO >*/
	    }
/*<          END DO >*/
	}
/*<       END IF >*/
    }
/*<       IF (JCH.EQ.38.OR.JCH.EQ.39) THEN >*/
    if (*jch == 38 || *jch == 39) {
/*<          W = W_DEU >*/
	w = w_deu__;
/*<          DO J=1,4 >*/
	for (j = 1; j <= 4; ++j) {
/*<             BETA(J) = -BETA_DEU(J) >*/
	    beta[j - 1] = -he3_1.beta_deu__[j - 1];
/*<          END DO >*/
	}
/*<          DO J=1,5 >*/
	for (j = 1; j <= 5; ++j) {
/*<             DO JJ=1,2 >*/
	    for (jj = 1; jj <= 2; ++jj) {
/*<                PN(J,JJ)=PN_DEU(J,JJ) >*/
		mom_nucleon__1.pn[j + jj * 5 - 6] = he3_1.pn_deu__[j + jj * 5 
			- 6];
/*<             END DO >*/
	    }
/*<          END DO >*/
	}
/*<       END IF >*/
    }
/*<       IF (JCH.EQ.40) THEN >*/
    if (*jch == 40) {
/*<          W = W_DIB >*/
	w = w_dib__;
/*<          DO J=1,4 >*/
	for (j = 1; j <= 4; ++j) {
/*<             BETA(J) = -BETA_DIB(J) >*/
	    beta[j - 1] = -he3_1.beta_dib__[j - 1];
/*<          END DO >*/
	}
/*<          DO J=1,5 >*/
	for (j = 1; j <= 5; ++j) {
/*<             DO JJ=1,3 >*/
	    for (jj = 1; jj <= 3; ++jj) {
/*<                PN(J,JJ)=PN_DIB(J,JJ) >*/
		mom_nucleon__1.pn[j + jj * 5 - 6] = he3_1.pn_dib__[j + jj * 5 
			- 6];
/*<             END DO >*/
	    }
/*<          END DO >*/
	}
/*<       END IF >*/
    }
/* *********************************************** */
/* *********************************************** */
/* *** KINEMATICS OF FINAL STATE */
/* *********************************************** */
/* *********************************************** */
/* ************************************* */
/* ***   First level */
/* ************************************* */
/*<       LEV = 1 >*/
    lev = 1;
/*<       JBR = 0 >*/
    jbr = 0;
/*<       IF((JCH.LE.0).OR.(JCH.GT.40)) THEN >*/
    if (*jch <= 0 || *jch > 40) {
/*<          WRITE(0,112) JCH >*/
	s_wsfe(&io___32);
	do_fio(&c__1, (char *)&(*jch), (ftnlen)sizeof(integer));
	e_wsfe();
/*<  112     FORMAT(5X,'JCH=',I5) >*/
/*<       ENDIF  >*/
    }
/* ************************************************************* */
/* *** Check on final masses */
/* ************************************************************* */
/*<       TOTMASS = 0. >*/
    totmass = 0.f;
/* *** Multipion channels */
/* *** Treat these channels as a special case, select mutiplicity here */
/*<       IF(JCH.EQ.27.OR.JCH.EQ.28) THEN >*/
    if (*jch == 27 || *jch == 28) {
/*<          S5_9=0.0 >*/
	s5_9__ = 0.f;
/*<          DO IVAR=1,5 >*/
	for (ivar = 1; ivar <= 5; ++ivar) {
/*<             S5_9=S5_9+TOP_GAMN(EGAM,IVAR+4) >*/
	    i__1 = ivar + 4;
	    s5_9__ += top_gamn__(&montecarl_1.egam, &i__1);
/*<             PROB_TOP(IVAR)=S5_9 >*/
	    prob_top__[ivar - 1] = s5_9__;
/*<          ENDDO >*/
	}
/*<          DO IVAR=1,5 >*/
	for (ivar = 1; ivar <= 5; ++ivar) {
/*<             PROB_TOP(IVAR)=PROB_TOP(IVAR)/PROB_TOP(5) >*/
	    prob_top__[ivar - 1] /= prob_top__[4];
/*<          ENDDO >*/
	}
/*<          call ranmar(vran,1) >*/
	ranmar_(vran, &c__1);
/*<          TOP_NUMBER=IDIHOTOMIA(vran(1),PROB_TOP,5)+4 >*/
	top_number__ = idihotomia_(vran, prob_top__, &c__5) + 4;
/*<          NP=TOP_NUMBER >*/
	genin_1.np = top_number__;
/*<          DO JM=1,NP-1           ! put pion masses >*/
	i__1 = genin_1.np - 1;
	for (jm = 1; jm <= i__1; ++jm) {
/*<             AMASS(JM) = RMASS(PARTIC(27,1,0,1)) >*/
	    genin_1.amass[jm - 1] = rmass_(channels_1.partic + 792, (ftnlen)
		    12);
/*<             AMASS1(JM) = AMASS(JM) >*/
	    amass1[jm - 1] = genin_1.amass[jm - 1];
/*<             TOTMASS = TOTMASS+AMASS(JM) >*/
	    totmass += genin_1.amass[jm - 1];
/*<          END DO >*/
	}
/*<          AMASS(NP)=RMASS(PARTIC(27,1,0,5)) ! put nucleon mass  >*/
	genin_1.amass[genin_1.np - 1] = rmass_(channels_1.partic + 31512, (
		ftnlen)12);
/*<          AMASS1(NP) = AMASS(NP) >*/
	amass1[genin_1.np - 1] = genin_1.amass[genin_1.np - 1];
/*<          TOTMASS = TOTMASS+AMASS(NP) >*/
	totmass += genin_1.amass[genin_1.np - 1];
/* ***  For other channels ... */
/*<       ELSE >*/
    } else {
/*<          NP = NPLEV(JCH,1) >*/
	genin_1.np = channels_1.nplev[*jch - 1];
/*<          DO JM=1,NP >*/
	i__1 = genin_1.np;
	for (jm = 1; jm <= i__1; ++jm) {
/*<             AMASS(JM) = RMASS(PARTIC(JCH,LEV,JBR,JM)) >*/
	    genin_1.amass[jm - 1] = rmass_(channels_1.partic + (*jch + (lev + 
		    (jbr + (jm << 2) << 2)) * 40 - 641) * 12, (ftnlen)12);
/*<             AMASS1(JM) = AMASS(JM) >*/
	    amass1[jm - 1] = genin_1.amass[jm - 1];
/*<             TOTMASS = TOTMASS+AMASS(JM) >*/
	    totmass += genin_1.amass[jm - 1];
/*<          END DO >*/
	}
/*<       ENDIF >*/
    }
/*<       IF (W.LT.TOTMASS+0.000001) THEN >*/
    if (w < totmass + 1e-6f) {
/*<          ITH = ITH+1 >*/
	++ith;
/*<          IF ((MOD(ITH,100000).EQ.0).AND.(ITH.GT.100000)) THEN >*/
	if (ith % 100000 == 0 && ith > 100000) {
/*         WRITE (6,*) 'REACTION IMPOSSIBLE ' */
/*         WRITE (6,*) 'PHOTON ENERGY UNDER THRESHOLD' */
/*<             WRITE (6,100) JCH,(AMASS(JJ),JJ=1,NP) >*/
	    s_wsfe(&io___40);
	    do_fio(&c__1, (char *)&(*jch), (ftnlen)sizeof(integer));
	    i__1 = genin_1.np;
	    for (jj = 1; jj <= i__1; ++jj) {
		do_fio(&c__1, (char *)&genin_1.amass[jj - 1], (ftnlen)sizeof(
			real));
	    }
	    e_wsfe();
/*<             WRITE (6,101) W,TOTMASS,EGAM >*/
	    s_wsfe(&io___41);
	    do_fio(&c__1, (char *)&w, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&totmass, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&montecarl_1.egam, (ftnlen)sizeof(real));
	    e_wsfe();
/*<  100        FORMAT(5X,'CHANNEL NUMBER=',I4/5X,'AMASS= ',18F10.6) >*/
/*<  101        FORMAT(5X,'W= ',F10.6,'  TOTMASS= ',F10.6,'  EGAM= ',F10.6) >*/
/*<          END IF >*/
	}
/*<          GOTO 11 >*/
	goto L11;
/*<       END IF >*/
    }
/* ****************************************************************************** */
/* ***   Choosing masses for target baryon, scattered baryon, spectator */
/* ***   Baryon, if any, and of scattered mesons */
/* ****************************************************************************** */
/* *** */
/* *** channels  1 - 14 e 19 - 22, 29, 30, 32, 33, 34, 35, 36 */
/* *** */
/*<        >*/
    if (*jch >= 1 && *jch <= 14 || *jch == 19 || *jch == 20 || *jch == 29 || *
	    jch == 30 || *jch >= 32 && *jch <= 36) {
/* **************************** */
/* *** Target */
/*<          RMB1 = masstarg >*/
	masses_1.rmb1 = masstarg;
/* **************************** */
/* *** Scattered baryon (nucleon) */
/*<          RMB2 = AMASS(2) >*/
	masses_1.rmb2 = genin_1.amass[1];
/* *** Scattered meson(pion) */
/*<          RMME = AMASS(1) >*/
	masses_1.rmme = genin_1.amass[0];
/* *** */
/* *** If scattered Baryon is delta */
/* *** */
/*<          IF (JCH.GE.5.AND.JCH.LE.10) THEN >*/
	if (*jch >= 5 && *jch <= 10) {
/*<             WW = DELTAMASS >*/
	    ww = masses_1.deltamass;
/*<             RMB2 = WW >*/
	    masses_1.rmb2 = ww;
/*<             IF (W.LT.(RMME+RMB2+0.000001)) GOTO 11 >*/
	    if (w < masses_1.rmme + masses_1.rmb2 + 1e-6f) {
		goto L11;
	    }
/*<          END IF >*/
	}
/* *** */
/* *** If scattered meson is rho */
/* *** */
/*<          IF ((JCH.GE.11.AND.JCH.LE.14).or.(jch.eq.35)) THEN >*/
	if (*jch >= 11 && *jch <= 14 || *jch == 35) {
/*<             WW = RHOMASS >*/
	    ww = masses_1.rhomass;
/*<             RMME = WW >*/
	    masses_1.rmme = ww;
/*<             IF (W.LT.(RMME+RMB2+0.000001)) GOTO 11 >*/
	    if (w < masses_1.rmme + masses_1.rmb2 + 1e-6f) {
		goto L11;
	    }
/*<          END IF >*/
	}
/*<       END IF >*/
    }
/* *** */
/* ***  If scattered meson is omega */
/* *** */
/*<       IF (JCH.EQ.21.OR.JCH.EQ.22) THEN >*/
    if (*jch == 21 || *jch == 22) {
/*<          RMB2 = AMASS(2) >*/
	masses_1.rmb2 = genin_1.amass[1];
/*<          WW = AMASS(1) >*/
	ww = genin_1.amass[0];
/*<          RMME = WW >*/
	masses_1.rmme = ww;
/*<          IF (W.LT.(RMME+RMB2+0.000001)) GOTO 11 >*/
	if (w < masses_1.rmme + masses_1.rmb2 + 1e-6f) {
	    goto L11;
	}
/*<       END IF >*/
    }
/* *** */
/* *** He3 reaction */
/* *** */
/*<       IF (JCH.EQ.40) WW=AMASS(2) >*/
    if (*jch == 40) {
	ww = genin_1.amass[1];
    }
/* ************************************************** */
/* *** Prova stampa */
/* ************************************************** */
/*      if (TARGFERMI.eq.'deuteron    ') then */
/*         amx=0. */
/*         write (*,35)amx,amx,amx,RMASS(targfermi),RMASS(targfermi) */
/*      else */
/*         amx=sqrt(pn_n_d(4,1)**2.- */
/*     +        (pn_n_d(1,1)**2.+pn_n_d(2,1)**2.+pn_n_d(3,1)**2.)) */
/*         write (*,37) (pn_n_d(j,1),j=1,4),amx */
/*         amx=sqrt(pn_n_d(4,2)**2.- */
/*     +        (pn_n_d(1,2)**2.+pn_n_d(2,2)**2.+pn_n_d(3,2)**2.)) */
/*         write (*,38) (pn_n_d(j,2),j=1,4),amx */
/*      endif */
/* 35   format (1x,'deut. ',5(1x,f8.5)) */
/* 37   format (1x,'N.int.',5(1x,f8.5)) */
/* 38   format (1x,'N. sp.',5(1x,f8.5)) */
/* ********************************************************************** */
/* *** Generating final state */
/* ********************************************************************** */
/*<  20    >*/
L20:
    if (*jch >= 15 && *jch <= 18 || *jch >= 23 && *jch <= 28 || *jch == 31 || 
	    *jch >= 37) {
/* *** */
/* *** Channels with 2 pions or more */
/* *** */
/*<          ECM = W >*/
	genin_1.ecm = w;
/* !    PRINT 113,JCH,W,(AMASS(IPRINT),IPRINT=1,NP) */
/*<  113     FORMAT (1X,'JCH=',I4,' W=',F9.5,' AMASS=',18F9.5) >*/
/* L113: */
/*<          CALL GENBOD >*/
	genbod_();
/*<       ELSE >*/
    } else {
/* *** */
/* *** Other channels */
/* *** */
/*<          CALL GEN_EVT(W,LEV,JCH) >*/
	gen_evt__(&w, &lev, jch);
/*<       END IF >*/
    }
/* ********************************************************************** */
/* *** CM final momenta */
/* ********************************************************************** */
/*<       DO IA=1,NP >*/
    i__1 = genin_1.np;
    for (ia = 1; ia <= i__1; ++ia) {
/*<          DO IB=1,5 >*/
	for (ib = 1; ib <= 5; ++ib) {
/*<             PCM1(IB,IA)=PCM(IB,IA) >*/
	    pcm1[ib + ia * 5 - 6] = genout_1.pcm[ib + ia * 5 - 6];
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/* ********************************************************************** */
/* *** Lorentz Transformations in LAB system */
/* ********************************************************************** */
/* ***      BETA -> -BETA  was assigned previously for Lorentz transformation */
/* ***      from CMS to Lab system */
/* ********************************************************************** */
/*< 	BETA( >*/
    beta[3] = 1.f / sqrt(1.f - (beta[0] * beta[0] + beta[1] * beta[1] + beta[
	    2] * beta[2]));
/* Transform to LAB SYS */
/*< 	DO J=1,NP >*/
    i__1 = genin_1.np;
    for (j = 1; j <= i__1; ++j) {
/*<             DO IL=1,4 >*/
	for (il = 1; il <= 4; ++il) {
/*<                P1(IL)=PCM1(IL,J) >*/
	    p1[il - 1] = pcm1[il + j * 5 - 6];
/*<             END DO >*/
	}
/*<             CALL GLOREN(BETA,P1,PL) >*/
	gloren_(beta, p1, pl);
/*<             DO IL=1,4 >*/
	for (il = 1; il <= 4; ++il) {
/*<                PLAB1(IL,J)=PL(IL) >*/
	    plab1[il + j * 5 - 6] = pl[il - 1];
/*<             END DO >*/
	}
/*<             PLAB1(5,J)=SQRT(PL(1)*PL(1)+PL(2)*PL(2)+PL(3)*PL(3)) >*/
	plab1[j * 5 - 1] = sqrt(pl[0] * pl[0] + pl[1] * pl[1] + pl[2] * pl[2])
		;
/*< 	END DO >*/
    }
/* ********************************************************************** */
/* ***   Computing t for two pions production (phase space) and */
/* ***   cheking for right values */
/* ********************************************************************** */
/*<       IF (JCH.GE.15.AND.JCH.LE.18) THEN >*/
    if (*jch >= 15 && *jch <= 18) {
/*<          DO IL=1,4 >*/
	for (il = 1; il <= 4; ++il) {
/*<             PN_T(IL) = PLAB1(IL,3)-PN(IL,1) >*/
	    pn_t__[il - 1] = plab1[il + 9] - mom_nucleon__1.pn[il - 1];
/*<          END DO >*/
	}
/*<          T = PN_T(1)*PN_T(1)+PN_T(2)*PN_T(2)+PN_T(3)*PN_T(3) >*/
	t = pn_t__[0] * pn_t__[0] + pn_t__[1] * pn_t__[1] + pn_t__[2] * 
		pn_t__[2];
/*<          T = PN_T(4)*PN_T(4)-T >*/
	t = pn_t__[3] * pn_t__[3] - t;
/*<          U = EXP(T_COEFF*T) >*/
	u = exp(t_coeff__ * t);
/*<          call ranmar(vran,1) >*/
	ranmar_(vran, &c__1);
/*<          YRNDM=vran(1) >*/
	yrndm = vran[0];
/*<          IF (YRNDM.GT.U) GOTO 20 >*/
	if (yrndm > u) {
	    goto L20;
	}
/*<       END IF >*/
    }
/* ********************************************************************** */
/* ********************************************************************** */
/* *** Second leve: decays */
/* ********************************************************************** */
/* ********************************************************************** */
/*<       IF(JDEC(JCH).EQ.0) THEN >*/
    if (channels_1.jdec[*jch - 1] == 0) {
/* ********************************************************************** */
/* ***   special channel on He3, Delta knockout */
/* Check conservation law for channel #37 */
/* in the lab system. */
/* For other cases it has to be done */
/* outside the generator body... */
/* ********************************************************************** */
/*<          IF (JCH.EQ.37) THEN >*/
	if (*jch == 37) {
/*<        >*/
	    delta[2] = mom_nucleon__1.pn[2] - (plab1[2] + plab1[7] + plab1[12]
		     + plab1[17] + plab1[22] + plab1[27]);
/*<        >*/
	    delta[3] = mom_nucleon__1.pn[3] - (plab1[3] + plab1[8] + plab1[13]
		     + plab1[18] + plab1[23] + plab1[28]);
/*<          END IF >*/
	}
/*<        >*/
	if (dabs(delta[0]) > 1e-4f || dabs(delta[1]) > 1e-4f || dabs(delta[2])
		 > 1e-4f || dabs(delta[3]) > 1e-4f) {
/*<             NERR = NERR+1 >*/
	    ++flag_gen__1.nerr;
/*<             GOTO 11 >*/
	    goto L11;
/*<          END IF >*/
	}
/* ********************************************************************** */
/*   Put angular distribution (C.M.S.) for the first particle in the */
/*    array for each channel */
/* ********************************************************************** */
/*<          IF(PCM(5,1).GT.1.E-07) THEN >*/
	if (genout_1.pcm[4] > 1e-7f) {
/*<             ARGUMENT=PCM(3,1)/PCM(5,1) >*/
	    argument = genout_1.pcm[2] / genout_1.pcm[4];
/*<             IF(ABS(ARGUMENT).LE.1.) THEN  >*/
	    if (dabs(argument) <= 1.f) {
/*               CALL HFILL(-JCH,ACOSDD(ARGUMENT),0.,1.) */
/*<             ELSEIF (ARGUMENT.LT.-1.) THEN >*/
	    } else if (argument < -1.f) {
/*               CALL HFILL(-JCH,180.,0.,1.) */
/*<             ELSEIF (ARGUMENT.GT.1.) THEN >*/
	    } else if (argument > 1.f) {
/*               CALL HFILL(-JCH,0.,0.,1.) */
/*<             ENDIF >*/
	    }
/*<          ENDIF >*/
	}
/*<       ELSE >*/
    } else {
/* ********************************************************************** */
/* *** DECAYS */
/* ********************************************************************** */
/* ******************************************** */
/* *** Calling the routine for extraction of kinematics */
/* *** of first level particles                       !M.Mirazita */
/* ******************************************** */
/* *** Prova per aggiunta canale 35 */
/* *** Finiva a jch=14 */
/* ******************************************** */
/*<          if ( ((jch.ge.4).and.(jch.le.14)).or.(jch.eq.35) ) then >*/
	if (*jch >= 4 && *jch <= 14 || *jch == 35) {
/*<             call prim_kine >*/
	    prim_kine__();
/*<          endif >*/
	}
/* ************************************************************************* */
/*<          IF(PCM(5,1).GT.1.E-07) THEN >*/
	if (genout_1.pcm[4] > 1e-7f) {
/*<             ARGUMENT=PCM(3,1)/PCM(5,1) >*/
	    argument = genout_1.pcm[2] / genout_1.pcm[4];
/*<             IF(ABS(ARGUMENT).LE.1.) THEN  >*/
	    if (dabs(argument) <= 1.f) {
/*               CALL HFILL(-JCH,ACOSDD(ARGUMENT),0.,1.) */
/*<             ELSEIF (ARGUMENT.LT.-1.) THEN >*/
	    } else if (argument < -1.f) {
/*               CALL HFILL(-JCH,180.,0.,1.) */
/*<             ELSEIF (ARGUMENT.GT.1.) THEN >*/
	    } else if (argument > 1.f) {
/*               CALL HFILL(-JCH,0.,0.,1.) */
/*<             ENDIF >*/
	    }
/*<          ENDIF >*/
	}
/* ********************************************************************** */
/* ***   Choosing Branching Ratio for second level channel definition */
/* ********************************************************************** */
/*<          JBR=1 >*/
	jbr = 1;
/*<          IF (NCH(JCH).GE.2) THEN >*/
	if (channels_1.nch[*jch - 1] >= 2) {
/*<             PROMUL(0) = 0. >*/
	    promul[0] = 0.f;
/*<             DO K=1,NCH(JCH) >*/
	    i__1 = channels_1.nch[*jch - 1];
	    for (k = 1; k <= i__1; ++k) {
/*<                PROMUL(K) = PROMUL(K-1)+BR(JCH,K) >*/
		promul[k] = promul[k - 1] + channels_1.br[*jch + k * 40 - 41];
/*<             END DO >*/
	    }
/* *** random */
/*<             call ranmar(vran,1) >*/
	    ranmar_(vran, &c__1);
/*<             Q=vran(1) >*/
	    q = vran[0];
/*<             DO K=1,NCH(JCH) >*/
	    i__1 = channels_1.nch[*jch - 1];
	    for (k = 1; k <= i__1; ++k) {
/*<                IF (Q.LE.PROMUL(K)) GOTO 66 >*/
		if (q <= promul[k]) {
		    goto L66;
		}
/*<             END DO >*/
	    }
/*<  66         JBR = K >*/
L66:
	    jbr = k;
/*<          END IF >*/
	}
/*<          NP=NPLEV(JCH,2) >*/
	genin_1.np = channels_1.nplev[*jch + 39];
/*<          JD=JDEC(JCH) >*/
	jd = channels_1.jdec[*jch - 1];
/* ********************************************************************** */
/* ***   Defining  BETA =-(P/ETOT) of decayng particle */
/* ********************************************************************** */
/*<          IF (JD.GT.6) THEN >*/
	if (jd > 6) {
/*<             WRITE(15,111) JD,JCH,JBR,NP >*/
	    s_wsfe(&io___54);
	    do_fio(&c__1, (char *)&jd, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*jch), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&jbr, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&genin_1.np, (ftnlen)sizeof(integer));
	    e_wsfe();
/*<  111        FORMAT(1X,' JD=',I6,' JCH=',I6,' JBR=',I6,' NP=',I6) >*/
/*<          ENDIF >*/
	}
/*<          DO K=1,3 >*/
	for (k = 1; k <= 3; ++k) {
/*<             BETA(K)=-PLAB1(K,JD)/PLAB1(4,JD) >*/
	    beta[k - 1] = -plab1[k + jd * 5 - 6] / plab1[jd * 5 - 2];
/*<          ENDDO >*/
	}
/*<          BETA(4)=0. >*/
	beta[3] = 0.f;
/*<          W=WW >*/
	w = ww;
/*<          TOTMASS = 0. >*/
	totmass = 0.f;
/*<          LEV=2 >*/
	lev = 2;
/*<          DO JM=1,NP >*/
	i__1 = genin_1.np;
	for (jm = 1; jm <= i__1; ++jm) {
/*<             AMASS(JM)=RMASS(PARTIC(JCH,LEV,JBR,JM)) >*/
	    genin_1.amass[jm - 1] = rmass_(channels_1.partic + (*jch + (lev + 
		    (jbr + (jm << 2) << 2)) * 40 - 641) * 12, (ftnlen)12);
/*<             AMASS2(JM) = AMASS(JM)         >*/
	    amass2[jm - 1] = genin_1.amass[jm - 1];
/*<             TOTMASS = TOTMASS+AMASS(JM) >*/
	    totmass += genin_1.amass[jm - 1];
/*<          ENDDO >*/
	}
/*<          IF (W.LT.(TOTMASS+0.000001)) GOTO 11 >*/
	if (w < totmass + 1e-6f) {
	    goto L11;
	}
/* ********************************************************************** */
/* ***   First step: decay with GENBOD */
/* ********************************************************************** */
/*<          ECM = W >*/
	genin_1.ecm = w;
/*<          CALL GENBOD >*/
	genbod_();
/*<          DO IA=1,NP >*/
	i__1 = genin_1.np;
	for (ia = 1; ia <= i__1; ++ia) {
/*<             DO IB=1,5 >*/
	    for (ib = 1; ib <= 5; ++ib) {
/*<                PCM2(IB,IA)=PCM(IB,IA) >*/
		pcm2[ib + ia * 5 - 6] = genout_1.pcm[ib + ia * 5 - 6];
/*<             ENDDO >*/
	    }
/*<          ENDDO >*/
	}
/* ********************************************************************** */
/* ***      BETA -> -BETA  was assigned previously for Lorentz transformation */
/* ***      from CMS to Lab system */
/* ********************************************************************** */
/*<        >*/
	beta[3] = 1.f / sqrt(1.f - (beta[0] * beta[0] + beta[1] * beta[1] + 
		beta[2] * beta[2]));
/*<          DO J=1,NP >*/
	i__1 = genin_1.np;
	for (j = 1; j <= i__1; ++j) {
/*<             DO IL=1,4 >*/
	    for (il = 1; il <= 4; ++il) {
/*<                P2(IL)=PCM2(IL,J) >*/
		p2[il - 1] = pcm2[il + j * 5 - 6];
/*<             ENDDO >*/
	    }
/*<             CALL GLOREN(BETA,P2,PL) >*/
	    gloren_(beta, p2, pl);
/*<             DO IL=1,4 >*/
	    for (il = 1; il <= 4; ++il) {
/*<                PLAB2(IL,J)=PL(IL) >*/
		plab2[il + j * 5 - 6] = pl[il - 1];
/*<             ENDDO >*/
	    }
/*<             PLAB2(5,J)=SQRT(PL(1)*PL(1)+PL(2)*PL(2)+PL(3)*PL(3)) >*/
	    plab2[j * 5 - 1] = sqrt(pl[0] * pl[0] + pl[1] * pl[1] + pl[2] * 
		    pl[2]);
/*<          ENDDO >*/
	}
/*<          IF (JCH.EQ.40) THEN    ! Check conservation law for channel #40 >*/
	if (*jch == 40) {
/* in the lab system. */
/* For other cases it has to be done */
/* outside the generator body... */
/*<             DELTA(1) = PLAB1(1,1)+PLAB2(1,1)+PLAB2(1,2)+PN(1,2)+PN(1,3) >*/
	    delta[0] = plab1[0] + plab2[0] + plab2[5] + mom_nucleon__1.pn[5] 
		    + mom_nucleon__1.pn[10];
/*<             DELTA(2) = PLAB1(2,1)+PLAB2(2,1)+PLAB2(2,2)+PN(2,2)+PN(2,3) >*/
	    delta[1] = plab1[1] + plab2[1] + plab2[6] + mom_nucleon__1.pn[6] 
		    + mom_nucleon__1.pn[11];
/*<             DELTA(3) = EGAM+PN(3,1)-(PLAB1(3,1)+PLAB2(3,1)+PLAB2(3,2)) >*/
	    delta[2] = montecarl_1.egam + mom_nucleon__1.pn[2] - (plab1[2] + 
		    plab2[2] + plab2[7]);
/*<             DELTA(4) = EGAM+PN(4,1)-(PLAB1(4,1)+PLAB2(4,1)+PLAB2(4,2)) >*/
	    delta[3] = montecarl_1.egam + mom_nucleon__1.pn[3] - (plab1[3] + 
		    plab2[3] + plab2[8]);
/*<          END IF >*/
	}
/*<        >*/
	if (dabs(delta[0]) > 1e-4f || dabs(delta[1]) > 1e-4f || dabs(delta[2])
		 > 1e-4f || dabs(delta[3]) > 1e-4f) {
/*<             NERR = NERR+1 >*/
	    ++flag_gen__1.nerr;
/*<             GOTO 11 >*/
	    goto L11;
/*<          END IF >*/
	}
/*<       END IF >*/
    }
/* ********************************************************************** */
/* ********************************************************************** */
/* *** Final momenta */
/* ********************************************************************** */
/* ********************************************************************** */
/*<       IF (TARGET.EQ.'deuteron') NUMPLUS = 1 >*/
    if (s_cmp(target_1.target, "deuteron", (ftnlen)12, (ftnlen)8) == 0) {
	numplus = 1;
    }
/*<       IF (TARGET.EQ.'He3') NUMPLUS = 2 >*/
    if (s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
	numplus = 2;
    }
/* ***    Canali 1,2,3,4,19,20,21,22,32 */
/*<        >*/
    if (*jch <= 4 || *jch == 19 || *jch == 20 || *jch == 21 || *jch == 22 || *
	    jch == 32) {
/*<          NUMPAR=NPLEV(JCH,1) >*/
	numpar = channels_1.nplev[*jch - 1];
/*<          NPART = NUMPAR+NUMPLUS >*/
	montecarl_1.npart = numpar + numplus;
/*<          P(1,1) = PLAB1(1,1) >*/
	montecarl_1.p[0] = plab1[0];
/*<          P(2,1) = PLAB1(2,1) >*/
	montecarl_1.p[1] = plab1[1];
/*<          P(3,1) = PLAB1(3,1) >*/
	montecarl_1.p[2] = plab1[2];
/*<          P(4,1) = AMASS1(1) >*/
	montecarl_1.p[3] = amass1[0];
/*<          I(1) = CODE(PARTIC(JCH,1,0,1)) >*/
	montecarl_1.i__[0] = code_(channels_1.partic + (*jch + 39) * 12, (
		ftnlen)12);
/*<          P(1,2) = PLAB1(1,2) >*/
	montecarl_1.p[4] = plab1[5];
/*<          P(2,2) = PLAB1(2,2) >*/
	montecarl_1.p[5] = plab1[6];
/*<          P(3,2) = PLAB1(3,2) >*/
	montecarl_1.p[6] = plab1[7];
/*<          P(4,2) = AMASS1(2) >*/
	montecarl_1.p[7] = amass1[1];
/*<          I(2) = CODE(PARTIC(JCH,1,0,2)) >*/
	montecarl_1.i__[1] = code_(channels_1.partic + (*jch + 679) * 12, (
		ftnlen)12);
/*<          IF (TARGET.EQ.'deuteron'.OR.TARGET.EQ.'He3') THEN >*/
	if (s_cmp(target_1.target, "deuteron", (ftnlen)12, (ftnlen)8) == 0 || 
		s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<             P(1,3) = PN(1,2) >*/
	    montecarl_1.p[8] = mom_nucleon__1.pn[5];
/*<             P(2,3) = PN(2,2) >*/
	    montecarl_1.p[9] = mom_nucleon__1.pn[6];
/*<             P(3,3) = PN(3,2) >*/
	    montecarl_1.p[10] = mom_nucleon__1.pn[7];
/*<             P(4,3) = RMASS(PARTIC(JCH,1,0,3)) >*/
	    montecarl_1.p[11] = rmass_(channels_1.partic + (*jch + 1319) * 12,
		     (ftnlen)12);
/*<             I(3) = CODE(PARTIC(JCH,1,0,3)) >*/
	    montecarl_1.i__[2] = code_(channels_1.partic + (*jch + 1319) * 12,
		     (ftnlen)12);
/*<             IF (TARGET.EQ.'He3') THEN >*/
	    if (s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<                P(1,4) = PN(1,3) >*/
		montecarl_1.p[12] = mom_nucleon__1.pn[10];
/*<                P(2,4) = PN(2,3) >*/
		montecarl_1.p[13] = mom_nucleon__1.pn[11];
/*<                P(3,4) = PN(3,3) >*/
		montecarl_1.p[14] = mom_nucleon__1.pn[12];
/*<                P(4,4) = RMASS(PARTIC(JCH,1,0,4)) >*/
		montecarl_1.p[15] = rmass_(channels_1.partic + (*jch + 1959) *
			 12, (ftnlen)12);
/*<                I(4) = CODE(PARTIC(JCH,1,0,4)) >*/
		montecarl_1.i__[3] = code_(channels_1.partic + (*jch + 1959) *
			 12, (ftnlen)12);
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<          GOTO 99 >*/
	goto L99;
/*<       END IF >*/
    }
/* ***    Canali 5,6,7,8,9,10 */
/*<       IF (JCH.GE.5.AND.JCH.LE.10) THEN >*/
    if (*jch >= 5 && *jch <= 10) {
/*<          NUMPAR=NPLEV(JCH,1)+NPLEV(JCH,2)-1 >*/
	numpar = channels_1.nplev[*jch - 1] + channels_1.nplev[*jch + 39] - 1;
/*<          NPART = NUMPAR+NUMPLUS >*/
	montecarl_1.npart = numpar + numplus;
/*<          P(1,1) = PLAB1(1,1) >*/
	montecarl_1.p[0] = plab1[0];
/*<          P(2,1) = PLAB1(2,1) >*/
	montecarl_1.p[1] = plab1[1];
/*<          P(3,1) = PLAB1(3,1) >*/
	montecarl_1.p[2] = plab1[2];
/*<          P(4,1) = AMASS1(1) >*/
	montecarl_1.p[3] = amass1[0];
/*<          I(1) = CODE(PARTIC(JCH,1,0,1)) >*/
	montecarl_1.i__[0] = code_(channels_1.partic + (*jch + 39) * 12, (
		ftnlen)12);
/*<          P(1,2) = PLAB2(1,1) >*/
	montecarl_1.p[4] = plab2[0];
/*<          P(2,2) = PLAB2(2,1) >*/
	montecarl_1.p[5] = plab2[1];
/*<          P(3,2) = PLAB2(3,1) >*/
	montecarl_1.p[6] = plab2[2];
/*<          P(4,2) = AMASS2(1) >*/
	montecarl_1.p[7] = amass2[0];
/*<          I(2) = CODE(PARTIC(JCH,2,JBR,1)) >*/
	montecarl_1.i__[1] = code_(channels_1.partic + (*jch + ((jbr + 4 << 2)
		 + 2) * 40 - 641) * 12, (ftnlen)12);
/*<          P(1,3) = PLAB2(1,2) >*/
	montecarl_1.p[8] = plab2[5];
/*<          P(2,3) = PLAB2(2,2) >*/
	montecarl_1.p[9] = plab2[6];
/*<          P(3,3) = PLAB2(3,2) >*/
	montecarl_1.p[10] = plab2[7];
/*<          P(4,3) = AMASS2(2) >*/
	montecarl_1.p[11] = amass2[1];
/*<          I(3) = CODE(PARTIC(JCH,2,JBR,2)) >*/
	montecarl_1.i__[2] = code_(channels_1.partic + (*jch + ((jbr + 8 << 2)
		 + 2) * 40 - 641) * 12, (ftnlen)12);
/*<          IF (TARGET.EQ.'deuteron'.OR.TARGET.EQ.'He3') THEN >*/
	if (s_cmp(target_1.target, "deuteron", (ftnlen)12, (ftnlen)8) == 0 || 
		s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<             P(1,4) = PN(1,2) >*/
	    montecarl_1.p[12] = mom_nucleon__1.pn[5];
/*<             P(2,4) = PN(2,2) >*/
	    montecarl_1.p[13] = mom_nucleon__1.pn[6];
/*<             P(3,4) = PN(3,2) >*/
	    montecarl_1.p[14] = mom_nucleon__1.pn[7];
/*<             P(4,4) = RMASS(PARTIC(JCH,1,0,3)) >*/
	    montecarl_1.p[15] = rmass_(channels_1.partic + (*jch + 1319) * 12,
		     (ftnlen)12);
/*<             I(4) = CODE(PARTIC(JCH,1,0,3)) >*/
	    montecarl_1.i__[3] = code_(channels_1.partic + (*jch + 1319) * 12,
		     (ftnlen)12);
/*<             IF (TARGET.EQ.'He3') THEN >*/
	    if (s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<                P(1,5) = PN(1,3) >*/
		montecarl_1.p[16] = mom_nucleon__1.pn[10];
/*<                P(2,5) = PN(2,3) >*/
		montecarl_1.p[17] = mom_nucleon__1.pn[11];
/*<                P(3,5) = PN(3,3) >*/
		montecarl_1.p[18] = mom_nucleon__1.pn[12];
/*<                P(4,5) = RMASS(PARTIC(JCH,1,0,4)) >*/
		montecarl_1.p[19] = rmass_(channels_1.partic + (*jch + 1959) *
			 12, (ftnlen)12);
/*<                I(5) = CODE(PARTIC(JCH,1,0,4)) >*/
		montecarl_1.i__[4] = code_(channels_1.partic + (*jch + 1959) *
			 12, (ftnlen)12);
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
    if (*jch >= 11 && *jch <= 14) {
/*<          NUMPAR=NPLEV(JCH,1)+NPLEV(JCH,2)-1 >*/
	numpar = channels_1.nplev[*jch - 1] + channels_1.nplev[*jch + 39] - 1;
/*<          NPART = NUMPAR+NUMPLUS >*/
	montecarl_1.npart = numpar + numplus;
/*<          P(1,1) = PLAB2(1,1) >*/
	montecarl_1.p[0] = plab2[0];
/*<          P(2,1) = PLAB2(2,1) >*/
	montecarl_1.p[1] = plab2[1];
/*<          P(3,1) = PLAB2(3,1) >*/
	montecarl_1.p[2] = plab2[2];
/*<          P(4,1) = AMASS2(1) >*/
	montecarl_1.p[3] = amass2[0];
/*<          I(1) = CODE(PARTIC(JCH,2,JBR,1)) >*/
	montecarl_1.i__[0] = code_(channels_1.partic + (*jch + ((jbr + 4 << 2)
		 + 2) * 40 - 641) * 12, (ftnlen)12);
/*<          P(1,2) = PLAB2(1,2) >*/
	montecarl_1.p[4] = plab2[5];
/*<          P(2,2) = PLAB2(2,2) >*/
	montecarl_1.p[5] = plab2[6];
/*<          P(3,2) = PLAB2(3,2) >*/
	montecarl_1.p[6] = plab2[7];
/*<          P(4,2) = AMASS2(2) >*/
	montecarl_1.p[7] = amass2[1];
/*<          I(2) = CODE(PARTIC(JCH,2,JBR,2)) >*/
	montecarl_1.i__[1] = code_(channels_1.partic + (*jch + ((jbr + 8 << 2)
		 + 2) * 40 - 641) * 12, (ftnlen)12);
/*<          P(1,3) = PLAB1(1,2) >*/
	montecarl_1.p[8] = plab1[5];
/*<          P(2,3) = PLAB1(2,2) >*/
	montecarl_1.p[9] = plab1[6];
/*<          P(3,3) = PLAB1(3,2) >*/
	montecarl_1.p[10] = plab1[7];
/*<          P(4,3) = AMASS1(2) >*/
	montecarl_1.p[11] = amass1[1];
/*<          I(3) = CODE(PARTIC(JCH,1,0,2)) >*/
	montecarl_1.i__[2] = code_(channels_1.partic + (*jch + 679) * 12, (
		ftnlen)12);
/*<          IF (TARGET.EQ.'deuteron'.OR.TARGET.EQ.'He3') THEN >*/
	if (s_cmp(target_1.target, "deuteron", (ftnlen)12, (ftnlen)8) == 0 || 
		s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<             P(1,4) = PN(1,2) >*/
	    montecarl_1.p[12] = mom_nucleon__1.pn[5];
/*<             P(2,4) = PN(2,2) >*/
	    montecarl_1.p[13] = mom_nucleon__1.pn[6];
/*<             P(3,4) = PN(3,2) >*/
	    montecarl_1.p[14] = mom_nucleon__1.pn[7];
/*<             P(4,4) = RMASS(PARTIC(JCH,1,0,3)) >*/
	    montecarl_1.p[15] = rmass_(channels_1.partic + (*jch + 1319) * 12,
		     (ftnlen)12);
/*<             I(4) = CODE(PARTIC(JCH,1,0,3)) >*/
	    montecarl_1.i__[3] = code_(channels_1.partic + (*jch + 1319) * 12,
		     (ftnlen)12);
/*<             IF (TARGET.EQ.'He3') THEN >*/
	    if (s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<                P(1,5) = PN(1,3) >*/
		montecarl_1.p[16] = mom_nucleon__1.pn[10];
/*<                P(2,5) = PN(2,3) >*/
		montecarl_1.p[17] = mom_nucleon__1.pn[11];
/*<                P(3,5) = PN(3,3) >*/
		montecarl_1.p[18] = mom_nucleon__1.pn[12];
/*<                P(4,5) = RMASS(PARTIC(JCH,1,0,4)) >*/
		montecarl_1.p[19] = rmass_(channels_1.partic + (*jch + 1959) *
			 12, (ftnlen)12);
/*<                I(5) = CODE(PARTIC(JCH,1,0,4)) >*/
		montecarl_1.i__[4] = code_(channels_1.partic + (*jch + 1959) *
			 12, (ftnlen)12);
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<          GOTO 99 >*/
	goto L99;
/*<       END IF >*/
    }
/* ***    Canali 15,16,17,18 */
/*<       IF (JCH.GE.15.AND.JCH.LE.18) THEN >*/
    if (*jch >= 15 && *jch <= 18) {
/*<          NUMPAR=NPLEV(JCH,1) >*/
	numpar = channels_1.nplev[*jch - 1];
/*<          NPART = NUMPAR+NUMPLUS >*/
	montecarl_1.npart = numpar + numplus;
/*<          P(1,1) = PLAB1(1,1) >*/
	montecarl_1.p[0] = plab1[0];
/*<          P(2,1) = PLAB1(2,1) >*/
	montecarl_1.p[1] = plab1[1];
/*<          P(3,1) = PLAB1(3,1) >*/
	montecarl_1.p[2] = plab1[2];
/*<          P(4,1) = AMASS1(1) >*/
	montecarl_1.p[3] = amass1[0];
/*<          I(1) = CODE(PARTIC(JCH,1,0,1)) >*/
	montecarl_1.i__[0] = code_(channels_1.partic + (*jch + 39) * 12, (
		ftnlen)12);
/*<          P(1,2) = PLAB1(1,2) >*/
	montecarl_1.p[4] = plab1[5];
/*<          P(2,2) = PLAB1(2,2) >*/
	montecarl_1.p[5] = plab1[6];
/*<          P(3,2) = PLAB1(3,2) >*/
	montecarl_1.p[6] = plab1[7];
/*<          P(4,2) = AMASS1(2) >*/
	montecarl_1.p[7] = amass1[1];
/*<          I(2) = CODE(PARTIC(JCH,1,0,2)) >*/
	montecarl_1.i__[1] = code_(channels_1.partic + (*jch + 679) * 12, (
		ftnlen)12);
/*<          P(1,3) = PLAB1(1,3) >*/
	montecarl_1.p[8] = plab1[10];
/*<          P(2,3) = PLAB1(2,3) >*/
	montecarl_1.p[9] = plab1[11];
/*<          P(3,3) = PLAB1(3,3) >*/
	montecarl_1.p[10] = plab1[12];
/*<          P(4,3) = AMASS1(3) >*/
	montecarl_1.p[11] = amass1[2];
/*<          I(3) = CODE(PARTIC(JCH,1,0,3)) >*/
	montecarl_1.i__[2] = code_(channels_1.partic + (*jch + 1319) * 12, (
		ftnlen)12);
/*<          IF (TARGET.EQ.'deuteron'.OR.TARGET.EQ.'He3') THEN >*/
	if (s_cmp(target_1.target, "deuteron", (ftnlen)12, (ftnlen)8) == 0 || 
		s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<             P(1,4) = PN(1,2) >*/
	    montecarl_1.p[12] = mom_nucleon__1.pn[5];
/*<             P(2,4) = PN(2,2) >*/
	    montecarl_1.p[13] = mom_nucleon__1.pn[6];
/*<             P(3,4) = PN(3,2) >*/
	    montecarl_1.p[14] = mom_nucleon__1.pn[7];
/*<             P(4,4) = RMASS(PARTIC(JCH,1,0,4)) >*/
	    montecarl_1.p[15] = rmass_(channels_1.partic + (*jch + 1959) * 12,
		     (ftnlen)12);
/*<             I(4) = CODE(PARTIC(JCH,1,0,4)) >*/
	    montecarl_1.i__[3] = code_(channels_1.partic + (*jch + 1959) * 12,
		     (ftnlen)12);
/*<             IF (TARGET.EQ.'He3') THEN >*/
	    if (s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<                P(1,5) = PN(1,3) >*/
		montecarl_1.p[16] = mom_nucleon__1.pn[10];
/*<                P(2,5) = PN(2,3) >*/
		montecarl_1.p[17] = mom_nucleon__1.pn[11];
/*<                P(3,5) = PN(3,3) >*/
		montecarl_1.p[18] = mom_nucleon__1.pn[12];
/*<                P(4,5) = RMASS(PARTIC(JCH,1,0,5)) >*/
		montecarl_1.p[19] = rmass_(channels_1.partic + (*jch + 2599) *
			 12, (ftnlen)12);
/*<                I(5) = CODE(PARTIC(JCH,1,0,5)) >*/
		montecarl_1.i__[4] = code_(channels_1.partic + (*jch + 2599) *
			 12, (ftnlen)12);
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<          GOTO 99 >*/
	goto L99;
/*<       END IF >*/
    }
/* ***    Canali 23,25 */
/*<       IF (JCH.EQ.23.OR.JCH.EQ.25) THEN >*/
    if (*jch == 23 || *jch == 25) {
/*<          NUMPAR=NPLEV(JCH,1) >*/
	numpar = channels_1.nplev[*jch - 1];
/*<          NPART = NUMPAR+NUMPLUS >*/
	montecarl_1.npart = numpar + numplus;
/*<          P(1,1) = PLAB1(1,1) >*/
	montecarl_1.p[0] = plab1[0];
/*<          P(2,1) = PLAB1(2,1) >*/
	montecarl_1.p[1] = plab1[1];
/*<          P(3,1) = PLAB1(3,1) >*/
	montecarl_1.p[2] = plab1[2];
/*<          P(4,1) = AMASS1(1) >*/
	montecarl_1.p[3] = amass1[0];
/*<          I(1) = CODE(PARTIC(JCH,1,0,1)) >*/
	montecarl_1.i__[0] = code_(channels_1.partic + (*jch + 39) * 12, (
		ftnlen)12);
/*<          P(1,2) = PLAB1(1,2) >*/
	montecarl_1.p[4] = plab1[5];
/*<          P(2,2) = PLAB1(2,2) >*/
	montecarl_1.p[5] = plab1[6];
/*<          P(3,2) = PLAB1(3,2) >*/
	montecarl_1.p[6] = plab1[7];
/*<          P(4,2) = AMASS1(2) >*/
	montecarl_1.p[7] = amass1[1];
/*<          I(2) = CODE(PARTIC(JCH,1,0,2)) >*/
	montecarl_1.i__[1] = code_(channels_1.partic + (*jch + 679) * 12, (
		ftnlen)12);
/*<          P(1,3) = PLAB1(1,3) >*/
	montecarl_1.p[8] = plab1[10];
/*<          P(2,3) = PLAB1(2,3) >*/
	montecarl_1.p[9] = plab1[11];
/*<          P(3,3) = PLAB1(3,3) >*/
	montecarl_1.p[10] = plab1[12];
/*<          P(4,3) = AMASS1(3) >*/
	montecarl_1.p[11] = amass1[2];
/*<          I(3) = CODE(PARTIC(JCH,1,0,3)) >*/
	montecarl_1.i__[2] = code_(channels_1.partic + (*jch + 1319) * 12, (
		ftnlen)12);
/*<          P(1,4) = PLAB1(1,4) >*/
	montecarl_1.p[12] = plab1[15];
/*<          P(2,4) = PLAB1(2,4) >*/
	montecarl_1.p[13] = plab1[16];
/*<          P(3,4) = PLAB1(3,4) >*/
	montecarl_1.p[14] = plab1[17];
/*<          P(4,4) = AMASS1(4) >*/
	montecarl_1.p[15] = amass1[3];
/*<          I(4) = CODE(PARTIC(JCH,1,0,4)) >*/
	montecarl_1.i__[3] = code_(channels_1.partic + (*jch + 1959) * 12, (
		ftnlen)12);
/*<          IF (TARGET.EQ.'deuteron'.OR.TARGET.EQ.'He3') THEN >*/
	if (s_cmp(target_1.target, "deuteron", (ftnlen)12, (ftnlen)8) == 0 || 
		s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<             P(1,5) = PN(1,2) >*/
	    montecarl_1.p[16] = mom_nucleon__1.pn[5];
/*<             P(2,5) = PN(2,2) >*/
	    montecarl_1.p[17] = mom_nucleon__1.pn[6];
/*<             P(3,5) = PN(3,2) >*/
	    montecarl_1.p[18] = mom_nucleon__1.pn[7];
/*<             P(4,5) = RMASS(PARTIC(JCH,1,0,5)) >*/
	    montecarl_1.p[19] = rmass_(channels_1.partic + (*jch + 2599) * 12,
		     (ftnlen)12);
/*<             I(5) = CODE(PARTIC(JCH,1,0,5)) >*/
	    montecarl_1.i__[4] = code_(channels_1.partic + (*jch + 2599) * 12,
		     (ftnlen)12);
/*<             IF (TARGET.EQ.'He3') THEN >*/
	    if (s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<                P(1,6) = PN(1,3) >*/
		montecarl_1.p[20] = mom_nucleon__1.pn[10];
/*<                P(2,6) = PN(2,3) >*/
		montecarl_1.p[21] = mom_nucleon__1.pn[11];
/*<                P(3,6) = PN(3,3) >*/
		montecarl_1.p[22] = mom_nucleon__1.pn[12];
/*<                P(4,6) = RMASS(PARTIC(JCH,1,0,6)) >*/
		montecarl_1.p[23] = rmass_(channels_1.partic + (*jch + 3239) *
			 12, (ftnlen)12);
/*<                I(6) = CODE(PARTIC(JCH,1,0,6)) >*/
		montecarl_1.i__[5] = code_(channels_1.partic + (*jch + 3239) *
			 12, (ftnlen)12);
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<          GOTO 99 >*/
	goto L99;
/*<       END IF >*/
    }
/* ***    Channels #24 and #26 */
/*<       IF (JCH.EQ.24.OR.JCH.EQ.26) THEN >*/
    if (*jch == 24 || *jch == 26) {
/*<          NUMPAR=NPLEV(JCH,1) >*/
	numpar = channels_1.nplev[*jch - 1];
/*<          NPART = NUMPAR+NUMPLUS >*/
	montecarl_1.npart = numpar + numplus;
/*<          IF(JCH.EQ.24) THEN >*/
	if (*jch == 24) {
/*<             CALL GAMN_CHARG(4,1,PART_NUMB) >*/
	    gamn_charg__(&c__4, &c__1, part_numb__);
/*<          ELSE >*/
	} else {
/*<             CALL GAMN_CHARG(4,0,PART_NUMB) >*/
	    gamn_charg__(&c__4, &c__0, part_numb__);
/*<          ENDIF           >*/
	}
/*<          P(1,1) = PLAB1(1,1) >*/
	montecarl_1.p[0] = plab1[0];
/*<          P(2,1) = PLAB1(2,1) >*/
	montecarl_1.p[1] = plab1[1];
/*<          P(3,1) = PLAB1(3,1) >*/
	montecarl_1.p[2] = plab1[2];
/*<          P(4,1) = AMASS1(1) >*/
	montecarl_1.p[3] = amass1[0];
/*<          I(1) = PART_NUMB(2)            >*/
	montecarl_1.i__[0] = part_numb__[1];
/*<          P(1,2) = PLAB1(1,2) >*/
	montecarl_1.p[4] = plab1[5];
/*<          P(2,2) = PLAB1(2,2) >*/
	montecarl_1.p[5] = plab1[6];
/*<          P(3,2) = PLAB1(3,2) >*/
	montecarl_1.p[6] = plab1[7];
/*<          P(4,2) = AMASS1(2) >*/
	montecarl_1.p[7] = amass1[1];
/*<          I(2) = PART_NUMB(3)                        >*/
	montecarl_1.i__[1] = part_numb__[2];
/*<          P(1,3) = PLAB1(1,3) >*/
	montecarl_1.p[8] = plab1[10];
/*<          P(2,3) = PLAB1(2,3) >*/
	montecarl_1.p[9] = plab1[11];
/*<          P(3,3) = PLAB1(3,3) >*/
	montecarl_1.p[10] = plab1[12];
/*<          P(4,3) = AMASS1(3) >*/
	montecarl_1.p[11] = amass1[2];
/*<          I(3) = PART_NUMB(4) >*/
	montecarl_1.i__[2] = part_numb__[3];
/*<          P(1,4) = PLAB1(1,4) >*/
	montecarl_1.p[12] = plab1[15];
/*<          P(2,4) = PLAB1(2,4) >*/
	montecarl_1.p[13] = plab1[16];
/*<          P(3,4) = PLAB1(3,4) >*/
	montecarl_1.p[14] = plab1[17];
/*<          P(4,4) = AMASS1(4) >*/
	montecarl_1.p[15] = amass1[3];
/*<          I(4) = PART_NUMB(1) >*/
	montecarl_1.i__[3] = part_numb__[0];
/*<          IF (TARGET.EQ.'deuteron'.OR.TARGET.EQ.'He3') THEN  >*/
	if (s_cmp(target_1.target, "deuteron", (ftnlen)12, (ftnlen)8) == 0 || 
		s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<             P(1,5) = PN(1,2) >*/
	    montecarl_1.p[16] = mom_nucleon__1.pn[5];
/*<             P(2,5) = PN(2,2) >*/
	    montecarl_1.p[17] = mom_nucleon__1.pn[6];
/*<             P(3,5) = PN(3,2) >*/
	    montecarl_1.p[18] = mom_nucleon__1.pn[7];
/*<             P(4,5) = RMASS(PARTIC(JCH,1,0,5)) >*/
	    montecarl_1.p[19] = rmass_(channels_1.partic + (*jch + 2599) * 12,
		     (ftnlen)12);
/*<             I(5) = CODE(PARTIC(JCH,1,0,5)) >*/
	    montecarl_1.i__[4] = code_(channels_1.partic + (*jch + 2599) * 12,
		     (ftnlen)12);
/*<             IF (TARGET.EQ.'He3') THEN  >*/
	    if (s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<                P(1,6) = PN(1,3) >*/
		montecarl_1.p[20] = mom_nucleon__1.pn[10];
/*<                P(2,6) = PN(2,3) >*/
		montecarl_1.p[21] = mom_nucleon__1.pn[11];
/*<                P(3,6) = PN(3,3) >*/
		montecarl_1.p[22] = mom_nucleon__1.pn[12];
/*<                P(4,6) = RMASS(PARTIC(JCH,1,0,6)) >*/
		montecarl_1.p[23] = rmass_(channels_1.partic + (*jch + 3239) *
			 12, (ftnlen)12);
/*<                I(6) = CODE(PARTIC(JCH,1,0,6)) >*/
		montecarl_1.i__[5] = code_(channels_1.partic + (*jch + 3239) *
			 12, (ftnlen)12);
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<          GOTO 99 >*/
	goto L99;
/*<       END IF >*/
    }
/* ***    Canali 27,28 */
/*<       IF (JCH.EQ.27.OR.JCH.EQ.28) THEN >*/
    if (*jch == 27 || *jch == 28) {
/*<          NUMPAR=TOP_NUMBER >*/
	numpar = top_number__;
/*<          NPART = NUMPAR+NUMPLUS >*/
	montecarl_1.npart = numpar + numplus;
/*<          IF(JCH.EQ.27) THEN >*/
	if (*jch == 27) {
/*<             CALL GAMN_CHARG(TOP_NUMBER,1,PART_NUMB) >*/
	    gamn_charg__(&top_number__, &c__1, part_numb__);
/*<          ELSE >*/
	} else {
/*<             CALL GAMN_CHARG(TOP_NUMBER,0,PART_NUMB) >*/
	    gamn_charg__(&top_number__, &c__0, part_numb__);
/*<          ENDIF >*/
	}
/*<         DO JVAR=1,TOP_NUMBER >*/
	i__1 = top_number__;
	for (jvar = 1; jvar <= i__1; ++jvar) {
/*<            P(1,JVAR) = PLAB1(1,JVAR) >*/
	    montecarl_1.p[(jvar << 2) - 4] = plab1[jvar * 5 - 5];
/*<            P(2,JVAR) = PLAB1(2,JVAR) >*/
	    montecarl_1.p[(jvar << 2) - 3] = plab1[jvar * 5 - 4];
/*<            P(3,JVAR) = PLAB1(3,JVAR) >*/
	    montecarl_1.p[(jvar << 2) - 2] = plab1[jvar * 5 - 3];
/*<            P(4,JVAR) = AMASS1(JVAR) >*/
	    montecarl_1.p[(jvar << 2) - 1] = amass1[jvar - 1];
/*<            IF(JVAR.EQ.TOP_NUMBER) THEN >*/
	    if (jvar == top_number__) {
/*<               I(JVAR) = PART_NUMB(1) >*/
		montecarl_1.i__[jvar - 1] = part_numb__[0];
/*<            ELSE >*/
	    } else {
/*<               I(JVAR) = PART_NUMB(JVAR+1) >*/
		montecarl_1.i__[jvar - 1] = part_numb__[jvar];
/*<            ENDIF >*/
	    }
/*<         ENDDO >*/
	}
/*<         IF (TARGET.EQ.'deuteron'.OR.TARGET.EQ.'He3') THEN >*/
	if (s_cmp(target_1.target, "deuteron", (ftnlen)12, (ftnlen)8) == 0 || 
		s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<            P(1,TOP_NUMBER+1) = PN(1,2) >*/
	    montecarl_1.p[(top_number__ + 1 << 2) - 4] = mom_nucleon__1.pn[5];
/*<            P(2,TOP_NUMBER+1) = PN(2,2) >*/
	    montecarl_1.p[(top_number__ + 1 << 2) - 3] = mom_nucleon__1.pn[6];
/*<            P(3,TOP_NUMBER+1) = PN(3,2) >*/
	    montecarl_1.p[(top_number__ + 1 << 2) - 2] = mom_nucleon__1.pn[7];
/*<            P(4,TOP_NUMBER+1) = RMASS(PARTIC(JCH,1,0,6)) >*/
	    montecarl_1.p[(top_number__ + 1 << 2) - 1] = rmass_(
		    channels_1.partic + (*jch + 3239) * 12, (ftnlen)12);
/*<            I(TOP_NUMBER+1) = CODE(PARTIC(JCH,1,0,6)) >*/
	    montecarl_1.i__[top_number__] = code_(channels_1.partic + (*jch + 
		    3239) * 12, (ftnlen)12);
/*<            IF (TARGET.EQ.'He3') THEN >*/
	    if (s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<               P(1,TOP_NUMBER+2) = PN(1,3) >*/
		montecarl_1.p[(top_number__ + 2 << 2) - 4] = 
			mom_nucleon__1.pn[10];
/*<               P(2,TOP_NUMBER+2) = PN(2,3) >*/
		montecarl_1.p[(top_number__ + 2 << 2) - 3] = 
			mom_nucleon__1.pn[11];
/*<               P(3,TOP_NUMBER+2) = PN(3,3) >*/
		montecarl_1.p[(top_number__ + 2 << 2) - 2] = 
			mom_nucleon__1.pn[12];
/*<               P(4,TOP_NUMBER+2) = RMASS(PARTIC(JCH,1,0,7)) >*/
		montecarl_1.p[(top_number__ + 2 << 2) - 1] = rmass_(
			channels_1.partic + (*jch + 3879) * 12, (ftnlen)12);
/*<               I(TOP_NUMBER+2) = CODE(PARTIC(JCH,1,0,7)) >*/
		montecarl_1.i__[top_number__ + 1] = code_(channels_1.partic + 
			(*jch + 3879) * 12, (ftnlen)12);
/*<            END IF >*/
	    }
/*<         END IF >*/
	}
/*<         GOTO 99 >*/
	goto L99;
/*<       END IF >*/
    }
/* ***    Canali 29,30 */
/*<       IF (JCH.EQ.29.OR.JCH.EQ.30) THEN >*/
    if (*jch == 29 || *jch == 30) {
/*<          NUMPAR=NPLEV(JCH,1) >*/
	numpar = channels_1.nplev[*jch - 1];
/*<          NPART = NUMPAR+NUMPLUS >*/
	montecarl_1.npart = numpar + numplus;
/*<          P(1,1) = PLAB1(1,1) >*/
	montecarl_1.p[0] = plab1[0];
/*<          P(2,1) = PLAB1(2,1) >*/
	montecarl_1.p[1] = plab1[1];
/*<          P(3,1) = PLAB1(3,1) >*/
	montecarl_1.p[2] = plab1[2];
/*<          P(4,1) = AMASS1(1) >*/
	montecarl_1.p[3] = amass1[0];
/*<          I(1) = CODE(PARTIC(JCH,1,0,1)) >*/
	montecarl_1.i__[0] = code_(channels_1.partic + (*jch + 39) * 12, (
		ftnlen)12);
/*<          P(1,2) = PLAB1(1,2) >*/
	montecarl_1.p[4] = plab1[5];
/*<          P(2,2) = PLAB1(2,2) >*/
	montecarl_1.p[5] = plab1[6];
/*<          P(3,2) = PLAB1(3,2) >*/
	montecarl_1.p[6] = plab1[7];
/*<          P(4,2) = AMASS1(2) >*/
	montecarl_1.p[7] = amass1[1];
/*<          I(2) = CODE(PARTIC(JCH,1,0,2)) >*/
	montecarl_1.i__[1] = code_(channels_1.partic + (*jch + 679) * 12, (
		ftnlen)12);
/*<          IF (TARGET.EQ.'deuteron'.OR.TARGET.EQ.'He3') THEN >*/
	if (s_cmp(target_1.target, "deuteron", (ftnlen)12, (ftnlen)8) == 0 || 
		s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<             P(1,3) = PN(1,2) >*/
	    montecarl_1.p[8] = mom_nucleon__1.pn[5];
/*<             P(2,3) = PN(2,2) >*/
	    montecarl_1.p[9] = mom_nucleon__1.pn[6];
/*<             P(3,3) = PN(3,2) >*/
	    montecarl_1.p[10] = mom_nucleon__1.pn[7];
/*<             P(4,3) = RMASS(PARTIC(JCH,1,0,3)) >*/
	    montecarl_1.p[11] = rmass_(channels_1.partic + (*jch + 1319) * 12,
		     (ftnlen)12);
/*<             I(3) = CODE(PARTIC(JCH,1,0,3)) >*/
	    montecarl_1.i__[2] = code_(channels_1.partic + (*jch + 1319) * 12,
		     (ftnlen)12);
/*<             IF (TARGET.EQ.'He3') THEN >*/
	    if (s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<                P(1,4) = PN(1,3) >*/
		montecarl_1.p[12] = mom_nucleon__1.pn[10];
/*<                P(2,4) = PN(2,3) >*/
		montecarl_1.p[13] = mom_nucleon__1.pn[11];
/*<                P(3,4) = PN(3,3) >*/
		montecarl_1.p[14] = mom_nucleon__1.pn[12];
/*<                P(4,4) = RMASS(PARTIC(JCH,1,0,4)) >*/
		montecarl_1.p[15] = rmass_(channels_1.partic + (*jch + 1959) *
			 12, (ftnlen)12);
/*<                I(4) = CODE(PARTIC(JCH,1,0,4)) >*/
		montecarl_1.i__[3] = code_(channels_1.partic + (*jch + 1959) *
			 12, (ftnlen)12);
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<          GOTO 99 >*/
	goto L99;
/*<       END IF >*/
    }
/* ***    Canale 31 */
/*<       IF (JCH.EQ.31) THEN >*/
    if (*jch == 31) {
/*<          NUMPAR=NPLEV(JCH,1) >*/
	numpar = channels_1.nplev[*jch - 1];
/*<          NPART = NUMPAR+NUMPLUS >*/
	montecarl_1.npart = numpar + numplus;
/*<          P(1,1) = PLAB1(1,1) >*/
	montecarl_1.p[0] = plab1[0];
/*<          P(2,1) = PLAB1(2,1) >*/
	montecarl_1.p[1] = plab1[1];
/*<          P(3,1) = PLAB1(3,1) >*/
	montecarl_1.p[2] = plab1[2];
/*<          P(4,1) = AMASS1(1) >*/
	montecarl_1.p[3] = amass1[0];
/*<          I(1) = CODE(PARTIC(JCH,1,0,1)) >*/
	montecarl_1.i__[0] = code_(channels_1.partic + (*jch + 39) * 12, (
		ftnlen)12);
/*<          P(1,2) = PLAB1(1,2) >*/
	montecarl_1.p[4] = plab1[5];
/*<          P(2,2) = PLAB1(2,2) >*/
	montecarl_1.p[5] = plab1[6];
/*<          P(3,2) = PLAB1(3,2) >*/
	montecarl_1.p[6] = plab1[7];
/*<          P(4,2) = AMASS1(2) >*/
	montecarl_1.p[7] = amass1[1];
/*<          I(2) = CODE(PARTIC(JCH,1,0,2)) >*/
	montecarl_1.i__[1] = code_(channels_1.partic + (*jch + 679) * 12, (
		ftnlen)12);
/*<          P(1,3) = PLAB1(1,3) >*/
	montecarl_1.p[8] = plab1[10];
/*<          P(2,3) = PLAB1(2,3) >*/
	montecarl_1.p[9] = plab1[11];
/*<          P(3,3) = PLAB1(3,3) >*/
	montecarl_1.p[10] = plab1[12];
/*<          P(4,3) = AMASS1(3) >*/
	montecarl_1.p[11] = amass1[2];
/*<          I(3) = CODE(PARTIC(JCH,1,0,3)) >*/
	montecarl_1.i__[2] = code_(channels_1.partic + (*jch + 1319) * 12, (
		ftnlen)12);
/*<          IF (TARGET.EQ.'deuteron'.OR.TARGET.EQ.'He3') THEN >*/
	if (s_cmp(target_1.target, "deuteron", (ftnlen)12, (ftnlen)8) == 0 || 
		s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<             P(1,4) = PN(1,2) >*/
	    montecarl_1.p[12] = mom_nucleon__1.pn[5];
/*<             P(2,4) = PN(2,2) >*/
	    montecarl_1.p[13] = mom_nucleon__1.pn[6];
/*<             P(3,4) = PN(3,2) >*/
	    montecarl_1.p[14] = mom_nucleon__1.pn[7];
/*<             P(4,4) = RMASS(PARTIC(JCH,1,0,4)) >*/
	    montecarl_1.p[15] = rmass_(channels_1.partic + (*jch + 1959) * 12,
		     (ftnlen)12);
/*<             I(4) = CODE(PARTIC(JCH,1,0,4)) >*/
	    montecarl_1.i__[3] = code_(channels_1.partic + (*jch + 1959) * 12,
		     (ftnlen)12);
/*<             IF (TARGET.EQ.'He3') THEN >*/
	    if (s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<                P(1,5) = PN(1,3) >*/
		montecarl_1.p[16] = mom_nucleon__1.pn[10];
/*<                P(2,5) = PN(2,3) >*/
		montecarl_1.p[17] = mom_nucleon__1.pn[11];
/*<                P(3,5) = PN(3,3) >*/
		montecarl_1.p[18] = mom_nucleon__1.pn[12];
/*<                P(4,5) = RMASS(PARTIC(JCH,1,0,5)) >*/
		montecarl_1.p[19] = rmass_(channels_1.partic + (*jch + 2599) *
			 12, (ftnlen)12);
/*<                I(5) = CODE(PARTIC(JCH,1,0,5)) >*/
		montecarl_1.i__[4] = code_(channels_1.partic + (*jch + 2599) *
			 12, (ftnlen)12);
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<          GOTO 99 >*/
	goto L99;
/*<       END IF >*/
    }
/* ***    Canale 33 */
/*<       IF (JCH.EQ.33) THEN >*/
    if (*jch == 33) {
/*<          NUMPAR=NPLEV(JCH,1) >*/
	numpar = channels_1.nplev[*jch - 1];
/*<          NPART = NUMPAR >*/
	montecarl_1.npart = numpar;
/*<          P(1,1) = PLAB1(1,1) >*/
	montecarl_1.p[0] = plab1[0];
/*<          P(2,1) = PLAB1(2,1) >*/
	montecarl_1.p[1] = plab1[1];
/*<          P(3,1) = PLAB1(3,1) >*/
	montecarl_1.p[2] = plab1[2];
/*<          P(4,1) = AMASS1(1) >*/
	montecarl_1.p[3] = amass1[0];
/*<          I(1) = CODE(PARTIC(JCH,1,0,1)) >*/
	montecarl_1.i__[0] = code_(channels_1.partic + (*jch + 39) * 12, (
		ftnlen)12);
/*<          P(1,2) = PLAB1(1,2) >*/
	montecarl_1.p[4] = plab1[5];
/*<          P(2,2) = PLAB1(2,2) >*/
	montecarl_1.p[5] = plab1[6];
/*<          P(3,2) = PLAB1(3,2) >*/
	montecarl_1.p[6] = plab1[7];
/*<          P(4,2) = AMASS1(2) >*/
	montecarl_1.p[7] = amass1[1];
/*<          I(2) = CODE(PARTIC(JCH,1,0,2)) >*/
	montecarl_1.i__[1] = code_(channels_1.partic + (*jch + 679) * 12, (
		ftnlen)12);
/*<          IF (TARGET.EQ.'He3') THEN >*/
	if (s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<             P(1,3) = PN(1,2) >*/
	    montecarl_1.p[8] = mom_nucleon__1.pn[5];
/*<             P(2,3) = PN(2,2) >*/
	    montecarl_1.p[9] = mom_nucleon__1.pn[6];
/*<             P(3,3) = PN(3,2) >*/
	    montecarl_1.p[10] = mom_nucleon__1.pn[7];
/*<             P(4,3) = RMASS(PARTIC(JCH,1,0,3)) >*/
	    montecarl_1.p[11] = rmass_(channels_1.partic + (*jch + 1319) * 12,
		     (ftnlen)12);
/*<             I(3) = CODE(PARTIC(JCH,1,0,3)) >*/
	    montecarl_1.i__[2] = code_(channels_1.partic + (*jch + 1319) * 12,
		     (ftnlen)12);
/*<          END IF >*/
	}
/*<          GOTO 99 >*/
	goto L99;
/*<       END IF >*/
    }
/* ***    Canale 34 */
/*<       IF (JCH.EQ.34) THEN >*/
    if (*jch == 34) {
/*<          NUMPAR=NPLEV(JCH,1) >*/
	numpar = channels_1.nplev[*jch - 1];
/*<          NPART = NUMPAR >*/
	montecarl_1.npart = numpar;
/*<          P(1,1) = PLAB1(1,1) >*/
	montecarl_1.p[0] = plab1[0];
/*<          P(2,1) = PLAB1(2,1) >*/
	montecarl_1.p[1] = plab1[1];
/*<          P(3,1) = PLAB1(3,1) >*/
	montecarl_1.p[2] = plab1[2];
/*<          P(4,1) = AMASS1(2) >*/
	montecarl_1.p[3] = amass1[1];
/*<          I(1) = CODE(PARTIC(JCH,1,0,1)) >*/
	montecarl_1.i__[0] = code_(channels_1.partic + (*jch + 39) * 12, (
		ftnlen)12);
/*<          P(1,2) = PLAB1(1,2) >*/
	montecarl_1.p[4] = plab1[5];
/*<          P(2,2) = PLAB1(2,2) >*/
	montecarl_1.p[5] = plab1[6];
/*<          P(3,2) = PLAB1(3,2) >*/
	montecarl_1.p[6] = plab1[7];
/*<          P(4,2) = AMASS1(2) >*/
	montecarl_1.p[7] = amass1[1];
/*<          I(2) = CODE(PARTIC(JCH,1,0,2)) >*/
	montecarl_1.i__[1] = code_(channels_1.partic + (*jch + 679) * 12, (
		ftnlen)12);
/*<          IF (TARGET.EQ.'He3') THEN >*/
	if (s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<             P(1,3) = PN(1,2) >*/
	    montecarl_1.p[8] = mom_nucleon__1.pn[5];
/*<             P(2,3) = PN(2,2) >*/
	    montecarl_1.p[9] = mom_nucleon__1.pn[6];
/*<             P(3,3) = PN(3,2) >*/
	    montecarl_1.p[10] = mom_nucleon__1.pn[7];
/*<             P(4,3) = RMASS(PARTIC(JCH,1,0,3)) >*/
	    montecarl_1.p[11] = rmass_(channels_1.partic + (*jch + 1319) * 12,
		     (ftnlen)12);
/*<             I(3) = CODE(PARTIC(JCH,1,0,3)) >*/
	    montecarl_1.i__[2] = code_(channels_1.partic + (*jch + 1319) * 12,
		     (ftnlen)12);
/*<          END IF >*/
	}
/*<          GOTO 99 >*/
	goto L99;
/*<       END IF >*/
    }
/*<       IF (JCH.eq.35) THEN >*/
    if (*jch == 35) {
/*<          NUMPAR=NPLEV(JCH,1)+NPLEV(JCH,2)-1          >*/
	numpar = channels_1.nplev[*jch - 1] + channels_1.nplev[*jch + 39] - 1;
/*<          NPART = NUMPAR >*/
	montecarl_1.npart = numpar;
/* *** non va aggiunto NUMPLUS perche' alla fine c'e' D e non p-n */
/* *** pione di decadimento della ro */
/*<          P(1,1) = PLAB2(1,1) >*/
	montecarl_1.p[0] = plab2[0];
/*<          P(2,1) = PLAB2(2,1) >*/
	montecarl_1.p[1] = plab2[1];
/*<          P(3,1) = PLAB2(3,1) >*/
	montecarl_1.p[2] = plab2[2];
/*<          P(4,1) = AMASS2(1) >*/
	montecarl_1.p[3] = amass2[0];
/*<          I(1) = CODE(PARTIC(JCH,2,JBR,1)) >*/
	montecarl_1.i__[0] = code_(channels_1.partic + (*jch + ((jbr + 4 << 2)
		 + 2) * 40 - 641) * 12, (ftnlen)12);
/* *** pione di decadimento della ro */
/*<          P(1,2) = PLAB2(1,2) >*/
	montecarl_1.p[4] = plab2[5];
/*<          P(2,2) = PLAB2(2,2) >*/
	montecarl_1.p[5] = plab2[6];
/*<          P(3,2) = PLAB2(3,2) >*/
	montecarl_1.p[6] = plab2[7];
/*<          P(4,2) = AMASS2(2) >*/
	montecarl_1.p[7] = amass2[1];
/*<          I(2) = CODE(PARTIC(JCH,2,JBR,2)) >*/
	montecarl_1.i__[1] = code_(channels_1.partic + (*jch + ((jbr + 8 << 2)
		 + 2) * 40 - 641) * 12, (ftnlen)12);
/* *** deutone diffuso */
/*<          P(1,3) = PLAB1(1,2) >*/
	montecarl_1.p[8] = plab1[5];
/*<          P(2,3) = PLAB1(2,2) >*/
	montecarl_1.p[9] = plab1[6];
/*<          P(3,3) = PLAB1(3,2) >*/
	montecarl_1.p[10] = plab1[7];
/*<          P(4,3) = AMASS1(2) >*/
	montecarl_1.p[11] = amass1[1];
/*<          I(3) = CODE(PARTIC(JCH,1,0,2)) >*/
	montecarl_1.i__[2] = code_(channels_1.partic + (*jch + 679) * 12, (
		ftnlen)12);
/* *** va aggiunto l'altro nucleone per He3 */
/*<          IF (TARGET.EQ.'He3') THEN >*/
	if (s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<             P(1,4) = PN(1,2) >*/
	    montecarl_1.p[12] = mom_nucleon__1.pn[5];
/*<             P(2,4) = PN(2,2) >*/
	    montecarl_1.p[13] = mom_nucleon__1.pn[6];
/*<             P(3,4) = PN(3,2) >*/
	    montecarl_1.p[14] = mom_nucleon__1.pn[7];
/*<             P(4,4) = RMASS(PARTIC(JCH,1,0,3)) >*/
	    montecarl_1.p[15] = rmass_(channels_1.partic + (*jch + 1319) * 12,
		     (ftnlen)12);
/*<             I(4) = CODE(PARTIC(JCH,1,0,3)) >*/
	    montecarl_1.i__[3] = code_(channels_1.partic + (*jch + 1319) * 12,
		     (ftnlen)12);
/*<          END IF >*/
	}
/*<          GOTO 99 >*/
	goto L99;
/*<       END IF >*/
    }
/* ***    Canale 36 */
/*<       IF (JCH.EQ.36) THEN >*/
    if (*jch == 36) {
/*<          NUMPAR=NPLEV(JCH,1) >*/
	numpar = channels_1.nplev[*jch - 1];
/*<          NPART = NUMPAR >*/
	montecarl_1.npart = numpar;
/*<          P(1,1) = PLAB1(1,1) >*/
	montecarl_1.p[0] = plab1[0];
/*<          P(2,1) = PLAB1(2,1) >*/
	montecarl_1.p[1] = plab1[1];
/*<          P(3,1) = PLAB1(3,1) >*/
	montecarl_1.p[2] = plab1[2];
/*<          P(4,1) = AMASS1(1) >*/
	montecarl_1.p[3] = amass1[0];
/*<          I(1) = CODE(PARTIC(JCH,1,0,1)) >*/
	montecarl_1.i__[0] = code_(channels_1.partic + (*jch + 39) * 12, (
		ftnlen)12);
/*<          P(1,2) = PLAB1(1,2) >*/
	montecarl_1.p[4] = plab1[5];
/*<          P(2,2) = PLAB1(2,2) >*/
	montecarl_1.p[5] = plab1[6];
/*<          P(3,2) = PLAB1(3,2) >*/
	montecarl_1.p[6] = plab1[7];
/*<          P(4,2) = AMASS1(2) >*/
	montecarl_1.p[7] = amass1[1];
/*<          I(2) = CODE(PARTIC(JCH,1,0,2)) >*/
	montecarl_1.i__[1] = code_(channels_1.partic + (*jch + 679) * 12, (
		ftnlen)12);
/*<          IF (TARGET.EQ.'He3') THEN >*/
	if (s_cmp(target_1.target, "He3", (ftnlen)12, (ftnlen)3) == 0) {
/*<             P(1,3) = PN(1,2) >*/
	    montecarl_1.p[8] = mom_nucleon__1.pn[5];
/*<             P(2,3) = PN(2,2) >*/
	    montecarl_1.p[9] = mom_nucleon__1.pn[6];
/*<             P(3,3) = PN(3,2) >*/
	    montecarl_1.p[10] = mom_nucleon__1.pn[7];
/*<             P(4,3) = RMASS(PARTIC(JCH,1,0,3)) >*/
	    montecarl_1.p[11] = rmass_(channels_1.partic + (*jch + 1319) * 12,
		     (ftnlen)12);
/*<             I(3) = CODE(PARTIC(JCH,1,0,3)) >*/
	    montecarl_1.i__[2] = code_(channels_1.partic + (*jch + 1319) * 12,
		     (ftnlen)12);
/*<          END IF >*/
	}
/*<          GOTO 99 >*/
	goto L99;
/*<       END IF >*/
    }
/* ***    Canale 37 */
/*<        IF (JCH.EQ.37) THEN >*/
    if (*jch == 37) {
/*<         NUMPAR=NPLEV(JCH,1) >*/
	numpar = channels_1.nplev[*jch - 1];
/*<         NPART = NUMPAR+NUMPLUS >*/
	montecarl_1.npart = numpar + numplus;
/*<         P(1,1) = PLAB1(1,1) >*/
	montecarl_1.p[0] = plab1[0];
/*<         P(2,1) = PLAB1(2,1) >*/
	montecarl_1.p[1] = plab1[1];
/*<         P(3,1) = PLAB1(3,1) >*/
	montecarl_1.p[2] = plab1[2];
/*<         P(4,1) = AMASS1(1) >*/
	montecarl_1.p[3] = amass1[0];
/*<         I(1) = CODE(PARTIC(JCH,1,0,1)) >*/
	montecarl_1.i__[0] = code_(channels_1.partic + (*jch + 39) * 12, (
		ftnlen)12);
/*<         P(1,2) = PLAB1(1,2) >*/
	montecarl_1.p[4] = plab1[5];
/*<         P(2,2) = PLAB1(2,2) >*/
	montecarl_1.p[5] = plab1[6];
/*<         P(3,2) = PLAB1(3,2) >*/
	montecarl_1.p[6] = plab1[7];
/*<         P(4,2) = AMASS1(2) >*/
	montecarl_1.p[7] = amass1[1];
/*<         I(2) = CODE(PARTIC(JCH,1,0,2)) >*/
	montecarl_1.i__[1] = code_(channels_1.partic + (*jch + 679) * 12, (
		ftnlen)12);
/*<         P(1,3) = PN(1,2) >*/
	montecarl_1.p[8] = mom_nucleon__1.pn[5];
/*<         P(2,3) = PN(2,2) >*/
	montecarl_1.p[9] = mom_nucleon__1.pn[6];
/*<         P(3,3) = PN(3,2) >*/
	montecarl_1.p[10] = mom_nucleon__1.pn[7];
/*<         P(4,3) = RMASS(PARTIC(JCH,1,0,3))          >*/
	montecarl_1.p[11] = rmass_(channels_1.partic + (*jch + 1319) * 12, (
		ftnlen)12);
/*<         I(3) = CODE(PARTIC(JCH,1,0,3)) >*/
	montecarl_1.i__[2] = code_(channels_1.partic + (*jch + 1319) * 12, (
		ftnlen)12);
/*<         P(1,4) = PN(1,3) >*/
	montecarl_1.p[12] = mom_nucleon__1.pn[10];
/*<         P(2,4) = PN(2,3) >*/
	montecarl_1.p[13] = mom_nucleon__1.pn[11];
/*<         P(3,4) = PN(3,3) >*/
	montecarl_1.p[14] = mom_nucleon__1.pn[12];
/*<         P(4,4) = RMASS(PARTIC(JCH,1,0,4))          >*/
	montecarl_1.p[15] = rmass_(channels_1.partic + (*jch + 1959) * 12, (
		ftnlen)12);
/*<         I(4) = CODE(PARTIC(JCH,1,0,4)) >*/
	montecarl_1.i__[3] = code_(channels_1.partic + (*jch + 1959) * 12, (
		ftnlen)12);
/*<         GOTO 99 >*/
	goto L99;
/*<        END IF >*/
    }
/* ***    Canali 38,39 */
/*<        IF (JCH.EQ.38.OR.JCH.EQ.39) THEN >*/
    if (*jch == 38 || *jch == 39) {
/*<         NUMPAR=NPLEV(JCH,1) >*/
	numpar = channels_1.nplev[*jch - 1];
/*<         NPART = NUMPAR+1 >*/
	montecarl_1.npart = numpar + 1;
/*<         IF (JCH.EQ.39) THEN >*/
	if (*jch == 39) {
/*<           KK=2 >*/
	    kk = 2;
/*<         ELSE >*/
	} else {
/*<           KK=1 >*/
	    kk = 1;
/*<         END IF >*/
	}
/*<         P(1,1) = PLAB1(1,1) >*/
	montecarl_1.p[0] = plab1[0];
/*<         P(2,1) = PLAB1(2,1) >*/
	montecarl_1.p[1] = plab1[1];
/*<         P(3,1) = PLAB1(3,1) >*/
	montecarl_1.p[2] = plab1[2];
/*<         P(4,1) = AMASS1(1) >*/
	montecarl_1.p[3] = amass1[0];
/*<         I(1) = CODE(PARTIC(JCH,1,0,1)) >*/
	montecarl_1.i__[0] = code_(channels_1.partic + (*jch + 39) * 12, (
		ftnlen)12);
/*<         IF (JCH.EQ.39) THEN >*/
	if (*jch == 39) {
/*<           P(1,2) = PLAB1(1,2) >*/
	    montecarl_1.p[4] = plab1[5];
/*<           P(2,2) = PLAB1(2,2) >*/
	    montecarl_1.p[5] = plab1[6];
/*<           P(3,2) = PLAB1(3,2) >*/
	    montecarl_1.p[6] = plab1[7];
/*<           P(4,2) = AMASS1(2) >*/
	    montecarl_1.p[7] = amass1[1];
/*<           I(2) = CODE(PARTIC(JCH,1,0,2)) >*/
	    montecarl_1.i__[1] = code_(channels_1.partic + (*jch + 679) * 12, 
		    (ftnlen)12);
/*<           P(1,3) = PLAB1(1,3) >*/
	    montecarl_1.p[8] = plab1[10];
/*<           P(2,3) = PLAB1(2,3) >*/
	    montecarl_1.p[9] = plab1[11];
/*<           P(3,3) = PLAB1(3,3) >*/
	    montecarl_1.p[10] = plab1[12];
/*<           P(4,3) = AMASS1(3) >*/
	    montecarl_1.p[11] = amass1[2];
/*<           I(3) = CODE(PARTIC(JCH,1,0,3)) >*/
	    montecarl_1.i__[2] = code_(channels_1.partic + (*jch + 1319) * 12,
		     (ftnlen)12);
/*<           P(1,4) = PN(1,2) >*/
	    montecarl_1.p[12] = mom_nucleon__1.pn[5];
/*<           P(2,4) = PN(2,2) >*/
	    montecarl_1.p[13] = mom_nucleon__1.pn[6];
/*<           P(3,4) = PN(3,2) >*/
	    montecarl_1.p[14] = mom_nucleon__1.pn[7];
/*<           P(4,4) = RMASS(PARTIC(JCH,1,0,4))          >*/
	    montecarl_1.p[15] = rmass_(channels_1.partic + (*jch + 1959) * 12,
		     (ftnlen)12);
/*<           I(4) = CODE(PARTIC(JCH,1,0,4)) >*/
	    montecarl_1.i__[3] = code_(channels_1.partic + (*jch + 1959) * 12,
		     (ftnlen)12);
/*<         ELSE >*/
	} else {
/*<           P(1,2) = PLAB1(1,2) >*/
	    montecarl_1.p[4] = plab1[5];
/*<           P(2,2) = PLAB1(2,2) >*/
	    montecarl_1.p[5] = plab1[6];
/*<           P(3,2) = PLAB1(3,2) >*/
	    montecarl_1.p[6] = plab1[7];
/*<           P(4,2) = AMASS1(2) >*/
	    montecarl_1.p[7] = amass1[1];
/*<           I(2) = CODE(PARTIC(JCH,1,0,2)) >*/
	    montecarl_1.i__[1] = code_(channels_1.partic + (*jch + 679) * 12, 
		    (ftnlen)12);
/*<           P(1,3) = PN(1,2) >*/
	    montecarl_1.p[8] = mom_nucleon__1.pn[5];
/*<           P(2,3) = PN(2,2) >*/
	    montecarl_1.p[9] = mom_nucleon__1.pn[6];
/*<           P(3,3) = PN(3,2) >*/
	    montecarl_1.p[10] = mom_nucleon__1.pn[7];
/*<           P(4,3) = RMASS(PARTIC(JCH,1,0,3))          >*/
	    montecarl_1.p[11] = rmass_(channels_1.partic + (*jch + 1319) * 12,
		     (ftnlen)12);
/*<           I(3) = CODE(PARTIC(JCH,1,0,3)) >*/
	    montecarl_1.i__[2] = code_(channels_1.partic + (*jch + 1319) * 12,
		     (ftnlen)12);
/*<         END IF >*/
	}
/*<         GOTO 99 >*/
	goto L99;
/*<        END IF >*/
    }
/* ***    Canale 40 */
/*<        IF (JCH.EQ.40) THEN >*/
    if (*jch == 40) {
/*<         NPART = 4 >*/
	montecarl_1.npart = 4;
/*<         P(1,1) = PLAB1(1,1) >*/
	montecarl_1.p[0] = plab1[0];
/*<         P(2,1) = PLAB1(2,1) >*/
	montecarl_1.p[1] = plab1[1];
/*<         P(3,1) = PLAB1(3,1) >*/
	montecarl_1.p[2] = plab1[2];
/*<         P(4,1) = AMASS1(1) >*/
	montecarl_1.p[3] = amass1[0];
/*<         I(1) = CODE(PARTIC(JCH,1,0,1)) >*/
	montecarl_1.i__[0] = code_(channels_1.partic + (*jch + 39) * 12, (
		ftnlen)12);
/*<         P(1,2) = PLAB2(1,1) >*/
	montecarl_1.p[4] = plab2[0];
/*<         P(2,2) = PLAB2(2,1) >*/
	montecarl_1.p[5] = plab2[1];
/*<         P(3,2) = PLAB2(3,1) >*/
	montecarl_1.p[6] = plab2[2];
/*<         P(4,2) = AMASS2(1) >*/
	montecarl_1.p[7] = amass2[0];
/*<         I(2) = CODE(PARTIC(JCH,2,1,1)) >*/
	montecarl_1.i__[1] = code_(channels_1.partic + (*jch + 239) * 12, (
		ftnlen)12);
/*<         P(1,3) = PLAB2(1,2) >*/
	montecarl_1.p[8] = plab2[5];
/*<         P(2,3) = PLAB2(2,2) >*/
	montecarl_1.p[9] = plab2[6];
/*<         P(3,3) = PLAB2(3,2) >*/
	montecarl_1.p[10] = plab2[7];
/*<         P(4,3) = AMASS2(2) >*/
	montecarl_1.p[11] = amass2[1];
/*<         I(3) = CODE(PARTIC(JCH,2,1,2)) >*/
	montecarl_1.i__[2] = code_(channels_1.partic + (*jch + 879) * 12, (
		ftnlen)12);
/*<         P(1,4) = PN(1,2) >*/
	montecarl_1.p[12] = mom_nucleon__1.pn[5];
/*<         P(2,4) = PN(2,2) >*/
	montecarl_1.p[13] = mom_nucleon__1.pn[6];
/*<         P(3,4) = PN(3,2) >*/
	montecarl_1.p[14] = mom_nucleon__1.pn[7];
/*<         P(4,4) = RMASS(PARTIC(JCH,1,0,3))          >*/
	montecarl_1.p[15] = rmass_(channels_1.partic + (*jch + 1319) * 12, (
		ftnlen)12);
/*<         I(4) = CODE(PARTIC(JCH,1,0,3)) >*/
	montecarl_1.i__[3] = code_(channels_1.partic + (*jch + 1319) * 12, (
		ftnlen)12);
/*<        END IF >*/
    }
/*<  99   CONTINUE >*/
L99:
/*      CALL HFILL(JCH,EGAM,0.,1.) */
/*<       RETURN >*/
    return 0;
/*<       END                                                !*** END EV_NUC >*/
} /* ev_nuc__ */

#undef plab1
#undef jch


/*< 	SUBROUTINE GLOREN(BETA,PA,PB) >*/
/* Subroutine */ int gloren_(real *beta, real *pa, real *pb)
{
    static real bpgam, betpa;

/* . */
/* .    ****************************************************************** */
/* .    *                                                                * */
/*     *       Routine to transform momentum and energy from the        * */
/*     *       Lorentz frame A to the Lorentz frame B                   * */
/*     *                                                                * */
/*     *       PA(1)                                                    * */
/*     *       PA(2)     Momentum components in frame A                 * */
/*     *       PA(3)                                                    * */
/*     *       PA(4)     Energy                                         * */
/*     *       PB(..)   same quantities in frame B                      * */
/*     *                                                                * */
/*     *       BETA(1)    Components of velocity of frame B             * */
/*     *       BETA(2)        as seen from frame A                      * */
/*     *       BETA(3)                                                  * */
/*     *       BETA(4)    1./SQRT(1.-BETA**2)                           * */
/* .    *                                                                * */
/* .    *       Author    M.Hansroul  *********                          * */
/* .    *                                                                * */
/* .    ****************************************************************** */
/* 	RIA : taken from geant321 */
/* . */
/*< 	DIMENSION BETA(4),PA(4),PB(4) >*/
/* . */
/* .    ------------------------------------------------------------------ */
/* . */
/*< 	BETPA  = BETA(1)*PA(1) + BETA(2)*PA(2) + BETA(3)*PA(3) >*/
    /* Parameter adjustments */
    --pb;
    --pa;
    --beta;

    /* Function Body */
    betpa = beta[1] * pa[1] + beta[2] * pa[2] + beta[3] * pa[3];
/*< 	BPGAM  = (BETPA * BETA(4)/(BETA(4) + 1.) - PA(4)) * BETA(4) >*/
    bpgam = (betpa * beta[4] / (beta[4] + 1.f) - pa[4]) * beta[4];
/*< 	PB(1) = PA(1) + BPGAM  * BETA(1) >*/
    pb[1] = pa[1] + bpgam * beta[1];
/*< 	PB(2) = PA(2) + BPGAM  * BETA(2) >*/
    pb[2] = pa[2] + bpgam * beta[2];
/*< 	PB(3) = PA(3) + BPGAM  * BETA(3) >*/
    pb[3] = pa[3] + bpgam * beta[3];
/*< 	PB(4) =(PA(4) - BETPA) * BETA(4) >*/
    pb[4] = (pa[4] - betpa) * beta[4];
/*< 	END >*/
    return 0;
} /* gloren_ */

