/* yn5i.F -- translated by f2c (version 20230428).
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

struct b_1_ {
    real b1[30], b2[30], ctm[30];
};

#define b_1 (*(struct b_1_ *) &b_)

struct tl_1_ {
    real tl[30];
    integer nl;
};

#define tl_1 (*(struct tl_1_ *) &tl_)

struct te_1_ {
    real te[30];
    integer ne;
};

#define te_1 (*(struct te_1_ *) &te_)

struct be_1_ {
    real be[30];
};

#define be_1 (*(struct be_1_ *) &be_)

struct coefa_1_ {
    real a1[72], a2[72], a3[72], a4[72], a5[72], a6[72], a7[32];
};

#define coefa_1 (*(struct coefa_1_ *) &coefa_)

struct coefbc_1_ {
    real b1[64], b2[64], c__[24];
};
struct coefbc_2_ {
    real bnkj[128]	/* was [4][4][8] */, ckj[24]	/* was [3][8] */;
};

#define coefbc_1 (*(struct coefbc_1_ *) &coefbc_)
#define coefbc_2 (*(struct coefbc_2_ *) &coefbc_)

struct {
    real w[8], fiks[8];
} cohelp_;

#define cohelp_1 cohelp_

/* Initialized data */

struct {
    real e_1[152];
    } coefbc_ = { .50278f, 3.1442f, -7.8172f, 8.1667f, .93482f, -10.59f, 
	    29.227f, -34.55f, -.096685f, 4.7335f, -14.298f, 17.685f, 
	    -.025041f, -.62478f, 2.0282f, -2.5895f, 1.1965f, -.82889f, 
	    1.0426f, -1.909f, .28703f, -4.9065f, 16.264f, -19.904f, -.24492f, 
	    2.9191f, -9.5776f, 11.938f, .037297f, -.422f, 1.3883f, -1.7476f, 
	    1.3508f, -4.3139f, 12.291f, -15.288f, -.20086f, 1.3641f, -3.403f, 
	    3.8559f, .012583f, -.083492f, .186f, -.20043f, -2.3628e-4f, 
	    .0013514f, -.0024324f, .0021906f, 1.2419f, -4.3633f, 13.743f, 
	    -18.592f, -.24404f, 1.3158f, -3.5691f, 4.3867f, .015693f, 
	    -.082579f, .21427f, -.25846f, -2.9386e-4f, .001406f, -.0033835f, 
	    .0038664f, .63054f, -3.7333f, 13.464f, -18.594f, 2.1801f, 1.5163f,
	     -16.38f, 27.944f, -1.2886f, -2.457f, 15.129f, -23.295f, .20915f, 
	    .52279f, -2.8687f, 4.2688f, .93363f, -1.8181f, 5.5157f, -8.5216f, 
	    1.7811f, -8.2927f, 20.607f, -20.827f, -1.5264f, 6.8433f, -16.067f,
	     16.845f, .27128f, -1.1944f, 2.7495f, -2.9045f, 1.9439f, -4.6268f,
	     9.7879f, -9.6074f, -.3464f, 1.1093f, -1.9313f, 1.7064f, .027054f,
	     -.11638f, .26969f, -.31853f, -6.6092e-4f, .0050728f, -.014995f, 
	    .019605f, 1.8693f, -5.5678f, 14.795f, -16.903f, -.49965f, 1.7874f,
	     -4.133f, 3.8393f, .046194f, -.18536f, .45315f, -.46273f, 
	    -.0013341f, .005771f, -.014554f, .015554f, .14509f, .4652f, 
	    -.033005f, .15376f, .27436f, -.014604f, .62959f, .17866f, 
	    -.0026216f, .8381f, .0086137f, .0032946f, .092852f, .53886f, 
	    -.054493f, .13032f, .40709f, -.028782f, .14909f, .38502f, 
	    -.012775f, .18024f, .33022f, -.0094491f };

struct {
    real e_1[90];
    } b_ = { 1.13f, 1.52f, 1.67f, 1.94f, 2.45f, 2.7f, 3.f, 3.54f, 4.14f, 
	    4.46f, 4.71f, 5.99f, 7.04f, 7.58f, 8.57f, 9.37f, 10.06f, 10.87f, 
	    11.12f, 12.27f, 13.55f, 14.73f, 15.32f, 15.09f, 15.33f, 17.5f, 
	    18.84f, 29.66f, 53.06f, 84.34f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
	     0.f, .18f, .21f, .23f, .3f, .35f, .38f, .429f, .468f, .503f, 
	    .544f, .56f, .61f, .677f, .736f, .77f, .8f, .851f, .875f, 1.05f, 
	    3.707f, 5.97f, 7.028f, -1.f, -1.f, -1.f, -.9f, -.8f, -.68f, -.53f,
	     -.37f, -.19f, -.02f, .11f, .21f, .27f, .3f, .34f, .4f, .44f, 
	    .48f, .5f, .54f, .57f, .6f, .61f, .63f, .65f, .67f, .72f, .75f, 
	    .906f, .913f };

struct {
    real e_1[30];
    integer e_2;
    } tl_ = { .02f, .04f, .05f, .063f, .083f, .11f, .137f, .175f, .226f, 
	    .288f, .334f, .425f, .499f, .538f, .608f, .664f, .713f, .771f, 
	    .83f, .916f, 1.011f, 1.1f, 1.144f, 1.189f, 1.279f, 1.379f, 1.572f,
	     2.203f, 4.84f, 7.115f, 30 };

struct {
    real e_1[30];
    integer e_2;
    } te_ = { .0437f, .116f, .143f, .167f, .192f, .214f, .231f, .249f, 2.68f, 
	    4.15f, 5.13f, 6.12f, 8.11f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 
	    0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 13 };

struct {
    real e_1[30];
    } be_ = { .838f, 1.284f, 1.695f, 1.575f, 1.989f, 2.172f, 2.388f, 1.908f, 
	    16.75f, 18.92f, 22.87f, 23.88f, 30.97f, 0.f, 0.f, 0.f, 0.f, 0.f, 
	    0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f };

struct {
    real e_1[464];
    } coefa_ = { 2.7404f, -9.6998f, 10.4f, 2.3882f, -7.5137f, 44.096f, 
	    -74.379f, 46.038f, 7.5479f, -39.274f, 64.835f, -41.609f, -1.8369f,
	     8.6911f, -13.06f, 7.188f, -30.853f, 106.24f, -129.39f, 54.339f, 
	    19.465f, -68.102f, 96.358f, -56.827f, -3.4831f, 12.341f, -18.592f,
	     12.024f, .18941f, -.6788f, 1.0665f, -.7291f, .10258f, -1.0542f, 
	    11.389f, -16.638f, -.49607f, 11.8f, -90.857f, 164.76f, 1.5437f, 
	    -33.769f, 251.92f, -450.71f, -1.2021f, 25.336f, -186.58f, 332.54f,
	     .15789f, 2.9671f, -5.5251f, 6.8925f, -7.0218f, -205.34f, 569.51f,
	     -898.58f, 134.96f, 4872.2f, -14674.f, 23924.f, -821.16f, 
	    -32586.f, 100980.f, -165530.f, .31531f, -7.4981f, 43.295f, 
	    -76.36f, -6.5373f, 193.07f, -1018.1f, 1742.6f, 46.864f, -1303.f, 
	    6729.1f, -11075.f, -95.192f, 2637.3f, -12857.f, 20294.f, -17.953f,
	     109.72f, -239.54f, 228.26f, 91.968f, -519.63f, 1126.6f, -1074.f, 
	    -132.7f, 741.12f, -1600.f, 1524.9f, 58.598f, -318.74f, 677.51f, 
	    -640.11f, .42169f, 147.05f, -653.35f, 915.07f, -3.5198f, -260.19f,
	     1225.f, -1748.1f, 3.6373f, 155.92f, -752.01f, 1079.6f, -.78041f, 
	    -30.563f, 147.95f, -212.5f, -.38288f, 3.7587f, -6.5144f, 6.774f, 
	    103.81f, -272.82f, 477.59f, -512.22f, -1788.2f, 4305.2f, -7931.4f,
	     9347.1f, 7147.5f, -3339.5f, -4139.2f, -4436.4f, .24991f, 32.028f,
	     -118.82f, 150.99f, -2.6994f, -460.45f, 1895.9f, -2519.f, 16.268f,
	     2138.4f, -9126.2f, 12431.f, -29.654f, -3182.3f, 13944.f, 
	    -19342.f, 3.9025f, -91.126f, 323.73f, -400.48f, -20.619f, 491.7f, 
	    -1715.5f, 2114.3f, 33.004f, -766.84f, 2700.3f, -3352.5f, -16.367f,
	     373.94f, -1320.2f, 1642.3f, 19.402f, -224.46f, 747.33f, -935.7f, 
	    -44.18f, 471.94f, -1485.6f, 1805.5f, 31.567f, -301.76f, 907.63f, 
	    -1077.3f, -6.8648f, 60.476f, -175.2f, 203.81f, .40693f, -4.1404f, 
	    14.044f, -17.265f, -3.6799f, 59.61f, -162.69f, 188.73f, 14.556f, 
	    -175.5f, 458.39f, -533.9f, -12.621f, 149.64f, -381.18f, 451.41f, 
	    -.47554f, 2.2641f, -12.528f, 24.647f, 5.162f, -9.9236f, 55.623f, 
	    -104.62f, -8.1117f, 19.315f, -84.255f, 139.08f, 3.5187f, -9.1783f,
	     34.95f, -51.243f, .48173f, 5.7726f, -13.745f, 27.125f, -4.4804f, 
	    -38.582f, 111.59f, -243.05f, 16.306f, 110.46f, -330.45f, 722.7f, 
	    -15.968f, -80.14f, 246.16f, -607.53f, -5.1646f, -6.0776f, 78.989f,
	     -107.05f, 21.871f, 56.915f, -401.59f, 512.15f, -27.993f, -94.67f,
	     569.28f, -696.21f, 11.587f, 45.998f, -245.6f, 284.52f, -53.067f, 
	    576.12f, -1543.8f, 164550.f, 147.5f, -1638.f, 4592.3f, -4994.9f, 
	    -134.36f, 1578.f, -4446.3f, 4902.2f, 40.253f, -488.6f, 1400.1f, 
	    -1560.6f, .14988f, 2.8753f, -5.3078f, 6.2233f, -5.9558f, -162.03f,
	     430.79f, -625.48f, 128.75f, 3140.2f, -7918.9f, 10983.f, -851.61f,
	     -18780.f, 44607.f, -58790.f, .53689f, -13.216f, 81.011f, 
	    -142.85f, -10.55f, 296.29f, -1695.7f, 2893.5f, 69.621f, -1924.5f, 
	    10620.f, -17468.f, -138.65f, 3928.1f, -20293.f, 32058.f, .65288f, 
	    .38977f, .84078f, .18893f, -4.3964f, 34.309f, -73.692f, 84.308f, 
	    14.889f, -143.8f, 312.27f, -350.14f, -15.658f, 171.6f, -372.12f, 
	    412.99f, .085591f, 5.039f, -13.782f, 14.661f, .054284f, -9.2324f, 
	    36.397f, -42.962f, -.051111f, 4.6003f, -20.534f, 27.731f, 
	    .0074514f, -.62529f, 2.9159f, -4.1101f, .071622f, 3.096f, 
	    -11.125f, 18.13f, .092581f, -3.2186f, 20.273f, -33.245f, 
	    -.051531f, .89886f, -7.5084f, 13.188f, .0058258f, -.0017288f, 
	    .70224f, -1.4854f, .0823f, .15854f, 3.7716f, -4.0562f, .010802f, 
	    -.33688f, 1.1727f, -.67476f, -.0021798f, .052166f, -.25816f, 
	    .32048f, 6.5764e-5f, -.0014711f, .0078209f, -.01058f, .11138f, 
	    .60396f, 3.0174f, -4.419f, -.017709f, .23015f, -1.8187f, 3.4518f, 
	    .0020977f, -.025458f, .21626f, -.40692f, -5.4799e-5f, 5.9111e-4f, 
	    -.0055552f, .010647f, .17288f, 7.108f, -17.961f, 16.403f, 
	    -.14504f, -13.032f, 41.781f, -40.799f, .04539f, 8.3515f, -30.26f, 
	    32.882f, -.0047961f, -1.4095f, 5.3505f, -6.0946f, .037596f, 
	    1.4331f, -3.135f, 6.4864f, .23827f, 1.8253f, 1.7648f, -16.735f, 
	    -.1541f, -1.5201f, -1.5692f, 17.185f, .025037f, .30588f, .3252f, 
	    -3.5277f, .12489f, 1.3573f, .82338f, -1.4595f, -.051577f, 
	    -.35778f, -1.169f, 1.8078f, .0074864f, .032888f, .23744f, 
	    -.39802f, -2.988e-4f, -7.5117e-4f, -.011402f, .019505f, .1847f, 
	    1.9269f, -3.2979f, 3.6843f, -.073932f, .27213f, 1.06f, -2.3354f, 
	    .018907f, -.056473f, -.16487f, .38426f, -9.2984e-4f, .0025506f, 
	    .0073052f, -.01722f, -1.0306f, 32.849f, -75.052f, 60.255f, 
	    7.9586f, -125.72f, 256.04f, -165.47f, -14.797f, 165.9f, -279.91f, 
	    113.33f, 8.2309f, -67.871f, 85.762f, 5.9727f, -237.22f, 968.9f, 
	    -1621.9f, 1363.7f, 658.f, -2694.1f, 4548.f, -3846.f, -606.53f, 
	    2498.3f, -4249.8f, 3613.6f, 186.04f, -769.33f, 1316.6f, -1124.2f }
	    ;


/* Table of constant values */

static integer c__1 = 1;

/*<         FUNCTION QB(X,T,B,N) >*/
doublereal qb_(real *x, real *t, real *b, integer *n)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3;

    /* Local variables */
    static real a, c__, d__;
    static integer lp;
    static real ph1, ph2, ph3, ps1, ps2, ps3, del, dela, delc, deld;

/*<         DIMENSION B(30),T(30) >*/
/*<         LP=1 >*/
    /* Parameter adjustments */
    --b;
    --t;

    /* Function Body */
    lp = 1;
/*<   18    IF(X-T(LP) < 0) THEN >*/
L18:
    if (*x - t[lp] < 0.f) {
/*<         GOTO 10 >*/
	goto L10;
/*<         ELSE IF (X-T(LP) == 0) THEN >*/
    } else if (*x - t[lp] == 0.f) {
/*<         GOTO 11 >*/
	goto L11;
/*<         ELSE IF (X-T(LP) > 0) THEN >*/
    } else if (*x - t[lp] > 0.f) {
/*<         GOTO 12 >*/
	goto L12;
/*<         END IF >*/
    }
/*<   10    IF(LP-1 < 0) THEN >*/
L10:
    if (lp - 1 < 0) {
/*<         GOTO 13 >*/
	goto L13;
/*<         ELSE IF (LP-1 == 0) THEN >*/
    } else if (lp - 1 == 0) {
/*<         GOTO 13 >*/
	goto L13;
/*<         ELSE IF (LP-1 > 0) THEN >*/
    } else if (lp - 1 > 0) {
/*<         GOTO 14 >*/
	goto L14;
/*<         END IF >*/
    }
/*<  13    QB=B(1) >*/
L13:
    ret_val = b[1];
/*<         RETURN >*/
    return ret_val;
/*<   14    IF(LP-(N-1) < 0) THEN >*/
L14:
    if (lp - (*n - 1) < 0) {
/*<         GOTO 17 >*/
	goto L17;
/*<         ELSE IF (LP-(N-1) == 0) THEN >*/
    } else if (lp - (*n - 1) == 0) {
/*<         GOTO 16 >*/
	goto L16;
/*<         ELSE IF (LP-(N-1) > 0) THEN >*/
    } else if (lp - (*n - 1) > 0) {
/*<         GOTO 16 >*/
	goto L16;
/*<         END IF >*/
    }
/*<   16    PH1=B(N-2) >*/
L16:
    ph1 = b[*n - 2];
/*<         PS1=T(N-2) >*/
    ps1 = t[*n - 2];
/*<         PH2=B(N-1) >*/
    ph2 = b[*n - 1];
/*<         PS2=T(N-1) >*/
    ps2 = t[*n - 1];
/*<         PH3=B(N) >*/
    ph3 = b[*n];
/*<         PS3=T(N) >*/
    ps3 = t[*n];
/*<         GO TO 15 >*/
    goto L15;
/*<   11    QB=B(LP) >*/
L11:
    ret_val = b[lp];
/*<         RETURN >*/
    return ret_val;
/*<   12    LP=LP+1 >*/
L12:
    ++lp;
/*<         IF(LP-N < 0) THEN >*/
    if (lp - *n < 0) {
/*<         GOTO 18 >*/
	goto L18;
/*<         ELSE IF ((LP-N) == 0) THEN >*/
    } else if (lp - *n == 0) {
/*<         GOTO 19 >*/
	goto L19;
/*<         ELSE IF ((LP-N) > 0) THEN >*/
    } else if (lp - *n > 0) {
/*<         GOTO 19 >*/
	goto L19;
/*<         END IF >*/
    }
/*<   19    QB=B(N) >*/
L19:
    ret_val = b[*n];
/*<         RETURN >*/
    return ret_val;
/*<   17    PH1=B(LP-1) >*/
L17:
    ph1 = b[lp - 1];
/*<         PS1=T(LP-1) >*/
    ps1 = t[lp - 1];
/*<         PH2=B(LP) >*/
    ph2 = b[lp];
/*<         PS2=T(LP) >*/
    ps2 = t[lp];
/*<         PH3=B(LP+1) >*/
    ph3 = b[lp + 1];
/*<         PS3=T(LP+1) >*/
    ps3 = t[lp + 1];
/*<         GO TO 15 >*/
    goto L15;
/*<   15    DEL=(PS2-PS3)*PS1**2+(PS3-PS1)*PS2**2+(PS1-PS2)*PS3**2 >*/
L15:
/* Computing 2nd power */
    r__1 = ps1;
/* Computing 2nd power */
    r__2 = ps2;
/* Computing 2nd power */
    r__3 = ps3;
    del = (ps2 - ps3) * (r__1 * r__1) + (ps3 - ps1) * (r__2 * r__2) + (ps1 - 
	    ps2) * (r__3 * r__3);
/*<         DELA=PH1*(PS2-PS3)+PH2*(PS3-PS1)+PH3*(PS1-PS2) >*/
    dela = ph1 * (ps2 - ps3) + ph2 * (ps3 - ps1) + ph3 * (ps1 - ps2);
/*<         DELD=(PH2-PH3)*PS1**2+(PH3-PH1)*PS2**2+(PH1-PH2)*PS3**2 >*/
/* Computing 2nd power */
    r__1 = ps1;
/* Computing 2nd power */
    r__2 = ps2;
/* Computing 2nd power */
    r__3 = ps3;
    deld = (ph2 - ph3) * (r__1 * r__1) + (ph3 - ph1) * (r__2 * r__2) + (ph1 - 
	    ph2) * (r__3 * r__3);
/*<        >*/
/* Computing 2nd power */
    r__1 = ps1;
/* Computing 2nd power */
    r__2 = ps2;
/* Computing 2nd power */
    r__3 = ps3;
    delc = (ps2 * ph3 - ps3 * ph2) * (r__1 * r__1) + (ps3 * ph1 - ps1 * ph3) *
	     (r__2 * r__2) + (ps1 * ph2 - ps2 * ph1) * (r__3 * r__3);
/*<         A=DELA/DEL >*/
    a = dela / del;
/*<         D=DELD/DEL >*/
    d__ = deld / del;
/*<         C=DELC/DEL >*/
    c__ = delc / del;
/*<         QB=A*X**2+D*X+C >*/
/* Computing 2nd power */
    r__1 = *x;
    ret_val = a * (r__1 * r__1) + d__ * *x + c__;
/*<         RETURN >*/
    return ret_val;
/*<         END >*/
} /* qb_ */

/*<         SUBROUTINE TC8 >*/
/* Subroutine */ int tc8_(void)
{
/*<         COMMON/B/B1(30),B2(30),CTM(30) >*/
/*<         COMMON/TL/TL(30),NL >*/
/*<         COMMON/TE/TE(30),NE >*/
/*<         COMMON/BE/BE(30) >*/
/*<        >*/
/*<        >*/
/*<        >*/
/*<        >*/
/*<        >*/
/*<        >*/
/*<         DATA NL/30/ >*/
/*<         DATA NE/13/ >*/
/*<         RETURN >*/
    return 0;
/*<         END >*/
} /* tc8_ */

/*<       FUNCTION FACTOR(M) >*/
doublereal factor_(integer *m)
{
    /* System generated locals */
    integer i__1;
    real ret_val;

    /* Local variables */
    static real f;
    static integer k;

/*<       F=1.0 >*/
    f = 1.f;
/*<       IF(M-2 < 0) THEN >*/
    if (*m - 2 < 0) {
/*<       GOTO 1 >*/
	goto L1;
/*<       ELSE IF (M-2 == 0) THEN >*/
    } else if (*m - 2 == 0) {
/*<       GOTO 2 >*/
	goto L2;
/*<       ELSE IF (M-2 > 0) THEN >*/
    } else if (*m - 2 > 0) {
/*<       GOTO 2 >*/
	goto L2;
/*<       END IF >*/
    }
/*<     1 FACTOR=F >*/
L1:
    ret_val = f;
/*<       RETURN >*/
    return ret_val;
/*<     2 DO 3 K=1,M >*/
L2:
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
/*<       F=FLOAT(K)*F >*/
	f = (real) k * f;
/*<     3 CONTINUE >*/
/* L3: */
    }
/*<       FACTOR=F >*/
    ret_val = f;
/*<       RETURN >*/
    return ret_val;
/*<       END >*/
} /* factor_ */

/*<       FUNCTION BINOM(N,J) >*/
doublereal binom_(integer *n, integer *j)
{
    /* System generated locals */
    integer i__1;
    real ret_val;

    /* Local variables */
    extern doublereal factor_(integer *);

/*<       BINOM=FACTOR(N)/(FACTOR(J)*FACTOR(N-J)) >*/
    i__1 = *n - *j;
    ret_val = factor_(n) / (factor_(j) * factor_(&i__1));
/*<       RETURN >*/
    return ret_val;
/*<       END >*/
} /* binom_ */

/*<       SUBROUTINE TC5 >*/
/* Subroutine */ int tc5_(void)
{
/*<       COMMON/COEFA/A1(72),A2(72),A3(72),A4(72),A5(72),A6(72),A7(32) >*/
/*<       COMMON/COEFBC/B1(64),B2(64),C(24) >*/
/*<        >*/
/*<        >*/
/*<        >*/
/*<        >*/
/*<        >*/
/*<        >*/
/*<        >*/
/*<        >*/
/*<        >*/
/*<        >*/
/*<          RETURN >*/
    return 0;
/*<          END >*/
} /* tc5_ */

/*<       SUBROUTINE CINEMA(PSTAR,V,P,CT,ST,CFI,SFI,T,CM) >*/
/* Subroutine */ int cinema_(real *pstar, real *v, real *p, real *ct, real *
	st, real *cfi, real *sfi, real *t, real *cm)
{
    /* System generated locals */
    real r__1, r__2, r__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real v2, pm, spv, temp1, temp2, temp3, temp4, tstar, pmstar;


/*     KINEMATIC BLOCK. BACKWARD LORENTS TRANSFORMATION: */

/*     PSTAR ---> -V ---> P */

/*     N.B.: SPECIAL CASE V=0 SHOULD BE TREATED SEPARATELY ! */

/*<       DIMENSION PSTAR(3),V(3),P(3) >*/
/*<       SPV = PSTAR(1)*V(1)+PSTAR(2)*V(2)+PSTAR(3)*V(3) >*/
    /* Parameter adjustments */
    --p;
    --v;
    --pstar;

    /* Function Body */
    spv = pstar[1] * v[1] + pstar[2] * v[2] + pstar[3] * v[3];
/*<       V2 = V(1)**2+V(2)**2+V(3)**2 >*/
/* Computing 2nd power */
    r__1 = v[1];
/* Computing 2nd power */
    r__2 = v[2];
/* Computing 2nd power */
    r__3 = v[3];
    v2 = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
/*<       TEMP1 = SQRT(1.-V2) >*/
    temp1 = sqrt(1.f - v2);
/*<       TEMP2 = SPV*(1./TEMP1-1.)/V2 >*/
    temp2 = spv * (1.f / temp1 - 1.f) / v2;
/*<       PMSTAR = SQRT(PSTAR(1)**2+PSTAR(2)**2+PSTAR(3)**2) >*/
/* Computing 2nd power */
    r__1 = pstar[1];
/* Computing 2nd power */
    r__2 = pstar[2];
/* Computing 2nd power */
    r__3 = pstar[3];
    pmstar = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
/*<       TSTAR=SQRT(PMSTAR**2+CM**2)-CM >*/
/* Computing 2nd power */
    r__1 = pmstar;
/* Computing 2nd power */
    r__2 = *cm;
    tstar = sqrt(r__1 * r__1 + r__2 * r__2) - *cm;
/*<       P(1) = PSTAR(1)+V(1)*TEMP2+V(1)*(TSTAR+CM)/TEMP1 >*/
    p[1] = pstar[1] + v[1] * temp2 + v[1] * (tstar + *cm) / temp1;
/*<       P(2) = PSTAR(2)+V(2)*TEMP2+V(2)*(TSTAR+CM)/TEMP1 >*/
    p[2] = pstar[2] + v[2] * temp2 + v[2] * (tstar + *cm) / temp1;
/*<       P(3) = PSTAR(3)+V(3)*TEMP2+V(3)*(TSTAR+CM)/TEMP1 >*/
    p[3] = pstar[3] + v[3] * temp2 + v[3] * (tstar + *cm) / temp1;
/*<       PM = SQRT(P(1)**2+P(2)**2+P(3)**2) >*/
/* Computing 2nd power */
    r__1 = p[1];
/* Computing 2nd power */
    r__2 = p[2];
/* Computing 2nd power */
    r__3 = p[3];
    pm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
/*<       IF (PM.GT.1.E-30) THEN >*/
    if (pm > 1e-30f) {
/*<           CT = P(3)/PM >*/
	*ct = p[3] / pm;
/*<           TEMP4 = 1.-CT**2 >*/
/* Computing 2nd power */
	r__1 = *ct;
	temp4 = 1.f - r__1 * r__1;
/*<           IF(TEMP4.GT.0.0) THEN >*/
	if (temp4 > 0.f) {
/*<               ST = SQRT(TEMP4) >*/
	    *st = sqrt(temp4);
/*<               TEMP3 = PM*ST   >*/
	    temp3 = pm * *st;
/*<               CFI=P(1)/TEMP3 >*/
	    *cfi = p[1] / temp3;
/*<               SFI=P(2)/TEMP3 >*/
	    *sfi = p[2] / temp3;
/*<           ELSE  ! FORWARD DIRECTION >*/
	} else {
/*<               CT=1. >*/
	    *ct = 1.f;
/*<               CFI=1. >*/
	    *cfi = 1.f;
/*<               SFI=0. >*/
	    *sfi = 0.f;
/*<               ST=0. >*/
	    *st = 0.f;
/*<           ENDIF >*/
	}
/*<       ELSE    ! BECAUSE PM TOO LOW >*/
    } else {
/*<           CT=1. >*/
	*ct = 1.f;
/*<           CFI=1. >*/
	*cfi = 1.f;
/*<           SFI=0. >*/
	*sfi = 0.f;
/*<           ST=0. >*/
	*st = 0.f;
/*<       ENDIF >*/
    }
/*<       T = SQRT(PM**2+CM**2)-CM >*/
/* Computing 2nd power */
    r__1 = pm;
/* Computing 2nd power */
    r__2 = *cm;
    *t = sqrt(r__1 * r__1 + r__2 * r__2) - *cm;
/*<       RETURN >*/
    return 0;
/*<          END >*/
} /* cinema_ */

/*<       SUBROUTINE CMS(P,V,PSTART,T,CM) >*/
/* Subroutine */ int cms_(real *p, real *v, real *pstart, real *t, real *cm)
{
    /* System generated locals */
    real r__1, r__2, r__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real v2, sv, spv, temp1, temp2;

/*<       DIMENSION P(3),V(3),PSTART(3) >*/
/*<       V2=V(1)**2+V(2)**2+V(3)**2 >*/
    /* Parameter adjustments */
    --pstart;
    --v;
    --p;

    /* Function Body */
/* Computing 2nd power */
    r__1 = v[1];
/* Computing 2nd power */
    r__2 = v[2];
/* Computing 2nd power */
    r__3 = v[3];
    v2 = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
/*<       TEMP1=SQRT(1.-V2) >*/
    temp1 = sqrt(1.f - v2);
/*<       SPV=P(1)*V(1)+P(2)*V(2)+P(3)*V(3) >*/
    spv = p[1] * v[1] + p[2] * v[2] + p[3] * v[3];
/*<       TEMP2=SPV/V2*(1./TEMP1-1.) >*/
    temp2 = spv / v2 * (1.f / temp1 - 1.f);
/*<       SV=(T+CM)/TEMP1 >*/
    sv = (*t + *cm) / temp1;
/*<       PSTART(1)=P(1)+V(1)*TEMP2+V(1)*SV >*/
    pstart[1] = p[1] + v[1] * temp2 + v[1] * sv;
/*<       PSTART(2)=P(2)+V(2)*TEMP2+V(2)*SV >*/
    pstart[2] = p[2] + v[2] * temp2 + v[2] * sv;
/*<       PSTART(3)=P(3)+V(3)*TEMP2+V(3)*SV >*/
    pstart[3] = p[3] + v[3] * temp2 + v[3] * sv;
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* cms_ */

/*<       FUNCTION FINT1 (FNAME,FLIM1,FLIM2,BS,RS) >*/
doublereal fint1_(E_fp fname, real *flim1, real *flim2, real *bs, real *rs)
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
    static integer it;
    static real fy, temp;


/*  GAUSS QUADRATURE METHOD FOR FUNCTION FNAME IN */
/*  THE REGION [FLIM1,FLIM2] */

/*<       COMMON/COHELP/W(8),FIKS(8) >*/
/*<       TEMP = 0. >*/
    temp = 0.f;
/*<       DO 30 IT=1,8 >*/
    for (it = 1; it <= 8; ++it) {
/*<       FY = ((FLIM1-FLIM2)*FIKS(IT)+FLIM1+FLIM2)/2. >*/
	fy = ((*flim1 - *flim2) * cohelp_1.fiks[it - 1] + *flim1 + *flim2) / 
		2.f;
/*<       TEMP = TEMP+W(IT)*FNAME(FY,BS,RS) >*/
	temp += cohelp_1.w[it - 1] * (*fname)(&fy, bs, rs);
/*<    30 CONTINUE >*/
/* L30: */
    }
/*<       FINT1 = TEMP*(FLIM1-FLIM2)/2. >*/
    ret_val = temp * (*flim1 - *flim2) / 2.f;
/*<       RETURN >*/
    return ret_val;
/*<          END >*/
} /* fint1_ */

/*<       FUNCTION PMOM (J,T) >*/
doublereal pmom_(integer *j, real *t)
{
    /* System generated locals */
    integer i__1, i__2;
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double pow_ri(real *, integer *), sqrt(doublereal);

    /* Local variables */
    static integer k, n;
    static real r1, s1, s2, s3, bnk[16]	/* was [4][4] */, vran[1];
    extern /* Subroutine */ int ranmar_(real *, integer *);

/*     BLOCK OF CALCULATION OF SECONDARY PARTICLES MOMENTUM. */
/*<       COMMON/COEFBC/BNKJ,CKJ >*/
/*<       DIMENSION BNKJ(4,4,8),CKJ(3,8),BNK(4,4) >*/
/*<       real vran(1) >*/
/*<       DO 10 K=1,4 >*/
    for (k = 1; k <= 4; ++k) {
/*<       DO N=1,4 >*/
	for (n = 1; n <= 4; ++n) {
/*<       BNK(N,K) = BNKJ(N,K,J) >*/
	    bnk[n + (k << 2) - 5] = coefbc_2.bnkj[n + (k + (*j << 2) << 2) - 
		    21];
/*<       END DO >*/
	}
/*<    10 CONTINUE >*/
/* L10: */
    }
/*<       S1 = 0. >*/
    s1 = 0.f;
/*<       call ranmar(vran,1) >*/
    ranmar_(vran, &c__1);
/*<       R1 = vran(1) >*/
    r1 = vran[0];
/*<       S2 = 0. >*/
    s2 = 0.f;
/*<       S3 = 0. >*/
    s3 = 0.f;
/*<       IF((T.GT.0.0).AND.(R1.GT.0.0)) THEN >*/
    if (*t > 0.f && r1 > 0.f) {
/*<          DO 11 N=1,4 >*/
	for (n = 1; n <= 4; ++n) {
/*<          DO K=1,4 >*/
	    for (k = 1; k <= 4; ++k) {
/*<        S1 = S1+BNK(N,K)*(T**(K-1))*(R1**(N-1)) >*/
		i__1 = k - 1;
		i__2 = n - 1;
		s1 += bnk[n + (k << 2) - 5] * pow_ri(t, &i__1) * pow_ri(&r1, &
			i__2);
/*<          END DO >*/
	    }
/*<    11    CONTINUE >*/
/* L11: */
	}
/*<       ELSE >*/
    } else {
/*<        S1 = BNK(1,1)+BNK(2,1)*R1+BNK(3,1)*R1**2+BNK(4,1)*R1**3 >*/
/* Computing 2nd power */
	r__1 = r1;
/* Computing 3rd power */
	r__2 = r1;
	s1 = bnk[0] + bnk[1] * r1 + bnk[2] * (r__1 * r__1) + bnk[3] * (r__2 * 
		(r__2 * r__2));
/*<       ENDIF >*/
    }
/*<       IF(T.GT.0.0) THEN >*/
    if (*t > 0.f) {
/*<          DO 12 N=1,4 >*/
	for (n = 1; n <= 4; ++n) {
/*<          DO K=1,4 >*/
	    for (k = 1; k <= 4; ++k) {
/*<            S2 = S2+BNK(N,K)*T**(K-1) >*/
		i__1 = k - 1;
		s2 += bnk[n + (k << 2) - 5] * pow_ri(t, &i__1);
/*<          END DO >*/
	    }
/*<    12    CONTINUE >*/
/* L12: */
	}
/*<         DO 13 K=1,3 >*/
	for (k = 1; k <= 3; ++k) {
/*<           S3 = S3+CKJ(K,J)*T**(K-1) >*/
	    i__1 = k - 1;
	    s3 += coefbc_2.ckj[k + *j * 3 - 4] * pow_ri(t, &i__1);
/*<    13   CONTINUE >*/
/* L13: */
	}
/*<       ELSE >*/
    } else {
/*<            S2 = BNK(1,1)+BNK(2,1)+BNK(3,1)+BNK(4,1) >*/
	s2 = bnk[0] + bnk[1] + bnk[2] + bnk[3];
/*<            S3 = CKJ(1,J) >*/
	s3 = coefbc_2.ckj[*j * 3 - 3];
/*<       ENDIF >*/
    }
/*<       PMOM = S3*SQRT(R1)*(S1+(1.-S2)*R1**4) >*/
/* Computing 4th power */
    r__1 = r1, r__1 *= r__1;
    ret_val = s3 * sqrt(r1) * (s1 + (1.f - s2) * (r__1 * r__1));
/*<       RETURN >*/
    return ret_val;
/*<          END >*/
} /* pmom_ */

