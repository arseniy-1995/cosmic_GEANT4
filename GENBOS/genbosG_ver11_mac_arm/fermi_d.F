CDECK  ID>, FERMI_D.
*CMZ :          13/04/94  18.02.52  by  Thomas Russew
*-- Author :
C***********************************************************************

      SUBROUTINE FERMI_D(W_N,W_DEU)

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                                                            
C FUNCTIONAL DESCRIPTION:                                         
C 
C    This subroutine computes invariant mass taking into account Fermi
C    motion. It defines also tetramomenta and beta of nucleons in CM
C
C    The momentum distribution is calculated by Vistor, wich uses 
C    Paris or Bonn potential
C
C     MODIFICHE
C
C      ******* by marco mirazita
C     29.04.99
C --> Aggiunto il file di include PROCESS.INC con le variabili per il
C     conteggio dei canali di reazione selezionati.
C     
C     05.05.99
C --> Eliminato il conteggio dei canali di reazione selezionati, che ora e'
C     fatto una volta per tutte nella routine CARDS (che sta in IN.F).
C --> Eliminata la parte in cui, nel caso siano selezionati sia canali su p
C     che su n, veniva scelto il nucleone bersaglio perche' era sbagliata.
C     Adesso la massa del nucleone bersaglio e' sempre "RMASS('proton   ')"
C     e quella dello spettatore "RMASS('neutron   ')", tanto i loro valori 
C     sono definiti uguali.
C --> Spostata l'estrazione dei numeri random subito dopo gli azzeramenti
C     iniziali.
C
C COMMON BLOCKS:                                               
C                                                             
C     /MONTECARL/
C     /FLAG_COMP/
C     /TARGET/
C     /H2/
C     /FLAG_PROC/
C     /MASSES/
C                                                             
C-----------------------------------------------------------------------

      IMPLICIT NONE

C   #include "process.inc"
C    #include "nt_kine.inc"

C     include "process.inc"

      INCLUDE 'process.inc'
      INCLUDE 'nt_kine.inc'




      INTEGER IBEAM,IELET,NCHAIN,ICHAIN(40),ICIBLE,IPOL,ILAM
      REAL ETAG,FWHM
      COMMON/FLAG_COMP/IELET,IBEAM,ICIBLE,ETAG,FWHM,NCHAIN,ICHAIN,IPOL,
     &   ILAM
      REAL EGAM, P(4,11)
      INTEGER I(11),NPART
      COMMON/MONTECARL/ EGAM,I,NPART,P
      DOUBLE PRECISION PI,TWOPI,PIBY2,DEGRAD,RADDEG,CLIGHT,BIG,EMASS
      DOUBLE PRECISION EMMU,PMASS,AVO
*
      PARAMETER (PI=3.14159265358979324D0)
      PARAMETER (TWOPI=6.28318530717958648D0)
      PARAMETER (PIBY2=1.57079632679489662D0)
      PARAMETER (DEGRAD=0.0174532925199432958D0)
      PARAMETER (RADDEG=57.2957795130823209D0)
      PARAMETER (CLIGHT=29979245800.D0)
      PARAMETER (BIG=10000000000.D0)
      PARAMETER (EMASS=0.0005109990615D0)
      PARAMETER (EMMU=0.105658387D0)
      PARAMETER (PMASS=0.9382723128D0)
      PARAMETER (AVO=0.60221367D0)
*
      INTEGER ITARG,IIDH
      CHARACTER*12 TARGET,TARGFERMI
      COMMON /TARGET/ ITARG,IIDH,TARGET,TARGFERMI
      REAL PN_N_D(5,3),PN_D(5,3),BETA_N_D(4),BETA_D(4)
      COMMON /H2/ PN_N_D,PN_D,BETA_N_D,BETA_D
      REAL RMB1,RMB2,RMME,DELTAMASS,DIBMASS,DEUMASS,RHOMASS
      COMMON /MASSES/ RMB1,RMB2,RMME,DELTAMASS,DIBMASS,DEUMASS,RHOMASS
      INTEGER ICP,ICN,J,JJ,K,KK,IDEU,IDP
      REAL W_N,W_DEU,ECM,Q,A,B,ETOT,S,RMASS,PP(0:100),
     &           DISTR_P(0:100),DISTR_B(0:100),P1,RND(5)
c      REAL RM(2)

      real mass_n_d,rm_fin,e_sp_fin,delta_e,totmass,totmass1,totmass2
      integer kv

      DATA DISTR_P / 0.0,6.7953602E-03,4.6238571E-02,
     &          0.1236217,0.2234159,
     &           0.3280089, 0.4262296, 0.5131546, 0.5876556, 0.6504297,
     &           0.7028700, 0.7465120, 0.7827951, 0.8129796, 0.8381307,
     &           0.8591338, 0.8767187, 0.8914840, 0.9039201, 0.9144295,
     &           0.9233416, 0.9309281, 0.9374120, 0.9429782, 0.9477784,
     &           0.9519395, 0.9555655, 0.9587440, 0.9615465, 0.9640333,
     &           0.9662549, 0.9682527, 0.9700612, 0.9717091, 0.9732212,
     &           0.9746166, 0.9759121, 0.9771215, 0.9782565, 0.9793263,
     &           0.9803385, 0.9812997, 0.9822153, 0.9830891, 0.9839256,
     &           0.9847271, 0.9854963, 0.9862349, 0.9869449, 0.9876275,
     &           0.9882837, 0.9889148, 0.9895214, 0.9901043, 0.9906641,
     &           0.9912014, 0.9917168, 0.9922106, 0.9926836, 0.9931361,
     &           0.9935685, 0.9939811, 0.9943746, 0.9947497, 0.9951068,
     &           0.9954458, 0.9957677, 0.9960729, 0.9963619, 0.9966350,
     &           0.9968930, 0.9971365, 0.9973655, 0.9975811, 0.9977834,
     &           0.9979731, 0.9981509, 0.9983172, 0.9984723, 0.9986169,
     &           0.9987513, 0.9988762, 0.9989921, 0.9990993, 0.9991982,
     &           0.9992898, 0.9993736, 0.9994510, 0.9995219, 0.9995865,
     &           0.9996458, 0.9996991, 0.9997482, 0.9997926, 0.9998326,
     &           0.9998688, 0.9999010, 0.9999300, 0.9999564, 0.9999794,
     &           1.000000/
      DATA DISTR_B / 0.0,7.0403279E-03,4.7886062E-02,
     &           0.1279474,0.2310553,
     &           0.3389300, 0.4400139, 0.5292494, 0.6055161, 0.6695801,
     &           0.7229201, 0.7671533, 0.8037899, 0.8341466, 0.8593346,
     &       0.8802744, 0.8977222, 0.9122964, 0.9245030, 0.9347540,
     &           0.9433876, 0.9506804, 0.9568596, 0.9621121, 0.9665918,
     &           0.9704258, 0.9737193, 0.9765596, 0.9790186, 0.9811566,
     &       0.9830235, 0.9846612, 0.9861036, 0.9873803, 0.9885156,
     &           0.9895294, 0.9904380, 0.9912565, 0.9919961, 0.9926671,
     &       0.9932781, 0.9938356, 0.9943457, 0.9948140, 0.9952440,
     &           0.9956399, 0.9960047, 0.9963413, 0.9966521, 0.9969389,
     &       0.9972037, 0.9974484, 0.9976741, 0.9978822, 0.9980739,
     &           0.9982503, 0.9984124, 0.9985616, 0.9986980, 0.9988228,
     &       0.9989371, 0.9990408, 0.9991356, 0.9992215, 0.9992995,
     &           0.9993696, 0.9994332, 0.9994900, 0.9995410, 0.9995866,
     &       0.9996274, 0.9996639, 0.9996958, 0.9997242, 0.9997495,
     &           0.9997714, 0.9997910, 0.9998078, 0.9998228, 0.9998358,
     &       0.9998475, 0.9998578, 0.9998671, 0.9998750, 0.9998825,
     &           0.9998896, 0.9998960, 0.9999021, 0.9999081, 0.9999138,
     &       0.9999200, 0.9999259, 0.9999323, 0.9999395, 0.9999465,
     &           0.9999539, 0.9999623, 0.9999706, 0.9999802, 0.9999896,
     &       1.000000/
      DATA IDP /1/         !Paris potential   (IDP=0 : Bonn potential)



C ***   Inizializzazioni
      DO J=0,99
          PP(J) = 0.01*J
      END DO
      ICP = 0
      ICN = 0
      IDEU = 0
      W_N = 0.
      W_DEU = 0.
      ECM = 0.
      DO J=1,5
        DO JJ=1,3
          PN_N_D(J,JJ) = 0.
          PN_D(J,JJ) = 0.
        END DO
      END DO

C *************************************
C ***   Target selection
C ***   NCHAIN_P = number of selected channels on proton
C ***   NCHAIN_N = number of selected channels on neutron
C *************************************

C *** Estrazione di 5 numeri random
 3    continue
 
      call ranmar(RND,5)

 1    continue

C ***   Computation of W_N_D and BETA_D for gamma-nucleon interactions
C ***   one nucleon is interacting, the other one is spectator
      A = ACOS(1.-2.*RND(2))
      B = 2.*PI*RND(3)
      Q = RND(4)
      DO J=0,99
         IF (IDP.EQ.1) THEN
            IF (Q.GE.DISTR_P(J).AND.Q.LT.DISTR_P(J+1)) GOTO 2
         ELSE
            IF (Q.GE.DISTR_B(J).AND.Q.LT.DISTR_B(J+1)) GOTO 2
         END IF
      END DO
 2    P1 = PP(J)+RND(5)*(PP(J+1)-PP(J))

C *** Check on nucleon mass inside deuteron
      totmass=deumass/2.
      if (abs(p1).gt.totmass) go to 3

C **************************************************
C *** Calculation of nucleon masses inside the deuterons
C *** (effective masses)
C **************************************************
      rm(1)=mass_n_d(p1)
      rm(2)=mass_n_d(p1)

C *** Target nucleon
      PN_N_D(1,2) = ABS(P1)*SIN(A)*COS(B)
      PN_N_D(2,2) = ABS(P1)*SIN(A)*SIN(B)
      PN_N_D(3,2) = ABS(P1)*COS(A)
      PN_N_D(5,2) = SQRT(PN_N_D(1,2)**2+PN_N_D(2,2)**2+PN_N_D(3,2)**2)
      PN_N_D(4,2) = SQRT(PN_N_D(5,2)**2+RM(2)**2)
C *** Spectator nucleon
      PN_N_D(1,1) = -PN_N_D(1,2)
      PN_N_D(2,1) = -PN_N_D(2,2)
      PN_N_D(3,1) = -PN_N_D(3,2)
      PN_N_D(5,1) = SQRT(PN_N_D(1,1)**2+PN_N_D(2,1)**2+PN_N_D(3,1)**2)
      PN_N_D(4,1) = SQRT(PN_N_D(5,1)**2+RM(1)**2)

C ***********************************************************
C *** Spectator nucleon has reduced mass. I should subtract
C *** from total energy for the reaction the mass variation of
C *** the spectator 
C ***********************************************************

      ETOT = EGAM+PN_N_D(4,1)

C *** Protone e neutrone hanno masse uguali
      rm_fin= RMASS('proton      ')
      E_sp_fin=sqrt(PN_N_D(5,1)**2+rm_fin**2.)
      delta_e=E_sp_fin-PN_N_D(4,1)
      etot=etot-delta_e

C *** CM  beta of fotone+interacting nucleon
      BETA_N_D(1)=PN_N_D(1,1)/ETOT
      BETA_N_D(2)=PN_N_D(2,1)/ETOT
      BETA_N_D(3)=(PN_N_D(3,1)+EGAM)/ETOT
      BETA_N_D(4)=0.

C *** Variabile per impulso di FErmi nella ntupla
C *** Fermi motion is computed in the frame in wich deuteron is at 
C *** rest, i.e. in LAB
      p_fermi=abs(P1)


C *************************************************************
C ***   S is Lorentz invariant. Scalar product has only third 
C ***   component if photon is along z
C *************************************************************
C ***
C *** Square of (gamma+target nucleon)
C ***
      S = RM(1)**2+2.*PN_N_D(4,1)*EGAM-2.*PN_N_D(3,1)*EGAM

C ***
C *** subtracting spectator mass variation
C ***
      S = S +delta_e**2.-2.*delta_e*(egam+PN_N_D(4,1))

C *********************************************
C *** Check on total energy for the reaction
      totmass1=rmass('proton      ')
      totmass2=rmass('pi-         ')
      totmass=(totmass1+totmass2)**2.
      if (s.lt.totmass+0.000001) go to 3
C *********************************************

      W_N = SQRT(S)

      RETURN
      END                                               !*** END FERMI_D
      
      real function mass_n_d(pf)

C ******************************
C *** This function computes reduced mass
C *** of a nucleon inside the deuteron
C *** from fermi motion so that energy
C *** is conserved
C ******************************

      rmn=rmass('deuteron    ')/2.

      x=sqrt(rmn*rmn-pf*pf)

      mass_n_d=x


      return
      end

