C***********************************************************************

      SUBROUTINE EV_NUC

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                                                           
C FUNCTIONAL DESCRIPTION:                                       
C                                                           
C    This subroutine generates the physical nuclear event, defines the
C    kinematics of selected process doing the proper LAB-CM transformation,
C    considering the decay channels and momentum distribution for processes
C    with bound nucleons.
C
C MODIFICATION DATE:
C 
C     02.05.95       by Igor Pshenichnov
C
C   ************     by Marco mirazita
C
C     23.04.99
C --> Aggiunto il file di include NT_KINE.INC con i dati delle particelle
C     dopo la prima fase di reazione (prima dei decadimenti).
C --> Aggiunta la chiamata alla routine PRIM_KINE che registra la 
C     cinematica delle particelle prima dei decadimenti.
C
C     29.04.99
C --> Aggiunto il file di include PROCESS.INC con le variabili per il
C     conteggio dei canali di reazione selezionati
C     
C     05.05.99  
C --> Eliminato He3. Se TARGET = HE3 ==> errore e il programma si ferma.
c --> Azzerata la variabile TARGFERMI all'inizio e ricalcolata dopo 
C     l'estrazione del canale di reazione, in base al canale estratto.
C --> Aggiunto il canale di reazione 35 gamma d --> rho0 d. Sono state 
C     modificate tutti i punti in cui vengono calcolate le varie grandezze
C     cinematiche per il nuovo canale. (Aggiunte anche alcune chiamate per
C     il canale 36 gamma d --> eta d, ma e' eventualmente da finire).
C
C   ************
C
C FORMAL PARAMETERS:                                              
C                                                            
C    NP          = NUMBER OF OUTGOING PARTICLES (2-18)
C    ECM         = TOTAL CM ENERGY (REAL)                            
C    AMASS       = MASS OF I-TH PARTICLE (REAL)
C    KGENEV      = CONSTANT CROSS SECTION (=1), FERMI DISTRIBUTION (=2)
C    PCM(1,I)    = PX OF I-TH PARTICLE             (CM-SYSTEM)
C    PCM(2,I)    = PY OF I-TH PARTICLE             (CM-SYSTEM)
C    PCM(3,I)    = PZ OF I-TH PARTICLE             (CM-SYSTEM)
C    PCM(4,I)    = ENERGY OF I-TH PARTICLE         (CM-SYSTEM)
C    PCM(5,I)    = |P| OF I-TH PARTICLE            (CM-SYSTEM)
C    WT          = WEIGHT OF THE EVENT
C    PLAB(1,I)   = PX OF I-TH PARTICLE             (LAB-SYSTEM)
C    PLAB(2,I)   = PY OF I-TH PARTICLE             (LAB-SYSTEM)
C    PLAB(3,I)   = PZ OF I-TH PARTICLE             (LAB-SYSTEM)
C    PLAB(4,I)   = ENERGY OF I-TH PARTICLE         (LAB-SYSTEM)
C    PLAB(5,I)   = |P| OF I-TH PARTICLE            (LAB-SYSTEM)
C    KINETIC ENERGY = PLAB(4,I) - AMASS(I)                            
C                                                           
C COMMON BLOCKS:                                              
C                                                            
C     /MONTECARL/                                              
C     /FLAG_COMP/                                                  
C     /TARGET/                                                    
C     /MASSES/                                                    
C     /NUCLEUS/                                                      
C     /NUCLEON/                                                      
C     /H2/                                                      
C     /HE3/                                                    
C     /FLAG_PROC/                                              
C     /CHANNELS/                                              
C     /GENIN/                                                 
C     /GENOUT/
C     /MOM_NUCLEON/
C     /LOGEG/
C-----------------------------------------------------------------------
      IMPLICIT NONE

#include "process.inc"

      INTEGER IBEAM,IELET,NCHAIN,ICHAIN(40),ICIBLE,IPOL,ILAM
      REAL ETAG,FWHM
      COMMON/FLAG_COMP/IELET,IBEAM,ICIBLE,ETAG,FWHM,NCHAIN,ICHAIN,IPOL,
     &   ILAM
      REAL EGAM, P(4,11)
      INTEGER I(11),NPART
      COMMON/MONTECARL/ EGAM,I,NPART,P
      INTEGER ITARG,IIDH
      CHARACTER*12 TARGET,TARGFERMI
      COMMON /TARGET/ ITARG,IIDH,TARGET,TARGFERMI
      CHARACTER*20 FILEN,FILEP
      COMMON/FNAME/FILEN,FILEP
      REAL PN_N(5,3),BETA_N(4)
      COMMON /NUCLEON/ PN_N,BETA_N
      REAL PN_N_D(5,3),PN_D(5,3),BETA_N_D(4),BETA_D(4)
      COMMON /H2/ PN_N_D,PN_D,BETA_N_D,BETA_D
      REAL PN_N_HE3(5,3),PN_DELTA(5,3),PN_DEU(5,3),PN_DIB(5,3),
     &     BETA_N_HE3(4),BETA_DELTA(4),BETA_DEU(4),BETA_DIB(4)
      COMMON /HE3/ PN_N_HE3,PN_DELTA,PN_DEU,PN_DIB,
     &     BETA_N_HE3,BETA_DELTA,BETA_DEU,BETA_DIB
      REAL RMB1,RMB2,RMME,DELTAMASS,DIBMASS,DEUMASS,RHOMASS
      COMMON /MASSES/ RMB1,RMB2,RMME,DELTAMASS,DIBMASS,DEUMASS,
     &     RHOMASS
      CHARACTER*12 PARTIC(40,0:3,0:3,8)
      INTEGER JDEC(40),NCH(40),NPLEV(40,3)
      REAL BR(40,3)
      COMMON /CHANNELS/ PARTIC,JDEC,NCH,NPLEV,BR
      INTEGER NP,KGENEV
      REAL ECM,AMASS,PCM,WT
      COMMON/GENIN/NP,ECM,AMASS(18),KGENEV
      COMMON/GENOUT/PCM(5,18),WT
      INTEGER NERR,LOOPV,NCASEG,LOOPC,NPROC,JCH
      COMMON /FLAG_GEN/ NERR,LOOPV,NCASEG,LOOPC,NPROC
      REAL RNDM,PROB_TOP(5),S5_9,TOP_GAMN
      INTEGER IPRINT,IVAR,JVAR,IDIHOTOMIA
      INTEGER*4 TOP_NUMBER, PART_NUMB(10)

      REAL MASSTARG,RMASS,W_N,W_DEU,W_DELTA,W_DIB,W,RMN,BETA(4),
     &     DELTA(4),
     &     PN(5,3),TOTMASS,WW,PCM1(5,18),PCM2(5,18),PN_T(4),T,U,Q,YRNDM,
     &     T_COEFF,PL(4),PROMUL(0:8),PLAB1(5,18),PLAB2(5,18),P1(4),P2(4)
      INTEGER A,Z,J,JJ,K,KK,LEV,JBR,JM,IA,IB,JL,IL,JD,
     &            NUMPLUS,NUMPAR,CODE,ITH
      COMMON /MOM_NUCLEON/PN
      LOGICAL EG,CMSF
      COMMON/LOGEG/EG,CMSF
      REAL AMASS1(18),AMASS2(18), ARGUMENT
      real acosdd

C *** common for ntuple variables
#include "nt_kine.inc"
      equivalence (plab1,pprim)
      equivalence (jch,j_channel)

      real vran(1),amx,mass_n_d,s


C ***   Inizialization
      KGENEV = 1
      DO J=1,11
         DO JJ=1,4
            P(JJ,J) = 0.
         END DO
         I(J) = 0
      END DO
      MASSTARG = 0.
      W_N = 0.
      W_DEU = 0.
      W_DELTA = 0.
      W_DIB = 0.
      W = 0.
      RMN = 0.
      JCH = 0
      DO J=1,4
          BETA(J) = 0.
          P1(J) = 0.
          P2(J) = 0.
          DELTA(J) = 0.
      END DO
      DO J=1,5
         DO JJ=1,3
            PN(J,JJ) = 0.
         END DO
         DO JJ=1,18
            PCM1(J,JJ) = 0.
            PCM2(J,JJ) = 0.
            PLAB1(J,JJ) = 0.
            PLAB2(J,JJ) = 0.
         END DO
      END DO
      LEV = 0
      JBR = 0
      TOTMASS = 0.
      RMB1 = 0.
      RMB2 = 0.
      RMME = 0.
      WW = 0.
      NUMPLUS = 0
      NUMPAR = 0
      DO J=0,8
         PROMUL(J) = 0.
      END DO
      ITH = 0
      T_COEFF = 3.
C **********************************
C *** Azzero anche TARGFERMI
C **********************************
      targfermi=' '
C *** Impulso di fermi
      p_fermi=0.


C ***  If the energy of gamma is below the pion threshold than keep
C ***  gamma in output array and nothing more (Target nucleon at rest)
      IF(EGAM.LE.0.151) THEN
         NPART  = 1
         P(1,1) = 0.
         P(2,1) = 0.
         P(3,1) = EGAM
         I(1)   = 1
         RETURN
      ENDIF

C ***   Definition of masses of decaying 
C ***   particles (Delta, Rho, Omega)
 11   DELTAMASS = RMASS('delta++     ')
      RHOMASS =   RMASS('rho0        ')
      DIBMASS =   RMASS('dibaryon    ')
      DEUMASS =   RMASS('deuteron    ')


C *************************************************************
C ***   Computing CM energy and beta, considering, if necessary,
C ***   the Fermi motion
C *************************************************************
      IF (TARGET.EQ.'deuteron    ') THEN
C *** Calcolo moto di fermi per i nucleoni
         if (nchain_p+nchain_n.ne.0) then
            CALL FERMI_D(W_N,W_DEU)
         endif
C *** Kinematical quantities for deuteron
         S = DEUMASS**2+2.*DEUMASS*EGAM
         W_DEU = SQRT(S)
         PN_D(1,1) = 0.
         PN_D(2,1) = 0.
         PN_D(3,1) = 0.
         PN_D(4,1) = DEUMASS
C *** Beta(CM) for (gamma+deuteron) system
         BETA_D(1) = 0.
         BETA_D(2) = 0.
         BETA_D(3) = EGAM/(EGAM+DEUMASS)
         BETA_D(4) = 0.
      ELSE IF (TARGET.EQ.'He3         ') THEN
         WRITE(15,*) '***********************************'
         WRITE(15,*) ' ERROR: Target He3 not implemented'
         WRITE(15,*) '***********************************'
         stop
c         CALL FERMI_HE3(W_N,W_DELTA,W_DEU,W_DIB)
      ELSE IF(TARGET.EQ.'proton      '.OR.TARGET.EQ.'neutron     ')THEN
C ****
c         RMN = MASSTARG
         RMN = rmass(target)
C ****
         PN_N(4,1) = RMN
         W_N = SQRT(RMN*RMN + 2.*RMN*EGAM)
         BETA_N(1) = 0.
         BETA_N(2) = 0.
         BETA_N(3) = EGAM/(EGAM+RMN)
         BETA_N(4) = 0.
      END IF

C *** Computing beta(CM) of (gamma+target) system  
      if (TARGET.EQ.'proton      '.OR.TARGET.EQ.'neutron     ') THEN
         do j=1,3
            beta_cm(j)=beta_n(j)
         enddo
      else if (TARGET.EQ.'deuteron    ') THEN
         do j=1,3
            beta_cm(j)=beta_d(j)
         enddo
      endif
      beta_cm(4)=1./SQRT(1.-(BETA_cm(1)*2.+BETA_cm(2)**2.
     &     +BETA_cm(3)**2.))


C *************************************************************
C ***   Selecting reaction channel
C *************************************************************
      IF (NCHAIN.EQ.1) THEN
         JCH = ICHAIN(1)
      ELSE
         CALL CHANNELEX(W_N,W_DELTA,W_DEU,W_DIB,JCH)
      END IF

C **********************************************************************
C *** L'impulso di Fermi e' calcolato in ogni caso.
C *** Se la reazione e' su d l'impulso di Fermi deve essere nullo
C **********************************************************************
      if (jch.ge.33) p_fermi=0.


C *****************************************************
C *** Calcolo TARGFERMI in base al canale estratto
C *** TARGFERMI is the mass of the interacting particle
C *****************************************************
C *** Protone
      if ( (jch.eq.1).or.(jch.eq.2).or.
     +     ((jch.ge.5).and.(jch.le.7)).or.
     +     (jch.eq.11).or.(jch.eq.12).or.
     +     (jch.eq.15).or.(jch.eq.16).or.
     +     (jch.eq.19).or.(jch.eq.21).or.
     +     (jch.eq.23).or.(jch.eq.24).or.
     +     (jch.eq.27).or.
     +     ((jch.ge.29).and.(jch.le.32))
     *                                   ) then

C *** Massa del bersaglio (incluso eventualmente il moto di Fermi)
         if (target.eq.'deuteron') then
            masstarg=rm(1)
         else
            MASSTARG = RMASS('proton      ')
         endif
         TARGFERMI = 'proton      '

C *** Neutrone
      else if ( (jch.eq.3).or.(jch.eq.4).or.
     +        ((jch.ge.8).and.(jch.le.10)).or.
     +        (jch.eq.13).or.(jch.eq.14).or.
     +        (jch.eq.17).or.(jch.eq.18).or.
     +        (jch.eq.20).or.(jch.eq.22).or.
     +        (jch.eq.25).or.(jch.eq.26).or.
     +        (jch.eq.28)
     *                                ) then

C *** Massa del bersaglio (incluso eventualmente il moto di Fermi)
         if (target.eq.'deuteron') then
            masstarg=rm(1)
         else
            MASSTARG = RMASS('neutron     ')
         endif
         TARGFERMI = 'neutron     '

C *** Deuterio
      else if ( (jch.ge.33).and.(jch.le.36)
     +                                ) then
         MASSTARG = RMASS('deuteron    ')
         TARGFERMI = 'deuteron    '

C *** He3
      else if ( (jch.ge.37).and.(jch.le.40)
     +                                ) then
         TARGFERMI = 'He3 '

      endif

C ***********************************************
C *** fine modifica
C ***********************************************


C ***********************************************
C ***    BETA calculation
C *** For bound nucleons (TARGET=deuteron) 
C *** photon+interacting nucleon must be considered
C ***********************************************

 12   IF (JCH.LE.32) THEN
         W = W_N
         DO J=1,4
            IF (TARGET.EQ.'deuteron') BETA(J) = -BETA_N_D(J)
c            IF (TARGET.EQ.'He3') BETA(J) = -BETA_N_HE3(J)
            IF (TARGET.EQ.'proton'.OR.TARGET.EQ.'neutron')
     &           BETA(J) = -BETA_N(J)
         END DO
         DO J=1,5
            DO JJ=1,3
               IF (TARGET.EQ.'deuteron') PN(J,JJ)=PN_N_D(J,JJ)
               IF (TARGET.EQ.'He3') PN(J,JJ)=PN_N_HE3(J,JJ)
               IF (TARGET.EQ.'proton'.OR.TARGET.EQ.'neutron')
     &              PN(J,JJ)=PN_N(J,JJ)
            END DO
         END DO
      END IF

C ***************************************************************
C ***    BETA calculation
C *** For diffractive reactions on deuteron, 
C *** photon+deuteron must be considered
C ***************************************************************
      IF ( (JCH.gE.33).and.(jch.le.36) ) THEN
         W = W_DEU
         DO J=1,4
            BETA(J) = -BETA_D(J)
         END DO
         DO J=1,5
            DO JJ=1,3
               PN(J,JJ)=PN_D(J,JJ)
            END DO
         END DO
      END IF

C ***************************************************************
C ***     CALCOLO DI BETA
C *** Le reazioni su He3 sono da implementare
C ***************************************************************
      IF (JCH.EQ.37) THEN
         W = W_DELTA
         DO J=1,4
            BETA(J) = -BETA_DELTA(J)
         END DO
         DO J=1,5
            DO JJ=1,3
               PN(J,JJ)=PN_DELTA(J,JJ)
            END DO
         END DO
      END IF

      IF (JCH.EQ.38.OR.JCH.EQ.39) THEN
         W = W_DEU
         DO J=1,4
            BETA(J) = -BETA_DEU(J)
         END DO
         DO J=1,5
            DO JJ=1,2
               PN(J,JJ)=PN_DEU(J,JJ)
            END DO
         END DO
      END IF

      IF (JCH.EQ.40) THEN
         W = W_DIB
         DO J=1,4
            BETA(J) = -BETA_DIB(J)
         END DO
         DO J=1,5
            DO JJ=1,3
               PN(J,JJ)=PN_DIB(J,JJ)
            END DO
         END DO
      END IF


C ***********************************************
C ***********************************************
C *** KINEMATICS OF FINAL STATE
C ***********************************************
C ***********************************************


C *************************************
C ***   First level
C *************************************
      LEV = 1
      JBR = 0

      IF((JCH.LE.0).OR.(JCH.GT.40)) THEN
         WRITE(0,112) JCH
 112     FORMAT(5X,'JCH=',I5)
      ENDIF 

C *************************************************************
C *** Check on final masses
C *************************************************************

      TOTMASS = 0.

C *** Multipion channels
C *** Treat these channels as a special case, select mutiplicity here

      IF(JCH.EQ.27.OR.JCH.EQ.28) THEN
         S5_9=0.0
         DO IVAR=1,5
            S5_9=S5_9+TOP_GAMN(EGAM,IVAR+4)
            PROB_TOP(IVAR)=S5_9
         ENDDO
         DO IVAR=1,5
            PROB_TOP(IVAR)=PROB_TOP(IVAR)/PROB_TOP(5)
         ENDDO
         call ranmar(vran,1)
         TOP_NUMBER=IDIHOTOMIA(vran(1),PROB_TOP,5)+4

         NP=TOP_NUMBER
         DO JM=1,NP-1           ! put pion masses
            AMASS(JM) = RMASS(PARTIC(27,1,0,1))
            AMASS1(JM) = AMASS(JM)
            TOTMASS = TOTMASS+AMASS(JM)
         END DO
         AMASS(NP)=RMASS(PARTIC(27,1,0,5)) ! put nucleon mass 
         AMASS1(NP) = AMASS(NP)
         TOTMASS = TOTMASS+AMASS(NP)

C ***  For other channels ...

      ELSE
         NP = NPLEV(JCH,1)
         DO JM=1,NP
            AMASS(JM) = RMASS(PARTIC(JCH,LEV,JBR,JM))
            AMASS1(JM) = AMASS(JM)
            TOTMASS = TOTMASS+AMASS(JM)
         END DO
      ENDIF


      IF (W.LT.TOTMASS+0.000001) THEN
         ITH = ITH+1
         IF ((MOD(ITH,100000).EQ.0).AND.(ITH.GT.100000)) THEN
       !         WRITE (6,*) 'REACTION IMPOSSIBLE '
       !         WRITE (6,*) 'PHOTON ENERGY UNDER THRESHOLD'
            WRITE (6,100) JCH,(AMASS(JJ),JJ=1,NP)
            WRITE (6,101) W,TOTMASS,EGAM
 100        FORMAT(5X,'CHANNEL NUMBER=',I4/5X,'AMASS= ',18F10.6)
 101        FORMAT(5X,'W= ',F10.6,'  TOTMASS= ',F10.6,'  EGAM= ',F10.6)
         END IF
         GOTO 11
      END IF


C ******************************************************************************
C ***   Choosing masses for target baryon, scattered baryon, spectator
C ***   Baryon, if any, and of scattered mesons 
C ******************************************************************************
C ***
C *** channels  1 - 14 e 19 - 22, 29, 30, 32, 33, 34, 35, 36
C ***
      IF ( ((JCH.GE.1).AND.(JCH.LE.14)) .OR.(JCH.EQ.19).OR.(JCH.EQ.20)
     &     .OR.(JCH.EQ.29).OR.(JCH.EQ.30).OR.( (JCH.GE.32)
     &     .AND.(JCH.LE.36)) ) THEN
 
C ****************************
C *** Target
         RMB1 = masstarg
C ****************************
C *** Scattered baryon (nucleon)
         RMB2 = AMASS(2)
C *** Scattered meson(pion)
         RMME = AMASS(1)

C ***
C *** If scattered Baryon is delta
C ***
         IF (JCH.GE.5.AND.JCH.LE.10) THEN
            WW = DELTAMASS
            RMB2 = WW
            IF (W.LT.(RMME+RMB2+0.000001)) GOTO 11
         END IF
C ***
C *** If scattered meson is rho
C ***
         IF ((JCH.GE.11.AND.JCH.LE.14).or.(jch.eq.35)) THEN
            WW = RHOMASS
            RMME = WW
            IF (W.LT.(RMME+RMB2+0.000001)) GOTO 11
         END IF
      END IF

C ***
C ***  If scattered meson is omega
C ***
      IF (JCH.EQ.21.OR.JCH.EQ.22) THEN
         RMB2 = AMASS(2)
         WW = AMASS(1)
         RMME = WW
         IF (W.LT.(RMME+RMB2+0.000001)) GOTO 11
      END IF

C ***
C *** He3 reaction
C***
      IF (JCH.EQ.40) WW=AMASS(2)

C **************************************************
C *** Prova stampa
C **************************************************
c      if (TARGFERMI.eq.'deuteron    ') then
c         amx=0.
c         write (*,35)amx,amx,amx,RMASS(targfermi),RMASS(targfermi)
c      else
c         amx=sqrt(pn_n_d(4,1)**2.-
c     +        (pn_n_d(1,1)**2.+pn_n_d(2,1)**2.+pn_n_d(3,1)**2.))
c         write (*,37) (pn_n_d(j,1),j=1,4),amx
c         amx=sqrt(pn_n_d(4,2)**2.-
c     +        (pn_n_d(1,2)**2.+pn_n_d(2,2)**2.+pn_n_d(3,2)**2.))
c         write (*,38) (pn_n_d(j,2),j=1,4),amx
c      endif
c 35   format (1x,'deut. ',5(1x,f8.5))
c 37   format (1x,'N.int.',5(1x,f8.5))
c 38   format (1x,'N. sp.',5(1x,f8.5))



C **********************************************************************
C *** Generating final state
C **********************************************************************
 20   IF ((JCH.GE.15.AND.JCH.LE.18).OR.(JCH.GE.23.AND.JCH.LE.28)
     &     .OR.JCH.EQ.31.OR.JCH.GE.37) THEN
C ***
C *** Channels with 2 pions or more
C ***
         ECM = W
C!    PRINT 113,JCH,W,(AMASS(IPRINT),IPRINT=1,NP)
 113     FORMAT (1X,'JCH=',I4,' W=',F9.5,' AMASS=',18F9.5)
         CALL GENBOD
      ELSE
C ***
C *** Other channels
C ***
         CALL GEN_EVT(W,LEV,JCH)
      END IF



C **********************************************************************
C *** CM final momenta 
C **********************************************************************
      DO IA=1,NP
         DO IB=1,5
            PCM1(IB,IA)=PCM(IB,IA)
         END DO
      END DO


C **********************************************************************
C *** Lorentz Transformations in LAB system
C **********************************************************************
C ***      BETA -> -BETA  was assigned previously for Lorentz transformation
C ***      from CMS to Lab system
C **********************************************************************
	BETA(4)  = 1./SQRT(1.-(BETA(1)*BETA(1)+BETA(2)*BETA(2)
     &    +BETA(3)*BETA(3)))

C Transform to LAB SYS
	DO J=1,NP
            DO IL=1,4
               P1(IL)=PCM1(IL,J)
            END DO
            CALL GLOREN(BETA,P1,PL)
            DO IL=1,4
               PLAB1(IL,J)=PL(IL)
            END DO
            PLAB1(5,J)=SQRT(PL(1)*PL(1)+PL(2)*PL(2)+PL(3)*PL(3))
	END DO


C **********************************************************************
C ***   Computing t for two pions production (phase space) and
C ***   cheking for right values 
C **********************************************************************
      IF (JCH.GE.15.AND.JCH.LE.18) THEN
         DO IL=1,4
            PN_T(IL) = PLAB1(IL,3)-PN(IL,1)
         END DO
         T = PN_T(1)*PN_T(1)+PN_T(2)*PN_T(2)+PN_T(3)*PN_T(3)
         T = PN_T(4)*PN_T(4)-T
         U = EXP(T_COEFF*T)
         call ranmar(vran,1)
         YRNDM=vran(1)

         IF (YRNDM.GT.U) GOTO 20
      END IF



C **********************************************************************
C **********************************************************************
C *** Second leve: decays
C **********************************************************************
C **********************************************************************





      IF(JDEC(JCH).EQ.0) THEN
      
C **********************************************************************
C ***   special channel on He3, Delta knockout
                              ! Check conservation law for channel #37
                              ! in the lab system.
                              ! For other cases it has to be done 
                              ! outside the generator body...
C **********************************************************************
         IF (JCH.EQ.37) THEN
            DELTA(3) = (PN(3,1))-(PLAB1(3,1)+PLAB1(3,2)+PLAB1(3,3)+
     &           PLAB1(3,4)+PLAB1(3,5)+PLAB1(3,6))
            DELTA(4) = (PN(4,1))-(PLAB1(4,1)+PLAB1(4,2)+PLAB1(4,3)+
     &           PLAB1(4,4)+PLAB1(4,5)+PLAB1(4,6))
         END IF
         IF (ABS(DELTA(1)).GT.1.E-04.OR.ABS(DELTA(2)).GT.1.E-04.OR.
     &        ABS(DELTA(3)).GT.1.E-04.OR.ABS(DELTA(4)).GT.1.E-04) THEN
            NERR = NERR+1
            GOTO 11
         END IF

C **********************************************************************
C   Put angular distribution (C.M.S.) for the first particle in the
C    array for each channel
C **********************************************************************
         IF(PCM(5,1).GT.1.E-07) THEN
            ARGUMENT=PCM(3,1)/PCM(5,1)
            IF(ABS(ARGUMENT).LE.1.) THEN 
!               CALL HFILL(-JCH,ACOSDD(ARGUMENT),0.,1.)
            ELSEIF (ARGUMENT.LT.-1.) THEN
!               CALL HFILL(-JCH,180.,0.,1.)
            ELSEIF (ARGUMENT.GT.1.) THEN
!               CALL HFILL(-JCH,0.,0.,1.)
            ENDIF
         ENDIF

      ELSE

C **********************************************************************
C *** DECAYS
C **********************************************************************

C ********************************************
C *** Calling the routine for extraction of kinematics 
C *** of first level particles                       !M.Mirazita
C ********************************************
C *** Prova per aggiunta canale 35
C *** Finiva a jch=14
C ********************************************

         if ( ((jch.ge.4).and.(jch.le.14)).or.(jch.eq.35) ) then
            call prim_kine
         endif


C *************************************************************************


         IF(PCM(5,1).GT.1.E-07) THEN
            ARGUMENT=PCM(3,1)/PCM(5,1)
            IF(ABS(ARGUMENT).LE.1.) THEN 
!               CALL HFILL(-JCH,ACOSDD(ARGUMENT),0.,1.)
            ELSEIF (ARGUMENT.LT.-1.) THEN
!               CALL HFILL(-JCH,180.,0.,1.)
            ELSEIF (ARGUMENT.GT.1.) THEN
!               CALL HFILL(-JCH,0.,0.,1.)
            ENDIF
         ENDIF


C **********************************************************************
C ***   Choosing Branching Ratio for second level channel definition
C **********************************************************************
         JBR=1
         IF (NCH(JCH).GE.2) THEN
            PROMUL(0) = 0.
            DO K=1,NCH(JCH)
               PROMUL(K) = PROMUL(K-1)+BR(JCH,K)
            END DO
C *** random
            call ranmar(vran,1)
            Q=vran(1)


            DO K=1,NCH(JCH)
               IF (Q.LE.PROMUL(K)) GOTO 66
            END DO
 66         JBR = K
         END IF
         NP=NPLEV(JCH,2)
         JD=JDEC(JCH)

C **********************************************************************
C ***   Defining  BETA =-(P/ETOT) of decayng particle
C **********************************************************************
         IF (JD.GT.6) THEN
            WRITE(15,111) JD,JCH,JBR,NP
 111        FORMAT(1X,' JD=',I6,' JCH=',I6,' JBR=',I6,' NP=',I6)
         ENDIF
         DO K=1,3
            BETA(K)=-PLAB1(K,JD)/PLAB1(4,JD)
         ENDDO
         BETA(4)=0.
         W=WW
         TOTMASS = 0.
         LEV=2
         DO JM=1,NP
            AMASS(JM)=RMASS(PARTIC(JCH,LEV,JBR,JM))
            AMASS2(JM) = AMASS(JM)        
            TOTMASS = TOTMASS+AMASS(JM)
         ENDDO

         IF (W.LT.(TOTMASS+0.000001)) GOTO 11

C **********************************************************************
C ***   First step: decay with GENBOD
C **********************************************************************
         ECM = W
         CALL GENBOD
         DO IA=1,NP
            DO IB=1,5
               PCM2(IB,IA)=PCM(IB,IA)
            ENDDO
         ENDDO

C **********************************************************************
C ***      BETA -> -BETA  was assigned previously for Lorentz transformation
C ***      from CMS to Lab system
C **********************************************************************
         BETA(4)  = 1./SQRT(1.-(BETA(1)*BETA(1)+BETA(2)*BETA(2)
     &        +BETA(3)*BETA(3)))
         DO J=1,NP
            DO IL=1,4
               P2(IL)=PCM2(IL,J)
            ENDDO
            CALL GLOREN(BETA,P2,PL)
            DO IL=1,4
               PLAB2(IL,J)=PL(IL)
            ENDDO
            PLAB2(5,J)=SQRT(PL(1)*PL(1)+PL(2)*PL(2)+PL(3)*PL(3))
         ENDDO
         
         IF (JCH.EQ.40) THEN    ! Check conservation law for channel #40
                           ! in the lab system.
                           ! For other cases it has to be done 
                           ! outside the generator body...
            DELTA(1) = PLAB1(1,1)+PLAB2(1,1)+PLAB2(1,2)+PN(1,2)+PN(1,3)
            DELTA(2) = PLAB1(2,1)+PLAB2(2,1)+PLAB2(2,2)+PN(2,2)+PN(2,3)
            DELTA(3) = EGAM+PN(3,1)-(PLAB1(3,1)+PLAB2(3,1)+PLAB2(3,2))
            DELTA(4) = EGAM+PN(4,1)-(PLAB1(4,1)+PLAB2(4,1)+PLAB2(4,2))
         END IF
         IF (ABS(DELTA(1)).GT.1.E-04.OR.ABS(DELTA(2)).GT.1.E-04.OR.
     &        ABS(DELTA(3)).GT.1.E-04.OR.ABS(DELTA(4)).GT.1.E-04) THEN
            NERR = NERR+1
            GOTO 11
         END IF
      END IF

C **********************************************************************
C **********************************************************************
C *** Final momenta
C **********************************************************************
C **********************************************************************

      IF (TARGET.EQ.'deuteron') NUMPLUS = 1
      IF (TARGET.EQ.'He3') NUMPLUS = 2

C ***    Canali 1,2,3,4,19,20,21,22,32
      IF (JCH.LE.4  .OR. JCH.EQ.19 .OR. JCH.EQ.20 .OR.
     #     JCH.EQ.21 .OR. JCH.EQ.22 .OR. JCH.EQ.32) THEN
         NUMPAR=NPLEV(JCH,1)
         NPART = NUMPAR+NUMPLUS
         P(1,1) = PLAB1(1,1)
         P(2,1) = PLAB1(2,1)
         P(3,1) = PLAB1(3,1)
         P(4,1) = AMASS1(1)
         I(1) = CODE(PARTIC(JCH,1,0,1))
         P(1,2) = PLAB1(1,2)
         P(2,2) = PLAB1(2,2)
         P(3,2) = PLAB1(3,2)
         P(4,2) = AMASS1(2)
         I(2) = CODE(PARTIC(JCH,1,0,2))
         IF (TARGET.EQ.'deuteron'.OR.TARGET.EQ.'He3') THEN
            P(1,3) = PN(1,2)
            P(2,3) = PN(2,2)
            P(3,3) = PN(3,2)
            P(4,3) = RMASS(PARTIC(JCH,1,0,3))
            I(3) = CODE(PARTIC(JCH,1,0,3))
            IF (TARGET.EQ.'He3') THEN
               P(1,4) = PN(1,3)
               P(2,4) = PN(2,3)
               P(3,4) = PN(3,3)
               P(4,4) = RMASS(PARTIC(JCH,1,0,4))
               I(4) = CODE(PARTIC(JCH,1,0,4))
            END IF
         END IF
         GOTO 99
      END IF

C ***    Canali 5,6,7,8,9,10
      IF (JCH.GE.5.AND.JCH.LE.10) THEN
         NUMPAR=NPLEV(JCH,1)+NPLEV(JCH,2)-1
         NPART = NUMPAR+NUMPLUS
         P(1,1) = PLAB1(1,1)
         P(2,1) = PLAB1(2,1)
         P(3,1) = PLAB1(3,1)
         P(4,1) = AMASS1(1)
         I(1) = CODE(PARTIC(JCH,1,0,1))
         P(1,2) = PLAB2(1,1)
         P(2,2) = PLAB2(2,1)
         P(3,2) = PLAB2(3,1)
         P(4,2) = AMASS2(1)
         I(2) = CODE(PARTIC(JCH,2,JBR,1))
         P(1,3) = PLAB2(1,2)
         P(2,3) = PLAB2(2,2)
         P(3,3) = PLAB2(3,2)
         P(4,3) = AMASS2(2)
         I(3) = CODE(PARTIC(JCH,2,JBR,2))
         IF (TARGET.EQ.'deuteron'.OR.TARGET.EQ.'He3') THEN
            P(1,4) = PN(1,2)
            P(2,4) = PN(2,2)
            P(3,4) = PN(3,2)
            P(4,4) = RMASS(PARTIC(JCH,1,0,3))
            I(4) = CODE(PARTIC(JCH,1,0,3))
            IF (TARGET.EQ.'He3') THEN
               P(1,5) = PN(1,3)
               P(2,5) = PN(2,3)
               P(3,5) = PN(3,3)
               P(4,5) = RMASS(PARTIC(JCH,1,0,4))
               I(5) = CODE(PARTIC(JCH,1,0,4))
            END IF
         END IF
         GOTO 99
      END IF

C ***    Canali 11,12,13,14
      IF (JCH.GE.11.AND.JCH.LE.14) THEN
         NUMPAR=NPLEV(JCH,1)+NPLEV(JCH,2)-1
         NPART = NUMPAR+NUMPLUS
         P(1,1) = PLAB2(1,1)
         P(2,1) = PLAB2(2,1)
         P(3,1) = PLAB2(3,1)
         P(4,1) = AMASS2(1)
         I(1) = CODE(PARTIC(JCH,2,JBR,1))
         P(1,2) = PLAB2(1,2)
         P(2,2) = PLAB2(2,2)
         P(3,2) = PLAB2(3,2)
         P(4,2) = AMASS2(2)
         I(2) = CODE(PARTIC(JCH,2,JBR,2))
         P(1,3) = PLAB1(1,2)
         P(2,3) = PLAB1(2,2)
         P(3,3) = PLAB1(3,2)
         P(4,3) = AMASS1(2)
         I(3) = CODE(PARTIC(JCH,1,0,2))
         IF (TARGET.EQ.'deuteron'.OR.TARGET.EQ.'He3') THEN
            P(1,4) = PN(1,2)
            P(2,4) = PN(2,2)
            P(3,4) = PN(3,2)
            P(4,4) = RMASS(PARTIC(JCH,1,0,3))
            I(4) = CODE(PARTIC(JCH,1,0,3))
            IF (TARGET.EQ.'He3') THEN
               P(1,5) = PN(1,3)
               P(2,5) = PN(2,3)
               P(3,5) = PN(3,3)
               P(4,5) = RMASS(PARTIC(JCH,1,0,4))
               I(5) = CODE(PARTIC(JCH,1,0,4))
            END IF
         END IF
         GOTO 99
      END IF

C ***    Canali 15,16,17,18
      IF (JCH.GE.15.AND.JCH.LE.18) THEN
         NUMPAR=NPLEV(JCH,1)
         NPART = NUMPAR+NUMPLUS
         P(1,1) = PLAB1(1,1)
         P(2,1) = PLAB1(2,1)
         P(3,1) = PLAB1(3,1)
         P(4,1) = AMASS1(1)
         I(1) = CODE(PARTIC(JCH,1,0,1))
         P(1,2) = PLAB1(1,2)
         P(2,2) = PLAB1(2,2)
         P(3,2) = PLAB1(3,2)
         P(4,2) = AMASS1(2)
         I(2) = CODE(PARTIC(JCH,1,0,2))
         P(1,3) = PLAB1(1,3)
         P(2,3) = PLAB1(2,3)
         P(3,3) = PLAB1(3,3)
         P(4,3) = AMASS1(3)
         I(3) = CODE(PARTIC(JCH,1,0,3))
         IF (TARGET.EQ.'deuteron'.OR.TARGET.EQ.'He3') THEN
            P(1,4) = PN(1,2)
            P(2,4) = PN(2,2)
            P(3,4) = PN(3,2)
            P(4,4) = RMASS(PARTIC(JCH,1,0,4))
            I(4) = CODE(PARTIC(JCH,1,0,4))
            IF (TARGET.EQ.'He3') THEN
               P(1,5) = PN(1,3)
               P(2,5) = PN(2,3)
               P(3,5) = PN(3,3)
               P(4,5) = RMASS(PARTIC(JCH,1,0,5))
               I(5) = CODE(PARTIC(JCH,1,0,5))
            END IF
         END IF
         GOTO 99
      END IF
      

C ***    Canali 23,25
      IF (JCH.EQ.23.OR.JCH.EQ.25) THEN
         NUMPAR=NPLEV(JCH,1)
         NPART = NUMPAR+NUMPLUS
         P(1,1) = PLAB1(1,1)
         P(2,1) = PLAB1(2,1)
         P(3,1) = PLAB1(3,1)
         P(4,1) = AMASS1(1)
         I(1) = CODE(PARTIC(JCH,1,0,1))
         P(1,2) = PLAB1(1,2)
         P(2,2) = PLAB1(2,2)
         P(3,2) = PLAB1(3,2)
         P(4,2) = AMASS1(2)
         I(2) = CODE(PARTIC(JCH,1,0,2))
         P(1,3) = PLAB1(1,3)
         P(2,3) = PLAB1(2,3)
         P(3,3) = PLAB1(3,3)
         P(4,3) = AMASS1(3)
         I(3) = CODE(PARTIC(JCH,1,0,3))
         P(1,4) = PLAB1(1,4)
         P(2,4) = PLAB1(2,4)
         P(3,4) = PLAB1(3,4)
         P(4,4) = AMASS1(4)
         I(4) = CODE(PARTIC(JCH,1,0,4))
         IF (TARGET.EQ.'deuteron'.OR.TARGET.EQ.'He3') THEN
            P(1,5) = PN(1,2)
            P(2,5) = PN(2,2)
            P(3,5) = PN(3,2)
            P(4,5) = RMASS(PARTIC(JCH,1,0,5))
            I(5) = CODE(PARTIC(JCH,1,0,5))
            IF (TARGET.EQ.'He3') THEN
               P(1,6) = PN(1,3)
               P(2,6) = PN(2,3)
               P(3,6) = PN(3,3)
               P(4,6) = RMASS(PARTIC(JCH,1,0,6))
               I(6) = CODE(PARTIC(JCH,1,0,6))
            END IF
         END IF
         GOTO 99
      END IF

C ***    Channels #24 and #26        
      IF (JCH.EQ.24.OR.JCH.EQ.26) THEN
         NUMPAR=NPLEV(JCH,1)
         NPART = NUMPAR+NUMPLUS
         IF(JCH.EQ.24) THEN
            CALL GAMN_CHARG(4,1,PART_NUMB)
         ELSE
            CALL GAMN_CHARG(4,0,PART_NUMB)
         ENDIF          
         P(1,1) = PLAB1(1,1)
         P(2,1) = PLAB1(2,1)
         P(3,1) = PLAB1(3,1)
         P(4,1) = AMASS1(1)
         I(1) = PART_NUMB(2)           
         P(1,2) = PLAB1(1,2)
         P(2,2) = PLAB1(2,2)
         P(3,2) = PLAB1(3,2)
         P(4,2) = AMASS1(2)
         I(2) = PART_NUMB(3)                       
         P(1,3) = PLAB1(1,3)
         P(2,3) = PLAB1(2,3)
         P(3,3) = PLAB1(3,3)
         P(4,3) = AMASS1(3)
         I(3) = PART_NUMB(4)
         P(1,4) = PLAB1(1,4)
         P(2,4) = PLAB1(2,4)
         P(3,4) = PLAB1(3,4)
         P(4,4) = AMASS1(4)
         I(4) = PART_NUMB(1)
         IF (TARGET.EQ.'deuteron'.OR.TARGET.EQ.'He3') THEN 
            P(1,5) = PN(1,2)
            P(2,5) = PN(2,2)
            P(3,5) = PN(3,2)
            P(4,5) = RMASS(PARTIC(JCH,1,0,5))
            I(5) = CODE(PARTIC(JCH,1,0,5))
            IF (TARGET.EQ.'He3') THEN 
               P(1,6) = PN(1,3)
               P(2,6) = PN(2,3)
               P(3,6) = PN(3,3)
               P(4,6) = RMASS(PARTIC(JCH,1,0,6))
               I(6) = CODE(PARTIC(JCH,1,0,6))
            END IF
         END IF
         GOTO 99
      END IF

C ***    Canali 27,28    
      IF (JCH.EQ.27.OR.JCH.EQ.28) THEN
         NUMPAR=TOP_NUMBER
         NPART = NUMPAR+NUMPLUS
         IF(JCH.EQ.27) THEN
            CALL GAMN_CHARG(TOP_NUMBER,1,PART_NUMB)
         ELSE
            CALL GAMN_CHARG(TOP_NUMBER,0,PART_NUMB)
         ENDIF
         
        DO JVAR=1,TOP_NUMBER
           P(1,JVAR) = PLAB1(1,JVAR)
           P(2,JVAR) = PLAB1(2,JVAR)
           P(3,JVAR) = PLAB1(3,JVAR)
           P(4,JVAR) = AMASS1(JVAR)
           IF(JVAR.EQ.TOP_NUMBER) THEN
              I(JVAR) = PART_NUMB(1)
           ELSE
              I(JVAR) = PART_NUMB(JVAR+1)
           ENDIF
        ENDDO
        
        IF (TARGET.EQ.'deuteron'.OR.TARGET.EQ.'He3') THEN
           P(1,TOP_NUMBER+1) = PN(1,2)
           P(2,TOP_NUMBER+1) = PN(2,2)
           P(3,TOP_NUMBER+1) = PN(3,2)
           P(4,TOP_NUMBER+1) = RMASS(PARTIC(JCH,1,0,6))
           I(TOP_NUMBER+1) = CODE(PARTIC(JCH,1,0,6))
           IF (TARGET.EQ.'He3') THEN
              P(1,TOP_NUMBER+2) = PN(1,3)
              P(2,TOP_NUMBER+2) = PN(2,3)
              P(3,TOP_NUMBER+2) = PN(3,3)
              P(4,TOP_NUMBER+2) = RMASS(PARTIC(JCH,1,0,7))
              I(TOP_NUMBER+2) = CODE(PARTIC(JCH,1,0,7))
           END IF
        END IF
        GOTO 99
      END IF
      
C ***    Canali 29,30
      IF (JCH.EQ.29.OR.JCH.EQ.30) THEN
         NUMPAR=NPLEV(JCH,1)
         NPART = NUMPAR+NUMPLUS
         P(1,1) = PLAB1(1,1)
         P(2,1) = PLAB1(2,1)
         P(3,1) = PLAB1(3,1)
         P(4,1) = AMASS1(1)
         I(1) = CODE(PARTIC(JCH,1,0,1))
         P(1,2) = PLAB1(1,2)
         P(2,2) = PLAB1(2,2)
         P(3,2) = PLAB1(3,2)
         P(4,2) = AMASS1(2)
         I(2) = CODE(PARTIC(JCH,1,0,2))
         IF (TARGET.EQ.'deuteron'.OR.TARGET.EQ.'He3') THEN
            P(1,3) = PN(1,2)
            P(2,3) = PN(2,2)
            P(3,3) = PN(3,2)
            P(4,3) = RMASS(PARTIC(JCH,1,0,3))
            I(3) = CODE(PARTIC(JCH,1,0,3))
            IF (TARGET.EQ.'He3') THEN
               P(1,4) = PN(1,3)
               P(2,4) = PN(2,3)
               P(3,4) = PN(3,3)
               P(4,4) = RMASS(PARTIC(JCH,1,0,4))
               I(4) = CODE(PARTIC(JCH,1,0,4))
            END IF
         END IF
         GOTO 99
      END IF

C ***    Canale 31
      IF (JCH.EQ.31) THEN
         NUMPAR=NPLEV(JCH,1)
         NPART = NUMPAR+NUMPLUS
         P(1,1) = PLAB1(1,1)
         P(2,1) = PLAB1(2,1)
         P(3,1) = PLAB1(3,1)
         P(4,1) = AMASS1(1)
         I(1) = CODE(PARTIC(JCH,1,0,1))
         P(1,2) = PLAB1(1,2)
         P(2,2) = PLAB1(2,2)
         P(3,2) = PLAB1(3,2)
         P(4,2) = AMASS1(2)
         I(2) = CODE(PARTIC(JCH,1,0,2))
         P(1,3) = PLAB1(1,3)
         P(2,3) = PLAB1(2,3)
         P(3,3) = PLAB1(3,3)
         P(4,3) = AMASS1(3)
         I(3) = CODE(PARTIC(JCH,1,0,3))
         IF (TARGET.EQ.'deuteron'.OR.TARGET.EQ.'He3') THEN
            P(1,4) = PN(1,2)
            P(2,4) = PN(2,2)
            P(3,4) = PN(3,2)
            P(4,4) = RMASS(PARTIC(JCH,1,0,4))
            I(4) = CODE(PARTIC(JCH,1,0,4))
            IF (TARGET.EQ.'He3') THEN
               P(1,5) = PN(1,3)
               P(2,5) = PN(2,3)
               P(3,5) = PN(3,3)
               P(4,5) = RMASS(PARTIC(JCH,1,0,5))
               I(5) = CODE(PARTIC(JCH,1,0,5))
            END IF
         END IF
         GOTO 99
      END IF

C ***    Canale 33
      IF (JCH.EQ.33) THEN
         NUMPAR=NPLEV(JCH,1)
         NPART = NUMPAR
         P(1,1) = PLAB1(1,1)
         P(2,1) = PLAB1(2,1)
         P(3,1) = PLAB1(3,1)
         P(4,1) = AMASS1(1)
         I(1) = CODE(PARTIC(JCH,1,0,1))
         P(1,2) = PLAB1(1,2)
         P(2,2) = PLAB1(2,2)
         P(3,2) = PLAB1(3,2)
         P(4,2) = AMASS1(2)
         I(2) = CODE(PARTIC(JCH,1,0,2))
         IF (TARGET.EQ.'He3') THEN
            P(1,3) = PN(1,2)
            P(2,3) = PN(2,2)
            P(3,3) = PN(3,2)
            P(4,3) = RMASS(PARTIC(JCH,1,0,3))
            I(3) = CODE(PARTIC(JCH,1,0,3))
         END IF
         GOTO 99
      END IF
      
C ***    Canale 34
      IF (JCH.EQ.34) THEN
         NUMPAR=NPLEV(JCH,1)
         NPART = NUMPAR
         P(1,1) = PLAB1(1,1)
         P(2,1) = PLAB1(2,1)
         P(3,1) = PLAB1(3,1)
         P(4,1) = AMASS1(2)
         I(1) = CODE(PARTIC(JCH,1,0,1))
         P(1,2) = PLAB1(1,2)
         P(2,2) = PLAB1(2,2)
         P(3,2) = PLAB1(3,2)
         P(4,2) = AMASS1(2)
         I(2) = CODE(PARTIC(JCH,1,0,2))
         IF (TARGET.EQ.'He3') THEN
            P(1,3) = PN(1,2)
            P(2,3) = PN(2,2)
            P(3,3) = PN(3,2)
            P(4,3) = RMASS(PARTIC(JCH,1,0,3))
            I(3) = CODE(PARTIC(JCH,1,0,3))
         END IF
         GOTO 99
      END IF
      

      IF (JCH.eq.35) THEN
         NUMPAR=NPLEV(JCH,1)+NPLEV(JCH,2)-1         
         NPART = NUMPAR
C *** non va aggiunto NUMPLUS perche' alla fine c'e' D e non p-n

C *** pione di decadimento della ro
         P(1,1) = PLAB2(1,1)
         P(2,1) = PLAB2(2,1)
         P(3,1) = PLAB2(3,1)
         P(4,1) = AMASS2(1)
         I(1) = CODE(PARTIC(JCH,2,JBR,1))
C *** pione di decadimento della ro
         P(1,2) = PLAB2(1,2)
         P(2,2) = PLAB2(2,2)
         P(3,2) = PLAB2(3,2)
         P(4,2) = AMASS2(2)
         I(2) = CODE(PARTIC(JCH,2,JBR,2))
C *** deutone diffuso
         P(1,3) = PLAB1(1,2)
         P(2,3) = PLAB1(2,2)
         P(3,3) = PLAB1(3,2)
         P(4,3) = AMASS1(2)
         I(3) = CODE(PARTIC(JCH,1,0,2))
C *** va aggiunto l'altro nucleone per He3
         IF (TARGET.EQ.'He3') THEN
            P(1,4) = PN(1,2)
            P(2,4) = PN(2,2)
            P(3,4) = PN(3,2)
            P(4,4) = RMASS(PARTIC(JCH,1,0,3))
            I(4) = CODE(PARTIC(JCH,1,0,3))
         END IF
         GOTO 99
      END IF


C ***    Canale 36
      IF (JCH.EQ.36) THEN
         NUMPAR=NPLEV(JCH,1)
         NPART = NUMPAR
         P(1,1) = PLAB1(1,1)
         P(2,1) = PLAB1(2,1)
         P(3,1) = PLAB1(3,1)
         P(4,1) = AMASS1(1)
         I(1) = CODE(PARTIC(JCH,1,0,1))
         P(1,2) = PLAB1(1,2)
         P(2,2) = PLAB1(2,2)
         P(3,2) = PLAB1(3,2)
         P(4,2) = AMASS1(2)
         I(2) = CODE(PARTIC(JCH,1,0,2))
         IF (TARGET.EQ.'He3') THEN
            P(1,3) = PN(1,2)
            P(2,3) = PN(2,2)
            P(3,3) = PN(3,2)
            P(4,3) = RMASS(PARTIC(JCH,1,0,3))
            I(3) = CODE(PARTIC(JCH,1,0,3))
         END IF
         GOTO 99
      END IF





C ***    Canale 37
       IF (JCH.EQ.37) THEN
        NUMPAR=NPLEV(JCH,1)
        NPART = NUMPAR+NUMPLUS
        P(1,1) = PLAB1(1,1)
        P(2,1) = PLAB1(2,1)
        P(3,1) = PLAB1(3,1)
        P(4,1) = AMASS1(1)
        I(1) = CODE(PARTIC(JCH,1,0,1))
        P(1,2) = PLAB1(1,2)
        P(2,2) = PLAB1(2,2)
        P(3,2) = PLAB1(3,2)
        P(4,2) = AMASS1(2)
        I(2) = CODE(PARTIC(JCH,1,0,2))
        P(1,3) = PN(1,2)
        P(2,3) = PN(2,2)
        P(3,3) = PN(3,2)
        P(4,3) = RMASS(PARTIC(JCH,1,0,3))         
        I(3) = CODE(PARTIC(JCH,1,0,3))
        P(1,4) = PN(1,3)
        P(2,4) = PN(2,3)
        P(3,4) = PN(3,3)
        P(4,4) = RMASS(PARTIC(JCH,1,0,4))         
        I(4) = CODE(PARTIC(JCH,1,0,4))
        GOTO 99
       END IF

C ***    Canali 38,39
       IF (JCH.EQ.38.OR.JCH.EQ.39) THEN
        NUMPAR=NPLEV(JCH,1)
        NPART = NUMPAR+1
        IF (JCH.EQ.39) THEN
          KK=2
        ELSE
          KK=1
        END IF
        P(1,1) = PLAB1(1,1)
        P(2,1) = PLAB1(2,1)
        P(3,1) = PLAB1(3,1)
        P(4,1) = AMASS1(1)
        I(1) = CODE(PARTIC(JCH,1,0,1))
        IF (JCH.EQ.39) THEN
          P(1,2) = PLAB1(1,2)
          P(2,2) = PLAB1(2,2)
          P(3,2) = PLAB1(3,2)
          P(4,2) = AMASS1(2)
          I(2) = CODE(PARTIC(JCH,1,0,2))
          P(1,3) = PLAB1(1,3)
          P(2,3) = PLAB1(2,3)
          P(3,3) = PLAB1(3,3)
          P(4,3) = AMASS1(3)
          I(3) = CODE(PARTIC(JCH,1,0,3))
          P(1,4) = PN(1,2)
          P(2,4) = PN(2,2)
          P(3,4) = PN(3,2)
          P(4,4) = RMASS(PARTIC(JCH,1,0,4))         
          I(4) = CODE(PARTIC(JCH,1,0,4))
        ELSE
          P(1,2) = PLAB1(1,2)
          P(2,2) = PLAB1(2,2)
          P(3,2) = PLAB1(3,2)
          P(4,2) = AMASS1(2)
          I(2) = CODE(PARTIC(JCH,1,0,2))
          P(1,3) = PN(1,2)
          P(2,3) = PN(2,2)
          P(3,3) = PN(3,2)
          P(4,3) = RMASS(PARTIC(JCH,1,0,3))         
          I(3) = CODE(PARTIC(JCH,1,0,3))
        END IF
        GOTO 99
       END IF

C ***    Canale 40
       IF (JCH.EQ.40) THEN
        NPART = 4
        P(1,1) = PLAB1(1,1)
        P(2,1) = PLAB1(2,1)
        P(3,1) = PLAB1(3,1)
        P(4,1) = AMASS1(1)
        I(1) = CODE(PARTIC(JCH,1,0,1))
        P(1,2) = PLAB2(1,1)
        P(2,2) = PLAB2(2,1)
        P(3,2) = PLAB2(3,1)
        P(4,2) = AMASS2(1)
        I(2) = CODE(PARTIC(JCH,2,1,1))
        P(1,3) = PLAB2(1,2)
        P(2,3) = PLAB2(2,2)
        P(3,3) = PLAB2(3,2)
        P(4,3) = AMASS2(2)
        I(3) = CODE(PARTIC(JCH,2,1,2))
        P(1,4) = PN(1,2)
        P(2,4) = PN(2,2)
        P(3,4) = PN(3,2)
        P(4,4) = RMASS(PARTIC(JCH,1,0,3))         
        I(4) = CODE(PARTIC(JCH,1,0,3))
       END IF


 99   CONTINUE
!      CALL HFILL(JCH,EGAM,0.,1.)
      RETURN
      END                                                !*** END EV_NUC

      

	SUBROUTINE GLOREN(BETA,PA,PB)
C.
C.    ******************************************************************
C.    *                                                                *
C     *       Routine to transform momentum and energy from the        *
C     *       Lorentz frame A to the Lorentz frame B                   *
C     *                                                                *
C     *       PA(1)                                                    *
C     *       PA(2)     Momentum components in frame A                 *
C     *       PA(3)                                                    *
C     *       PA(4)     Energy                                         *
C     *       PB(..)   same quantities in frame B                      *
C     *                                                                *
C     *       BETA(1)    Components of velocity of frame B             *
C     *       BETA(2)        as seen from frame A                      *
C     *       BETA(3)                                                  *
C     *       BETA(4)    1./SQRT(1.-BETA**2)                           *
C.    *                                                                *
C.    *       Author    M.Hansroul  *********                          *
C.    *                                                                *
C.    ******************************************************************
C	RIA : taken from geant321
C.
	DIMENSION BETA(4),PA(4),PB(4)
C.
C.    ------------------------------------------------------------------
C.
	BETPA  = BETA(1)*PA(1) + BETA(2)*PA(2) + BETA(3)*PA(3)
	BPGAM  = (BETPA * BETA(4)/(BETA(4) + 1.) - PA(4)) * BETA(4)
	PB(1) = PA(1) + BPGAM  * BETA(1)
	PB(2) = PA(2) + BPGAM  * BETA(2)
	PB(3) = PA(3) + BPGAM  * BETA(3)
	PB(4) =(PA(4) - BETPA) * BETA(4)
	END
 
