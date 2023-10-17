
      SUBROUTINE CHANNELEX(E_N,E_DELTA,E_DEU,E_DIB,JCH)

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                                                            
C FUNCTIONAL DESCRIPTION:                                         
C                                                            
C    Select a channel to simulate according to the tabulated cross sections.
C
C    PROB(i)= probability for i-th selected channel
C
C    NCHAIN_P = number of selected reaction channel on proton
C    NCHAIN_N = number of selected reaction channel on neutron
C    NCHAIN_D = number of selected reaction channel on deuteron
c
C MODIFICATION DATE:
C 
C     12.05.95       by Igor Pshenichnov
C
C  ************   by marco mirazita
C
C     23.04.99
C --> Aggiunto l'azzeramento del vettore PROB(0:40), altrimenti non viene
C     calcolata correttamente la probabilita' dei vari canali di reazione.
C     Nella versione che girava su VMS questo probabilmente non era un 
C     problema perche' il vettore era comunque azzerato dal compilatore.
C     29.04.99
C --> Aggiunto il file di include PROCESS.INC con le variabili per il
C     conteggio dei canali di reazione selezionati
C --> Aggiunto il caso del canale 35 gamma d --> d rho0 (ed anche il canale
C     36 gamma d --> d eta, che non e' implementato completamente).
C --> Modificata la parte che interpola le sezioni d'urto totali per i
C     canali su d (era stata in precedenza modificata anche da federico).
C     05.05.99
C --> Modificata tutta la parte di estrazione del canale di reazione, che
C     aveva dei giri strani. Adesso e' possibile selezionare un numero
C     qualsiasi di canali di reazione.
C --> Modificato il calcolo di NSUM (numero di canali totali) e NSTART 
C     (indice del canale di partenza del conteggio), che ora sono calcolati
C     in base al numero di canali selezionati per ciascun tipo di bersaglio
C     (p, n, d, He3).
C --> Eliminata la variabile TARGFERMI. Adesso la selezione dei canali di 
C     reazione e' fatta solo in base a TARGET. I canali su p (n) vengono 
C     considerati solo se TARGET e' p (n), d oppure He3.
C
C COMMON BLOCKS:                                               
C                                                             
C     /MONTECARL/                                               
C     /FLAG_COMP/                                                     
C     /FLAG_PROC/                                               
C     /TARGET/                                                       
C                                                             
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     include "process.inc"

      INCLUDE 'process.inc'

      INTEGER IBEAM,IELET,NCHAIN,ICHAIN(40),ICIBLE,IPOL,ILAM,IARM
      REAL ETAG,FWHM
      COMMON/FLAG_COMP/IELET,IBEAM,ICIBLE,ETAG,FWHM,NCHAIN,ICHAIN,IPOL,
     &   ILAM
      REAL EGAM, P(4,11)
      INTEGER I(11),NPART
      COMMON/MONTECARL/ EGAM,I,NPART,P
      INTEGER ITARG,IIDH
      CHARACTER*12 TARGET,TARGFERMI
      COMMON /TARGET/ ITARG,IIDH,TARGET,TARGFERMI
      REAL WR(263),SIGR(40,263)
      COMMON/SIGTO/WR,SIGR
      INTEGER IPRL
      COMMON/IPRL/IPRL
      INTEGER NPRO(40),NDEU,NHE3,NSUM,NSTART,M,N,ICHANNEL,K,JCH,J
      INTEGER IPOS,IDIHOTOMIA
      REAL SIGTOTRIF(2),XX(2),SIGREL(40),PROB(0:40),
     &           E_N,E_DELTA,E_DEU,E_DIB,ECM,RNDA
      REAL RNDM,ALI

      real vran(1)

C ***  Initialization
      DO J=1,40
         NPRO(J) = 0
      END DO
      DO J=1,40
         SIGREL(J) = 0.
      END DO
      NDEU = 0
      NHE3 = 0
      NSUM = 0
      NSTART = 0
C *****************************************
C *** Initialization of the probabilities of reaction channels.
      do j=0,40
         prob(j)=0.
      enddo
C *****************************************
C ***   Different channels separation
C ************************************
C *** Se non ci sono canali selezionati, li prendo tutti
      IF (NCHAIN.EQ.0) THEN
         DO J=1,40
            ICHAIN(J) = NPRO(J)
         END DO
      END IF
C ************************************

C ************************************
C *** Find energy bin for proton or neutron
C ************************************
      IPOS = IDIHOTOMIA(E_N,WR,263)

C ************************************
C *** Free or bound proton
C ********************************************
      IF ( (TARGet.EQ.'proton      ').or.
     +     (TARGet.EQ.'deuteron    ').or.
     +     (TARGet.EQ.'He3      ') ) THEN
         IF (NCHAIN.EQ.0) then
            NSUM = nsum + nprop_tot
         else
            NSUM = nsum + NCHAIN_p
         endif
         DO J=nstart+1,nsum
            NPRO(J) = NPROP(J-nstart)
         END DO

         DO N=nstart+1,nsum

            ICHANNEL = NPRO(N)

C *** Linear interpolation in the energy bin
            SIGTOTRIF(1) = SIGR(ICHANNEL,IPOS-1)
            XX(1) = WR(IPOS-1)
C *** IPOS deve essere < 264
            if (ipos.le.263) then
               SIGTOTRIF(2) = SIGR(ICHANNEL,IPOS)
               XX(2) = WR(IPOS)
            else
               SIGTOTRIF(2) = SIGR(ICHANNEL,263)
               XX(2) = WR(263)
            endif


c            SIGTOTRIF(2) = SIGR(ICHANNEL,IPOS)
c            XX(2) = WR(IPOS)
            SIGREL(n) = ALI(E_N,XX,SIGTOTRIF,2)
c            write (15,*) 'N=',n,' xx=',xx
c            write (15,*) 'N=',n,' sigmarif=',sigtotrif
c            write (15,*) '  W=',e_n,' sigma=',sigrel(n)
         END DO

         nstart=nstart+nchain_p
      END IF

C ************************************
C *** Free or bound neutron
C ********************************************
      IF ( (TARGet.EQ.'neutron      ').or.
     +     (TARGet.EQ.'deuteron    ').or.
     +     (TARGet.EQ.'He3      ') ) THEN
         IF (NCHAIN.EQ.0) then
            NSUM = nsum + npron_tot
         else
            NSUM = nsum + NCHAIN_n
         endif
         DO J=nstart+1,nsum
            NPRO(J) = NPRON(J-nstart)
         END DO

c         write (15,*) 'NEUTRONE'
         DO N=nstart+1,nsum

            ICHANNEL = NPRO(N)

C *** Linear interpolation in the energy bin
            SIGTOTRIF(1) = SIGR(ICHANNEL,IPOS-1)
            XX(1) = WR(IPOS-1)

C *** IPOS deve essere < 264
            if (ipos.le.263) then
               SIGTOTRIF(2) = SIGR(ICHANNEL,IPOS)
               XX(2) = WR(IPOS)
            else
               SIGTOTRIF(2) = SIGR(ICHANNEL,263)
               XX(2) = WR(263)
            endif

c            SIGTOTRIF(2) = SIGR(ICHANNEL,IPOS)
c            XX(2) = WR(IPOS)
            SIGREL(n) = ALI(E_N,XX,SIGTOTRIF,2)
c            write (15,*) 'N=',n,' xx=',xx
c            write (15,*) 'N=',n,' sigmarif=',sigtotrif
c            write (15,*) '  W=',e_n,' sigma=',sigrel(n)

         END DO


         nstart=nstart+nchain_n
      END IF


C ********************************************
C *** Speciali deuteron channels
C ********************************************
      IF (TARGET.EQ.'deuteron'.AND.NCHAIN_D.GT.0) THEN
          IF (NCHAIN.EQ.0) then
             NSUM = nsum + nprod_tot
          else
             NSUM = nsum + NCHAIN_d
          endif
C *** Aggiungo a NPRO i canali su deutone
         DO J=nstart+1,nsum
            NPRO(J) = NPROD(J-nstart)
         END DO


C *** Select energy bin
         IPOS = IDIHOTOMIA(e_deu,WR,263) ! FIND AN ENERGY INTERVAL 
                                       ! TO WORK ON

C *** Select channel cross sections
         DO M=NSTART+1,Nsum

            ICHANNEL = NPRO(m)

C *** Linear interpolation in the energy bin
            SIGTOTRIF(1) = SIGR(ICHANNEL,IPOS-1)
            XX(1) = WR(IPOS-1)

C *** IPOS deve essere < 264
            if (ipos.le.263) then
               SIGTOTRIF(2) = SIGR(ICHANNEL,IPOS)
               XX(2) = WR(IPOS)
            else
               SIGTOTRIF(2) = SIGR(ICHANNEL,263)
               XX(2) = WR(263)
            endif

            SIGREL(M) = ALI(E_Deu,XX,SIGTOTRIF,2)
         END DO
         nstart=nstart+nchain_n

      END IF

C ********************************************
C *** Special He3 channels - not implemented
C ********************************************

      IF (TARGET.EQ.'He3'.AND.NCHAIN_HE3.GT.0) THEN
          IF (NCHAIN.EQ.0) then
             NSUM = nsum + nprohe3_tot
          else
             NSUM = nsum + NCHAIN_he3
          endif
C *** Aggiungo a NPRO i canali du He3
         DO J=nstart+1,nsum
            NPRO(J) = NPROHE3(J-nstart)
         END DO

         DO M=NSTART+1,NSum

 4          ICHANNEL = NPRO(m)

C *** Select target energy 
            IF (ICHANNEL.EQ.37) ECM = E_DELTA
            IF (ICHANNEL.EQ.38.OR.ICHANNEL.EQ.39) ECM = E_DEU
            IF (ICHANNEL.EQ.40) ECM = E_DIB

            IPOS = IDIHOTOMIA(Ecm,WR,263) ! FIND AN ENERGY INTERVAL 
                                       ! TO WORK ON

            SIGTOTRIF(1) = SIGR(ICHANNEL,IPOS-1)
            XX(1) = WR(IPOS-1)
            SIGTOTRIF(2) = SIGR(ICHANNEL,IPOS)
            XX(2) = WR(IPOS)
            SIGREL(M) = ALI(Ecm,XX,SIGTOTRIF,2) ! MAKE A LINEAR INTERPOLATION
                                                ! IN THIS INTERVAL ONLY
         END DO
         nstart=nstart+nchain_he3
      END IF


      DO N=1,nsum
         PROB(N) = PROB(N-1)+SIGREL(N)
      END DO

      DO N=1,Nsum
         PROB(N) = PROB(N)/PROB(nsum)
      END DO

      call ranmar(vran,1)
      rnda=vran(1)

      DO J=1,NSum
         IF (RNDA.GE.PROB(J-1).AND.RNDA.LT.PROB(J)) THEN
            JCH = npro(j)
            ICHANNEL = JCH
            RETURN
         END IF
      END DO


      RETURN
      END                       !*** END CHANNELEX

C RIA : moved from graph1.F

            INTEGER FUNCTION IDIHOTOMIA(TRY,ARRAY,NA)
C
C                FIND THE NUMBER ID SATISFIED:
C                ARRAY(ID-1).LT.TRY.LE.ARRAY(ID)
C                
C

         REAL*4 ARRAY(NA)
         IF(TRY.LE.ARRAY(1)) THEN
           IDIHOTOMIA=1
         ELSE IF(TRY.GT.ARRAY(NA)) THEN
           IDIHOTOMIA=NA+1
         ELSE
                             NBEGIN=1
                             NEND=NA
                      NRANGE=(NEND-NBEGIN)/2+NBEGIN
                 DO WHILE ((NEND-NBEGIN).GT.1)
                    IF(TRY.GT.ARRAY(NRANGE)) THEN
                        NBEGIN=NRANGE
                    ELSE
                        NEND=NRANGE
                    END IF
                    NRANGE=(NEND-NBEGIN)/2+NBEGIN
                 END DO
         IDIHOTOMIA=NEND
         END IF
         RETURN 
         END         


                 REAL*4 FUNCTION ALI(X,T,B,N)
C  
C          MAKE A LINEAR INTERPOLATION OF THE B(N) ARRAY 
C                     AT X ABSCISSA VALUE
C
                 REAL*4 B(N),T(N)

        IK=IDIHOTOMIA(X,T,N)
        IF (IK.EQ.1) THEN
           ALI=B(1)
           RETURN
        ELSE IF (IK.GT.N) THEN
           ALI=B(N)
           RETURN
        ELSE 
        ALI=((B(IK)-B(IK-1))/(T(IK)-T(IK-1)))*(X-T(IK-1))+B(IK-1)
        END IF
        RETURN
        END
