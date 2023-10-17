	SUBROUTINE GENBOS_REACTIONS(NREAC, IREAC)
	implicit none
	INTEGER NREAC, IREAC(40)
      
	INTEGER IBEAM,IELET,NCHAIN,ICHAIN(40),ICIBLE,IPOL,ILAM
	REAL ETAG,FWHM,the_min,the_max,dphi
	COMMON/FLAG_COMP/IELET,IBEAM,ICIBLE,ETAG,FWHM,NCHAIN,ICHAIN,IPOL,
     &  	ILAM
	INTEGER I

	IF (NREAC.LT.40) NREAC=40
	DO I=1,40
	  IF (I.LE.NREAC) THEN
		ICHAIN(I) = IREAC(I)
	  ELSE 
		ICHAIN(I) = 0
	  ENDIF
	ENDDO
	
	END

C-----
	SUBROUTINE GENBOS_TARG(TG)
	implicit none
	REAL TG
      
	INTEGER IBEAM,IELET,NCHAIN,ICHAIN(40),ICIBLE,IPOL,ILAM
	REAL ETAG,FWHM
	COMMON/FLAG_COMP/IELET,IBEAM,ICIBLE,ETAG,FWHM,NCHAIN,ICHAIN,IPOL,
     &  	ILAM
	

	ICIBLE = TG
	
	END

	SUBROUTINE GENBOS_BEAM(BE, EG, EW)
	implicit none
	REAL BE, EG, EW
      
	INTEGER IBEAM,IELET,NCHAIN,ICHAIN(40),ICIBLE,IPOL,ILAM
	REAL ETAG,FWHM
	COMMON/FLAG_COMP/IELET,IBEAM,ICIBLE,ETAG,FWHM,NCHAIN,ICHAIN,IPOL,
     &  	ILAM
	

	IBEAM = BE
	ETAG = EG
	FWHM = EW
	
	END

C-----
	SUBROUTINE GENBOS_START
	implicit none
	INTEGER IX,IA,IB
	DATA IX,IA,IB/76762,1000,0/
	
	INTEGER IBEAM,IELET,NCHAIN,ICHAIN(40),ICIBLE,IPOL,ILAM
	REAL ETAG,FWHM,the_min,the_max,dphi
	COMMON/FLAG_COMP/IELET,IBEAM,ICIBLE,ETAG,FWHM,NCHAIN,ICHAIN,IPOL,
     &  	ILAM
	CHARACTER*50 NCANAL(40)
	COMMON/CLABEL/NCANAL
	INTEGER NEV
	COMMON/COUNT/NEV
     
	INTEGER K
     
C *** open log file
	open (unit=15,FILE='genbos.log',STATUS='unknown')

        CALL CARDS
        CALL GEN_INIT

C *** genbos inizialization
        
C *** random generator inizialization
        CALL RMARIN(IX,IA,IB)
        
        NEV=0;

	write (15,*) 'Random Seed:',ix
	write (15,*)

	if (ibeam.eq.0) then
	   write (15,62) etag,fwhm
	else if (ibeam.eq.2) then
	   write (15,63) etag+fwhm,etag-fwhm
	else if (ibeam.eq.3) then
	   write (15,64) etag,fwhm
	endif

	if (icible.eq.0) then
	   write (15,65) ' Neutron target'
	else if (icible.eq.1) then
	   write (15,65) ' Proton target'
	else if (icible.eq.2) then
	   write (15,65) ' Deuteron target'
	endif

	write (15,65) ' Laboratory system for outgoing particles'


 62	format (' Gaussian spectrum: E=',F10.3,' +/- ',F10.3,/)
 63	format (' Bremsstrahlun spectrum: Emax=',F10.3,' Emin=',F10.3,/)
 64	format (' Uniform spectrum: E=',F10.3,' +/- ',F10.3,/)
 65	format (a,/)
 66	format (' channel n.',i3,3x,a)

	
	END

C-----
	SUBROUTINE GENBOS_STOP
	implicit none
	INTEGER NEV
	COMMON/COUNTER/NEV

	WRITE(15,10) NEV
10	FORMAT (' TOTAL EVENTS PROCESSED : ',I10)
	CLOSE (15)
	
	END
	
C-----
	SUBROUTINE GENBOS_EVENT(EFOT,NREAC,NP,ID,PX,PY,PZ)
	implicit none
	REAL EFOT, PX(11),PY(11),PZ(11)
	INTEGER NREAC,NP,ID(11)
	
	COMMON/COUNTER/NEV
	INTEGER NEV
        COMMON/REACCH/J_CHANNEL
        INTEGER J_CHANNEL
	COMMON/MONTECARL/ EGAM,I,NPART,P
	REAL EGAM, P(4,11)
	INTEGER I(11),NPART

	INTEGER K

	CALL MONTECARLO
	
	NREAC = J_CHANNEL
	NP = NPART
	EFOT = EGAM
	
	IF (NP.EQ.0) RETURN
	
	DO K=1,NPART
	  ID(K) = I(K)
	  PX(K) = P(1,K)
	  PY(K) = P(2,K)
	  PZ(K) = P(3,K)	
	ENDDO
	
	NEV = NEV+1
	
	END

C-----
	SUBROUTINE GENBOS_RAND(IX)
	implicit none
	INTEGER IX
	INTEGER IA,IB
	DATA IA,IB/1000,0/
C *** random generator inizialization
        CALL RMARIN(IX,IA,IB)
        
	write (15,*) 'New Random Seed:',ix
	write (15,*)

	END
	
C***********************************************************************

      SUBROUTINE CARDS  
      
#include "process.inc"
#include "nt_kine.inc"

      INTEGER ITARG,IIDH
      CHARACTER*12 TARGET,TARGFERMI
      COMMON /TARGET/ ITARG,IIDH,TARGET,TARGFERMI
      INTEGER IBEAM,IELET,NCHAIN,ICHAIN(40),ICIBLE,IPOL,ILAM
      REAL ETAG,FWHM
      COMMON/FLAG_COMP/IELET,IBEAM,ICIBLE,ETAG,FWHM,NCHAIN,ICHAIN,IPOL,
     &   ILAM
      COMMON/NSTATIS/NSTATIS      
      LOGICAL EG,CMSF
      COMMON/LOGEG/EG,CMSF

	DATA ICIBLE,IBEAM /2,2/
	DATA ETAG,FWHM /0.8,0.2/
	DATA ICHAIN /1,2,3,4,5,6,7,8,9,10,15,16,17,18,33,34,24*0/ 
	
	DATA NSTATIS,IRUN /500,0/


	ielet=0
	EG=.TRUE.    ! NUMBERING AS IN INC (.FALSE.) OR AS IN GEANT (.TRUE.)
	CMSF=.FALSE. ! All the momenta of particles are in CMS (.TRUE.)
                    ! or in LAB system (.FALSE.)

       IF (ICIBLE.EQ.0) THEN   ! NEUTRON AS A TARGET
         TARGET = 'neutron     '
       ENDIF

       IF (ICIBLE.EQ.1) THEN   ! PROTON AS A TARGET
         TARGET = 'proton      '
       ENDIF

       IF (ICIBLE.EQ.2) THEN   ! DEUTERON
         TARGET = 'deuteron    '
       ENDIF

 
C
C ***      Count processes selected for simulation
C
      NCHAIN = 0

      nchain_p=0
      nchain_n=0
      nchain_d=0
      nchain_he3=0

      ILLA  = 0
      DO I = 1,40
         IF (ICHAIN(I).GT.0 ) NCHAIN = NCHAIN + 1
         IF (ICHAIN(I).EQ.29) ILLA  = 1
C *********************************
C *** Proton
         if ( (ichain(i).eq.1).or.(ichain(i).eq.2).or.
     +        ((ichain(i).ge.5).and.(ichain(i).le.7)).or.
     +        (ichain(i).eq.11).or.(ichain(i).eq.12).or.
     +        (ichain(i).eq.15).or.(ichain(i).eq.16).or.
     +        (ichain(i).eq.19).or.(ichain(i).eq.21).or.
     +        (ichain(i).eq.23).or.(ichain(i).eq.24).or.
     +        (ichain(i).eq.27).or.
     +        ((ichain(i).ge.29).and.(ichain(i).le.32))
     *                                ) then
            nchain_p=nchain_p+1
            nprop(nchain_p)=ichain(i)

C *** Neutron
         else if ( (ichain(i).eq.3).or.(ichain(i).eq.4).or.
     +        ((ichain(i).ge.8).and.(ichain(i).le.10)).or.
     +        (ichain(i).eq.13).or.(ichain(i).eq.14).or.
     +        (ichain(i).eq.17).or.(ichain(i).eq.18).or.
     +        (ichain(i).eq.20).or.(ichain(i).eq.22).or.
     +        (ichain(i).eq.25).or.(ichain(i).eq.26).or.
     +        (ichain(i).eq.28)
     *                                ) then
            nchain_n=nchain_n+1
            npron(nchain_n)=ichain(i)

C *** Deuteron
         else if ( (ichain(i).ge.33).and.(ichain(i).le.36)
     +                                ) then
            nchain_d=nchain_d+1
            nprod(nchain_d)=ichain(i)
C *** He3
         else if ( (ichain(i).ge.37).and.(ichain(i).le.40)
     +                                ) then
            nchain_he3=nchain_he3+1
            nprohe3(nchain_he3)=ichain(i)

         endif
      ENDDO
C *** Check on selected channels
      if (icible.ne.2) nchain_d=0
      if (icible.ne.3) nchain_he3=0
      if (icible.eq.0) nchain_p=0
      if (icible.eq.1) nchain_n=0

C      IF (ILLA.NE.1.OR.NCHAIN.GT.1) ILAM  = 0
      ILAM = 0

      RETURN
      END


C***********************************************************************


      SUBROUTINE GEN_INIT
C+
C FUNCTIONAL DESCRIPTION:
C 
C ***   Initialisation of event generator
C
C MODIFICATION DATE:
C 
C     19.05.95      by Igor Pshenichnov
C 
C-
      IMPLICIT NONE
      INTEGER ITARG,IIDH
      CHARACTER*12 TARGET,TARGFERMI
      COMMON /TARGET/ ITARG,IIDH,TARGET,TARGFERMI
      INTEGER IBEAM,IELET,NCHAIN,ICHAIN(40),ICIBLE,IPOL,ILAM
      REAL ETAG,FWHM
      COMMON/FLAG_COMP/IELET,IBEAM,ICIBLE,ETAG,FWHM,NCHAIN,ICHAIN,IPOL,
     &   ILAM
      CHARACTER*20 FILEN,FILEP
      COMMON/FNAME/FILEN,FILEP
      INTEGER JCH,INW,INTH,IE
      REAL XSECT(40,50,0:18)
      COMMON /SIGMA/ XSECT,IE

      IF (IELET.EQ.0) THEN

C ***   Put the names of particles for each channel by the next call
         CALL CHANNELDEF

C ***   Total cross section for multiplication with energy spectrum,
C ***   Verification of reaction channel
         CALL TOTALSECT
      ELSE
         WRITE(15,*) '    Electromagnetic events'
      ENDIF

2323  FORMAT (1X)
   7  FORMAT (20X,5G10.3)
   8  FORMAT (7G10.3)
      RETURN
      END                                           !*** END GEN_INIT



C***********************************************************************


C***********************************************************************

      SUBROUTINE CHANNELDEF

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                                                    
C ABSTRACT:                                             
C                                                     
C     This subroutine defines particles for different reaction channels,
C     branching ratios and final states from decays
C                                                     
C
C     PARTIC(JCH,LEV,JBR,IP) = particle (char)
C
C     JCH = channel number
C     LEV = level number (0 = target particle)
C     JBR = number of decay channel
C     IP  = progressiv particle number
C
C     JDEC(JCH) = number of decaying particle at first level 
C                 for channel JCH
C     NCH(LEV) = number of decay channel for level LEV
C                                                   
C     BR(JCH,NCH(LEV)) = branching ratio for reaction JCH and
C                        decay decay channel NCH(LEV)
C     NPLEV(JCH,LEV) = number of produced particles at level LEV
C                      from reaction channel JCH 
C
C     List of GEANT
C                                                   
C      protone		14                                 
C      neutrone		13                                 
C      deuteron		45                                 
C      pai +		8                                 
C      pai 0		7                                 
C      pai -		9                                 
C      eta		17                                 
C      lambda		18                                 
C      K +		11                                 
C      K 0 short	16                                 
C      K 0 long		10                                 
C      sigma 0		20                                 
C                                                   
C      delta ++		101                                 
C      delta +		102                                 
C      delta 0		103                                 
C      delta -		104                                 
C      rho +		58                                 
C      rho 0		57                                 
C      rho -		59                                 
C      omega		60                                 
C                                                   
C      He3		49                                
C                                                   
C COMMON BLOCKS:                                       
C                                                     
C     /CHANNELS/                                       
C     /CLABEL/
C                                                     
C     MODIFICATIONS
C
C   **********    by marco mirazita
C
C     05.05.99
C --> Aggiunto il canale di reazione 35 gamma d --> rho0 d e 36 
C     gamma d --> rho0 eta.
C-----------------------------------------------------------------------

      IMPLICIT NONE
      CHARACTER*12 PARTIC(40,0:3,0:3,8)
      INTEGER JDEC(40),NCH(40),NPLEV(40,3)
      REAL BR(40,3)
      COMMON /CHANNELS/ PARTIC,JDEC,NCH,NPLEV,BR

      INTEGER I,J,K,L,J1

      CHARACTER*50 NCANAL(40)
      COMMON/CLABEL/NCANAL

      NCANAL(1)='gamma + p --->  pi+  n'
      NCANAL(2)='gamma + p --->  pi0  p'
      NCANAL(3)='gamma + n --->  pi-  p'
      NCANAL(4)='gamma + n --->  pi0  n'
      NCANAL(5)='gamma + p --->  del++ pi-'
      NCANAL(6)='gamma + p --->  del+  pi0'
      NCANAL(7)='gamma + p --->  del0  pi+'
      NCANAL(8)='gamma + n --->  del+  pi-'
      NCANAL(9)='gamma + n --->  del0  pi0'
      NCANAL(10)='gamma + n --->  del-  pi+'
      NCANAL(11)='gamma + p --->  rho0  p'
      NCANAL(12)='gamma + p --->  rho+  n'
      NCANAL(13)='gamma + n --->  rho-  p'
      NCANAL(14)='gamma + n --->  rho0  n'
      NCANAL(15)='gamma + p --->  pi+  pi-  p'
      NCANAL(16)='gamma + p --->  pi0  pi+  n'
      NCANAL(17)='gamma + n --->  pi+  pi-  n'
      NCANAL(18)='gamma + n --->  pi0  pi-  p'
      NCANAL(19)='gamma + p --->  eta  p'
      NCANAL(20)='gamma + n --->  eta  n'
      NCANAL(21)='gamma + p --->  ome  p'
      NCANAL(22)='gamma + n --->  ome  n'
      NCANAL(23)='gamma + p --->  pi+  pi-  pi0  p'
      NCANAL(24)='gamma + p ---> 3pi0 p. pi+ 2pi0 n. 2pi+ pi- n'
      NCANAL(25)='gamma + n --->  pi+  pi-  pi0  n'
      NCANAL(26)='gamma + n ---> 3pi0 n. pi- 2pi0 p. pi+ 2pi- p'
      NCANAL(27)='gamma + p --->  (4-8)pi N'
      NCANAL(28)='gamma + n --->  (4-8)pi N'
      NCANAL(29)='gamma + p --->  lam  K+'
      NCANAL(30)='gamma + p --->  sig  K+'
      NCANAL(31)='gamma + p --->  lam  K0  pi+'
      NCANAL(32)='gamma + p --->  gamma + p'
      NCANAL(33)='gamma + D --->  pi0  D'
      NCANAL(34)='gamma + D --->  p  n'
      NCANAL(35)='gamma + D --->  rho0  D'
      NCANAL(36)='gamma + D --->  eta  D'
      NCANAL(37)='gamma + He3 --->  delta++ knockout'
      NCANAL(38)='gamma + He3 --->  pi0  D  p'
      NCANAL(39)='gamma + He3 --->  pi+  pi-  D  p'
      NCANAL(40)='gamma + He3 --->  pi-  dib  p'

C ***   Inizializzazioni
        DO I=1,40
            JDEC(I) = 0
            NCH(I) = 0
                    DO J1=1,3
                        NPLEV(I,J1) = 0
                        BR(I,J1) = 0.
                    ENDDO
            DO J=0,3
                DO K=0,3
                    DO L=1,8
                        PARTIC(I,J,K,L) = ' '
                    END DO
                END DO
            END DO
        END DO

C ***   Canale 1
      PARTIC(1,0,0,1)='proton'
      PARTIC(1,1,0,1)='pi+'
      PARTIC(1,1,0,2)='neutron'
      PARTIC(1,1,0,3)='neutron'
      PARTIC(1,1,0,4)='proton'
      JDEC(1)=0
      NCH(1)=0
      NPLEV(1,1)=2

C ***   Canale 2
      PARTIC(2,0,0,1)='proton'
      PARTIC(2,1,0,1)='pi0'
      PARTIC(2,1,0,2)='proton'
      PARTIC(2,1,0,3)='neutron'
      PARTIC(2,1,0,4)='proton'
      JDEC(2)=0
      NCH(2)=0
      NPLEV(2,1)=2

C ***   Canale 3
      PARTIC(3,0,0,1)='neutron'
      PARTIC(3,1,0,1)='pi-'
      PARTIC(3,1,0,2)='proton'
      PARTIC(3,1,0,3)='proton'
      PARTIC(3,1,0,4)='proton'
      JDEC(3)=0
      NCH(3)=0
      NPLEV(3,1)=2

C ***   Canale 4
      PARTIC(4,0,0,1)='neutron'
      PARTIC(4,1,0,1)='pi0'
      PARTIC(4,1,0,2)='neutron'
      PARTIC(4,1,0,3)='proton'
      PARTIC(4,1,0,4)='proton'
      JDEC(4)=0
      NCH(4)=0
      NPLEV(4,1)=2

C ***   Canale 5
      PARTIC(5,0,0,1)='proton'
      PARTIC(5,1,0,1)='pi-'
      PARTIC(5,1,0,2)='delta++'
      PARTIC(5,1,0,3)='neutron'
      PARTIC(5,1,0,4)='proton'
      JDEC(5)=2
      NCH(5)=1
      BR(5,1)=1.
      PARTIC(5,2,1,1)='pi+'
      PARTIC(5,2,1,2)='proton'
      NPLEV(5,1)=2
      NPLEV(5,2)=2

C ***   Canale 6
      PARTIC(6,0,0,1)='proton'
      PARTIC(6,1,0,1)='pi0'
      PARTIC(6,1,0,2)='delta+'
      PARTIC(6,1,0,3)='neutron'
      PARTIC(6,1,0,4)='proton'
      JDEC(6)=2
      NCH(6)=2
      BR(6,1)=0.33
      BR(6,2)=0.67
      PARTIC(6,2,1,1)='pi+'
      PARTIC(6,2,1,2)='neutron'
      PARTIC(6,2,2,1)='pi0'
      PARTIC(6,2,2,2)='proton'
      NPLEV(6,1)=2
      NPLEV(6,2)=2

C ***   Canale 7
      PARTIC(7,0,0,1)='proton'
      PARTIC(7,1,0,1)='pi+'
      PARTIC(7,1,0,2)='delta0'
      PARTIC(7,1,0,3)='neutron'
      PARTIC(7,1,0,4)='proton'
      JDEC(7)=2
      NCH(7)=2
      BR(7,1)=0.33
      BR(7,2)=0.67
      PARTIC(7,2,1,1)='pi-'
      PARTIC(7,2,1,2)='proton'
      PARTIC(7,2,2,1)='pi0'
      PARTIC(7,2,2,2)='neutron'
      NPLEV(7,1)=2
      NPLEV(7,2)=2

C ***   Canale 8
      PARTIC(8,0,0,1)='neutron'
      PARTIC(8,1,0,1)='pi-'
      PARTIC(8,1,0,2)='delta+'
      PARTIC(8,1,0,3)='proton'
      PARTIC(8,1,0,4)='proton'
      JDEC(8)=2
      NCH(8)=2
      BR(8,1)=0.33
      BR(8,2)=0.67
      PARTIC(8,2,1,1)='pi+'
      PARTIC(8,2,1,2)='neutron'
      PARTIC(8,2,2,1)='pi0'
      PARTIC(8,2,2,2)='proton'
      NPLEV(8,1)=2
      NPLEV(8,2)=2

C ***   Canale 9
      PARTIC(9,0,0,1)='neutron'
      PARTIC(9,1,0,1)='pi0'
      PARTIC(9,1,0,2)='delta0'
      PARTIC(9,1,0,3)='proton'
      PARTIC(9,1,0,4)='proton'
      JDEC(9)=2
      NCH(9)=2
      BR(9,1)=0.33
      BR(9,2)=0.67
      PARTIC(9,2,1,1)='pi-'
      PARTIC(9,2,1,2)='proton'
      PARTIC(9,2,2,1)='pi0'
      PARTIC(9,2,2,2)='neutron'
      NPLEV(9,1)=2
      NPLEV(9,2)=2

C ***   Canale 10
      PARTIC(10,0,0,1)='neutron'
      PARTIC(10,1,0,1)='pi+'
      PARTIC(10,1,0,2)='delta-'
      PARTIC(10,1,0,3)='proton'
      PARTIC(10,1,0,4)='proton'
      JDEC(10)=2
      NCH(10)=1
      BR(10,1)=1.
      PARTIC(10,2,1,1)='pi-'
      PARTIC(10,2,1,2)='neutron'
      NPLEV(10,1)=2
      NPLEV(10,2)=2

C ***   Canale 11
      PARTIC(11,0,0,1)='proton'
      PARTIC(11,1,0,1)='rho0'
      PARTIC(11,1,0,2)='proton'
      PARTIC(11,1,0,3)='neutron'
      PARTIC(11,1,0,4)='proton'
      JDEC(11)=1
      NCH(11)=1
      BR(11,1)=1.
      PARTIC(11,2,1,1)='pi+'
      PARTIC(11,2,1,2)='pi-'
      NPLEV(11,1)=2
      NPLEV(11,2)=2

C ***   Canale 12
      PARTIC(12,0,0,1)='proton'
      PARTIC(12,1,0,1)='rho+'
      PARTIC(12,1,0,2)='neutron'
      PARTIC(12,1,0,3)='neutron'
      PARTIC(12,1,0,4)='proton'
      JDEC(12)=1
      NCH(12)=1
      BR(12,1)=1.
      PARTIC(12,2,1,1)='pi+'
      PARTIC(12,2,1,2)='pi0'
      NPLEV(12,1)=2
      NPLEV(12,2)=2

C ***   Canale 13
      PARTIC(13,0,0,1)='neutron'
      PARTIC(13,1,0,1)='rho-'
      PARTIC(13,1,0,2)='proton'
      PARTIC(13,1,0,3)='proton'
      PARTIC(13,1,0,4)='proton'
      JDEC(13)=1
      NCH(13)=1
      BR(13,1)=1.
      PARTIC(13,2,1,1)='pi-'
      PARTIC(13,2,1,2)='pi0'
      NPLEV(13,1)=2
      NPLEV(13,2)=2

C ***   Canale 14
      PARTIC(14,0,0,1)='neutron'
      PARTIC(14,1,0,1)='rho0'
      PARTIC(14,1,0,2)='neutron'
      PARTIC(14,1,0,3)='proton'
      PARTIC(14,1,0,4)='proton'
      JDEC(14)=1
      NCH(14)=1
      BR(14,1)=1.
      PARTIC(14,2,1,1)='pi+'
      PARTIC(14,2,1,2)='pi-'
      NPLEV(14,1)=2
      NPLEV(14,2)=2

C ***   Canale 15
      PARTIC(15,0,0,1)='proton'
      PARTIC(15,1,0,1)='pi+'
      PARTIC(15,1,0,2)='pi-'
      PARTIC(15,1,0,3)='proton'
      PARTIC(15,1,0,4)='neutron'
      PARTIC(15,1,0,5)='proton'
      JDEC(15)=0
      NCH(15)=0
      NPLEV(15,1)=3

C ***   Canale 16
      PARTIC(16,0,0,1)='proton'
      PARTIC(16,1,0,1)='pi+'
      PARTIC(16,1,0,2)='pi0'
      PARTIC(16,1,0,3)='neutron'
      PARTIC(16,1,0,4)='neutron'
      PARTIC(16,1,0,5)='proton'
      JDEC(16)=0
      NCH(16)=0
      NPLEV(16,1)=3

C ***   Canale 17
      PARTIC(17,0,0,1)='neutron'
      PARTIC(17,1,0,1)='pi+'
      PARTIC(17,1,0,2)='pi-'
      PARTIC(17,1,0,3)='neutron'
      PARTIC(17,1,0,4)='proton'
      PARTIC(17,1,0,5)='proton'
      JDEC(17)=0
      NCH(17)=0
      NPLEV(17,1)=3

C ***   Canale 18
      PARTIC(18,0,0,1)='neutron'
      PARTIC(18,1,0,1)='pi-'
      PARTIC(18,1,0,2)='pi0'
      PARTIC(18,1,0,3)='proton'
      PARTIC(18,1,0,4)='proton'
      PARTIC(18,1,0,5)='proton'
      JDEC(18)=0
      NCH(18)=0
      NPLEV(18,1)=3

C ***   Canale 19
      PARTIC(19,0,0,1)='proton'
      PARTIC(19,1,0,1)='eta'
      PARTIC(19,1,0,2)='proton'
      PARTIC(19,1,0,3)='neutron'
      PARTIC(19,1,0,4)='proton'
      JDEC(19)=0
      NCH(19)=0
      NPLEV(19,1)=2

C ***   Canale 20
      PARTIC(20,0,0,1)='neutron'
      PARTIC(20,1,0,1)='eta'
      PARTIC(20,1,0,2)='neutron'
      PARTIC(20,1,0,3)='proton'
      PARTIC(20,1,0,4)='proton'
      JDEC(20)=0
      NCH(20)=0
      NPLEV(20,1)=2

C ***   Canale 21
      PARTIC(21,0,0,1)='proton'
      PARTIC(21,1,0,1)='omega'
      PARTIC(21,1,0,2)='proton'
      PARTIC(21,1,0,3)='neutron'
      PARTIC(21,1,0,4)='proton'
      JDEC(21)=0
      NCH(21)=0
      NPLEV(21,1)=2

C ***   Canale 22
      PARTIC(22,0,0,1)='neutron'
      PARTIC(22,1,0,1)='omega'
      PARTIC(22,1,0,2)='neutron'
      PARTIC(22,1,0,3)='proton'
      PARTIC(22,1,0,4)='proton'
      JDEC(22)=0
      NCH(22)=0
      NPLEV(22,1)=2

C ***   Canale 23
      PARTIC(23,0,0,1)='proton'
      PARTIC(23,1,0,1)='pi+'
      PARTIC(23,1,0,2)='pi-'
      PARTIC(23,1,0,3)='pi0'
      PARTIC(23,1,0,4)='proton'
      PARTIC(23,1,0,5)='neutron'
      PARTIC(23,1,0,6)='proton'
      JDEC(23)=0
      NCH(23)=0
      NPLEV(23,1)=4

C ***   Channel #24 includes p 3pi0, n pi+ 2pi0, n 2pi+ pi-,
C ***   n 2pi+ pi- presented here as a placeholder only
      PARTIC(24,0,0,1)='proton'
      PARTIC(24,1,0,1)='pi+'
      PARTIC(24,1,0,2)='pi+'
      PARTIC(24,1,0,3)='pi-'
      PARTIC(24,1,0,4)='neutron'
      PARTIC(24,1,0,5)='neutron'
      PARTIC(24,1,0,6)='proton'
      JDEC(24)=0
      NCH(24)=0
      NPLEV(24,1)=4

C ***   Canale 25
      PARTIC(25,0,0,1)='neutron'
      PARTIC(25,1,0,1)='pi+'
      PARTIC(25,1,0,2)='pi-'
      PARTIC(25,1,0,3)='pi0'
      PARTIC(25,1,0,4)='neutron'
      PARTIC(25,1,0,5)='proton'
      PARTIC(25,1,0,6)='proton'
      JDEC(25)=0
      NCH(25)=0
      NPLEV(25,1)=4

C ***   Channel #26 includes n 3pi0, p pi- 2pi0, p pi+ 2pi-,
C ***   p pi+ 2pi- presented here as a placeholder only
      PARTIC(26,0,0,1)='neutron'
      PARTIC(26,1,0,1)='pi+'
      PARTIC(26,1,0,2)='pi-'
      PARTIC(26,1,0,3)='pi-'
      PARTIC(26,1,0,4)='proton'
      PARTIC(26,1,0,5)='proton'
      PARTIC(26,1,0,6)='proton'
      JDEC(26)=0
      NCH(26)=0
      NPLEV(26,1)=4

C ***   Channel #27 includes all the states (4-8)pi N in gamma p,
C ***   p 2pi+ 2pi- presented here as a placeholder only
      PARTIC(27,0,0,1)='proton' 
      PARTIC(27,1,0,1)='pi+'
      PARTIC(27,1,0,2)='pi+'
      PARTIC(27,1,0,3)='pi-'
      PARTIC(27,1,0,4)='pi-'
      PARTIC(27,1,0,5)='proton'
      PARTIC(27,1,0,6)='neutron'
      PARTIC(27,1,0,7)='proton'
      JDEC(27)=0
      NCH(27)=0
      NPLEV(27,1)=5

C ***   Channel #28 includes all the states (4-8)pi N in gamma n,
C ***   n 2pi+ 2pi- presented here as a placeholder only
      PARTIC(28,0,0,1)='neutron'
      PARTIC(28,1,0,1)='pi+'
      PARTIC(28,1,0,2)='pi+'
      PARTIC(28,1,0,3)='pi-'
      PARTIC(28,1,0,4)='pi-'
      PARTIC(28,1,0,5)='neutron'
      PARTIC(28,1,0,6)='proton'
      PARTIC(28,1,0,7)='proton'
      JDEC(28)=0
      NCH(28)=0
      NPLEV(28,1)=5

C ***   Canale 29
      PARTIC(29,0,0,1)='proton'
      PARTIC(29,1,0,1)='kappa+'
      PARTIC(29,1,0,2)='lambda'
      PARTIC(29,1,0,3)='neutron'
      PARTIC(29,1,0,4)='proton'
      JDEC(29)=0
      NCH(29)=0
      NPLEV(29,1)=2

C ***   Canale 30
      PARTIC(30,0,0,1)='proton'
      PARTIC(30,1,0,1)='kappa+'
      PARTIC(30,1,0,2)='sigma0'
      PARTIC(30,1,0,3)='neutron'
      PARTIC(30,1,0,4)='proton'
      JDEC(30)=0
      NCH(30)=0
      NPLEV(30,1)=2

C ***   Canale 31
      PARTIC(31,0,0,1)='proton'
      PARTIC(31,1,0,1)='pi+'
      PARTIC(31,1,0,2)='kappas'      !ATTENZIONE E' SHORT !!!!!!!
      PARTIC(31,1,0,3)='lambda'
      PARTIC(31,1,0,4)='neutron'
      PARTIC(31,1,0,5)='proton'
      JDEC(31)=0
      NCH(31)=0
      NPLEV(31,1)=3

C ***      CANALE COMPTON

C ***   Canale 32
      PARTIC(32,0,0,1)='proton'
      PARTIC(32,1,0,1)='gamma'
      PARTIC(32,1,0,2)='proton'
      PARTIC(32,1,0,3)='neutron'
      PARTIC(32,1,0,4)='proton'
      JDEC(32)=0
      NCH(32)=0
      NPLEV(32,1)=2

C ***   CANALI SPECIALI SU DEUTONE

C ***   Canale 33
      PARTIC(33,0,0,1)='deuteron'
      PARTIC(33,1,0,1)='pi0'
      PARTIC(33,1,0,2)='deuteron'
      PARTIC(33,1,0,3)='proton'
      JDEC(33)=0
      NCH(33)=0
      NPLEV(33,1)=2

C ***     Canale 34
      PARTIC(34,0,0,1)='deuteron'
      PARTIC(34,1,0,1)='proton'
      PARTIC(34,1,0,2)='neutron'
      PARTIC(34,1,0,3)='proton'
      JDEC(34)=0
      NCH(34)=0
      NPLEV(34,1)=2

C ***     Canale 35
C *** aggiunto per prova
      PARTIC(35,0,0,1)='deuteron'
      PARTIC(35,1,0,1)='rho0'
      PARTIC(35,1,0,2)='deuteron'
      PARTIC(35,1,0,3)='proton'
      JDEC(35)=1
      NCH(35)=1
      BR(35,1)=1.
      NPLEV(35,1)=2
      PARTIC(35,2,1,1)='pi+'
      PARTIC(35,2,1,2)='pi-'
      NPLEV(35,2)=2

C ***    Canale 36
C *** Aggiunto per prova
      PARTIC(36,0,0,1)='deuteron'
      PARTIC(36,1,0,1)='eta'
      PARTIC(36,1,0,2)='deuteron'
      PARTIC(36,1,0,3)='proton'
      JDEC(36)=0
      NCH(36)=0
      NPLEV(36,1)=2


C ***   CANALI SPECIALI SU HE3

C ***   Canale 37
        PARTIC(37,0,0,1)='delta++'
        PARTIC(37,1,0,1)='pi+'
        PARTIC(37,1,0,2)='proton'
        PARTIC(37,1,0,3)='neutron'
        PARTIC(37,1,0,4)='neutron'
        JDEC(37)=0
        NCH(37)=0
        NPLEV(37,1)=2

C ***   Canale 38
        PARTIC(38,0,0,1)='deuteron'
        PARTIC(38,1,0,1)='pi0'
        PARTIC(38,1,0,2)='deuteron'
        PARTIC(38,1,0,3)='proton'
        JDEC(38)=0
        NCH(38)=0
        NPLEV(38,1)=2

C ***   Canale 39
        PARTIC(39,0,0,1)='deuteron'
        PARTIC(39,1,0,1)='pi+'
        PARTIC(39,1,0,2)='pi-'
        PARTIC(39,1,0,3)='deuteron'
        PARTIC(39,1,0,4)='proton'
        JDEC(39)=0
        NCH(39)=0
        NPLEV(39,1)=3

C ***   Canale 40
        PARTIC(40,0,0,1)='deuteron'
        PARTIC(40,1,0,1)='pi-'
        PARTIC(40,1,0,2)='dibaryon'
        PARTIC(40,1,0,3)='proton'
        JDEC(40)=2
        NCH(40)=1
        BR(40,1)=1.
        PARTIC(40,2,1,1)='proton'
        PARTIC(40,2,1,2)='proton'
        NPLEV(40,1)=2
        NPLEV(40,2)=2

      RETURN
      END                                            !*** END CHANNELDEF


C***********************************************************************

      SUBROUTINE TOTALSECT

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C FUNCTIONAL DESCRIPTION:
C
C    Sum of total cross sections of reactions tabulated in BLOCK DATA and 
C    ordered for Monte Carlo simulation.
C
C MODIFICATION DATE:
C 
C     19.05.95      by Igor Pshenichnov
C
C   **************   by marco mirazita
C
C     23.04.99
C --> Aggiunto l'azzeramento del vettore NPRO, altrimenti i valori non 
C     definiti dopo sono numeri casuali.
C
C     05.05.99
C --> Modificata la definizione dei valori di NPRO per i nuovi canali
C     su d.
C-----------------------------------------------------------------------

#include "process.inc"

      INTEGER IBEAM,IELET,NCHAIN,ICHAIN(40),ICIBLE,IPOL,ILAM
      REAL ETAG,FWHM
      COMMON/FLAG_COMP/IELET,IBEAM,ICIBLE,ETAG,FWHM,NCHAIN,ICHAIN,IPOL,
     &   ILAM
      REAL WR(263),SIGR(40,263)
      COMMON/SIGTO/WR,SIGR
      INTEGER ITARG,IIDH
      CHARACTER*12 TARGET,TARGFERMI
      COMMON /TARGET/ ITARG,IIDH,TARGET,TARGFERMI
      REAL HMAX,SIGTOTRIF(263),SIGSUM(263)
      COMMON/CR_TOT/SIGSUM
      INTEGER NPRO(40)

      CHARACTER*50 NCANAL(40)
      COMMON/CLABEL/NCANAL

C *********************************
C *** azzeramento
	do j=1,40
           npro(j)=0
	enddo
C *********************************

      DO J = 1,263
         SIGSUM(J) = 0.
      ENDDO

      NSUM = NCHAIN

      IF (TARGET.EQ.'neutron     ') THEN
         DO J=1,14
            NPRO(J) = NPRON(J)
         END DO
         IF (NCHAIN.EQ.0) NSUM = 20
      ELSE IF (TARGET.EQ.'proton      ') THEN
         DO J=1,18
            NPRO(J) = NPROP(J)
         END DO
         IF (NCHAIN.EQ.0) NSUM = 20
      ELSE IF (TARGET.EQ.'deuteron    ') THEN
         DO J=1,18
            NPRO(J)    = NPROP(J)
            NPRO(J+18) = NPRON(J)
         END DO

         NPRO(33) = NPROD(1)
         NPRO(34) = NPROD(2)
         NPRO(35) = NPROD(3)
         NPRO(36) = NPROD(4)

C ***********************************************
         IF (NCHAIN.EQ.0) NSUM = 40
      ELSE IF (TARGET.EQ.'He3         ') THEN
         DO J=1,18
            NPRO(J)    = NPROP(J)
            NPRO(J+18) = NPRON(J)
         END DO
         NPRO(37) = NPROHE3(1)
         NPRO(38) = NPROHE3(2)
         NPRO(39) = NPROHE3(3)
         NPRO(40) = NPROHE3(4)
         IF (NCHAIN.EQ.0) NSUM = 40
      END IF

      IF (NCHAIN.EQ.0) THEN
        DO J=1,40
          ICHAIN(J) = NPRO(J)
        END DO
      END IF


      WRITE(15,*) 'Reaction channels:'
      DO M=1,NSUM

        DO N=1,40
          IF(ICHAIN(M).EQ.NPRO(N).AND.NPRO(N).NE.0) GO TO 1
        END DO

        IF (N.EQ.41) THEN
           WRITE(15,*) 'Channel',ICHAIN(M),' not defined for target ',TARGET
           GOTO 2
        ENDIF

1       ICHANNEL = NPRO(N)

        IF (ICHANNEL.GE.33) IIDH = 1

        IF(TARGET.EQ.'He3         ') THEN
            IF((ICHANNEL.EQ.1.OR.ICHANNEL.EQ.2.OR.ICHANNEL.EQ.5.OR.
     &         ICHANNEL.EQ.6.OR.ICHANNEL.EQ.11.OR.ICHANNEL.EQ.12.OR.
     &         ICHANNEL.EQ.15.OR.ICHANNEL.EQ.19.OR.ICHANNEL.EQ.21.OR.
     &         ICHANNEL.EQ.23.OR.ICHANNEL.EQ.24.OR.ICHANNEL.EQ.27.OR.
     &        (ICHANNEL.GE.29.AND.ICHANNEL.LE.32))) THEN
                  DO K=1,263
                     SIGTOTRIF(K) = SIGR(ICHANNEL,K)*2.
                  ENDDO
            ELSEIF (ICHANNEL.EQ.7) THEN
                  DO K=1,263
                     SIGTOTRIF(K) = SIGR(5,K)*.67
                  ENDDO
            ELSEIF (ICHANNEL.EQ.16.OR.ICHANNEL.EQ.18) THEN
                  DO K=1,263
                     SIGTOTRIF(K) = SIGR(15,K)
                  ENDDO
            ELSEIF (ICHANNEL.EQ.17) THEN
                  DO K=1,263
                     SIGTOTRIF(K) = SIGR(15,K)*2.
                  END DO
            ELSE
                  DO K=1,263
                     SIGTOTRIF(K) = SIGR(ICHANNEL,K)
                  ENDDO
            ENDIF
        ELSE  ! DEUTERON AND NUCLEONS 
          DO K=1,263
            SIGTOTRIF(K) = SIGR(ICHANNEL,K)
          ENDDO
        ENDIF


        DO K=1,263
            SIGSUM(K) = SIGSUM(K) + SIGTOTRIF(K)
        END DO
        IF (ICHANNEL.NE.0) WRITE(15,*) '  ',ICHANNEL,NCANAL(ICHANNEL)
2     END DO

	DO K=1,263
	IF (SIGSUM(K).GE.1.E-8) GOTO 3
	ENDDO
	STOP ' PROGRAM STOPPED, NO POSSIBLE REACTION'

3	CONTINUE
      WRITE(15,*) 
      RETURN
      END

