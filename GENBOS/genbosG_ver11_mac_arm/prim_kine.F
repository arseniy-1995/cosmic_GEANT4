C***********************************************************************

      SUBROUTINE prim_kine

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                                                           
C FUNCTIONAL DESCRIPTION:
C
C Extracts kinematics of outgoing particles after first level
C and before decays. Channels 5->10 (delta), 11->14 and 35 (rho)
C are involved.
C
C MODIFICATION DATE:
C 
C     05.05.99       by marco mirazita
C
C                                                          
C FORMAL PARAMETERS:                                              
C                                                            
C    Npart_1    = NUMBER OF OUTGOING PARTICLES (2-18)
C    RMASS      = MASS OF I-TH spectator PARTICLE (REAL)
C    PN(1,I)    = PX OF I-TH spectator PARTICLE
C    PN(2,I)    = PY OF I-TH spectator PARTICLE
C    PN(3,I)    = PZ OF I-TH spectator PARTICLE
C    PN(4,I)    = ENERGY OF I-TH spectator PARTICLE
C    PN(5,I)    = |P| OF I-TH spectator PARTICLE
C
C    PPRIM(1,I)   = PX OF I-TH primary partICLE
C    PPRIM(2,I)   = PY OF I-TH primary partICLE
C    PPRIM(3,I)   = PZ OF I-TH primary partICLE
C    PPRIM(4,I)   = ENERGY OF I-TH primary partICLE
C
C    P_1(1,I)  =  PX OF I-TH primary partICLE
C    P_1(2,I)  =  PY OF I-TH primary partICLE
C    P_1(3,I)  =  PZ OF I-TH primary partICLE
C    P_1(4,I)  =  MASS OF I-TH primary partICLE
C
C COMMON BLOCKS:
C     /TARGET/                                                    
C     /CHANNELS/                                              
C     /MOM_NUCLEON/
C-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER ITARG,IIDH
      CHARACTER*12 TARGET,TARGFERMI
      COMMON /TARGET/ ITARG,IIDH,TARGET,TARGFERMI

      CHARACTER*12 PARTIC(40,0:3,0:3,8)
      INTEGER JDEC(40),NCH(40),NPLEV(40,3)
      REAL BR(40,3)
      COMMON /CHANNELS/ PARTIC,JDEC,NCH,NPLEV,BR

      COMMON /MOM_NUCLEON/PN
      real pn(5,3)

      real rmass,xm
      integer j,jj,jch,code

C *** ntuple variables
C     #include "nt_kine.inc"
      INCLUDE 'nt_kine.inc'
      jch=j_channel


C *** initialization
      npart_1=0
      do j=1,4
         do jj=1,4
            P_1(jj,j)=0.
         enddo
         I_1(j)=0
      enddo


C ***    Canali 5,6,7,8,9,10
      IF (JCH.GE.5.AND.JCH.LE.10) THEN
         p_1(1,1) = pprim(1,1)
         p_1(2,1) = pprim(2,1)
         p_1(3,1) = pprim(3,1)
c         p_1(4,1) = pprim(4,1)
         p_1(4,1) = sqrt( pprim(4,1)**2.- (pprim(1,1)**2.+
     +        pprim(2,1)**2.+pprim(3,1)**2.) )
c         write (*,'(a)') PARTIC(JCH,1,0,1)
         I_1(1) = CODE(PARTIC(JCH,1,0,1))
         
         p_1(1,2) = pprim(1,2)
         p_1(2,2) = pprim(2,2)
         p_1(3,2) = pprim(3,2)
c         p_1(4,2) = pprim(4,2)
         p_1(4,2) = sqrt( pprim(4,2)**2.- (pprim(1,2)**2.+
     +        pprim(2,2)**2.+pprim(3,2)**2.) )
c        write (*,'(a)') PARTIC(JCH,1,0,2)
         I_1(2) = CODE(PARTIC(JCH,1,0,2))

         NPART_1 = 2
         IF (TARGET.EQ.'deuteron'.OR.TARGET.EQ.'He3') THEN
            p_1(1,3) = PN(1,2)
            p_1(2,3) = PN(2,2)
            p_1(3,3) = PN(3,2)
            p_1(4,3) = RMASS(PARTIC(JCH,1,0,3))
            I_1(3) = CODE(PARTIC(JCH,1,0,3))
            NPART_1 = NPART_1+1
            IF (TARGET.EQ.'He3') THEN
               p_1(1,4) = PN(1,3)
               p_1(2,4) = PN(2,3)
               p_1(3,4) = PN(3,3)
               p_1(4,4) = RMASS(PARTIC(JCH,1,0,4))
               I_1(4) = CODE(PARTIC(JCH,1,0,4))
               NPART_1 = NPART_1+1
            END IF
         END IF
         GOTO 99
      END IF

C ***    Canali 11,12,13,14
      IF (JCH.GE.11.AND.JCH.LE.14) THEN
         p_1(1,1) = pprim(1,1)
         p_1(2,1) = pprim(2,1)
         p_1(3,1) = pprim(3,1)
c         p_1(4,1) = pprim(4,1)
         p_1(4,1) = sqrt( pprim(4,1)**2.- (pprim(1,1)**2.+
     +        pprim(2,1)**2.+pprim(3,1)**2.) )
c        write (*,'(a)') PARTIC(JCH,1,0,1)
         I_1(1) = CODE(PARTIC(JCH,1,0,1))
         
         p_1(1,2) = pprim(1,2)
         p_1(2,2) = pprim(2,2)
         p_1(3,2) = pprim(3,2)
c         p_1(4,2) = pprim(4,2)
         p_1(4,2) = sqrt( pprim(4,2)**2.- (pprim(1,2)**2.+
     +        pprim(2,2)**2.+pprim(3,2)**2.) )
c        write (*,'(a)') PARTIC(JCH,1,0,2)
         I_1(2) = CODE(PARTIC(JCH,1,0,2))

         NPART_1 = 2
         IF (TARGET.EQ.'deuteron'.OR.TARGET.EQ.'He3') THEN
            p_1(1,3) = PN(1,2)
            p_1(2,3) = PN(2,2)
            p_1(3,3) = PN(3,2)
            p_1(4,3) = RMASS(PARTIC(JCH,1,0,3))
            I_1(3) = CODE(PARTIC(JCH,1,0,3))
            NPART_1 = NPART_1+1
            IF (TARGET.EQ.'He3') THEN
               p_1(1,4) = PN(1,3)
               p_1(2,4) = PN(2,3)
               p_1(3,4) = PN(3,3)
               p_1(4,4) = RMASS(PARTIC(JCH,1,0,4))
               I_1(4) = CODE(PARTIC(JCH,1,0,4))
               NPART_1 = NPART_1+1
            END IF
         END IF
         GOTO 99
      END IF


C ***    Canale 35
      IF (JCH.EQ.35) THEN
         p_1(1,1) = pprim(1,1)
         p_1(2,1) = pprim(2,1)
         p_1(3,1) = pprim(3,1)
c         p_1(4,1) = pprim(4,1)
         p_1(4,1) = sqrt( pprim(4,1)**2.- (pprim(1,1)**2.+
     +        pprim(2,1)**2.+pprim(3,1)**2.) )
c     write (*,'(a)') PARTIC(JCH,1,0,1)
         I_1(1) = CODE(PARTIC(JCH,1,0,1))
         
         p_1(1,2) = pprim(1,2)
         p_1(2,2) = pprim(2,2)
         p_1(3,2) = pprim(3,2)
c         p_1(4,2) = pprim(4,2)
         p_1(4,2) = sqrt( pprim(4,2)**2.- (pprim(1,2)**2.+
     +        pprim(2,2)**2.+pprim(3,2)**2.) )
c     write (*,'(a)') PARTIC(JCH,1,0,2)
         I_1(2) = CODE(PARTIC(JCH,1,0,2))

         NPART_1 = 2
         IF (TARGET.EQ.'He3') THEN
            p_1(1,4) = PN(1,3)
            p_1(2,4) = PN(2,3)
            p_1(3,4) = PN(3,3)
            p_1(4,4) = RMASS(PARTIC(JCH,1,0,4))
            I_1(4) = CODE(PARTIC(JCH,1,0,4))
            NPART_1 = NPART_1+1
         END IF
         GOTO 99
      END IF


 99   continue

      RETURN
      END                       !*** END prim_kine
