


	REAL*4 FUNCTION TOP_GAMN(KIN_ENERGY,TOP_NUMBER)
C+
C 
C FUNCTIONAL DESCRIPTION:	
C 
C    The function presents a calculation of gamma-
C    nucleon TOPOLOGICAL cross sections by means of approximation
C    of cross sections of certain channel and SU(2) statistical model
C
C CREATION DATE: 
C 
C     18.03.94       by Igor Pshenichnov 
C 
C FORMAL PARAMETERS:
C 
C    KIN_ENERGY  - projectile kinetic energy (GeV) 
C 
C    TOP_NUMBER  - topology (number of particles in the final
C                  state (n=3,...,9)).
C 
C COMMON BLOCKS:
C
C    none 
C 
C FUNCTION VALUE:
C 
C    cross-section value (mb)
C 
C-
	REAL*4 KIN_ENERGY
        INTEGER*4 TOP_NUMBER
        
                IF(TOP_NUMBER.EQ.3) THEN
                      TOP_GAMN=GAMN(KIN_ENERGY,6)/0.7
                ELSEIF(TOP_NUMBER.EQ.4) THEN
                      TOP_GAMN=GAMN(KIN_ENERGY,9)/0.75
                ELSEIF(TOP_NUMBER.EQ.5) THEN
                      TOP_GAMN=GAMN(KIN_ENERGY,14)/0.2667
                ELSEIF(TOP_NUMBER.EQ.6) THEN
                      TOP_GAMN=GAMN(KIN_ENERGY,19)/0.4381
                ELSEIF(TOP_NUMBER.EQ.7) THEN
                      TOP_GAMN=GAMN(KIN_ENERGY,26)/0.125
                ELSEIF(TOP_NUMBER.EQ.8) THEN
                      TOP_GAMN=GAMN(KIN_ENERGY,33)/0.2755
                ELSEIF(TOP_NUMBER.EQ.9) THEN
                      TOP_GAMN=GAMN(KIN_ENERGY,42)/0.05614
                ELSE
                   WRITE(15,1) TOP_NUMBER
 1                 FORMAT(' A WRONG TOPOLOGY NUMBER: ',I4)
                   TOP_GAMN=0.0
                   RETURN
                ENDIF
 
        RETURN
	END



	REAL*4 FUNCTION GAMN(KIN_ENERGY,TAB_NUMBER)
C+
C 
C FUNCTIONAL DESCRIPTION:	
C 
C    The function presents a calculation of approximation gamma-
C    nucleon partial cross-sections from experiments:
C    S.I.Alekhin et al., CERN-HERA 87-01
C
C CREATION DATE: 
C 
C     16.03.94       by Igor Pshenichnov 
C 
C FORMAL PARAMETERS:
C 
C    KIN_ENERGY  - projectile kinetic energy (GeV) 
C 
C    TAB_NUMBER  - reaction number (see table)
C 
C COMMON BLOCKS:
C 
C    /LAGERR/PL(0:20) - Lagerr polinom values
C
C    /RE_NAME/RE_NAME - reaction notation in output
C 
C FUNCTION VALUE:
C 
C    cross-section value (mb)
C 
C-
	REAL*4 PL(0:20),KIN_ENERGY,TK(7),A(0:4,7)
        REAL*4 ALFA,C_UNI,C_DPI
        INTEGER*4 TAB_NUMBER,RE_NUMBER,MAX_POWER(7)
        CHARACTER*45 RE_NAME,ST_NAMES(7)

        COMMON /LAGERR/PL
	COMMON /RE_NAME/RE_NAME
C ...............................................
	DATA ST_NAMES/' 6: g p ==> p pi+ pi-',
     #                ' 9: g p ==> p pi+ pi- pi0',
     #                '14: g p ==> p pi+ pi+ pi- pi-',
     #                '19: g p ==> p pi+ pi+ pi- pi- pi0',
     #                '26: g p ==> p pi+ pi+ pi+ pi- pi- pi-',
     #                '33: g p ==> p pi+ pi+ pi+ pi- pi- pi- pi0',     
     #                '42: g p ==> p pi+ pi+ pi+ pi+ pi- pi- pi- pi-'/        
     
		DATA TK/  321.00,   506.00,  727.00,   952.00,
     #                   1215.00,  1481.00, 1788.00/

C               .....................   6  ...........................
	DATA A/ 0.3317900,-0.0821800, 0.0997600, 0.0001300, 0.0000000,
C               .....................   9  ...........................
     #          0.0537900, 0.1350400,-0.0317300, 0.0000000, 0.0000000,
C               .....................  14  ...........................
     #         -0.0009800, 0.0983250,-0.0475448, 0.0000000, 0.0000000,
C               .....................  19  ...........................
     #          0.2056800,-0.0629900, 0.0000000, 0.0000000, 0.0000000,
C               .....................  26  ...........................
     #          0.0619000,-0.0192100, 0.0000000, 0.0000000, 0.0000000,
C               .....................  33  ...........................
     #          0.1113700,-0.0409400, 0.0000000, 0.0000000, 0.0000000,
C               .....................  42  ...........................
     #          0.0336780,-0.0130250, 0.0000000, 0.0000000, 0.0000000/

		DATA MAX_POWER/  3,  2,  2,  1,  1,  1,  1/    

                      ALFA=2.0
                      C=1.0

                IF(TAB_NUMBER.EQ.6) THEN
                      RE_NUMBER=1
                ELSEIF(TAB_NUMBER.EQ.9) THEN
                      RE_NUMBER=2
                ELSEIF(TAB_NUMBER.EQ.14) THEN
                      RE_NUMBER=3
                ELSEIF(TAB_NUMBER.EQ.19) THEN
                      RE_NUMBER=4
                ELSEIF(TAB_NUMBER.EQ.26) THEN
                      RE_NUMBER=5
                ELSEIF(TAB_NUMBER.EQ.33) THEN
                      RE_NUMBER=6
                ELSEIF(TAB_NUMBER.EQ.42) THEN
                      RE_NUMBER=7
                ELSE
                   WRITE(15,100) TAB_NUMBER
 100               FORMAT('   THIS IS NOT A BASE REACTION: ',I4)
                   GAMN=0.0
                   RETURN
                ENDIF
 
                RE_NAME=ST_NAMES(RE_NUMBER)
                IF(KIN_ENERGY*1000.0.LE.TK(RE_NUMBER)) THEN
                     GAMN=0.0
                     RETURN
                ENDIF

        X=C*ALOG(KIN_ENERGY*1000./TK(RE_NUMBER))   
	ALV=AL(MAX_POWER(RE_NUMBER),ALFA,X)
                      F=0.0
        DO 1 N=0,MAX_POWER(RE_NUMBER)
        F=F+A(N,RE_NUMBER)*PL(N)
 1      CONTINUE
        F=F*(X**(ALFA/2.0))/EXP(X/2.0)
        GAMN=F**2 
	RETURN 
	END


	REAL*4 FUNCTION AL(N,ALFA,X)
C+
C 
C FUNCTIONAL DESCRIPTION:	
C
C The recurrent calculation of LAGERR's polinoms up to order=20 
C (See for details: Handbook of mathematical functions. Ed. by 
C M.Abramowitz and I.A.Stegun)
C
C CREATION DATE: 
C 
C     08.02.93       by Igor Pshenichnov 
C 
C FORMAL PARAMETERS:
C 
C    N    -  polinom's order
C    ALFA -  parameter
C    X    -  point to be calculated 
C 
C COMMON BLOCKS:
C 
C    /LAGERR/PLAGERR(0:20) array contais polinom's values from 
C                      zero oder to N 
C 
C FUNCTION VALUE:
C 
C    AL   -  N-oder value of LAGERR's polinom    
C 
C-
	COMMON /LAGERR/PLAGERR(0:20)

        DO INDEX=0,20
          PLAGERR(INDEX)=0.0
        END DO
	 
	PLAGERR(0) = 1.
        AL = PLAGERR(0)
	IF (N.EQ.0) RETURN
C
	PLAGERR(1) = ALFA+1-X
        AL = PLAGERR(1)
        IF (N.EQ.1) RETURN
C       
        DO 1 M=1,N-1 
               PLAGERR(M+1) = ((2*M+ALFA+1-X)*PLAGERR(M)-
     #                             (M+ALFA)*PLAGERR(M-1))/(M+1)
  1     CONTINUE
        AL = PLAGERR(N)
        RETURN
	END



        SUBROUTINE GAMN_CHARG(TOP_NUMBER,CHARGE,PART_NUMB)
C+
C 
C FUNCTIONAL DESCRIPTION:	
C 
C    This subroutine is a generator of charges of outgoing particles.
C    It contains SU(2) coefficients table for N+mpi
C 
C CREATION DATE: 
C 
C     14.03.94       by Igor Pshenichnov
C
C FORMAL PARAMETERS:
C 
C 
C    TOP_NUMBER  - topology (number of particles in the final
C                  state (m=3,...,9)).
C 
C    CHARGE      - charge of the taget nucleon   
C 
C    PART_NUMB   - outgoing particle numbers, nucleon comes first
C    
C 
C COMMON BLOCKS:
C
C    /LOGEG/EG,CMSF   - governs the mode of operation: 
C
C   EG:
C   PARTICLE NUMBERS ARE AS IN INC (.FALSE.) OR AS IN GEANT (.TRUE.)
C   
C   CMSF:
C   All the momenta of particles are in CMS (.TRUE.)
C   or in LAB system (.FALSE.)
C    
C-
        
        INTEGER*4 TOP_NUMBER, PART_NUMB(10), NI, NF, NUMPI(3)
        INTEGER*4 NUMNUC(2),CHARGE,IRENUM(7)
        INTEGER*2 CHARG_DAT(4,42)
        REAL*4 SU2(42), PROBAB(9)
        LOGICAL EG,CMSF
        COMMON/LOGEG/EG,CMSF ! Insure compartibility with the INC model.

        real vran(1)


	DATA CHARG_DAT/
C --------- notation is : N(p=1, n=0) kpi+ lpi- mpi0 --------
C ________________________ 3 PARTICLES ______________________
     # 1, 0, 0, 2,   1, 1, 1, 0,   0, 1, 0, 1,        
C ________________________ 4 PARTICLES ______________________
     # 1, 0, 0, 3,   1, 1, 1, 1,   0, 1, 0, 2,   0, 2, 1, 0,
C ________________________ 5 PARTICLES ______________________
     # 1, 0, 0, 4,   1, 1, 1, 2,   1, 2, 2, 0,   0, 1, 0, 3,
     # 0, 2, 1, 1,
C ________________________ 6 PARTICLES ______________________
     # 1, 0, 0, 5,   1, 1, 1, 3,   1, 2, 2, 1,   0, 1, 0, 4,
     # 0, 2, 1, 2,   0, 3, 2, 0, 
C ________________________ 7 PARTICLES ______________________
     # 1, 0, 0, 6,   1, 1, 1, 4,   1, 2, 2, 2,   1, 3, 3, 0,
     # 0, 1, 0, 5,   0, 2, 1, 3,   0, 3, 2, 1,
C ________________________ 8 PARTICLES ______________________
     # 1, 0, 0, 7,   1, 1, 1, 5,   1, 2, 2, 3,   1, 3, 3, 1,
     # 0, 1, 0, 6,   0, 2, 1, 4,   0, 3, 2, 2,   0, 4, 3, 0,
C ________________________ 9 PARTICLES ______________________                  
     # 1, 0, 0, 8,   1, 1, 1, 6,   1, 2, 2, 4,   1, 3, 3, 2,
     # 1, 4, 4, 0,   0, 1, 0, 7,   0, 2, 1, 5,   0, 3, 2, 3,
     # 0, 4, 3, 1/

	DATA SU2/
C ________________________ 3 PARTICLES ______________________
     # 0.10000000,   0.70000000,   0.20000000,
C ________________________ 4 PARTICLES ______________________
     # 0.01000000,   0.75000000,   0.04000000,   0.20000000,
C ________________________ 5 PARTICLES ______________________
     # 0.03333000,   0.36670000,   0.26670000,   0.06667000,
     # 0.26670000,
C ________________________ 6 PARTICLES ______________________
     # 0.00476200,   0.22380000,   0.43810000,   0.02857000,
     # 0.20950000,   0.09524000,
C ________________________ 7 PARTICLES ______________________
     # 0.00476000,   0.11550000,   0.42140000,   0.12500000,
     # 0.01190000,   0.14290000,   0.17860000,
C ________________________ 8 PARTICLES ______________________
     # 0.00061000,   0.05632000,   0.33420000,   0.27550000,
     # 0.00488400,   0.08425000,   0.20150000,   0.04273000,
C ________________________ 9 PARTICLES ______________________
     # 0.00061050,   0.02624000,   0.22600000,   0.35770000,
     # 0.05614000,   0.00383100,   0.09195000,   0.35630000,
     # 0.21460000/


        IF (EG) THEN
          NUMNUC(1)=14  ! proton number
          NUMNUC(2)=13  ! neutron number
          NUMPI(1)=8
          NUMPI(2)=9
          NUMPI(3)=7
          SU2(4)=0.090902
          SU2(5)=0.0  ! since it is simulating via different branch
          SU2(6)=0.363552
          SU2(7)=0.545464
        ELSE   ! "as was"
          NUMNUC(1)=37  ! proton number
          NUMNUC(2)=38  ! neutron number
          NUMPI(1)=1
          NUMPI(2)=2
          NUMPI(3)=7
          SU2(4)=0.01
          SU2(5)=0.75
          SU2(6)=0.04
          SU2(7)=0.2     
        ENDIF

        DO I=1,10
          PART_NUMB(I)=0.0
        ENDDO

C ____ DETERMINE POSITIONS OF CERTAIN MULTIPLICITY IN THE TABLE ___

        NF=TOP_NUMBER*(TOP_NUMBER+1)/2+1
        NI=NF-TOP_NUMBER+1
        SUM=0.
            DO I=1,TOP_NUMBER
                  SUM=SUM+SU2(NI+I-5)
                  PROBAB(I)=SUM
            ENDDO
        call ranmar(vran,1)
        TEMP1=vran(1)
        IABSOL=IDIHOTOMIA(TEMP1,PROBAB,TOP_NUMBER)+NI-1
 
        IF(CHARG_DAT(1,IABSOL-4).EQ.1) THEN
            PART_NUMB(1)=NUMNUC(1)
        ELSEIF(CHARG_DAT(1,IABSOL-4).EQ.0) THEN
            PART_NUMB(1)=NUMNUC(2)
        ENDIF

        IPOSITION=1

        DO K=1,3

        IF (CHARG_DAT(K+1,IABSOL-4).GT.0) THEN
           DO I=1,CHARG_DAT(K+1,IABSOL-4)
               IPOSITION=IPOSITION+1
               PART_NUMB(IPOSITION)=NUMPI(K)
           ENDDO
        ENDIF

        ENDDO
   
                  
      IF(CHARGE.EQ.0) THEN   ! CASE OF GAMMA + NEUTRON
      I=1
      DO WHILE (PART_NUMB(I).GT.0)
        IF(PART_NUMB(I).EQ.NUMNUC(1)) THEN
              PART_NUMB(I)=NUMNUC(2)
        ELSEIF(PART_NUMB(I).EQ.NUMNUC(2)) THEN
              PART_NUMB(I)=NUMNUC(1)
        ELSEIF(PART_NUMB(I).EQ.NUMPI(2)) THEN
              PART_NUMB(I)=NUMPI(1)
        ELSEIF(PART_NUMB(I).EQ.NUMPI(1)) THEN
              PART_NUMB(I)=NUMPI(2)
        ENDIF
         I=I+1
      ENDDO
      ENDIF
        
        RETURN
	END





