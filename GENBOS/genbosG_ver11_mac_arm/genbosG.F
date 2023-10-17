	PROGRAM GENBOS
	implicit none

C ******************************************************************
C
C ====  RIA 2020 -- GEANT4
C

	INTEGER NSTATIS
	COMMON/NSTATIS/NSTATIS


C *** General of the event
	common/evnt/nev,efot,nreac
	integer nev,nreac
	real efot

C *** final state particles
	common/finst/np,idg(11),mass(11),imp(11)
     +		    ,px(11),py(11),pz(11)
	integer np,idg
	real mass,imp,px,py,pz

c *** local variables
	integer ido,iii,k
	real xx,rr,pmom

	INTEGER cut1
	REAL phi,the,deg

	CALL GENBOS_START
	CALL GENBOS_RAND(76762)

	DO  IDO=1,NSTATIS
C *** event generation
	   CALL GENBOS_EVENT(EFOT,NREAC,NP,IDG,PX,PY,PZ)
	   
	   WRITE(0,10) EFOT,NREAC,NP
	   
	ENDDO

10	FORMAT ('EFOT=',F5.3,' NREAC=',I3,' NP=',I3)

	CALL GENBOS_STOP
		
	STOP 
	END

