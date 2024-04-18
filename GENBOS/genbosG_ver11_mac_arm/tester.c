#include <stdio.h>
#include <stdlib.h>

#include "genbos.h"

int main(int argc, char **argv){

  int nev,nreac,np,id[11];
  float efot, px[11],py[11],pz[11];
  int i,k;
  int target=1;


	int ix=76762;
	
	genbos_targ_(&target);
	genbos_start_();
	genbos_rand_(&ix);
	
	for(k=0;k<500;k++){
		genbos_event_(&efot, &nreac, &np, id, px, py, pz);
		printf("%5.3f %2d %2d |",efot,nreac,np);
		for(i=0;i<np;i++){
			printf(" %2d % 5.3f % 5.3f % 5.3f |",id[i],px[i],py[i],pz[i]);
		}
		printf("\n");		
	}
	
	
	genbos_stop_();
	
	
	return 0;
	

}

