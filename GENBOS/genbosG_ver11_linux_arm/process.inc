      INTEGER NPROP(20),NPRON(20),NPROD(4),NPROHE3(4)
      COMMON /FLAG_PROC/ NPROP,NPRON,NPROD,NPROHE3

	integer nprop_tot,npron_tot,nprod_tot,nprohe3_tot
	common/proc_tot/nprop_tot,npron_tot,nprod_tot,nprohe3_tot

	integer nchain_p,nchain_n,nchain_d,nchain_he3
	common/proc_sel/nchain_p,nchain_n,nchain_d,nchain_he3
