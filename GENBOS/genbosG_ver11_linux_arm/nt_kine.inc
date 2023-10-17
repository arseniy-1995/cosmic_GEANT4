C *** canale di reazione
	common/reacch/j_channel
	integer j_channel

C *** common per le variabili della cinematica di primo livello
        common/kine_1/I_1,npart_1,pprim,P_1
        real P_1(4,11),pprim(5,18)
        integer I_1(4),npart_1

C *** Impulso di fermi
	common/fermi_motion/p_fermi,rm
	real p_fermi,rm(2)

C *** Beta(CM) of (gamma+target) system
	common/betacm/beta_cm
	real beta_cm(4)
