//
// Created by Арсений Юрченко on 19.10.2022.
//



#ifndef CosmicGenbos_h
#define CosmicGenbos_h 1

extern"C" {

void genbos_reactions_(int* nreac, int *ireac);
//	nreac -- number of considered reactions
//	ireac[nreac] -- list of considered reactions

void genbos_targ_(int* tg);
//	tg -- target nucleus:  0-neutron 1-proton 2-deuteron

void genbos_beam_(int* be, float* eg, float* ew);
//	be -- photon beam spectrum: 0-gaussian, 2-bremsstrahlung, 3-uniform
//	eg -- Egamma mean
//	ew -- Egamma width
//  -->  egmin = eg - ew/2  ;   egmax = eg + ew/2



void genbos_start_(int *RN);

void genbos_stop_();

void genbos_event_(float* efot, int* nreac, int* np, int *id, float *px, float* py, float* pz);
//	efot - photon energy, GeV
//	nreac - reaction index
//	np - number of generated particles
//	id[np] - list of particle's indexes
//	px[np],py[np],pz[np] - lists of momentum components

void genbos_rand_(int *ix);
//	ix -- random seed


}


#endif //CosmicGenbos_h
