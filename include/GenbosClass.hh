//
// Created by Арсений Юрченко on 29.10.2022.
//

#ifndef COSMIC_GENBOSCLASS_HH
#define COSMIC_GENBOSCLASS_HH

#include "Genbos.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"

//#define GENBOSCLASSON

namespace Cosmic {

#ifdef GENBOSCLASSON
    class GenbosClass {


    public:
        GenbosClass(G4int *RN);

        ~GenbosClass();

        void GenbosGenerator();

        void SetEgMin(G4double val_min) {
            G4float EgMin = val_min / GeV;
            G4int n = 2;
            //genbos_beam_(&n, &EgMin, &few);
            fbe = n;
            feg = EgMin;

            // genbos_beam_(&n, &EgMin, &EgMax);

            G4float EgMean = (EgMax + EgMin) / 2.;
            G4float EgWidht = (EgMax - EgMin);
            genbos_beam_(&n, &EgMean, &EgWidht);

        }

        void SetEgMax(G4double val_max) {

            G4float EgMax = val_max / GeV;
            G4int n = 2;
           // genbos_beam_(&n, &feg, &EgMax);
            fbe = n;
            few = EgMax;

            // genbos_beam_(&n, &EgMin, &EgMax);

            G4float EgMean = (EgMax + EgMin) / 2.;
            G4float EgWidht = (EgMax - EgMin);
            genbos_beam_(&n, &EgMean, &EgWidht);
        }

        void SetEgMinMax(G4double val_min, G4double val_max) {
            G4float EgMin = val_min / GeV;
            G4float EgMax = val_max / GeV;
            G4int n = 2;
           // genbos_beam_(&n, &EgMin, &EgMax);
            fbe = n;
            feg = EgMin;
            few = EgMax;

            // genbos_beam_(&n, &EgMin, &EgMax);

            G4float EgMean = (EgMax + EgMin) / 2.;
            G4float EgWidht = (EgMax - EgMin);
            genbos_beam_(&n, &EgMean, &EgWidht);
        }

        void SetRandom(G4int ix) {
            genbos_rand_(&ix);
            fix = ix;
        }

        void SetMode(G4int val_mode) {
            int nreac, ireac[40];

            //   if (val_mode == 34) { // это p+n
            memset(ireac, 0, sizeof(ireac));
            nreac = 1;
            ireac[0] = val_mode;
            // G4AutoLock lock(&aMutex);
            genbos_reactions_(&nreac, ireac);
            //  genbos_start_(&FileNum);

            //   }
        }

        G4float GetEfot() { return fefot; }

        G4int GetNreac() { return fnreac; }

        G4int GetNp() { return fnp; }

        G4int *GetId() { return fid; }

        G4float *GetPx() { return fpx; }

        G4float *GetPy() { return fpy; }

        G4float *GetPz() { return fpz; }

        G4int GetIx() { return fix; }

        G4int GetBe() { return fbe; }

        G4int GetRN() { return fRN; }


    private:

        // for genbos_reactions_(int* nreac, int *ireac)
        G4int fnreac; //	nreac -- number of considered reactions
        G4int *fireac; //	ireac[nreac] -- list of considered reactions

        // for genbos_targ_(int* tg)
        G4int ftg; //	tg -- target nucleus:  0-neutron 1-proton 2-deuteron

        // for genbos_beam_(int* be, float* eg, float* ew)
        G4int fbe; //	be -- photon beam spectrum: 0-gaussian, 2-bremsstrahlung, 3-uniform
        G4float feg; //	egmin -- Egamma min
        G4float few; //	egmax -- Egamma max

        // genbos_start_(int *RN)
        G4int fRN;

        // for genbos_stop_()


        // for genbos_event_(float* efot, int* nreac, int* np, int *id, float *px, float* py, float* pz)

        G4float fefot; //	efot - photon energy, GeV
        // int fnreac; //	nreac - reaction index
        G4int fnp; //	np - number of generated particles
        G4int *fid = new G4int[11]; //	id[np] - list of particle's indexes
        G4float *fpx = new G4float[11]; //	px[np],py[np],pz[np] - lists of momentum components
        G4float *fpy = new G4float[11];
        G4float *fpz = new G4float[11];

        // for genbos_rand_(int *ix)

        G4int fix; //	ix -- random seed


    };
#endif // GENBOSCLASSON

}

#endif //COSMIC_GENBOSCLASS_HH
