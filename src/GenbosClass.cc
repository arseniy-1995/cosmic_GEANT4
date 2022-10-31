//
// Created by Арсений Юрченко on 29.10.2022.
//

#include "GenbosClass.hh"
//#include "Genbos.hh"

namespace Cosmic {

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    GenbosClass::GenbosClass(int *RN) {

        genbos_start_(RN);
        fRN = *RN;
        // memset(fid, 0, sizeof(fid));
        // memset(fpx, 0, sizeof(fpx));
        // memset(fpy, 0, sizeof(fpy));
        // memset(fpz, 0, sizeof(fpz));
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    GenbosClass::~GenbosClass() {

        delete fid;
        delete fpx;
        delete fpy;
        delete fpz;
        genbos_stop_();

    }

    void GenbosClass::GenbosGenerator() {

        genbos_event_(&fefot, &fnreac, &fnp, fid, fpx, fpy, fpz);
    }

}
