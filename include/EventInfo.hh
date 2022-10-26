//
// Created by Арсений Юрченко on 25.10.2022.
//

#ifndef COSMIC_EVENTINFO_HH
#define COSMIC_EVENTINFO_HH

#include "G4VUserEventInformation.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class EventInfo : public G4VUserEventInformation {
public:
    EventInfo();

    ~EventInfo();

    void Print() const {};

    // Set Methods
    void SetEgamma(G4double e) { fEgamma = e; };

    void SetXbeam(G4double val) { fXbeam = val; }

    void SetYbeam(G4double val) { fYbeam = val; }

    void SetPzz(G4double val) { fPzz = val; }

    void SetProCM(G4double val) { ProCM = val; }

    void SetNreac(G4int val) { nreac = val; }

    void SetNp(G4int val) { np = val; }

    void SetEntry(G4int val) { entry = val; }

    // Get Methods

    G4double GetEgamma() { return fEgamma; }

    G4double GetProCM() { return ProCM; }

    G4double GetXbeam() { return fXbeam; }

    G4double GetYbeam() { return fYbeam; }

    G4double GetPzz() { return fPzz; }

    G4int GetNreac() { return nreac; }

    G4int GetNp() { return np; }

    G4int GetEntry() { return entry; }

private:
    G4double fEgamma = 0.;
    G4double fXbeam = 0., fYbeam = 0.;
    G4double fPzz = 0.;

    G4int nreac = 0, np = 0, entry = 0;

    G4double ProCM = 0.;

};


#endif //COSMIC_EVENTINFO_HH
