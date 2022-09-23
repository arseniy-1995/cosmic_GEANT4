//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file EventAction.hh
/// \brief Definition of the Cosmic::EventAction class

#ifndef CosmicEventAction_h
#define CosmicEventAction_h 1

#include "Constants.hh"
#include "G4UserEventAction.hh"

#include "PlasticHit.hh"
#include "HadronCalorimeterHit.hh"

#include "globals.hh"

#include <vector>
#include <array>

const G4int NofLayers_plastic_fat_nsys1 = 6;
const G4int NofLayers_plastic_fat_nsys2 = 8;
const G4int NofLayers_plastic_thin_nsys1 = 2;
const G4int NofLayers_plastic_thin_nsys2 = 2;


namespace Cosmic
{

/// Event action class
///
/// In EndOfEventAction(), it prints the accumulated quantities of the energy
/// deposit and track lengths of charged particles in Absober and Gap layers
/// stored in the hits collections.

class EventAction : public G4UserEventAction
{
public:
  EventAction();
  ~EventAction() override;

  void  BeginOfEventAction(const G4Event* event) override;
  void    EndOfEventAction(const G4Event* event) override;

    std::vector<G4double> &GetPlasticFatEdep(G4double nsys = 1) {
        if (nsys == 2) return fPlastic_fat_nsys2Edep;
        return fPlastic_fat_nsys1Edep;
    }

    std::vector<G4double> &GetPlasticFatTrackLength(G4double nsys = 1) {
        if (nsys == 2) return fPlastic_fat_nsys2TrackLength;
        return fPlastic_fat_nsys1TrackLength;
    }

    std::vector<G4double> &GetPlasticFatToF(G4double nsys = 1) {
        if (nsys == 2) return fPlastic_fat_nsys2ToF;
        return fPlastic_fat_nsys1ToF;
    }

    std::vector<G4double> &GetPlasticThinEdep(G4double nsys = 1) {
        if (nsys == 2) return fPlastic_thin_nsys2Edep;
        return fPlastic_thin_nsys1Edep;
    }

    std::vector<G4double> &GetPlasticThinTrackLength(G4double nsys = 1) {
        if (nsys == 2) return fPlastic_thin_nsys2TrackLength;
        return fPlastic_thin_nsys1TrackLength;
    }

    std::vector<G4double> &GetPlasticThinToF(G4double nsys = 1) {
        if (nsys == 2) return fPlastic_thin_nsys2ToF;
        return fPlastic_thin_nsys1ToF;
    }


private:
    // methods
    PlasticHitsCollection *GetHitsCollection(G4int hcID,const G4Event *event) const;

        void PrintEventStatistics(G4double aplasticEdep, G4double plasticTrackLength) const;


    std::vector<G4double> fPlastic_fat_nsys1Edep {std::vector<G4double>(fNofLayers_plastic_fat_nsys1 + 1,0.0)};
    std::vector<G4double> fPlastic_fat_nsys2Edep {std::vector<G4double>(fNofLayers_plastic_fat_nsys2 + 1,0.0)};
    std::vector<G4double> fPlastic_fat_nsys1TrackLength {std::vector<G4double>(fNofLayers_plastic_fat_nsys1 + 1,0.0)};
    std::vector<G4double> fPlastic_fat_nsys2TrackLength {std::vector<G4double>(fNofLayers_plastic_fat_nsys2 + 1,0.0)};
    std::vector<G4double> fPlastic_fat_nsys1ToF {std::vector<G4double>(fNofLayers_plastic_fat_nsys1 + 1,0.0)};
    std::vector<G4double> fPlastic_fat_nsys2ToF {std::vector<G4double>(fNofLayers_plastic_fat_nsys2 + 1,0.0)};

    std::vector<G4double> fPlastic_thin_nsys1Edep {std::vector<G4double>(fNofLayers_plastic_thin_nsys1 + 1,0.0)};
    std::vector<G4double> fPlastic_thin_nsys2Edep {std::vector<G4double>(fNofLayers_plastic_thin_nsys2 + 1,0.0)};
    std::vector<G4double> fPlastic_thin_nsys1TrackLength {std::vector<G4double>(fNofLayers_plastic_thin_nsys1 + 1,0.0)};
    std::vector<G4double> fPlastic_thin_nsys2TrackLength {std::vector<G4double>(fNofLayers_plastic_thin_nsys2 + 1,0.0)};
    std::vector<G4double> fPlastic_thin_nsys1ToF {std::vector<G4double>(fNofLayers_plastic_thin_nsys1 + 1,0.0)};
    std::vector<G4double> fPlastic_thin_nsys2ToF {std::vector<G4double>(fNofLayers_plastic_thin_nsys2 + 1,0.0)};

  // data members

  G4int fplastic_fat_nsys1HCID = -1;
  G4int fplastic_fat_nsys2HCID = -1;
  G4int fplastic_thin_nsys1HCID = -1;
  G4int fplastic_thin_nsys2HCID = -1;

};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


