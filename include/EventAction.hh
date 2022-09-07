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
/// \brief Definition of the Cosmic_sim::EventAction class

#ifndef Cosmic_simEventAction_h
#define Cosmic_simEventAction_h 1

#include "G4UserEventAction.hh"

#include "PlasticHit.hh"

#include "globals.hh"

namespace Cosmic_sim
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

private:
  // methods
  PlasticHitsCollection* GetHitsCollection(G4int hcID,
                                            const G4Event* event) const;
//  void PrintEventStatistics(G4double absoEdep, G4double absoTrackLength, G4double gapEdep, G4double gapTrackLength) const;

    void PrintEventStatistics(G4double aplastic_fat_nsys1Edep, G4double plastic_fat_nsys1TrackLength) const;



  // data members
 // G4int fAbsHCID = -1;
 // G4int fGapHCID = -1;
  G4int fplastic_fat_nsys1HCID = -1;


};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


