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
/// \file CalorimeterSD.hh
/// \brief Definition of the Cosmic::CalorimeterSD class

#ifndef CosmicPlasticSD_h
#define CosmicPlasticSD_h 1

#include "G4VSensitiveDetector.hh"


#include "G4RunManager.hh"
#include "DetectorConstruction.hh"

#include <vector>

class DetectorConstruction;
class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

#include "PlasticHit.hh"


namespace Cosmic
{

/// Calorimeter sensitive detector class
///
/// In Initialize(), it creates one hit for each calorimeter layer and one more
/// hit for accounting the total quantities in all layers.
///
/// The values are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step.

class PlasticSD : public G4VSensitiveDetector
{
  public:
  //  PlasticSD(const G4String& name, const G4String& hitsCollectionName, G4int nofCells);
 // PlasticSD(const G4String& name, const G4String &hitsCollectionName, DetectorConstruction*);

      PlasticSD(const G4String &name, const G4String &hitsCollectionName,
                G4int nofLayers /*, /*Cosmic::#1#DetectorConstruction* detector*/);
      ~PlasticSD() override;

      // methods from base class
      void Initialize(G4HCofThisEvent *hitCollection) override;
      G4bool ProcessHits(G4Step *step, G4TouchableHistory *ROhist) override;

      void EndOfEvent(G4HCofThisEvent *hitCollection) override;

      void SetDiscrThres(G4double val) { discr_threshold = val; };

  private:
    PlasticHitsCollection *fHitsCollection = nullptr;
    G4int fHCID = -1;

    G4int fNofLayers = 0;

    DetectorConstruction *Detector;

    // DetectorConstruction* Detector;
    // G4int*                   HitID;

    G4double discr_threshold = 0.1 * MeV;
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

