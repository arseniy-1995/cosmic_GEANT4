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
/// \file EventAction.cc
/// \brief Implementation of the Cosmic_sim::EventAction class

#include "EventAction.hh"
#include "PlasticSD.hh"
#include "PlasticHit.hh"

#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

namespace Cosmic_sim
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    PlasticHitsCollection*
EventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  auto hitsCollection
    = static_cast<PlasticHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));

  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID;
    G4Exception("EventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }

  return hitsCollection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//void EventAction::PrintEventStatistics(G4double absoEdep, G4double absoTrackLength, G4double gapEdep, G4double gapTrackLength) const

    void EventAction::PrintEventStatistics( G4double plastic_fat_nsys1Edep, G4double plastic_fat_nsys1TrackLength) const
    {
    /*
  // print event statistics
  G4cout
     << "   Absorber: total energy: "
     << std::setw(7) << G4BestUnit(absoEdep, "Energy")
     << "       total track length: "
     << std::setw(7) << G4BestUnit(absoTrackLength, "Length")
     << G4endl
     << "        Gap: total energy: "
     << std::setw(7) << G4BestUnit(gapEdep, "Energy")
     << "       total track length: "
     << std::setw(7) << G4BestUnit(gapTrackLength, "Length")
     << G4endl;

     */

    // print event statistics
    G4cout
            << "   plastic_fat_nsys1: total energy: "
            << std::setw(7) << G4BestUnit(plastic_fat_nsys1Edep, "Energy")
            << "       total track length: "
            << std::setw(7) << G4BestUnit(plastic_fat_nsys1TrackLength, "Length")
            << G4endl;



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
    /*
  // Get hits collections IDs (only once)
  if ( fAbsHCID == -1 ) {
    fAbsHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("AbsorberHitsCollection");
    fGapHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("GapHitsCollection");
  }
*/

    // Get hits collections IDs (only once)
    if ( fplastic_fat_nsys1HCID == -1 ) {
        fplastic_fat_nsys1HCID
                = G4SDManager::GetSDMpointer()->GetCollectionID("plastic_fat_nsys1HitsCollection");

    }


/*
  // Get hits collections
  auto absoHC = GetHitsCollection(fAbsHCID, event);
  auto gapHC = GetHitsCollection(fGapHCID, event);

  // Get hit with total values
  auto absoHit = (*absoHC)[absoHC->entries()-1];
  auto gapHit = (*gapHC)[gapHC->entries()-1];

  */


    // Get hits collections
    auto plastic_fat_nsys1HC = GetHitsCollection(fplastic_fat_nsys1HCID, event);

    // Get hit with total values
    auto plastic_fat_nsys1Hit = (*plastic_fat_nsys1HC)[plastic_fat_nsys1HC->entries()-1];


  // Print per event (modulo n)
  //
  auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;

    PrintEventStatistics(
    //  absoHit->GetEdep(), absoHit->GetTrackLength(),
    //  gapHit->GetEdep(), gapHit->GetTrackLength());

            plastic_fat_nsys1Hit->GetEdep(), plastic_fat_nsys1Hit->GetTrackLength());
  }

  // Fill histograms, ntuple
  //

  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  /*
  // fill histograms
  analysisManager->FillH1(0, absoHit->GetEdep());
  analysisManager->FillH1(1, gapHit->GetEdep());
  analysisManager->FillH1(2, absoHit->GetTrackLength());
  analysisManager->FillH1(3, gapHit->GetTrackLength());

  // fill ntuple
  analysisManager->FillNtupleDColumn(0, absoHit->GetEdep());
  analysisManager->FillNtupleDColumn(1, gapHit->GetEdep());
  analysisManager->FillNtupleDColumn(2, absoHit->GetTrackLength());
  analysisManager->FillNtupleDColumn(3, gapHit->GetTrackLength());
  analysisManager->AddNtupleRow();

   */


    // fill histograms
    analysisManager->FillH1(0, plastic_fat_nsys1Hit->GetEdep());
    analysisManager->FillH1(1, plastic_fat_nsys1Hit->GetTrackLength());


    // fill ntuple
    analysisManager->FillNtupleDColumn(0, plastic_fat_nsys1Hit->GetEdep());
    analysisManager->FillNtupleDColumn(1, plastic_fat_nsys1Hit->GetTrackLength());
    analysisManager->AddNtupleRow();


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
