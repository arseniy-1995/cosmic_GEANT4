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
/// \file RunAction.cc
/// \brief Implementation of the Cosmic::RunAction class

#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "EventAction.hh"


namespace Cosmic
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(EventAction* eventAction)
: fEventAction(eventAction)
{
  // set printing event number per each event
 // G4RunManager::GetRunManager()->SetPrintProgress(1);

  // Create analysis manager
  // The choice of the output format is done via the specified
  // file extension.
 auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType("root");
    // If the filename extension is not provided, the default file type (root)
    // will be used for all files specified without extension.


  // Create directories
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
    // Note: merging ntuples is available only with Root output

    analysisManager->SetFileName("Cosmic");

  // Book histograms, ntuple
  //

  // Creating histograms


    analysisManager->CreateH1("Eplastic_fat_nsys1", "Edep in Plastic_fat_nsys1", 100, 0., 800 * MeV);
    analysisManager->CreateH1("Lplastic_fat_nsys1", "trackL in Plastic_fat_nsys1", 100, 0., 1 * m);
    analysisManager->CreateH1("Eplastic_fat_nsys2", "Edep in Plastic_fat_nsys2", 100, 0., 800 * MeV);
    analysisManager->CreateH1("Lplastic_fat_nsys2", "trackL in Plastic_fat_nsys2", 100, 0., 1 * m);

    analysisManager->CreateH1("Eplastic_thin_nsys1", "Edep in Plastic_thin_nsys1", 100, 0., 800 * MeV);
    analysisManager->CreateH1("Lplastic_thin_nsys1", "trackL in Plastic_thin_nsys1", 100, 0., 1 * m);
    analysisManager->CreateH1("Eplastic_thin_nsys2", "Edep in Plastic_thin_nsys2", 100, 0., 800 * MeV);
    analysisManager->CreateH1("Lplastic_thin_nsys2", "trackL in Plastic_thin_nsys2", 100, 0., 1 * m);


    // Creating ntuple
    //
    if ( fEventAction ) {
        analysisManager->CreateNtuple("Cosmic", "Edep and TrackL");

        analysisManager->CreateNtupleDColumn("Eplastic_fat_nsys1"); // column Id = 0
        analysisManager->CreateNtupleDColumn("Lplastic_fat_nsys1"); // column Id = 1
        analysisManager->CreateNtupleDColumn("Eplastic_fat_nsys2"); // column Id = 2
        analysisManager->CreateNtupleDColumn("Lplastic_fat_nsys2"); // column Id = 3

      //  analysisManager->CreateNtupleDColumn("Eplastic_fat_nsys1_vector", fEventAction->GetPlasticEdep()); // column Id = 4
      //  analysisManager->CreateNtupleDColumn("Eplastic_fat_nsys2_vector"); // column Id = 5

        analysisManager->CreateNtupleDColumn("Eplastic_fat_nsys1Vector", fEventAction->GetPlasticFatEdep(1)); // column Id = 4
        analysisManager->CreateNtupleDColumn("Eplastic_fat_nsys2Vector", fEventAction->GetPlasticFatEdep(2)); // column Id = 5
        analysisManager->CreateNtupleDColumn("Lplastic_fat_nsys1Vector", fEventAction->GetPlasticFatTrackLength(1)); // column Id = 6
        analysisManager->CreateNtupleDColumn("Lplastic_fat_nsys2Vector", fEventAction->GetPlasticFatTrackLength(2)); // column Id = 7
        analysisManager->CreateNtupleDColumn("ToFplastic_fat_nsys1Vector", fEventAction->GetPlasticFatToF(1)); // column Id = 8
        analysisManager->CreateNtupleDColumn("ToFplastic_fat_nsys2Vector", fEventAction->GetPlasticFatToF(2)); // column Id = 9

        analysisManager->CreateNtupleDColumn("Eplastic_thin_nsys1Vector", fEventAction->GetPlasticThinEdep(1)); // column Id = 10
        analysisManager->CreateNtupleDColumn("Eplastic_thin_nsys2Vector", fEventAction->GetPlasticThinEdep(2)); // column Id = 11
        analysisManager->CreateNtupleDColumn("Lplastic_thin_nsys1Vector", fEventAction->GetPlasticThinTrackLength(1)); // column Id = 12
        analysisManager->CreateNtupleDColumn("Lplastic_thin_nsys2Vector", fEventAction->GetPlasticThinTrackLength(2)); // column Id = 13
        analysisManager->CreateNtupleDColumn("ToFplastic_thin_nsys1Vector", fEventAction->GetPlasticThinToF(1)); // column Id = 14
        analysisManager->CreateNtupleDColumn("ToFplastic_thin_nsys2Vector", fEventAction->GetPlasticThinToF(2)); // column Id = 15

        analysisManager->FinishNtuple(0);
    }

// Set ntuple output file
    analysisManager->SetNtupleFileName(0, "Cosmicntuple");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

    // Reset histograms from previous run
   analysisManager->Reset();

  // Open an output file
  //
 // G4String fileName = "Cosmic.root";
  // Other supported output types:
  // G4String fileName = "Cosmic.csv";
  // G4String fileName = "Cosmic.hdf5";
  // G4String fileName = "Cosmic.xml";
 // analysisManager->OpenFile(fileName);

    // The default file name is set in RunAction::RunAction(),
    // it can be overwritten in a macro
    analysisManager->OpenFile();
  G4cout << "Using " << analysisManager->GetType() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  //
  auto analysisManager = G4AnalysisManager::Instance();



  if ( analysisManager->GetH1(0) ) {
    G4cout << G4endl << " ----> print histograms statistic ";
    if(isMaster) {
      G4cout << "for the entire run " << G4endl << G4endl;
    }
    else {
      G4cout << "for the local thread " << G4endl << G4endl;
    }


      G4cout << " E : mean = "
             << G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy")
             << " rms = "
             << G4BestUnit(analysisManager->GetH1(0)->rms(),  "Energy") << G4endl;


      G4cout << " L : mean = "
             << G4BestUnit(analysisManager->GetH1(1)->mean(), "Length")
             << " rms = "
             << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Length") << G4endl;



  }



  // save histograms & ntuple
  //
  analysisManager->Write();
 analysisManager->CloseFile(false);

    // Keep content of histos so that they are plotted.
    // The content will be reset at start of the next run.
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
