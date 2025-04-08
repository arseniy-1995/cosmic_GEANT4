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

#include "Constants.hh"
#include "EventAction.hh"
#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Timer.hh"
#include "G4UnitsTable.hh"
#include "HistoManager.hh"

#ifdef GENBOS

#include "Genbos.hh"
#include "G4AutoLock.hh"

#endif

//namespace {
//    G4Mutex aMutex = G4MUTEX_INITIALIZER;
//
//}

namespace Cosmic {

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    RunAction::RunAction(EventAction *eventAction)
            : fEventAction(eventAction), FileNum(0) {
        timer = new G4Timer;

        fHistoManager = new HistoManager(fEventAction); // задание гистограмм и ntuples

        // G4AutoLock lock(&aMutex);
        //genbos_start_(&FileNum);

        // set printing event number per each event
        // G4RunManager::GetRunManager()->SetPrintProgress(1);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
    delete timer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void RunAction::BeginOfRunAction(const G4Run *aRun) {

        G4cerr << "### Run " << aRun->GetRunID() << " start." << G4endl;
        timer->Start();

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
       // if (analysisManager->IsActive()) {
            analysisManager->OpenFile();
        //}
  G4cout << "Using " << analysisManager->GetType() << G4endl;
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void RunAction::EndOfRunAction(const G4Run *aRun) {
        timer->Stop();

        // print histogram statistics
        //
        auto analysisManager = G4AnalysisManager::Instance();


        if (analysisManager->GetH1(0)) {
            G4cout << G4endl << " ----> print histograms statistic ";
            if (isMaster) {
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
       // if (analysisManager->IsActive()) {
          //  analysisManager->Write();
         //   analysisManager->CloseFile(false);
        // }
        // its true fir new versin GEANT

        analysisManager->Write();
        analysisManager->CloseFile(false);
        // Keep content of histos so that they are plotted.
        // The content will be reset at start of the next run.

        G4cerr << "number of event = " << aRun->GetNumberOfEvent() << " " << *timer << G4endl;


    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
