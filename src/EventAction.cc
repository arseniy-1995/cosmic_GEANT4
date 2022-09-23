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
/// \brief Implementation of the Cosmic::EventAction class

#include "Constants.hh"
#include "EventAction.hh"

#include "PlasticSD.hh"
#include "HadronCalorimeterSD.hh"

#include "PlasticHit.hh"
#include "HadronCalorimeterHit.hh"

#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>


//using std::array;
//using std::vector;

namespace {

// Utility function which finds a hit collection with the given Id
// and print warnings if not found
    G4VHitsCollection *GetHC(const G4Event *event, G4int collId) {
        auto hce = event->GetHCofThisEvent();
        if (!hce) {
            G4ExceptionDescription msg;
            msg << "No hits collection of this event found." << G4endl;
            G4Exception("EventAction::EndOfEventAction()",
                        "Code001", JustWarning, msg);
            return nullptr;
        }

        auto hc = hce->GetHC(collId);
        if (!hc) {
            G4ExceptionDescription msg;
            msg << "Hits collection " << collId << " of this event not found." << G4endl;
            G4Exception("EventAction::EndOfEventAction()",
                        "Code001", JustWarning, msg);
        }
        return hc;
    }

}


namespace Cosmic {

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    EventAction::EventAction() {


        G4RunManager::GetRunManager()->SetPrintProgress(1);
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    EventAction::~EventAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    PlasticHitsCollection *
    EventAction::GetHitsCollection(G4int hcID, const G4Event *event) const {
        auto hitsCollection
                = static_cast<PlasticHitsCollection *>(
                        event->GetHCofThisEvent()->GetHC(hcID));

        if (!hitsCollection) {
            G4ExceptionDescription msg;
            msg << "Cannot access hitsCollection ID " << hcID;
            G4Exception("EventAction::GetHitsCollection()",
                        "MyCode0003", FatalException, msg);
        }

        return hitsCollection;
    }





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//void EventAction::PrintEventStatistics(G4double absoEdep, G4double absoTrackLength, G4double gapEdep, G4double gapTrackLength) const



    void EventAction::PrintEventStatistics(G4double plasticEdep, G4double plasticTrackLength) const {

        // print event statistics
        G4cout
                << "   plastic_fat_nsys1: total energy: "
                << std::setw(7) << G4BestUnit(plasticEdep, "Energy")
                << "       total track length: "
                << std::setw(7) << G4BestUnit(plasticTrackLength, "Length")
                << G4endl;


    }



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void EventAction::BeginOfEventAction(const G4Event * /*event*/) {

        // Find hit collections and histogram Ids by names (just once)
        // and save them in the data members of this class

        /*
        if (fHodHCID[0] == -1) {
            auto sdManager = G4SDManager::GetSDMpointer();
            auto analysisManager = G4AnalysisManager::Instance();

            // hits collections names
            array<G4String, kDim> hHCName
                    = {{ "hodoscope1/hodoscopeColl", "hodoscope2/hodoscopeColl" }};
            array<G4String, kDim> dHCName
                    = {{ "chamber1/driftChamberColl", "chamber2/driftChamberColl" }};
            array<G4String, kDim> cHCName
                    = {{ "EMcalorimeter/EMcalorimeterColl", "HadCalorimeter/HadCalorimeterColl" }};

            // histograms names
            array<array<G4String, kDim>, kDim> histoName
                    = {{ {{ "Chamber1", "Chamber2" }}, {{ "Chamber1 XY", "Chamber2 XY" }} }};

            for (G4int iDet = 0; iDet < kDim; ++iDet) {
                // hit collections IDs
                fHodHCID[iDet]   = sdManager->GetCollectionID(hHCName[iDet]);
                fDriftHCID[iDet] = sdManager->GetCollectionID(dHCName[iDet]);
                fCalHCID[iDet]   = sdManager->GetCollectionID(cHCName[iDet]);
                // histograms IDs
                fDriftHistoID[kH1][iDet] = analysisManager->GetH1Id(histoName[kH1][iDet]);
                fDriftHistoID[kH2][iDet] = analysisManager->GetH2Id(histoName[kH2][iDet]);
            }
        }
        */

        // Get hits collections IDs (only once)
        //   if (fplasticHCID == -1) {
        //      auto SDmanp = G4SDManager::GetSDMpointer();
        //      auto analysisManager = G4AnalysisManager::Instance();
        //      fplasticHCID = SDmanp->GetCollectionID("plasticHitsCollection");

        //  }


    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void EventAction::EndOfEventAction(const G4Event *event) {

        //G4cerr << "........" << G4endl;

        // Get hits collections IDs (only once)
        if (fplastic_fat_nsys1HCID == -1) {
            G4SDManager *SDmanp = G4SDManager::GetSDMpointer();
            fplastic_fat_nsys1HCID = SDmanp->GetCollectionID("plastic_fat_nsys1HitsCollection");
        }
        if (fplastic_fat_nsys2HCID == -1) {
            G4SDManager *SDmanp = G4SDManager::GetSDMpointer();
            fplastic_fat_nsys2HCID = SDmanp->GetCollectionID("plastic_fat_nsys2HitsCollection");
        }
        if (fplastic_thin_nsys1HCID == -1) {
            G4SDManager *SDmanp = G4SDManager::GetSDMpointer();
            fplastic_thin_nsys1HCID = SDmanp->GetCollectionID("plastic_thin_nsys1HitsCollection");
        }
        if (fplastic_thin_nsys2HCID == -1) {
            G4SDManager *SDmanp = G4SDManager::GetSDMpointer();
            fplastic_thin_nsys2HCID = SDmanp->GetCollectionID("plastic_thin_nsys2HitsCollection");
        }

        //
        // Fill histograms & ntuple
        //

        // Get analysis manager
        auto analysisManager = G4AnalysisManager::Instance();

        // Get hits collections
       // auto plasticHC = GetHitsCollection(fplasticHCID, event);


      // vector<G4double> totalPlasticEdep[10];
     //   vector<G4double> totalPlasticTrackLength[10];

        // Get hits collections
     //   auto plasticHC = GetHC(event, fplasticHCID);
        auto plastic_fat_nsys1HC = GetHitsCollection(fplastic_fat_nsys1HCID, event);
        auto plastic_fat_nsys2HC = GetHitsCollection(fplastic_fat_nsys2HCID, event);
        auto plastic_thin_nsys1HC = GetHitsCollection(fplastic_thin_nsys1HCID, event);
        auto plastic_thin_nsys2HC = GetHitsCollection(fplastic_thin_nsys2HCID, event);
        if ( ! plastic_fat_nsys1HC ) return;
        if ( ! plastic_fat_nsys2HC ) return;
        if ( ! plastic_thin_nsys1HC ) return;
        if ( ! plastic_thin_nsys2HC ) return;
        // Get hit with total values

        // нулевой индекс массива это полный Хит, далее это конкретные пластики
        PlasticHit *plastic_fat_nsys1Hit[fNofLayers_plastic_fat_nsys1 + 1];
        PlasticHit *plastic_fat_nsys2Hit[fNofLayers_plastic_fat_nsys2 + 1];
        PlasticHit *plastic_thin_nsys1Hit[fNofLayers_plastic_thin_nsys1 + 1];
        PlasticHit *plastic_thin_nsys2Hit[fNofLayers_plastic_thin_nsys2 + 1];

        plastic_fat_nsys1Hit[0] = (*plastic_fat_nsys1HC)[plastic_fat_nsys1HC->entries() - 1];
        plastic_fat_nsys2Hit[0] = (*plastic_fat_nsys2HC)[plastic_fat_nsys2HC->entries() - 1];
        plastic_thin_nsys1Hit[0] = (*plastic_thin_nsys1HC)[plastic_thin_nsys1HC->entries() - 1];
        plastic_thin_nsys2Hit[0] = (*plastic_thin_nsys2HC)[plastic_thin_nsys2HC->entries() - 1];

        if (plastic_fat_nsys1Hit[0]->GetEdep() > 0)
        {
            fPlastic_fat_nsys1Edep[0] = plastic_fat_nsys1Hit[0]->GetEdep() / MeV;
            fPlastic_fat_nsys1TrackLength[0] = plastic_fat_nsys1Hit[0]->GetTrackLength() / cm;
            fPlastic_fat_nsys1ToF[0] = plastic_fat_nsys1Hit[0]->GetToF() / ns;
        }
        if (plastic_fat_nsys2Hit[0]->GetEdep() > 0)
        {
            fPlastic_fat_nsys2Edep[0] = plastic_fat_nsys2Hit[0]->GetEdep() / MeV;
            fPlastic_fat_nsys2TrackLength[0] = plastic_fat_nsys2Hit[0]->GetTrackLength() / cm;
            fPlastic_fat_nsys2ToF[0] = plastic_fat_nsys2Hit[0]->GetToF() / ns;
        }
        if (plastic_thin_nsys1Hit[0]->GetEdep() > 0)
        {
            fPlastic_thin_nsys1Edep[0] = plastic_thin_nsys1Hit[0]->GetEdep() / MeV;
            fPlastic_thin_nsys1TrackLength[0] = plastic_thin_nsys1Hit[0]->GetTrackLength() / cm;
            fPlastic_thin_nsys1ToF[0] = plastic_thin_nsys1Hit[0]->GetToF() / ns;
        }
        if (plastic_thin_nsys2Hit[0]->GetEdep() > 0)
        {
            fPlastic_thin_nsys2Edep[0] = plastic_thin_nsys2Hit[0]->GetEdep() / MeV;
            fPlastic_thin_nsys2TrackLength[0] = plastic_thin_nsys2Hit[0]->GetTrackLength() / cm;
            fPlastic_thin_nsys2ToF[0] = plastic_thin_nsys2Hit[0]->GetToF() / ns;
        }

        G4double plastic_fat_threshold = 1.0*MeV; // порог записи в файл
        G4double plastic_thin_threshold = 0.1 * MeV;

        for (G4int i = 1; i <= fNofLayers_plastic_fat_nsys1; i++) {
            plastic_fat_nsys1Hit[i] = (*plastic_fat_nsys1HC)[i - 1];
            if (plastic_fat_nsys1Hit[i]->GetEdep() > plastic_fat_threshold) {
                fPlastic_fat_nsys1Edep[i] = plastic_fat_nsys1Hit[i]->GetEdep() / MeV;
                fPlastic_fat_nsys1TrackLength[i] = plastic_fat_nsys1Hit[i]->GetTrackLength() / cm;
                fPlastic_fat_nsys1ToF[i] = plastic_fat_nsys1Hit[i]->GetToF() / ns;
            }
        }

        for (G4int i = 1; i <= fNofLayers_plastic_fat_nsys2; i++) {
            plastic_fat_nsys2Hit[i] = (*plastic_fat_nsys2HC)[i - 1];
            if (plastic_fat_nsys2Hit[i]->GetEdep() > plastic_fat_threshold) {
                fPlastic_fat_nsys2Edep[i] = plastic_fat_nsys2Hit[i]->GetEdep() / MeV;
                fPlastic_fat_nsys2TrackLength[i] = plastic_fat_nsys2Hit[i]->GetTrackLength() / cm;
                fPlastic_fat_nsys2ToF[i] = plastic_fat_nsys2Hit[i]->GetToF() / ns;
            }
        }

        for(G4int i=1;i<=fNofLayers_plastic_thin_nsys1;i++) {
            plastic_thin_nsys1Hit[i] = (*plastic_thin_nsys1HC)[i-1];
            if (plastic_thin_nsys1Hit[i]->GetEdep()>plastic_thin_threshold){
                fPlastic_thin_nsys1Edep[i]= plastic_thin_nsys1Hit[i]->GetEdep() / MeV;
                fPlastic_thin_nsys1TrackLength[i]= plastic_thin_nsys1Hit[i]->GetTrackLength() / cm;
                fPlastic_thin_nsys1ToF[i] = plastic_thin_nsys1Hit[i]->GetToF() / ns;
            }
        }
        for(G4int i=1;i<=fNofLayers_plastic_thin_nsys2;i++) {
            plastic_thin_nsys2Hit[i] = (*plastic_thin_nsys2HC)[i-1];
            if (plastic_thin_nsys2Hit[i]->GetEdep()>plastic_thin_threshold) {
                fPlastic_thin_nsys2Edep[i]= plastic_thin_nsys2Hit[i]->GetEdep() / MeV;
                fPlastic_thin_nsys2TrackLength[i]= plastic_thin_nsys2Hit[i]->GetTrackLength() / cm;
                fPlastic_thin_nsys2ToF[i] = plastic_thin_nsys2Hit[i]->GetToF() / ns;
            }
        }


   //     totalPlasticEdep = 0.;
     //   totalPlasticTrackLength =0.0;
      //  for (unsigned long i = 0; i < plasticHC->GetSize(); ++i) {
      //      G4double edep = 0.;
            // The EM and Had calorimeter hits are of different types

        //        auto plasticHit = static_cast<PlasticHit*>(plasticHC->GetHit(i));
        //        edep = plasticHit->GetEdep();

         //   if ( edep > 0. ) {

        //        totalPlasticEdep[i] += edep;
        //    }
          //  fPlasticEdep[i] = edep;

           // fPlasticEdep[]
     //   }


      //  for(G4int i=1;i<=fNofLayers_plastic_fat_nsys1;i++) {
       //     G4cerr <<"!!!"  <<G4endl;
       //     G4cerr <<Edepplastic_fat_nsys1[i]  <<G4endl;
      //  }

       // auto plasticHit = static_cast<PlasticHit*>(plasticHC - 1);


        // Print per event (modulo n)
        //

        auto eventID = event->GetEventID();
        auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();


        if ((printModulo > 0) && (eventID % printModulo == 0)) {
            G4cout << "---> End of event: " << eventID << G4endl;

           PrintEventStatistics(
                   plastic_fat_nsys1Hit[0]->GetEdep(), plastic_fat_nsys1Hit[0]->GetTrackLength()
                   );


        }



// РАЧЕК

/*

        G4bool verb=false;//FALSE;

        G4int nvx = event->GetNumberOfPrimaryVertex();
        if(nvx<2) return;

        G4HCofThisEvent* HCE = event->GetHCofThisEvent();
        PlasticHitsCollection* DHC = 0;
        G4int n_hit = 0;
//  G4double blkE,blkT,blkL,totE=0;
        G4int i,k,n,nch,ip;

        DHC = 0;
        if(HCE) DHC = (PlasticHitsCollection*)(HCE->GetHC(fplasticHCID));
        else {
            if(verb)	G4cout << " HCE!" << G4endl;
            else	return;
        }	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (DHC==NULL)return;
        n_hit = DHC->entries();
        if(verb) G4cout << "n_hit=" << n_hit << G4endl;




        PlasticHit *hit;



        for(i=0;i<10;i++) {
            hit=(*DHC)[i];
            G4cerr<<hit->GetEdep()<<G4endl;

        }

*/

        // Fill histograms, ntuple
        //

        // get analysis manager
       //auto analysisManager = G4AnalysisManager::Instance();


        // fill histograms

        if(plastic_fat_nsys1Hit[0]->GetEdep()>plastic_fat_threshold) {
            analysisManager->FillH1(0, plastic_fat_nsys1Hit[0]->GetEdep() / MeV);
            analysisManager->FillH1(1, plastic_fat_nsys1Hit[0]->GetTrackLength() / cm);
        }

        if(plastic_fat_nsys2Hit[0]->GetEdep()>plastic_fat_threshold) {
            analysisManager->FillH1(2, plastic_fat_nsys2Hit[0]->GetEdep() / MeV);
            analysisManager->FillH1(3, plastic_fat_nsys2Hit[0]->GetTrackLength() / cm);
        }

        // fill ntuple


       // for (unsigned int i = 0; i<plasticHC->GetSize(); ++i) {
       //     auto hit = static_cast<PlasticHit*>(plasticHC->GetHit(i));
            // columns 0, 1
        //    analysisManager->FillNtupleDColumn(0, hit->GetEdep());
       //     analysisManager->FillNtupleDColumn(1, hit->GetTrackLength());
     //   }
     //   if(plastic_fat_nsys1Hit[0]->GetEdep()/MeV>1. && plastic_fat_nsys2Hit[0]->GetEdep()/MeV>1.) {
            analysisManager->FillNtupleDColumn(0, plastic_fat_nsys1Hit[0]->GetEdep() / MeV);
            analysisManager->FillNtupleDColumn(1, plastic_fat_nsys1Hit[0]->GetTrackLength() / cm);
            analysisManager->FillNtupleDColumn(2, plastic_fat_nsys2Hit[0]->GetEdep() / MeV);
            analysisManager->FillNtupleDColumn(3, plastic_fat_nsys2Hit[0]->GetTrackLength() / cm);

           // analysisManager->FillNtupleDColumn(4, *Edepplastic_fat_nsys1);
           // analysisManager->FillNtupleDColumn(5, *Edepplastic_fat_nsys2);

            analysisManager->AddNtupleRow(0);
      //  }


    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
