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
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

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
        G4SDManager *SDmanp = G4SDManager::GetSDMpointer();

        if (fplastic_fat_nsys1HCID == -1) fplastic_fat_nsys1HCID = SDmanp->GetCollectionID("plastic_fat_nsys1HitsCollection");
        if (fplastic_fat_nsys2HCID == -1) fplastic_fat_nsys2HCID = SDmanp->GetCollectionID("plastic_fat_nsys2HitsCollection");
        if (fplastic_thin_nsys1HCID == -1) fplastic_thin_nsys1HCID = SDmanp->GetCollectionID("plastic_thin_nsys1HitsCollection");
        if (fplastic_thin_nsys2HCID == -1) fplastic_thin_nsys2HCID = SDmanp->GetCollectionID("plastic_thin_nsys2HitsCollection");
        if (fplastic_LQ_nsys1HCID == -1) fplastic_LQ_nsys1HCID = SDmanp->GetCollectionID("plastic_LQ_nsys1HitsCollection");
        if (fplastic_LQ_nsys2HCID == -1) fplastic_LQ_nsys2HCID = SDmanp->GetCollectionID("plastic_LQ_nsys2HitsCollection");

        if (fHadronCalorimeter_nsys1HCID == -1) fHadronCalorimeter_nsys1HCID = SDmanp->GetCollectionID("hadron_calorimeter_nsys1HitsCollection");
        if (fHadronCalorimeter_nsys2HCID == -1) fHadronCalorimeter_nsys2HCID = SDmanp->GetCollectionID("hadron_calorimeter_nsys2HitsCollection");


        if (fW_Chamber_nsys1HCID == -1) fW_Chamber_nsys1HCID = SDmanp->GetCollectionID("W_chamber_nsys1HitsCollection");
        if (fW_Chamber_nsys2HCID == -1) fW_Chamber_nsys2HCID = SDmanp->GetCollectionID("W_chamber_nsys2HitsCollection");
        if (fV_Chamber_nsys1HCID == -1) fV_Chamber_nsys1HCID = SDmanp->GetCollectionID("V_chamber_nsys1HitsCollection");
        if (fV_Chamber_nsys2HCID == -1) fV_Chamber_nsys2HCID = SDmanp->GetCollectionID("V_chamber_nsys2HitsCollection");

        auto event_info = event->GetUserInformation();

        // if (event_info == NULL) return;
        //  Eg=info->GetEgamma();
        //  nr=info->GetNreac();
        // nev=info->GetEntry();
        number_vertex = event->GetNumberOfPrimaryVertex();
        auto primary_vertex = event->GetPrimaryVertex();

        vertex_x = primary_vertex->GetX0();
        vertex_y = primary_vertex->GetY0();
        vertex_z = primary_vertex->GetZ0();
        vertex_index = primary_vertex->GetPrimary()->GetPDGcode();
        vertex_mass = primary_vertex->GetPrimary()->GetMass();
        vertex_energy = primary_vertex->GetPrimary()->GetKineticEnergy();
        vertex_momentum = primary_vertex->GetPrimary()->GetMomentum().mag();
        vertex_theta = primary_vertex->GetPrimary()->GetMomentum().theta();
        vertex_phi = primary_vertex->GetPrimary()->GetMomentum().phi();


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
        auto plastic_LQ_nsys1HC = GetHitsCollection(fplastic_LQ_nsys1HCID, event);
        auto plastic_LQ_nsys2HC = GetHitsCollection(fplastic_LQ_nsys2HCID, event);

        auto HadronCalorimeter_nsys1HC = GetHC(event, fHadronCalorimeter_nsys1HCID);
        auto HadronCalorimeter_nsys2HC = GetHC(event, fHadronCalorimeter_nsys2HCID);

        auto W_Chamber_nsys1HC = GetHC(event, fW_Chamber_nsys1HCID);
        auto W_Chamber_nsys2HC = GetHC(event, fW_Chamber_nsys2HCID);
        auto V_Chamber_nsys1HC = GetHC(event, fV_Chamber_nsys1HCID);
        auto V_Chamber_nsys2HC = GetHC(event, fV_Chamber_nsys2HCID);

        if (!plastic_fat_nsys1HC) return;
        if (!plastic_fat_nsys2HC) return;
        if (!plastic_thin_nsys1HC) return;
        if (!plastic_thin_nsys2HC) return;
        if (!plastic_LQ_nsys1HC) return;
        if (!plastic_LQ_nsys2HC) return;

        if (!HadronCalorimeter_nsys1HC) return;
        if (!HadronCalorimeter_nsys2HC) return;

        if (!W_Chamber_nsys1HC) return;
        if (!W_Chamber_nsys2HC) return;
        if (!V_Chamber_nsys1HC) return;
        if (!V_Chamber_nsys2HC) return;

        // Get hit with total values



        // нулевой индекс массива это полный Хит, далее это конкретные пластики
        PlasticHit *plastic_fat_nsys1Hit[fNofLayers_plastic_fat_nsys1 + 1];
        PlasticHit *plastic_fat_nsys2Hit[fNofLayers_plastic_fat_nsys2 + 1];
        PlasticHit *plastic_thin_nsys1Hit[fNofLayers_plastic_thin_nsys1 + 1];
        PlasticHit *plastic_thin_nsys2Hit[fNofLayers_plastic_thin_nsys2 + 1];
        PlasticHit *plastic_LQ_nsys1Hit[fNofLayers_plastic_LQ_nsys1 + 1];
        PlasticHit *plastic_LQ_nsys2Hit[fNofLayers_plastic_LQ_nsys2 + 1];


        HadronCalorimeterHit *HadronCalorimeter_nsys1Hit[fNofLayers_HadrtonCalorimeter_nsys1 + 1];
        HadronCalorimeterHit *HadronCalorimeter_nsys2Hit[fNofLayers_HadrtonCalorimeter_nsys2 + 1];

        plastic_fat_nsys1Hit[0] = (*plastic_fat_nsys1HC)[plastic_fat_nsys1HC->entries() - 1];
        plastic_fat_nsys2Hit[0] = (*plastic_fat_nsys2HC)[plastic_fat_nsys2HC->entries() - 1];
        plastic_thin_nsys1Hit[0] = (*plastic_thin_nsys1HC)[plastic_thin_nsys1HC->entries() - 1];
        plastic_thin_nsys2Hit[0] = (*plastic_thin_nsys2HC)[plastic_thin_nsys2HC->entries() - 1];
        plastic_LQ_nsys1Hit[0] = (*plastic_LQ_nsys1HC)[plastic_LQ_nsys1HC->entries() - 1];
        plastic_LQ_nsys2Hit[0] = (*plastic_LQ_nsys2HC)[plastic_LQ_nsys2HC->entries() - 1];

        //  HadronCalorimeter_nsys1HC->GetHit(0)->

        //   HadronCalorimeter_nsys1Hit[0] = static_cast<HadronCalorimeterHit*>(HadronCalorimeter_nsys1HC->GetHit(0));
        //   HadronCalorimeter_nsys2Hit[0] = static_cast<HadronCalorimeterHit*>(HadronCalorimeter_nsys2HC->GetHit(0));


        //    G4cerr <<"!!!!" << HadronCalorimeter_nsys1HC->GetSize()<<G4endl;
        //    G4cerr <<"!!!!" << HadronCalorimeter_nsys2HC->GetSize()<<G4endl;


        //auto HadronCalorimeter_nsys1Hit[fNofLayers_HadrtonCalorimeter_nsys1 + 1] = static_cast<HadronCalorimeterHit*>(HadronCalorimeter_nsys1HC);

        for (G4int i = 0; i <= fNofLayers_plastic_fat_nsys1; i++) {
            if (i >= 1) plastic_fat_nsys1Hit[i] = (*plastic_fat_nsys1HC)[i - 1];
            if (plastic_fat_nsys1Hit[i]->GetEdep() > plastic_fat_threshold) {
                fPlastic_fat_nsys1Edep[i] = plastic_fat_nsys1Hit[i]->GetEdep() / MeV;
                fPlastic_fat_nsys1LO[i] = plastic_fat_nsys1Hit[i]->GetLO() / MeV;
                fPlastic_fat_nsys1A1[i] = plastic_fat_nsys1Hit[i]->GetA1() / MeV;
                fPlastic_fat_nsys1A2[i] = plastic_fat_nsys1Hit[i]->GetA2() / MeV;
                fPlastic_fat_nsys1T1[i] = plastic_fat_nsys1Hit[i]->GetT1() / ns;
                fPlastic_fat_nsys1T2[i] = plastic_fat_nsys1Hit[i]->GetT2() / ns;
                fPlastic_fat_nsys1TrackLength[i] = plastic_fat_nsys1Hit[i]->GetTrackLength() / cm;
                fPlastic_fat_nsys1ToF[i] = plastic_fat_nsys1Hit[i]->GetToF() / ns;

                fPlastic_fat_nsys1XPos[i] = plastic_fat_nsys1Hit[i]->GetLocalPos().x() / cm;
                fPlastic_fat_nsys1YPos[i] = plastic_fat_nsys1Hit[i]->GetLocalPos().y() / cm;
                fPlastic_fat_nsys1ZPos[i] = plastic_fat_nsys1Hit[i]->GetLocalPos().z() / cm;
            }
        }

        for (G4int i = 0; i <= fNofLayers_plastic_fat_nsys2; i++) {
            if (i >= 1) plastic_fat_nsys2Hit[i] = (*plastic_fat_nsys2HC)[i - 1];
            if (plastic_fat_nsys2Hit[i]->GetEdep() > plastic_fat_threshold) {
                fPlastic_fat_nsys2Edep[i] = plastic_fat_nsys2Hit[i]->GetEdep() / MeV;
                fPlastic_fat_nsys2LO[i] = plastic_fat_nsys2Hit[i]->GetLO() / MeV;
                fPlastic_fat_nsys2A1[i] = plastic_fat_nsys2Hit[i]->GetA1() / MeV;
                fPlastic_fat_nsys2A2[i] = plastic_fat_nsys2Hit[i]->GetA2() / MeV;
                fPlastic_fat_nsys2T1[i] = plastic_fat_nsys2Hit[i]->GetT1() / ns;
                fPlastic_fat_nsys2T2[i] = plastic_fat_nsys2Hit[i]->GetT2() / ns;
                fPlastic_fat_nsys2TrackLength[i] = plastic_fat_nsys2Hit[i]->GetTrackLength() / cm;
                fPlastic_fat_nsys2ToF[i] = plastic_fat_nsys2Hit[i]->GetToF() / ns;

                fPlastic_fat_nsys2XPos[i] = plastic_fat_nsys2Hit[i]->GetLocalPos().x() / cm;
                fPlastic_fat_nsys2YPos[i] = plastic_fat_nsys2Hit[i]->GetLocalPos().y() / cm;
                fPlastic_fat_nsys2ZPos[i] = plastic_fat_nsys2Hit[i]->GetLocalPos().z() / cm;
            }
        }
/////////////////////

        for (G4int i = 0; i <= fNofLayers_plastic_thin_nsys1; i++) {
            if (i >= 1) plastic_thin_nsys1Hit[i] = (*plastic_thin_nsys1HC)[i - 1];
            if (plastic_thin_nsys1Hit[i]->GetEdep() > plastic_thin_threshold) {
                fPlastic_thin_nsys1Edep[i] = plastic_thin_nsys1Hit[i]->GetEdep() / MeV;
                fPlastic_thin_nsys1LO[i] = plastic_thin_nsys1Hit[i]->GetLO() / MeV;
                fPlastic_thin_nsys1A1[i] = plastic_thin_nsys1Hit[i]->GetA1() / MeV;
                fPlastic_thin_nsys1T1[i] = plastic_thin_nsys1Hit[i]->GetT1() / ns;
                fPlastic_thin_nsys1TrackLength[i] = plastic_thin_nsys1Hit[i]->GetTrackLength() / cm;
                fPlastic_thin_nsys1ToF[i] = plastic_thin_nsys1Hit[i]->GetToF() / ns;

                fPlastic_thin_nsys1XPos[i] = plastic_thin_nsys1Hit[i]->GetLocalPos().x() / cm;
                fPlastic_thin_nsys1YPos[i] = plastic_thin_nsys1Hit[i]->GetLocalPos().y() / cm;
                fPlastic_thin_nsys1ZPos[i] = plastic_thin_nsys1Hit[i]->GetLocalPos().z() / cm;
            }
        }
        for (G4int i = 0; i <= fNofLayers_plastic_thin_nsys2; i++) {
            if (i >= 1) plastic_thin_nsys2Hit[i] = (*plastic_thin_nsys2HC)[i - 1];
            if (plastic_thin_nsys2Hit[i]->GetEdep() > plastic_thin_threshold) {
                fPlastic_thin_nsys2Edep[i] = plastic_thin_nsys2Hit[i]->GetEdep() / MeV;
                fPlastic_thin_nsys2LO[i] = plastic_thin_nsys2Hit[i]->GetLO() / MeV;
                fPlastic_thin_nsys2A1[i] = plastic_thin_nsys2Hit[i]->GetA1() / MeV;
                fPlastic_thin_nsys2T1[i] = plastic_thin_nsys2Hit[i]->GetT1() / ns;
                fPlastic_thin_nsys2TrackLength[i] = plastic_thin_nsys2Hit[i]->GetTrackLength() / cm;
                fPlastic_thin_nsys2ToF[i] = plastic_thin_nsys2Hit[i]->GetToF() / ns;

                fPlastic_thin_nsys2XPos[i] = plastic_thin_nsys2Hit[i]->GetLocalPos().x() / cm;
                fPlastic_thin_nsys2YPos[i] = plastic_thin_nsys2Hit[i]->GetLocalPos().y() / cm;
                fPlastic_thin_nsys2ZPos[i] = plastic_thin_nsys2Hit[i]->GetLocalPos().z() / cm;
            }
        }

//////////////////

        for (G4int i = 0; i <= fNofLayers_plastic_LQ_nsys1; i++) {
            if (i >= 1) plastic_LQ_nsys1Hit[i] = (*plastic_LQ_nsys1HC)[i - 1];
            if (plastic_LQ_nsys1Hit[i]->GetEdep() > plastic_LQ_threshold) {
                fPlastic_LQ_nsys1Edep[i] = plastic_LQ_nsys1Hit[i]->GetEdep() / MeV;
                fPlastic_LQ_nsys1LO[i] = plastic_LQ_nsys1Hit[i]->GetLO() / MeV;
                fPlastic_LQ_nsys1A1[i] = plastic_LQ_nsys1Hit[i]->GetA1() / MeV;
                fPlastic_LQ_nsys1T1[i] = plastic_LQ_nsys1Hit[i]->GetT1() / ns;
                fPlastic_LQ_nsys1TrackLength[i] = plastic_LQ_nsys1Hit[i]->GetTrackLength() / cm;
                fPlastic_LQ_nsys1ToF[i] = plastic_LQ_nsys1Hit[i]->GetToF() / ns;

                fPlastic_LQ_nsys1XPos[i] = plastic_LQ_nsys1Hit[i]->GetLocalPos().x() / cm;
                fPlastic_LQ_nsys1YPos[i] = plastic_LQ_nsys1Hit[i]->GetLocalPos().y() / cm;
                fPlastic_LQ_nsys1ZPos[i] = plastic_LQ_nsys1Hit[i]->GetLocalPos().z() / cm;
            }
        }
        for (G4int i = 0; i <= fNofLayers_plastic_LQ_nsys2; i++) {
            if (i >= 1) plastic_LQ_nsys2Hit[i] = (*plastic_LQ_nsys2HC)[i - 1];
            if (plastic_LQ_nsys2Hit[i]->GetEdep() > plastic_LQ_threshold) {
                fPlastic_LQ_nsys2Edep[i] = plastic_LQ_nsys2Hit[i]->GetEdep() / MeV;
                fPlastic_LQ_nsys2LO[i] = plastic_LQ_nsys2Hit[i]->GetLO() / MeV;
                fPlastic_LQ_nsys2A1[i] = plastic_LQ_nsys2Hit[i]->GetA1() / MeV;
                fPlastic_LQ_nsys2T1[i] = plastic_LQ_nsys2Hit[i]->GetT1() / ns;
                fPlastic_LQ_nsys2TrackLength[i] = plastic_LQ_nsys2Hit[i]->GetTrackLength() / cm;
                fPlastic_LQ_nsys2ToF[i] = plastic_LQ_nsys2Hit[i]->GetToF() / ns;

                fPlastic_LQ_nsys2XPos[i] = plastic_LQ_nsys2Hit[i]->GetLocalPos().x() / cm;
                fPlastic_LQ_nsys2YPos[i] = plastic_LQ_nsys2Hit[i]->GetLocalPos().y() / cm;
                fPlastic_LQ_nsys2ZPos[i] = plastic_LQ_nsys2Hit[i]->GetLocalPos().z() / cm;
            }
        }
        
///////////////////////
        G4int index = 0;
        G4int k = 0;

        for (G4int i = 0; i < (G4int) HadronCalorimeter_nsys1HC->GetSize(); i++) {
            auto HadronCalorimeter_nsys1Hit_ = static_cast<HadronCalorimeterHit *>(HadronCalorimeter_nsys1HC->GetHit(
                    i));
            if (HadronCalorimeter_nsys1Hit_->GetEdep() > HadronCalorimeter_threshold) {

                // index = i -  AC_IND - DE_IND - HCX_IND;
                index = i;

                if (index >= 0 && index < N_HCX) {
                    k = index % NX_BARS;
                    k = (k % N_UNITS) * 2 + (k / N_UNITS);    // re-numbering bars in a layer

                    fHCX_nsys1AL[index] = k;
                    fHCX_nsys1N[index]++;
                    fHCX_nsys1Edep[index] = HadronCalorimeter_nsys1Hit_->GetEdep() / MeV;
                    fHCX_nsys1LO[index] = HadronCalorimeter_nsys1Hit_->GetLO() / MeV;
                    fHCX_nsys1A1[index] = HadronCalorimeter_nsys1Hit_->GetA1() / MeV;
                    fHCX_nsys1T1[index] = HadronCalorimeter_nsys1Hit_->GetT1() / ns;
                    fHCX_nsys1TrackLength[index] = HadronCalorimeter_nsys1Hit_->GetTrackLength() / cm;
                    fHCX_nsys1ToF[index] = HadronCalorimeter_nsys1Hit_->GetToF() / ns;

                    fHCX_nsys1XPos[index] = HadronCalorimeter_nsys1Hit_->GetLocalPos().x() / cm;
                    fHCX_nsys1YPos[index] = HadronCalorimeter_nsys1Hit_->GetLocalPos().y() / cm;
                    fHCX_nsys1ZPos[index] = HadronCalorimeter_nsys1Hit_->GetLocalPos().z() / cm;

                    continue;
                }

                index = index - HCZ_IND;
                if (index >= 0 && index < N_HCZ) {
                    k = index % NZ_BARS;
                    k = (k % N_UNITS) * 2 + (k / N_UNITS);    // re-numbering bars in a layer

                    fHCZ_nsys1AL[index] = k;
                    fHCZ_nsys1N[index]++;
                    fHCZ_nsys1Edep[index] = HadronCalorimeter_nsys1Hit_->GetEdep() / MeV;
                    fHCZ_nsys1LO[index] = HadronCalorimeter_nsys1Hit_->GetLO() / MeV;
                    fHCZ_nsys1A1[index] = HadronCalorimeter_nsys1Hit_->GetA1() / MeV;
                    fHCZ_nsys1T1[index] = HadronCalorimeter_nsys1Hit_->GetT1() / ns;
                    fHCZ_nsys1TrackLength[index] = HadronCalorimeter_nsys1Hit_->GetTrackLength() / cm;
                    fHCZ_nsys1ToF[index] = HadronCalorimeter_nsys1Hit_->GetToF() / ns;

                    fHCZ_nsys1XPos[index] = HadronCalorimeter_nsys1Hit_->GetLocalPos().x() / cm;
                    fHCZ_nsys1YPos[index] = HadronCalorimeter_nsys1Hit_->GetLocalPos().y() / cm;
                    fHCZ_nsys1ZPos[index] = HadronCalorimeter_nsys1Hit_->GetLocalPos().z() / cm;

                    continue;

                }
            }
        }


        index = 0;
        k = 0;

        for (G4int i = 0; i < (G4int) HadronCalorimeter_nsys2HC->GetSize(); i++) {
            auto HadronCalorimeter_nsys2Hit_ = static_cast<HadronCalorimeterHit *>(HadronCalorimeter_nsys2HC->GetHit(
                    i));
            if (HadronCalorimeter_nsys2Hit_->GetEdep() > HadronCalorimeter_threshold) {

                //  index = i -  AC_IND - DE_IND - HCX_IND;
                index = i;

                if (index >= 0 && index < N_HCX) {
                    k = index % NX_BARS;
                    k = (k % N_UNITS) * 2 + (k / N_UNITS);    // re-numbering bars in a layer

                    fHCX_nsys2AL[index] = k;
                    fHCX_nsys2N[index]++;
                    fHCX_nsys2Edep[index] = HadronCalorimeter_nsys2Hit_->GetEdep() / MeV;
                    fHCX_nsys2LO[index] = HadronCalorimeter_nsys2Hit_->GetLO() / MeV;
                    fHCX_nsys2A1[index] = HadronCalorimeter_nsys2Hit_->GetA1() / MeV;
                    fHCX_nsys2T1[index] = HadronCalorimeter_nsys2Hit_->GetT1() / ns;
                    fHCX_nsys2TrackLength[index] = HadronCalorimeter_nsys2Hit_->GetTrackLength() / cm;
                    fHCX_nsys2ToF[index] = HadronCalorimeter_nsys2Hit_->GetToF() / ns;

                    fHCX_nsys2XPos[index] = HadronCalorimeter_nsys2Hit_->GetLocalPos().x() / cm;
                    fHCX_nsys2YPos[index] = HadronCalorimeter_nsys2Hit_->GetLocalPos().y() / cm;
                    fHCX_nsys2ZPos[index] = HadronCalorimeter_nsys2Hit_->GetLocalPos().z() / cm;

                    continue;
                }

                index = index - HCZ_IND;
                if (index >= 0 && index < N_HCZ) {
                    k = index % NZ_BARS;
                    k = (k % N_UNITS) * 2 + (k / N_UNITS);    // re-numbering bars in a layer

                    fHCZ_nsys2AL[index] = k;
                    fHCZ_nsys2N[index]++;
                    fHCZ_nsys2Edep[i] = HadronCalorimeter_nsys2Hit_->GetEdep() / MeV;
                    fHCZ_nsys2LO[index] = HadronCalorimeter_nsys2Hit_->GetLO() / MeV;
                    fHCZ_nsys2A1[index] = HadronCalorimeter_nsys2Hit_->GetA1() / MeV;
                    fHCZ_nsys2T1[index] = HadronCalorimeter_nsys2Hit_->GetT1() / ns;
                    fHCZ_nsys2TrackLength[index] = HadronCalorimeter_nsys2Hit_->GetTrackLength() / cm;
                    fHCZ_nsys2ToF[index] = HadronCalorimeter_nsys2Hit_->GetToF() / ns;

                    fHCZ_nsys2XPos[index] = HadronCalorimeter_nsys2Hit_->GetLocalPos().x() / cm;
                    fHCZ_nsys2YPos[index] = HadronCalorimeter_nsys2Hit_->GetLocalPos().y() / cm;
                    fHCZ_nsys2ZPos[index] = HadronCalorimeter_nsys2Hit_->GetLocalPos().z() / cm;

                    continue;

                }
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


        // Fill histograms, ntuple
        //

        // get analysis manager
        //auto analysisManager = G4AnalysisManager::Instance();


        // fill histograms

        if (plastic_fat_nsys1Hit[0]->GetEdep() > plastic_fat_threshold) {
            analysisManager->FillH1(0, plastic_fat_nsys1Hit[0]->GetEdep() / MeV);
            analysisManager->FillH1(1, plastic_fat_nsys1Hit[0]->GetTrackLength() / cm);
        }

        if (plastic_fat_nsys2Hit[0]->GetEdep() > plastic_fat_threshold) {
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
        analysisManager->FillNtupleIColumn(0, 1);

        analysisManager->FillNtupleDColumn(1, plastic_fat_nsys1Hit[0]->GetEdep() / MeV);
        analysisManager->FillNtupleDColumn(2, plastic_fat_nsys1Hit[0]->GetTrackLength() / cm);
        analysisManager->FillNtupleDColumn(3, plastic_fat_nsys2Hit[0]->GetEdep() / MeV);
        analysisManager->FillNtupleDColumn(4, plastic_fat_nsys2Hit[0]->GetTrackLength() / cm);


        analysisManager->FillNtupleDColumn(5, vertex_x / cm);
        analysisManager->FillNtupleDColumn(6, vertex_y / cm);
        analysisManager->FillNtupleDColumn(7, vertex_z / cm);
        analysisManager->FillNtupleIColumn(8, vertex_index);
        analysisManager->FillNtupleIColumn(9, vertex_mass / MeV);
        analysisManager->FillNtupleDColumn(10, vertex_energy / MeV);
        analysisManager->FillNtupleDColumn(11, vertex_momentum / MeV);
        analysisManager->FillNtupleDColumn(12, vertex_theta / degree);
        analysisManager->FillNtupleDColumn(13, vertex_phi / degree);


        analysisManager->AddNtupleRow(0);
        //  }


    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
