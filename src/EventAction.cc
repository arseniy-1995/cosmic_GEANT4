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
#include "EventInfo.hh"

#include "PlasticSD.hh"
#include "HadronCalorimeterSD.hh"
#include "ChamberSD.hh"

#include "PlasticHit.hh"
#include "HadronCalorimeterHit.hh"
#include "ChamberHit.hh"

#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"

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

    EventAction::EventAction() :
            drawFlag("all")
    //drawFlag("charged")
    //drawFlag("neutral")
    {

        G4RunManager::GetRunManager()->SetPrintProgress(1);

        G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
        G4ParticleDefinition *particle = particleTable->FindParticle("proton");
        base_code = particle->GetPDGEncoding();
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

        if (fplastic_fat_nsys1HCID == -1)
            fplastic_fat_nsys1HCID = SDmanp->GetCollectionID("plastic_fat_nsys1HitsCollection");
        if (fplastic_fat_nsys2HCID == -1)
            fplastic_fat_nsys2HCID = SDmanp->GetCollectionID("plastic_fat_nsys2HitsCollection");
        if (fplastic_thin_nsys1HCID == -1)
            fplastic_thin_nsys1HCID = SDmanp->GetCollectionID("plastic_thin_nsys1HitsCollection");
        if (fplastic_thin_nsys2HCID == -1)
            fplastic_thin_nsys2HCID = SDmanp->GetCollectionID("plastic_thin_nsys2HitsCollection");
        if (fplastic_LQ_nsys1HCID == -1)
            fplastic_LQ_nsys1HCID = SDmanp->GetCollectionID("plastic_LQ_nsys1HitsCollection");
        if (fplastic_LQ_nsys2HCID == -1)
            fplastic_LQ_nsys2HCID = SDmanp->GetCollectionID("plastic_LQ_nsys2HitsCollection");

        if (fHadronCalorimeter_nsys1HCID == -1)
            fHadronCalorimeter_nsys1HCID = SDmanp->GetCollectionID("hadron_calorimeter_nsys1HitsCollection");
        if (fHadronCalorimeter_nsys2HCID == -1)
            fHadronCalorimeter_nsys2HCID = SDmanp->GetCollectionID("hadron_calorimeter_nsys2HitsCollection");


        if (fW_Chamber_nsys1HCID == -1) fW_Chamber_nsys1HCID = SDmanp->GetCollectionID("W_Chamber_nsys1HitsCollection");
        if (fW_Chamber_nsys2HCID == -1) fW_Chamber_nsys2HCID = SDmanp->GetCollectionID("W_Chamber_nsys2HitsCollection");

#if defined(VCARM1) && defined(RUN21)
        if (fV_Chamber_nsys1HCID == -1) fV_Chamber_nsys1HCID = SDmanp->GetCollectionID("V_Chamber_nsys1HitsCollection");
#endif
#if defined(VCARM2) && defined(RUN21)
        if (fV_Chamber_nsys2HCID == -1) fV_Chamber_nsys2HCID = SDmanp->GetCollectionID("V_Chamber_nsys2HitsCollection");
#endif

        auto event_info = (EventInfo *) event->GetUserInformation();
        // EventInfo* event_info =(EventInfo*)event->GetUserInformation();

        if (event_info == NULL) return;


        G4int number_vertex = event->GetNumberOfPrimaryVertex();
        // if(number_vertex<2) return;
//TODO только для GENBOS

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

        vertex_Pzz = event_info->GetPzz();
        vertex_Energy_gamma = event_info->GetEgamma();
        vertex_number_reaction = event_info->GetNreac();
        vertex_number_event = event_info->GetEntry();

        G4int i_vertex, k_vertex, n_vertex, nch_vertex, ip_vertex;
        G4PrimaryVertex *primary_vertex_with_index;

        for (ip_vertex = 0; ip_vertex < number_vertex; ip_vertex++) {
            primary_vertex_with_index = event->GetPrimaryVertex(ip_vertex);
            if (primary_vertex_with_index->GetPrimary()->GetPDGcode() == base_code) break; // proton ?
        }

        if (ip_vertex == number_vertex) primary_vertex_with_index = event->GetPrimaryVertex(0);

        G4ThreeVector vertex_with_index;
        vertex_number_vertex = 0;
        if (ip_vertex < number_vertex) {
            primary_vertex_with_index = event->GetPrimaryVertex(ip_vertex);
            index_vertex[vertex_number_vertex] = primary_vertex_with_index->GetPrimary()->GetPDGcode();
            energy_vertex[vertex_number_vertex] = primary_vertex_with_index->GetPrimary()->GetKineticEnergy();
            vertex_with_index = primary_vertex_with_index->GetPrimary()->GetMomentum();
            theta_vertex[vertex_number_vertex] = vertex_with_index.theta();
            phi_vertex[vertex_number_vertex] = vertex_with_index.phi();

            vertex_x_vector[vertex_number_vertex] = primary_vertex_with_index->GetX0() / cm;
            vertex_y_vector[vertex_number_vertex] = primary_vertex_with_index->GetY0() / cm;
            vertex_z_vector[vertex_number_vertex] = primary_vertex_with_index->GetZ0() / cm;
            vertex_index_vector[vertex_number_vertex] = primary_vertex_with_index->GetPrimary()->GetPDGcode();
            vertex_mass_vector[vertex_number_vertex] = primary_vertex_with_index->GetPrimary()->GetMass() / MeV;
            vertex_energy_vector[vertex_number_vertex] =
                    primary_vertex_with_index->GetPrimary()->GetKineticEnergy() / MeV;
            vertex_momentum_vector[vertex_number_vertex] =
                    primary_vertex_with_index->GetPrimary()->GetMomentum().mag() / MeV;
            vertex_theta_vector[vertex_number_vertex] = vertex_with_index.theta() / degree;
            vertex_phi_vector[vertex_number_vertex] = vertex_with_index.phi() / degree;

            vertex_number_vertex++;

        }

        for (i_vertex = 0; i_vertex < number_vertex; i_vertex++) {
            if (i_vertex == ip_vertex)continue;
            primary_vertex_with_index = event->GetPrimaryVertex(i_vertex);
            index_vertex[vertex_number_vertex] = primary_vertex_with_index->GetPrimary()->GetPDGcode();
            energy_vertex[vertex_number_vertex] = primary_vertex_with_index->GetPrimary()->GetKineticEnergy();
            vertex_with_index = primary_vertex_with_index->GetPrimary()->GetMomentum();
            theta_vertex[vertex_number_vertex] = vertex_with_index.theta();
            phi_vertex[vertex_number_vertex] = vertex_with_index.phi();

            vertex_x_vector[vertex_number_vertex] = primary_vertex_with_index->GetX0() / cm;
            vertex_y_vector[vertex_number_vertex] = primary_vertex_with_index->GetY0() / cm;
            vertex_z_vector[vertex_number_vertex] = primary_vertex_with_index->GetZ0() / cm;
            vertex_index_vector[vertex_number_vertex] = primary_vertex_with_index->GetPrimary()->GetPDGcode();
            vertex_mass_vector[vertex_number_vertex] = primary_vertex_with_index->GetPrimary()->GetMass() / MeV;
            vertex_energy_vector[vertex_number_vertex] =
                    primary_vertex_with_index->GetPrimary()->GetKineticEnergy() / MeV;
            vertex_momentum_vector[vertex_number_vertex] =
                    primary_vertex_with_index->GetPrimary()->GetMomentum().mag() / MeV;
            vertex_theta_vector[vertex_number_vertex] = vertex_with_index.theta() / degree;
            vertex_phi_vector[vertex_number_vertex] = vertex_with_index.phi() / degree;


            vertex_number_vertex++;
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
        auto plastic_LQ_nsys1HC = GetHitsCollection(fplastic_LQ_nsys1HCID, event);
        auto plastic_LQ_nsys2HC = GetHitsCollection(fplastic_LQ_nsys2HCID, event);

        auto HadronCalorimeter_nsys1HC = GetHC(event, fHadronCalorimeter_nsys1HCID);
        auto HadronCalorimeter_nsys2HC = GetHC(event, fHadronCalorimeter_nsys2HCID);

        auto W_Chamber_nsys1HC = GetHC(event, fW_Chamber_nsys1HCID);
        auto W_Chamber_nsys2HC = GetHC(event, fW_Chamber_nsys2HCID);

#if defined(VCARM1) && defined(RUN21)
        auto V_Chamber_nsys1HC = GetHC(event, fV_Chamber_nsys1HCID);
        auto V_Chamber_nsys2HC = GetHC(event, fV_Chamber_nsys2HCID);
#endif
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

#if defined(VCARM1) && defined(RUN21)
        if (!V_Chamber_nsys1HC) return;
        if (!V_Chamber_nsys2HC) return;
#endif
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

        ChamberHit *W_Chamber_nsys1Hit[fNofLayers_W_Chamber_nsys1 + 1];
        ChamberHit *W_Chamber_nsys2Hit[fNofLayers_W_Chamber_nsys2 + 1];
        ChamberHit *V_Chamber_nsys1Hit[fNofLayers_V_Chamber_nsys1 + 1];
        ChamberHit *V_Chamber_nsys2Hit[fNofLayers_V_Chamber_nsys2 + 1];

        plastic_fat_nsys1Hit[0] = (*plastic_fat_nsys1HC)[plastic_fat_nsys1HC->entries() - 1];
        plastic_fat_nsys2Hit[0] = (*plastic_fat_nsys2HC)[plastic_fat_nsys2HC->entries() - 1];
        plastic_thin_nsys1Hit[0] = (*plastic_thin_nsys1HC)[plastic_thin_nsys1HC->entries() - 1];
        plastic_thin_nsys2Hit[0] = (*plastic_thin_nsys2HC)[plastic_thin_nsys2HC->entries() - 1];
        plastic_LQ_nsys1Hit[0] = (*plastic_LQ_nsys1HC)[plastic_LQ_nsys1HC->entries() - 1];
        plastic_LQ_nsys2Hit[0] = (*plastic_LQ_nsys2HC)[plastic_LQ_nsys2HC->entries() - 1];

        //
        // For trigger
        G4int hcdep[2] = {0, 0}; // адронный калориметр
        G4int wcdep[2] = {0, 0}; // камеры
        G4int de_fat_dep[2] = {0, 0}; // толстые счетчики
        G4int de_thin_dep[2] = {0, 0}; // тонкие счетчики
        G4int de_LQ_dep[2] = {0, 0}; // счетчики LQ

        //  HadronCalorimeter_nsys1HC->GetHit(0)->

        //   HadronCalorimeter_nsys1Hit[0] = static_cast<HadronCalorimeterHit*>(HadronCalorimeter_nsys1HC->GetHit(0));
        //   HadronCalorimeter_nsys2Hit[0] = static_cast<HadronCalorimeterHit*>(HadronCalorimeter_nsys2HC->GetHit(0));


        //    G4cerr <<"!!!!" << HadronCalorimeter_nsys1HC->GetSize()<<G4endl;
        //    G4cerr <<"!!!!" << HadronCalorimeter_nsys2HC->GetSize()<<G4endl;


        //auto HadronCalorimeter_nsys1Hit[fNofLayers_HadrtonCalorimeter_nsys1 + 1] = static_cast<HadronCalorimeterHit*>(HadronCalorimeter_nsys1HC);

        for (G4int i = 0; i <= fNofLayers_plastic_fat_nsys1; i++) {
            if (i >= 1) plastic_fat_nsys1Hit[i] = (*plastic_fat_nsys1HC)[i - 1];
            if (plastic_fat_nsys1Hit[i]->GetEdep() > plastic_fat_threshold) {
                fPlastic_fatTrackID[0][i] = plastic_fat_nsys1Hit[i]->GetTrackID();
                fPlastic_fatEdep[0][i] = plastic_fat_nsys1Hit[i]->GetEdep() / MeV;
                fPlastic_fatLO[0][i] = plastic_fat_nsys1Hit[i]->GetLO() / MeV;
                fPlastic_fatA1[0][i] = plastic_fat_nsys1Hit[i]->GetA1() / MeV;
                fPlastic_fatA2[0][i] = plastic_fat_nsys1Hit[i]->GetA2() / MeV;
                fPlastic_fatT1[0][i] = plastic_fat_nsys1Hit[i]->GetT1() / ns;
                fPlastic_fatT2[0][i] = plastic_fat_nsys1Hit[i]->GetT2() / ns;
                fPlastic_fatTrackLength[0][i] = plastic_fat_nsys1Hit[i]->GetTrackLength() / cm;
                fPlastic_fatToF[0][i] = plastic_fat_nsys1Hit[i]->GetToF() / ns;

                fPlastic_fatXPos[0][i] = plastic_fat_nsys1Hit[i]->GetLocalPos().x() / cm;
                fPlastic_fatYPos[0][i] = plastic_fat_nsys1Hit[i]->GetLocalPos().y() / cm;
                fPlastic_fatZPos[0][i] = plastic_fat_nsys1Hit[i]->GetLocalPos().z() / cm;
                fPlastic_fat_global_XPos[0][i] = plastic_fat_nsys1Hit[i]->GetWorldPos().x() / cm;
                fPlastic_fat_global_YPos[0][i] = plastic_fat_nsys1Hit[i]->GetWorldPos().y() / cm;
                fPlastic_fat_global_ZPos[0][i] = plastic_fat_nsys1Hit[i]->GetWorldPos().z() / cm;
                fPlastic_fatTheta[0][i] = plastic_fat_nsys1Hit[i]->GetPosTheta() / degree;
                fPlastic_fatPhi[0][i] = plastic_fat_nsys1Hit[i]->GetPosPhi() / degree;

                de_fat_dep[0] = 1;
            }
        }

        for (G4int i = 0; i <= fNofLayers_plastic_fat_nsys2; i++) {
            if (i >= 1) plastic_fat_nsys2Hit[i] = (*plastic_fat_nsys2HC)[i - 1];
            if (plastic_fat_nsys2Hit[i]->GetEdep() > plastic_fat_threshold) {
                fPlastic_fatTrackID[1][i] = plastic_fat_nsys2Hit[i]->GetTrackID();
                fPlastic_fatEdep[1][i] = plastic_fat_nsys2Hit[i]->GetEdep() / MeV;
                fPlastic_fatLO[1][i] = plastic_fat_nsys2Hit[i]->GetLO() / MeV;
                fPlastic_fatA1[1][i] = plastic_fat_nsys2Hit[i]->GetA1() / MeV;
                fPlastic_fatA2[1][i] = plastic_fat_nsys2Hit[i]->GetA2() / MeV;
                fPlastic_fatT1[1][i] = plastic_fat_nsys2Hit[i]->GetT1() / ns;
                fPlastic_fatT2[1][i] = plastic_fat_nsys2Hit[i]->GetT2() / ns;
                fPlastic_fatTrackLength[1][i] = plastic_fat_nsys2Hit[i]->GetTrackLength() / cm;
                fPlastic_fatToF[1][i] = plastic_fat_nsys2Hit[i]->GetToF() / ns;

                fPlastic_fatXPos[1][i] = plastic_fat_nsys2Hit[i]->GetLocalPos().x() / cm;
                fPlastic_fatYPos[1][i] = plastic_fat_nsys2Hit[i]->GetLocalPos().y() / cm;
                fPlastic_fatZPos[1][i] = plastic_fat_nsys2Hit[i]->GetLocalPos().z() / cm;
                fPlastic_fat_global_XPos[1][i] = plastic_fat_nsys2Hit[i]->GetWorldPos().x() / cm;
                fPlastic_fat_global_YPos[1][i] = plastic_fat_nsys2Hit[i]->GetWorldPos().y() / cm;
                fPlastic_fat_global_ZPos[1][i] = plastic_fat_nsys2Hit[i]->GetWorldPos().z() / cm;
                fPlastic_fatTheta[1][i] = plastic_fat_nsys2Hit[i]->GetPosTheta() / degree;
                fPlastic_fatPhi[1][i] = plastic_fat_nsys2Hit[i]->GetPosPhi() / degree;

                de_fat_dep[1] = 1;
            }
        }
/////////////////////

        for (G4int i = 0; i <= fNofLayers_plastic_thin_nsys1; i++) {
            if (i >= 1) plastic_thin_nsys1Hit[i] = (*plastic_thin_nsys1HC)[i - 1];
            if (plastic_thin_nsys1Hit[i]->GetEdep() > plastic_thin_threshold) {
                fPlastic_thinTrackID[0][i] = plastic_thin_nsys1Hit[i]->GetTrackID();
                fPlastic_thinEdep[0][i] = plastic_thin_nsys1Hit[i]->GetEdep() / MeV;
                fPlastic_thinLO[0][i] = plastic_thin_nsys1Hit[i]->GetLO() / MeV;
                fPlastic_thinA1[0][i] = plastic_thin_nsys1Hit[i]->GetA1() / MeV;
                fPlastic_thinT1[0][i] = plastic_thin_nsys1Hit[i]->GetT1() / ns;
                fPlastic_thinTrackLength[0][i] = plastic_thin_nsys1Hit[i]->GetTrackLength() / cm;
                fPlastic_thinToF[0][i] = plastic_thin_nsys1Hit[i]->GetToF() / ns;

                fPlastic_thinXPos[0][i] = plastic_thin_nsys1Hit[i]->GetLocalPos().x() / cm;
                fPlastic_thinYPos[0][i] = plastic_thin_nsys1Hit[i]->GetLocalPos().y() / cm;
                fPlastic_thinZPos[0][i] = plastic_thin_nsys1Hit[i]->GetLocalPos().z() / cm;
                fPlastic_thin_global_XPos[0][i] = plastic_thin_nsys1Hit[i]->GetWorldPos().x() / cm;
                fPlastic_thin_global_YPos[0][i] = plastic_thin_nsys1Hit[i]->GetWorldPos().y() / cm;
                fPlastic_thin_global_ZPos[0][i] = plastic_thin_nsys1Hit[i]->GetWorldPos().z() / cm;
                fPlastic_thinTheta[0][i] = plastic_thin_nsys1Hit[i]->GetPosTheta() / degree;
                fPlastic_thinPhi[0][i] = plastic_thin_nsys1Hit[i]->GetPosPhi() / degree;

                de_thin_dep[0] = 1;
            }
        }
        for (G4int i = 0; i <= fNofLayers_plastic_thin_nsys2; i++) {
            if (i >= 1) plastic_thin_nsys2Hit[i] = (*plastic_thin_nsys2HC)[i - 1];
            if (plastic_thin_nsys2Hit[i]->GetEdep() > plastic_thin_threshold) {
                fPlastic_thinTrackID[1][i] = plastic_thin_nsys2Hit[i]->GetTrackID();
                fPlastic_thinEdep[1][i] = plastic_thin_nsys2Hit[i]->GetEdep() / MeV;
                fPlastic_thinLO[1][i] = plastic_thin_nsys2Hit[i]->GetLO() / MeV;
                fPlastic_thinA1[1][i] = plastic_thin_nsys2Hit[i]->GetA1() / MeV;
                fPlastic_thinT1[1][i] = plastic_thin_nsys2Hit[i]->GetT1() / ns;
                fPlastic_thinTrackLength[1][i] = plastic_thin_nsys2Hit[i]->GetTrackLength() / cm;
                fPlastic_thinToF[1][i] = plastic_thin_nsys2Hit[i]->GetToF() / ns;

                fPlastic_thinXPos[1][i] = plastic_thin_nsys2Hit[i]->GetLocalPos().x() / cm;
                fPlastic_thinYPos[1][i] = plastic_thin_nsys2Hit[i]->GetLocalPos().y() / cm;
                fPlastic_thinZPos[1][i] = plastic_thin_nsys2Hit[i]->GetLocalPos().z() / cm;
                fPlastic_thin_global_XPos[1][i] = plastic_thin_nsys2Hit[i]->GetWorldPos().x() / cm;
                fPlastic_thin_global_YPos[1][i] = plastic_thin_nsys2Hit[i]->GetWorldPos().y() / cm;
                fPlastic_thin_global_ZPos[1][i] = plastic_thin_nsys2Hit[i]->GetWorldPos().z() / cm;
                fPlastic_thinTheta[1][i] = plastic_thin_nsys2Hit[i]->GetPosTheta() / degree;
                fPlastic_thinPhi[1][i] = plastic_thin_nsys2Hit[i]->GetPosPhi() / degree;

                de_thin_dep[1] = 1;
            }
        }

//////////////////

        for (G4int i = 0; i <= fNofLayers_plastic_LQ_nsys1; i++) {
            if (i >= 1) plastic_LQ_nsys1Hit[i] = (*plastic_LQ_nsys1HC)[i - 1];
            if (plastic_LQ_nsys1Hit[i]->GetEdep() > plastic_LQ_threshold) {
                fPlastic_LQTrackID[0][i] = plastic_LQ_nsys1Hit[i]->GetTrackID();
                fPlastic_LQ_Edep[0][i] = plastic_LQ_nsys1Hit[i]->GetEdep() / MeV;
                fPlastic_LQ_LO[0][i] = plastic_LQ_nsys1Hit[i]->GetLO() / MeV;
                fPlastic_LQ_A1[0][i] = plastic_LQ_nsys1Hit[i]->GetA1() / MeV;
                fPlastic_LQ_T1[0][i] = plastic_LQ_nsys1Hit[i]->GetT1() / ns;
                fPlastic_LQ_TrackLength[0][i] = plastic_LQ_nsys1Hit[i]->GetTrackLength() / cm;
                fPlastic_LQ_ToF[0][i] = plastic_LQ_nsys1Hit[i]->GetToF() / ns;

                fPlastic_LQ_XPos[0][i] = plastic_LQ_nsys1Hit[i]->GetLocalPos().x() / cm;
                fPlastic_LQ_YPos[0][i] = plastic_LQ_nsys1Hit[i]->GetLocalPos().y() / cm;
                fPlastic_LQ_ZPos[0][i] = plastic_LQ_nsys1Hit[i]->GetLocalPos().z() / cm;
                fPlastic_LQ_global_XPos[0][i] = plastic_LQ_nsys1Hit[i]->GetWorldPos().x() / cm;
                fPlastic_LQ_global_YPos[0][i] = plastic_LQ_nsys1Hit[i]->GetWorldPos().y() / cm;
                fPlastic_LQ_global_ZPos[0][i] = plastic_LQ_nsys1Hit[i]->GetWorldPos().z() / cm;
                fPlastic_LQ_Theta[0][i] = plastic_LQ_nsys1Hit[i]->GetPosTheta() / degree;
                fPlastic_LQ_Phi[0][i] = plastic_LQ_nsys1Hit[i]->GetPosPhi() / degree;

                de_LQ_dep[0] = 1;
            }
        }
        for (G4int i = 0; i <= fNofLayers_plastic_LQ_nsys2; i++) {
            if (i >= 1) plastic_LQ_nsys2Hit[i] = (*plastic_LQ_nsys2HC)[i - 1];
            if (plastic_LQ_nsys2Hit[i]->GetEdep() > plastic_LQ_threshold) {
                fPlastic_LQTrackID[1][i] = plastic_LQ_nsys2Hit[i]->GetTrackID();
                fPlastic_LQ_Edep[1][i] = plastic_LQ_nsys2Hit[i]->GetEdep() / MeV;
                fPlastic_LQ_LO[1][i] = plastic_LQ_nsys2Hit[i]->GetLO() / MeV;
                fPlastic_LQ_A1[1][i] = plastic_LQ_nsys2Hit[i]->GetA1() / MeV;
                fPlastic_LQ_T1[1][i] = plastic_LQ_nsys2Hit[i]->GetT1() / ns;
                fPlastic_LQ_TrackLength[1][i] = plastic_LQ_nsys2Hit[i]->GetTrackLength() / cm;
                fPlastic_LQ_ToF[1][i] = plastic_LQ_nsys2Hit[i]->GetToF() / ns;

                fPlastic_LQ_XPos[1][i] = plastic_LQ_nsys2Hit[i]->GetLocalPos().x() / cm;
                fPlastic_LQ_YPos[1][i] = plastic_LQ_nsys2Hit[i]->GetLocalPos().y() / cm;
                fPlastic_LQ_ZPos[1][i] = plastic_LQ_nsys2Hit[i]->GetLocalPos().z() / cm;
                fPlastic_LQ_global_XPos[1][i] = plastic_LQ_nsys2Hit[i]->GetWorldPos().x() / cm;
                fPlastic_LQ_global_YPos[1][i] = plastic_LQ_nsys2Hit[i]->GetWorldPos().y() / cm;
                fPlastic_LQ_global_ZPos[1][i] = plastic_LQ_nsys2Hit[i]->GetWorldPos().z() / cm;
                fPlastic_LQ_Theta[1][i] = plastic_LQ_nsys2Hit[i]->GetPosTheta() / degree;
                fPlastic_LQ_Phi[1][i] = plastic_LQ_nsys2Hit[i]->GetPosPhi() / degree;

                de_LQ_dep[1] = 1;
            }
        }

/////////////////////// КАЛОРИМЕТР
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
                    fHCX_AL[0][index] = k;
                    fHCX_N[0][index]++;
                    fHCX_Edep[0][index] = HadronCalorimeter_nsys1Hit_->GetEdep() / MeV;
                    fHCX_LO[0][index] = HadronCalorimeter_nsys1Hit_->GetLO() / MeV;
                    fHCX_A[0][index] += HadronCalorimeter_nsys1Hit_->GetLO() / MeV;
                    fHCX_A1[0][index] = HadronCalorimeter_nsys1Hit_->GetA1() / MeV;
                    fHCX_A2[0][index] = HadronCalorimeter_nsys1Hit_->GetA2() / MeV;
                    fHCX_T1[0][index] = HadronCalorimeter_nsys1Hit_->GetT1() / ns;
                    fHCX_T2[0][index] = HadronCalorimeter_nsys1Hit_->GetT2() / ns;
                    fHCX_TrackLength[0][index] = HadronCalorimeter_nsys1Hit_->GetTrackLength() / cm;
                    fHCX_ToF[0][index] = HadronCalorimeter_nsys1Hit_->GetToF() / ns;

                    fHCX_XPos[0][index] = HadronCalorimeter_nsys1Hit_->GetLocalPos().x() / cm;
                    fHCX_YPos[0][index] = HadronCalorimeter_nsys1Hit_->GetLocalPos().y() / cm;
                    fHCX_ZPos[0][index] = HadronCalorimeter_nsys1Hit_->GetLocalPos().z() / cm;
                    hcdep[0] |= 1;
                    continue;
                }

                index = index - HCZ_IND;
                if (index >= 0 && index < N_HCZ) {
                    k = index % NZ_BARS;
                    k = (k % N_UNITS) * 2 + (k / N_UNITS);    // re-numbering bars in a layer

                    fHCZ_AL[0][index] = k;
                    fHCZ_N[0][index]++;
                    fHCZ_Edep[0][index] = HadronCalorimeter_nsys1Hit_->GetEdep() / MeV;
                    fHCZ_LO[0][index] = HadronCalorimeter_nsys1Hit_->GetLO() / MeV;
                    fHCZ_A[0][index] += HadronCalorimeter_nsys1Hit_->GetLO() / MeV;
                    fHCZ_A1[0][index] = HadronCalorimeter_nsys1Hit_->GetA1() / MeV;
                    fHCZ_A2[0][index] = HadronCalorimeter_nsys1Hit_->GetA2() / MeV;
                    fHCZ_T1[0][index] = HadronCalorimeter_nsys1Hit_->GetT1() / ns;
                    fHCZ_T2[0][index] = HadronCalorimeter_nsys1Hit_->GetT2() / ns;

                    fHCZ_TrackLength[0][index] = HadronCalorimeter_nsys1Hit_->GetTrackLength() / cm;
                    fHCZ_ToF[0][index] = HadronCalorimeter_nsys1Hit_->GetToF() / ns;

                    fHCZ_XPos[0][index] = HadronCalorimeter_nsys1Hit_->GetLocalPos().x() / cm;
                    fHCZ_YPos[0][index] = HadronCalorimeter_nsys1Hit_->GetLocalPos().y() / cm;
                    fHCZ_ZPos[0][index] = HadronCalorimeter_nsys1Hit_->GetLocalPos().z() / cm;
                    hcdep[0] |= 2;
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

                    fHCX_AL[1][index] = k;
                    fHCX_N[1][index]++;
                    fHCX_Edep[1][index] = HadronCalorimeter_nsys2Hit_->GetEdep() / MeV;
                    fHCX_LO[1][index] = HadronCalorimeter_nsys2Hit_->GetLO() / MeV;
                    fHCX_A[1][index] += HadronCalorimeter_nsys2Hit_->GetLO() / MeV;
                    fHCX_A1[1][index] = HadronCalorimeter_nsys2Hit_->GetA1() / MeV;
                    fHCX_A2[1][index] = HadronCalorimeter_nsys2Hit_->GetA2() / MeV;
                    fHCX_T1[1][index] = HadronCalorimeter_nsys2Hit_->GetT1() / ns;
                    fHCX_T2[1][index] = HadronCalorimeter_nsys2Hit_->GetT2() / ns;

                    fHCX_TrackLength[1][index] = HadronCalorimeter_nsys2Hit_->GetTrackLength() / cm;
                    fHCX_ToF[1][index] = HadronCalorimeter_nsys2Hit_->GetToF() / ns;

                    fHCX_XPos[1][index] = HadronCalorimeter_nsys2Hit_->GetLocalPos().x() / cm;
                    fHCX_YPos[1][index] = HadronCalorimeter_nsys2Hit_->GetLocalPos().y() / cm;
                    fHCX_ZPos[1][index] = HadronCalorimeter_nsys2Hit_->GetLocalPos().z() / cm;
                    hcdep[1] |= 1;
                    continue;
                }
                index = index - HCZ_IND;
                if (index >= 0 && index < N_HCZ) {
                    k = index % NZ_BARS;
                    k = (k % N_UNITS) * 2 + (k / N_UNITS);    // re-numbering bars in a layer

                    fHCZ_AL[1][index] = k;
                    fHCZ_N[1][index]++;
                    fHCZ_Edep[1][i] = HadronCalorimeter_nsys2Hit_->GetEdep() / MeV;
                    fHCZ_LO[1][index] = HadronCalorimeter_nsys2Hit_->GetLO() / MeV;
                    fHCZ_A[1][index] += HadronCalorimeter_nsys2Hit_->GetLO() / MeV;
                    fHCZ_A1[1][index] = HadronCalorimeter_nsys2Hit_->GetA1() / MeV;
                    fHCZ_A2[1][index] = HadronCalorimeter_nsys2Hit_->GetA2() / MeV;
                    fHCZ_T1[1][index] = HadronCalorimeter_nsys2Hit_->GetT1() / ns;
                    fHCZ_T2[1][index] = HadronCalorimeter_nsys2Hit_->GetT2() / ns;

                    fHCZ_TrackLength[1][index] = HadronCalorimeter_nsys2Hit_->GetTrackLength() / cm;
                    fHCZ_ToF[1][index] = HadronCalorimeter_nsys2Hit_->GetToF() / ns;

                    fHCZ_XPos[1][index] = HadronCalorimeter_nsys2Hit_->GetLocalPos().x() / cm;
                    fHCZ_YPos[1][index] = HadronCalorimeter_nsys2Hit_->GetLocalPos().y() / cm;
                    fHCZ_ZPos[1][index] = HadronCalorimeter_nsys2Hit_->GetLocalPos().z() / cm;
                    hcdep[1] |= 2;
                    continue;
                }
            }
        }

        /// ВЕРШИННЫЕ КАМЕРЫ

#if defined(VCARM1) && defined(RUN21)

        for (G4int i = 0; i < (G4int) V_Chamber_nsys1HC->GetSize(); i++) {
            auto V_Chamber_nsys1Hit_ = static_cast<ChamberHit *>(V_Chamber_nsys1HC->GetHit(i));
            if (V_Chamber_nsys1Hit_->GetEdep() > VChamber_threshold) {
                index = i;
                if (index >= 0 && index < NVC_WRS) {
                    fVC_N[0][index]++;
                    fVC_XPos[0][index] = V_Chamber_nsys1Hit_->GetLocalPos().x() / cm;
                    fVC_YPos[0][index] = V_Chamber_nsys1Hit_->GetLocalPos().y() / cm;
                    fVC_ZPos[0][index] = V_Chamber_nsys1Hit_->GetLocalPos().z() / cm;
                    fVC_RPos[0][index] = V_Chamber_nsys1Hit_->GetRho().z() / cm;
                    fVC_Theta[0][index] = V_Chamber_nsys1Hit_->GetPosTheta() / degree;
                    fVC_Phi[0][index] = V_Chamber_nsys1Hit_->GetPosPhi() / degree;
                    fVC_Mass[0][index] = V_Chamber_nsys1Hit_->GetMass() / MeV;
                    fVC_KineticEnergy[0][index] = V_Chamber_nsys1Hit_->GetKineticEnergy() / MeV;
                    wcdep[0] |= 1;
                    continue;
                }
            }
        }

        for (G4int i = 0; i < (G4int) V_Chamber_nsys2HC->GetSize(); i++) {
            auto V_Chamber_nsys2Hit_ = static_cast<ChamberHit *>(V_Chamber_nsys2HC->GetHit(i));
            if (V_Chamber_nsys2Hit_->GetEdep() > VChamber_threshold) {
                index = i;
                if (index >= 0 && index < NVC_WRS) {
                    fVC_N[1][index]++;
                    fVC_XPos[1][index] = V_Chamber_nsys2Hit_->GetLocalPos().x() / cm;
                    fVC_YPos[1][index] = V_Chamber_nsys2Hit_->GetLocalPos().y() / cm;
                    fVC_ZPos[1][index] = V_Chamber_nsys2Hit_->GetLocalPos().z() / cm;
                    fVC_RPos[1][index] = V_Chamber_nsys2Hit_->GetRho().z() / cm;
                    fVC_Theta[1][index] = V_Chamber_nsys2Hit_->GetPosTheta() / degree;
                    fVC_Phi[1][index] = V_Chamber_nsys2Hit_->GetPosPhi() / degree;
                    fVC_Mass[1][index] = V_Chamber_nsys2Hit_->GetMass() / MeV;
                    fVC_KineticEnergy[1][index] = V_Chamber_nsys2Hit_->GetKineticEnergy() / MeV;
                    wcdep[1] |= 1;
                    continue;
                }
            }
        }

#endif

        /////////////////////// ДРЕЙФОВЫЕ КАМЕРЫ

        for (G4int i = 0; i < (G4int) W_Chamber_nsys1HC->GetSize(); i++) {
            auto W_Chamber_nsys1Hit_ = static_cast<ChamberHit *>(W_Chamber_nsys1HC->GetHit(i));
            if (W_Chamber_nsys1Hit_->GetEdep() > WChamber_threshold) {
                index = i;
                if (index >= 0 && index < NW1_WRS) {
                    fWa_N[0][index]++;
                    fWa_XPos[0][index] = W_Chamber_nsys1Hit_->GetLocalPos().x() / cm;
                    fWa_YPos[0][index] = W_Chamber_nsys1Hit_->GetLocalPos().y() / cm;
                    fWa_ZPos[0][index] = W_Chamber_nsys1Hit_->GetLocalPos().z() / cm;
                    fWa_RPos[0][index] = W_Chamber_nsys1Hit_->GetRho().z() / cm;
                    fWa_Theta[0][index] = W_Chamber_nsys1Hit_->GetPosTheta() / degree;
                    fWa_Phi[0][index] = W_Chamber_nsys1Hit_->GetPosPhi() / degree;
                    fWa_Mass[0][index] = W_Chamber_nsys1Hit_->GetMass() / MeV;
                    fWa_KineticEnergy[0][index] = W_Chamber_nsys1Hit_->GetKineticEnergy() / MeV;
                    wcdep[0] |= 2;
                    continue;
                }
                index = index - WC2_IND;
                if (index >= 0 && index < NW2_WRS) {
                    fWb_N[0][index]++;
                    fWb_XPos[0][index] = W_Chamber_nsys1Hit_->GetLocalPos().x() / cm;
                    fWb_YPos[0][index] = W_Chamber_nsys1Hit_->GetLocalPos().y() / cm;
                    fWb_ZPos[0][index] = W_Chamber_nsys1Hit_->GetLocalPos().z() / cm;
                    fWb_RPos[0][index] = W_Chamber_nsys1Hit_->GetRho().x() / cm;
                    fWb_Theta[0][index] = W_Chamber_nsys1Hit_->GetPosTheta() / degree;
                    fWb_Phi[0][index] = W_Chamber_nsys1Hit_->GetPosPhi() / degree;
                    fWb_Mass[0][index] = W_Chamber_nsys1Hit_->GetMass() / MeV;
                    fWb_KineticEnergy[0][index] = W_Chamber_nsys1Hit_->GetKineticEnergy() / MeV;
                    wcdep[0] |= 4;
                    continue;
                }
                index = index - WC3_IND;
                if (index >= 0 && index < NW3_WRS) {
                    fWc_N[0][index]++;
                    fWc_XPos[0][index] = W_Chamber_nsys1Hit_->GetLocalPos().x() / cm;
                    fWc_YPos[0][index] = W_Chamber_nsys1Hit_->GetLocalPos().y() / cm;
                    fWc_ZPos[0][index] = W_Chamber_nsys1Hit_->GetLocalPos().z() / cm;
                    fWc_RPos[0][index] = W_Chamber_nsys1Hit_->GetRho().z() / cm;
                    fWc_Theta[0][index] = W_Chamber_nsys1Hit_->GetPosTheta() / degree;
                    fWc_Phi[0][index] = W_Chamber_nsys1Hit_->GetPosPhi() / degree;
                    fWc_Mass[0][index] = W_Chamber_nsys1Hit_->GetMass() / MeV;
                    fWc_KineticEnergy[0][index] = W_Chamber_nsys1Hit_->GetKineticEnergy() / MeV;
                    wcdep[0] |= 8;
                    continue;
                }
            }
        }


        for (G4int i = 0; i < (G4int) W_Chamber_nsys2HC->GetSize(); i++) {
            auto W_Chamber_nsys2Hit_ = static_cast<ChamberHit *>(W_Chamber_nsys2HC->GetHit(i));
            if (W_Chamber_nsys2Hit_->GetEdep() > WChamber_threshold) {
                index = i;
                if (index >= 0 && index < NW1_WRS) {
                    fWa_N[1][index]++;
                    fWa_XPos[1][index] = W_Chamber_nsys2Hit_->GetLocalPos().x() / cm;
                    fWa_YPos[1][index] = W_Chamber_nsys2Hit_->GetLocalPos().y() / cm;
                    fWa_ZPos[1][index] = W_Chamber_nsys2Hit_->GetLocalPos().z() / cm;
                    fWa_RPos[1][index] = W_Chamber_nsys2Hit_->GetRho().z() / cm;
                    fWa_Theta[1][index] = W_Chamber_nsys2Hit_->GetPosTheta() / degree;
                    fWa_Phi[1][index] = W_Chamber_nsys2Hit_->GetPosPhi() / degree;
                    fWa_Mass[1][index] = W_Chamber_nsys2Hit_->GetMass() / MeV;
                    fWa_KineticEnergy[1][index] = W_Chamber_nsys2Hit_->GetKineticEnergy() / MeV;
                    wcdep[1] |= 2;
                    continue;
                }
                index = index - WC2_IND;
                if (index >= 0 && index < NW2_WRS) {
                    fWb_N[1][index]++;
                    fWb_XPos[1][index] = W_Chamber_nsys2Hit_->GetLocalPos().x() / cm;
                    fWb_YPos[1][index] = W_Chamber_nsys2Hit_->GetLocalPos().y() / cm;
                    fWb_ZPos[1][index] = W_Chamber_nsys2Hit_->GetLocalPos().z() / cm;
                    fWb_RPos[1][index] = W_Chamber_nsys2Hit_->GetRho().x() / cm;
                    fWb_Theta[1][index] = W_Chamber_nsys2Hit_->GetPosTheta() / degree;
                    fWb_Phi[1][index] = W_Chamber_nsys2Hit_->GetPosPhi() / degree;
                    fWb_Mass[1][index] = W_Chamber_nsys2Hit_->GetMass() / MeV;
                    fWb_KineticEnergy[1][index] = W_Chamber_nsys2Hit_->GetKineticEnergy() / MeV;
                    wcdep[1] |= 4;
                    continue;
                }
                index = index - WC3_IND;
                if (index >= 0 && index < NW3_WRS) {
                    fWc_N[1][index]++;
                    fWc_XPos[1][index] = W_Chamber_nsys2Hit_->GetLocalPos().x() / cm;
                    fWc_YPos[1][index] = W_Chamber_nsys2Hit_->GetLocalPos().y() / cm;
                    fWc_ZPos[1][index] = W_Chamber_nsys2Hit_->GetLocalPos().z() / cm;
                    fWc_RPos[1][index] = W_Chamber_nsys2Hit_->GetRho().z() / cm;
                    fWc_Theta[1][index] = W_Chamber_nsys2Hit_->GetPosTheta() / degree;
                    fWc_Phi[1][index] = W_Chamber_nsys2Hit_->GetPosPhi() / degree;
                    fWc_Mass[1][index] = W_Chamber_nsys2Hit_->GetMass() / MeV;
                    fWc_KineticEnergy[1][index] = W_Chamber_nsys2Hit_->GetKineticEnergy() / MeV;
                    wcdep[1] |= 8;
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

        if (vertex_energy > 0) analysisManager->FillH1(8, vertex_energy / MeV);
        if (vertex_Energy_gamma > 0) analysisManager->FillH1(9, vertex_Energy_gamma / MeV);

        // fill ntuple


        // for (unsigned int i = 0; i<plasticHC->GetSize(); ++i) {
        //     auto hit = static_cast<PlasticHit*>(plasticHC->GetHit(i));
        // columns 0, 1
        //    analysisManager->FillNtupleDColumn(0, hit->GetEdep());
        //     analysisManager->FillNtupleDColumn(1, hit->GetTrackLength());
        //   }
        //   if(plastic_fat_nsys1Hit[0]->GetEdep()/MeV>1. && plastic_fat_nsys2Hit[0]->GetEdep()/MeV>1.) {
        analysisManager->FillNtupleIColumn(0, vertex_number_event);

        analysisManager->FillNtupleDColumn(1, plastic_fat_nsys1Hit[0]->GetEdep() / MeV);
        analysisManager->FillNtupleDColumn(2, plastic_fat_nsys1Hit[0]->GetTrackLength() / cm);
        analysisManager->FillNtupleDColumn(3, plastic_fat_nsys2Hit[0]->GetEdep() / MeV);
        analysisManager->FillNtupleDColumn(4, plastic_fat_nsys2Hit[0]->GetTrackLength() / cm);


        //    analysisManager->FillNtupleDColumn(5, vertex_x / cm);
        //    analysisManager->FillNtupleDColumn(6, vertex_y / cm);
        //    analysisManager->FillNtupleDColumn(7, vertex_z / cm);
        //    analysisManager->FillNtupleIColumn(8, vertex_index);
        //   analysisManager->FillNtupleIColumn(9, vertex_mass / MeV);
        //   analysisManager->FillNtupleDColumn(10, vertex_energy / MeV);
        //   analysisManager->FillNtupleDColumn(11, vertex_momentum / MeV);
        //   analysisManager->FillNtupleDColumn(12, vertex_theta / degree);
        //   analysisManager->FillNtupleDColumn(13, vertex_phi / degree);

        analysisManager->FillNtupleDColumn(5, vertex_Pzz);
        analysisManager->FillNtupleDColumn(6, vertex_Energy_gamma / MeV);
        analysisManager->FillNtupleDColumn(7, vertex_number_reaction);
        analysisManager->FillNtupleDColumn(8, vertex_number_event);



/// Trigger записи в файл
        if (analysisManager) {
            if (1
                && ((hcdep[0] == 3 && hcdep[1] == 3)    // hit in both HC arms
                    || ((hcdep[0] == 3) && de_fat_dep[1]) || ((hcdep[1] == 3) && de_fat_dep[0]))
                //     && ((ahx[0]>HC_THRES && ahz[0]>HC_THRES)		// HC hit in arm 1
                //	|| (ahx[1]>HC_THRES && ahz[1]>HC_THRES))	// HC hit in arm 2
                && (wcdep[0] >= 14 || wcdep[1] >= 14)    // track at least in one arm (VX not mandatory)
                // || (( de_thin_dep[0] && de_LQ_dep[1]) || (de_thin_dep[1] && de_LQ_dep[0])) // for LQ
                //  || (( de_thin_dep[0] && de_fat_dep[0]) && (de_thin_dep[1] && de_fat_dep[1])) // for cosmic
                    ) {
//G4cout << "Taken!"<<G4endl;
                //  TO->Fill();
                analysisManager->AddNtupleRow(0);
//G4cout<<" way="<<way[0]<<G4endl;
            } else {
//G4cout 	<< "Rejected! : ahx="<<ahx[0]<<","<<ahx[1]<<"  ahz="<<ahz[0]<<","<<ahz[1]
//	<< " wcdep="<<wcdep[0]<<","<<wcdep[1]<<G4endl;

            }
        }

        // Это без триггера, все. Если нет тригера по энергии, то пишется нулями
        // analysisManager->AddNtupleRow(0);

        //  }


        // extract the trajectories and draw them

        if (G4VVisManager::GetConcreteInstance()) {
            G4TrajectoryContainer *trajectoryContainer = event->GetTrajectoryContainer();
            G4int n_trajectories = 0;
            if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

            for (G4int i = 0; i < n_trajectories; i++) {
                G4Trajectory *trj = (G4Trajectory *)
                        ((*(event->GetTrajectoryContainer()))[i]);
                if (drawFlag == "all") trj->DrawTrajectory();
                else if ((drawFlag == "charged") && (trj->GetCharge() != 0.))
                    trj->DrawTrajectory();
                else if ((drawFlag == "neutral") && (trj->GetCharge() == 0.))
                    trj->DrawTrajectory();
            }
        }


    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
