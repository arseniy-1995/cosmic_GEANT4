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
/// \file optical/OpNovice2/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "EventAction.hh"
#include "G4UnitsTable.hh"
#include "Constants.hh"

namespace Cosmic
{
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    HistoManager::HistoManager(EventAction *eventAction)
      : fFileName("Cosmic"), fEventAction(eventAction)
    {
        Book();
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    HistoManager::~HistoManager() {}

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void HistoManager::Book() {



        //
        //
        // // auto analysisManager = G4AnalysisManager::Instance();
        //
        //
        // // Creating ntuple
        // //
        //
        // analysisMan->SetFirstNtupleId(1);
        // auto fCurrentNtupleId1 = analysisMan->CreateNtuple("TO", "Edep and TrackL");
        // auto fCurrentNtupleId2 = analysisMan->CreateNtuple("TO_light", "Edep and TrackL");
        //
        // //G4cerr << "!!!" << fCurrentNtupleId1 << " " << fCurrentNtupleId2 << std::endl;
        //

        // // G4cerr << "CurrentNtupleId = " << fCurrentNtupleId << G4endl;




// Create or get analysis manager
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

        analysisManager->SetFileName(fFileName);
         //analysisMan->SetActivation(true);  // enable inactivation of histograms

        // Book histograms, ntuple
        //

        // Creating histograms

        //
        // // Define histograms
        //
        //
        // Create all histograms as inactivated
        // for(G4int i = 0; i < analysisMan->GetNofH1s(); ++i)
        // {
        //   analysisMan->SetH1Activation(i, true);
        // }

        analysisManager->CreateH1("Eplastic_fat_nsys1", "Edep in Plastic_fat_nsys1", 100, 0., 800 * MeV); //0
        analysisManager->CreateH1("Lplastic_fat_nsys1", "trackL in Plastic_fat_nsys1", 100, 0., 1 * m); //1
        analysisManager->CreateH1("Eplastic_fat_nsys2", "Edep in Plastic_fat_nsys2", 100, 0., 800 * MeV); //2
        analysisManager->CreateH1("Lplastic_fat_nsys2", "trackL in Plastic_fat_nsys2", 100, 0., 1 * m); //3

        analysisManager->CreateH1("Eplastic_thin_nsys1", "Edep in Plastic_thin_nsys1", 100, 0., 800 * MeV); //4
        analysisManager->CreateH1("Lplastic_thin_nsys1", "trackL in Plastic_thin_nsys1", 100, 0., 1 * m); //5
        analysisManager->CreateH1("Eplastic_thin_nsys2", "Edep in Plastic_thin_nsys2", 100, 0., 800 * MeV); //6
        analysisManager->CreateH1("Lplastic_thin_nsys2", "trackL in Plastic_thin_nsys2", 100, 0., 1 * m); //7

        analysisManager->CreateH1("EVertex0", "Initial Energy Vertex", 2000, 0., 10000 * MeV); //8
        analysisManager->CreateH1("Egamma0", "Initial Energy Gamma", 2000, 0., 10000 * MeV); //9

        analysisManager->CreateH1("XVertex0", "Initial X Vertex", 1000, -50. * cm, 50. * cm); //10
        analysisManager->CreateH1("YVertex0", "Initial Y Vertex", 1000, -50. * cm, 50. * cm); //11
        analysisManager->CreateH1("ZVertex0", "Initial Z Vertex", 1000, -50. * cm, 50. * cm); //12

        analysisManager->CreateH1("ThetaVertex0", "Initial Theta Vertex", 1800, 0, 180); //13
        analysisManager->CreateH1("PhiVertex0", "Initial Phi Vertex", 3600, -180, 180); //14


        // Creating ntuple
        //

        G4int fCurrentNtupleId = 0;

        if (fEventAction) {

             analysisManager->SetFirstNtupleId(0);
            fCurrentNtupleId = analysisManager->CreateNtuple("TO", "Edep and TrackL");


            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"EventNumber"); // column Id = 0

            analysisManager->CreateNtupleDColumn(fCurrentNtupleId,"Eplastic_fat_nsys1_sum"); // column Id = 1
            analysisManager->CreateNtupleDColumn(fCurrentNtupleId,"Lplastic_fat_nsys1_sum"); // column Id = 2
            analysisManager->CreateNtupleDColumn(fCurrentNtupleId,"Eplastic_fat_nsys2_sum"); // column Id = 3
            analysisManager->CreateNtupleDColumn(fCurrentNtupleId,"Lplastic_fat_nsys2_sum"); // column Id = 4


            //  analysisManager->CreateNtupleDColumn(fCurrentNtupleId,"XVertex"); // column Id = 5  // Позиция вершины (точка генерации)
            //  analysisManager->CreateNtupleDColumn(fCurrentNtupleId,"YVertex"); // column Id = 6
            //  analysisManager->CreateNtupleDColumn(fCurrentNtupleId,"ZVertex"); // column Id = 7

            //  analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"IndexVertex"); // column Id = 8 // Индекс частицы по PDG
            //  analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"MassVertex"); // column Id = 9 // Индекс частицы по PDG
            //   analysisManager->CreateNtupleDColumn(fCurrentNtupleId,"EnergyVertex"); // column Id = 10 // Кинетическая энергия
            //  analysisManager->CreateNtupleDColumn(fCurrentNtupleId,"MomentumVertex"); // column Id = 11 // Импульс
            //   analysisManager->CreateNtupleDColumn(fCurrentNtupleId,"ThetaVertex"); // column Id = 12 // theta
            //  analysisManager->CreateNtupleDColumn(fCurrentNtupleId,"PhiVertex"); // column Id = 13 // phi

            analysisManager->CreateNtupleDColumn(fCurrentNtupleId,"Pzz"); // column Id = 5 // Pzz
            analysisManager->CreateNtupleDColumn(fCurrentNtupleId,"Egamma"); // column Id = 6 // Energy gamma
            analysisManager->CreateNtupleDColumn(fCurrentNtupleId,"NumReact"); // column Id = 7 // Number reaction
            analysisManager->CreateNtupleDColumn(fCurrentNtupleId,"NumEvnt"); // column Id = 8 // Number event

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"XVertex", fEventAction->GetVertexXYZ(1));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"YVertex", fEventAction->GetVertexXYZ(2));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ZVertex", fEventAction->GetVertexXYZ(3));
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"IndexVertex", fEventAction->GetVertexIndex());
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"MassVertex", fEventAction->GetVertexMass());
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"EnergyVertex", fEventAction->GetVertexEnergy());
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"MomentumVertex", fEventAction->GetVertexMomentum());
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ThetaVertex", fEventAction->GetVertexTheta());
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"PhiVertex", fEventAction->GetVertexPhi());

     //   TO->Branch("nev",&nev,"nev/I");
     //   TO->Branch("nr",&nr,"nr/I");
     //   TO->Branch("eg",&Eg,"eg/F");
     //   TO->Branch("vx",&vx,"vx/F");
    //    TO->Branch("vy",&vy,"vy/F");
    //    TO->Branch("vz",&vz,"vz/F");

            //   TO->Branch("nv",&nv,"nv/I");
            //    TO->Branch("iv",iv,"iv[nv]/F");
            //    TO->Branch("ev",ev,"ev[nv]/F");
            //    TO->Branch("tv",tv,"tv[nv]/F");
            //    TO->Branch("fv",fv,"fv[nv]/F");


            //  analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Eplastic_fat_nsys1_vector", fEventAction->GetPlasticEdep()); // column Id = 4

            //  analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Eplastic_fat_nsys2_vector"); // column Id = 5

#ifndef isGenLQ

#ifdef PF1_FAT

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"TrackIDplastic_fat_nsys1", fEventAction->GetTrackIDPlasticFat(1));

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Eplastic_fat_nsys1", fEventAction->GetPlasticFatEdep(
                                                     1)); // column Id = 4 // Энерговыделение


            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"LOplastic_fat_nsys1",
                                                 fEventAction->GetPlasticFatLO(1)); // column Id = 6 // Световыход


            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A1plastic_fat_nsys1", fEventAction->GetPlasticFatA1(
                                                     1)); // column Id = 8 // Амплутуды с противоположных торцов

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A2plastic_fat_nsys1",
                                                 fEventAction->GetPlasticFatA2(1)); // column Id = 10 //


            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T1plastic_fat_nsys1", fEventAction->GetPlasticFatT1(1));
            // column Id = 8 // Времена с противоположных торцов
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T2plastic_fat_nsys1", fEventAction->GetPlasticFatT2(1));
            // column Id = 10 //

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Lplastic_fat_nsys1", fEventAction->GetPlasticFatTrackLength(1));
            // column Id = 8 // Длина пробега


            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ToFplastic_fat_nsys1", fEventAction->GetPlasticFatToF(1));
            // column Id = 10 // TimeOfFlight


            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xplastic_fat_nsys1", fEventAction->GetPlasticFatPos(1,
                                                     1)); // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yplastic_fat_nsys1",
                                                 fEventAction->GetPlasticFatPos(1, 2)); // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zplastic_fat_nsys1",
                                                 fEventAction->GetPlasticFatPos(1, 3)); // column Id = 4 //


            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xglob_plastic_fat_nsys1", fEventAction->GetPlasticFatGlobalPos(1,
                                                     1)); // column Id = 4 // Глобальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yglob_plastic_fat_nsys1",
                                                 fEventAction->GetPlasticFatGlobalPos(1, 2)); // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zglob_plastic_fat_nsys1",
                                                 fEventAction->GetPlasticFatGlobalPos(1, 3)); // column Id = 4 //


            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Thetaplastic_fat_nsys1",
                                                 fEventAction->GetPlasticFatAngle(1, 1)); // column Id = 4 // Углы
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Phiplastic_fat_nsys1",
                                                 fEventAction->GetPlasticFatAngle(1, 2)); // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ThetaGlobplastic_fat_nsys1",
                                                 fEventAction->GetPlasticFatAngleGlob(1, 1)); // column Id = 4 // Углы
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"PhiGlobplastic_fat_nsys1",
                                                 fEventAction->GetPlasticFatAngleGlob(1, 2)); // column Id = 4 //

#endif

#ifdef PF2_FAT

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"TrackIDplastic_fat_nsys2", fEventAction->GetTrackIDPlasticFat(2));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"LOplastic_fat_nsys2",
                                                 fEventAction->GetPlasticFatLO(2)); // column Id = 7
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Eplastic_fat_nsys2",
                                                 fEventAction->GetPlasticFatEdep(2)); // column Id = 5

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Lplastic_fat_nsys2", fEventAction->GetPlasticFatTrackLength(2));
            // column Id = 9
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A1plastic_fat_nsys2",
                                                 fEventAction->GetPlasticFatA1(2)); // column Id = 9
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A2plastic_fat_nsys2",
                                                 fEventAction->GetPlasticFatA2(2)); // column Id = 11
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T1plastic_fat_nsys2", fEventAction->GetPlasticFatT1(2));
            // column Id = 9
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T2plastic_fat_nsys2", fEventAction->GetPlasticFatT2(2));
            // column Id = 11


            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ToFplastic_fat_nsys2",
                                                 fEventAction->GetPlasticFatToF(2)); // column Id = 11

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xplastic_fat_nsys2", fEventAction->GetPlasticFatPos(2,
                                                     1)); // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yplastic_fat_nsys2",
                                                 fEventAction->GetPlasticFatPos(2, 2)); // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zplastic_fat_nsys2",
                                                 fEventAction->GetPlasticFatPos(2, 3)); // column Id = 4 //

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xglob_plastic_fat_nsys2", fEventAction->GetPlasticFatGlobalPos(2,
                                                     1)); // column Id = 4 // Глобальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yglob_plastic_fat_nsys2",
                                                 fEventAction->GetPlasticFatGlobalPos(2, 2)); // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zglob_plastic_fat_nsys2",
                                                 fEventAction->GetPlasticFatGlobalPos(2, 3)); // column Id = 4 //


            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Thetaplastic_fat_nsys2", fEventAction->GetPlasticFatAngle(2,
                                                     1)); // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Phiplastic_fat_nsys2",
                                                 fEventAction->GetPlasticFatAngle(2, 2)); // column Id = 4 //

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ThetaGlobplastic_fat_nsys2", fEventAction->GetPlasticFatAngleGlob(2,
                                                     1)); // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"PhiGlobplastic_fat_nsys2",
                                                 fEventAction->GetPlasticFatAngleGlob(2, 2)); // column Id = 4 //

#endif

            /////////////////////// Для Адронного Калориметра

#ifdef HADCAL1
            //nsys1
            // X-bars
            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"TrackID_HCX_nsys1", fEventAction->Get_HCX_TrackID(1));
            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"N_HCX_nsys1",
                                                 fEventAction->Get_HCX_N(1)); // column Id = 4 // Число срабатываний

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"NSum_HCX_nsys1",
                                                 fEventAction->Get_HCX_NSum(1)); // column Id = 4 // Число срабатываний

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"NSum2_HCX_nsys1",
                                                 fEventAction->Get_HCX_NSum2(1)); // column Id = 4 // Число срабатываний

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"AL_HCX_nsys1",
                                                 fEventAction->Get_HCX_AL(1)); // column Id = 4 // Номер слоя


            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"E_HCX_nsys1",
                                                 fEventAction->Get_HCX_Edep(1)); // column Id = 4 // Энерговыделение
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"LO_HCX_nsys1",
                                                 fEventAction->Get_HCX_LO(1)); // column Id = 6 // Световыход

        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A1_HCX_nsys1", fEventAction->Get_HCX_A1(
                1)); // column Id = 8 // Амплутуды с противоположных торцов
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A2_HCX_nsys1", fEventAction->Get_HCX_A2(1)); // column Id = 10 //

        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T1_HCX_nsys1", fEventAction->Get_HCX_T1(
                1)); // column Id = 8 // Времена с противоположных торцов
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T2_HCX_nsys1", fEventAction->Get_HCX_T2(1)); // column Id = 10 //

        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"L_HCX_nsys1",
                                             fEventAction->Get_HCX_TrackLength(1)); // column Id = 8 // Длина пробега

        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ToF_HCX_nsys1",
                                             fEventAction->Get_HCX_ToF(1)); // column Id = 10 // Световыход

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ESum_HCX_nsys1",
                                                 fEventAction->Get_HCX_EdepSum(1)); // column Id = 4 // Энерговыделение
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"LOSum_HCX_nsys1",
                                                 fEventAction->Get_HCX_LOSum(1)); // column Id = 6 // Световыход

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A1Sum_HCX_nsys1", fEventAction->Get_HCX_A1Sum(
                                                     1)); // column Id = 8 // Амплутуды с противоположных торцов
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A2Sum_HCX_nsys1", fEventAction->Get_HCX_A2Sum(1));
            // column Id = 10 //

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T1Sum_HCX_nsys1", fEventAction->Get_HCX_T1Sum(
                                                     1)); // column Id = 8 // Времена с противоположных торцов
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T2Sum_HCX_nsys1", fEventAction->Get_HCX_T2Sum(1));
            // column Id = 10 //

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"LSum_HCX_nsys1",
                                                 fEventAction->Get_HCX_TrackLengthSum(1));
            // column Id = 8 // Длина пробега

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ToFSum_HCX_nsys1",
                                                 fEventAction->Get_HCX_ToFSum(1)); // column Id = 10 // Световыход


            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"X_HCX_nsys1", fEventAction->Get_HCX_Pos(1, 1));
            // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Y_HCX_nsys1", fEventAction->Get_HCX_Pos(1, 2)); // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Z_HCX_nsys1", fEventAction->Get_HCX_Pos(1, 3)); // column Id = 4 //

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xglob_HCX_nsys1", fEventAction->Get_HCX_GlobalPos(1, 1));
            // column Id = 4 // Глобальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yglob_HCX_nsys1", fEventAction->Get_HCX_GlobalPos(1, 2));
            // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zglob_HCX_nsys1", fEventAction->Get_HCX_GlobalPos(1, 3));
            // column Id = 4 //


            //Z-bars
            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"TrackID_HCZ_nsys1", fEventAction->Get_HCZ_TrackID(1));

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"N_HCZ_nsys1",
                                                 fEventAction->Get_HCZ_N(1)); // column Id = 4 // Число срабатываний

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"NSum_HCZ_nsys1",
                                                 fEventAction->Get_HCZ_NSum(1)); // column Id = 4 // Число срабатываний

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"NSum2_HCZ_nsys1",
                                                 fEventAction->Get_HCZ_NSum2(1)); // column Id = 4 // Число срабатываний

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"AL_HCZ_nsys1",
                                                 fEventAction->Get_HCZ_AL(1)); // column Id = 4 // Номер слоя

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"E_HCZ_nsys1",
                                                 fEventAction->Get_HCZ_Edep(1)); // column Id = 4 // Энерговыделение
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"LO_HCZ_nsys1",
                                                 fEventAction->Get_HCZ_LO(1)); // column Id = 6 // Световыход

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A1_HCZ_nsys1", fEventAction->Get_HCZ_A1(
                                                     1)); // column Id = 8 // Амплутуды с противоположных торцов
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A2_HCZ_nsys1", fEventAction->Get_HCZ_A2(1)); // column Id = 10 //

        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T1_HCZ_nsys1", fEventAction->Get_HCZ_T1(
                1)); // column Id = 8 // Времена с противоположных торцов
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T2_HCZ_nsys1", fEventAction->Get_HCZ_T2(1)); // column Id = 10 //

        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"L_HCZ_nsys1",
                                             fEventAction->Get_HCZ_TrackLength(1)); // column Id = 8 // Длина пробега

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ToF_HCZ_nsys1",
                                                 fEventAction->Get_HCZ_ToF(1)); // column Id = 10 // Световыход


            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ESum_HCZ_nsys1",
                                                 fEventAction->Get_HCZ_EdepSum(1)); // column Id = 4 // Энерговыделение
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"LOSum_HCZ_nsys1",
                                                 fEventAction->Get_HCZ_LOSum(1)); // column Id = 6 // Световыход

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A1Sum_HCZ_nsys1", fEventAction->Get_HCZ_A1Sum(
                                                     1)); // column Id = 8 // Амплутуды с противоположных торцов
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A2Sum_HCZ_nsys1", fEventAction->Get_HCZ_A2Sum(1));
            // column Id = 10 //

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T1Sum_HCZ_nsys1", fEventAction->Get_HCZ_T1Sum(
                                                     1)); // column Id = 8 // Времена с противоположных торцов
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T2Sum_HCZ_nsys1", fEventAction->Get_HCZ_T2Sum(1));
            // column Id = 10 //

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"LSum_HCZ_nsys1",
                                                 fEventAction->Get_HCZ_TrackLengthSum(1));
            // column Id = 8 // Длина пробега

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ToFSum_HCZ_nsys1",
                                                 fEventAction->Get_HCZ_ToFSum(1)); // column Id = 10 // Световыход


            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"X_HCZ_nsys1", fEventAction->Get_HCZ_Pos(1, 1));
            // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Y_HCZ_nsys1", fEventAction->Get_HCZ_Pos(1, 2)); // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Z_HCZ_nsys1", fEventAction->Get_HCZ_Pos(1, 3)); // column Id = 4 //

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xglob_HCZ_nsys1", fEventAction->Get_HCZ_GlobalPos(1, 1));
            // column Id = 4 // Глобальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yglob_HCZ_nsys1", fEventAction->Get_HCZ_GlobalPos(1, 2));
            // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zglob_HCZ_nsys1", fEventAction->Get_HCZ_GlobalPos(1, 3));
            // column Id = 4 //


#endif

#ifdef HADCAL2

            //nsys2

            // X-bars
            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"TrackID_HCX_nsys2", fEventAction->Get_HCX_TrackID(2));
            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"N_HCX_nsys2",
                                                 fEventAction->Get_HCX_N(2)); // column Id = 4 // Число срабатываний

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"NSum_HCX_nsys2",
                                                 fEventAction->Get_HCX_NSum(2)); // column Id = 4 // Число срабатываний

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"NSum2_HCX_nsys2",
                                                 fEventAction->Get_HCX_NSum2(2)); // column Id = 4 // Число срабатываний

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"AL_HCX_nsys2",
                                                 fEventAction->Get_HCX_AL(2)); // column Id = 4 // Номер слоя

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"E_HCX_nsys2",
                                                 fEventAction->Get_HCX_Edep(2)); // column Id = 4 // Энерговыделение
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"LO_HCX_nsys2",
                                             fEventAction->Get_HCX_LO(2)); // column Id = 6 // Световыход

        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A1_HCX_nsys2", fEventAction->Get_HCX_A1(
                2)); // column Id = 8 // Амплутуды с противоположных торцов
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A2_HCX_nsys2", fEventAction->Get_HCX_A2(2)); // column Id = 10 //

        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T1_HCX_nsys2", fEventAction->Get_HCX_T1(
                2)); // column Id = 8 // Времена с противоположных торцов
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T2_HCX_nsys2", fEventAction->Get_HCX_T2(2)); // column Id = 10 //

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"L_HCX_nsys2",
                                                 fEventAction->Get_HCX_TrackLength(2));
            // column Id = 8 // Длина пробега

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ToF_HCX_nsys2",
                                                 fEventAction->Get_HCX_ToF(2)); // column Id = 10 // Световыход


            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ESum_HCX_nsys2",
                                                 fEventAction->Get_HCX_EdepSum(2)); // column Id = 4 // Энерговыделение
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"LOSum_HCX_nsys2",
                                                 fEventAction->Get_HCX_LOSum(2)); // column Id = 6 // Световыход

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A1Sum_HCX_nsys2", fEventAction->Get_HCX_A1Sum(
                                                     2)); // column Id = 8 // Амплутуды с противоположных торцов
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A2Sum_HCX_nsys2", fEventAction->Get_HCX_A2Sum(2));
            // column Id = 10 //

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T1Sum_HCX_nsys2", fEventAction->Get_HCX_T1Sum(
                                                     2)); // column Id = 8 // Времена с противоположных торцов
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T2Sum_HCX_nsys2", fEventAction->Get_HCX_T2Sum(2));
            // column Id = 10 //

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"LSum_HCX_nsys2",
                                                 fEventAction->Get_HCX_TrackLengthSum(2));
            // column Id = 8 // Длина пробега

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ToFSum_HCX_nsys2",
                                                 fEventAction->Get_HCX_ToFSum(2)); // column Id = 10 // Световыход


            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"X_HCX_nsys2", fEventAction->Get_HCX_Pos(2, 1));
            // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Y_HCX_nsys2", fEventAction->Get_HCX_Pos(2, 2)); // column Id = 4 //
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Z_HCX_nsys2", fEventAction->Get_HCX_Pos(2, 3)); // column Id = 4 //


            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xglob_HCX_nsys2", fEventAction->Get_HCX_GlobalPos(2, 1));
            // column Id = 4 // Глобальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yglob_HCX_nsys2", fEventAction->Get_HCX_GlobalPos(2, 2));
            // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zglob_HCX_nsys2", fEventAction->Get_HCX_GlobalPos(2, 3));
            // column Id = 4 //



//Z-bars

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"TrackID_HCZ_nsys2", fEventAction->Get_HCZ_TrackID(2));
            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"N_HCZ_nsys2",
                                                 fEventAction->Get_HCZ_N(2)); // column Id = 4 // Число срабатываний
            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"NSum_HCZ_nsys2",
                                                 fEventAction->Get_HCZ_NSum(2)); // column Id = 4 // Число срабатываний

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"NSum2_HCZ_nsys2",
                                                 fEventAction->Get_HCZ_NSum2(2)); // column Id = 4 // Число срабатываний

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"AL_HCZ_nsys2",
                                                 fEventAction->Get_HCZ_AL(2)); // column Id = 4 // Номер слоя

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"E_HCZ_nsys2",
                                                 fEventAction->Get_HCZ_Edep(2)); // column Id = 4 // Энерговыделение
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"LO_HCZ_nsys2",
                                                 fEventAction->Get_HCZ_LO(2)); // column Id = 6 // Световыход

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A1_HCZ_nsys2", fEventAction->Get_HCZ_A1(
                                                     2)); // column Id = 8 // Амплутуды с противоположных торцов
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A2_HCZ_nsys2", fEventAction->Get_HCZ_A2(2)); // column Id = 10 //

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T1_HCZ_nsys2", fEventAction->Get_HCZ_T1(
                                                     2)); // column Id = 8 // Времена с противоположных торцов
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T2_HCZ_nsys2", fEventAction->Get_HCZ_T2(2)); // column Id = 10 //

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"L_HCZ_nsys2",
                                                 fEventAction->Get_HCZ_TrackLength(2));
            // column Id = 8 // Длина пробега

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ToF_HCZ_nsys2",
                                                 fEventAction->Get_HCZ_ToF(2)); // column Id = 10 // Световыход


            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ESum_HCZ_nsys2",
                                                 fEventAction->Get_HCZ_EdepSum(2)); // column Id = 4 // Энерговыделение
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"LOSum_HCZ_nsys2",
                                                 fEventAction->Get_HCZ_LOSum(2)); // column Id = 6 // Световыход

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A1Sum_HCZ_nsys2", fEventAction->Get_HCZ_A1Sum(
                                                     2)); // column Id = 8 // Амплутуды с противоположных торцов
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A2Sum_HCZ_nsys2", fEventAction->Get_HCZ_A2Sum(2));
            // column Id = 10 //

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T1Sum_HCZ_nsys2", fEventAction->Get_HCZ_T1Sum(
                                                     2)); // column Id = 8 // Времена с противоположных торцов
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T2Sum_HCZ_nsys2", fEventAction->Get_HCZ_T2Sum(2));
            // column Id = 10 //

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"LSum_HCZ_nsys2",
                                                 fEventAction->Get_HCZ_TrackLengthSum(2));
            // column Id = 8 // Длина пробега

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ToFSum_HCZ_nsys2",
                                                 fEventAction->Get_HCZ_ToFSum(2)); // column Id = 10 // Световыход


            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"X_HCZ_nsys2", fEventAction->Get_HCZ_Pos(2, 1)); // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Y_HCZ_nsys2", fEventAction->Get_HCZ_Pos(2, 2)); // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Z_HCZ_nsys2", fEventAction->Get_HCZ_Pos(2, 3)); // column Id = 4 //

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xglob_HCZ_nsys2", fEventAction->Get_HCZ_GlobalPos(2, 1));
            // column Id = 4 // Глобальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yglob_HCZ_nsys2", fEventAction->Get_HCZ_GlobalPos(2, 2));
            // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zglob_HCZ_nsys2", fEventAction->Get_HCZ_GlobalPos(2, 3)); // column Id = 4 //




#endif

#endif// notdef isGenLQ

            // ТОНКИЕ ПЛАСТИКИ

/////////////////////////////

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"TrackIDplastic_thin_nsys1", fEventAction->GetTrackIDPlasticThin(1));
            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"TrackIDplastic_thin_nsys2", fEventAction->GetTrackIDPlasticThin(2));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Eplastic_thin_nsys1",
                                                 fEventAction->GetPlasticThinEdep(1)); // column Id = 12
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Eplastic_thin_nsys2",
                                                 fEventAction->GetPlasticThinEdep(2)); // column Id = 13

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"LOplastic_thin_nsys1",
                                                 fEventAction->GetPlasticThinLO(1)); // column Id = 14
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"LOplastic_thin_nsys2",
                                                 fEventAction->GetPlasticThinLO(2)); // column Id = 15

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Lplastic_thin_nsys1",
                                                 fEventAction->GetPlasticThinTrackLength(1)); // column Id = 16
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Lplastic_thin_nsys2",
                                                 fEventAction->GetPlasticThinTrackLength(2)); // column Id = 17

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A1plastic_thin_nsys1", fEventAction->GetPlasticThinA1(
                    1)); // column Id = 8 // Амплутуды с противоположных торцов
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A1plastic_thin_nsys2", fEventAction->GetPlasticThinA1(2)); // column Id = 9

        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T1plastic_thin_nsys1", fEventAction->GetPlasticThinT1(1)); // column Id = 8 // Времена с противоположных торцов
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T1plastic_thin_nsys2", fEventAction->GetPlasticThinT1(2)); // column Id = 9

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ToFplastic_thin_nsys1",
                                                 fEventAction->GetPlasticThinToF(1)); // column Id = 18
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ToFplastic_thin_nsys2",
                                                 fEventAction->GetPlasticThinToF(2)); // column Id = 19

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xplastic_thin_nsys1", fEventAction->GetPlasticThinPos(1,
                                                                                                        1)); // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yplastic_thin_nsys1",
                                                 fEventAction->GetPlasticThinPos(1, 2)); // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zplastic_thin_nsys1",
                                                 fEventAction->GetPlasticThinPos(1, 3)); // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xplastic_thin_nsys2", fEventAction->GetPlasticThinPos(2,
                                                                                                        1)); // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yplastic_thin_nsys2",
                                                 fEventAction->GetPlasticThinPos(2, 2)); // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zplastic_thin_nsys2",
                                                 fEventAction->GetPlasticThinPos(2, 3)); // column Id = 4 //


            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xglob_plastic_thin_nsys1", fEventAction->GetPlasticThinGlobalPos(1,
                                                                                                                   1)); // column Id = 4 // Глобальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yglob_plastic_thin_nsys1",
                                                 fEventAction->GetPlasticThinGlobalPos(1, 2)); // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zglob_plastic_thin_nsys1",
                                                 fEventAction->GetPlasticThinGlobalPos(1, 3)); // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xglob_plastic_thin_nsys2", fEventAction->GetPlasticThinGlobalPos(2,
                                                                                                                   1)); // column Id = 4 // Глобальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yglob_plastic_thin_nsys2",
                                                 fEventAction->GetPlasticThinGlobalPos(2, 2)); // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zglob_plastic_thin_nsys2",
                                                 fEventAction->GetPlasticThinGlobalPos(2, 3)); // column Id = 4 //

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Thetaplastic_thin_nsys1",
                                                 fEventAction->GetPlasticThinAngle(1, 1)); // column Id = 4 // Углы
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Phiplastic_thin_nsys1",
                                                 fEventAction->GetPlasticThinAngle(1, 2)); // column Id = 4 //

                        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ThetaGlobplastic_thin_nsys1",
                                                             fEventAction->GetPlasticThinAngleGlob(1, 1));
            // column Id = 4 // Углы
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"PhiGlobplastic_thin_nsys1",
                                                 fEventAction->GetPlasticThinAngleGlob(1, 2)); // column Id = 4 //

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Thetaplastic_thin_nsys2", fEventAction->GetPlasticThinAngle(2,
                                                     1)); // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Phiplastic_thin_nsys2",
                                                 fEventAction->GetPlasticThinAngle(2, 2)); // column Id = 4 //

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ThetaGlobplastic_thin_nsys2", fEventAction->GetPlasticThinAngleGlob(
                                                     2,
                                                     1)); // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"PhiGlobplastic_thin_nsys2",
                                                 fEventAction->GetPlasticThinAngleGlob(2, 2)); // column Id = 4 //

            // Для пластиков LQ-поляриметра
            /////////////////////////////
            // здесь индексы системы перепутаны!!! 1<->2

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"TrackIDplastic_LQ_nsys1", fEventAction->GetTrackIDPlasticLQ(2));
            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"TrackIDplastic_LQ_nsys2", fEventAction->GetTrackIDPlasticLQ(1));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Eplastic_LQ_nsys1",
                                                 fEventAction->GetPlasticLQEdep(2)); // column Id = 12
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Eplastic_LQ_nsys2",
                                                 fEventAction->GetPlasticLQEdep(1)); // column Id = 13

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"LOplastic_LQ_nsys1",
                                                 fEventAction->GetPlasticLQLO(2)); // column Id = 14
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"LOplastic_LQ_nsys2",
                                                 fEventAction->GetPlasticLQLO(1)); // column Id = 15

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Lplastic_LQ_nsys1",
                                                 fEventAction->GetPlasticLQTrackLength(2)); // column Id = 16
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Lplastic_LQ_nsys2",
                                                 fEventAction->GetPlasticLQTrackLength(1)); // column Id = 17

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A1plastic_LQ_nsys1", fEventAction->GetPlasticLQA1(
                                                     2)); // column Id = 8 // Амплутуды с противоположных торцов
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A2plastic_LQ_nsys1", fEventAction->GetPlasticLQA2(2));
            // column Id = 9
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A1plastic_LQ_nsys2", fEventAction->GetPlasticLQA1(1));
            // column Id = 9
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"A2plastic_LQ_nsys2", fEventAction->GetPlasticLQA2(1));

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T1plastic_LQ_nsys1", fEventAction->GetPlasticLQT1(2));
            // column Id = 8 // Времена с противоположных торцов
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T2plastic_LQ_nsys1", fEventAction->GetPlasticLQT2(2));
            // column Id = 9
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T1plastic_LQ_nsys2", fEventAction->GetPlasticLQT1(1));
            // column Id = 9
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"T2plastic_LQ_nsys2", fEventAction->GetPlasticLQT2(1));
            // column Id = 9

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ToFplastic_LQ_nsys1",
                                                 fEventAction->GetPlasticLQToF(2)); // column Id = 18
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ToFplastic_LQ_nsys2",
                                                 fEventAction->GetPlasticLQToF(1)); // column Id = 19

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xplastic_LQ_nsys1", fEventAction->GetPlasticLQPos(2,
                                                     1)); // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yplastic_LQ_nsys1",
                                                 fEventAction->GetPlasticLQPos(2, 2)); // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zplastic_LQ_nsys1",
                                                 fEventAction->GetPlasticLQPos(2, 3)); // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xplastic_LQ_nsys2", fEventAction->GetPlasticLQPos(1,
                                                     1)); // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yplastic_LQ_nsys2",
                                                 fEventAction->GetPlasticLQPos(1, 2)); // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zplastic_LQ_nsys2",
                                                 fEventAction->GetPlasticLQPos(1, 3)); // column Id = 4 //


            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xglob_plastic_LQ_nsys1", fEventAction->GetPlasticLQGlobalPos(2,
                                                     1)); // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yglob_plastic_LQ_nsys1",
                                                 fEventAction->GetPlasticLQGlobalPos(2, 2)); // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zglob_plastic_LQ_nsys1",
                                                 fEventAction->GetPlasticLQGlobalPos(2, 3)); // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xglob_plastic_LQ_nsys2", fEventAction->GetPlasticLQGlobalPos(1,
                                                     1)); // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yglob_plastic_LQ_nsys2",
                                                 fEventAction->GetPlasticLQGlobalPos(1, 2)); // column Id = 4 //
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zglob_plastic_LQ_nsys2",
                                                 fEventAction->GetPlasticLQGlobalPos(1, 3)); // column Id = 4 //


            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Thetaplastic_LQ_nsys1",
                                                 fEventAction->GetPlasticLQAngle(2, 1)); // column Id = 4 // Углы
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Phiplastic_LQ_nsys1",
                                                 fEventAction->GetPlasticLQAngle(2, 2)); // column Id = 4 //

                    analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ThetaGlobplastic_LQ_nsys1",
                                                 fEventAction->GetPlasticLQAngleGlob(2, 1)); // column Id = 4 // Углы
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"PhiGlobplastic_LQ_nsys1",
                                                 fEventAction->GetPlasticLQAngleGlob(2, 2)); // column Id = 4 //

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Thetaplastic_LQ_nsys2", fEventAction->GetPlasticLQAngle(1,
                                                     1)); // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Phiplastic_LQ_nsys2",
                                                 fEventAction->GetPlasticLQAngle(1, 2)); // column Id = 4 //

                        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"ThetaGlobplastic_LQ_nsys2",
                                                             fEventAction->GetPlasticLQAngleGlob(1,
                                                                 1)); // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"PhiGlobplastic_LQ_nsys2",
                                                 fEventAction->GetPlasticLQAngleGlob(1, 2)); // column Id = 4 //


            /// Для трековых камер
#if defined(DCARM1)
            // Камера A

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"TrackIDwa_nsys1", fEventAction->Get_W_TrackID(1, 1));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Ewa_nsys1", fEventAction->Get_W_Edep(1, 1));
            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"Nwa_nsys1", fEventAction->Get_W_N(1, 1));
            // column Id = 4 // Число срабатываний
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xwa_nsys1", fEventAction->Get_W_Pos(1, 1, 1));
            // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Ywa_nsys1", fEventAction->Get_W_Pos(1, 2, 1));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zwa_nsys1", fEventAction->Get_W_Pos(1, 3, 1));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Rwa_nsys1", fEventAction->Get_W_Pos(1, 4, 1));


            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xglobwa_nsys1", fEventAction->Get_W_GlobalPos(1, 1, 1));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yglobwa_nsys1", fEventAction->Get_W_GlobalPos(1, 2, 1));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zglobwa_nsys1", fEventAction->Get_W_GlobalPos(1, 3, 1));


            // Камера B
            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"TrackIDwb_nsys1", fEventAction->Get_W_TrackID(1, 2));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Ewb_nsys1", fEventAction->Get_W_Edep(1, 2));
            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"Nwb_nsys1", fEventAction->Get_W_N(1, 2));
            // column Id = 4 // Число срабатываний
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xwb_nsys1", fEventAction->Get_W_Pos(1, 1, 2));
            // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Ywb_nsys1", fEventAction->Get_W_Pos(1, 2, 2));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zwb_nsys1", fEventAction->Get_W_Pos(1, 3, 2));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Rwb_nsys1", fEventAction->Get_W_Pos(1, 4, 2));


            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xglobwb_nsys1", fEventAction->Get_W_GlobalPos(1, 1, 2));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yglobwb_nsys1", fEventAction->Get_W_GlobalPos(1, 2, 2));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zglobwb_nsys1", fEventAction->Get_W_GlobalPos(1, 3, 2));


            // Камера C

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"TrackIDwc_nsys1", fEventAction->Get_W_TrackID(1, 3));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Ewc_nsys1", fEventAction->Get_W_Edep(1, 3));
            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"Nwc_nsys1", fEventAction->Get_W_N(1, 3));
            // column Id = 4 // Число срабатываний
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xwc_nsys1", fEventAction->Get_W_Pos(1, 1, 3));
            // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Ywc_nsys1", fEventAction->Get_W_Pos(1, 2, 3));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zwc_nsys1", fEventAction->Get_W_Pos(1, 3, 3));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Rwc_nsys1", fEventAction->Get_W_Pos(1, 4, 3));

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xglobwc_nsys1", fEventAction->Get_W_GlobalPos(1, 1, 3));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yglobwc_nsys1", fEventAction->Get_W_GlobalPos(1, 2, 3));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zglobwc_nsys1", fEventAction->Get_W_GlobalPos(1, 3, 3));



#endif


#if defined(DCARM2)

            // Камера A

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"TrackIDwa_nsys2", fEventAction->Get_W_TrackID(2, 1));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Ewa_nsys2", fEventAction->Get_W_Edep(2, 1));
            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"Nwa_nsys2", fEventAction->Get_W_N(2));
            // column Id = 4 // Число срабатываний
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xwa_nsys2", fEventAction->Get_W_Pos(2, 1, 1));
            // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Ywa_nsys2", fEventAction->Get_W_Pos(2, 2, 1));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zwa_nsys2", fEventAction->Get_W_Pos(2, 3, 1));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Rwa_nsys2", fEventAction->Get_W_Pos(2, 4, 1));

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xglobwa_nsys2", fEventAction->Get_W_GlobalPos(2, 1, 1));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yglobwa_nsys2", fEventAction->Get_W_GlobalPos(2, 2, 1));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zglobwa_nsys2", fEventAction->Get_W_GlobalPos(2, 3, 1));


            // Камера B

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"TrackIDwb_nsys2", fEventAction->Get_W_TrackID(2, 2));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Ewb_nsys2", fEventAction->Get_W_Edep(2, 2));
            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"Nwb_nsys2", fEventAction->Get_W_N(2, 2));
            // column Id = 4 // Число срабатываний
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xwb_nsys2", fEventAction->Get_W_Pos(2, 1, 2));
            // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Ywb_nsys2", fEventAction->Get_W_Pos(2, 2, 2));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zwb_nsys2", fEventAction->Get_W_Pos(2, 3, 2));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Rwb_nsys2", fEventAction->Get_W_Pos(2, 4, 2));

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xglobwb_nsys2", fEventAction->Get_W_GlobalPos(2, 1, 2));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yglobwb_nsys2", fEventAction->Get_W_GlobalPos(2, 2, 2));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zglobwb_nsys2", fEventAction->Get_W_GlobalPos(2, 3, 2));


            // Камера C

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"TrackIDwc_nsys2", fEventAction->Get_W_TrackID(2, 3));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Ewc_nsys2", fEventAction->Get_W_Edep(2, 3));
            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"Nwc_nsys2", fEventAction->Get_W_N(2, 3));
            // column Id = 4 // Число срабатываний
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xwc_nsys2", fEventAction->Get_W_Pos(2, 1, 3));
            // column Id = 4 // Локальная точка
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Ywc_nsys2", fEventAction->Get_W_Pos(2, 2, 3));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zwc_nsys2", fEventAction->Get_W_Pos(2, 3, 3));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Rwc_nsys2", fEventAction->Get_W_Pos(2, 4, 3));

            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xglobwc_nsys2", fEventAction->Get_W_GlobalPos(2, 1, 3));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yglobwc_nsys2", fEventAction->Get_W_GlobalPos(2, 2, 3));
            analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zglobwc_nsys2", fEventAction->Get_W_GlobalPos(2, 3, 3));


#endif

            /// Для вершинных камер

#if defined(VCARM1) && defined(RUN21)

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"TrackIDvc_nsys1", fEventAction->Get_VC_TrackID(1));
            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"Evb_nsys1", fEventAction->Get_VC_Edep(1));
        analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"Nvc_nsys1",
                                             fEventAction->Get_VC_N(1)); // column Id = 4 // Число срабатываний
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xvc_nsys1",
                                             fEventAction->Get_VC_Pos(1, 1)); // column Id = 4 // Локальная точка
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yvc_nsys1", fEventAction->Get_VC_Pos(1, 2));
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zvc_nsys1", fEventAction->Get_VC_Pos(1, 3));
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Rvc_nsys1", fEventAction->Get_VC_Pos(1, 4));
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Massvc_nsys1", fEventAction->Get_VC_Mass(1));
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Energyvc_nsys1", fEventAction->Get_VC_KineticEnergy(1));
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Thetavc_nsys1", fEventAction->Get_VC_Theta(1));
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Phivc_nsys1", fEventAction->Get_VC_Phi(1));

        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xglobvc_nsys1", fEventAction->Get_VC_GlobalPos(1, 1));
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yglobvc_nsys1", fEventAction->Get_VC_GlobalPos(1, 2));
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zglobvc_nsys1", fEventAction->Get_VC_GlobalPos(1, 3));




#endif
#if defined(VCARM2) && defined(RUN21)

            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"TrackIDvc_nsys2", fEventAction->Get_VC_TrackID(2));
            analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"Evb_nsys2", fEventAction->Get_VC_Edep(2));
        analysisManager->CreateNtupleIColumn(fCurrentNtupleId,"Nvc_nsys2",
                                             fEventAction->Get_VC_N(2)); // column Id = 4 // Число срабатываний
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xvc_nsys2",
                                             fEventAction->Get_VC_Pos(2, 1)); // column Id = 4 // Локальная точка
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yvc_nsys2", fEventAction->Get_VC_Pos(2, 2));
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zvc_nsys2", fEventAction->Get_VC_Pos(2, 3));
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Rvc_nsys2", fEventAction->Get_VC_Pos(2, 4));
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Massvc_nsys2", fEventAction->Get_VC_Mass(2));
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Energyvc_nsys2", fEventAction->Get_VC_KineticEnergy(2));
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Thetavc_nsys2", fEventAction->Get_VC_Theta(2));
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Phivc_nsys2", fEventAction->Get_VC_Phi(2));

        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Xglobvc_nsys2", fEventAction->Get_VC_GlobalPos(2, 1));
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Yglobvc_nsys2", fEventAction->Get_VC_GlobalPos(2, 2));
        analysisManager->CreateNtupleFColumn(fCurrentNtupleId,"Zglobvc_nsys2", fEventAction->Get_VC_GlobalPos(2, 3));



#endif

        /////////////////

        analysisManager->FinishNtuple(fCurrentNtupleId);
    }

// Set ntuple output file
    analysisManager->SetNtupleFileName(fCurrentNtupleId, "Cosmicntuple");





    }
}