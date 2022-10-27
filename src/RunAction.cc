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
#include "Genbos.hh"


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

        analysisManager->CreateNtupleIColumn("EventNumber"); // column Id = 0

        analysisManager->CreateNtupleDColumn("Eplastic_fat_nsys1_sum"); // column Id = 1
        analysisManager->CreateNtupleDColumn("Lplastic_fat_nsys1_sum"); // column Id = 2
        analysisManager->CreateNtupleDColumn("Eplastic_fat_nsys2_sum"); // column Id = 3
        analysisManager->CreateNtupleDColumn("Lplastic_fat_nsys2_sum"); // column Id = 4

      //  analysisManager->CreateNtupleDColumn("XVertex"); // column Id = 5  // Позиция вершины (точка генерации)
      //  analysisManager->CreateNtupleDColumn("YVertex"); // column Id = 6
      //  analysisManager->CreateNtupleDColumn("ZVertex"); // column Id = 7

      //  analysisManager->CreateNtupleIColumn("IndexVertex"); // column Id = 8 // Индекс частицы по PDG
      //  analysisManager->CreateNtupleIColumn("MassVertex"); // column Id = 9 // Индекс частицы по PDG
     //   analysisManager->CreateNtupleDColumn("EnergyVertex"); // column Id = 10 // Кинетическая энергия
      //  analysisManager->CreateNtupleDColumn("MomentumVertex"); // column Id = 11 // Импульс
     //   analysisManager->CreateNtupleDColumn("ThetaVertex"); // column Id = 12 // theta
      //  analysisManager->CreateNtupleDColumn("PhiVertex"); // column Id = 13 // phi

        analysisManager->CreateNtupleDColumn("Pzz"); // column Id = 5 // Pzz
        analysisManager->CreateNtupleDColumn("Egamma"); // column Id = 6 // Energy gamma
        analysisManager->CreateNtupleDColumn("NumReact"); // column Id = 7 // Number reaction
        analysisManager->CreateNtupleDColumn("NumEvnt"); // column Id = 8 // Number event

        analysisManager->CreateNtupleFColumn("XVertex", fEventAction->GetVertexXYZ(1));
        analysisManager->CreateNtupleFColumn("YVertex", fEventAction->GetVertexXYZ(2));
        analysisManager->CreateNtupleFColumn("ZVertex", fEventAction->GetVertexXYZ(3));
        analysisManager->CreateNtupleFColumn("IndexVertex", fEventAction->GetVertexIndex());
        analysisManager->CreateNtupleFColumn("MassVertex", fEventAction->GetVertexMass());
        analysisManager->CreateNtupleFColumn("EnergyVertex", fEventAction->GetVertexEnergy());
        analysisManager->CreateNtupleFColumn("MomentumVertex", fEventAction->GetVertexMomentum());
        analysisManager->CreateNtupleFColumn("ThetaVertex", fEventAction->GetVertexTheta());
        analysisManager->CreateNtupleFColumn("PhiVertex", fEventAction->GetVertexPhi());

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


      //  analysisManager->CreateNtupleFColumn("Eplastic_fat_nsys1_vector", fEventAction->GetPlasticEdep()); // column Id = 4
      //  analysisManager->CreateNtupleFColumn("Eplastic_fat_nsys2_vector"); // column Id = 5

        analysisManager->CreateNtupleFColumn("Eplastic_fat_nsys1", fEventAction->GetPlasticFatEdep(1)); // column Id = 4 // Энерговыделение
        analysisManager->CreateNtupleFColumn("Eplastic_fat_nsys2", fEventAction->GetPlasticFatEdep(2)); // column Id = 5

        analysisManager->CreateNtupleFColumn("LOplastic_fat_nsys1", fEventAction->GetPlasticFatLO(1)); // column Id = 6 // Световыход
        analysisManager->CreateNtupleFColumn("LOplastic_fat_nsys2", fEventAction->GetPlasticFatLO(2)); // column Id = 7

        analysisManager->CreateNtupleFColumn("A1plastic_fat_nsys1", fEventAction->GetPlasticFatA1(1)); // column Id = 8 // Амплутуды с противоположных торцов
        analysisManager->CreateNtupleFColumn("A1plastic_fat_nsys2", fEventAction->GetPlasticFatA1(2)); // column Id = 9
        analysisManager->CreateNtupleFColumn("A2plastic_fat_nsys1", fEventAction->GetPlasticFatA2(1)); // column Id = 10 //
        analysisManager->CreateNtupleFColumn("A2plastic_fat_nsys2", fEventAction->GetPlasticFatA2(2)); // column Id = 11

        analysisManager->CreateNtupleFColumn("T1plastic_fat_nsys1", fEventAction->GetPlasticFatT1(1)); // column Id = 8 // Времена с противоположных торцов
        analysisManager->CreateNtupleFColumn("T1plastic_fat_nsys2", fEventAction->GetPlasticFatT1(2)); // column Id = 9
        analysisManager->CreateNtupleFColumn("T2plastic_fat_nsys1", fEventAction->GetPlasticFatT2(1)); // column Id = 10 //
        analysisManager->CreateNtupleFColumn("T2plastic_fat_nsys2", fEventAction->GetPlasticFatT2(2)); // column Id = 11

        analysisManager->CreateNtupleFColumn("Lplastic_fat_nsys1", fEventAction->GetPlasticFatTrackLength(1)); // column Id = 8 // Длина пробега
        analysisManager->CreateNtupleFColumn("Lplastic_fat_nsys2", fEventAction->GetPlasticFatTrackLength(2)); // column Id = 9

        analysisManager->CreateNtupleFColumn("ToFplastic_fat_nsys1", fEventAction->GetPlasticFatToF(1)); // column Id = 10 // TimeOfFlight
        analysisManager->CreateNtupleFColumn("ToFplastic_fat_nsys2", fEventAction->GetPlasticFatToF(2)); // column Id = 11

        analysisManager->CreateNtupleFColumn("Xplastic_fat_nsys1", fEventAction->GetPlasticFatPos(1, 1)); // column Id = 4 // Локальная точка
        analysisManager->CreateNtupleFColumn("Yplastic_fat_nsys1", fEventAction->GetPlasticFatPos(1, 2)); // column Id = 4 //
        analysisManager->CreateNtupleFColumn("Zplastic_fat_nsys1", fEventAction->GetPlasticFatPos(1, 3)); // column Id = 4 //
        analysisManager->CreateNtupleFColumn("Xplastic_fat_nsys2", fEventAction->GetPlasticFatPos(2, 1)); // column Id = 4 // Локальная точка
        analysisManager->CreateNtupleFColumn("Yplastic_fat_nsys2", fEventAction->GetPlasticFatPos(2, 2)); // column Id = 4 //
        analysisManager->CreateNtupleFColumn("Zplastic_fat_nsys2", fEventAction->GetPlasticFatPos(2, 3)); // column Id = 4 //


        analysisManager->CreateNtupleFColumn("Thetaplastic_fat_nsys1", fEventAction->GetPlasticFatAngle(1, 1)); // column Id = 4 // Углы
        analysisManager->CreateNtupleFColumn("Phiplastic_fat_nsys1", fEventAction->GetPlasticFatAngle(1, 2)); // column Id = 4 //

        analysisManager->CreateNtupleFColumn("Thetaplastic_fat_nsys2", fEventAction->GetPlasticFatAngle(2, 1)); // column Id = 4 // Локальная точка
        analysisManager->CreateNtupleFColumn("Phiplastic_fat_nsys2", fEventAction->GetPlasticFatAngle(2, 2)); // column Id = 4 //


    /////////////////////// Для Адронного Калориметра

    //nsys1
    // X-bars
        analysisManager->CreateNtupleIColumn("N_HCX_nsys1", fEventAction->Get_HCX_N(1)); // column Id = 4 // Число срабатываний
        analysisManager->CreateNtupleIColumn("AL_HCX_nsys1", fEventAction->Get_HCX_AL(1)); // column Id = 4 // Номер слоя

        analysisManager->CreateNtupleFColumn("E_HCX_nsys1", fEventAction->Get_HCX_Edep(1)); // column Id = 4 // Энерговыделение
        analysisManager->CreateNtupleFColumn("LO_HCX_nsys1", fEventAction->Get_HCX_LO(1)); // column Id = 6 // Световыход

        analysisManager->CreateNtupleFColumn("A1_HCX_nsys1", fEventAction->Get_HCX_A1(1)); // column Id = 8 // Амплутуды с противоположных торцов
        analysisManager->CreateNtupleFColumn("A2_HCX_nsys1", fEventAction->Get_HCX_A2(1)); // column Id = 10 //

        analysisManager->CreateNtupleFColumn("T1_HCX_nsys1", fEventAction->Get_HCX_T1(1)); // column Id = 8 // Времена с противоположных торцов
        analysisManager->CreateNtupleFColumn("T2_HCX_nsys1", fEventAction->Get_HCX_T2(1)); // column Id = 10 //

        analysisManager->CreateNtupleFColumn("L_HCX_nsys1", fEventAction->Get_HCX_TrackLength(1)); // column Id = 8 // Длина пробега

        analysisManager->CreateNtupleFColumn("ToF_HCX_nsys1", fEventAction->Get_HCX_ToF(1)); // column Id = 10 // Световыход


        analysisManager->CreateNtupleFColumn("X_HCX_nsys1", fEventAction->Get_HCX_Pos(1, 1)); // column Id = 4 // Локальная точка
        analysisManager->CreateNtupleFColumn("Y_HCX_nsys1", fEventAction->Get_HCX_Pos(1, 2)); // column Id = 4 //
        analysisManager->CreateNtupleFColumn("Z_HCX_nsys1", fEventAction->Get_HCX_Pos(1, 3)); // column Id = 4 //

 //Z-bars

        analysisManager->CreateNtupleIColumn("N_HCZ_nsys1", fEventAction->Get_HCZ_N(1)); // column Id = 4 // Число срабатываний
        analysisManager->CreateNtupleIColumn("AL_HCZ_nsys1", fEventAction->Get_HCZ_AL(1)); // column Id = 4 // Номер слоя

        analysisManager->CreateNtupleFColumn("E_HCZ_nsys1", fEventAction->Get_HCZ_Edep(1)); // column Id = 4 // Энерговыделение
        analysisManager->CreateNtupleFColumn("LO_HCZ_nsys1", fEventAction->Get_HCZ_LO(1)); // column Id = 6 // Световыход

        analysisManager->CreateNtupleFColumn("A1_HCZ_nsys1", fEventAction->Get_HCZ_A1(1)); // column Id = 8 // Амплутуды с противоположных торцов
        analysisManager->CreateNtupleFColumn("A2_HCZ_nsys1", fEventAction->Get_HCZ_A2(1)); // column Id = 10 //

        analysisManager->CreateNtupleFColumn("T1_HCZ_nsys1", fEventAction->Get_HCZ_T1(1)); // column Id = 8 // Времена с противоположных торцов
        analysisManager->CreateNtupleFColumn("T2_HCZ_nsys1", fEventAction->Get_HCZ_T2(1)); // column Id = 10 //

        analysisManager->CreateNtupleFColumn("L_HCZ_nsys1", fEventAction->Get_HCZ_TrackLength(1)); // column Id = 8 // Длина пробега

        analysisManager->CreateNtupleFColumn("ToF_HCZ_nsys1", fEventAction->Get_HCZ_ToF(1)); // column Id = 10 // Световыход


        analysisManager->CreateNtupleFColumn("X_HCZ_nsys1", fEventAction->Get_HCZ_Pos(1, 1)); // column Id = 4 // Локальная точка
        analysisManager->CreateNtupleFColumn("Y_HCZ_nsys1", fEventAction->Get_HCZ_Pos(1, 2)); // column Id = 4 //
        analysisManager->CreateNtupleFColumn("Z_HCZ_nsys1", fEventAction->Get_HCZ_Pos(1, 3)); // column Id = 4 //
 
//nsys2

        // X-bars
        analysisManager->CreateNtupleIColumn("N_HCX_nsys2", fEventAction->Get_HCX_N(2)); // column Id = 4 // Число срабатываний
        analysisManager->CreateNtupleIColumn("AL_HCX_nsys2", fEventAction->Get_HCX_AL(2)); // column Id = 4 // Номер слоя

        analysisManager->CreateNtupleFColumn("E_HCX_nsys2", fEventAction->Get_HCX_Edep(2)); // column Id = 4 // Энерговыделение
        analysisManager->CreateNtupleFColumn("LO_HCX_nsys2", fEventAction->Get_HCX_LO(2)); // column Id = 6 // Световыход

        analysisManager->CreateNtupleFColumn("A1_HCX_nsys2", fEventAction->Get_HCX_A1(2)); // column Id = 8 // Амплутуды с противоположных торцов
        analysisManager->CreateNtupleFColumn("A2_HCX_nsys2", fEventAction->Get_HCX_A2(2)); // column Id = 10 //

        analysisManager->CreateNtupleFColumn("T1_HCX_nsys2", fEventAction->Get_HCX_T1(2)); // column Id = 8 // Времена с противоположных торцов
        analysisManager->CreateNtupleFColumn("T2_HCX_nsys2", fEventAction->Get_HCX_T2(2)); // column Id = 10 //

        analysisManager->CreateNtupleFColumn("Lp_HCX_nsys2", fEventAction->Get_HCX_TrackLength(2)); // column Id = 8 // Длина пробега

        analysisManager->CreateNtupleFColumn("ToF_HCX_nsys2", fEventAction->Get_HCX_ToF(2)); // column Id = 10 // Световыход


        analysisManager->CreateNtupleFColumn("X_HCX_nsys2", fEventAction->Get_HCX_Pos(2, 1)); // column Id = 4 // Локальная точка
        analysisManager->CreateNtupleFColumn("Y_HCX_nsys2", fEventAction->Get_HCX_Pos(2, 2)); // column Id = 4 //
        analysisManager->CreateNtupleFColumn("Z_HCX_nsys2", fEventAction->Get_HCX_Pos(2, 3)); // column Id = 4 //

//Z-bars

        analysisManager->CreateNtupleIColumn("N_HCZ_nsys2", fEventAction->Get_HCZ_N(2)); // column Id = 4 // Число срабатываний
        analysisManager->CreateNtupleIColumn("AL_HCZ_nsys2", fEventAction->Get_HCZ_AL(2)); // column Id = 4 // Номер слоя

        analysisManager->CreateNtupleFColumn("E_HCZ_nsys2", fEventAction->Get_HCZ_Edep(2)); // column Id = 4 // Энерговыделение
        analysisManager->CreateNtupleFColumn("LO_HCZ_nsys2", fEventAction->Get_HCZ_LO(2)); // column Id = 6 // Световыход

        analysisManager->CreateNtupleFColumn("A1_HCZ_nsys2", fEventAction->Get_HCZ_A1(2)); // column Id = 8 // Амплутуды с противоположных торцов
        analysisManager->CreateNtupleFColumn("A2_HCZ_nsys2", fEventAction->Get_HCZ_A2(2)); // column Id = 10 //

        analysisManager->CreateNtupleFColumn("T1_HCZ_nsys2", fEventAction->Get_HCZ_T1(2)); // column Id = 8 // Времена с противоположных торцов
        analysisManager->CreateNtupleFColumn("T2_HCZ_nsys2", fEventAction->Get_HCZ_T2(2)); // column Id = 10 //

        analysisManager->CreateNtupleFColumn("Lp_HCZ_nsys2", fEventAction->Get_HCZ_TrackLength(2)); // column Id = 8 // Длина пробега

        analysisManager->CreateNtupleFColumn("ToF_HCZ_nsys2", fEventAction->Get_HCZ_ToF(2)); // column Id = 10 // Световыход


        analysisManager->CreateNtupleFColumn("X_HCZ_nsys2", fEventAction->Get_HCZ_Pos(2, 1)); // column Id = 4 // Локальная точка
        analysisManager->CreateNtupleFColumn("Y_HCZ_nsys2", fEventAction->Get_HCZ_Pos(2, 2)); // column Id = 4 //
        analysisManager->CreateNtupleFColumn("Z_HCZ_nsys2", fEventAction->Get_HCZ_Pos(2, 3)); // column Id = 4 //



/////////////////////////////

        analysisManager->CreateNtupleFColumn("Eplastic_thin_nsys1", fEventAction->GetPlasticThinEdep(1)); // column Id = 12
        analysisManager->CreateNtupleFColumn("Eplastic_thin_nsys2", fEventAction->GetPlasticThinEdep(2)); // column Id = 13

        analysisManager->CreateNtupleFColumn("LOplastic_thin_nsys1", fEventAction->GetPlasticThinLO(1)); // column Id = 14
        analysisManager->CreateNtupleFColumn("LOplastic_thin_nsys2", fEventAction->GetPlasticThinLO(2)); // column Id = 15

        analysisManager->CreateNtupleFColumn("Lplastic_thin_nsys1", fEventAction->GetPlasticThinTrackLength(1)); // column Id = 16
        analysisManager->CreateNtupleFColumn("Lplastic_thin_nsys2", fEventAction->GetPlasticThinTrackLength(2)); // column Id = 17

        analysisManager->CreateNtupleFColumn("A1plastic_thin_nsys1", fEventAction->GetPlasticThinA1(1)); // column Id = 8 // Амплутуды с противоположных торцов
        analysisManager->CreateNtupleFColumn("A1plastic_thin_nsys2", fEventAction->GetPlasticThinA1(2)); // column Id = 9

        analysisManager->CreateNtupleFColumn("T1plastic_thin_nsys1", fEventAction->GetPlasticThinT1(1)); // column Id = 8 // Времена с противоположных торцов
        analysisManager->CreateNtupleFColumn("T1plastic_thin_nsys2", fEventAction->GetPlasticThinT1(2)); // column Id = 9

        analysisManager->CreateNtupleFColumn("ToFplastic_thin_nsys1", fEventAction->GetPlasticThinToF(1)); // column Id = 18
        analysisManager->CreateNtupleFColumn("ToFplastic_thin_nsys2", fEventAction->GetPlasticThinToF(2)); // column Id = 19

        analysisManager->CreateNtupleFColumn("Xplastic_thin_nsys1", fEventAction->GetPlasticThinPos(1, 1)); // column Id = 4 // Локальная точка
        analysisManager->CreateNtupleFColumn("Yplastic_thin_nsys1", fEventAction->GetPlasticThinPos(1, 2)); // column Id = 4 //
        analysisManager->CreateNtupleFColumn("Zplastic_thin_nsys1", fEventAction->GetPlasticThinPos(1, 3)); // column Id = 4 //
        analysisManager->CreateNtupleFColumn("Xplastic_thin_nsys2", fEventAction->GetPlasticThinPos(2, 1)); // column Id = 4 // Локальная точка
        analysisManager->CreateNtupleFColumn("Yplastic_thin_nsys2", fEventAction->GetPlasticThinPos(2, 2)); // column Id = 4 //
        analysisManager->CreateNtupleFColumn("Zplastic_thin_nsys2", fEventAction->GetPlasticThinPos(2, 3)); // column Id = 4 //

        analysisManager->CreateNtupleFColumn("Thetaplastic_thin_nsys1", fEventAction->GetPlasticThinAngle(1, 1)); // column Id = 4 // Углы
        analysisManager->CreateNtupleFColumn("Phiplastic_thin_nsys1", fEventAction->GetPlasticThinAngle(1, 2)); // column Id = 4 //

        analysisManager->CreateNtupleFColumn("Thetaplastic_thin_nsys2", fEventAction->GetPlasticThinAngle(2, 1)); // column Id = 4 // Локальная точка
        analysisManager->CreateNtupleFColumn("Phiplastic_thin_nsys2", fEventAction->GetPlasticThinAngle(2, 2)); // column Id = 4 //

// Для пластиков LQ-поляриметра
/////////////////////////////

        analysisManager->CreateNtupleFColumn("Eplastic_LQ_nsys1", fEventAction->GetPlasticLQEdep(1)); // column Id = 12
        analysisManager->CreateNtupleFColumn("Eplastic_LQ_nsys2", fEventAction->GetPlasticLQEdep(2)); // column Id = 13

        analysisManager->CreateNtupleFColumn("LOplastic_LQ_nsys1", fEventAction->GetPlasticLQLO(1)); // column Id = 14
        analysisManager->CreateNtupleFColumn("LOplastic_LQ_nsys2", fEventAction->GetPlasticLQLO(2)); // column Id = 15

        analysisManager->CreateNtupleFColumn("Lplastic_LQ_nsys1", fEventAction->GetPlasticLQTrackLength(1)); // column Id = 16
        analysisManager->CreateNtupleFColumn("Lplastic_LQ_nsys2", fEventAction->GetPlasticLQTrackLength(2)); // column Id = 17

        analysisManager->CreateNtupleFColumn("A1plastic_LQ_nsys1", fEventAction->GetPlasticLQA1(1)); // column Id = 8 // Амплутуды с противоположных торцов
        analysisManager->CreateNtupleFColumn("A2plastic_LQ_nsys1", fEventAction->GetPlasticLQA2(1)); // column Id = 9
        analysisManager->CreateNtupleFColumn("A1plastic_LQ_nsys2", fEventAction->GetPlasticLQA1(2)); // column Id = 9
        analysisManager->CreateNtupleFColumn("A2plastic_LQ_nsys2", fEventAction->GetPlasticLQA2(2));

        analysisManager->CreateNtupleFColumn("T1plastic_LQ_nsys1", fEventAction->GetPlasticLQT1(1)); // column Id = 8 // Времена с противоположных торцов
        analysisManager->CreateNtupleFColumn("T2plastic_LQ_nsys1", fEventAction->GetPlasticLQT2(1)); // column Id = 9
        analysisManager->CreateNtupleFColumn("T1plastic_LQ_nsys2", fEventAction->GetPlasticLQT1(2)); // column Id = 9
        analysisManager->CreateNtupleFColumn("T2plastic_LQ_nsys2", fEventAction->GetPlasticLQT2(2)); // column Id = 9

        analysisManager->CreateNtupleFColumn("ToFplastic_LQ_nsys1", fEventAction->GetPlasticLQToF(1)); // column Id = 18
        analysisManager->CreateNtupleFColumn("ToFplastic_LQ_nsys2", fEventAction->GetPlasticLQToF(2)); // column Id = 19

        analysisManager->CreateNtupleFColumn("Xplastic_LQ_nsys1", fEventAction->GetPlasticLQPos(1, 1)); // column Id = 4 // Локальная точка
        analysisManager->CreateNtupleFColumn("Yplastic_LQ_nsys1", fEventAction->GetPlasticLQPos(1, 2)); // column Id = 4 //
        analysisManager->CreateNtupleFColumn("Zplastic_LQ_nsys1", fEventAction->GetPlasticLQPos(1, 3)); // column Id = 4 //
        analysisManager->CreateNtupleFColumn("Xplastic_LQ_nsys2", fEventAction->GetPlasticLQPos(2, 1)); // column Id = 4 // Локальная точка
        analysisManager->CreateNtupleFColumn("Yplastic_LQ_nsys2", fEventAction->GetPlasticLQPos(2, 2)); // column Id = 4 //
        analysisManager->CreateNtupleFColumn("Zplastic_LQ_nsys2", fEventAction->GetPlasticLQPos(2, 3)); // column Id = 4 //

        analysisManager->CreateNtupleFColumn("Thetaplastic_LQ_nsys1", fEventAction->GetPlasticLQAngle(1, 1)); // column Id = 4 // Углы
        analysisManager->CreateNtupleFColumn("Phiplastic_LQ_nsys1", fEventAction->GetPlasticLQAngle(1, 2)); // column Id = 4 //

        analysisManager->CreateNtupleFColumn("Thetaplastic_LQ_nsys2", fEventAction->GetPlasticLQAngle(2, 1)); // column Id = 4 // Локальная точка
        analysisManager->CreateNtupleFColumn("Phiplastic_LQ_nsys2", fEventAction->GetPlasticLQAngle(2, 2)); // column Id = 4 //


        /// Для трековых камер

        // Камера A
        analysisManager->CreateNtupleIColumn("nwa_nsys1", fEventAction->Get_W_N(1,1)); // column Id = 4 // Число срабатываний
        analysisManager->CreateNtupleFColumn("Xwa_nsys1", fEventAction->Get_W_Pos(1,1,1)); // column Id = 4 // Локальная точка
        analysisManager->CreateNtupleFColumn("Ywa_nsys1", fEventAction->Get_W_Pos(1,2,1));
        analysisManager->CreateNtupleFColumn("Zwa_nsys1", fEventAction->Get_W_Pos(1,3,1));
        analysisManager->CreateNtupleFColumn("Rwa_nsys1", fEventAction->Get_W_Pos(1,4,1));

        analysisManager->CreateNtupleIColumn("nwa_nsys2", fEventAction->Get_W_N(2)); // column Id = 4 // Число срабатываний
        analysisManager->CreateNtupleFColumn("Xwa_nsys2", fEventAction->Get_W_Pos(2,1,1)); // column Id = 4 // Локальная точка
        analysisManager->CreateNtupleFColumn("Ywa_nsys2", fEventAction->Get_W_Pos(2,2,1));
        analysisManager->CreateNtupleFColumn("Zwa_nsys2", fEventAction->Get_W_Pos(2,3,1));
        analysisManager->CreateNtupleFColumn("Rwa_nsys2", fEventAction->Get_W_Pos(2,4,1));

        // Камера B
        analysisManager->CreateNtupleIColumn("nwb_nsys1", fEventAction->Get_W_N(1,2)); // column Id = 4 // Число срабатываний
        analysisManager->CreateNtupleFColumn("Xwb_nsys1", fEventAction->Get_W_Pos(1,1,2)); // column Id = 4 // Локальная точка
        analysisManager->CreateNtupleFColumn("Ywb_nsys1", fEventAction->Get_W_Pos(1,2,2));
        analysisManager->CreateNtupleFColumn("Zwb_nsys1", fEventAction->Get_W_Pos(1,3,2));
        analysisManager->CreateNtupleFColumn("Rwb_nsys1", fEventAction->Get_W_Pos(1,4,2));



        analysisManager->CreateNtupleIColumn("nwb_nsys2", fEventAction->Get_W_N(2,2)); // column Id = 4 // Число срабатываний
        analysisManager->CreateNtupleFColumn("Xwb_nsys2", fEventAction->Get_W_Pos(2,1,2)); // column Id = 4 // Локальная точка
        analysisManager->CreateNtupleFColumn("Ywb_nsys2", fEventAction->Get_W_Pos(2,2,2));
        analysisManager->CreateNtupleFColumn("Zwb_nsys2", fEventAction->Get_W_Pos(2,3,2));
        analysisManager->CreateNtupleFColumn("Rwb_nsys2", fEventAction->Get_W_Pos(2,4,2));


        // Камера C
        analysisManager->CreateNtupleIColumn("nwc_nsys1", fEventAction->Get_W_N(1,3)); // column Id = 4 // Число срабатываний
        analysisManager->CreateNtupleFColumn("Xwc_nsys1", fEventAction->Get_W_Pos(1,1,3)); // column Id = 4 // Локальная точка
        analysisManager->CreateNtupleFColumn("Ywc_nsys1", fEventAction->Get_W_Pos(1,2,3));
        analysisManager->CreateNtupleFColumn("Zwc_nsys1", fEventAction->Get_W_Pos(1,3,3));
        analysisManager->CreateNtupleFColumn("Rwc_nsys1", fEventAction->Get_W_Pos(1,4,3));

        analysisManager->CreateNtupleIColumn("nwc_nsys2", fEventAction->Get_W_N(2,3)); // column Id = 4 // Число срабатываний
        analysisManager->CreateNtupleFColumn("Xwc_nsys2", fEventAction->Get_W_Pos(2,1,3)); // column Id = 4 // Локальная точка
        analysisManager->CreateNtupleFColumn("Ywc_nsys2", fEventAction->Get_W_Pos(2,2,3));
        analysisManager->CreateNtupleFColumn("Zwc_nsys2", fEventAction->Get_W_Pos(2,3,3));
        analysisManager->CreateNtupleFColumn("Rwc_nsys2", fEventAction->Get_W_Pos(2,4,3));


        /// Для вершинных камер

        analysisManager->CreateNtupleIColumn("nvc_nsys1", fEventAction->Get_VC_N(1)); // column Id = 4 // Число срабатываний
        analysisManager->CreateNtupleFColumn("Xvc_nsys1", fEventAction->Get_VC_Pos(1,1)); // column Id = 4 // Локальная точка
        analysisManager->CreateNtupleFColumn("Yvc_nsys1", fEventAction->Get_VC_Pos(1,2));
        analysisManager->CreateNtupleFColumn("Zvc_nsys1", fEventAction->Get_VC_Pos(1,3));
        analysisManager->CreateNtupleFColumn("Rvc_nsys1", fEventAction->Get_VC_Pos(1,4));

        analysisManager->CreateNtupleIColumn("nvc_nsys2", fEventAction->Get_VC_N(2)); // column Id = 4 // Число срабатываний
        analysisManager->CreateNtupleFColumn("Xvc_nsys2", fEventAction->Get_VC_Pos(2,1)); // column Id = 4 // Локальная точка
        analysisManager->CreateNtupleFColumn("Yvc_nsys2", fEventAction->Get_VC_Pos(2,2));
        analysisManager->CreateNtupleFColumn("Zvc_nsys2", fEventAction->Get_VC_Pos(2,3));
        analysisManager->CreateNtupleFColumn("Rvc_nsys2", fEventAction->Get_VC_Pos(2,4));

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
