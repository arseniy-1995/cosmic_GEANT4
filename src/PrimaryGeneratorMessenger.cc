//
// Created by Арсений Юрченко on 19.10.2022.
//
// Этот класс позволяет устанавливать параметры с помощью макросов

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
/// \file PrimaryGeneratorMessenger.cc
/// \brief Implementation of the Cosmic::PrimaryGeneratorMessenger class

#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"

//#include "G4ParticleGunMessenger.hh"
#include "G4ParticleGun.hh"
#include "G4Geantino.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleTable.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "globals.hh"
//#include <iostream.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace Cosmic {

    //  PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction *Gun) :
    //         GeneratorAction(Gun) {

    PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction *GeneratorAction)
            : fGeneratorAction(GeneratorAction) {


        GenbosBoolCmd = new G4UIcmdWithAnInteger("/gun/GenbosBool", this);
        GenbosBoolCmd->SetGuidance("Genbos Generator on/off");
        GenbosBoolCmd->SetGuidance("  Choice : off(default), on");
        GenbosBoolCmd->SetParameterName("GenbosBool", false);
        GenbosBoolCmd->SetRange("GenbosBool>=0");
        GenbosBoolCmd->SetDefaultValue(0);

        RndmCmd = new G4UIcmdWithAString("/gun/random", this);
        RndmCmd->SetGuidance("Randomize using current time.");
        RndmCmd->SetGuidance("  Choice : off(default), on");
        RndmCmd->SetParameterName("choice", true);
        RndmCmd->SetDefaultValue("off");
        RndmCmd->SetCandidates("off on");
//  RndmCmd->AvailableForStates(PreInit,Idle);

        CountCmd = new G4UIcmdWithAString("/gun/count", this);
        CountCmd->SetGuidance("Display event count on terminal.");
        CountCmd->SetGuidance("  Choice : off(default), on");
        CountCmd->SetParameterName("choice", true);
        CountCmd->SetDefaultValue("off");
        CountCmd->SetCandidates("off on");

        CStepCmd = new G4UIcmdWithAnInteger("/gun/CStep", this);
        CStepCmd->SetGuidance("Count Step");
        CStepCmd->SetParameterName("CStep", false);
        CStepCmd->SetRange("CStep>=5");
        CStepCmd->SetDefaultValue(100);

        ModeCmd = new G4UIcmdWithAnInteger("/gun/mode", this);
        ModeCmd->SetGuidance("0: all gamma+d; 1: elastic ed; 34: pn");
        ModeCmd->SetParameterName("mode", false);
        //  ModeCmd->SetRange("");
        ModeCmd->SetDefaultValue(0);


        TargetTypeCmd = new G4UIcmdWithAnInteger("/gun/TargetType", this);
        TargetTypeCmd->SetGuidance("0: neutron; 1: proton; 2: deuteron");
        TargetTypeCmd->SetParameterName("mode", false);
        //  TargetTypeCmd->SetRange("");
        TargetTypeCmd->SetDefaultValue(2);


        BeamSpectrumCmd = new G4UIcmdWithAnInteger("/gun/BeamSpectrum", this);
        BeamSpectrumCmd->SetGuidance("0: gaussian; 2: bremsstrahlung; 3: uniform");
        BeamSpectrumCmd->SetParameterName("mode", false);
        //  BeamSpectrumCmd->SetRange("");
        BeamSpectrumCmd->SetDefaultValue(2);


        VertexCmd = new G4UIcmdWithAString("/gun/vertex", this);
        VertexCmd->SetGuidance("Generate vertex in the cell?");
        VertexCmd->SetGuidance("  Choice : off(default), on");
        VertexCmd->SetParameterName("choice", true);
        VertexCmd->SetDefaultValue("off");
        VertexCmd->SetCandidates("off on");

        EntryCmd = new G4UIcmdWithAnInteger("/gun/FileNum", this);
        EntryCmd->SetGuidance("File/Process index");
        EntryCmd->SetParameterName("FileNum", false);
        EntryCmd->SetRange("FileNum>=0");
        EntryCmd->SetDefaultValue(0);

        EgminCmd = new G4UIcmdWithADoubleAndUnit("/gun/EgMin", this);
        EgminCmd->SetGuidance("Set min Egamma");
        EgminCmd->SetParameterName("EgMin", false);
        EgminCmd->SetRange("EgMin>100.");
        EgminCmd->SetUnitCategory("Energy");

        EgmaxCmd = new G4UIcmdWithADoubleAndUnit("/gun/EgMax", this);
        EgmaxCmd->SetGuidance("Set max Egamma");
        EgmaxCmd->SetParameterName("EgMax", false);
        EgmaxCmd->SetRange("EgMax>200.");
        EgmaxCmd->SetUnitCategory("Energy");

    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger() {
        delete GenbosBoolCmd;
        delete RndmCmd;
        delete CountCmd;
        delete CStepCmd;
        delete ModeCmd;
        delete TargetTypeCmd;
        delete BeamSpectrumCmd;
        delete VertexCmd;
        delete EntryCmd;
        delete EgminCmd;
        delete EgmaxCmd;
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand *command, G4String newValue) {

#ifdef GENBOS

        if (command == GenbosBoolCmd) {
            fGeneratorAction->SetGenbosBool(GenbosBoolCmd->GetNewIntValue(newValue));
        }
        if (command == RndmCmd) {
            fGeneratorAction->SetRndmFlag(newValue);
            fGeneratorAction->DoRandomize();
        }
        if (command == CountCmd) { fGeneratorAction->SetCountFlag(newValue); }
        if (command == CStepCmd)
        {
            fGeneratorAction->SetCStep(CStepCmd->GetNewIntValue(newValue));
        }
        if (command == ModeCmd)
        {
            fGeneratorAction->SetMode(ModeCmd->GetNewIntValue(newValue));
        }
        if (command == TargetTypeCmd)
        {
            fGeneratorAction->SetTargetType(TargetTypeCmd->GetNewIntValue(newValue));
        }
        if (command == BeamSpectrumCmd)
        {
            fGeneratorAction->SetBeamSpectrum(BeamSpectrumCmd->GetNewIntValue(newValue)); }
        if (command == VertexCmd) { fGeneratorAction->SetVertexFlag(newValue); }
        if (command == EntryCmd) { fGeneratorAction->SetFileNum(EntryCmd->GetNewIntValue(newValue)); }
        if (command == EgminCmd) { fGeneratorAction->SetEgMin(EgminCmd->GetNewDoubleValue(newValue)); }
        if (command == EgmaxCmd) { fGeneratorAction->SetEgMax(EgmaxCmd->GetNewDoubleValue(newValue)); }

#endif
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}