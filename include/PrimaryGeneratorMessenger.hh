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

#ifndef CosmicPrimaryGeneratorMessenger_h
#define CosmicPrimaryGeneratorMessenger_h 1

//#include "PrimaryGeneratorAction.hh"
#include "G4UImessenger.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"

class G4ParticleGun;
class G4ParticleTable;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


namespace Cosmic {

    class PrimaryGeneratorAction;

    class PrimaryGeneratorMessenger : public G4UImessenger {
    public:
       // PrimaryGeneratorMessenger(PrimaryGeneratorAction *);
        PrimaryGeneratorMessenger(PrimaryGeneratorAction * GeneratorAction);

        ~PrimaryGeneratorMessenger();

        void SetNewValue(G4UIcommand * command, G4String newValues);

    private:
        G4ParticleGun * fParticleGun;
        G4ParticleTable * particleTable;
        PrimaryGeneratorAction *fGeneratorAction;

    private: //commands

        G4UIcmdWithAnInteger *GenbosBoolCmd;

        G4UIcmdWithAString *RndmCmd;
        G4UIcmdWithAString *CountCmd;
        G4UIcmdWithAnInteger *CStepCmd;
        G4UIcmdWithAnInteger *ModeCmd;
        G4UIcmdWithAString *VertexCmd;
        G4UIcmdWithAnInteger *EntryCmd;

        G4UIcmdWithADoubleAndUnit *EgminCmd;
        G4UIcmdWithADoubleAndUnit *EgmaxCmd;
    };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
#endif //CosmicPrimaryGeneratorMessenger_h
