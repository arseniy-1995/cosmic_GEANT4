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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the Cosmic::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

namespace Cosmic
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(): G4VUserPrimaryGeneratorAction(), fParticleGun(0)
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

    fRandomDirection = true;

  // default particle kinematic
  //
  auto particleDefinition
    = G4ParticleTable::GetParticleTable()->FindParticle("mu-");
  fParticleGun->SetParticleDefinition(particleDefinition);



  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0,-1.0,0.0));
  fParticleGun->SetParticleEnergy(2.*GeV);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume
  // from G4LogicalVolumeStore
  //
  G4double worldZHalfLength = 0.;
  auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

  // Check that the world volume has box shape
  G4Box* worldBox = nullptr;
  if (  worldLV ) {
    worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  }

  if ( worldBox ) {
    worldZHalfLength = worldBox->GetZHalfLength();
  }
  else  {
    G4ExceptionDescription msg;
    msg << "World volume of box shape not found." << G4endl;
    msg << "Perhaps you have changed geometry." << G4endl;
    msg << "The gun will be place in the center.";
    G4Exception("PrimaryGeneratorAction::GeneratePrimaries()",
      "MyCode0002", JustWarning, msg);
  }
//fRandomDirection=false;

    if(fRandomDirection) {

        G4double x_initial_pos, y_initial_pos = 4.5*m, z_initial_pos;
        G4double vx, vy, vz;
        G4double theta, psi;
        G4double initial_theta = 0.0, final_theta = 30.*M_PI/180. ;
        G4double initial_momentum = 10.0, final_momentum = 10000.0; // в МэВ
        G4double max_f = 1.0;
        G4double momentum, kinetic_energy;

        random_Neumann_theta_energy_sin(initial_theta, final_theta, initial_momentum,
                                        initial_momentum, max_f, theta,
                                        momentum, kinetic_energy);
        psi = 2.0 * M_PI * G4UniformRand();

        vx = sin(theta) * cos(psi);
        vz = sin(theta) * sin(psi);
        vy = -cos(theta);

        x_initial_pos=4.0*m*G4UniformRand() -2.0*m;
        z_initial_pos=4.0*m*G4UniformRand() -2.0*m;

        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(vx,vy,vz));
        fParticleGun->SetParticleEnergy(kinetic_energy);

        fParticleGun->SetParticlePosition(G4ThreeVector(x_initial_pos, y_initial_pos, z_initial_pos));

    }else{

        // Set gun position
        fParticleGun
                ->SetParticlePosition(G4ThreeVector(0., 4.5 * m, 0.));
    }





  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
