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
/// \file PrimaryGeneratorAction.hh
/// \brief Definition of the PrimaryGeneratorAction class

#ifndef CosmicPrimaryGeneratorAction_h
#define CosmicPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include <Randomize.hh>
#include "G4ThreeVector.hh"
#include "Constants.hh"

#ifdef GENBOS

#include "Genbos.hh"

#endif

class G4ParticleGun;

class G4Event;

namespace Cosmic {

/// The primary generator action class with particle gum.
///
/// It defines a single particle which hits the calorimeter
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class
/// (see the macros provided with this example).

    class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
    public:
        PrimaryGeneratorAction();

        ~PrimaryGeneratorAction() override;

        void GeneratePrimaries(G4Event *event) override;
        void GenerateCosmic(G4Event* event); // for cosmic-generator
        void ShowParticleTable();

        // set methods
        void SetRandomFlag(G4bool value);

        ///// FOR GENBOS

        void GenerateGenbos(G4Event* event); // for pn-generator GENBOS
        void SetFileNum(G4int val){FileNum=val;}
        inline void SetEgMin(G4double val){
            EgMin=val/GeV;G4int n=2;
            genbos_beam_(&n,&EgMin, &EgMax);
        }
        inline void SetEgMax(G4double val){
            EgMax=val/GeV;G4int n=2;
            genbos_beam_(&n,&EgMin, &EgMax);
        }

        inline void SetEgMin(G4double val_min, G4double val_max){
            EgMin=val_min/GeV;EgMax=val_max/GeV;
            G4int n=2;
            genbos_beam_(&n,&EgMin, &EgMax);
        }
        //////


    private:
        G4ParticleGun *fParticleGun = nullptr; // G4 particle gun

        G4bool fRandomDirection;

        G4double m_mu = 105.658374524; // в МэВ
// другое распрределение по углу и энергии из работы Шебалина
        G4double f_theta_energy_sin(G4double theta, G4double momentum) {
            G4double temp = 0.0;
            G4double k = 2.26;
            G4double p0 = 500.0; //МэВ

            G4double a = 2.86;
            G4double b = 1.54;
            G4double sigma = (a - cos(theta)) / b;
            //  long double A = 1.0 / 0.158363; // константа интегрирования 0..pi/2
            G4double A = 1.0;
            // здесь без нормировки
            temp = A * pow(cos(theta), k) * exp(-pow(log(momentum / p0), 2.0) / (2.0 * pow(sigma, 2.0))) * sin(theta);
            return temp;
        }

        void
        random_Neumann_theta_energy_sin(G4double initial_theta, G4double final_theta, G4double initial_momentum,
                                        G4double final_momentum, G4double max_f, G4double &theta,
                                        G4double &momentum, G4double &kinetic_energy);


        ///// FOR GENBOS

    private:
        G4ThreeVector GenVertex(G4double, G4double, G4double, G4double); //  метод генерации точки вылета из пучка по гауссу

        G4int FileNum;
        G4double Egamma;
        G4float EgMin;
        G4float EgMax;

        G4String part_name[64];

        void PrepareNames();

        G4float efot; //photon energy, GeV
        G4int nreac; // reaction index
        G4float vx;
        G4float vy;
        G4float vz;
        G4int np; // number of generated particles
        G4int idg[11];   //[np] // list of particle's indexes
        G4float cx[11];   //[np] // lists of momentum components
        G4float cy[11];   //[np]
        G4float cz[11];   //[np]

        //////


    };

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
