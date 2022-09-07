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

class G4ParticleGun;
class G4Event;

namespace Cosmic
{

/// The primary generator action class with particle gum.
///
/// It defines a single particle which hits the calorimeter
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class
/// (see the macros provided with this example).

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event* event) override;

  // set methods
  void SetRandomFlag(G4bool value);

private:
  G4ParticleGun* fParticleGun = nullptr; // G4 particle gun

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
                                    G4double &momentum, G4double &kinetic_energy) {

        G4double theta_temp = 0.0;
        G4double momentum_temp = 0.0;
        G4double f_temp = 0.0;
        G4double r1, r2, r3;

        G4bool flaq = false;

        while (flaq == false) {

            r1 = G4UniformRand();
            r2 = G4UniformRand();
            r3 = G4UniformRand();

            theta_temp = initial_theta + r1 * (final_theta - initial_theta);
            momentum_temp = initial_momentum + r2 * (final_momentum - initial_momentum);
            f_temp = max_f * r3;

            if (f_temp <= f_theta_energy_sin(theta_temp, momentum_temp)) {
                flaq = true;
            } else {
                flaq = false;
            }
        }

        theta = theta_temp;

        G4double gamma = sqrt(pow(momentum_temp / m_mu, 2.0) + 1.0);
        momentum = momentum_temp;
        kinetic_energy = m_mu * (gamma - 1.0); // возвращаем не импульс, а кинетическую энергию

    }

};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
