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
#include "G4ios.hh"

#ifdef GENBOS

#include "Genbos.hh"

#endif

namespace Cosmic {

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    PrimaryGeneratorAction::PrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(), fParticleGun(0) {
        G4int nofParticles = 1;
        fParticleGun = new G4ParticleGun(nofParticles);

        fRandomDirection = true;

        // default particle kinematic
        //

        auto particleTable = G4ParticleTable::GetParticleTable();

        auto particleDefinition = particleTable->FindParticle("mu-");
        fParticleGun->SetParticleDefinition(particleDefinition);


        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, -1.0, 0.0));
        fParticleGun->SetParticleEnergy(2. * GeV);

       // ShowParticleTable();

    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    PrimaryGeneratorAction::~PrimaryGeneratorAction() {
        delete fParticleGun;
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {
        // This function is called at the begining of event

        // In order to avoid dependence of PrimaryGeneratorAction
        // on DetectorConstruction class we get world volume
        // from G4LogicalVolumeStore
        //
        G4double worldZHalfLength = 0.;
        auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

        // Check that the world volume has box shape
        G4Box *worldBox = nullptr;
        if (worldLV) {
            worldBox = dynamic_cast<G4Box *>(worldLV->GetSolid());
        }

        if (worldBox) {
            worldZHalfLength = worldBox->GetZHalfLength();
        } else {
            G4ExceptionDescription msg;
            msg << "World volume of box shape not found." << G4endl;
            msg << "Perhaps you have changed geometry." << G4endl;
            msg << "The gun will be place in the center.";
            G4Exception("PrimaryGeneratorAction::GeneratePrimaries()",
                        "MyCode0002", JustWarning, msg);
        }

        GenerateCosmic(anEvent);

    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


    void PrimaryGeneratorAction::ShowParticleTable() {

#define SW(N)    " "<<std::setw( N )
        G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
        G4cout << std::setprecision(1) << std::setiosflags(std::ios::fixed);
        for (G4int ii = 0; ii < particleTable->entries(); ii++) {
            G4cout << SW(3) << ii;
            G4cout << SW(10) << particleTable->GetParticle(ii)->GetParticleName();
            G4cout << SW(6) << particleTable->GetParticle(ii)->GetPDGMass();
            G4cout << SW(2) << particleTable->GetParticle(ii)->GetPDGCharge();
            G4cout << SW(3) << particleTable->GetParticle(ii)->GetPDGSpin();
            G4cout << SW(4) << particleTable->GetParticle(ii)->GetPDGEncoding();
            G4cout << G4endl;
        }
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void PrimaryGeneratorAction::random_Neumann_theta_energy_sin(G4double initial_theta, G4double final_theta,
                                                                 G4double initial_momentum,
                                                                 G4double final_momentum, G4double max_f,
                                                                 G4double &theta,
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

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void PrimaryGeneratorAction::GenerateCosmic(G4Event *event) {

        G4ThreeVector initial_pos(0.0, 0.0, 0.0);
        G4ThreeVector initial_momentum_direction(0.0, 0.0, 0.0);

        //fRandomDirection=false;

        if (fRandomDirection) {

            G4double x_initial_pos, y_initial_pos = 4.5 * m, z_initial_pos;
            G4double x_initial = -0.5 * m, y_initial = 4.5 * m, z_initial = 0.5 * m;
            G4double x_final = 0.5 * m, y_final = 4.5 * m, z_final = 1.0 * m;
            G4double vx, vy, vz;
            G4double theta, psi;
            G4double initial_theta = 0.0, final_theta = 30. * M_PI / 180.;
            G4double initial_momentum = 10.0, final_momentum = 10000.0; // в МэВ
            G4double max_f = 1.0;
            G4double momentum, kinetic_energy;

            random_Neumann_theta_energy_sin(initial_theta, final_theta, initial_momentum,
                                            final_momentum, max_f, theta,
                                            momentum, kinetic_energy);
            psi = 2.0 * M_PI * G4UniformRand();

            vx = sin(theta) * cos(psi);
            vz = sin(theta) * sin(psi);
            vy = -cos(theta);

            x_initial_pos = x_initial + (x_final - x_initial) * G4UniformRand();
            z_initial_pos = z_initial + (z_final - z_initial) * G4UniformRand();


            initial_momentum_direction = G4ThreeVector(vx, vy, vz);

            //  G4cerr<<"momentum = "<< momentum<< "kinetic_energy ="<< kinetic_energy << " theta = "<< theta* 180./M_PI << G4endl;

            fParticleGun->SetParticleEnergy(kinetic_energy * MeV);
            //  fParticleGun->SetParticleMomentum(momentum*MeV);
            initial_pos = G4ThreeVector(x_initial_pos, y_initial_pos, z_initial_pos);
            //   initial_pos = G4ThreeVector(0., 4.5 * m, 1.0 * m - 200 * mm);

        } else {

            // Set gun position
            initial_momentum_direction = G4ThreeVector(0.0, -1.0, 0.0);
            initial_pos = G4ThreeVector(0., 4.5 * m, 1.0 * m + 100 * mm);
        }

        fParticleGun->SetParticleMomentumDirection(initial_momentum_direction);
        fParticleGun->SetParticlePosition(initial_pos);
        fParticleGun->GeneratePrimaryVertex(event);

    }

    void PrimaryGeneratorAction::GenerateGenbos(G4Event *event) {


        G4ThreeVector vertex;
        G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
        G4ParticleDefinition *particle;
        G4ParticleMomentum Momentum;

        G4int qq = 0;

        while (qq != 3) {
            genbos_event_(&efot, &nreac, &np, idg, cx, cy, cz);
            for (G4int i = 0; i < np; i++) {
                Momentum = G4ParticleMomentum(cx[i], cy[i], cz[i]);
                if (Momentum.theta() < THDN) continue;
                if (fabs(Momentum.phi() - PHIC) < DPHI) qq |= 1;
                if (fabs(Momentum.phi() + PHIC) < DPHI) qq |= 2;
            }
        }

// Vertex position  ------------

        G4double Xbeam = 0., Ybeam = 0.;
        vertex = GenVertex(Xbeam, Ybeam, Xsigma_beam, Ysigma_beam);

        Egamma = efot * GeV;

//G4cout<<"**** Eg="<<Egamma<<"  np="<<np<<G4endl;

        for (int i = 0; i < np; i++) {
            G4double mass, energy;
            Momentum = G4ParticleMomentum(cx[i], cy[i], cz[i]);
            Momentum *= GeV;

            particle = particleTable->FindParticle(part_name[idg[i]]);

            mass = particle->GetPDGMass();
            energy = sqrt(Momentum.mag2() + pow(mass, 2.0)) - mass;
            fParticleGun->SetParticleEnergy(energy);
            fParticleGun->SetParticleMomentumDirection(Momentum);

//     particleGun = new G4ParticleGun(1);
//     particleGun->SetParticleMomentum(Mom);

            fParticleGun->SetParticleDefinition(particle);
            fParticleGun->SetParticlePosition(vertex);
            fParticleGun->GeneratePrimaryVertex(event);

//     delete particleGun;
        }

        /*
        pn2020EventInfo* info = new pn2020EventInfo();
//   pn2020EventInfo* info =(pn2020EventInfo*)anEvent->GetUserInformation();
        info->SetEgamma(Egamma);
        info->SetNreac(nreac);
        info->SetNp(np);
        info->SetEntry(FileNum);
        event->SetUserInformation(info);

         */
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    G4ThreeVector PrimaryGeneratorAction::GenVertex(G4double X, G4double Y, G4double sigmaX, G4double sigmaY) {
        static G4ThreeVector position;
        G4double pos_x, pos_y, pos_z, pos_t;
#ifdef XY_RAND
        pos_x = G4RandGauss::shoot(X, sigmaX);
        pos_y = G4RandGauss::shoot(Y, sigmaY);
#else
        pos_x = X;
        pos_y = Y;
        sigmaX=0;
        sigmaY=0;
#endif
        do {
            pos_z = G4UniformRand();
            pos_t = G4UniformRand();
        } while (pos_t > (1.0 - pos_z));
        pos_z *= (G4UniformRand() > 0.5) ? cell_z_size_2 : -cell_z_size_2;
        // pos_z *= cm;

        position = G4ThreeVector(pos_x, pos_y, pos_z);
        return position;
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void PrimaryGeneratorAction::PrepareNames() {

//       PDG-id                 g, e+, e-, mu+, mu-, pi0, pi+, pi-, K0L, K+, K-,
        G4int partid[64] = {1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12,
//                 n, p, K0S,eta,Lambda,sigma+,
                            13, 14, 16, 17, 18, 19,
//               sigma0,sigma-,omega,rho+,rho0,rho-
                            20, 21, 60, 58, 57, 59,
//               deuteron
                            45, -1};
        const char partname[64][12] =
                {"gamma", "e+", "e-", "mu+", "mu-", "pi0", "pi+", "pi-", "kaon0L", "kaon+", "kaon-",
                 "neutron", "proton", "kaon0S", "eta", "lambda", "sigma+",
                 "sigma0", "sigma-", "omega", "rho+", "rho0", "rho-", "deuteron"};

        for (G4int i = 0; i < 64; i++) {
            part_name[i] = "";
        }
        G4int i = 0;
        while (partid[i] >= 0 && partid[i] < 64) {
            part_name[partid[i]] = partname[i];
            i++;
        }

    }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


}
