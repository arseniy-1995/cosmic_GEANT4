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
//#include "PrimaryGeneratorMessenger.hh"
#include "EventInfo.hh"

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
#include <unistd.h>
#include <fcntl.h>
#include "G4MTRunManager.hh"
#include "G4AutoLock.hh"
#endif

namespace Cosmic {


    namespace	{
        G4Mutex	aMutex	=	G4MUTEX_INITIALIZER;
    }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    PrimaryGeneratorAction::PrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(), fParticleGun(0),
                                                       GenbosBool(1), cstep(100), countFlag("off"), rndmFlag("off"),
                                                       vertexFlag("off"), Mode(0),
                                                       FileNum(0),
                                                    //   FileNum(G4Threading::G4GetThreadId()),
                                                       EgMin(400*MeV), EgMax(650*MeV) {

       // G4cerr<<"!!!!!!!!"<<G4Threading::G4GetThreadId()<< std::endl;
        G4int nofParticles = 1;
        fParticleGun = new G4ParticleGun(nofParticles);

        gunMessenger = new PrimaryGeneratorMessenger(this); // Вызов конструктора для класса пользовательского ввода

        fRandomDirection = true;

        // default particle kinematic
        //

        auto particleTable = G4ParticleTable::GetParticleTable();

        auto particleDefinition = particleTable->FindParticle("mu-");
        fParticleGun->SetParticleDefinition(particleDefinition);


        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, -1.0, 0.0));
        fParticleGun->SetParticleEnergy(2. * GeV);

        // ShowParticleTable();

        if (GenbosBool == 1) {
            G4AutoLock lock(&aMutex);
           	genbos_start_(&FileNum);
           // lock.unlock();
               PrepareNames();
        }
      //  genbos_start_(&FileNum);


       // PrepareNames();
       // genbos_start_(&FileNum);
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    PrimaryGeneratorAction::~PrimaryGeneratorAction() {
        //if (GenbosBool == 1) {
            G4AutoLock lock(&aMutex);
      //  }
       //genbos_stop_();
      // lock.unlock();
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

        if (GenbosBool == 0) {
           //   GenerateCosmic(anEvent);
            //GenerateLowQ_method1(anEvent);
           // GenerateLowQ_method2(anEvent);
            GenerateProton(anEvent);
           // GenerateGenbos(anEvent);
        } else if (GenbosBool == 1) {

            G4AutoLock lock(&aMutex);
          //  G4cout<<"!!!"<< std::endl;
            //genbos_start_(&FileNum);
            GenerateGenbos(anEvent);
         //   genbos_stop_();
           // lock.unlock();
        }

        if (countFlag == "on") {
            G4int nev = anEvent->GetEventID() + 1;
            if (nev % cstep == 0)
                fprintf(stderr, "\r event #%05d/%d \n    ", nev, FileNum);
        }

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

    void PrimaryGeneratorAction::random_Neumann_LQ_method1(G4double initial_theta_e, G4double final_theta_e,
                                                   G4double initial_zz_cell, G4double final_zz_cell,
                                                   G4double Pzz1, G4double Pzz2, G4double max_f,
                                                   G4double &theta_electron, G4double &phi_electron, G4double &energy_electron,
                                                   G4double &theta_deuteron, G4double &phi_deuteron, G4double &energy_deuteron,
                                                   G4double &Pzz,
                                                   G4double &xx_cell, G4double &yy_cell, G4double &zz_cell) {

        G4double theta_e_temp = 0.0;
        G4double phi_e_temp = 0.0;
        G4double xx_cell_temp =0.0, yy_cell_temp=0.0, zz_cell_temp = 0.0;
        G4double Pzz_temp = 0.0;
        G4double f_temp = 0.0;
        G4double r1 = 0., r2 = 0., r3 = 0., r4 = 0., r5 = 0.;



        G4double f_dsdo =0.0;
        G4bool flaq = false;

        while (flaq == false) {
            r1 = G4UniformRand();
            r2 = G4UniformRand();
            r3 = G4UniformRand();
            r4 = G4UniformRand();
            r5 = G4UniformRand();

            theta_e_temp = initial_theta_e + r1 * (final_theta_e - initial_theta_e);
            xx_cell_temp = 0.1*G4RandGauss::shoot(meanX_beam, Xsigma_beam) / mm;
            yy_cell_temp = 0.1*G4RandGauss::shoot(meanY_beam, Ysigma_beam) / mm;
            zz_cell_temp = initial_zz_cell + r2 * (final_zz_cell - initial_zz_cell);
            Pzz_temp = (r3 > 0.5) ? Pzz1 : Pzz2;
            phi_e_temp = 2.0 * M_PI * r4;
            f_temp = max_f * r5;
            f_dsdo= dsdo_LQ(theta_e_temp, Pzz_temp, zz_cell_temp);
         //   G4cerr << "!!! dsdo= " << f_dsdo <<" rnd = " << f_temp <<std::endl;
            if (f_temp <= f_dsdo) {
                flaq = true;
            } else {
                flaq = false;
            }
        }

        // Returns
        theta_electron = theta_e_temp;
        G4double rc = 1. + 2. * Ebeam * pow(sin(theta_electron / 2.), 2.) / Md; // фактор отдачи
        energy_deuteron = (1. - 1. / rc) * Ebeam; // кин энергия дейтрона
        energy_electron = Ebeam / rc; // энергия электрона после рассеяния
        theta_deuteron = atan(sin(theta_electron) / (rc - cos(theta_electron)));
        phi_electron = phi_e_temp;
        phi_deuteron = M_PI - phi_electron;
        zz_cell = zz_cell_temp;
        xx_cell = xx_cell_temp;
        yy_cell = yy_cell_temp;
    }



    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void PrimaryGeneratorAction::random_Neumann_LQ_method2(G4double initial_x_counter_e, G4double final_x_counter_e,
                                                           G4double initial_y_counter_e, G4double final_y_counter_e,
                                                           G4double initial_zz_cell, G4double final_zz_cell,
                                                           G4double Pzz1, G4double Pzz2, G4double max_f,
                                                           G4double &theta_electron, G4double &phi_electron, G4double &energy_electron,
                                                           G4double &theta_deuteron, G4double &phi_deuteron, G4double &energy_deuteron,
                                                           G4double &Pzz,
                                                           G4double &xx_cell, G4double &yy_cell, G4double &zz_cell) {

        G4double theta_e_temp = 0.0, phi_e_temp = 0.0;
        G4double x_counter_e_temp = 0.0;
        G4double y_counter_e_temp = 0.0;
        G4double xx_cell_temp = 0.0, yy_cell_temp = 0.0, zz_cell_temp = 0.0;
        G4double Pzz_temp = 0.0;
        G4double f_temp = 0.0;
        G4double r1 = 0., r2 = 0., r3 = 0., r4 = 0., r5 = 0.;


        G4double f_dsdo =0.0;
        G4bool flaq = false;

        while (flaq == false) {

            r1 = G4UniformRand();
            r2 = G4UniformRand();
            r3 = G4UniformRand();
            r4 = G4UniformRand();
            r5 = G4UniformRand();

            x_counter_e_temp = initial_x_counter_e + r1 * (final_x_counter_e - initial_x_counter_e);
            y_counter_e_temp = initial_x_counter_e + r2 * (final_y_counter_e - initial_y_counter_e);
            xx_cell_temp = 0.1*G4RandGauss::shoot(meanX_beam, Xsigma_beam) / mm;
            yy_cell_temp = 0.1*G4RandGauss::shoot(meanY_beam, Ysigma_beam) / mm;
            zz_cell_temp = initial_zz_cell + r3 * (final_zz_cell - initial_zz_cell);
            Pzz_temp = (r4 > 0.5) ? Pzz1 : Pzz2;
            f_temp = max_f * r5;

            G4double h = x_counter_e_temp * cos(theta0_counter) + R0_counter * sin(theta0_counter);// высота точки на счетчике (координата x) в лабораторной системе отсчета
            //   printf ("h = %f \n", h);

            //   p6 = zz_cell + R0_counter / cos(theta0_counter);                       // расстояние от точки на мишени до точки пересечения счетчика с осью Z

            G4double p6 = -zz_cell_temp + R0_counter / cos(theta0_counter);                       // расстояние от точки на мишени до точки пересечения счетчика с осью Z

            G4double p1 = sqrt(pow(p6 - h * tan(theta0_counter), 2.) + pow(h, 2.));         // расстояние от точки мишени до точки на счетчике в плоскость xz
            G4double p2 = sqrt(pow(p1, 2.) + pow(y_counter_e_temp, 2.));                           //радиус из точки мишени на точку dxdy счетчика
            G4double R_ds = p2;
            G4double p3 = zz_cell_temp * cos(theta0_counter) + R0_counter;

            G4double p4 = sqrt(pow(y_counter_e_temp, 2.) + pow(x_counter_e_temp - zz_cell_temp * sin(theta0_counter), 2.));
            G4double p5 = sqrt(pow(y_counter_e_temp, 2.) + pow(x_counter_e_temp + R0_counter * tan(theta0_counter), 2.));

            G4double cos_theta = (pow(p2, 2.) + pow(p6, 2.) - pow(p5, 2.)) / (2. * p2 * p6);// по теореме косинусов угол между векторами
            G4double cos_ds = (pow(p2, 2.) + pow(p3, 2.) - pow(p4, 2.)) / (2. * p2 * p3);

            //    dom_method_2 = d_x_counter * d_y_counter * (1./cos_ds )/ pow(R_ds, 2.);

            //G4double dom_method_2 = d_zz_cell * d_x_counter * d_y_counter * cos_ds / pow(R_ds, 2.);
            // G4double dom_without_z_cell_distr_method_2 =  d_x_counter * d_y_counter * cos_ds / pow(R_ds, 2.);

            //    dom_method_2 = d_x_counter * d_y_counter / pow(R_ds, 2.);


            //  printf ("dx = %e \n", d_x_counter);
            //  printf ("dy = %e \n", d_y_counter);
            //  printf ("dom = %e \n", dom);

            theta_e_temp = acos(cos_theta);
            //  if (y_counter>0) phi_e = atan(fabs(y_counter) / h);
            //  if (y_counter<0) phi_e = -atan(fabs(y_counter) / h);
            // printf ("h = %f y_counter = %f theta_e = %f \n", h,y_counter, phi_e/D);

            ///// внимание, именно так
            phi_e_temp = atan(y_counter_e_temp / h);

            f_dsdo= dsdo_LQ(theta_e_temp, Pzz_temp, zz_cell_temp);
            //   G4cerr << "!!! dsdo= " << f_dsdo <<" rnd = " << f_temp <<std::endl;
            if (f_temp <= f_dsdo) {
                flaq = true;
            } else {
                flaq = false;
            }
        }

        // Returns
        theta_electron = theta_e_temp;
        G4double rc = 1. + 2. * Ebeam * pow(sin(theta_electron / 2.), 2.) / Md; // фактор отдачи
        energy_deuteron = (1. - 1. / rc) * Ebeam; // кин энергия дейтрона
        energy_electron = Ebeam / rc; // энергия электрона после рассеяния
        theta_deuteron = atan(sin(theta_electron) / (rc - cos(theta_electron)));
        phi_electron = (G4UniformRand() > 0.5) ? phi_e_temp : phi_e_temp +
                                                              M_PI; // возможность, чтобы электроны летели в нижний счетчик тоже
        phi_deuteron = M_PI - phi_electron;
        xx_cell = xx_cell_temp;
        yy_cell = yy_cell_temp;
        zz_cell = zz_cell_temp;
        Pzz = Pzz_temp;

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

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


    void PrimaryGeneratorAction::GenerateLowQ_method1(G4Event *event) {

        G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();

        G4double initial_theta_e = 30. * M_PI / 180., final_theta_e = 90. * M_PI / 180.;
        G4double initial_zz_cell = -l_zz_cell / 2., final_zz_cell = l_zz_cell / 2.; // в cm
        G4double max_f = 1.0;
        G4double momentum, kinetic_energy;
        G4double Pzz1 = 1.0, r = -2.0;
        G4double Pzz2 = r * Pzz1;
        G4double Pzz =0.0;

        G4double theta_electron = 0.0, theta_deuteron = 0.0;
        G4double phi_electron = 0.0, phi_deuteron = 0.0;
        G4double energy_electron = 0.0, energy_deuteron = 0.0;
        G4double xx_cell = 0.0, yy_cell = 0.0, zz_cell = 0.0;

        random_Neumann_LQ_method1(initial_theta_e, final_theta_e, initial_zz_cell,
                          final_zz_cell, Pzz1, Pzz2, max_f,
                          theta_electron, phi_electron, energy_electron,
                          theta_deuteron, phi_deuteron, energy_deuteron,
                          Pzz, xx_cell, yy_cell, zz_cell);

     //   G4cerr << "!!! electron " << theta_electron <<"   " << phi_electron << "   "<< energy_electron <<std::endl;
     //   G4cerr << "!!! deuteron " << theta_deuteron <<"   " << phi_deuteron << "   "<< energy_deuteron <<std::endl;
     //   G4cerr << "!!! cell " << xx_cell <<"   " << yy_cell << "   "<< zz_cell <<std::endl;
     //   G4cerr <<std::endl;

        //  theta_electron = 50. * M_PI / 180.;
        //  theta_deuteron = 20. * M_PI / 180.;
        //  phi_electron = 50. * M_PI / 180.;
        //  phi_deuteron = 20. * M_PI / 180.;
        //  energy_electron = 60;
        //  energy_deuteron = 500;

        p4vector electron, deuteron, p0;
        G4double md = 1875.63, ra = 3.14159 / 180.;

        p0.e() = md + Ebeam;
        p0.z() = Ebeam;
        // electron.e()=Ebeam-ed;
        electron.e() = energy_electron;
        electron.x() = electron.e() * sin(theta_electron) * sin(phi_electron);
        electron.y() = electron.e() * sin(theta_electron) * cos(phi_electron);
        electron.z() = electron.e() * cos(theta_electron);
        electron.theta() = electron.thetar();
        electron.phi() = electron.phir();

        // deuteron=p0-e;
        deuteron.e() = energy_deuteron;
        deuteron.x() = deuteron.e() * sin(theta_deuteron) * sin(phi_deuteron);
        deuteron.y() = deuteron.e() * sin(theta_deuteron) * cos(phi_deuteron);
        deuteron.z() = deuteron.e() * cos(theta_deuteron);
        deuteron.theta() = deuteron.thetar();
        deuteron.phi() = deuteron.phir();


        G4ParticleDefinition *particle; //for electron
        particle = particleTable->FindParticle("e-");
        fParticleGun->SetParticleDefinition(particle);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(electron.x(), electron.y(), electron.z()));
        // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta_electron)*cos(phi_electron),sin(theta_electron)*sin(phi_electron),cos(theta_electron)));
        fParticleGun->SetParticleEnergy(energy_electron * MeV);
        fParticleGun->SetParticlePosition(G4ThreeVector(xx_cell * cm, yy_cell * cm, zz_cell * cm));
        fParticleGun->GeneratePrimaryVertex(event);


        particle = particleTable->FindParticle("deuteron"); //for deuteron
        fParticleGun->SetParticleDefinition(particle);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(deuteron.x(), deuteron.y(), deuteron.z()));
        // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta_deuteron)*cos(phi_deuteron),sin(theta_deuteron)*sin(phi_deuteron),cos(theta_deuteron)));

        fParticleGun->SetParticleEnergy(energy_deuteron * MeV);
        fParticleGun->SetParticlePosition(G4ThreeVector(xx_cell * cm, yy_cell * cm, zz_cell * cm));
        fParticleGun->GeneratePrimaryVertex(event);

        EventInfo* info = new EventInfo();
//   EventInfo* info =(EventInfo*)anEvent->GetUserInformation();
        //info->SetEgamma(Egamma);
        // info->SetNreac(nreac);
        //info->SetNp(np);
        //info->SetEntry(FileNum);
        info->SetPzz(Pzz);
        event->SetUserInformation(info);

    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void PrimaryGeneratorAction::GenerateLowQ_method2(G4Event *event) {

        G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();


        G4double x_counter_initial = -l_theta_counter, x_counter_final = l_theta_counter;
        G4double y_counter_initial = -l_phi_counter, y_counter_final = l_phi_counter;

        G4double initial_zz_cell = -l_zz_cell / 2., final_zz_cell = l_zz_cell / 2.; // в cm
        G4double max_f = 1.0;
        G4double momentum, kinetic_energy;
        G4double Pzz = 0.0;

        G4double theta_electron = 0.0, theta_deuteron = 0.0;
        G4double phi_electron = 0.0, phi_deuteron = 0.0;
        G4double energy_electron = 0.0, energy_deuteron = 0.0;
        G4double xx_cell = 0.0, yy_cell = 0.0, zz_cell = 0.0;

        random_Neumann_LQ_method2(x_counter_initial, x_counter_final,
                                  y_counter_initial, y_counter_final,
                                  initial_zz_cell, final_zz_cell, Pzz1, Pzz2, max_f,
                                  theta_electron, phi_electron, energy_electron,
                                  theta_deuteron, phi_deuteron, energy_deuteron,
                                  Pzz, xx_cell, yy_cell, zz_cell);

        //   G4cerr << "!!! electron " << theta_electron <<"   " << phi_electron << "   "<< energy_electron <<std::endl;
        //   G4cerr << "!!! deuteron " << theta_deuteron <<"   " << phi_deuteron << "   "<< energy_deuteron <<std::endl;
        //   G4cerr << "!!! cell " << xx_cell <<"   " << yy_cell << "   "<< zz_cell <<std::endl;
        //   G4cerr <<std::endl;

        //  theta_electron = 50. * M_PI / 180.;
        //  theta_deuteron = 20. * M_PI / 180.;
        //  phi_electron = 50. * M_PI / 180.;
        //  phi_deuteron = 20. * M_PI / 180.;
        //  energy_electron = 60;
        //  energy_deuteron = 500;

        p4vector electron, deuteron, p0;
        G4double md = 1875.63, ra = 3.14159 / 180.;

        p0.e() = md + Ebeam;
        p0.z() = Ebeam;
        // electron.e()=Ebeam-ed;
        electron.e() = energy_electron;
        electron.x() = electron.e() * sin(theta_electron) * sin(phi_electron);
        electron.y() = electron.e() * sin(theta_electron) * cos(phi_electron);
        electron.z() = electron.e() * cos(theta_electron);
        electron.theta() = electron.thetar();
        electron.phi() = electron.phir();

        // deuteron=p0-e;
        deuteron.e() = energy_deuteron;
        deuteron.x() = deuteron.e() * sin(theta_deuteron) * sin(phi_deuteron);
        deuteron.y() = deuteron.e() * sin(theta_deuteron) * cos(phi_deuteron);
        deuteron.z() = deuteron.e() * cos(theta_deuteron);
        deuteron.theta() = deuteron.thetar();
        deuteron.phi() = deuteron.phir();


        G4ParticleDefinition *particle; //for electron
        particle = particleTable->FindParticle("e-");
        fParticleGun->SetParticleDefinition(particle);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(electron.x(), electron.y(), electron.z()));
        // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta_electron)*cos(phi_electron),sin(theta_electron)*sin(phi_electron),cos(theta_electron)));
        fParticleGun->SetParticleEnergy(energy_electron * MeV);
        fParticleGun->SetParticlePosition(G4ThreeVector(xx_cell * cm, yy_cell * cm, zz_cell * cm));
        fParticleGun->GeneratePrimaryVertex(event);


        particle = particleTable->FindParticle("deuteron"); //for deuteron
        fParticleGun->SetParticleDefinition(particle);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(deuteron.x(), deuteron.y(), deuteron.z()));
        // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta_deuteron)*cos(phi_deuteron),sin(theta_deuteron)*sin(phi_deuteron),cos(theta_deuteron)));

        fParticleGun->SetParticleEnergy(energy_deuteron * MeV);
        fParticleGun->SetParticlePosition(G4ThreeVector(xx_cell * cm, yy_cell * cm, zz_cell * cm));
        fParticleGun->GeneratePrimaryVertex(event);

          // G4cerr << "!!! energy " <<energy_electron <<"   " << energy_deuteron <<std::endl;
          // G4cerr <<std::endl;

        EventInfo* info = new EventInfo();
//   EventInfo* info =(EventInfo*)anEvent->GetUserInformation();
        //info->SetEgamma(Egamma);
       // info->SetNreac(nreac);
        //info->SetNp(np);
        //info->SetEntry(FileNum);
        info->SetPzz(Pzz);
        event->SetUserInformation(info);

    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void PrimaryGeneratorAction::GenerateProton(G4Event *event) {
        G4ThreeVector vertex;
        G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
        G4ParticleDefinition *particle;
        G4ParticleMomentum Momentum;
        G4double energy1, theta1, phi1;
        G4double energy2, theta2, phi2;

        G4double Xbeam = 0., Ybeam = 0.;
        vertex = GenVertex(Xbeam, Ybeam, 0.7*mm, 0.3*mm);

        particle = particleTable->FindParticle("proton");

        energy1 = (10.+450.*G4UniformRand())*MeV;
        theta1 = M_PI*G4UniformRand();
        phi1 = 2.0*M_PI*G4UniformRand();

        // Нужно генерить две частицы для срабатывания тригера
        energy2 = (10.+450.*G4UniformRand())*MeV;
        theta2 = M_PI*G4UniformRand();
        phi2 =M_PI-phi1;

        auto v_direction1 = G4ThreeVector (sin(theta1) * sin(phi1), sin(theta1) * cos(phi1), cos(theta1));
        auto v_direction2 = G4ThreeVector (sin(theta2) * sin(phi2), sin(theta2) * cos(phi2), cos(theta2));

        fParticleGun->SetParticleEnergy(energy1);
        //fParticleGun->SetParticleMomentumDirection(v_direction1);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(2.*G4UniformRand()-1,1.,G4UniformRand()));
        fParticleGun->SetParticleDefinition(particle);
        fParticleGun->SetParticlePosition(vertex);
        fParticleGun->GeneratePrimaryVertex(event);

        fParticleGun->SetParticleEnergy(energy2);
      //  fParticleGun->SetParticleMomentumDirection(v_direction2);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(2.*G4UniformRand()-1,-1.,G4UniformRand()));
        fParticleGun->SetParticleDefinition(particle);
        fParticleGun->SetParticlePosition(vertex);
        fParticleGun->GeneratePrimaryVertex(event);

        EventInfo* info = new EventInfo();
//   pn2020EventInfo* info =(pn2020EventInfo*)anEvent->GetUserInformation();
        info->SetEgamma(500.0);
        info->SetNreac(1);
        info->SetNp(2);
        info->SetEntry(FileNum);
        event->SetUserInformation(info);

    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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


        EventInfo* info = new EventInfo();
//   EventInfo* info =(EventInfo*)anEvent->GetUserInformation();
        info->SetEgamma(Egamma);
        info->SetNreac(nreac);
        info->SetNp(np);
        info->SetEntry(FileNum);
        info->SetPzz((G4UniformRand() > 0.5) ? Pzz1 : Pzz2);
        event->SetUserInformation(info);


    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void PrimaryGeneratorAction::DoRandomize() {
        // Привязка генератора к системному времени (через устройство)
        // Этот метод работает только в Линукс

        if (rndmFlag == "on") {
// Randomizer ;-)
            G4AutoLock lock(&aMutex);
            G4int i;
            long prand;

            i = open("/dev/urandom", O_RDONLY);
            if (i < 0) prand = time(NULL);
            else {
                read(i, &prand, sizeof(long));
                close(i);
            }
            prand &= 0xFFFFFF;
            i = (G4int) prand;
            CLHEP::HepRandom::setTheSeed(prand);
            G4cout << "\n/\\/\\/\\ Randomizied !  prand = " << prand << " /\\/\\/\\" << G4endl;;

            genbos_rand_(&i);
// ---------------------
        }

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

    void PrimaryGeneratorAction::SetMode(G4int val) {
        int nreac, ireac[40];
        Mode = val;
        G4AutoLock lock(&aMutex);
        if (Mode == 34) { // это p+n
            memset(ireac, 0, sizeof(ireac));
            nreac = 1;
            ireac[0] = 34;
           // G4AutoLock lock(&aMutex);
            genbos_reactions_(&nreac, ireac);
            genbos_start_(&FileNum);
            
        }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
}
