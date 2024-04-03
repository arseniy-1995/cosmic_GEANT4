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


    namespace {
        G4Mutex aMutex = G4MUTEX_INITIALIZER;

    }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    PrimaryGeneratorAction::PrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(), fParticleGun(0),
                                                       GenbosBool(0), cstep(100), countFlag("off"), rndmFlag("off"),
                                                       vertexFlag("off"), Mode(0),
                                                       FileNum(0),
            // FileNum(G4Threading::G4GetThreadId()),
                                                       EgMin(400 * MeV), EgMax(650 * MeV) {

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


        //G4AutoLock lock(&aMutex);

        if (GenbosBool == 1) {
            G4AutoLock lock(&aMutex);



           // lock.lock();
#ifdef GENBOS
           // lock.lock();
            genbos_start_(&FileNum);

           // SetMode(34);
#endif //GENBOS

            G4int n = 2;
#ifdef GENBOS
            genbos_beam_(&n, &EgMin, &EgMax);

#endif //GENBOS

            // lock.unlock();
            PrepareNames();

              //fGenbosClass = new GenbosClass(&FileNum);

            // G4Threading::G4GetThreadId();
        }
        //  genbos_start_(&FileNum);

        // PrepareNames();
        // genbos_start_(&FileNum);
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    PrimaryGeneratorAction::~PrimaryGeneratorAction() {
        if (GenbosBool == 1) {
            G4AutoLock lock(&aMutex);
        }

#ifdef GENBOS
        genbos_stop_();
#endif
        // lock.unlock();

        // delete fGenbosClass;
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
#ifdef isGenCosmic
            GenerateCosmic(anEvent);
#endif

#ifdef isGenLQ
            //GenerateLowQ_ed_method1(anEvent);
            GenerateLowQ_ed_method2(anEvent); ///!!!!

#endif

#ifdef isGenPNpair
            GenerateProtonNeutron(anEvent);
            // GenerateProtonNeutron_rachek(anEvent);
#endif

#ifdef isGenPPpair
             GenerateProton(anEvent);
#endif
            //  GenerateNeutron(anEvent);

            // GenerateDeuteronPi0(anEvent);
           // GenerateGamma(anEvent);
           //  GenerateGenbos(anEvent);
        } else if (GenbosBool == 1) {
            //genbos_start_(&FileNum);
            G4AutoLock lock(&aMutex);
            //  G4cout<<"!!!"<< std::endl;

           // lock.unlock();
            //genbos_start_(&FileNum);
            //  fGenbosClass = new GenbosClass(&FileNum);

            //genbos_start_(&FileNum);

            GenerateGenbos(anEvent);
            // delete fGenbosClass;
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


    void PrimaryGeneratorAction::ShowParticleTable()
    {

        // #define SW(N)    " "<<std::setw( N )
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
                                                           G4double &theta_electron, G4double &phi_electron,
                                                           G4double &energy_electron,
                                                           G4double &theta_deuteron, G4double &phi_deuteron,
                                                           G4double &energy_deuteron,
                                                           G4double &Pzz,
                                                           G4double &xx_cell, G4double &yy_cell, G4double &zz_cell) {

        G4double theta_e_temp = 0.0;
        G4double phi_e_temp = 0.0;
        G4double xx_cell_temp = 0.0, yy_cell_temp = 0.0, zz_cell_temp = 0.0;
        G4double Pzz_temp = 0.0;
        G4double f_temp = 0.0;
        G4double r1 = 0., r2 = 0., r3 = 0., r4 = 0., r5 = 0.;


        G4double f_dsdo = 0.0;
        G4bool flaq = false;

        while (flaq == false) {
            r1 = G4UniformRand();
            r2 = G4UniformRand();
            r3 = G4UniformRand();
            r4 = G4UniformRand();
            r5 = G4UniformRand();

            theta_e_temp = initial_theta_e + r1 * (final_theta_e - initial_theta_e);
            xx_cell_temp = 0.1 * G4RandGauss::shoot(meanX_beam, Xsigma_beam) / mm;
            yy_cell_temp = 0.1 * G4RandGauss::shoot(meanY_beam, Ysigma_beam) / mm;
            zz_cell_temp = initial_zz_cell + r2 * (final_zz_cell - initial_zz_cell);
            Pzz_temp = (r3 > 0.5) ? Pzz1 : Pzz2;
            phi_e_temp = 2.0 * M_PI * r4;
            f_temp = max_f * r5;
            f_dsdo = dsdo_LQ(theta_e_temp, Pzz_temp, zz_cell_temp);
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
                                                           G4double &theta_electron, G4double &phi_electron,
                                                           G4double &energy_electron,
                                                           G4double &theta_deuteron, G4double &phi_deuteron,
                                                           G4double &energy_deuteron,
                                                           G4double &Pzz,
                                                           G4double &xx_cell, G4double &yy_cell, G4double &zz_cell) {

        G4double theta_e_temp = 0.0, phi_e_temp = 0.0;
        G4double x_counter_e_temp = 0.0;
        G4double y_counter_e_temp = 0.0;
        G4double xx_cell_temp = 0.0, yy_cell_temp = 0.0, zz_cell_temp = 0.0;
        G4double Pzz_temp = 0.0;
        G4double f_temp = 0.0;
        G4double r1 = 0., r2 = 0., r3 = 0., r4 = 0., r5 = 0.;


        G4double f_dsdo = 0.0;
        G4bool flaq = false;

        while (flaq == false) {

            r1 = G4UniformRand();
            r2 = G4UniformRand();
            r3 = G4UniformRand();
            r4 = G4UniformRand();
            r5 = G4UniformRand();

            x_counter_e_temp = initial_x_counter_e + r1 * (final_x_counter_e - initial_x_counter_e);
            y_counter_e_temp = initial_x_counter_e + r2 * (final_y_counter_e - initial_y_counter_e);
            xx_cell_temp = 0.1 * G4RandGauss::shoot(meanX_beam, Xsigma_beam) / mm; // в см
            yy_cell_temp = 0.1 * G4RandGauss::shoot(meanY_beam, Ysigma_beam) / mm;
            zz_cell_temp = initial_zz_cell + r3 * (final_zz_cell - initial_zz_cell);
            Pzz_temp = (r4 > 0.5) ? Pzz1 : Pzz2;
            f_temp = max_f * r5;

            G4double h = x_counter_e_temp * cos(theta0_counter) + R0_counter *
                                                                  sin(theta0_counter);// высота точки на счетчике (координата x) в лабораторной системе отсчета
            //   printf ("h = %f \n", h);

            //   p6 = zz_cell + R0_counter / cos(theta0_counter);                       // расстояние от точки на мишени до точки пересечения счетчика с осью Z

            G4double p6 = -zz_cell_temp + R0_counter /
                                          cos(theta0_counter);                       // расстояние от точки на мишени до точки пересечения счетчика с осью Z

            G4double p1 = sqrt(pow(p6 - h * tan(theta0_counter), 2.) +
                               pow(h, 2.));         // расстояние от точки мишени до точки на счетчике в плоскость xz
            G4double p2 = sqrt(pow(p1, 2.) + pow(y_counter_e_temp,
                                                 2.));                           //радиус из точки мишени на точку dxdy счетчика
            G4double R_ds = p2;
            G4double p3 = zz_cell_temp * cos(theta0_counter) + R0_counter;

            G4double p4 = sqrt(
                    pow(y_counter_e_temp, 2.) + pow(x_counter_e_temp - zz_cell_temp * sin(theta0_counter), 2.));
            G4double p5 = sqrt(
                    pow(y_counter_e_temp, 2.) + pow(x_counter_e_temp + R0_counter * tan(theta0_counter), 2.));

            G4double cos_theta = (pow(p2, 2.) + pow(p6, 2.) - pow(p5, 2.)) /
                                 (2. * p2 * p6);// по теореме косинусов угол между векторами
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

            f_dsdo = dsdo_LQ(theta_e_temp, Pzz_temp, zz_cell_temp);
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
        phi_electron = (G4UniformRand() > 0.5)
                           ? phi_e_temp
                           : phi_e_temp +
                           M_PI; // возможность, чтобы электроны летели в нижний счетчик тоже

        // phi_electron = phi_e_temp + M_PI;

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

            G4double x_initial_pos, y_initial_pos = 12.0 * m + 4.5 * m, z_initial_pos;
            G4double x_initial = -1.0 * m, y_initial = 12.0 * m + 4.5 * m, z_initial = -0.8 * m + 0.8 * m;
            G4double x_final = 1.0 * m, y_final = 12.0 * m + 4.5 * m, z_final = 0.8 * m + 0.8 * m;
            G4double vx, vy, vz;
            G4double theta, psi;
            G4double initial_theta = 0. * M_PI / 180., final_theta = 10. * M_PI / 180.;
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

        auto particleTable = G4ParticleTable::GetParticleTable();

        auto particleDefinition = particleTable->FindParticle("mu-");
        fParticleGun->SetParticleDefinition(particleDefinition);


        fParticleGun->SetParticleMomentumDirection(initial_momentum_direction);
        fParticleGun->SetParticlePosition(initial_pos);
        fParticleGun->GeneratePrimaryVertex(event);
        // fParticleGun->GeneratePrimaryVertex(event); // две частицы из одной точки

        EventInfo *info = new EventInfo();
//   pn2020EventInfo* info =(pn2020EventInfo*)anEvent->GetUserInformation();
        info->SetEgamma(500. * MeV);
        info->SetNreac(1);
        info->SetNp(2);
        info->SetEntry(FileNum);
        event->SetUserInformation(info);

    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// метод 1 - это просто по углам, без учета ячейки
void PrimaryGeneratorAction::GenerateLowQ_ed_method1(G4Event* event)
{
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

    G4double initial_theta_e = 30. * M_PI / 180., final_theta_e = 90. * M_PI / 180.;
    G4double initial_zz_cell = -l_zz_cell / 2., final_zz_cell = l_zz_cell / 2.; // в cm
    G4double max_f = 1.0;
    G4double momentum, kinetic_energy;
    G4double Pzz1 = 1.0, r = -2.0;
        G4double Pzz2 = r * Pzz1;
        G4double Pzz = 0.0;

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

        EventInfo *info = new EventInfo();
//   EventInfo* info =(EventInfo*)anEvent->GetUserInformation();
        //info->SetEgamma(Egamma);
        // info->SetNreac(nreac);
        //info->SetNp(np);
        //info->SetEntry(FileNum);
        info->SetPzz(Pzz);
        event->SetUserInformation(info);

    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// метод 2 - это по поверхности счетчика, с учетам ячейки
// упрогое ed-рассеяние
void PrimaryGeneratorAction::GenerateLowQ_ed_method2(G4Event* event)
{
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();


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


    // G4cout << "theta = " << deuteron.thetar() * 180./ M_PI << " " << deuteron.theta() * 180./ M_PI << " " << theta_deuteron * 180./ M_PI << G4endl;
    // G4cout << "phi = " << deuteron.phir() * 180./ M_PI << " " << deuteron.phi() * 180./ M_PI << " " << phi_deuteron * 180./ M_PI << G4endl;

    G4ParticleDefinition* particle; //for electron
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

    EventInfo *info = new EventInfo();
//   EventInfo* info =(EventInfo*)anEvent->GetUserInformation();
        //info->SetEgamma(Egamma);
        // info->SetNreac(nreac);
        //info->SetNp(np);
        //info->SetEntry(FileNum);
        info->SetPzz(Pzz);
        event->SetUserInformation(info);

    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// упрогое ep-рассеяние
void PrimaryGeneratorAction::GenerateLowQ_ep_method2(G4Event* event)
{
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();


    G4double x_counter_initial = -l_theta_counter, x_counter_final = l_theta_counter;
    G4double y_counter_initial = -l_phi_counter, y_counter_final = l_phi_counter;

    G4double initial_zz_cell = -l_zz_cell / 2., final_zz_cell = l_zz_cell / 2.; // в cm
    G4double max_f = 1.0;
    G4double momentum, kinetic_energy;
    G4double Pzz = 0.0;

    G4double theta_electron = 0.0, theta_proton = 0.0;
    G4double phi_electron = 0.0, phi_proton = 0.0;
    G4double energy_electron = 0.0, energy_proton = 0.0;
    G4double xx_cell = 0.0, yy_cell = 0.0, zz_cell = 0.0;

    // TODO
    // random_Neumann_LQ_ep_method2(x_counter_initial, x_counter_final,
    //                              y_counter_initial, y_counter_final,
    //                              initial_zz_cell, final_zz_cell, Pzz1, Pzz2, max_f,
    //                              theta_electron, phi_electron, energy_electron,
    //                              theta_proton, phi_proton, energy_proton,
    //                              Pzz, xx_cell, yy_cell, zz_cell);

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

    p4vector electron, proton, p0;
    G4double mp = 938.2720881629, ra = 3.14159 / 180.;

    p0.e() = mp + Ebeam;
    p0.z() = Ebeam;
    // electron.e()=Ebeam-ed;
    electron.e() = energy_electron;
        electron.x() = electron.e() * sin(theta_electron) * sin(phi_electron);
        electron.y() = electron.e() * sin(theta_electron) * cos(phi_electron);
        electron.z() = electron.e() * cos(theta_electron);
        electron.theta() = electron.thetar();
        electron.phi() = electron.phir();

        // proton=p0-e;
        proton.e() = energy_proton;
        proton.x() = proton.e() * sin(theta_proton) * sin(phi_proton);
        proton.y() = proton.e() * sin(theta_proton) * cos(phi_proton);
        proton.z() = proton.e() * cos(theta_proton);
        proton.theta() = proton.thetar();
        proton.phi() = proton.phir();


        G4ParticleDefinition* particle; //for electron
        particle = particleTable->FindParticle("e-");
        fParticleGun->SetParticleDefinition(particle);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(electron.x(), electron.y(), electron.z()));
        // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta_electron)*cos(phi_electron),sin(theta_electron)*sin(phi_electron),cos(theta_electron)));
        fParticleGun->SetParticleEnergy(energy_electron * MeV);
        fParticleGun->SetParticlePosition(G4ThreeVector(xx_cell * cm, yy_cell * cm, zz_cell * cm));
        fParticleGun->GeneratePrimaryVertex(event);


        particle = particleTable->FindParticle("proton"); //for proton
        fParticleGun->SetParticleDefinition(particle);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(proton.x(), proton.y(), proton.z()));
        // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta_proton)*cos(phi_proton),sin(theta_proton)*sin(phi_proton),cos(theta_proton)));

        fParticleGun->SetParticleEnergy(energy_proton * MeV);
        fParticleGun->SetParticlePosition(G4ThreeVector(xx_cell * cm, yy_cell * cm, zz_cell * cm));
        fParticleGun->GeneratePrimaryVertex(event);

        // G4cerr << "!!! energy " <<energy_electron <<"   " << energy_proton <<std::endl;
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

    // квази-упрогое ep-рассеяние
    void PrimaryGeneratorAction::GenerateLowQ_ep_quasi_elastic_method2(G4Event* event)
    {
        G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();


        G4double x_counter_initial = -l_theta_counter, x_counter_final = l_theta_counter;
        G4double y_counter_initial = -l_phi_counter, y_counter_final = l_phi_counter;

        G4double initial_zz_cell = -l_zz_cell / 2., final_zz_cell = l_zz_cell / 2.; // в cm
        G4double max_f = 1.0;
        G4double momentum, kinetic_energy;
        G4double Pzz = 0.0;

        G4double theta_electron = 0.0, theta_proton = 0.0;
        G4double phi_electron = 0.0, phi_proton = 0.0;
        G4double energy_electron = 0.0, energy_proton = 0.0;
        G4double xx_cell = 0.0, yy_cell = 0.0, zz_cell = 0.0;


        ///

        // это интеграл плотности вероятности по импульсу с шагом 10 МэВ?
        // см. fermi_d.F из GENBOS

        G4double distr_paris_potential[101] = {
            0.0, 6.7953602E-03, 4.6238571E-02, 0.1236217, 0.2234159,
            0.3280089, 0.4262296, 0.5131546, 0.5876556, 0.6504297,
            0.7028700, 0.7465120, 0.7827951, 0.8129796, 0.8381307,
            0.8591338, 0.8767187, 0.8914840, 0.9039201, 0.9144295,
            0.9233416, 0.9309281, 0.9374120, 0.9429782, 0.9477784,
            0.9519395, 0.9555655, 0.9587440, 0.9615465, 0.9640333,
            0.9662549, 0.9682527, 0.9700612, 0.9717091, 0.9732212,
            0.9746166, 0.9759121, 0.9771215, 0.9782565, 0.9793263,
            0.9803385, 0.9812997, 0.9822153, 0.9830891, 0.9839256,
            0.9847271, 0.9854963, 0.9862349, 0.9869449, 0.9876275,
            0.9882837, 0.9889148, 0.9895214, 0.9901043, 0.9906641,
            0.9912014, 0.9917168, 0.9922106, 0.9926836, 0.9931361,
            0.9935685, 0.9939811, 0.9943746, 0.9947497, 0.9951068,
            0.9954458, 0.9957677, 0.9960729, 0.9963619, 0.9966350,
            0.9968930, 0.9971365, 0.9973655, 0.9975811, 0.9977834,
            0.9979731, 0.9981509, 0.9983172, 0.9984723, 0.9986169,
            0.9987513, 0.9988762, 0.9989921, 0.9990993, 0.9991982,
            0.9992898, 0.9993736, 0.9994510, 0.9995219, 0.9995865,
            0.9996458, 0.9996991, 0.9997482, 0.9997926, 0.9998326,
            0.9998688, 0.9999010, 0.9999300, 0.9999564, 0.9999794,
            1.000000
        };

        G4double distr_bonn_potential[101] = {
            0.0, 7.0403279E-03, 4.7886062E-02, 0.1279474, 0.2310553,
            0.3389300, 0.4400139, 0.5292494, 0.6055161, 0.6695801,
            0.7229201, 0.7671533, 0.8037899, 0.8341466, 0.8593346,
            0.8802744, 0.8977222, 0.9122964, 0.9245030, 0.9347540,
            0.9433876, 0.9506804, 0.9568596, 0.9621121, 0.9665918,
            0.9704258, 0.9737193, 0.9765596, 0.9790186, 0.9811566,
            0.9830235, 0.9846612, 0.9861036, 0.9873803, 0.9885156,
            0.9895294, 0.9904380, 0.9912565, 0.9919961, 0.9926671,
            0.9932781, 0.9938356, 0.9943457, 0.9948140, 0.9952440,
            0.9956399, 0.9960047, 0.9963413, 0.9966521, 0.9969389,
            0.9972037, 0.9974484, 0.9976741, 0.9978822, 0.9980739,
            0.9982503, 0.9984124, 0.9985616, 0.9986980, 0.9988228,
            0.9989371, 0.9990408, 0.9991356, 0.9992215, 0.9992995,
            0.9993696, 0.9994332, 0.9994900, 0.9995410, 0.9995866,
            0.9996274, 0.9996639, 0.9996958, 0.9997242, 0.9997495,
            0.9997714, 0.9997910, 0.9998078, 0.9998228, 0.9998358,
            0.9998475, 0.9998578, 0.9998671, 0.9998750, 0.9998825,
            0.9998896, 0.9998960, 0.9999021, 0.9999081, 0.9999138,
            0.9999200, 0.9999259, 0.9999323, 0.9999395, 0.9999465,
            0.9999539, 0.9999623, 0.9999706, 0.9999802, 0.9999896,
            1.000000
        };


        // TODO
        // random_Neumann_LQ_ep_quasi_elastic_method2(x_counter_initial, x_counter_final,
        //                                            y_counter_initial, y_counter_final,
        //                                            initial_zz_cell, final_zz_cell, Pzz1, Pzz2, max_f,
        //                                            theta_electron, phi_electron, energy_electron,
        //                                            theta_proton, phi_proton, energy_proton,
        //                                            Pzz, xx_cell, yy_cell, zz_cell);

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

        p4vector electron, proton, p0;

        G4double mp = 938.2720881629, ra = 3.14159 / 180.;

        p0.e() = mp + Ebeam;
        p0.z() = Ebeam;
        // electron.e()=Ebeam-ed;
        electron.e() = energy_electron;
        electron.x() = electron.e() * sin(theta_electron) * sin(phi_electron);
        electron.y() = electron.e() * sin(theta_electron) * cos(phi_electron);
        electron.z() = electron.e() * cos(theta_electron);
        electron.theta() = electron.thetar();
        electron.phi() = electron.phir();

        // deuteron=p0-e;
        proton.e() = energy_proton;
        proton.x() = proton.e() * sin(theta_proton) * sin(phi_proton);
        proton.y() = proton.e() * sin(theta_proton) * cos(phi_proton);
        proton.z() = proton.e() * cos(theta_proton);
        proton.theta() = proton.thetar();
        proton.phi() = proton.phir();


        G4ParticleDefinition* particle; //for electron
        particle = particleTable->FindParticle("e-");
        fParticleGun->SetParticleDefinition(particle);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(electron.x(), electron.y(), electron.z()));
        // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta_electron)*cos(phi_electron),sin(theta_electron)*sin(phi_electron),cos(theta_electron)));
        fParticleGun->SetParticleEnergy(energy_electron * MeV);
        fParticleGun->SetParticlePosition(G4ThreeVector(xx_cell * cm, yy_cell * cm, zz_cell * cm));
        fParticleGun->GeneratePrimaryVertex(event);


        particle = particleTable->FindParticle("proton"); //for proton
        fParticleGun->SetParticleDefinition(particle);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(proton.x(), proton.y(), proton.z()));
        // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta_proton)*cos(phi_proton),sin(theta_proton)*sin(phi_proton),cos(theta_proton)));

        fParticleGun->SetParticleEnergy(energy_proton * MeV);
        fParticleGun->SetParticlePosition(G4ThreeVector(xx_cell * cm, yy_cell * cm, zz_cell * cm));
        fParticleGun->GeneratePrimaryVertex(event);

        // G4cerr << "!!! energy " <<energy_electron <<"   " << energy_proton <<std::endl;
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
        G4ParticleDefinition *particle1, *particle2;
        G4ParticleMomentum Momentum;
        G4double energy1, theta1, phi1;
        G4double energy2, theta2, phi2;

        G4double Xbeam = 0., Ybeam = 0.;
        vertex = GenVertex(Xbeam, Ybeam, Xsigma_beam, Ysigma_beam);

        particle1 = particleTable->FindParticle("proton");
        particle2 = particleTable->FindParticle("proton");

        const G4double ra = M_PI / 180.;

        G4double max_phi = 40. * ra;
        G4double min_theta = 35. * ra;
        G4double max_theta = 120. * ra;

        G4ThreeVector v_direction1(0,0,1);
        G4ThreeVector v_direction2(0,0,1);

        if(1)
        {
            do
            {
                energy1 = (1. + 459. * G4UniformRand()) * MeV; // от 1 до 460 МэВ
                //energy1 = 1000 * MeV;
                //theta1 = M_PI * G4UniformRand();
                theta1 = acos(1.0 - 2.0 * G4UniformRand());
                phi1 = 2.0 * M_PI * G4UniformRand();

                // Нужно генерить две частицы для срабатывания тригера
                energy2 = (1. + 459. * G4UniformRand()) * MeV;
                //energy2 = 1000 * MeV;
                // theta2 = M_PI * G4UniformRand();
                theta2 = acos(1.0 - 2.0 * G4UniformRand());
                //theta2= theta1;
                //theta2= M_PI - theta1;
                phi2 = M_PI + phi1;

                //v_direction1 = G4ThreeVector(sin(theta1) * cos(phi1), sin(theta1) * sin(phi1), cos(theta1));
                //v_direction2 = G4ThreeVector(sin(theta2) * cos(phi2), sin(theta2) * sin(phi2), cos(theta2));

                v_direction1.setTheta(theta1);
                v_direction1.setPhi(phi1);


                v_direction2.setTheta(theta2);
                v_direction2.setPhi(phi2);

                v_direction1 = G4ThreeVector(2. * G4UniformRand() - 1, 1., G4UniformRand());
                v_direction2 = G4ThreeVector(2. * G4UniformRand() - 1, -1., G4UniformRand());

            }while (

            ((v_direction1.theta() < min_theta || v_direction1.theta() > max_theta) ||
                            (v_direction2.theta() < min_theta || v_direction2.theta() > max_theta))
                            ||
                not(
                   ((abs(v_direction1.phi() - (-90.) * ra) < max_phi) && v_direction1.phi() < 0. && (abs(v_direction2.phi() - (90.) * ra) < max_phi) && v_direction2.phi() > 0.) // протон вниз - нейтрон вверх
                   or
                   ((abs(v_direction1.phi() - (90.) * ra) < max_phi) && v_direction1.phi() > 0. && (abs(v_direction2.phi() - (-90.) * ra) < max_phi) && v_direction2.phi() < 0.) // протон вверх - нейтрон вниз
                )
                );
        }

        fParticleGun->SetParticleEnergy(energy1);
         fParticleGun->SetParticleMomentumDirection(v_direction1);
        //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(2. * G4UniformRand() - 1, 1., G4UniformRand()));
        fParticleGun->SetParticleDefinition(particle1);
        fParticleGun->SetParticlePosition(vertex);
        fParticleGun->GeneratePrimaryVertex(event);

        fParticleGun->SetParticleEnergy(energy2);
           fParticleGun->SetParticleMomentumDirection(v_direction2);
        //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(2. * G4UniformRand() - 1, -1., G4UniformRand()));
        fParticleGun->SetParticleDefinition(particle2);
        fParticleGun->SetParticlePosition(vertex);
        fParticleGun->GeneratePrimaryVertex(event);

        EventInfo *info = new EventInfo();
//   pn2020EventInfo* info =(pn2020EventInfo*)anEvent->GetUserInformation();
        info->SetEgamma(500.0);
        info->SetNreac(1);
        info->SetNp(2);
        info->SetEntry(FileNum);
        event->SetUserInformation(info);

    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void PrimaryGeneratorAction::GenerateNeutron(G4Event *event) {
        G4ThreeVector vertex;
        G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
        G4ParticleDefinition *particle;
        G4ParticleMomentum Momentum;
        G4double energy1, theta1, phi1;
        G4double energy2, theta2, phi2;

        G4double Xbeam = 0., Ybeam = 0.;
        vertex = GenVertex(Xbeam, Ybeam, Xsigma_beam, Ysigma_beam);

        particle = particleTable->FindParticle("neutron");

        energy1 = (10. + 450. * G4UniformRand()) * MeV;
        // theta1 = M_PI * G4UniformRand();
        theta1 = acos(1.0 - 2.0 * G4UniformRand());
        phi1 = 2.0 * M_PI * G4UniformRand();

        // Нужно генерить две частицы для срабатывания тригера
        energy2 = (10. + 450. * G4UniformRand()) * MeV;
        // theta2 = M_PI * G4UniformRand();
        theta2 = acos(1.0 - 2.0 * G4UniformRand());
        phi2 = M_PI - phi1;

        auto v_direction1 = G4ThreeVector(sin(theta1) * cos(phi1), sin(theta1) * sin(phi1), cos(theta1));
        auto v_direction2 = G4ThreeVector(sin(theta2) * cos(phi2), sin(theta2) * sin(phi2), cos(theta2));

        fParticleGun->SetParticleEnergy(energy1);
        //fParticleGun->SetParticleMomentumDirection(v_direction1);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(2. * G4UniformRand() - 1., 1., G4UniformRand()));
        fParticleGun->SetParticleDefinition(particle);
        fParticleGun->SetParticlePosition(vertex);
        fParticleGun->GeneratePrimaryVertex(event);

        fParticleGun->SetParticleEnergy(energy2);
        //  fParticleGun->SetParticleMomentumDirection(v_direction2);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(2. * G4UniformRand() - 1., -1., G4UniformRand()));
        fParticleGun->SetParticleDefinition(particle);
        fParticleGun->SetParticlePosition(vertex);
        fParticleGun->GeneratePrimaryVertex(event);

        EventInfo *info = new EventInfo();
//   pn2020EventInfo* info =(pn2020EventInfo*)anEvent->GetUserInformation();
        info->SetEgamma(500.0);
        info->SetNreac(1);
        info->SetNp(2);
        info->SetEntry(FileNum);
        event->SetUserInformation(info);

    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void PrimaryGeneratorAction::GenerateProtonNeutron(G4Event *event) {
        G4ThreeVector vertex;
        G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
        G4String particleName;
        G4ParticleDefinition *particle1, *particle2;
        G4ParticleMomentum Momentum;
        G4double energy1, theta1, phi1;
        G4double energy2, theta2, phi2;

        G4double Xbeam = 0., Ybeam = 0.;
        vertex = GenVertex(Xbeam, Ybeam, Xsigma_beam, Ysigma_beam);

        particle1 = particleTable->FindParticle("proton");
        particle2 = particleTable->FindParticle("neutron");

        // Нужно генерить две частицы для срабатывания тригера


        G4double x0 = 0.0, y0 = 0.0, z0 = 0.0;
        const G4double ra = M_PI / 180.;
        //  G4double costet = 0.0;
        G4double theta = 0.0, phi = 0.0, M2 = 0.0;
        G4int n_det = 0;

        p4vector q, d, n, p, pn;
        p4vector p_pn; // продольная составляющая


        G4double max_phi = 40. * ra;
        G4double min_theta = 35. * ra;
        G4double max_theta = 120. * ra;


        // G4double max_phi = 5. * ra;
        // G4double min_theta = 63. * ra;
        // G4double max_theta = 75. * ra;

        do
        {
            n_det++;
            n.mm() = Mn * Mn;
            p_pn.mm() = p.mm() = Mp * Mp;
            d.e() = Md; // дейтрон покоиться, развал
            d.mm() = Md * Md;
            q.e() = virtual_photon_random(); //generate energy virtual photon
            q.z() = q.e();
            pn = q + d;
            M2 = sqrt(pn.mmr());
            p_pn.e() = (M2 * M2 + Mp * Mp - Mn * Mn) / 2. / M2;
            // costet = 2. * rand() / 2147483647. - 1.;
            // phi=360.*ra*rand()/2147483647.;
            //  theta = M_PI * G4UniformRand();
            // theta = acos(2.*G4UniformRand()-1);
            theta = acos(1.0 - 2.0 * G4UniformRand());
            //phi = 2.0 * M_PI * G4UniformRand();
            phi = G4UniformRand() > 0.5
                      ? max_phi * (2. * G4UniformRand() - 1) + 90. * ra
                      : max_phi * (2. * G4UniformRand() - 1) - 90. * ra;
            // phi = max_phi * (2.*G4UniformRand()-1) + 90. * ra;
            p_pn.x() = sqrt((p_pn.e() * p_pn.e() - Mp * Mp)) * sin(theta) * cos(phi);
            p_pn.y() = sqrt((p_pn.e() * p_pn.e() - Mp * Mp)) * sin(theta) * sin(phi);
            p_pn.z() = sqrt(p_pn.e() * p_pn.e() - Mp * Mp) * cos(theta);
            p = prodol(pn.p(), pn.e(), M2, p_pn);

            p.theta() = p.thetar();
            p.phi() = p.phir();
            n = pn - p;
            n.theta() = n.thetar();
            n.phi() = n.phir();

        } while (
            (((p.e() - Mp) < 20. || (p.e() - Mp) > 1000.) ||
                ((n.e() - Mn) < 5. || (n.e() - Mn) > 1000.) ||
                // p.theta() < 30. * ra || p.theta() > 110. * ra ||
                // n.theta() < 30. * ra || n.theta() > 110. * ra
                (p.theta() < min_theta || p.theta() > max_theta) ||
                (n.theta() < min_theta || n.theta() > max_theta)) ||

            // (((abs(p.phi()-(-90.) * ra) > 45.*ra) && p.phi() > 0) || ((abs(n.phi()-(90.) * ra) > 45.*ra) && n.phi() < 0))  // протон вниз - нейтрон вверх
            //   || (((abs(p.phi()-(90.) * ra) > 45.*ra) && p.phi() < 0) || ((abs(n.phi()-(-90.) * ra) > 45.*ra) && n.phi() > 0)) // протон вверх - нейтрон вниз

            //(abs(p.phi()-(-90.) * ra) > 45.*ra)

            //(abs(p.phi()-(90.) * ra) > 45.*ra && abs(p.phi()-(-90.) * ra) > 45.*ra)
            // (p.phi() < (-90. - 35.) * ra && p.phi() > (-90. + 35.) * ra && p.phi() < 0 && n.phi() < (90. -  35.) * ra && n.phi() > (90. + 35.) * ra && n.phi() > 0)


            not(
               ((abs(p.phi() - (-90.) * ra) < max_phi) && p.phi() < 0. && (abs(n.phi() - (90.) * ra) < max_phi) && n
                   .phi() > 0.) // протон вниз - нейтрон вверх
               or
               ((abs(p.phi() - (90.) * ra) < max_phi) && p.phi() > 0. && (abs(n.phi() - (-90.) * ra) < max_phi) && n
                    .phi() < 0.) // протон вверх - нейтрон вниз
            )

        );
        //n.fi()<60.*ra||n.fi()>120.*ra );

        //G4cout<<"n_det="<<n_det<<"\t"<<"Eq="<<q.e()<<G4endl;
        // G4cout<<"Ep="<<p.e()-mp<<" tetp="<<p.tet()/ra<<" fip="<<p.fi()/ra<<G4endl;
        // G4cout<<"En="<<n.e()-mn<<" tetn="<<n.tet()/ra<<" fin="<<n.fi()/ra<<G4endl;
//cout<<q.e()+md<<" "<<p.e()+n.e()<<"\t\t"<<p.x()<<" "<<n.x()<<endl;


        // default particle kinematic

        particle1 = particleTable->FindParticle(particleName = "proton");
        fParticleGun->SetParticleDefinition(particle1);
        //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,-1.*G4UniformRand(),G4UniformRand()));
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(p.x(), p.y(), p.z()));
        //G4double ep = (50.+200.*G4UniformRand())*MeV;
        fParticleGun->SetParticleEnergy((p.e() - Mp) * MeV);
        //  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
        fParticleGun->SetParticlePosition(vertex);
        fParticleGun->GeneratePrimaryVertex(event);


        particle2 = particleTable->FindParticle(particleName = "neutron");
        fParticleGun->SetParticleDefinition(particle2);
        // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,G4UniformRand(),G4UniformRand()));
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(n.x(), n.y(), n.z()));
        //G4double en = (20.+200.*G4UniformRand())*MeV;
        fParticleGun->SetParticleEnergy((n.e() - Mn) * MeV);
        //  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
        fParticleGun->SetParticlePosition(vertex);
        fParticleGun->GeneratePrimaryVertex(event);


        // вручную точеная частица


        /*particle1 = particleTable->FindParticle(particleName = "proton");
       fParticleGun->SetParticleDefinition(particle1);
       fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., -1., 0.));
       fParticleGun->SetParticleEnergy(600. * MeV);
       fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
       fParticleGun->GeneratePrimaryVertex(event);


       particle2 = particleTable->FindParticle(particleName = "neutron");
       fParticleGun->SetParticleDefinition(particle2);
       fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 1., 0.));
       fParticleGun->SetParticleEnergy(600. * MeV);
       fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
       fParticleGun->GeneratePrimaryVertex(event);*/

        /*particle1 = particleTable->FindParticle(particleName = "mu-");
      fParticleGun->SetParticleDefinition(particle1);
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., -1., 0.));
      fParticleGun->SetParticleEnergy(1000. * MeV);
      fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.5*m));
      fParticleGun->GeneratePrimaryVertex(event);


      particle2 = particleTable->FindParticle(particleName = "mu-");
      fParticleGun->SetParticleDefinition(particle2);
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 1., 0.));
      fParticleGun->SetParticleEnergy(1000. * MeV);
      fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.5*m));
      fParticleGun->GeneratePrimaryVertex(event);*/


        // auto vertex2 =G4ThreeVector(1*(2*G4UniformRand()-1)*m, 0., 1*(2*G4UniformRand()-1)*m +0.7*m);
        //   particle1 = particleTable->FindParticle(particleName = "mu-");
        // fParticleGun->SetParticleDefinition(particle1);
        // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., -1., 0.));
        // fParticleGun->SetParticleEnergy(1000. * MeV);
        // fParticleGun->SetParticlePosition(vertex2 );
        // fParticleGun->GeneratePrimaryVertex(event);
        //
        //
        // particle2 = particleTable->FindParticle(particleName = "mu-");
        // fParticleGun->SetParticleDefinition(particle2);
        // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 1., 0.));
        // fParticleGun->SetParticleEnergy(1000. * MeV);
        // fParticleGun->SetParticlePosition(vertex2 );
        // fParticleGun->GeneratePrimaryVertex(event);


        EventInfo* info = new EventInfo();
        //   pn2020EventInfo* info =(pn2020EventInfo*)anEvent->GetUserInformation();
        info->SetEgamma(q.e());
        info->SetNreac(1);
        info->SetNp(2);
        info->SetEntry(FileNum);
        event->SetUserInformation(info);

    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void PrimaryGeneratorAction::GenerateProtonNeutron_rachek(G4Event *event) {
        G4ThreeVector vertex;
        G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
        G4String particleName;
        G4ParticleDefinition *d_particle, *p_particle, *n_particle;
        G4ParticleMomentum Momentum;
        G4double p_energy, p_theta, p_phi;
        G4double n_energy, n_theta, n_phi;
        G4double p_theta_cm; // угол протона в с.ц.м.

        G4double gamma_energy;

        G4double Xbeam = 0., Ybeam = 0.;
        vertex = GenVertex(Xbeam, Ybeam, Xsigma_beam, Ysigma_beam);

        // generate particles
        d_particle = particleTable->FindParticle("deuteron");
        G4double d_mass = d_particle->GetPDGMass();
        p_particle = particleTable->FindParticle("proton");
        G4double p_mass = p_particle->GetPDGMass();
        n_particle = particleTable->FindParticle("neutron");
        G4double n_mass = n_particle->GetPDGMass();

        // Нужно генерить две частицы для срабатывания тригера


        // default particle kinematic


        gamma_energy = virtual_photon_random(); //generate energy virtual photon
        // Range of Egamma and p_theta_CM  -  uniform distribution here

        // gamma_energy = (10.0 + 690.0 * G4UniformRand()) * MeV;    // from 10 to 700 MeV
        p_theta_cm = (20.0 + 110.0 * G4UniformRand()) * deg;    // from 20 to 140 degree

        // two-body photodesintegration kinematics


        G4LorentzVector deuteron4;
        G4LorentzVector proton4;
        G4LorentzVector neutron4;

        G4LorentzVector photon4(Egamma, G4ThreeVector(0.0, 0.0, gamma_energy));
        G4LorentzVector target(d_mass, G4ThreeVector());
        G4LorentzVector dg = photon4 + target;
        G4double dg_mass = dg.m();
        G4double beta = dg.beta();

        G4double pcm = sqrt((dg_mass * dg_mass - (p_mass + n_mass) * (p_mass + n_mass)) *
                            (dg_mass * dg_mass - (p_mass - n_mass) * (p_mass - n_mass))) /
                       (2.0 * dg_mass);
        G4double epcm = (dg_mass * dg_mass + p_mass * p_mass - n_mass * n_mass) /
                        (2.0 * dg_mass);

        proton4 = G4LorentzVector(epcm, G4ThreeVector(0.0, 0.0, pcm));

        proton4.setTheta(p_theta_cm);
        proton4.boostZ(beta);


        p_phi = 270.0 + (G4UniformRand() * 2.0 - 1.0) * 40.0;    // phi proton
        proton4.setPhi(p_phi * deg);
        p_energy = proton4.e() - p_mass;

        neutron4 = dg - proton4;
        n_energy = neutron4.e() - n_mass;


        G4ThreeVector p_direction(proton4.vect().unit());
        G4ThreeVector n_direction(neutron4.vect().unit());

        fParticleGun->SetParticleDefinition(p_particle);
        fParticleGun->SetParticleMomentumDirection(p_direction);
        fParticleGun->SetParticleEnergy(p_energy);
        fParticleGun->SetParticlePosition(vertex);
        fParticleGun->GeneratePrimaryVertex(event);


        fParticleGun->SetParticleDefinition(n_particle);
        fParticleGun->SetParticleMomentumDirection(n_direction);
        fParticleGun->SetParticleEnergy(n_energy);
        fParticleGun->SetParticlePosition(vertex);
        fParticleGun->GeneratePrimaryVertex(event);


        EventInfo *info = new EventInfo();
//   pn2020EventInfo* info =(pn2020EventInfo*)anEvent->GetUserInformation();
        info->SetEgamma(gamma_energy);
        info->SetNreac(1);
        info->SetNp(2);
        info->SetEntry(FileNum);
        event->SetUserInformation(info);

    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    void PrimaryGeneratorAction::GenerateDeuteronPi0(G4Event *event) {
        G4ThreeVector vertex;
        G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
        G4String particleName;
        G4ParticleDefinition *d_particle, *pi_particle;

        G4double pi_energy, d_energy;

        G4double d_theta, d_phi;


        G4double gamma_energy;

        G4double Xbeam = 0., Ybeam = 0.;
        vertex = GenVertex(Xbeam, Ybeam, Xsigma_beam, Ysigma_beam);

        // generate particles
        d_particle = particleTable->FindParticle("deuteron");
        G4double d_mass = d_particle->GetPDGMass();
        pi_particle = particleTable->FindParticle("pi0");
        G4double pi_mass = pi_particle->GetPDGMass();

        // Нужно генерить две частицы для срабатывания тригера


        // default particle kinematic


        //  gamma_energy = virtual_photon_random(); //generate energy virtual photon
        // Range of Egamma and p_theta_CM  -  uniform distribution here

        // gamma_energy = (10.0 + 690.0 * G4UniformRand()) * MeV;    // from 10 to 700 MeV
        gamma_energy = 400. * MeV; //generate energy virtual photon

        //G4double pi_theta_cm = (20.0 + 110.0 * G4UniformRand()) * deg;    // from 20 to 140 degree
        G4double pi_theta_cm = 120.0 * deg;


        G4LorentzVector photon4(Egamma, G4ThreeVector(0.0, 0.0, gamma_energy));
        G4LorentzVector target(d_mass, G4ThreeVector());
        G4LorentzVector dg = photon4 + target;
        G4double dg_mass = dg.m();
        G4double beta = dg.beta();


        G4double pcm = sqrt((dg_mass * dg_mass - (d_mass + pi_mass) * (d_mass + pi_mass)) *
                            (dg_mass * dg_mass - (d_mass - pi_mass) * (d_mass - pi_mass))) /
                       (2.0 * dg_mass);
        G4double edcm = (dg_mass * dg_mass + d_mass * d_mass - pi_mass * pi_mass) /
                        (2.0 * dg_mass);


        G4LorentzVector deuteron4(edcm, G4ThreeVector(0.0, 0.0, pcm));

        deuteron4.setTheta(M_PI - pi_theta_cm);
        deuteron4.boostZ(beta);
        d_theta = deuteron4.theta() * deg;

        d_phi = 270.0 + (G4UniformRand() * 2. - 1.) * 40.0;
        deuteron4.setPhi(d_phi * deg);
        d_energy = deuteron4.e() - d_mass;

        G4LorentzVector pion4 = dg - deuteron4;
        pi_energy = pion4.e() - pi_mass;

        G4ThreeVector pi_direction(pion4.vect().unit());
        G4ThreeVector d_direction(deuteron4.vect().unit());

        fParticleGun->SetParticleDefinition(d_particle);
        fParticleGun->SetParticleEnergy(d_energy);
        fParticleGun->SetParticlePosition(vertex);
        fParticleGun->SetParticleMomentumDirection(d_direction);
        fParticleGun->GeneratePrimaryVertex(event);

        fParticleGun->SetParticleDefinition(pi_particle);
        fParticleGun->SetParticlePosition(vertex);
        fParticleGun->SetParticleMomentumDirection(pi_direction);
        fParticleGun->SetParticleEnergy(pi_energy);
        fParticleGun->GeneratePrimaryVertex(event);


        EventInfo *info = new EventInfo();
//   pn2020EventInfo* info =(pn2020EventInfo*)anEvent->GetUserInformation();
        info->SetEgamma(gamma_energy);
        info->SetNreac(1);
        info->SetNp(2);
        info->SetEntry(FileNum);
        event->SetUserInformation(info);

    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void PrimaryGeneratorAction::GenerateGamma(G4Event *event) {
        G4ThreeVector vertex;
        G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
        G4String particleName;
        G4ParticleDefinition *gamma_particle;


        G4double gamma_energy1, gamma_theta1, gamma_phi1;
        G4double gamma_energy2, gamma_theta2, gamma_phi2;


        gamma_energy1 = (200. + 200. * G4UniformRand()) * MeV; // 200..400 MeV
        // theta1 = M_PI * G4UniformRand();

      //  gamma_theta1 = acos(1.0 - 2.0 * G4UniformRand());
     //   gamma_phi1 = 2.0 * M_PI * G4UniformRand();

        // Нужно генерить две частицы для срабатывания тригера
        gamma_energy2 = (200. + 200. * G4UniformRand()) * MeV; // 200..400 MeV
        // theta2 = M_PI * G4UniformRand();

      //  gamma_theta2 = acos(1.0 - 2.0 * G4UniformRand());
       // gamma_phi2 = M_PI - gamma_phi1;

      //  auto v_direction1 = G4ThreeVector(sin(gamma_theta1) * cos(gamma_phi1), sin(gamma_theta1) * sin(gamma_phi1),
      //                                    cos(gamma_theta1));
     //   auto v_direction2 = G4ThreeVector(sin(gamma_theta2) * cos(gamma_phi2), sin(gamma_theta2) * sin(gamma_phi2),
      //                                    cos(gamma_theta2));


        G4double Xbeam = 0., Ybeam = 0.;
        vertex = GenVertex(Xbeam, Ybeam, Xsigma_beam, Ysigma_beam);

        // generate particles
        gamma_particle = particleTable->FindParticle("gamma");


        // Нужно генерить две частицы для срабатывания тригера

        // default particle kinematic


        // gamma_energy = virtual_photon_random(); //generate energy virtual photon
        // Range of Egamma and p_theta_CM  -  uniform distribution here

        // gamma_energy = (10.0 + 690.0 * G4UniformRand()) * MeV;    // from 10 to 700 MeV
        //gamma_energy1 = 300. * MeV; //generate energy virtual photon
       // gamma_energy2 = 300. * MeV

      //  G4ThreeVector gamma_direction1(0, 1, 0);
       //  gamma_theta1 = 130.0;//+20.0*G4UniformRand();
      // gamma_theta1 = 30.0 +100.0*G4UniformRand(); // 30...130
      //  gamma_direction1.setTheta(gamma_theta1 * deg);
     //  gamma_phi1 = -30.0 + 60.0*G4UniformRand(); // -30...30
      //  gamma_direction1.setPhi(gamma_phi1 * deg);

       // G4ThreeVector gamma_direction2(0, -1, 0);
       // gamma_theta2 = 130.0;//+20.0*G4UniformRand();
       // gamma_theta2 = 30.0 + 100.0*G4UniformRand(); // 30...130
       // gamma_theta2 = 180. - gamma_theta1;
       // gamma_direction2.setTheta(gamma_theta2 * deg);
       // gamma_phi2 = 180. - gamma_phi1;
       // gamma_phi2 = -30.0 + 60.0*G4UniformRand(); // -30...30 вниз
       // gamma_direction2.setPhi(gamma_phi2 * deg);



//        auto v_direction1 = G4ThreeVector(sin(gamma_theta1) * cos(gamma_phi1), sin(gamma_theta1) * sin(gamma_phi1),
//                                            cos(gamma_theta1));
//           auto v_direction2 = G4ThreeVector(sin(gamma_theta2) * cos(gamma_phi2), sin(gamma_theta2) * sin(gamma_phi2),
//                                           cos(gamma_theta2));



        auto gamma_direction1 = G4ThreeVector(sin(gamma_theta1 * deg) * cos(gamma_phi1* deg), 1.,
                                          cos(gamma_theta1* deg));
        auto gamma_direction2 = G4ThreeVector(sin(gamma_theta2* deg) * cos(gamma_phi2* deg), -1.,
                                          cos(gamma_theta2* deg));

        fParticleGun->SetParticleDefinition(gamma_particle);
       //  fParticleGun->SetParticleMomentumDirection(gamma_direction1);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(2. * G4UniformRand() - 1., 1., G4UniformRand()));
        fParticleGun->SetParticleEnergy(gamma_energy1);
        fParticleGun->SetParticlePosition(vertex);
        fParticleGun->GeneratePrimaryVertex(event);


        fParticleGun->SetParticleDefinition(gamma_particle);
        // fParticleGun->SetParticleMomentumDirection(gamma_direction2);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(2. * G4UniformRand() - 1., -1., G4UniformRand()));
        fParticleGun->SetParticleEnergy(gamma_energy2);
        fParticleGun->SetParticlePosition(vertex);
        fParticleGun->GeneratePrimaryVertex(event);


        EventInfo *info = new EventInfo();
//   pn2020EventInfo* info =(pn2020EventInfo*)anEvent->GetUserInformation();
        info->SetEgamma(0.5 * (gamma_energy1 + gamma_energy2));
        info->SetNreac(1);
        info->SetNp(2);
        info->SetEntry(FileNum);
        event->SetUserInformation(info);

    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......





    void PrimaryGeneratorAction::GenerateGenbos(G4Event *event) {
        //  fGenbosClass = new GenbosClass(&FileNum);

        G4ThreeVector vertex;
        G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
        G4ParticleDefinition *particle;
        G4ParticleMomentum Momentum;

        G4int qq = 0;

        while (qq != 3) {

#ifdef GENBOS
            genbos_event_(&efot, &nreac, &np, idg, cx, cy, cz);
#endif //GENBOS


            /*

           //fGenbosClass->SetEgMinMax(400.*MeV,650.*MeV);

           fGenbosClass->GenbosGenerator();

           efot = fGenbosClass->GetEfot();
           nreac = fGenbosClass->GetNreac();
           np = fGenbosClass->GetNp();
           idg = fGenbosClass->GetId();
           cx = fGenbosClass->GetPx();
           cy = fGenbosClass->GetPy();
           cz = fGenbosClass->GetPz();
           */

            for (G4int i = 0; i < np; i++) {

                // G4cerr<<"cx[i] = " <<cx[i]<<" cy[i] = " <<cy[i]<< " cz[i] = " <<cz[i] <<std::endl;
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


        EventInfo *info = new EventInfo();
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
            G4int i = 0;
            long prand = 1;

#ifdef GENBOS
            i = open("/dev/urandom", O_RDONLY);
           // i = 5;
#endif //GENBOS


            if (i < 0) prand = time(NULL);
            else {
#ifdef GENBOS
                read(i, &prand, sizeof(long));
                close(i);
#endif //GENBOS

            }
            prand &= 0xFFFFFF;
            i = (G4int) prand;
            CLHEP::HepRandom::setTheSeed(prand);
            G4cout << "\n/\\/\\/\\ Randomizied !  prand = " << prand << " /\\/\\/\\" << G4endl;;
#ifdef GENBOS
            genbos_rand_(&i);
#endif //GENBOS


            //fGenbosClass->SetRandom(i);
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
        //   G4AutoLock lock(&aMutex);
        if (Mode == 34)
        {
            // это p+n
            memset(ireac, 0, sizeof(ireac));
            nreac = 1;
            ireac[0] = 34;
            //fGenbosClass->SetMode(Mode);
        }

        if (Mode == 100)
        {
            // это тип 3 и 8 и 18 только реакции с двумя протонами и пи-мезоном(и)
            memset(ireac, 0, sizeof(ireac));
            nreac = 3;
            ireac[0] = 3;
            ireac[1] = 8;
            ireac[1] = 18;
            //fGenbosClass->SetMode(Mode);
        }


         G4AutoLock lock(&aMutex);
#ifdef GENBOS
            genbos_reactions_(&nreac, ireac);
            genbos_start_(&FileNum);
#endif //GENBOS
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    G4double PrimaryGeneratorAction::virtual_photon_distrib(G4double Egamma) {
        G4double Eelectron = 2000.0;
        G4double fl1, fl2, f1, f2, sl;
        G4double x = Egamma / Eelectron;
        G4double Egamma_distr_temp;

        // Это код из Томска
        //  fl1 = 2. * Eelectron * (Eelectron - Egamma) / Me / (2. * Eelectron - Egamma);
        //  fl2 = (2. * Eelectron - Egamma) / Egamma;
        //  f1 = 2. * (1. - x + pow(x, 2.0) / 2.);
        //  f2 = pow(x, 2.0) / 2.;
        //  sl = x - 1.;
        //  Egamma_distr_temp = (alpha_em * (f1 * log10(fabs(fl1)) + f2 * log10(fabs(fl2)) + sl) / M_PI);

        //G4double theta_e_max = 5.*M_PI/180.;
        G4double theta_e_max = 84.0 * pow(10.0, -3); //84 mrad
        Eelectron = 800.0;
        Egamma_distr_temp = (2.0 * alpha_em / M_PI) *
                            ((1. - x + pow(x, 2.0) / 2.) * log(fabs((Eelectron - Egamma) * theta_e_max / (Me * x))) -
                             (1. - x));

        return Egamma_distr_temp;
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    G4double PrimaryGeneratorAction::virtual_photon_random() {
        G4double r1, r2, Egamma, y;

        G4double Egamma_min = 10., Egamma_max = 700.;
        do {
            // r1 = rand() / 2147483647.;
            // r2 = rand() / 2147483647.;
            r1 = G4UniformRand();
            r2 = G4UniformRand();
            //  Egamma = 145. * exp(r1 * log(10.)); // Энергия виртуального фотона в диапозоне 145...1450 МэВ
            //  y = 0.033 * r2;

            Egamma = Egamma_min * exp(r1 * log(Egamma_max /
                                               Egamma_min)); // Энергия виртуального фотона в диапозоне Egamma_min...Egamma_max
            y = 0.065 * r2; // максимум 0,065

        } while (y > virtual_photon_distrib(Egamma));
//cerr<<y<<"\t"<<dal(om)<<endl;
        return Egamma;
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
}
