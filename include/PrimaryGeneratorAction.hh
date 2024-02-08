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
#include "PrimaryGeneratorMessenger.hh"
#include "globals.hh"
#include <Randomize.hh>
#include "G4ThreeVector.hh"
#include "Constants.hh"
#include "p4vector.hh"

//#include "GenbosClass.hh"

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

        void GenerateCosmic(G4Event *event); // for cosmic-generator
        void GenerateLowQ_ed_method1(G4Event* event); // for ed-generator (elastic, polarized) for LQ-polarimeter
        void GenerateLowQ_ed_method2(G4Event* event);
        void GenerateLowQ_ep_method2(G4Event* event); // for ep-generator (elastic, unpolarized) for LQ-polarimeter
        void GenerateLowQ_ep_quasi_elastic_method2(G4Event* event);
        // for ep-generator (quasi-elastic, unpolarized) for LQ-polarimeter


        void GenerateProton(G4Event *event); // for pp-generator
        void GenerateNeutron(G4Event *event); // for nn-generator
        void GenerateProtonNeutron(G4Event *event); // for pn-generator
        void GenerateProtonNeutron_rachek(G4Event *event); // for pn-generator
        void GenerateDeuteronPi0(G4Event *event); // for dpi0-generator
        void GenerateGamma(G4Event *event); // for gamma-generator




        G4double virtual_photon_random();

        G4double virtual_photon_distrib(G4double);

        // set methods
        void SetRandomFlag(G4bool value);

        ///// FOR GENBOS
    public:
        void GenerateGenbos(G4Event *event); // for pn-generator GENBOS
        inline G4double GetEgamma() { return Egamma; };

        void SetRndmFlag(G4String val) { rndmFlag = val; }

        void SetCountFlag(G4String val) { countFlag = val; }

        void SetVertexFlag(G4String val) { vertexFlag = val; }

        G4bool GetCountFlag() { return (countFlag == "on"); };

        G4bool GetVertexFlag() { return (vertexFlag == "on"); };

        void SetCStep(G4int val) { cstep = val; }

        G4int GetCStep() { return cstep; }

        void DoRandomize();

        void SetMode(G4int val);

        void SetFileNum(G4int val) { FileNum = val; }

        void SetGenbosBool(G4int val) { GenbosBool = val;
        //    genbos_start_(&FileNum);
         //   PrepareNames();
        }

        void SetEgMin(G4double val) {
            EgMin = val / GeV;
            G4int n = 2;
#ifdef GENBOS
            genbos_beam_(&n, &EgMin, &EgMax);
#endif //GENBOS
        }

        void SetEgMax(G4double val) {
            EgMax = val / GeV;
            G4int n = 2;
#ifdef GENBOS
            genbos_beam_(&n, &EgMin, &EgMax);
#endif //GENBOS


        }

        void SetEgMinMax(G4double val_min, G4double val_max) {
            EgMin = val_min / GeV;
            EgMax = val_max / GeV;
            G4int n = 2;

#ifdef GENBOS
            genbos_beam_(&n, &EgMin, &EgMax);
#endif //GENBOS

        }
        //////


    private:
        G4ParticleGun *fParticleGun = nullptr; // G4 particle gun
        PrimaryGeneratorMessenger *gunMessenger; //messenger of this class
        void ShowParticleTable();

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


        p4vector prodol(G4double p0, G4double w0, G4double m0, p4vector p)
        {
            G4double pe, px, py, pz;
            px=p.x();
            py=p.y();
            pz=w0*p.z()/m0+p0*p.e()/m0;
            pe=w0*p.e()/m0+p0*p.z()/m0;
            return p4vector(pe, px, py, pz);
        }

        ///// FOR GENBOS

    private:
        G4ThreeVector
        GenVertex(G4double, G4double, G4double, G4double); //  метод генерации точки вылета из пучка по гауссу

        G4int GenbosBool;

        // GenbosClass *fGenbosClass = nullptr;

        G4int cstep;
        G4String countFlag;
        G4String rndmFlag;
        G4String vertexFlag;
        G4int Mode;
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
        G4int *idg = new G4int[11];   //[np] // list of particle's indexes
        G4float *cx = new G4float[11];   //[np] // lists of momentum components
        G4float *cy = new G4float[11];   //[np]
        G4float *cz = new G4float[11];   //[np]



        //////

        //////// FOR LQ-generator


        // формфакторы дейтрона изйна  таблицы


        G4double Gq_deuteron_spline(G4double Q2) {
            G4double Q2_temp =
                    Q2 / pow(10., 6.); // в функцию передается квадрат переданного импульса в МэВ, переводиться в ГэВ
            G4double a = -2.10393;
            G4double b = 1.40109;
            G4double c = 0.0462036;
            G4double temp = a + b / (c + Q2_temp); // сплайн гиперболой
            return temp;
        }


        G4double Ge_deuteron_spline(G4double Q2) {
            G4double Q2_temp =
                    Q2 / pow(10., 6.); // в функцию передается квадрат переданного импульса в МэВ, переводиться в ГэВ
            G4double a = -0.112158;
            G4double b = 0.0587308;
            G4double c = 0.0443795;
            G4double temp = a + b / (c + Q2_temp); // сплайн гиперболой
            return temp;
        }

        G4double Gm_deuteron_spline(G4double Q2) {
            G4double Q2_temp =
                    Q2 / pow(10., 6.); // в функцию передается квадрат переданного импульса в МэВ, переводиться в ГэВ
            G4double a = -0.132418;
            G4double b = 0.0829012;
            G4double c = -0.00188002;
            G4double temp = a + b / Q2_temp + c / pow(Q2_temp, 2.); // сплайн гиперболой
            return temp;
        }

        G4double assym_calc(G4double a, G4double b, G4double r = -2.) {
            // a - число событий с pol=1 (P_zz=+1), b - число событий с pol=2 (P_zz=-2)
            G4double assym_temp = (a - b) / (b - r * a);
            //  cout<< pf1t[1]<<endl;
            return assym_temp;
        }

        G4double assym_err_calc(G4double a, G4double b, G4double r = -2.) {

            //  float err_temp = 3. * sqrt(a * b * (a + b)) / pow(2. * a + b,2.); // Получено по квадратичному сложению ошибок
            G4double err_temp = fabs(r - 1.) * sqrt(a * b * (a + b)) / pow(b - r * a, 2.);
            return err_temp;
        }



        G4double dsdo_LQ(G4double theta_e, G4double Pzz, G4double zz_cell) {
            G4int zp0 = T2M_zanulenie[0], zp1 = T2M_zanulenie[1], zp2 = T2M_zanulenie[2];
            // зануление отдельных компонент
            G4double coeff_target =1.;

            if (zz_cell >= 0) coeff_target = (l_zz_cell / 2. - zz_cell) / pow(l_zz_cell / 2., 2.);// коэффициент для учета труегольного распределения частиц по мишени
            if (zz_cell < 0) coeff_target = (l_zz_cell / 2. + zz_cell) / pow(l_zz_cell / 2., 2.);

            G4double rc = 1. + (2. * Ebeam / Md) * pow(sin(theta_e / 2.), 2.); // ?????? ??????
            G4double Q2 = 4. * pow(Ebeam, 2.) * pow(sin(theta_e / 2.), 2.) / rc; // ???????? ????????
            G4double Qf2 = Q2 / pow(Fm, 2.);

            G4double energy_e = Ebeam / rc;

            G4double nu = Q2 / (4. * Md * Md);
            G4double epsilon = (1. + nu) * pow(tan(theta_e / 2.), 2.);

            G4double Gc = Ge_deuteron_spline(Q2); // Формфакторы
            G4double Gm = Gm_deuteron_spline(Q2);
            G4double Gq = Gq_deuteron_spline(Q2);

            G4double A = pow(Gc, 2.0) + (8. / 9.) * nu * nu * pow(Gq, 2.0) +
                         (2. / 3.) * nu * pow(Gm, 2.0); // Структурные функции
            G4double B = (4. / 3.) * nu * (1. + nu) * pow(Gm, 2.);
            G4double S = A + B * pow(tan(theta_e / 2.), 2.);

            G4double Mt = (pow(alpha_em * Fm, 2.) * 10000. / 4.) * pow(cos(theta_e / 2.), 2.) * S / (rc * pow(sin(theta_e / 2.), 4.) * pow(Ebeam, 2.));
            // Mott cross-Section , mkb


            G4double T20 = -(4. / 3.) * sqrt(2.) * (nu / S) * (Gc * Gq + nu / 3. * pow(Gq, 2.) +
                                                               (.5 + epsilon) * pow(Gm, 2.) /
                                                               4.);// ????????-????????????? ???????????
            G4double T21 = -(Gq * Gm / S) * 2. * nu * sqrt(nu * (1. + epsilon) / 3.);
            G4double T22 = -pow(Gm, 2.) * nu / (2. * sqrt(3.) * S);

            G4double q3md = sqrt(pow(Ebeam, 2.) + pow(energy_e, 2.) - 2. * Ebeam * energy_e * cos(theta_e));        // вектор передачи от начального электрона к конечному
            G4double cos_theta_otdachi = (pow(q3md, 2.) + pow(Ebeam, 2.) - pow(energy_e, 2.)) / (2. * q3md * Ebeam);// косинус угла между векторами по теореме косинусов
            G4double theta_otdachi = acos(cos_theta_otdachi);                                                 // угол дейтрона отдачи для формулы поляризованного сечения



            G4double d20 = (3. * pow(cos(theta_otdachi), 2.) - 1.) / 2.; // это правильно
            G4double d21 = -sqrt(3. / 2.) * sin(2. * theta_otdachi);
            G4double d22 = sqrt(3. / 2.) * pow(sin(theta_otdachi), 2.);

            G4double dsdo_temp = sin(theta_e) * (1. + (Pzz / sqrt(2)) * (d20 * T20 * zp0 + d21 * T21 * zp1 + d22 * T22 * zp2)) * Mt * coeff_target;

            return dsdo_temp;

        }

        void random_Neumann_LQ_method1(G4double initial_theta_e, G4double final_theta_e,
                               G4double initial_zz_cell, G4double final_zz_cell,
                               G4double Pzz1, G4double Pzz2, G4double max_f,
                               G4double &theta_electron, G4double &phi_electron, G4double &energy_electron,
                               G4double &theta_deuteron, G4double &phi_deuteron, G4double &energy_deuteron,
                               G4double &Pzz,
                               G4double &xx_cell, G4double &yy_cell, G4double &zz_cell);

        void random_Neumann_LQ_method2(G4double initial_x_counter_e, G4double final_x_counter_e,
                                       G4double initial_y_counter_e, G4double final_y_counter_e,
                                       G4double initial_zz_cell, G4double final_zz_cell,
                                       G4double Pzz1, G4double Pzz2, G4double max_f,
                                       G4double &theta_electron, G4double &phi_electron, G4double &energy_electron,
                                       G4double &theta_deuteron, G4double &phi_deuteron, G4double &energy_deuteron,
                                       G4double &Pzz,
                                       G4double &xx_cell, G4double &yy_cell, G4double &zz_cell);

        void random_Neumann_LQ_ep_method2(G4double initial_x_counter_e, G4double final_x_counter_e,
                                       G4double initial_y_counter_e, G4double final_y_counter_e,
                                       G4double initial_zz_cell, G4double final_zz_cell,
                                       G4double Pzz1, G4double Pzz2, G4double max_f,
                                       G4double &theta_electron, G4double &phi_electron, G4double &energy_electron,
                                       G4double &theta_deuteron, G4double &phi_deuteron, G4double &energy_deuteron,
                                       G4double &Pzz,
                                       G4double &xx_cell, G4double &yy_cell, G4double &zz_cell);

        void random_Neumann_LQ_ep_quasi_elastic_method2(G4double initial_x_counter_e, G4double final_x_counter_e,
                                       G4double initial_y_counter_e, G4double final_y_counter_e,
                                       G4double initial_zz_cell, G4double final_zz_cell,
                                       G4double Pzz1, G4double Pzz2, G4double max_f,
                                       G4double &theta_electron, G4double &phi_electron, G4double &energy_electron,
                                       G4double &theta_deuteron, G4double &phi_deuteron, G4double &energy_deuteron,
                                       G4double &Pzz,
                                       G4double &xx_cell, G4double &yy_cell, G4double &zz_cell);


    };

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
