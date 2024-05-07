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

#include <Randomize.hh>
#include "Constants.hh"
#include "G4ThreeVector.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "globals.hh"
#include "p4vector.hh"
// #include <TVector3.h>
#include <G4LorentzVector.hh>
#include <G4ThreeVector.hh>
#include <TFoam.h>
#include <TRandom3.h>
#include "Riostream.h"
#include "TFoamIntegrand.h"
#include "TMath.h"

// #include "GenbosClass.hh"

#ifdef GENBOS

#include "Genbos.hh"

#endif

using namespace TMath;
using namespace CLHEP;

class TFDISTR_LQ : public TFoamIntegrand
{
public:
    TFDISTR_LQ(){};
    Double_t Density(int nDim, Double_t *Xarg)
    {
        // Integrand for mFOAM

        // Double_t dsdo_LQ_ed_elastic_TFoam(Int_t nDim, Double_t *Xarg) {

        Double_t theta_e = Xarg[0];
        Double_t Pzz = Xarg[1];
        Double_t zz_cell = Xarg[2];

        Double_t dsdo_temp = 1;
        // dsdo_LQ_ed_elastic(theta_e, Pzz, zz_cell);

        return dsdo_temp;

        // }

        // return 1.;
    }
};


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
    ///


    class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
    {
    public:
        PrimaryGeneratorAction();

        ~PrimaryGeneratorAction() override;

        void GeneratePrimaries(G4Event *event) override;

        void GenerateCosmic(G4Event *event); // for cosmic-generator
        void GenerateLowQ_ed_method1(G4Event *event); // for ed-generator (elastic, polarized) for LQ-polarimeter
        void GenerateLowQ_ed_method2(G4Event *event);
        void GenerateLowQ_ep_method2(G4Event *event); // for ep-generator (elastic, unpolarized) for LQ-polarimeter
        void GenerateLowQ_ep_plus_ed_method2(G4Event *event); // for ep- (elastic or quasi-elastic, unpolarized) +
                                                              // ed-generator (elastic, polarized)  for LQ-polarimeter
        // for ep-generator (, unpolarized) for LQ-polarimeter


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

        void SetGenbosBool(G4int val)
        {
            GenbosBool = val;
            //    genbos_start_(&FileNum);
            //   PrepareNames();
        }

        void SetEgMin(G4double val)
        {
            EgMin = val / GeV;
            G4int n = 2;
#ifdef GENBOS

            // genbos_beam_(&n, &EgMin, &EgMax);

            G4float EgMean = (EgMax + EgMin) / 2.;
            G4float EgWidht = (EgMax - EgMin) / 2.;
            genbos_beam_(&n, &EgMean, &EgWidht);

#endif // GENBOS
        }

        void SetEgMax(G4double val)
        {
            EgMax = val / GeV;
            G4int n = 2;
#ifdef GENBOS

            // genbos_beam_(&n, &EgMin, &EgMax);

            G4float EgMean = (EgMax + EgMin) / 2.;
            G4float EgWidht = (EgMax - EgMin) / 2.;
            genbos_beam_(&n, &EgMean, &EgWidht);
#endif // GENBOS
        }

        void SetEgMinMax(G4double val_min, G4double val_max)
        {
            EgMin = val_min / GeV;
            EgMax = val_max / GeV;
            G4int n = 2;

#ifdef GENBOS

            // genbos_beam_(&n, &EgMin, &EgMax);

            G4float EgMean = (EgMax + EgMin) / 2.;
            G4float EgWidht = (EgMax - EgMin) / 2.;
            genbos_beam_(&n, &EgMean, &EgWidht);

#endif // GENBOS
        }
        //////


    private:
        G4ParticleGun *fParticleGun = nullptr; // G4 particle gun
        PrimaryGeneratorMessenger *gunMessenger; // messenger of this class
        void ShowParticleTable();

        G4bool fRandomDirection;

        G4double m_mu = 105.658374524; // в МэВ
        // другое распрределение по углу и энергии из работы Шебалина
        G4double f_theta_energy_sin(G4double theta, G4double momentum)
        {
            G4double temp = 0.0;
            G4double k = 2.26;
            G4double p0 = 500.0; // МэВ

            G4double a = 2.86;
            G4double b = 1.54;
            G4double sigma = (a - cos(theta)) / b;
            //  long double A = 1.0 / 0.158363; // константа интегрирования 0..pi/2
            G4double A = 1.0;
            // здесь без нормировки
            temp = A * pow(cos(theta), k) * exp(-pow(log(momentum / p0), 2.0) / (2.0 * pow(sigma, 2.0))) * sin(theta);
            return temp;
        }

        void random_Neumann_theta_energy_sin(G4double initial_theta, G4double final_theta, G4double initial_momentum,
                                             G4double final_momentum, G4double max_f, G4double &theta,
                                             G4double &momentum, G4double &kinetic_energy);


        p4vector prodol(G4double p0, G4double w0, G4double m0, p4vector p)
        {
            G4double pe, px, py, pz;
            px = p.x();
            py = p.y();
            pz = w0 * p.z() / m0 + p0 * p.e() / m0;
            pe = w0 * p.e() / m0 + p0 * p.z() / m0;
            return p4vector(pe, px, py, pz);
        }

        ///// FOR GENBOS

    private:
        G4ThreeVector GenVertex(G4double, G4double, G4double,
                                G4double); //  метод генерации точки вылета из пучка по гауссу

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

        G4float efot; // photon energy, GeV
        G4int nreac; // reaction index
        G4float vx;
        G4float vy;
        G4float vz;
        G4int np; // number of generated particles
        G4int *idg = new G4int[11]; //[np] // list of particle's indexes
        G4float *cx = new G4float[11]; //[np] // lists of momentum components
        G4float *cy = new G4float[11]; //[np]
        G4float *cz = new G4float[11]; //[np]


        TFDISTR_LQ *density_dsdo_LQ_ed_elastic_TFoam;
        TFoam *Foam;

        //////


        //////// FOR LQ-generator

        bool is_quasi_elastic_pd = true;

        // формфакторы дейтрона изйна  таблицы


        G4double Gq_deuteron_spline(G4double Q2)
        {
            G4double Q2_temp =
                Q2 / pow(10., 6.); // в функцию передается квадрат переданного импульса в МэВ, переводиться в ГэВ
            G4double a = -2.10393;
            G4double b = 1.40109;
            G4double c = 0.0462036;
            G4double temp = a + b / (c + Q2_temp); // сплайн гиперболой
            return temp;
        }


        G4double Ge_deuteron_spline(G4double Q2)
        {
            G4double Q2_temp =
                Q2 / pow(10., 6.); // в функцию передается квадрат переданного импульса в МэВ, переводиться в ГэВ
            G4double a = -0.112158;
            G4double b = 0.0587308;
            G4double c = 0.0443795;
            G4double temp = a + b / (c + Q2_temp); // сплайн гиперболой
            return temp;
        }

        G4double Gm_deuteron_spline(G4double Q2)
        {
            G4double Q2_temp =
                Q2 / pow(10., 6.); // в функцию передается квадрат переданного импульса в МэВ, переводиться в ГэВ
            G4double a = -0.132418;
            G4double b = 0.0829012;
            G4double c = -0.00188002;
            G4double temp = a + b / Q2_temp + c / pow(Q2_temp, 2.); // сплайн гиперболой
            return temp;
        }


        /////////


        // формфакторы протона из таблицы
        G4double Gdipole(G4double Q2)
        {
            G4double Q2_temp =
                Q2 / pow(10., 6.); // в функцию передается квадрат переданного импульса в МэВ, переводиться в ГэВ

            G4double M2_dipole = 0.71; // GeV^2
            G4double temp = 1. / pow(1. + Q2_temp / M2_dipole, 2.0); // сплайн гиперболой
            return temp;
        }


        G4double Ge_proton_spline(G4double Q2)
        {

            G4double temp = Gdipole(Q2); // сплайн гиперболой
            return temp;
        }

        G4double Gm_proton_spline(G4double Q2)
        {

            G4double temp = mu_proton * Ge_proton_spline(Q2); // сплайн гиперболой
            return temp;
        }


        //-------------------------
        // GEp from Kelly Parameterization
        G4double Ge_proton_spline_kelly(G4double Q2)
        {

           // Q2 = Q2 / 1e6; // в функцию передается квадрат переданного импульса в МэВ, переводиться в ГэВ

            G4double gep_par[4] = {-.299, 11.11, 14.11, 15.7};
            G4double a1 = gep_par[0];
            G4double b1 = gep_par[1];
            G4double b2 = gep_par[2];
            G4double b3 = gep_par[3];

            G4double tau = Q2 / (4. * Mp * Mp);

            G4double GE = (1. + a1 * tau) / (1. + b1 * tau + b2 * pow(tau, 2.) + b3 * pow(tau, 3.));
            return GE;
        }

        //-------------------------
        // GMp from Kelly Parameterization
        G4double Gm_proton_spline_kelly(G4double Q2)
        {

           // Q2 = Q2 / 1e6; // в функцию передается квадрат переданного импульса в МэВ, переводиться в ГэВ

            G4double gmp_par[4] = {.081, 11.15, 18.45, 5.31};
            G4double a1 = gmp_par[0];
            G4double b1 = gmp_par[1];
            G4double b2 = gmp_par[2];
            G4double b3 = gmp_par[3];

            G4double tau = Q2 / (4. * Mp * Mp);

            G4double GM = mu_proton * (1. + a1 * tau) / (1. + b1 * tau + b2 * pow(tau, 2.) + b3 * pow(tau, 3.));
            return GM;
        }

        G4double dsdo_LQ_ed_elastic(G4double theta_e, G4double Pzz, G4double zz_cell)
        {

            G4int zp0 = T2M_zanulenie[0], zp1 = T2M_zanulenie[1], zp2 = T2M_zanulenie[2];
            // зануление отдельных компонент
            G4double coeff_target = 1.;

            if (zz_cell >= 0)
                coeff_target = (l_zz_cell / 2. - zz_cell) /
                    pow(l_zz_cell / 2., 2.); // коэффициент для учета труегольного распределения частиц по мишени
            if (zz_cell < 0)
                coeff_target = (l_zz_cell / 2. + zz_cell) / pow(l_zz_cell / 2., 2.);

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

            G4double gamma_electron = Ebeam / Me;
            G4double beta_electron = sqrt(1. - 1. / pow(gamma_electron, 2.0));

            // pow(cos(theta_e / 2.), 2.)

            G4double Mt = (pow(alpha_em * Fm, 2.) * 10000. / 4.) *
                (1. - pow(beta_electron, 2.) * pow(sin(theta_e / 2.), 2.)) * S /
                (rc * pow(sin(theta_e / 2.), 4.) * pow(Ebeam, 2.));
            // Mott cross-Section , mkb

            // G4cout << "deuteron Q2 [Gev^2] = " << Q2 / 1e6  << " S = " << S << G4endl;

            G4double T20 = -(4. / 3.) * sqrt(2.) * (nu / S) *
                (Gc * Gq + nu / 3. * pow(Gq, 2.) +
                 (.5 + epsilon) * pow(Gm, 2.) / 4.); // ????????-????????????? ???????????
            G4double T21 = -(Gq * Gm / S) * 2. * nu * sqrt(nu * (1. + epsilon) / 3.);
            G4double T22 = -pow(Gm, 2.) * nu / (2. * sqrt(3.) * S);

            G4double q3md =
                sqrt(pow(Ebeam, 2.) + pow(energy_e, 2.) -
                     2. * Ebeam * energy_e * cos(theta_e)); // вектор передачи от начального электрона к конечному
            G4double cos_theta_otdachi = (pow(q3md, 2.) + pow(Ebeam, 2.) - pow(energy_e, 2.)) /
                (2. * q3md * Ebeam); // косинус угла между векторами по теореме косинусов
            G4double theta_otdachi =
                acos(cos_theta_otdachi); // угол дейтрона отдачи для формулы поляризованного сечения


            G4double d20 = (3. * pow(cos(theta_otdachi), 2.) - 1.) / 2.; // это правильно
            G4double d21 = -sqrt(3. / 2.) * sin(2. * theta_otdachi);
            G4double d22 = sqrt(3. / 2.) * pow(sin(theta_otdachi), 2.);

            G4double dsdo_temp = sin(theta_e) *
                (1. + (Pzz / sqrt(2)) * (d20 * T20 * zp0 + d21 * T21 * zp1 + d22 * T22 * zp2)) * Mt * coeff_target;

            if (dsdo_temp > 1.0)
            {
                // G4cout << "dsdo_temp = " << dsdo_temp << G4endl;
                return -1;
            }

            return dsdo_temp;
        }


        G4double dsdo_LQ_ep_elastic(G4double theta_e, G4double zz_cell)
        {


            G4double Mass = 0.;

            // Mass = Md;
            Mass = Mp;

            G4double coeff_target = 1.;

            if (zz_cell >= 0)
                coeff_target = (l_zz_cell / 2. - zz_cell) /
                    pow(l_zz_cell / 2., 2.); // коэффициент для учета труегольного распределения частиц по мишени
            if (zz_cell < 0)
                coeff_target = (l_zz_cell / 2. + zz_cell) / pow(l_zz_cell / 2., 2.);

            G4double rc = 1. + (2. * Ebeam / Mass) * pow(sin(theta_e / 2.), 2.); // фактор отдачи
            // G4double Q2 = 4. * pow(Ebeam, 2.) * pow(sin(theta_e / 2.), 2.) / rc; // ???????? ????????


            G4double energy_e = Ebeam / rc; // после рассеяния

            G4double Q2 = 4. * Ebeam * energy_e * pow(sin(theta_e / 2.), 2.);
            G4double Qf2 = Q2 / pow(Fm, 2.);
            G4double nu = Q2 / (4. * Mass * Mass);
            G4double epsilon = (1. + nu) * pow(tan(theta_e / 2.), 2.);

            // G4double Ge = Ge_proton_spline(Q2); // Формфакторы
            // G4double Gm = Gm_proton_spline(Q2);


            G4double Ge = Ge_proton_spline_kelly(Q2); // Формфакторы
            G4double Gm = Gm_proton_spline_kelly(Q2);

            G4double A = (nu * pow(Gm, 2.0) + pow(Ge, 2.0)) / (1.0 + nu); // Структурные функции

            G4double B = 2.0 * nu * pow(Gm, 2.);

            G4double S = A + B * pow(tan(theta_e / 2.), 2.);

            G4double gamma_electron = Ebeam / Me;
            G4double beta_electron = sqrt(1. - 1. / pow(gamma_electron, 2.0));

            // pow(cos(theta_e / 2.), 2.)

            G4double Mt = (pow(alpha_em * Fm, 2.) * 10000. / 4.) *
                (1. - pow(beta_electron, 2.) * pow(sin(theta_e / 2.), 2.)) * S /
                (rc * pow(sin(theta_e / 2.), 4.) * pow(Ebeam, 2.));
            // Mott cross-Section , mkb


            // G4cout << "Mt = " << Mt << G4endl;

            //  G4cout << "Q2 [Gev^2] = " << Q2 / 1e6 << " Ge = " << Ge  << " Gm = " << Gm << " S = " << S << G4endl;

            //   G4cout << "proton Q2 [Gev^2] = " << Q2 / 1e6  << " S = " << S << G4endl;

            G4double q3md =
                sqrt(pow(Ebeam, 2.) + pow(energy_e, 2.) -
                     2. * Ebeam * energy_e * cos(theta_e)); // вектор передачи от начального электрона к конечному
            G4double cos_theta_otdachi = (pow(q3md, 2.) + pow(Ebeam, 2.) - pow(energy_e, 2.)) /
                (2. * q3md * Ebeam); // косинус угла между векторами по теореме косинусов
            G4double theta_otdachi =
                acos(cos_theta_otdachi); // угол дейтрона отдачи для формулы поляризованного сечения


            G4double dsdo_temp = sin(theta_e) * Mt * coeff_target;

            if (dsdo_temp > 10.0)
            {
                // G4cout << "dsdo_temp = " << dsdo_temp << G4endl;
                return -1;
            }

            return dsdo_temp;
        }


        G4double dsdo_LQ_ep_quasi_elastic(G4LorentzVector &Pelectron_initial, G4LorentzVector &Pelectron_final,
                                          G4LorentzVector &Pproton_initial, G4LorentzVector &Pproton_final,
                                          G4double theta_e, G4double phi_e, G4double zz_cell)
        {


            // это 4-вектора в системе где протон покоиться

            G4LorentzVector Pelectron_rest_initial(0, 0, 0, 0);
            G4LorentzVector Pelectron_rest_final(0, 0, 0, 0);
            G4LorentzVector Pproton_rest_initial(0, 0, 0, 0);
            G4LorentzVector Pproton_rest_final(0, 0, 0, 0);

            /////

            G4double Eelectron_initial = Ebeam; // кин энергия электрона в ЛСО

            G4double Pe_initial_z = sqrt(Eelectron_initial * (Eelectron_initial + 2. * Me)); // импульс электрона по z
            G4ThreeVector Pe_vec_initial(0., 0., Pe_initial_z);
            // Pelectron_initial.setVect(Pe_vec_initial);
            // Pelectron_initial.setE(Eelectron_initial + Me); // полная энергия
            Pelectron_initial.setVectMag(Pe_vec_initial, Me);

            // `добавление испульса протона в дейтроне
            // перейдем в систему отсчета движищегося протона

            ///

            // это интеграл плотности вероятности по импульсу с шагом 10 МэВ?
            // см. fermi_d.F из GENBOS

            G4double distr_paris_potential[101] = {
                0.0,       6.7953602E-03, 4.6238571E-02, 0.1236217, 0.2234159, 0.3280089, 0.4262296, 0.5131546,
                0.5876556, 0.6504297,     0.7028700,     0.7465120, 0.7827951, 0.8129796, 0.8381307, 0.8591338,
                0.8767187, 0.8914840,     0.9039201,     0.9144295, 0.9233416, 0.9309281, 0.9374120, 0.9429782,
                0.9477784, 0.9519395,     0.9555655,     0.9587440, 0.9615465, 0.9640333, 0.9662549, 0.9682527,
                0.9700612, 0.9717091,     0.9732212,     0.9746166, 0.9759121, 0.9771215, 0.9782565, 0.9793263,
                0.9803385, 0.9812997,     0.9822153,     0.9830891, 0.9839256, 0.9847271, 0.9854963, 0.9862349,
                0.9869449, 0.9876275,     0.9882837,     0.9889148, 0.9895214, 0.9901043, 0.9906641, 0.9912014,
                0.9917168, 0.9922106,     0.9926836,     0.9931361, 0.9935685, 0.9939811, 0.9943746, 0.9947497,
                0.9951068, 0.9954458,     0.9957677,     0.9960729, 0.9963619, 0.9966350, 0.9968930, 0.9971365,
                0.9973655, 0.9975811,     0.9977834,     0.9979731, 0.9981509, 0.9983172, 0.9984723, 0.9986169,
                0.9987513, 0.9988762,     0.9989921,     0.9990993, 0.9991982, 0.9992898, 0.9993736, 0.9994510,
                0.9995219, 0.9995865,     0.9996458,     0.9996991, 0.9997482, 0.9997926, 0.9998326, 0.9998688,
                0.9999010, 0.9999300,     0.9999564,     0.9999794, 1.000000};

            G4double distr_bonn_potential[101] = {
                0.0,       7.0403279E-03, 4.7886062E-02, 0.1279474, 0.2310553, 0.3389300, 0.4400139, 0.5292494,
                0.6055161, 0.6695801,     0.7229201,     0.7671533, 0.8037899, 0.8341466, 0.8593346, 0.8802744,
                0.8977222, 0.9122964,     0.9245030,     0.9347540, 0.9433876, 0.9506804, 0.9568596, 0.9621121,
                0.9665918, 0.9704258,     0.9737193,     0.9765596, 0.9790186, 0.9811566, 0.9830235, 0.9846612,
                0.9861036, 0.9873803,     0.9885156,     0.9895294, 0.9904380, 0.9912565, 0.9919961, 0.9926671,
                0.9932781, 0.9938356,     0.9943457,     0.9948140, 0.9952440, 0.9956399, 0.9960047, 0.9963413,
                0.9966521, 0.9969389,     0.9972037,     0.9974484, 0.9976741, 0.9978822, 0.9980739, 0.9982503,
                0.9984124, 0.9985616,     0.9986980,     0.9988228, 0.9989371, 0.9990408, 0.9991356, 0.9992215,
                0.9992995, 0.9993696,     0.9994332,     0.9994900, 0.9995410, 0.9995866, 0.9996274, 0.9996639,
                0.9996958, 0.9997242,     0.9997495,     0.9997714, 0.9997910, 0.9998078, 0.9998228, 0.9998358,
                0.9998475, 0.9998578,     0.9998671,     0.9998750, 0.9998825, 0.9998896, 0.9998960, 0.9999021,
                0.9999081, 0.9999138,     0.9999200,     0.9999259, 0.9999323, 0.9999395, 0.9999465, 0.9999539,
                0.9999623, 0.9999706,     0.9999802,     0.9999896, 1.000000};

            G4bool is_distr_paris = true;

            G4double Qrand = 0.1;

            G4double Pnucleus = 0.;

            G4double Pfermi = 0.;

            while (1)
            {

                Qrand = G4UniformRand();

                for (G4int i = 0; i < 100; i++)
                {
                    if (is_distr_paris)
                    {
                        if (Qrand >= distr_paris_potential[i] and Qrand < distr_paris_potential[i + 1])
                        {
                            Pnucleus = 10. * (i + G4UniformRand() * 1.);

                            //  G4cout << Pnucleus << G4endl;
                        }
                    }
                    else
                    {
                        if (Qrand >= distr_bonn_potential[i] and Qrand < distr_bonn_potential[i + 1])
                        {
                            Pnucleus = 10. * (i + G4UniformRand() * 1.);
                        }
                    }
                }

                if (abs(Pnucleus) <= Md / 2.)
                    break;
                else
                    continue;
            }


            if (is_quasi_elastic_pd == false)
            {
                Pnucleus = 0.;
            }

            /*******************************
            *** This function computes reduced mass
            *** of a nucleon inside the deuteron
            *** from fermi motion so that energy
            *** is conserved
            *******************************/


            G4double red_mass_p = sqrt(Md / 2. * Md / 2. - Pnucleus * Pnucleus); // недостающая масса в ядре reduce
            G4double red_mass_n = sqrt(Md / 2. * Md / 2. - Pnucleus * Pnucleus);

            // red_mass_p = Mp;
            //  red_mass_n = Mn;

            G4double Arnd = acos(1. - 2. * G4UniformRand());
            G4double Brnd = 2. * M_PI * G4UniformRand();

            // *** Variabile per impulso di FErmi nella ntupla
            // *** Fermi motion is computed in the frame in wich deuteron is at
            // *** rest, i.e. in LAB
            Pfermi = abs(Pnucleus);

            // *** Target nucleon // мишень протон

            /*PN(1,I)    = PX OF I-TH spectator PARTICLE
            PN(2,I)    = PY OF I-TH spectator PARTICLE
            PN(3,I)    = PZ OF I-TH spectator PARTICLE
            PN(4,I)    = ENERGY OF I-TH spectator PARTICLE
            PN(5,I)    = |P| OF I-TH spectator PARTICLE*/


            G4LorentzVector PN_proton_in_deuteron;
            G4ThreeVector PN_proton_in_deuteron_3vec(0, 0, 1);
            PN_proton_in_deuteron_3vec.setMag(Pfermi);
            PN_proton_in_deuteron_3vec.setTheta(Arnd);
            PN_proton_in_deuteron_3vec.setPhi(Brnd);
            PN_proton_in_deuteron.setVectMag(PN_proton_in_deuteron_3vec, red_mass_p);


            Pproton_initial.setVectMag(PN_proton_in_deuteron_3vec, red_mass_p);
            // Pproton_initial.setVect(PN_proton_in_deuteron_3vec);
            // Pproton_initial.setE(sqrt(red_mass_p*red_mass_p+Pfermi*Pfermi));
            //
            // Pproton_initial.setE(sqrt(Mp*Mp+Pfermi*Pfermi));


            //  PN_proton_in_deuteron->setE(Pfermi*Pfermi + red_mass_p*red_mass_p);
            // PN_proton_in_deuteron->setPx( Pfermi*sin(Arnd)*cos(Brnd));
            // PN_proton_in_deuteron->setPy(Pfermi*sin(Arnd)*cos(Brnd));
            // PN_proton_in_deuteron->setPz(Pfermi*cos(Arnd));
            // G4cout << PN_proton_in_deuteron->getV().mag() << " !!! " << Pfermi << G4endl;
            // G4cout << PN_proton_in_deuteron->getX() << " !!! " << Pfermi*sin(Arnd)*cos(Brnd) << G4endl;
            // G4cout << PN_proton_in_deuteron->getY() << " !!! " <<Pfermi*sin(Arnd)*cos(Brnd) << G4endl;
            // G4cout << PN_proton_in_deuteron->getZ() << " !!! " <<Pfermi*cos(Arnd) << G4endl;
            // G4cout << PN_proton_in_deuteron->getT() << " !!! " <<sqrt(Pfermi*Pfermi+red_mass_p*red_mass_p) << G4endl;
            //   PN_N_D(5,2) = sqrt(PN_N_D(1,2)**2+PN_N_D(2,2)**2+PN_N_D(3,2)**2);
            //  PN_N_D(4,2) = sqrt(PN_N_D(5,2)**2+RM(2)**2);

            // *** Spectator nucleon // наюлюдаемый спектатор нейтрон

            // PN_N_D(1,1) = -PN_N_D(1,2);
            // PN_N_D(2,1) = -PN_N_D(2,2);
            // PN_N_D(3,1) = -PN_N_D(3,2);
            // PN_N_D(5,1) = sqrt(PN_N_D(1,1)**2+PN_N_D(2,1)**2+PN_N_D(3,1)**2);
            // PN_N_D(4,1) = sqrt(PN_N_D(5,1)**2+RM(1)**2);

            G4LorentzVector PN_neutron_in_deuteron;
            G4ThreeVector PN_neutron_in_deuteron_3vec(0, 0, 1);
            PN_neutron_in_deuteron_3vec = -PN_proton_in_deuteron_3vec;

            PN_neutron_in_deuteron.setVectM(PN_neutron_in_deuteron_3vec, red_mass_n);

            // PN_neutron_in_deuteron->boost()


            // G4cout << PN_neutron_in_deuteron->getV().mag() << " !!! " << Pfermi << G4endl;
            //  G4cout << PN_neutron_in_deuteron->getX() << " !!! " << Pfermi*sin(Arnd)*cos(Brnd) << G4endl;
            // G4cout << PN_neutron_in_deuteron->getY() << " !!! " <<Pfermi*sin(Arnd)*cos(Brnd) << G4endl;
            // G4cout << PN_neutron_in_deuteron->getZ() << " !!! " <<Pfermi*cos(Arnd) << G4endl;
            // G4cout << PN_neutron_in_deuteron->getT() << " !!! " <<sqrt(Pfermi*Pfermi+red_mass_p*red_mass_p) <<
            // G4endl;


            /*
            ***********************************************************
            *** Spectator nucleon has reduced mass. I should subtract
            *** from total energy for the reaction the mass variation of
            *** the spectator
            ***********************************************************
            */

            G4double Egamma = 0.;
            G4double Eneutron_in_deu = PN_neutron_in_deuteron.e();
            G4double ETotal = Egamma + Eneutron_in_deu;

            // *** Protone e neutrone hanno masse uguali
            G4double rm_fin = Mp;
            G4double Pneutron_in_deu = PN_neutron_in_deuteron.getV().mag();
            G4double E_sp_fin = sqrt(Pneutron_in_deu * Pneutron_in_deu + rm_fin * rm_fin);
            G4double delta_e = E_sp_fin - Eneutron_in_deu;
            ETotal = ETotal - delta_e;

            // *** CM  beta of fotone+interacting nucleon

            G4ThreeVector Beta_N_D(0, 0, 0);
            Beta_N_D.setX(PN_neutron_in_deuteron.x() / ETotal);
            Beta_N_D.setY(PN_neutron_in_deuteron.y() / ETotal);
            Beta_N_D.setZ((PN_neutron_in_deuteron.z() + Egamma) / ETotal);

            //       BETA_N_D(2)=PN_N_D(2,1)/ETOT
            //       BETA_N_D(3)=(PN_N_D(3,1)+EGAM)/ETOT
            //       BETA_N_D(4)=0.


            /**************************************************************
             ***   S is Lorentz invariant. Scalar product has only third
             ***   component if photon is along z
             *************************************************************
             ***
             *** Square of (gamma+target nucleon)
             ****/

            G4double Sinv =
                red_mass_p * red_mass_p + 2. * Eneutron_in_deu * Egamma - 2. * PN_neutron_in_deuteron.z() * Egamma;

            /****
            *** subtracting spectator mass variation
            ****/

            Sinv = Sinv + delta_e * delta_e - 2. * delta_e * (Egamma + Eneutron_in_deu);


            //////

            // Pelectron_initial, Pproton_initial `вычислены в ЛО где протон движется с импульсом Ферми, перейдем в
            // систему где протон покоиться

            //  G4cout << " !!! P is LAB e_i = " << Pelectron_initial << " p_i = " << Pproton_initial << G4endl;

            G4ThreeVector boost_vec = Pproton_initial.boostVector();

            //  boost_vec = G4ThreeVector(0,0,0);

            // G4cout << "boost = " << boost_vec << G4endl;

            G4double gamma = Pproton_initial.gamma();

            G4double gamma_z = 1. / sqrt(1. - pow(boost_vec.z(), 2.));

            //   G4cout << "beta_vector = " << boost_vec << " gamma_z = " << gamma_z << " gamma = " << gamma <<G4endl;

            Pelectron_rest_initial = Pelectron_initial;
            Pproton_rest_initial = Pproton_initial;

            Pelectron_rest_initial.boost(-boost_vec);
            Pproton_rest_initial.boost(-boost_vec);

            //  G4cout << " !!! P is REST e_i = " << Pelectron_rest_initial << " p_i = " << Pproton_rest_initial <<
            //  G4endl;

            // делаем поворот, что электрон был направлен по Z

            G4double angle_rot1 = Pelectron_rest_initial.getV().getPhi();
            Pelectron_rest_initial.rotateZ(-angle_rot1);
            G4double angle_rot2 = Pelectron_rest_initial.getV().getTheta();
            Pelectron_rest_initial.rotateY(-angle_rot2);

            Pproton_rest_initial.rotateZ(-angle_rot1);
            Pproton_rest_initial.rotateY(-angle_rot2);

            // G4cout << " !!! P is REST after rot e_i = " << Pelectron_rest_initial << " p_i = " <<
            // Pproton_rest_initial << G4endl;

            ///////

            G4double Mass = 0.;

            // Mass = Md;
            Mass = Mp;

            G4double coeff_target = 1.;

            // сокращение линейных размеров
            zz_cell /= gamma_z;
            G4double l_zz_cell_temp = l_zz_cell / gamma_z;

            if (zz_cell >= 0)
                coeff_target = (l_zz_cell_temp / 2. - zz_cell) /
                    pow(l_zz_cell_temp / 2., 2.); // коэффициент для учета труегольного распределения частиц по мишени
            if (zz_cell < 0)
                coeff_target = (l_zz_cell_temp / 2. + zz_cell) / pow(l_zz_cell_temp / 2., 2.);


            Eelectron_initial = Pelectron_rest_initial.e() - Me; // кин = полная - масса

            G4double rc = 1. + (2. * Eelectron_initial / Mass) * pow(sin(theta_e / 2.), 2.); // фактор отдачи
            // G4double Q2 = 4. * pow(Ebeam, 2.) * pow(sin(theta_e / 2.), 2.) / rc; // ???????? ????????


            G4double Eelectron_final = Eelectron_initial / rc; // энергия электрона после рассеяния

            // кинематика конечных состояний

            // в системе покоя протона

            G4double Pe_final_z = sqrt(Eelectron_final * (Eelectron_final + 2. * Me)); // импульс электрона модуль
            G4ThreeVector Pe_vec_final(0., 0., 1.);
            Pe_vec_final.setMag(Pe_final_z);
            Pe_vec_final.setTheta(theta_e);
            Pe_vec_final.setPhi(phi_e);

            Pelectron_rest_final.setVectMag(Pe_vec_final, Me);


            G4double Eproton_final = (1. - 1. / rc) * Eelectron_initial; // кин энергия протона
            G4double theta_proton_final = atan(sin(theta_e) / (rc - cos(theta_e)));

            // phi_electron = phi_e_temp + M_PI;

            G4double phi_proton_final = M_PI - phi_e;

            G4double Pp_final_z = sqrt(Eproton_final * (Eproton_final + 2. * Mp)); // импульс электрона модуль
            G4ThreeVector Pp_vec_final(0., 0., 1.);

            Pp_vec_final.setMag(Pp_final_z);
            Pp_vec_final.setTheta(theta_proton_final);
            Pp_vec_final.setPhi(phi_proton_final);

            Pproton_rest_final.setVectMag(Pp_vec_final, Mp);

            // обратный поворот

            Pelectron_rest_final.rotateY(angle_rot2);
            Pelectron_rest_final.rotateZ(angle_rot1);

            Pproton_rest_final.rotateY(angle_rot2);
            Pproton_rest_final.rotateZ(angle_rot1);


            Pelectron_rest_initial.rotateY(angle_rot2);
            Pelectron_rest_initial.rotateZ(angle_rot1);

            Pproton_rest_initial.rotateY(angle_rot2);
            Pproton_rest_initial.rotateZ(angle_rot1);

            // в ЛСО

            // обратный буст

            Pelectron_final = Pelectron_rest_final;
            Pproton_final = Pproton_rest_final;

            Pelectron_final.boost(boost_vec);
            Pproton_final.boost(boost_vec);

            ////


            //  G4cout << " !!! P is LAB e_f = " << Pelectron_final << " p_f = " << Pproton_final << G4endl;
            //  G4cout << " !!! P is REST e_f = " << Pelectron_rest_final << " p_f = " << Pproton_rest_final << G4endl;


            G4double Q2 = 4. * Eelectron_initial * Eelectron_final * pow(sin(theta_e / 2.), 2.);

            Q2 = -(Pelectron_rest_final - Pelectron_rest_initial).mag2(); // это точно из 4-импульса

            G4double Qf2 = Q2 / pow(Fm, 2.);
            G4double nu = Q2 / (4. * Mass * Mass);
            G4double epsilon = (1. + nu) * pow(tan(theta_e / 2.), 2.);

            // G4double Ge = Ge_proton_spline(Q2); // Формфакторы
            // G4double Gm = Gm_proton_spline(Q2);

            G4double Ge = Ge_proton_spline_kelly(Q2); // Формфакторы
            G4double Gm = Gm_proton_spline_kelly(Q2);

         //   G4cout << "sqrt(Q2) = " << sqrt(Q2) << " Q2 = " << Q2/1e6 << " GeV^2 " << " Ge = " << Ge_proton_spline(Q2) << " " <<  Ge_proton_spline_kelly(Q2) << " Gm = "
         //   << Gm_proton_spline(Q2) << " " << Gm_proton_spline_kelly(Q2) << G4endl;

           //   G4cout << "Q2 = " << Q2/1e6 << " GeV^2 " << " Ge/Gd = " <<  Ge_proton_spline_kelly(Q2)/Ge_proton_spline(Q2)  << " Gm/(mu*Gd) = "
           //   << Gm_proton_spline_kelly(Q2) / Gm_proton_spline(Q2) << G4endl;


            G4double A = (nu * pow(Gm, 2.0) + pow(Ge, 2.0)) / (1.0 + nu); // Структурные функции

            G4double B = 2.0 * nu * pow(Gm, 2.);

            G4double S = A + B * pow(tan(theta_e / 2.), 2.);

            G4double gamma_electron = Eelectron_initial / Me;
            G4double beta_electron = sqrt(1. - 1. / pow(gamma_electron, 2.0));

            // pow(cos(theta_e / 2.), 2.)

            // G4double Mt = (pow(alpha_em * Fm, 2.) * 10000. / 4.) * (1. - pow(beta_electron,2.) * pow(sin(theta_e
            // / 2.), 2.))
            // * S / (rc * pow(sin(theta_e / 2.), 4.) * pow(Eelectron_initial, 2.)); Mott cross-Section , mkb

            //  G4cout << "Mt = " << Mt << G4endl;

            G4double Mt = (pow(alpha_em * Fm, 2.) * 10000. * 4.) *
                (1. - pow(beta_electron, 2.) * pow(sin(theta_e / 2.), 2.)) * S * pow(Eelectron_initial, 2.) /
                (pow(rc, 3.) * pow(Q2, 2.));


            //  G4cout << "Mt = " << Mt << G4endl;

            //  G4cout << "Q2 [Gev^2] = " << Q2 / 1e6 << " Ge = " << Ge  << " Gm = " << Gm << " S = " << S << G4endl;

            //   G4cout << "proton Q2 [Gev^2] = " << Q2 / 1e6  << " S = " << S << G4endl;

            G4double q3md = sqrt(pow(Ebeam, 2.) + pow(Eelectron_final, 2.) -
                                 2. * Eelectron_initial * Eelectron_final *
                                     cos(theta_e)); // вектор передачи от начального электрона к конечному
            G4double cos_theta_otdachi = (pow(q3md, 2.) + pow(Eelectron_initial, 2.) - pow(Eelectron_final, 2.)) /
                (2. * q3md * Eelectron_initial); // косинус угла между векторами по теореме косинусов
            G4double theta_otdachi =
                acos(cos_theta_otdachi); // угол дейтрона отдачи для формулы поляризованного сечения

            ////


            //  G4cout << "sqrt(Q2) = " << sqrt(Q2) << " " << -(Pelectron_rest_final-Pelectron_rest_initial).mag() << "
            //  " << -(Pelectron_final-Pelectron_initial).mag() << G4endl;

            G4double dsdo_temp = sin(theta_e) * Mt * coeff_target;

            G4double factor_cross_sec_to_lo = 1.; // для перевода дифф сечения из покоящийся СО протона в ЛАБ

            // factor_cross_sec_to_lo = 4.*sin(theta_e/2.);

            G4double gamma_m = Me / Mp;

            // factor_cross_sec_to_lo = pow(1. + gamma_m*gamma_m + 2.*gamma_m*cos(theta_e), 3./2.) / abs(1. + gamma_m *
            // cos(theta_e));

            // factor_cross_sec_to_lo = 4. * sin(theta_e / 2.);

            // G4cout << "factor = " << factor_cross_sec_to_lo << G4endl;

            // dsdo_temp *= factor_cross_sec_to_lo;

            // G4cout << "dsdo_temp = " << dsdo_temp << G4endl;

            if (dsdo_temp > 10.0)
            {
                //G4cout << "dsdo_temp = " << dsdo_temp << G4endl;
                return -1;
            }

            return dsdo_temp;
        }


        /////////


        G4double assym_calc(G4double a, G4double b, G4double r = -2.)
        {
            // a - число событий с pol=1 (P_zz=+1), b - число событий с pol=2 (P_zz=-2)
            G4double assym_temp = (a - b) / (b - r * a);
            //  cout<< pf1t[1]<<endl;
            return assym_temp;
        }

        G4double assym_err_calc(G4double a, G4double b, G4double r = -2.)
        {

            //  float err_temp = 3. * sqrt(a * b * (a + b)) / pow(2. * a + b,2.); // Получено по квадратичному сложению
            //  ошибок
            G4double err_temp = fabs(r - 1.) * sqrt(a * b * (a + b)) / pow(b - r * a, 2.);
            return err_temp;
        }


        void random_Neumann_LQ_method1(G4double initial_theta_e, G4double final_theta_e, G4double initial_zz_cell,
                                       G4double final_zz_cell, G4double Pzz1, G4double Pzz2, G4double max_f,
                                       G4double &theta_electron, G4double &phi_electron, G4double &energy_electron,
                                       G4double &theta_deuteron, G4double &phi_deuteron, G4double &energy_deuteron,
                                       G4double &Pzz, G4double &xx_cell, G4double &yy_cell, G4double &zz_cell,
                                       G4double &dsdo);

        void random_Neumann_LQ_method2(G4double initial_x_counter_e, G4double final_x_counter_e,
                                       G4double initial_y_counter_e, G4double final_y_counter_e,
                                       G4double initial_zz_cell, G4double final_zz_cell, G4double Pzz1, G4double Pzz2,
                                       G4double max_f, G4double &theta_electron, G4double &phi_electron,
                                       G4double &energy_electron, G4double &theta_deuteron, G4double &phi_deuteron,
                                       G4double &energy_deuteron, G4double &Pzz, G4double &xx_cell, G4double &yy_cell,
                                       G4double &zz_cell, G4double &dsdo);


        void random_Neumann_LQ_ep_method2(G4double initial_x_counter_e, G4double final_x_counter_e,
                                          G4double initial_y_counter_e, G4double final_y_counter_e,
                                          G4double initial_zz_cell, G4double final_zz_cell, G4double max_f,
                                          G4double& theta_electron, G4double& phi_electron, G4double& energy_electron,
                                          G4double& theta_proton, G4double& phi_proton, G4double& energy_proton,
                                          G4double& xx_cell, G4double& yy_cell, G4double& zz_cell, G4double &dsdo);

        void random_Neumann_LQ_ep_plus_ed_method2(G4double initial_x_counter_e, G4double final_x_counter_e,
                                  G4double initial_y_counter_e, G4double final_y_counter_e,
                                  G4double initial_zz_cell, G4double final_zz_cell, G4double *max_f,
                                  G4double& theta_electron, G4double& phi_electron, G4double& energy_electron,
                                  G4double& theta_proton, G4double& phi_proton, G4double& energy_proton,
                                  G4double& theta_deuteron, G4double& phi_deuteron, G4double& energy_deuteron,
                                  G4double& Pzz,
                                  G4double& xx_cell, G4double& yy_cell, G4double& zz_cell, G4double &dsdo, G4int &type_cross_section);

        void random_Neumann_LQ_ep_quasi_elastic_method2(G4double initial_x_counter_e, G4double final_x_counter_e,
                                                        G4double initial_y_counter_e, G4double final_y_counter_e,
                                                        G4double initial_zz_cell, G4double final_zz_cell, G4double max_f,
                                                        G4double& theta_electron, G4double& phi_electron,
                                                        G4double& energy_electron,
                                                        G4double& theta_proton, G4double& phi_proton,
                                                        G4double& energy_proton,
                                                        G4double& xx_cell, G4double& yy_cell, G4double& zz_cell, G4double &dsdo);

    };






}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
