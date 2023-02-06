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
/// \file PlasticHit.cc
/// \brief Implementation of the Cosmic::PlasticHit class

#include "PlasticHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include "Constants.hh"

#include <iomanip>

namespace Cosmic {

    G4ThreadLocal G4Allocator<PlasticHit> *PlasticHitAllocator = nullptr;
    //   G4ThreadLocal G4Allocator<PlasticHit>* PlasticHitAllocator;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    PlasticHit::PlasticHit() {

        // G4cout<<"??????"<<G4endl;

        //  fEdep = A1 = A2 = LO = 0.;
        //  ToF = -1.*ns;
        //  fLocalPos = G4ThreeVector();
        //  pde=0.;
        // Trig = false;
        //  blkN = -1;
        // Vpos = G4ThreeVector();
        //  Vrot=G4RotationMatrix();
        //  Nprim=-1;
        //  rhoX=rhoZ=1000.;
    }

    PlasticHit::PlasticHit(G4int layerID)
            : fLayerID(layerID) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
    PlasticHit::PlasticHit(G4double l, G4double a, G4double thr, G4ThreeVector vp, const G4RotationMatrix *rr):
            halflength(l),absorbtion(a),threshold(thr)
    {
        fEdep = A1 = A2 = LO = 0.;
        ToF = -1.*ns;
        fPos = G4ThreeVector(); pde=0.;
        Trig = false;
        blkN = -1;
        Nprim=-1;
        Vpos = vp;
        if(rr) Vrot = *rr; else Vrot=G4RotationMatrix();
        rhoX=rhoZ=1000.;
        CalcRho(fPos);

    }
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
    PlasticHit::PlasticHit(G4ThreeVector vp,const G4RotationMatrix *rr):
            halflength(1.*mm),absorbtion(1000.*m),threshold(0.)
    {
        fEdep = A1 = A2 = LO = 0.;
        ToF = -1.*ns;
        Pos = G4ThreeVector();
        pde=0.;
        Trig = false;
        blkN = -1;
        Nprim= -1;
        Vpos = vp;
        if(rr) Vrot = *rr; else Vrot=G4RotationMatrix();
        rhoX=rhoZ=1000.;
        CalcRho(Pos);


    }

*/

    PlasticHit::~PlasticHit() {}

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// CalorHit::CalorHit(const CalorHit& right)
//   : G4VHit()
// {
//   fEdep        = right.fEdep;
//   fTrackLength = right.fTrackLength;
// }

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//const PlasticHit &PlasticHit::operator=(const PlasticHit &right) {
//        fEdep = right.fEdep;
//        fTrackLength = right.fTrackLength;


    //       halflength = right.halflength;
    //       absorbtion = right.absorbtion;
    //       threshold = right.threshold;
    //      fLocalPos = right.fLocalPos;
    //      pde = right.pde;
    //      Vpos = right.Vpos;
    //      Vrot = right.Vrot;
    //      LO = right.LO;
    //      blkN = right.blkN;
    //      Nprim = right.Nprim;
    //       A1 = right.A1;
    //       A2 = right.A2;
    //       ToF = right.ToF;
    //      Trig = right.Trig;

    //       return *this;
    // }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



    G4bool PlasticHit::operator==(const PlasticHit &right) const {

        return (this == &right) ? true : false;

    }



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    const std::map<G4String, G4AttDef> *PlasticHit::GetAttDefs() const {
        G4bool isNew;
        auto store = G4AttDefStore::GetInstance("EmCalorimeterHit", isNew);

        if (isNew) {
            (*store)["HitType"]
                    = G4AttDef("HitType", "Hit Type", "Physics", "", "G4String");

            (*store)["ID"]
                    = G4AttDef("ID", "ID", "Physics", "", "G4int");

            (*store)["Energy"]
                    = G4AttDef("Energy", "Energy Deposited", "Physics", "G4BestUnit",
                               "G4double");

            (*store)["Pos"]
                    = G4AttDef("Pos", "Position", "Physics", "G4BestUnit",
                               "G4ThreeVector");

            (*store)["LVol"]
                    = G4AttDef("LVol", "Logical Volume", "Physics", "", "G4String");
        }
        return store;
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    std::vector<G4AttValue> *PlasticHit::CreateAttValues() const {
        auto values = new std::vector<G4AttValue>;

        values
                ->push_back(G4AttValue("HitType", "EmCalorimeterHit", ""));
        values
                ->push_back(G4AttValue("ID", G4UIcommand::ConvertToString(fLayerID), ""));
        values
                ->push_back(G4AttValue("Energy", G4BestUnit(fEdep, "Energy"), ""));
        values
                ->push_back(G4AttValue("Pos", G4BestUnit(fLocalPos, "Length"), ""));

        if (fPLogV)
            values->push_back(G4AttValue("LVol", fPLogV->GetName(), ""));
        else
            values->push_back(G4AttValue("LVol", " ", ""));

        return values;
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void PlasticHit::Print() {
        G4cout
                << "Edep: "
                << std::setw(7) << G4BestUnit(fEdep, "Energy")
                << " track length: "
                << std::setw(7) << G4BestUnit(fTrackLength, "Length")
                << G4endl;
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void PlasticHit::CalcRho(G4ThreeVector pos) {
        G4ThreeVector r = pos - fVPos;
        G4double rx = sqrt(r.x() * r.x() + r.y() * r.y());
        G4double rz = sqrt(r.z() * r.z() + r.y() * r.y());
        G4double ry = sqrt(r.x() * r.x() + r.z() * r.z());

        if (rx < fRhoX)fRhoX = rx;
        if (rz < fRhoZ)fRhoZ = rz;
        if (ry < fRhoY)fRhoY = ry;

    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void PlasticHit::AddWorldPos(G4ThreeVector pos) {
        fWorldPos = pos;
        CalcRho(pos);
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void PlasticHit::AddLocalPos(G4ThreeVector pos) {
        fLocalPos = pos;
        CalcRho(pos);
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PlasticHit::AddLO(G4double de, G4ThreeVector pos, G4ThreeVector delta, G4double velocity, G4double ToF) {
        G4double dlo = 0.;
        G4double dA1 = 0., dA2 = 0.;
        G4double dT1 = 0., dT2 = 0.;
        G4double dx = delta.mag(); // модуль вектора
      //  G4ThreeVector posit = pos - fVPos;
       // posit.transform(fRot);
       // G4double xpos = posit.getX();
        G4double xpos = fLocalPos.getX();

//  G4cout << "  => " << posit << " xx = " << xx << G4endl;

        if (dx > 0.0) {
            G4double a = (de / MeV) / (dx / cm * DENSITY_LO);
            //G4double k1 = 0.001, k2 = 0.000009;
            G4double k1 = 0.011, k2 = 0.000009;
            G4double f = 1.0 + k1 * a + k2 * pow(a, 2.);
            dlo = de / f;
//G4cout << " ___ a="<<a<<" f="<<f<<" lo="<<lo<<G4endl;
            fLO += dlo;

            G4double halflength = fHalfLength.x();
            //  G4cerr << "!!! "<<halflength/cm<<G4endl;

            G4double absorbtion = 150.0 * cm;
            G4double n_refract = 1.58;
            G4double velocity_light = c_light / n_refract;

            dA1 = dlo * exp(-(halflength - xpos) / absorbtion);
            dA2 = dlo * exp(-(halflength + xpos) / absorbtion);
            dT1 = (halflength - xpos) / velocity_light + ToF; // распространение света
            dT2 = (halflength + xpos) / velocity_light + ToF;
            fA1 += dA1;
            fA2 += dA2;
            fT1 = dT1;
            fT2 = dT2;

        }

    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
}
