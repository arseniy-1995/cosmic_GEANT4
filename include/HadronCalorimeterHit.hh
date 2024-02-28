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
/// \file HadronCalorimeterHit.hh
/// \brief Definition of the Cosmic::HadronCalorimeterHit class

#ifndef CosmicHadronCalorimeterHit_h
#define CosmicHadronCalorimeterHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"
#include "G4LogicalVolume.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4RotationMatrix.hh"

//#include "PlasticHit.hh"

class G4AttDef;
class G4AttValue;

//class PlasticHit;

namespace Cosmic
{

/// Calorimeter hit class
///
/// It defines data members to store the the energy deposit and track lengths
/// of charged particles in a selected volume:
/// - fEdep, fTrackLength

class HadronCalorimeterHit : public G4VHit /*public PlasticHit*/
{
  public:
    HadronCalorimeterHit();
    HadronCalorimeterHit(G4int layerxID, G4int layeryID, G4int layerzID);
    //HadronCalorimeterHit(G4double, G4double, G4double, G4ThreeVector, const G4RotationMatrix*);
    //HadronCalorimeterHit(G4ThreeVector, const G4RotationMatrix*);
    HadronCalorimeterHit(const HadronCalorimeterHit& right) = default;
    ~HadronCalorimeterHit() override;

    // operators
   // const HadronCalorimeterHit& operator=(const HadronCalorimeterHit& right);
    HadronCalorimeterHit& operator=(const HadronCalorimeterHit &right) = default;
    G4bool operator==(const HadronCalorimeterHit& right) const;

    inline void* operator new(size_t);
    inline void  operator delete(void *aHit);

    // methods from base class
    void Draw()  override{}
    const std::map<G4String,G4AttDef>* GetAttDefs() const override;
    std::vector<G4AttValue>* CreateAttValues() const override;
    void Print() override;

    // methods to handle data

    // set, get and add methods


// void AddEdep(G4double de);
        // void AddTrackLength(G4double dl);
        //  void AddLO(G4double de, G4ThreeVector pos, G4ThreeVector delta); // световыход

        G4RotationMatrix GetVRot() const;
    G4int GetBlkN();
    // G4double GetToF() const;
        G4bool GetTrig() const;
        G4int GetNprim() const;

        // set, get and add methods
        inline void SetBlkN(G4int n)	{ blkN = n;};

        inline void SetTrig(G4bool v)	{Trig=v;};
        inline void SetNprim(G4int v)	{Nprim=v;};


        inline void SetEdep(G4double de) { fEdep = de; }
        void AddEdep(G4double de) { fEdep += de; }
        inline G4double GetEdep() const { return fEdep; }

        inline void SetLO(G4double lo) { fLO = lo; }
        void AddLO(G4double de, G4ThreeVector pos, G4ThreeVector delta, G4double velocity, G4double ToF);
        // void AddLO(G4double lo) { fLO += lo; }
        inline G4double GetLO() const { return fLO; }

        inline void SetA1(G4double a1) { fA1 = a1; }
        //  void AddA1(G4double a1) { fA1 += a1; }
        inline G4double GetA1() const { return fA1; }

        inline void SetA2(G4double a2) { fA2 = a2; }
        //  void AddA2(G4double a2) { fA2 += a2; }
        inline G4double GetA2() const { return fA2; }


        inline void SetT1(G4double t1) { fT1 = t1; }
        //  void AddT1(G4double t1) { fT1 += t1; }
        inline G4double GetT1() const { return fT1; }

        inline void SetT2(G4double t2) { fT2 = t2; }
        //  void AddT2(G4double t2) { fT2 += t2; }
        inline G4double GetT2() const { return fT2; }


        inline void SetTrackLength(G4double dl) { fTrackLength = dl; }
        void AddTrackLength(G4double dl) { fTrackLength += dl; }
        inline G4double GetTrackLength() const { return fTrackLength; }


        inline void SetToF(G4double tof) { fToF = tof; }
        void AddToF(G4double tof) { fToF = tof; }
        inline G4double GetToF() const { return fToF; }


        inline void SetWorldPos(G4ThreeVector pos) { fWorldPos = pos; }
        void AddWorldPos(G4ThreeVector pos);
        inline G4ThreeVector GetWorldPos() const { return fWorldPos; }

        inline void SetLocalPos(G4ThreeVector pos) { fLocalPos = pos; }
        void AddLocalPos(G4ThreeVector pos);
        inline G4ThreeVector GetLocalPos() const { return fLocalPos; }

        inline void SetVPos(G4ThreeVector pos)	{ fVPos =pos;}
        //void AddVPos(G4ThreeVector pos) { fVPos = pos;}
        inline G4ThreeVector GetVPos() const { return fVPos; }

        inline void SetRot(G4RotationMatrix rmat) { fRot = rmat; }
        // void AddRot(G4ThreeVector rmat) { fRot = rmat;}
        inline G4RotationMatrix GetRot() const { return fRot; }

        inline void SetHalfLength(G4ThreeVector xyz)	{ fHalfLength =xyz;}
        //void AddHalfLength(G4ThreeVector x) { fHalfLength_X =x}
        inline G4ThreeVector GetHalfLength() const { return fHalfLength; }

        inline void SetAbsorbtion(G4double absorb){ fAbsorbtion  = absorb; }
        //void AddAbsorbtion(G4ThreeVector absorb) { fAbsorbtion = absorb}
        inline G4double GetAbsorbtion() const { return fAbsorbtion; }

        inline G4double GetRhoX() const	{ return fRhoX; };
        inline G4double GetRhoY() const	{ return fRhoY; };
        inline G4double GetRhoZ() const	{ return fRhoZ; };
        inline void SetRho(G4ThreeVector rho)	{
            fRhoX = rho.x();
            fRhoY = rho.y();
            fRhoZ = rho.z();
        }

    void SetLayerID(G4int x, G4int y, G4int z) { fLayerxID = x;fLayeryID = y;fLayerzID = z; }
    G4int GetLayerxID() const { return fLayerxID; }
    G4int GetLayeryID() const { return fLayeryID; }
    G4int GetLayerzID() const { return fLayerzID; }

    void SetLogV(G4LogicalVolume* val) { fPLogV = val; }
    const G4LogicalVolume* GetLogV() const { return fPLogV; }

  private:
    G4int fLayerxID = -1;
    G4int fLayeryID = -1;
    G4int fLayerzID = -1;
    G4double fToF = 0.; // Time of Flight
    G4double fEdep = 0.;        ///< Energy deposit in the sensitive volume
    G4double fLO = 0.;        ///< Light Output in the sensitive volume
    G4double fTrackLength = 0.; ///< Track length in the  sensitive volume
    G4double fA1 = 0.;  // Амплитуды с двух концов счетчиков
    G4double fA2 = 0.;
    G4double fT1 = 0.;  // Времена с двух концов счетчиков
    G4double fT2 = 0.;
    G4ThreeVector fHalfLength = G4ThreeVector();
    G4double fAbsorbtion  = 0.;
    G4ThreeVector fLocalPos = G4ThreeVector(); // точка взаимодействия в координатах детектора
    G4ThreeVector fWorldPos = G4ThreeVector();
    G4ThreeVector fVPos = G4ThreeVector(); // вершина (точка генерации)
    G4RotationMatrix fRot = G4RotationMatrix();
    G4double fRhoX = 0.,fRhoZ = 0., fRhoY = 0.;
    const G4LogicalVolume* fPLogV = nullptr;

    void CalcRho(G4ThreeVector v);

    G4double threshold = 0.;
    G4double pde = 0.;

    G4int blkN = 0;

    G4bool Trig;
    G4int Nprim = 0;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using HadronCalorimeterHitsCollection = G4THitsCollection<HadronCalorimeterHit>;

extern G4ThreadLocal G4Allocator<HadronCalorimeterHit>* HadronCalorimeterHitAllocator;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    inline void *HadronCalorimeterHit::operator new(size_t) {
        if (!HadronCalorimeterHitAllocator) {
            HadronCalorimeterHitAllocator = new G4Allocator<HadronCalorimeterHit>;
        }
        void *hit;
        hit = (void *) HadronCalorimeterHitAllocator->MallocSingle();
        return hit;
    }

    inline void HadronCalorimeterHit::operator delete(void *hit) {
        if (!HadronCalorimeterHitAllocator) {
            HadronCalorimeterHitAllocator = new G4Allocator<HadronCalorimeterHit>;
        }
        HadronCalorimeterHitAllocator->FreeSingle((HadronCalorimeterHit *) hit);
    }


    inline G4int HadronCalorimeterHit::GetBlkN() { return blkN; }


    inline G4bool HadronCalorimeterHit::GetTrig() const { return Trig; };

    inline G4int HadronCalorimeterHit::GetNprim() const { return Nprim; };


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
