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

class G4AttDef;
class G4AttValue;

#define DENSITY		1.032	// scintill. density - for light output calculation

namespace Cosmic
{

/// Calorimeter hit class
///
/// It defines data members to store the the energy deposit and track lengths
/// of charged particles in a selected volume:
/// - fEdep, fTrackLength

class HadronCalorimeterHit : public G4VHit
{
  public:
    HadronCalorimeterHit();
    HadronCalorimeterHit(G4int layerID);
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
   // void AddEdep(G4double de);
   // void AddTrackLength(G4double dl);
  //  void AddLO(G4double de, G4ThreeVector pos, G4ThreeVector delta); // световыход

    // get methods
 //   G4double GetEdep() const;
 //   G4double GetTrackLength() const;
    G4double GetLO() const;
    G4double GetA1() const;
    G4double GetA2() const;
  //  G4ThreeVector GetPos() const;
    G4ThreeVector GetVPos() const;
    G4RotationMatrix GetVRot() const;
    G4int GetBlkN() const;
   // G4double GetToF() const;
    G4bool GetTrig() const;
    G4int GetNprim() const;

    // set methods
    inline void SetBlkN(G4int n)	{ blkN = n;};
   // inline void SetEdep(G4double v)	{ fEdep=v; };
    inline void SetLO(G4double v)	{ LO=v; };
    inline void SetA1(G4double v)	{ A1=v; };
    inline void SetA2(G4double v)	{ A2=v; };
    inline void SetHalfLength(G4double v){ halflength=v; };
    inline void SetAbsorbtion(G4double v){ absorbtion=v; };

    inline void SetVPos(G4ThreeVector v)	{ Vpos =v;} ;
    inline void SetVRot(G4RotationMatrix v)	{ Vrot =v;} ;
   // inline void SetToF(G4double t)	{ToF = t;};
    inline void SetTrig(G4bool v)	{Trig=v;};
    inline void SetNprim(G4int v)	{Nprim=v;};



    void SetEdep(G4double de) { fEdep = de; }
    void AddEdep(G4double de) { fEdep += de; }
    G4double GetEdep() const { return fEdep; }

    void SetTrackLength(G4double dl) { fTrackLength = dl; }
    void AddTrackLength(G4double dl) { fTrackLength += dl; }
    G4double GetTrackLength() const { return fTrackLength; }


    void SetToF(G4double tof) { fToF = tof; }
    void AddToF(G4double tof) { fToF = tof; }
    G4double GetToF() const { return fToF; }


    void SetPos(G4ThreeVector xyz) { fLocalPos = xyz; }
    G4ThreeVector GetPos() const { return fLocalPos; }

    void SetRot(G4RotationMatrix rmat) { fRot = rmat; }
    G4RotationMatrix GetRot() const { return fRot; }

    void SetLayerID(G4int z) { fLayerID = z; }
    G4int GetLayerID() const { return fLayerID; }

    void SetLogV(G4LogicalVolume* val) { fPLogV = val; }
    const G4LogicalVolume* GetLogV() const { return fPLogV; }

  private:
    G4int fLayerID = -1;
    G4double fToF = 0.;
    G4double fEdep = 0.;        ///< Energy deposit in the sensitive volume
    G4double fTrackLength = 0.; ///< Track length in the  sensitive volume
    G4ThreeVector fLocalPos = {0.,0.0,0.0};
    G4ThreeVector fWorldPos = {0.,0.0,0.0};
    G4RotationMatrix fRot;
    const G4LogicalVolume* fPLogV = nullptr;

    void CalcRho(G4ThreeVector v);

    G4double halflength = 0.;
    G4double absorbtion = 0.;

    G4double threshold = 0.;
    G4double pde = 0.;


    G4ThreeVector Vpos = {0.,0.0,0.0};
    G4RotationMatrix Vrot;
    G4double LO = 0.;
    G4int blkN = 0;
    G4double A1 = 0.;
    G4double A2 = 0.;
    G4bool Trig;
    G4int Nprim = 0;
    G4double rhoX = 0.,rhoZ = 0.;

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

   // inline void HadronCalorimeterHit::AddEdep(G4double de) {
   //     fEdep += de;

  //  }
  //  inline void HadronCalorimeterHit::AddTrackLength(G4double dl) {
  //      fTrackLength += dl;
  //  }

   // inline G4double HadronCalorimeterHit::GetEdep() const { return fEdep; }

  //  inline G4double HadronCalorimeterHit::GetTrackLength() const { return fTrackLength; }

    inline G4double HadronCalorimeterHit::GetLO() const { return LO; }

    inline G4double HadronCalorimeterHit::GetA1() const { return A1; }

    inline G4double HadronCalorimeterHit::GetA2() const { return A2; }

    //inline G4ThreeVector HadronCalorimeterHit::GetPos() const { return fLocalPos; }

    inline G4ThreeVector HadronCalorimeterHit::GetVPos() const { return Vpos; }

    inline G4RotationMatrix HadronCalorimeterHit::GetVRot() const { return Vrot; }

    inline G4int HadronCalorimeterHit::GetBlkN() const { return blkN; }

  //  inline G4double HadronCalorimeterHit::GetToF() const { return ToF; }

    inline G4bool HadronCalorimeterHit::GetTrig() const { return Trig; };

    inline G4int HadronCalorimeterHit::GetNprim() const { return Nprim; };

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
