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
/// \file PlasticHit.hh
/// \brief Definition of the Cosmic::PlasticHit class

#ifndef CosmicPlasticHit_h
#define CosmicPlasticHit_h 1

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


namespace Cosmic
{

/// Calorimeter hit class
///
/// It defines data members to store the the energy deposit and track lengths
/// of charged particles in a selected volume:
/// - fEdep, fTrackLength

class PlasticHit : public G4VHit
{
  public:
    PlasticHit();
    PlasticHit(G4int layerID);
    //PlasticHit(G4double, G4double, G4double, G4ThreeVector, const G4RotationMatrix*);
    //PlasticHit(G4ThreeVector, const G4RotationMatrix*);
    PlasticHit(const PlasticHit& right) = default;
    ~PlasticHit() override;

    // operators
   // const PlasticHit& operator=(const PlasticHit& right);
    PlasticHit& operator=(const PlasticHit &right) = default;
    G4bool operator==(const PlasticHit& right) const;

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

    G4RotationMatrix GetVRot() const;
    G4int GetBlkN() const;

    // G4double GetToF() const;
    G4bool GetTrig() const;
    G4double GetThreshold() const;


    G4int GetNprim() const;

    // set, get and add methods
    inline void SetBlkN(G4int n) { blkN = n; };

    inline void SetTrig(G4bool val) { flaq_is_trig_tof = val; };
    inline void SetThreshold(G4double val) { threshold = val; };

    inline void SetNprim(G4int v) { Nprim = v; };

    inline void SetTrackID(G4int TrackID) { fTrackID = TrackID; }

    void AddTrackID(G4int TrackID) { fTrackID = TrackID; }

    inline G4int GetTrackID() const { return fTrackID; }

    inline void SetEdep(G4double de) { fEdep = de; }

    void AddEdep(G4double de)
    {
        fEdep += de;
        fdEdep = de;
    }

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

    inline void SetPosTheta(G4double theta) { fPosTheta = theta; }
    inline G4double GetPosTheta() const { return fPosTheta; }

    inline void SetPosThetaGlob(G4double thetaGlob) { fPosThetaGlob = thetaGlob; }
    inline G4double GetPosThetaGlob() const { return fPosThetaGlob; }

    inline void SetPosPhi(G4double phi) { fPosPhi = phi; }
    inline G4double GetPosPhi() const { return fPosPhi; }

    inline void SetPosPhiGlob(G4double phiGlob) { fPosPhiGlob = phiGlob; }
    inline G4double GetPosPhiGlob() const { return fPosPhiGlob; }

    inline void SetToF(G4double tof) { fToF = tof; }
    void AddToF(G4double tof)
    {
        // fToF = tof;


        //  G4cout << threshold / MeV << G4endl;

        // G4float ToF_temp = 0.;


        // if (fLO > plastic_threshold) flaq_is_trig_tof = true;
        if (fLO > threshold && flaq_is_trig_tof == true)
        {
            fToF = tof;
            flaq_is_trig_tof = false;
        }

        // G4cout << flaq_is_trig_tof << G4endl;
    }
    inline G4double GetToF() const { return fToF; }


    inline void SetWorldPos(G4ThreeVector pos) { fWorldPos = pos; }
    void AddWorldPos(G4ThreeVector pos);
    inline G4ThreeVector GetWorldPos() const
    {
        //return fWorldPos;
        return fLO>0. ? fWorldPos / fLO: fWorldPos;
    }

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

    void SetLayerID(G4int z) { fLayerID = z; }
    G4int GetLayerID() const { return fLayerID; }

    void SetLogV(G4LogicalVolume* val) { fPLogV = val; }
    const G4LogicalVolume* GetLogV() const { return fPLogV; }


private:
    G4int fLayerID = -1;
    G4int fTrackID = -1;
    G4double fToF = 0.; // Time of Flight
    G4double fEdep = 0.;        ///< Energy deposit in the sensitive volume
    G4double fdEdep = 0.;
    G4double fLO = 0.;        ///< Light Output in the sensitive volume
    G4double fdLO = 0.;
    G4double fTrackLength = 0.; ///< Track length in the  sensitive volume
    G4double fA1 = 0.;  // Амплитуды с двух концов счетчиков
    G4double fA2 = 0.;
    G4double fT1 = 0.;  // Времена с двух концов счетчиков
    G4double fT2 = 0.;
    G4ThreeVector fHalfLength = G4ThreeVector();
    G4double fAbsorbtion = 0.;
    G4double fPosTheta = 0., fPosPhi = 0.; // углы по направлению движения/импульса
    G4double fPosThetaGlob = 0., fPosPhiGlob = 0.; // углы в глобальной системе координат
    G4ThreeVector fLocalPos = G4ThreeVector(); // точка взаимодействия в координатах детектора
    G4ThreeVector fWorldPos = G4ThreeVector();
    G4ThreeVector fVPos = G4ThreeVector(); // вершина (точка генерации)
    G4RotationMatrix fRot = G4RotationMatrix();
    G4double fRhoX = 0., fRhoZ = 0., fRhoY = 0.;
    const G4LogicalVolume *fPLogV = nullptr;

    void CalcRho(G4ThreeVector v);

    G4double threshold = 0.1 * MeV;
    G4double pde = 0.;



    G4int blkN = 0;

    G4bool Trig;
    G4int Nprim = 0;

    G4bool flaq_is_trig_tof = true;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using PlasticHitsCollection = G4THitsCollection<PlasticHit>;

extern G4ThreadLocal G4Allocator<PlasticHit>* PlasticHitAllocator;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    inline void *PlasticHit::operator new(size_t) {
        if (!PlasticHitAllocator) {
            PlasticHitAllocator = new G4Allocator<PlasticHit>;
        }
        void *hit;
        hit = (void *) PlasticHitAllocator->MallocSingle();
        return hit;
    }

    inline void PlasticHit::operator delete(void *hit) {
        if (!PlasticHitAllocator) {
            PlasticHitAllocator = new G4Allocator<PlasticHit>;
        }
        PlasticHitAllocator->FreeSingle((PlasticHit *)hit);
    }

    // inline void PlasticHit::AddEdep(G4double de) {
   //     fEdep += de;

    //  }
    //  inline void PlasticHit::AddTrackLength(G4double dl) {
    //      fTrackLength += dl;
    //  }


    inline G4int PlasticHit::GetBlkN() const { return blkN; }


    inline G4bool PlasticHit::GetTrig() const { return flaq_is_trig_tof; };

    inline G4double PlasticHit::GetThreshold() const { return threshold; };

    inline G4int PlasticHit::GetNprim() const { return Nprim; };

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
