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
/// \brief Definition of the Cosmic_sim::PlasticHit class

#ifndef Cosmic_simPlasticHit_h
#define Cosmic_simPlasticHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"

namespace Cosmic_sim
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
    PlasticHit(const PlasticHit&) = default;
    ~PlasticHit() override;

    // operators
    PlasticHit& operator=(const PlasticHit&) = default;
    G4bool operator==(const PlasticHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    void Draw()  override{}
    void Print() override;

    // methods to handle data
    void Add(G4double de, G4double dl);

    // get methods
    G4double GetEdep() const;
    G4double GetTrackLength() const;

  private:
    G4double fEdep = 0.;        ///< Energy deposit in the sensitive volume
    G4double fTrackLength = 0.; ///< Track length in the  sensitive volume
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using PlasticHitsCollection = G4THitsCollection<PlasticHit>;

extern G4ThreadLocal G4Allocator<PlasticHit>* PlasticHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* PlasticHit::operator new(size_t)
{
  if (!PlasticHitAllocator) {
      PlasticHitAllocator = new G4Allocator<PlasticHit>;
  }
  void *hit;
  hit = (void *) PlasticHitAllocator->MallocSingle();
  return hit;
}

inline void PlasticHit::operator delete(void *hit)
{
  if (!PlasticHitAllocator) {
      PlasticHitAllocator = new G4Allocator<PlasticHit>;
  }
    PlasticHitAllocator->FreeSingle((PlasticHit*) hit);
}

inline void PlasticHit::Add(G4double de, G4double dl) {
  fEdep += de;
  fTrackLength += dl;
}

inline G4double PlasticHit::GetEdep() const {
  return fEdep;
}

inline G4double PlasticHit::GetTrackLength() const {
  return fTrackLength;
}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
