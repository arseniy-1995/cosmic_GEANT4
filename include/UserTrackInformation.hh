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
/// \file optical/OpNovice2/include/TrackInformation.hh
/// \brief Definition of the TrackInformation class
//
//
//

#ifndef UserTrackInformation_h
#define UserTrackInformation_h 1

#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4VUserTrackInformation.hh"
#include "globals.hh"

#include "Constants.hh"

class G4Track;


enum TrackStatus // перечесления
{
    active = 1,
#if defined(isRUN23) || defined(isRUN21)
    hitSiPM = 2,
#endif


#ifdef isRUN25
    hitPMT = 2,
#endif
#ifdef isWLSFiber
    hitSiPM = 2,
#endif

#ifdef isFEU187
    hitPMT = 2,
#endif

#ifdef TRIGCOUNTsmall_MCP_PMT_BC_422
    hitPMT_TRIG = 2,
#endif

    absorbed = 4,
    boundaryAbsorbed = 8,
    hitSphere = 16,
    inactive = 14
};

/*TrackStatus:
  active: still being tracked
  hitSiPM: stopped by being detected in a SiPM
  absorbed: stopped by being absorbed with G4OpAbsorbtion
  boundaryAbsorbed: stopped by being aborbed with G4OpAbsorbtion
  hitSphere: track hit the sphere at some point
  inactive: track is stopped for some reason
   -This is the sum of all stopped flags so can be used to remove stopped flags

*/

enum TrackStatusWLS
{
    undefined = 0,
    left = 1, // track is going -z
    right = 2, // track is going +z
    defined = 3, // left or right flag is on (Can't be Set)
    EscapedFromSide = 4, // photon escaped through the side of the fiber
    EscapedFromReadOut = 8, // photon escaped through read-out end of the fiber
    ReflectedAtMirror = 16, // photon has been reflected by mirror at far end
    ReflectedAtReadOut = 32, // photon has been reflected at the read-out end
    murderee = 64, // photon is artificially killed
    InsideOfFiber = 128, // Flag is on if the photon is inside of fiber
    OutsideOfFiber = 256 // Flag is on if the photon is outside of fiber
};


class UserTrackInformation : public G4VUserTrackInformation
{
public:
    UserTrackInformation();
    UserTrackInformation(const G4Track *aTrack);
    UserTrackInformation(const UserTrackInformation *aTrackInfo);
    ~UserTrackInformation();

    inline void *operator new(size_t);
    inline void operator delete(void *aTrackInfo);

    UserTrackInformation &operator=(const UserTrackInformation &right);

    void SetSourceTrackInformation(const G4Track *aTrack);
    inline virtual void Print() const;

public:
    inline G4bool GetIsFirstTankX() const { return fFirstTankX; }
    inline void SetIsFirstTankX(G4bool b) { fFirstTankX = b; }


    // Sets the track status to s (does not check validity of flags)
    void SetTrackStatusFlags(int s) { fStatus = s; }
    // Does a smart add of track status flags (disabling old flags that conflict)
    // If s conflicts with itself it will not be detected
    void AddTrackStatusFlag(int s);

    int GetTrackStatus() const { return fStatus; }

    void IncReflections() { ++fReflections; }
    G4int GetReflectionCount() const { return fReflections; }

    void SetForceDrawTrajectory(G4bool b) { fForcedraw = b; }
    G4bool GetForceDrawTrajectory() { return fForcedraw; }

    /////////// for WLS

    const G4ThreeVector &GetExitPosition() const { return fExitPosition; }
    void SetExitPosition(const G4ThreeVector &pos) { fExitPosition = pos; }

    // Try adding a status flag and return if it is successful or not
    // Cannot Add Undefine or a flag that conflicts with another flag
    // Return true if the addition of flag is successful, false otherwise
    G4bool AddStatusFlagWLS(TrackStatusWLS s);

    // Check if a certain flag is on
    G4bool IsStatus(TrackStatusWLS s) { return s == undefined ? !(fStatus &= defined) : fStatus & s; }


private:
    G4bool fFirstTankX;
    G4int fStatus = active;
    G4int fReflections;
    G4bool fForcedraw = false;
    G4ThreeVector fExitPosition;
};

extern G4ThreadLocal G4Allocator<UserTrackInformation> *aTrackInformationAllocator;

inline void *UserTrackInformation::operator new(size_t)
{
    if (!aTrackInformationAllocator)
        aTrackInformationAllocator = new G4Allocator<UserTrackInformation>;
    return (void *)aTrackInformationAllocator->MallocSingle();
}

inline void UserTrackInformation::operator delete(void *aTrackInfo)
{
    aTrackInformationAllocator->FreeSingle((UserTrackInformation *)aTrackInfo);
}

#endif
