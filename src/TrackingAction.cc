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
/// \file optical/OpNovice2/src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
//

#include "TrackingAction.hh"

#include "UserTrackInformation.hh"

#include "G4Track.hh"
#include "G4TrackingManager.hh"


#include "Constants.hh"
#include "DetectorConstruction.hh"
#include "Trajectory.hh"

#include "G4OpticalPhoton.hh"
#include "G4Track.hh"
#include "G4TrackingManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
TrackingAction::TrackingAction()
        : G4UserTrackingAction() {}
*/


void TrackingAction::PreUserTrackingAction(const G4Track *aTrack)
{

    // Let this be up to the user via vis.mac
    //  fpTrackingManager->SetStoreTrajectory(true);

    // Use custom trajectory class
    fpTrackingManager->SetTrajectory(new Trajectory(aTrack));

    // This user track information is only relevant to the photons
    fpTrackingManager->SetUserTrackInformation(new UserTrackInformation);

    /*
    // Create trajectory only for track in tracking region
    TrackInformation *trackInfo =
            (TrackInformation *) (aTrack->GetUserInformation());

    if (!trackInfo) {
        trackInfo = new TrackInformation(aTrack);
        trackInfo->SetIsFirstTankX(true);
        aTrack->SetUserInformation(trackInfo);
    }

    trackInfo->SetIsFirstTankX(true);
    */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TrackingAction::PostUserTrackingAction(const G4Track *aTrack)
{


    // G4cerr << "!!!" << G4endl;
    /*
        G4TrackVector *secondaries = fpTrackingManager->GimmeSecondaries(); // только вторичные

        if (secondaries) {
            TrackInformation *info = (TrackInformation *) (aTrack->GetUserInformation());
            size_t nSeco = secondaries->size();
            if (nSeco > 0) {
                for (size_t i = 0; i < nSeco; ++i) {
                    TrackInformation *infoNew = new TrackInformation(info);
                    (*secondaries)[i]->SetUserInformation(infoNew);
                }
            }
        }
    */

    auto trajectory = (Trajectory *)fpTrackingManager->GimmeTrajectory();
    auto trackInformation = (UserTrackInformation *)aTrack->GetUserInformation();

    // Let's choose to draw only the photons that hit the sphere and a pmt
    if (aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
    {
        const G4VProcess *creator = aTrack->GetCreatorProcess();
        if (creator && creator->GetProcessName() == "OpWLS")
        {
            trajectory->WLS();
            trajectory->SetDrawTrajectory(true);
        }

        // trajectory->SetDrawTrajectory(true);
        // trajectory->SetDrawTrajectory_noeff(true); // рисование непопавших в фотокатод другим цветом

#ifdef isWLSFiber
        // trajectory->SetDrawTrajectory(true);
        // trajectory->SetDrawTrajectory_noeff(true); // рисование непопавших в фотокатод другим цветом
#endif

        if (1 /*DetectorConstruction::GetSphereOn()*/)
        {
            if (

                (trackInformation->GetTrackStatus()
#if defined(isRUN23) || defined(isRUN21)
                 & hitSiPM
#endif
#ifndef isWLSFiber
#if defined(isRUN25) || defined(isFEU187)
                 & hitPMT
#endif
#endif
#ifdef isWLSFiber
                 & hitSiPM
#endif

                 )

                && (trackInformation->GetTrackStatus()
#if defined(isRUN25) || defined(isFEU187)
                    & hitPMT
#endif
#ifdef isWLSFiber
                    & hitSiPM
#endif
                    )) // hitSphere
            {
                trajectory->SetDrawTrajectory(true);
                // trajectory->SetDrawTrajectory_noeff(true); // рисование непопавших в фотокатод другим цветом
                trajectory->SetDrawTrajectory_noeff(false);
            }
        }
        else
        {
            if (trackInformation->GetTrackStatus()

#if defined(isRUN23) || defined(isRUN21)
                & hitSiPM
#endif
#if defined(isRUN25) || defined(isFEU187)
                & hitPMT
#endif

            )
                trajectory->SetDrawTrajectory(true);
            // trajectory->SetDrawTrajectory(true);
            trajectory->SetDrawTrajectory_noeff(true); // рисование непопавших в фотокатод другим цветом
                                                       // 1;
        }

        trajectory->SetDrawTrajectory(true);
        // trajectory->SetDrawTrajectory_noeff(true); // рисование непопавших в фотокатод другим цветом
        trajectory->SetDrawTrajectory_noeff(false);
    }
    // draw all other (not optical photon) trajectories
    else
    {

        trajectory->SetDrawTrajectory(true);
    }


    if (trackInformation->GetForceDrawTrajectory())
    {
        trajectory->SetDrawTrajectory(true);
    }


    trajectory->SetDrawTrajectory(true); // все траектории
    trajectory->SetDrawTrajectory_noeff(false);
    // trajectory->SetDrawTrajectory(false);
}
