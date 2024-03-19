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
/// \file ChamberSD.cc
/// \brief Implementation of the Cosmic::ChamberSD class

#include "ChamberSD.hh"
#include "DetectorConstruction.hh"
#include "ChamberHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"

#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4Box.hh"



namespace Cosmic
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  //  ChamberSD::ChamberSD(const G4String &name, const G4String &hitsCollectionName, G4int nofCells)
    //ChamberSD::ChamberSD(const G4String &name,const G4String &hitsCollectionName, DetectorConstruction* det)
    ChamberSD::ChamberSD(const G4String &name,const G4String &hitsCollectionName, G4int nsystem )
    :G4VSensitiveDetector(name), fNSystem(nsystem)
  {
        collectionName.insert(hitsCollectionName);
        //HitID = new G4int[NHITS];

    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChamberSD::~ChamberSD()
{
   // delete [] HitID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChamberSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection
    = new ChamberHitsCollection(SensitiveDetectorName, collectionName[0]);

    // Add this collection in hce
    auto hcID
            = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    hce->AddHitsCollection( hcID, fHitsCollection );


    // for(G4int j=0;j<NHITS;j++) HitID[j] = -1;

  //  hce=NULL;

  // Add this collection in hce
   // if (fHCID<0) {
   //     fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  //  }
  //  hce->AddHitsCollection( fHCID, fHitsCollection );

  // Create hits
  // fNofCells for cells + one more for total sums

    G4int fNofCells = (NW1_WRS + NW2_WRS + NW3_WRS + NVC_WRS) + 100;
   // G4cerr << " fNofCells = " << fNofCells << G4endl;

    for (G4int i = 0; i < fNofCells + 1; i++) {
        fHitsCollection->insert(new ChamberHit(0, 0, 0));
    }


 /*
    // fill calorimeter hits with zero energy deposition
    for (auto yID = 0; yID < N_LAYERS + 1; yID++) {
        for (auto xID = 0; xID < NX_BARS + 1; xID++) {
            fHitsCollection->insert(new ChamberHit(xID, yID, 0));
        }
        for (auto zID = 0; zID < NZ_BARS + 1; zID++) {
            fHitsCollection->insert(new ChamberHit(0, yID, zID));
        }

    }

  */

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool ChamberSD::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{

   ChamberHit* hit;
  // energy deposit
  auto edep = aStep->GetTotalEnergyDeposit();
  if (edep == 0.) return true;

  G4Track *aTrack=aStep->GetTrack();
  auto preStepPoint = aStep->GetPreStepPoint();
  auto postStepPoint = aStep->GetPostStepPoint();
 // G4cout<<"Position:" << preStepPoint->GetPosition()<<G4endl;

//  aTrack->SetTrackStatus(fStopAndKill);


    auto tof = aTrack->GetGlobalTime(); // время
    auto posit = aTrack->GetPosition(); // позиция
    G4ThreeVector posit_local(0, 0, 0);
    auto dx = aStep->GetDeltaPosition();
    auto velosity = preStepPoint->GetVelocity(); // скорость

    // step length
    G4double stepLength = 0.;
    if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0.) {
        stepLength = aStep->GetStepLength();
    }

   G4int TrackID = aTrack->GetTrackID();

   G4bool did = (TrackID == 0) or (TrackID == 1); // position only for original particle
   //  G4bool did = (TrackID ==0); // position only for original particle

   // if ( edep==0. && stepLength == 0. ) return false;

   auto mass = preStepPoint->GetMass();
   auto kinetic_energy = preStepPoint->GetKineticEnergy();
    //auto theta = preStepPoint->GetPosition().getTheta();
    //auto phi = preStepPoint->GetPosition().getPhi();

    G4double theta = 0., phi = 0.;

   //  if (did) theta = aTrack->GetPosition().getTheta(); // position only for original particle
   //  if (did) phi = aTrack->GetPosition().getPhi(); // position only for original particle


   // theta = aTrack->GetPosition().getTheta();
   //phi = aTrack->GetPosition().getPhi();
   // theta = aTrack->GetMomentum().getTheta();
   //  phi = aTrack->GetMomentum().getPhi();

   theta = aTrack->GetMomentumDirection().getTheta();
   phi = aTrack->GetMomentumDirection().getPhi();

   auto touchable = preStepPoint->GetTouchable();
    auto copyNo = touchable->GetCopyNumber();

    auto motherPhysical = touchable->GetVolume(1); // mother
    auto copyNo_mother = motherPhysical->GetCopyNo();

    auto physicalVol = touchable->GetVolume();
    auto logicalVol = preStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
    auto copyNo_phys = physicalVol->GetCopyNo();
    auto positionDetector = physicalVol->GetTranslation();
    auto rotation = touchable->GetRotation();


    auto box = (G4Box *) touchable->GetSolid();
    auto halflength = G4ThreeVector(box->GetXHalfLength(), box->GetYHalfLength(), box->GetZHalfLength());


    auto Vpos = touchable->GetTranslation();
    posit_local = posit - Vpos;
    posit_local.transform(*rotation);

    // if(did) theta = Vpos.getTheta();
    //if(did) phi = Vpos.getPhi();

    // G4cout<<"Detector position= "<< positionDetector<<G4endl;
    // G4cout<<"Detector position= "<< positionDetector<<G4endl;

    // Get calorimeter cell id
    auto layerNumber = touchable->GetReplicaNumber(1);
    auto hitID = layerNumber; // начинается с нуля
    auto hitTime = preStepPoint->GetGlobalTime();

    auto evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

   // G4cout <<"event_number= "<< evt <<" layerNumber= " << layerNumber << " copyNo_mother= " <<copyNo_mother<<" copyNo_phys= "<<copyNo_phys<< " copyNo="<< copyNo<< G4endl;


    G4int CB=0;
    G4int ni = touchable->GetHistoryDepth();
    for(G4int i=0;i<ni;i++){		// determines element label
        G4int k=touchable->GetReplicaNumber(i);
        if(k>0) CB += k;
    }


   // auto rowNo = touchable->GetCopyNumber(2);
  //  auto columnNo = touchable->GetCopyNumber(3);
  //  auto hitID = kNofHadRows*columnNo+rowNo;

   // G4cerr  <<" layerNumber= " << touchable->GetReplicaNumber(1) <<  " copyNo="<< touchable->GetCopyNumber(1)<< G4endl;
   // G4cerr << " CB= " << CB<<G4endl;

    //G4cerr  <<" layerNumber0= " << touchable->GetReplicaNumber(0) <<" layerNumber1= " << touchable->GetReplicaNumber(1)<<" layerNumber2= " << touchable->GetReplicaNumber(2)<<" layerNumber3= " << touchable->GetReplicaNumber(3)<<" layerNumber4= " << touchable->GetReplicaNumber(4)<<G4endl;
   // G4cerr << " CB= " << CB - ARM2_IND<<G4endl;


    hitID = CB;

    hit = (*fHitsCollection)[hitID];


    if ( ! hit ) {
        G4ExceptionDescription msg;
        msg << "Cannot access hit " << layerNumber;
        G4Exception("ChamberSD::ProcessHits()",
                    "MyCode0004", FatalException, msg);
    }


    if (!(hit->GetLogV())) {
        // fill volume information
        hit->SetLogV(physicalVol->GetLogicalVolume());
        G4AffineTransform transform = touchable->GetHistory()->GetTopTransform();
        transform.Invert();
        hit->SetRot(transform.NetRotation());
        hit->SetLocalPos(transform.NetTranslation());
    }



    // check if it is first touch
    //if (hit->GetColumnID()<0) {
     //   hit->SetColumnID(columnNo);
     //   hit->SetRowID(rowNo);
     //   auto depth = touchable->GetHistory()->GetDepth();
      //  auto transform = touchable->GetHistory()->GetTransform(depth-2);
      //  transform.Invert();
     //   hit->SetRot(transform.NetRotation());
      //  hit->SetPos(transform.NetTranslation());
 //   }

    // Get hit for total accounting

    auto hitTotal
            = (*fHitsCollection)[fHitsCollection->entries()-1];


    // add energy deposition
    // Add values


   hit->AddTrackID(TrackID);
   hit->AddEdep(edep);
   hit->AddLO(edep, posit, dx, velosity, tof);
   hit->AddTrackLength(stepLength);
   hit->AddToF(tof);
   hit->AddWorldPos(posit);
   hit->AddLocalPos(posit_local);
   hit->SetHalfLength(halflength);
   hit->SetPosTheta(theta);
   hit->SetPosPhi(phi);
   hit->SetMass(mass);
   hit->SetKineticEnergy(kinetic_energy);
   hit->SetBlkN(CB);

   hitTotal->AddTrackID(TrackID);
   hitTotal->AddEdep(edep);
   hitTotal->AddLO(edep, posit, dx, velosity, tof);
   hitTotal->AddTrackLength(stepLength);
   hitTotal->AddToF(tof);
   hitTotal->AddWorldPos(posit);
   hitTotal->AddLocalPos(posit_local);
   hitTotal->SetHalfLength(halflength);
   hitTotal->SetPosTheta(theta);
   hitTotal->SetPosPhi(phi);
   hitTotal->SetMass(mass);
   hitTotal->SetKineticEnergy(kinetic_energy);
   hitTotal->SetBlkN(CB);


   ROhist = NULL;
   return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChamberSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) {
     auto nofHits = fHitsCollection->entries();
     G4cout
       << G4endl
       << "-------->Hits Collection: in this event they are " << nofHits
       << " hits in the tracker chambers: " << G4endl;
     for ( std::size_t i=0; i<nofHits; ++i ) (*fHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
