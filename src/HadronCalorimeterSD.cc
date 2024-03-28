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
/// \file HadronCalorimeterSD.cc
/// \brief Implementation of the Cosmic::HadronCalorimeterSD class

#include "HadronCalorimeterSD.hh"
#include "DetectorConstruction.hh"
#include "HadronCalorimeterHit.hh"

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

  //  HadronCalorimeterSD::HadronCalorimeterSD(const G4String &name, const G4String &hitsCollectionName, G4int nofCells)
    //HadronCalorimeterSD::HadronCalorimeterSD(const G4String &name,const G4String &hitsCollectionName, DetectorConstruction* det)
    HadronCalorimeterSD::HadronCalorimeterSD(const G4String &name,const G4String &hitsCollectionName, G4int nsystem )
    :G4VSensitiveDetector(name), fNSystem(nsystem)
  {
        collectionName.insert(hitsCollectionName);
        //HitID = new G4int[NHITS];

    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HadronCalorimeterSD::~HadronCalorimeterSD()
{
   // delete [] HitID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HadronCalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection
    = new HadronCalorimeterHitsCollection(SensitiveDetectorName, collectionName[0]);

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

  G4int fNofCells = (N_HCX+N_HCZ) + 100;

  for (G4int i=0; i<fNofCells+1; i++ ) {
    fHitsCollection->insert(new HadronCalorimeterHit(0,0,0));
  }


 /*
    // fill calorimeter hits with zero energy deposition
    for (auto yID = 0; yID < N_LAYERS + 1; yID++) {
        for (auto xID = 0; xID < NX_BARS + 1; xID++) {
            fHitsCollection->insert(new HadronCalorimeterHit(xID, yID, 0));
        }
        for (auto zID = 0; zID < NZ_BARS + 1; zID++) {
            fHitsCollection->insert(new HadronCalorimeterHit(0, yID, zID));
        }

    }

  */

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool HadronCalorimeterSD::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{

   HadronCalorimeterHit* hit;
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
  if ( aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    stepLength = aStep->GetStepLength();
  }

   G4int TrackID = aTrack->GetTrackID();

   // if ( edep==0. && stepLength == 0. ) return false;


   // auto theta = aTrack->GetPosition().getTheta();
   // auto phi = aTrack->GetPosition().getPhi();
   // auto theta = aTrack->GetMomentum().getTheta();
   // auto phi = aTrack->GetMomentum().getPhi();

    auto theta = aTrack->GetMomentumDirection().getTheta();
    auto phi = aTrack->GetMomentumDirection().getPhi();


   auto touchable = preStepPoint->GetTouchable();
   G4TouchableHistory* theTouchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
   auto copyNo = touchable->GetCopyNumber();

   //  G4cout << copyNo <<G4endl;

   auto motherPhysical = touchable->GetVolume(1); // mother
   auto copyNo_mother = motherPhysical->GetCopyNo();

   //  G4cout << copyNo_mother <<G4endl;

   // G4cout << copyNo_mother + copyNo <<G4endl;

   auto physicalVol = touchable->GetVolume();
   auto logicalVol = preStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
   auto copyNo_phys = physicalVol->GetCopyNo();
   auto positionDetector = physicalVol->GetTranslation();
   auto rotation = touchable->GetRotation();


   auto box = (G4Box*)touchable->GetSolid();
   auto halflength = G4ThreeVector(box->GetXHalfLength(), box->GetYHalfLength(), box->GetZHalfLength());


   auto Vpos = touchable->GetTranslation(); // положение центра детектора
   posit_local = posit - Vpos;
   posit_local.transform(*rotation);
   // G4cout<<"Detector position= "<< positionDetector<<G4endl;
   // G4cout<<"Detector position= "<< positionDetector<<G4endl;

    // Get calorimeter cell id
    auto layerNumber = touchable->GetReplicaNumber(1);
    auto hitID = layerNumber; // начинается с нуля
    auto hitTime = preStepPoint->GetGlobalTime();

    auto evt=G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

   // G4cout <<"event_number= "<< evt <<" layerNumber= " << layerNumber << " copyNo_mother= " <<copyNo_mother<<" copyNo_phys= "<<copyNo_phys<< " copyNo="<< copyNo<< G4endl;

   // G4cout << copyNo_mother + copyNo << " " << hitID <<G4endl;

   G4int CB = 0;
   G4int ni = theTouchable->GetHistoryDepth();
   for (G4int i = 0; i < ni; i++)
   {
       // determines element label
       G4int k = theTouchable->GetReplicaNumber(i);
       if (k > 0) CB += k;
   }


   // auto rowNo = touchable->GetCopyNumber(2);
  //  auto columnNo = touchable->GetCopyNumber(3);
  //  auto hitID = kNofHadRows*columnNo+rowNo;

  //  G4cerr  <<" layerNumber= " << touchable->GetReplicaNumber(1) <<  " copyNo="<< touchable->GetCopyNumber(1)<< G4endl;
    //G4cerr << " CB= " << CB<<G4endl;

    //G4cerr  <<" layerNumber0= " << touchable->GetReplicaNumber(0) <<" layerNumber1= " << touchable->GetReplicaNumber(1)<<" layerNumber2= " << touchable->GetReplicaNumber(2)<<" layerNumber3= " << touchable->GetReplicaNumber(3)<<" layerNumber4= " << touchable->GetReplicaNumber(4)<<G4endl;
   // G4cerr << " CB= " << CB - ARM2_IND<<G4endl;




    hitID = CB;

   G4int number_X = 0, number_Z = 0;
   G4int index = hitID;

   index = index - HCX_IND;
   if (index >= 0 && index < N_HCX)
   {
       number_X = index % NX_BARS;
       number_X = (index % N_UNITS) * 2 + (index / N_UNITS); // re-numbering bars in a layer
   }

   index = index - HCZ_IND;
   if (index >= 0 && index < N_HCZ)
   {
       number_Z = index % NZ_BARS;
       number_Z = (number_Z % N_UNITS) * 2 + (number_Z / N_UNITS); // re-numbering bars in a layer
   }

   //   G4cout << copyNo_mother + copyNo << " " << hitID <<G4endl;

   //     G4cout << hitID << " " << number_X << " " << number_Z <<G4endl;

   hit = (*fHitsCollection)[hitID];


   if (!hit)
   {
       G4ExceptionDescription msg;
       msg << "Cannot access hit " << layerNumber;
       G4Exception("HadronCalorimeterSD::ProcessHits()",
                   "MyCode0004", FatalException, msg);
   }

   if (!(hit->GetLogV()))
   {
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

    //   G4cout << CB <<G4endl;

    hit->SetHalfLength(halflength);
   hit->AddTrackID(TrackID);
   hit->AddEdep(edep);
   hit->AddLO(edep, posit, dx, velosity, tof);
   hit->AddTrackLength(stepLength);
   hit->AddToF(tof);
   hit->AddWorldPos(posit);
   hit->AddLocalPos(posit_local);
   hit->SetBlkN(CB);
   hit->SetVPos(Vpos);

   hitTotal->SetHalfLength(halflength);
   hitTotal->AddTrackID(TrackID);
   hitTotal->AddEdep(edep);
   hitTotal->AddLO(edep, posit, dx, velosity, tof);
   hitTotal->AddTrackLength(stepLength);
   hitTotal->AddToF(tof);
   hitTotal->AddWorldPos(posit);
   hitTotal->AddLocalPos(posit_local);
   hitTotal->SetBlkN(CB);
   hitTotal->SetVPos(Vpos);


   ROhist = NULL;
   return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HadronCalorimeterSD::EndOfEvent(G4HCofThisEvent*)
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
