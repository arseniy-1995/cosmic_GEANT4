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


//#define NHITS 20

namespace Cosmic
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  //  HadronCalorimeterSD::HadronCalorimeterSD(const G4String &name, const G4String &hitsCollectionName, G4int nofCells)
    //HadronCalorimeterSD::HadronCalorimeterSD(const G4String &name,const G4String &hitsCollectionName, DetectorConstruction* det)
    HadronCalorimeterSD::HadronCalorimeterSD(const G4String &name,const G4String &hitsCollectionName, G4int nofLayers )
    :G4VSensitiveDetector(name), fNofLayers(nofLayers)
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
  for (G4int i=0; i<fNofLayers+1; i++ ) {
    fHitsCollection->insert(new HadronCalorimeterHit(i));
  }
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
  auto posit=aTrack->GetPosition(); // позиция
  auto dx = aStep->GetDeltaPosition();


  // step length
  G4double stepLength = 0.;
  if ( aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    stepLength = aStep->GetStepLength();
  }

    // if ( edep==0. && stepLength == 0. ) return false;




    auto touchable = preStepPoint->GetTouchable();
    auto copyNo = touchable->GetCopyNumber();

    auto motherPhysical = touchable->GetVolume(1); // mother
    auto copyNo_mother = motherPhysical->GetCopyNo();

    auto physicalVol = touchable->GetVolume();
    auto logicalVol = preStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
    auto copyNo_phys = physicalVol->GetCopyNo();
    auto positionDetector = physicalVol->GetTranslation();

   // G4cout<<"Detector position= "<< positionDetector<<G4endl;

    // Get calorimeter cell id
    auto layerNumber = touchable->GetReplicaNumber(1);
    auto hitID = layerNumber; // начинается с нуля
    auto hitTime = preStepPoint->GetGlobalTime();

    auto evt=G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

   // G4cout <<"event_number= "<< evt <<" layerNumber= " << layerNumber << " copyNo_mother= " <<copyNo_mother<<" copyNo_phys= "<<copyNo_phys<< " copyNo="<< copyNo<< G4endl;

  //  G4cerr  <<" layerNumber= " << touchable->GetReplicaNumber(1) <<  " copyNo="<< touchable->GetCopyNumber(1)<< G4endl;


    hit = (*fHitsCollection)[hitID];


    if ( ! hit ) {
        G4ExceptionDescription msg;
        msg << "Cannot access hit " << layerNumber;
        G4Exception("HadronCalorimeterSD::ProcessHits()",
                    "MyCode0004", FatalException, msg);
    }

    if (!(hit->GetLogV())) {
        // fill volume information
        hit->SetLogV(physicalVol->GetLogicalVolume());
        G4AffineTransform transform = touchable->GetHistory()->GetTopTransform();
        transform.Invert();
        hit->SetRot(transform.NetRotation());
        hit->SetPos(transform.NetTranslation());
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
    hit->AddEdep(edep);
    hit->AddTrackLength(stepLength);
    hit->AddToF(tof);

    hitTotal->AddEdep(edep);
    hitTotal->AddTrackLength(stepLength);
    hitTotal->AddToF(tof);

/*
    G4int CB=0;
    G4int ni = touchable->GetHistoryDepth();
    for(G4int i=0;i<ni;i++){		// determines element label
        G4int k=touchable->GetReplicaNumber(i);
        if(k>0) CB += k;
    }
*/
/*

  G4cout<<"!!!!!"<<layerNumber<<G4endl;
 //   G4cout<<"eeee"<<layerNumber<<G4endl;
   // touchable->GetCopyNumber()
  // Get hit accounting data for this cell
  auto hit = (*fHitsCollection)[layerNumber];


  // Get hit for total accounting



 */

/*

    if (HitID[CB]==-1){
        G4Box *box = (G4Box*) touchable->GetSolid();
        G4double halflength = box->GetXHalfLength();
        G4double attenuation_length = 0.0, discr_threshold =0.0;
        //G4double attenuation_length = Detector->GetAttenuL(CB);
        //G4double discr_threshold = Detector->GetDiscrThr(CB);
        G4ThreeVector tv = touchable->GetTranslation();
        const G4RotationMatrix *rr = touchable->GetRotation();
        hit= new HadronCalorimeterHit(halflength,attenuation_length,discr_threshold,tv,rr);
        hit->SetBlkN(CB);
        hcHit->Add(edep, stepLength);
      //  hit->AddLO(edep, posit, dx);
      // hit->AddPos(edep,posit,prtn,tof);
        HitID[CB] = fHitsCollection->insert(hit) - 1;

//if((CB%ARM2_IND)== (HCZ_IND+5)) G4cout<<"=== CB="<<CB<<"  "<<CB/ARM2_IND<<"  edep="<<edep/keV<<G4endl;

    }else{
        hcHit = (*fHitsCollection)[HitID[CB]];
        G4ThreeVector tv = touchable->GetTranslation();
//	if(tv.y()==hcHit->GetVPos().y()){	// exclude similar cells in different layers of superlayer
        hcHit->Add(edep,stepLength);
      //  hcHit->AddLO(edep, posit, dx);
       // hcHit->AddPos(edep,posit,prtn,tof);
//	}
    }

    */

  //  fHitsCollection->insert(hit);
    ROhist=NULL;
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
