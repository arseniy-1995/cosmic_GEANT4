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
/// \file PlasticSD.cc
/// \brief Implementation of the Cosmic::PlasticSD class

#include "PlasticSD.hh"
#include "DetectorConstruction.hh"
#include "PlasticHit.hh"

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

  //  PlasticSD::PlasticSD(const G4String &name, const G4String &hitsCollectionName, G4int nofCells)
    //PlasticSD::PlasticSD(const G4String &name,const G4String &hitsCollectionName, DetectorConstruction* det)
    PlasticSD::PlasticSD(const G4String &name,const G4String &hitsCollectionName, G4int nofLayers )
    :G4VSensitiveDetector(name), fNofLayers(nofLayers) // Это список инициализации
  {
        collectionName.insert(hitsCollectionName);
        //HitID = new G4int[NHITS];
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PlasticSD::~PlasticSD()
{
   // delete [] HitID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PlasticSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection
    = new PlasticHitsCollection(SensitiveDetectorName, collectionName[0]);

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
    fHitsCollection->insert(new PlasticHit(i));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool PlasticSD::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{

    PlasticHit *hit;
    PlasticHit *hitTotal;

    // energy deposit
    auto edep = aStep->GetTotalEnergyDeposit();
    if (edep == 0.) return true;

    G4Track *aTrack = aStep->GetTrack();


    G4int TrackID = aTrack->GetTrackID();

    auto preStepPoint = aStep->GetPreStepPoint();
    auto postStepPoint = aStep->GetPostStepPoint();
    // G4cout<<"Position:" << preStepPoint->GetPosition()<<G4endl;

//  aTrack->SetTrackStatus(fStopAndKill);

    auto tof = aTrack->GetGlobalTime(); // время
    auto posit = aTrack->GetPosition(); // позиция (вектор)
    G4ThreeVector posit_local(0, 0, 0);
  auto dx = aStep->GetDeltaPosition(); // смещение (вектор)
  auto velosity = preStepPoint->GetVelocity(); // скорость


    // auto theta = aTrack->GetPosition().getTheta();
    //auto phi = aTrack->GetPosition().getPhi();
    // auto theta = aTrack->GetMomentum().getTheta();
    //  auto phi = aTrack->GetMomentum().getPhi();

    auto theta = aTrack->GetMomentumDirection().getTheta();
    auto phi = aTrack->GetMomentumDirection().getPhi();

    auto thetaGlob = aTrack->GetPosition().getTheta();
    auto phiGlob = aTrack->GetPosition().getPhi();

    //   G4cout<< aTrack->GetMomentum().getTheta()/degree << " "<< aTrack->GetMomentumDirection().getTheta()/degree << " "<< aTrack->GetPosition().getTheta()/degree<<G4endl;


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
    auto rotation = touchable->GetRotation();


    auto box = (G4Box*) touchable->GetSolid();
    auto halflength = G4ThreeVector (box->GetXHalfLength(),box->GetYHalfLength(),box->GetZHalfLength());



    auto Vpos = touchable->GetTranslation();
    posit_local =posit-Vpos;
    posit_local.transform(*rotation);
    // G4cout<<"Detector position= "<< positionDetector<<G4endl;

    // Get calorimeter cell id
    auto layerNumber = touchable->GetReplicaNumber(1);
    auto hitID = layerNumber; // начинается с нуля
    auto hitTime = preStepPoint->GetGlobalTime();

    auto evt=G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

   // G4cout <<"event_number= "<< evt <<" layerNumber= " << layerNumber << " copyNo_mother= " <<copyNo_mother<<" copyNo_phys= "<<copyNo_phys<< " copyNo="<< copyNo<< G4endl;

  //  G4cerr  <<" layerNumber= " << touchable->GetReplicaNumber(1) <<  " copyNo="<< touchable->GetCopyNumber(1)<< G4endl;


    hit = (*fHitsCollection)[hitID];

    // Get hit for total accounting
    hitTotal = (*fHitsCollection)[fHitsCollection->entries()-1];

    if ( ! hit ) {
        G4ExceptionDescription msg;
        msg << "Cannot access hit " << layerNumber;
        G4Exception("PlasticSD::ProcessHits()",
                    "MyCode0004", FatalException, msg);
    }

    if (!(hit->GetLogV())) {
        // fill volume information
        hit->SetLogV(physicalVol->GetLogicalVolume());
        G4AffineTransform transform = touchable->GetHistory()->GetTopTransform();
        transform.Invert();
        hit->SetRot(transform.NetRotation());
        //hit->SetWorldPos(transform.NetTranslation());
      //  hitTotal->SetRot(transform.NetRotation());
       // hitTotal->SetWorldPos(transform.NetTranslation());
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
    hit->SetPosThetaGlob(thetaGlob);
    hit->SetPosPhiGlob(phiGlob);

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
    hitTotal->SetPosThetaGlob(thetaGlob);
    hitTotal->SetPosPhiGlob(phiGlob);

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
        hit= new PlasticHit(halflength,attenuation_length,discr_threshold,tv,rr);
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

void PlasticSD::EndOfEvent(G4HCofThisEvent*)
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
