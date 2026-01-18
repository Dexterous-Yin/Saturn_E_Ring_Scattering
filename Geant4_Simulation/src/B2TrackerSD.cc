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
/// \file B2TrackerSD.cc
/// \brief Implementation of the B2TrackerSD class

#include "B2TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4TrackVector.hh"
#include "G4TrackStatus.hh"
#include "G4UnitsTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "B2SteppingAction.hh" // 添加SteppingAction的头文件

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2TrackerSD::B2TrackerSD(const G4String& name,
                         const G4String& hitsCollectionName) 
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2TrackerSD::~B2TrackerSD() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2TrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection 
    = new B2TrackerHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce

  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool B2TrackerSD::ProcessHits(G4Step* aStep, 
                                     G4TouchableHistory*)
{
  // parent ID
  G4int parent_id = aStep->GetTrack()->GetParentID();
  if (parent_id!=0) return false;

  G4String PostVolumeName = aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName();
  // G4cout << PostVolumeName << G4endl;
  if (std::strcmp(PostVolumeName,"World") != 0) return false;

  auto newHit = new B2TrackerHit();

  // newHit->SetParentID(aStep->GetTrack()->GetParentID());
  // newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
  // newHit->SetPos (aStep->GetPreStepPoint()->GetPosition());
  // newHit->SetP(aStep->GetPreStepPoint()->GetMomentum());
  // newHit->SetPdirection(aStep->GetPreStepPoint()->GetMomentumDirection());
  G4ThreeVector testP = aStep->GetPreStepPoint()->GetMomentumDirection();
  newHit->SetPdir(testP[2]);
  newHit->SetE(aStep->GetPreStepPoint()->GetKineticEnergy());
  
  // 获取当前粒子的trackID
  G4int trackID = aStep->GetTrack()->GetTrackID();
  
  // 检查电子是否经历了制动辐射过程，通过全局跟踪映射查询
  G4bool hasEBrem = B2SteppingAction::HasTrackExperiencedEBrem(trackID);
  // G4cout << "Zefan: stored process is " << hasEBrem << G4endl;
  G4bool hasEIoni = B2SteppingAction::HasTrackExperiencedEIoni(trackID);
  
  // 设置eBrem标志
  newHit->SetHasEBrem(hasEBrem);
  newHit->SetHasEIoni(hasEIoni);

  fHitsCollection->insert( newHit );
  
  // G4ThreeVector testP = aStep->GetPreStepPoint()->GetMomentumDirection();
  // if (testP[2]<0) {
  //   newHit->Print();
  // }
  // newHit->Print();

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2TrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) { 
     G4int nofHits = fHitsCollection->entries();
     G4cout << G4endl
            << "-------->Hits Collection: in this event they are " << nofHits 
            << " hits" << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
