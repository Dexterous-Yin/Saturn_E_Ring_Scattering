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
/// \file B2EventAction.cc
/// \brief Implementation of the B2EventAction class

#include "B2EventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "HistoManager.hh"
#include "B2TrackerHit.hh"
#include "B2RunAction.hh"
#include "B2SteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2EventAction::B2EventAction(B2RunAction* run, HistoManager* histo)
: G4UserEventAction(),
  fRunAct(run), fHistoManager(histo)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2EventAction::~B2EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2EventAction::BeginOfEventAction(const G4Event*)
{
  // 清除SteppingAction中的trackID-eBrem映射，为新事件做准备
  B2SteppingAction::ClearTrackEBremMap();
  B2SteppingAction::ClearTrackEIoniMap();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2EventAction::EndOfEventAction(const G4Event* event)
{

  // // periodic printing
  // G4int eventID = event->GetEventID();
  // if ( eventID < 100 || eventID % 100 == 0) {
  //   G4cout << ">>> Event: " << eventID  << G4endl;
  //   G4VHitsCollection* hc = event->GetHCofThisEvent()->GetHC(0);
  //   G4cout << "    "  
  //          << hc->GetSize() << " hits stored in this event" << G4endl;
  // }
  
  // Result Output

  G4VHitsCollection* hc = event->GetHCofThisEvent()->GetHC(0);
  G4int hc_num = hc->GetSize();
  if (hc_num>1) {
    G4cout << "Zefan: There are more than one hits in this event" << G4endl;
    G4cout << hc_num << " hits stored in this event" << G4endl;
  }
  if (hc_num==0){
    G4cout << "Zefan: Here" << G4endl;
    G4cout << hc_num << " hits stored in this event" << G4endl;
    return;
  }
  
  // auto hit = static_cast<B2TrackerHit*>(hc->GetHit(0));
  B2TrackerHitsCollection* myhc = (B2TrackerHitsCollection*)hc;
  auto hit = (*myhc)[0];
  // G4int ParentID = hit->GetParentID();
  // G4int TrackID = hit->GetTrackID();
  // G4ThreeVector Position = hit->GetPos();
  // G4ThreeVector P = hit->GetP();
  // G4ThreeVector Pdir = hit->GetPdirection();
  G4double Pdir = hit->GetPdir();
  G4double E = hit->GetE();
  G4bool hasEBrem = hit->GetHasEBrem();
  G4bool hasEIoni = hit->GetHasEIoni();
  
  //fill histograms
  //
  fHistoManager->FillHisto(0, E);
  
  //fill ntuple
  //
  // fHistoManager->FillNtuple(ParentID, TrackID, Position[0], Position[1], Position[2], P[0], P[1], P[2], Pdir[0], Pdir[1], Pdir[2], E);
  fHistoManager->FillNtuple(Pdir, E, hasEBrem, hasEIoni);
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
