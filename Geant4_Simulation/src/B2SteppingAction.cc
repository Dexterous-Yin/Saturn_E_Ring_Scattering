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
/// \file B2SteppingAction.cc
/// \brief Implementation of the B2::SteppingAction class

#include "B2SteppingAction.hh"
#include "B2EventAction.hh"
#include "B2aDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4VProcess.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// 初始化静态线程局部变量
G4ThreadLocal std::map<G4int, G4bool>* B2SteppingAction::fTrackEBremMap = nullptr;
G4ThreadLocal std::map<G4int, G4bool>* B2SteppingAction::fTrackEIoniMap = nullptr;

B2SteppingAction::B2SteppingAction(B2EventAction* eventAction)
: fEventAction(eventAction)
{
  // 初始化fTrackEBremMap
  if (!fTrackEBremMap) {
    fTrackEBremMap = new std::map<G4int, G4bool>();
  }
    if (!fTrackEIoniMap) {
    fTrackEIoniMap = new std::map<G4int, G4bool>();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2SteppingAction::UserSteppingAction(const G4Step* step)
{
  // 获取步骤的物理过程
  const G4StepPoint* postStepPoint = step->GetPostStepPoint();
  const G4VProcess* process = postStepPoint ? postStepPoint->GetProcessDefinedStep() : nullptr;
  
  // 获取粒子ID
  G4int trackID = step->GetTrack()->GetTrackID();
  
  // 检查是否是物理过程
  if (process) {
    G4String processName = process->GetProcessName();
    
    // 检查是否是电子制动辐射过程 (eBrem)
    if (!HasTrackExperiencedEBrem(trackID) && 
        (processName == "eBrem" || processName == "eBrems" || processName == "eBremsstrahlung")) {
      if (!fTrackEBremMap) {
        fTrackEBremMap = new std::map<G4int, G4bool>();
      }
      // 记录该粒子经历了eBrem过程
      (*fTrackEBremMap)[trackID] = true;
      
      // 输出调试信息
      G4cout << "Track ID " << trackID 
             << " experienced eBrem process at energy " 
             << step->GetPreStepPoint()->GetKineticEnergy()/MeV 
             << " MeV" << G4endl;
    }
    
    // 检查是否是电子电离过程 (eIoni)
    if (!HasTrackExperiencedEIoni(trackID) && 
        (processName == "eIoni" || processName == "electron ionisation")) {
      if (!fTrackEIoniMap) {
        fTrackEIoniMap = new std::map<G4int, G4bool>();
      }
      // 记录该粒子经历了eIoni过程
      (*fTrackEIoniMap)[trackID] = true;
      
      // 输出调试信息
      G4cout << "Track ID " << trackID 
             << " experienced eIoni process at energy " 
             << step->GetPreStepPoint()->GetKineticEnergy()/MeV 
             << " MeV" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2SteppingAction::ClearTrackEBremMap()
{
  // 清空map，通常在每个事件开始前调用
  if (fTrackEBremMap) {
    fTrackEBremMap->clear();
  }
}

void B2SteppingAction::ClearTrackEIoniMap()
{
  // 清空map，通常在每个事件开始前调用
  if (fTrackEIoniMap) {
    fTrackEIoniMap->clear();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool B2SteppingAction::HasTrackExperiencedEBrem(G4int trackID)
{
  // 检查指定trackID是否经历了eBrem过程
  if (!fTrackEBremMap) {
    return false;
  }
  
  auto it = fTrackEBremMap->find(trackID);
  return (it != fTrackEBremMap->end() && it->second);
}

G4bool B2SteppingAction::HasTrackExperiencedEIoni(G4int trackID)
{
  // 检查指定trackID是否经历了eIoni过程
  if (!fTrackEIoniMap) {
    return false;
  }
  
  auto it = fTrackEIoniMap->find(trackID);
  return (it != fTrackEIoniMap->end() && it->second);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
