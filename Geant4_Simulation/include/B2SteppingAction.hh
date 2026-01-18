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
/// \file B2/include/B2SteppingAction.hh
/// \brief Definition of the B2::SteppingAction class

#ifndef B2SteppingAction_h
#define B2SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <map>

class G4LogicalVolume;

class B2EventAction;

/// Stepping action class

class B2SteppingAction : public G4UserSteppingAction
{
  public:
    B2SteppingAction(B2EventAction* eventAction);
    ~B2SteppingAction() override = default;

    // method from the base class
    void UserSteppingAction(const G4Step*) override;
    
    // 清除事件开始时的静态数据
    static void ClearTrackEBremMap();
    static void ClearTrackEIoniMap();
    
    // 检查粒子是否经历了eBrem过程
    static G4bool HasTrackExperiencedEBrem(G4int trackID);
    static G4bool HasTrackExperiencedEIoni(G4int trackID);
    
  private:
    B2EventAction* fEventAction = nullptr;
    G4LogicalVolume* fScoringVolume = nullptr;
    
    // 使用G4ThreadLocal确保线程安全
    static G4ThreadLocal std::map<G4int, G4bool>* fTrackEBremMap;
    static G4ThreadLocal std::map<G4int, G4bool>* fTrackEIoniMap;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
