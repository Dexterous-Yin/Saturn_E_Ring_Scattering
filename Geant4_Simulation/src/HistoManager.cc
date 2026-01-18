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
/// \file analysis/AnaEx01/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
{
  // Create or get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  if ( ! fFactoryOn ) {
    //
    analysisManager->SetVerboseLevel(1);
    // Only merge in MT mode to avoid warning when running in Sequential mode
  #ifdef G4MULTITHREADED
    analysisManager->SetNtupleMerging(true);
  #endif

    // Create directories
    analysisManager->SetHistoDirectoryName("histo");
    analysisManager->SetNtupleDirectoryName("ntuple");
  }

  // Open an output file
  //
  G4bool fileOpen = analysisManager->OpenFile("Record");
  if (! fileOpen) {
    G4cerr << "\n---> HistoManager::Book(): cannot open "
           << analysisManager->GetFileName() << G4endl;
    return;
  }

  if ( ! fFactoryOn ) {
    // Create histograms.
    // Histogram ids are generated automatically starting from 0.
    // The start value can be changed by:
    // analysisManager->SetFirstHistoId(1);

    // id = 0
    analysisManager->CreateH1("E","Final Energy (MeV)", 10, 0*MeV, 10*MeV);

    // Create ntuples.
    // Ntuples ids are generated automatically starting from 0.
    // The start value can be changed by:
    // analysisManager->SetFirstMtupleId(1);

    // Create 1st ntuple (id = 0)
    analysisManager->CreateNtuple("Ntuple1", "Result");
    // analysisManager->CreateNtupleIColumn("ParentID"); // column Id = 0
    // analysisManager->CreateNtupleIColumn("TrackID"); // column Id = 1
    // analysisManager->CreateNtupleDColumn("X"); // column Id = 2
    // analysisManager->CreateNtupleDColumn("Y"); // column Id = 3
    // analysisManager->CreateNtupleDColumn("Z"); // column Id = 4
    // analysisManager->CreateNtupleDColumn("Px"); // column Id = 5
    // analysisManager->CreateNtupleDColumn("Py"); // column Id = 6
    // analysisManager->CreateNtupleDColumn("Pz"); // column Id = 7
    // analysisManager->CreateNtupleDColumn("Pdirx"); // column Id = 8
    // analysisManager->CreateNtupleDColumn("Pdiry"); // column Id = 9
    // analysisManager->CreateNtupleDColumn("Pdirz"); // column Id = 10
    analysisManager->CreateNtupleDColumn("Pdir"); // column Id = 0
    analysisManager->CreateNtupleDColumn("E"); // column Id = 1
    analysisManager->CreateNtupleIColumn("HasEBrem"); // column Id = 2, 布尔值以整数形式存储
    analysisManager->CreateNtupleIColumn("HasEIoni"); // column Id = 3, 布尔值以整数形式存储
    analysisManager->FinishNtuple();

    fFactoryOn = true;
  }

  G4cout << "\n----> Output file is open in "
         << analysisManager->GetFileName() << "."
         << analysisManager->GetFileType() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Save()
{
  if (! fFactoryOn) { return; }

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();

  G4cout << "\n----> Histograms and ntuples are saved\n" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillHisto(G4int ih, G4double xbin, G4double weight)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillH1(ih, xbin, weight);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Normalize(G4int ih, G4double fac)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  auto h1 = analysisManager->GetH1(ih);
  if (h1 != nullptr) { h1->scale(fac);
}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillNtuple(G4double Pdir, G4double E, G4bool hasEBrem, G4bool hasEIoni)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // Fill 1st ntuple ( id = 0)
  // analysisManager->FillNtupleIColumn(0, 0, ParentID);
  // analysisManager->FillNtupleIColumn(0, 1, TrackID);
  // analysisManager->FillNtupleDColumn(0, 2, X);
  // analysisManager->FillNtupleDColumn(0, 3, Y);
  // analysisManager->FillNtupleDColumn(0, 4, Z);
  // analysisManager->FillNtupleDColumn(0, 5, Px);
  // analysisManager->FillNtupleDColumn(0, 6, Py);
  // analysisManager->FillNtupleDColumn(0, 7, Pz);
  // analysisManager->FillNtupleDColumn(0, 8, Pdirx);
  // analysisManager->FillNtupleDColumn(0, 9, Pdiry);
  // analysisManager->FillNtupleDColumn(0, 10, Pdirz);
  // analysisManager->FillNtupleDColumn(0, 11, E);
  analysisManager->FillNtupleDColumn(0, 0, Pdir);
  analysisManager->FillNtupleDColumn(0, 1, E);
  analysisManager->FillNtupleIColumn(0, 2, hasEBrem ? 1 : 0); // 将布尔值转换为整数存储
  analysisManager->FillNtupleIColumn(0, 3, hasEIoni ? 1 : 0); // 将布尔值转换为整数存储
  analysisManager->AddNtupleRow(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::PrintStatistic()
{
  if (! fFactoryOn) { return; }

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  G4cout << "\n ----> print histograms statistic \n" << G4endl;
  for ( G4int i=0; i<analysisManager->GetNofH1s(); ++i ) {
    G4String name = analysisManager->GetH1Name(i);
    auto h1 = analysisManager->GetH1(i);

    G4String unitCategory;
    if (name[0U] == 'E' ) { unitCategory = "Energy"; }
    if (name[0U] == 'L' ) { unitCategory = "Length"; }
         // we use an explicit unsigned int type for operator [] argument
         // to avoid problems with windows compiler

    G4cout << name
           << ": mean = " << G4BestUnit(h1->mean(), unitCategory)
           << " rms = " << G4BestUnit(h1->rms(), unitCategory )
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
