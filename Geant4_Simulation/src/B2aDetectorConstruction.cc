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
/// \file B2aDetectorConstruction.cc
/// \brief Implementation of the B2aDetectorConstruction class
 
#include "B2aDetectorConstruction.hh"
#include "B2aDetectorMessenger.hh"
#include "B2TrackerSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4ThreadLocal 
G4GlobalMagFieldMessenger* B2aDetectorConstruction::fMagFieldMessenger = 0;

B2aDetectorConstruction::B2aDetectorConstruction()
:G4VUserDetectorConstruction(), 
 fNbOfChambers(0),
 fLogicTarget(NULL), fLogicChamber(NULL), 
 fTargetMaterial(NULL), fChamberMaterial(NULL), 
 fStepLimit(NULL),
 fCheckOverlaps(true)
{

  fTargetThickness = 10.*um;
  ComputeGeomParameters();

  fMessenger = new B2aDetectorMessenger(this);

  fNbOfChambers = 1;
  fLogicChamber = new G4LogicalVolume*[fNbOfChambers];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
B2aDetectorConstruction::~B2aDetectorConstruction()
{
  delete [] fLogicChamber; 
  delete fStepLimit;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* B2aDetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::DefineMaterials()
{
  // Material definition 

  G4NistManager* nistManager = G4NistManager::Instance();

  // Space defined using NIST Manager
  nistManager->FindOrBuildMaterial("G4_Galactic");
  
  G4double ice_density = 0.917*g/cm3;  // 冰的密度
  G4Element* elH = nistManager->FindOrBuildElement("H");
  G4Element* elO = nistManager->FindOrBuildElement("O");
  G4Material* ice = new G4Material("Ice_Crystal", ice_density, 2);
  ice->AddElement(elH, 2);
  ice->AddElement(elO, 1);
  fTargetMaterial = ice;

  // 或者直接使用水并注释说这是简化模型
  // fTargetMaterial = nistManager->FindOrBuildMaterial("G4_WATER");

  // Meterial gas defined using NIST Manager
  fChamberMaterial = nistManager->FindOrBuildMaterial("G4_Galactic");

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::ComputeGeomParameters()
{
  if(nullptr != fPhysiWorld) { ChangeGeometry(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B2aDetectorConstruction::DefineVolumes()
{
  G4Material* galactic  = G4Material::GetMaterial("G4_Galactic");
  //     
  // World
  //
  G4double world_sizeXY = 200*um; //200*mm;
  G4double world_sizeZ  = 200*um; //200*mm;
  //
  // Target
  // Size of Target
  // G4double Target_Z = 18.66*0.001/19.32*cm;
  G4double Target_Z = fTargetThickness; //17.72*um;
  G4double Target_XY = fTargetThickness; //1./4.*world_sizeXY;
  // //
  // // Tracker
  // //
  // G4ThreeVector TrackerPos = G4ThreeVector(0*cm, 0*cm, 0*cm);
  // G4double Tracker_XY = world_sizeXY;
  // G4double Tracker_Z = world_sizeZ/2.;

  // Chamber
  // Sphere shape
  G4double Chamber_MaxR = world_sizeXY/2.-10.*um; //world_sizeXY/2.-10.*mm;
  G4double Chamber_MinR = Chamber_MaxR-5.*um; //Chamber_MaxR-5.*mm;

  // Definitions of Solids, Logical Volumes, Physical Volumes

  // World

  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(world_sizeZ);

  G4cout << "Computed tolerance = "
         << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
         << " mm" << G4endl;

  G4Box* worldS
    = new G4Box("world",                                    //its name
                world_sizeXY/2,world_sizeXY/2,world_sizeZ/2); //its size
  G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                 worldS,   //its solid
                 galactic,      //its material
                 "World"); //its name
  
  //  Must place the World Physical volume unrotated at (0,0,0).
  // 
  fPhysiWorld
    = new G4PVPlacement(
                 0,               // no rotation
                 G4ThreeVector(), // at (0,0,0)
                 worldLV,         // its logical volume
                 "World",         // its name
                 0,               // its mother  volume
                 false,           // no boolean operations
                 0,               // copy number
                 fCheckOverlaps); // checking overlaps 

  // Target

  G4ThreeVector positionTarget = G4ThreeVector(0,0,-Target_Z/2);
  fSolidTarget
    = new G4Sphere("target", 0, Target_Z/2, 0.*deg, 360.*deg, 0.*deg, 180.*deg);
  
  fLogicTarget
    = new G4LogicalVolume(fSolidTarget, fTargetMaterial,"Target",0,0,0);
  fPhysiTarget = new G4PVPlacement(0,               // no rotation
                    positionTarget,  // at (x,y,z)
                    fLogicTarget,    // its logical volume
                    "Target",        // its name
                    worldLV,         // its mother volume
                    false,           // no boolean operations
                    0,               // copy number
                    fCheckOverlaps); // checking overlaps 

  G4cout << "Target is " << Target_Z/2/um << " um radius of "
         << fTargetMaterial->GetName() << G4endl;

  // // Tracker
 
  // G4ThreeVector positionTracker = G4ThreeVector(0,0,Tracker_Z/2);

  // G4Box* trackerS
  //   = new G4Box("tracker", Tracker_XY/2, Tracker_XY/2, Tracker_Z/2);
  // G4LogicalVolume* trackerLV
  //   = new G4LogicalVolume(trackerS, galactic, "Tracker",0,0,0);  
  // new G4PVPlacement(0,               // no rotation
  //                   positionTracker, // at (x,y,z)
  //                   trackerLV,       // its logical volume
  //                   "Tracker",       // its name
  //                   worldLV,         // its mother  volume
  //                   false,           // no boolean operations
  //                   0,               // copy number
  //                   fCheckOverlaps); // checking overlaps 

  // Visualization attributes

  G4VisAttributes* boxVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  G4VisAttributes* targetVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  G4VisAttributes* chamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0,0.2));

  worldLV      ->SetVisAttributes(boxVisAtt);
  fLogicTarget ->SetVisAttributes(targetVisAtt);
  // trackerLV    ->SetVisAttributes(boxVisAtt);

  // Tracker segments

  G4cout << "There are " << fNbOfChambers << " chambers in the tracker region. "
         << G4endl
         << "The chambers are of " << fChamberMaterial->GetName() << G4endl;

  for (G4int copyNo=0; copyNo<fNbOfChambers; copyNo++) {

      G4Sphere* chamberS
        = new G4Sphere("Chamber_solid", Chamber_MinR, Chamber_MaxR, 0.*deg, 360.*deg, 0.*deg, 180.*deg);

      fLogicChamber[copyNo] =
              new G4LogicalVolume(chamberS,fChamberMaterial,"Chamber_LV",0,0,0);

      fLogicChamber[copyNo]->SetVisAttributes(chamberVisAtt);

      new G4PVPlacement(0,                            // no rotation
                        G4ThreeVector(0,0,0), // at (x,y,z)
                        fLogicChamber[copyNo],        // its logical volume
                        "Chamber_PV",                 // its name
                        worldLV,                    // its mother  volume
                        false,                        // no boolean operations
                        copyNo,                       // copy number
                        fCheckOverlaps);              // checking overlaps 

  }

  // Example of User Limits
  //
  // Below is an example of how to set tracking constraints in a given
  // logical volume
  //
  // Sets a max step length in the tracker region, with G4StepLimiter

  // G4double maxStep = 0.5*(Chamber_MaxR-Chamber_MinR);
  // fStepLimit = new G4UserLimits(maxStep);
  // fLogicChamber[0]->SetUserLimits(fStepLimit);
 
  /// Set additional contraints on the track, with G4UserSpecialCuts
  ///
  /// G4double maxLength = 2*trackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
  /// trackerLV->SetUserLimits(new G4UserLimits(maxStep,
  ///                                           maxLength,
  ///                                           maxTime,
  ///                                           minEkin));

  // Always return the physical world

  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::PrintGeomParameters()
{
  G4cout << "Zefan - The height of the world is " 
         << G4BestUnit(fTargetThickness,"Length") << G4endl;  
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void B2aDetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  G4String trackerChamberSDname = "B2/TrackerChamberSD";
  B2TrackerSD* aTrackerSD = new B2TrackerSD(trackerChamberSDname,
                                            "TrackerHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  // Setting aTrackerSD to all logical volumes with the same name 
  // of "Chamber_LV".
  SetSensitiveDetector("Chamber_LV", aTrackerSD, true);

  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void B2aDetectorConstruction::SetTargetMaterial(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial = 
              nistManager->FindOrBuildMaterial(materialName);

  if (fTargetMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fTargetMaterial = pttoMaterial;
        if (fLogicTarget) fLogicTarget->SetMaterial(fTargetMaterial);
        G4cout 
          << G4endl 
          << "----> The target is made of " << materialName << G4endl;
     } else {
        G4cout 
          << G4endl 
          << "-->  WARNING from SetTargetMaterial : "
          << materialName << " not found" << G4endl;
     }
  }
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::SetChamberMaterial(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial =
              nistManager->FindOrBuildMaterial(materialName);

  if (fChamberMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fChamberMaterial = pttoMaterial;
        for (G4int copyNo=0; copyNo<fNbOfChambers; copyNo++) {
            if (fLogicChamber[copyNo]) fLogicChamber[copyNo]->
                                               SetMaterial(fChamberMaterial);
        }
        G4cout 
          << G4endl 
          << "----> The chambers are made of " << materialName << G4endl;
     } else {
        G4cout 
          << G4endl 
          << "-->  WARNING from SetChamberMaterial : "
          << materialName << " not found" << G4endl;
     }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::SetTargetThickness(G4double val)
{
  fTargetThickness = val;
  ComputeGeomParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
  fCheckOverlaps = checkOverlaps;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::ChangeGeometry()
{
  fSolidTarget->SetOuterRadius(fTargetThickness*0.5);
}
