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
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
//#include "G4GlobalMagFieldMessenger.hh"
#include "G4MagneticField.hh"
#include "G4AutoDelete.hh"
#include "G4ThreeVector.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4DecayTable.hh"
#include "G4NistManager.hh"
#include "G4Sphere.hh"
#include "G4Decay.hh"
#include "G4RadioactiveDecay.hh"
#include "G4UniformMagField.hh"
#include "G4SystemOfUnits.hh"
#include "G4FieldManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4ThreadLocal 
G4GlobalMagFieldMessenger* B2aDetectorConstruction::fMagFieldMessenger = 0;

//G4double fMagneticFieldStrength; //追加した
//G4ThreeVector fMagneticFieldDirection; //追加した

B2aDetectorConstruction::B2aDetectorConstruction()
:G4VUserDetectorConstruction(), 
 fNbOfChambers(0),
 fLogicTarget(NULL), fLogicChamber(NULL), 
 fTargetMaterial(NULL), fChamberMaterial(NULL), 
 fStepLimit(NULL),
 fCheckOverlaps(true)
 //fMagneticFieldStrength(7.0*tesla),//追加した
 //fMagneticFieldDirection(1.0, 0.0, 0.0)//追加した
{
  fMessenger = new B2aDetectorMessenger(this);

  fNbOfChambers = 5;
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
  
  //G4ThreeVector fieldValue = fMagneticFieldStrength * fMagneticFieldDirection;//追加した
  //fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);//追加した
  //fMagFieldMessenger->SetVerboseLevel(1);//追加した
  

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// ...



void B2aDetectorConstruction::DefineMaterials()
{
  // Material definition 

  G4NistManager* nistManager = G4NistManager::Instance();

  //ストロンチウム線源？
  //G4Material* strontiumMaterial = nistManager->FindOrBuildMaterial("G4_STRONTIUM");

  // Air defined using NIST Manager
  nistManager->FindOrBuildMaterial("G4_AIR");
  
  //コリメーターのマテリアルの定義
  fCollimatorMaterial  = nistManager->FindOrBuildMaterial("G4_Pb");
  
  // Lead defined using NIST Manager
  fTargetMaterial  = nistManager->FindOrBuildMaterial("G4_AIR");

  // Xenon gas defined using NIST Manager
  fChamberMaterial = nistManager->FindOrBuildMaterial("G4_POLYSTYRENE");

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B2aDetectorConstruction::DefineVolumes()
{
  G4Material* air  = G4Material::GetMaterial("G4_AIR");

  // Sizes of the principal geometrical components (solids)
  
  //G4double chamberSpacing = 80*cm; // from chamber center to center!

  G4double trackerWidth = 3.7*cm; // width of the chambers
  G4double targetWidth = 2.5*cm; // width of the target

  G4double targetLength = 0.5*cm; // full length of Target
  //G4double magnetLength = 4.0*cm;
  G4double trackerLength = 2*cm;

  G4double trackerRadius  = trackerWidth/2;   // Radius of Target
  G4double targetRadius  = targetWidth/2;   // Radius of Target


  G4double collimatorRmin = 0.1*cm; 
  G4double collimatorRmax = targetRadius;
  G4double collimatorDz = 0.2*cm;

  //G4double worldLength = 1.2 * (targetLength + trackerLength + magnetLength);
  G4double worldLength = 1.2 * (targetRadius*2 + trackerWidth + collimatorDz);
  G4double worldWidth = 1.2 * (targetRadius*3 + trackerLength);


  

  
  //targetLength = 0.5*targetLength;             // Half length of the Target  
  //G4double trackerSize   = 0.5*trackerLength;  // Half length of the Tracker

  // Definitions of Solids, Logical Volumes, Physical Volumes

  //ストロンチウム線源？
  
  // G4double strontiumRadius = 0.5*cm; // 適切な半径を設定
  // G4double strontiumZPosition = 1.0*cm; // 適切なZ座標を設定

  // G4Sphere* strontiumS = new G4Sphere("strontium", 0., strontiumRadius, 0., 360.*deg, 0., 180.*deg);
  // G4LogicalVolume* strontiumLV = new G4LogicalVolume(strontiumS, strontiumMaterial, "Strontium");
  // new G4PVPlacement(0, G4ThreeVector(0, 0, strontiumZPosition), strontiumLV, "Strontium", worldLV, false, 0, fCheckOverlaps);
  // G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  // G4ParticleDefinition* betaParticle = particleTable->FindParticle("e-");

  // G4DecayTable* strontiumDecayTable = new G4DecayTable();
  // G4Decay* strontiumDecay = new G4Decay();
  // strontiumDecay->SetDecayMode(G4RadioactiveDecay::BetaMinusDecay);
  // strontiumDecayTable->Insert(strontiumDecay);
  // strontiumLV->SetDecayTable(strontiumDecayTable);


  // World

  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(worldLength);

  G4cout << "Computed tolerance = "
         << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
         << " mm" << G4endl;

  G4Box* worldS
    = new G4Box("world",                                    //its name
                worldWidth/2,worldWidth/2,worldLength/2); //its size
  G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                 worldS,   //its solid
                 air,      //its material
                 "World"); //its name
  
  //  Must place the World Physical volume unrotated at (0,0,0).
  // 
  G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                 0,               // no rotation
                 G4ThreeVector(), // at (0,0,0)
                 worldLV,         // its logical volume
                 "World",         // its name
                 0,               // its mother  volume
                 false,           // no boolean operations
                 0,               // copy number
                 fCheckOverlaps); // checking overlaps 

  // コリメーター
  
  G4Tubs* collimatorS
   = new G4Tubs("collimator",collimatorRmin,collimatorRmax,collimatorDz, 0.*deg, 360.*deg);
  fLogicCollimator
    = new G4LogicalVolume(collimatorS, fCollimatorMaterial,"Collimator",0,0,0);
  G4ThreeVector positionCollimator = G4ThreeVector(targetRadius - (targetRadius*3 + trackerLength)/2,0,collimatorDz/2 - (targetRadius*2 + trackerWidth + collimatorDz)/2);

  new G4PVPlacement(0,               // no rotation
                  positionCollimator,  // at (x,y,z)
                  fLogicCollimator,    // its logical volume
                  "Collimator",        // its name
                  worldLV,         // its mother volume
                  false,           // no boolean operations
                  0,               // copy number
                  fCheckOverlaps); // checking overlaps 



  // Target
  
  G4ThreeVector positionTarget = G4ThreeVector(targetRadius - (targetRadius*3 + trackerLength)/2,0,targetRadius*2 + collimatorDz - (targetRadius*2 + trackerWidth + collimatorDz)/2);

  G4Tubs* targetS
    = new G4Tubs("target",0.,targetRadius,targetLength/2,0.*deg,360.*deg);
  fLogicTarget
    = new G4LogicalVolume(targetS, fTargetMaterial,"Target",0,0,0);
  
  G4RotationMatrix* rotationTarget = new G4RotationMatrix();
  rotationTarget->rotateX(90.*deg);

  new G4PVPlacement(rotationTarget,               // no rotation
                    positionTarget,  // at (x,y,z)
                    fLogicTarget,    // its logical volume
                    "Target",        // its name
                    worldLV,         // its mother volume
                    false,           // no boolean operations
                    0,               // copy number
                    fCheckOverlaps); // checking overlaps 

  G4cout << "Target is " << 2*targetLength/cm << " cm of "
         << fTargetMaterial->GetName() << G4endl;

  // Tracker 
 
  //G4ThreeVector positionTracker = G4ThreeVector(0,0,magnetLength + targetLength + trackerLength/2 -(magnetLength+targetLength+trackerLength)/2 );
  G4ThreeVector positionTracker = G4ThreeVector(targetRadius*3 + trackerLength/2 - (targetRadius*3 + trackerLength)/2,0,targetRadius*2 + trackerWidth/2 + collimatorDz - (targetRadius*2 + trackerWidth + collimatorDz)/2);

  G4Tubs* trackerS
    = new G4Tubs("tracker",0,trackerRadius,trackerLength/2, 0.*deg, 360.*deg);
  
  G4LogicalVolume* trackerLV
    = new G4LogicalVolume(trackerS, fChamberMaterial, "Tracker",0,0,0);  

  G4RotationMatrix* rotationTracker = new G4RotationMatrix();
  rotationTracker->rotateY(90.*deg);

  new G4PVPlacement(rotationTracker,               // no rotation
                    positionTracker, // at (x,y,z)
                    trackerLV,       // its logical volume
                    "Tracker",       // its name
                    worldLV,         // its mother  volume
                    false,           // no boolean operations
                    0,               // copy number
                    fCheckOverlaps); // checking overlaps 

  // Visualization attributes

  G4VisAttributes* boxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  G4VisAttributes* chamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));

  worldLV      ->SetVisAttributes(boxVisAtt);
  fLogicTarget ->SetVisAttributes(boxVisAtt);
  trackerLV    ->SetVisAttributes(boxVisAtt);

  // Tracker segments

  

  // Example of User Limits
  //
  // Below is an example of how to set tracking constraints in a given
  // logical volume
  //
  // Sets a max step length in the tracker region, with G4StepLimiter

  //G4double maxStep = 0.5*chamberWidth;
  //fStepLimit = new G4UserLimits(maxStep);
  //trackerLV->SetUserLimits(fStepLimit);
 
  /// Set additional contraints on the track, with G4UserSpecialCuts
  ///
  /// G4double maxLength = 2*trackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
  /// trackerLV->SetUserLimits(new G4UserLimits(maxStep,
  ///                                           maxLength,
  ///                                           maxTime,
  ///                                           minEkin));

  // Always return the physical world

  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 

void B2aDetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  G4String trackerChamberSDname = "B2/TrackerChamberSD";
  B2TrackerSD* aTrackerSD = new B2TrackerSD(trackerChamberSDname,
                                            "TrackerHitsCollection"); //HitsCollectionの場所
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  // Setting aTrackerSD to all logical volumes with the same name 
  // of "Chamber_LV".
  SetSensitiveDetector("Tracker", aTrackerSD, true);

  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.

  //一様磁場オブジェクトの作成
  G4MagneticField* magField = new	G4UniformMagField(G4ThreeVector(0, 0.3*tesla,0.));
  //新たにG4FieldManagerを作る　ー　コンストラクタにlocal fieldを表現するG4MagneticFieldオブジェクトを渡す
  G4FieldManager*	localFieldMgr	=	new	G4FieldManager(magField);	
  fLogicTarget->SetFieldManager(localFieldMgr, true);




  G4ThreeVector fieldValue = G4ThreeVector();
  //fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue); 前回のエラー
  //fMagFieldMessenger->SetVerboseLevel(1);
  
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

void B2aDetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
  fCheckOverlaps = checkOverlaps;
}  
