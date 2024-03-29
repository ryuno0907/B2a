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
/// \file B2PrimaryGeneratorAction.cc
/// \brief Implementation of the B2PrimaryGeneratorAction class

#include "B2PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"
#include "G4Geantino.hh"
#include "Randomize.hh"
#include "parameter.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2PrimaryGeneratorAction::B2PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction()
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

   
  //G4ParticleDefinition* particleDefinition 
  // = G4ParticleTable::GetParticleTable()->FindParticle("e-");

  //fParticleGun->SetParticleDefinition(particleDefinition);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
  fParticleGun->SetParticleEnergy(0*MeV);
  //fParticleGun -> SetParticlePosition(G4ThreeVector(p_targetRadius*cm - (p_targetRadius*3.0 + p_trackerLength)*cm/2.0,0,- (p_targetRadius*2.0 + p_trackerWidth + p_collimatorDz)*cm/2.0));
  

  //β線源の配置

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2PrimaryGeneratorAction::~B2PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume
  // from G4LogicalVolumeStore.

  G4double worldZHalfLength = 0;
  G4LogicalVolume* worldLV
    = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
  
  G4Box* worldBox = NULL;
  if ( worldLV ) worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  if ( worldBox ) worldZHalfLength = worldBox->GetZHalfLength();
  else  {
    G4cerr << "World volume of box not found." << G4endl;
    G4cerr << "Perhaps you have changed geometry." << G4endl;
    G4cerr << "The gun will be place in the center." << G4endl;
  }

  //fParticleGun -> SetParticlePosition(G4ThreeVector(0., 0., 0.));
  
  
  if (fParticleGun->GetParticleDefinition()==G4Geantino::Geantino()){
    G4int Z = 38;
    G4int A = 90;
    G4double ionCharge = 0*eplus;
    G4double excitEnergy = 0*keV;

    

    G4ParticleDefinition* ionDefinition 
      = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
    fParticleGun->SetParticleDefinition(ionDefinition);
    fParticleGun->SetParticleCharge(ionCharge);
  }
  
  

  // Note that this particular case of starting a primary particle on the world boundary
  // requires shooting in a direction towards inside the world.
  fParticleGun->SetParticlePosition(G4ThreeVector(p_targetRadius*cm - (p_targetRadius*3 + p_trackerLength)*cm/2,0,0 - (p_targetRadius*2 + p_trackerWidth + p_collimatorDz)*cm/2));
  

  


  fParticleGun->GeneratePrimaryVertex(anEvent);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
