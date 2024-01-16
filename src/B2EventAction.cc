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
#include "B2TrackerHit.hh"

#include "G4SystemOfUnits.hh"
#include "g4root.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2EventAction::B2EventAction()
: G4UserEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2EventAction::~B2EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2EventAction::BeginOfEventAction(const G4Event*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2EventAction::EndOfEventAction(const G4Event* event)
{
  // get number of stored trajectories

  G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

  // periodic printing

  G4int eventID = event->GetEventID();
  if ( eventID < 100 || eventID % 100 == 0) {
    G4cout << ">>> Event: " << eventID  << G4endl;
    if ( trajectoryContainer ) {
      G4cout << "    " << n_trajectories
             << " trajectories stored in this event." << G4endl;
    }
    G4VHitsCollection* hc = event->GetHCofThisEvent()->GetHC(0);
    G4cout << "    "  
           << hc->GetSize() << " hits stored in this event" << G4endl;
  }
  // fill Ntuple
  auto analysisManager = G4AnalysisManager::Instance();
  G4HCofThisEvent* HCTE = static_cast<G4HCofThisEvent*>(event->GetHCofThisEvent());

  G4SDManager* sdManager = G4SDManager::GetSDMpointer();
  G4int CID = sdManager->GetCollectionID("TrackerHitsCollection");
  B2TrackerHitsCollection* HC = static_cast<B2TrackerHitsCollection*>(HCTE->GetHC(CID));

  G4int nhits = HC -> entries();
//欲しい情報をColumnにとってこよう！
  for(G4int i=0; i<nhits; i++){

    B2TrackerHit* hit = (*HC)[i];

    G4ThreeVector position = hit->GetPos();
    G4ThreeVector momentum = hit->GetMom(); //追加した
    //G4double momentum = hit->GetMom();

    analysisManager->FillNtupleDColumn(0, position.x()/cm);
    analysisManager->FillNtupleDColumn(1, position.y()/cm);
    analysisManager->FillNtupleDColumn(2, position.z()/cm);
    analysisManager->FillNtupleDColumn(3, momentum.x()*MeV); //追加した
    analysisManager->FillNtupleDColumn(4, momentum.y()*MeV);
    analysisManager->FillNtupleDColumn(5, momentum.z()*MeV);
    //analysisManager->FillNtupleDColumn(6, trackID*MeV);
    //analysisManager->FillNtupleDColumn(3, momentum()*MeV);
  

    analysisManager->AddNtupleRow();
    
  }

}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
