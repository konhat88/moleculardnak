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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publications:
// Med. Phys. 45 (2018) e722-e739
// Phys. Med. 31 (2015) 861-874
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file TrackingAction.cc
/// \brief Implementation of the TrackingAction class

#include "TrackingAction.hh"
#include "AnalysisManager.hh"

#include "G4Alpha.hh"
#include "G4AnalysisManager.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4PhysicalConstants.hh"
#include "G4Proton.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(AnalysisManager* man)
  : fAnalysisManager(man) 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  G4double flagParticle = -1.;
  G4double x, y, z, dirx, diry, dirz;

  // Particle identification

  // The following method avoids the usage of string comparison
  G4ParticleDefinition* partDef = aTrack->GetDynamicParticle()->GetDefinition();
  G4DNAGenericIonsManager* instance = G4DNAGenericIonsManager::Instance();

  if (partDef == G4Gamma::GammaDefinition())
     flagParticle = 0; 
  else if (partDef == G4Electron::ElectronDefinition())
     flagParticle = 1; 
  else if (partDef == G4Proton::ProtonDefinition())
     flagParticle = 2; 
  else if (partDef == G4Alpha::AlphaDefinition())
     flagParticle = 4; 

  else if (partDef && partDef->GetAtomicNumber() == 89 && partDef->GetAtomicMass() == 225){
     flagParticle = 20; //Ac225
     G4cout << "ParticleName: " << partDef->GetParticleName() << G4endl;}
  else if (partDef && partDef->GetAtomicNumber() == 87 && partDef->GetAtomicMass() == 221){
     flagParticle = 21; //Fr221
     G4cout << "ParticleName: " <<  partDef->GetParticleName() << G4endl;}
  else if (partDef && partDef->GetAtomicNumber() == 85 && partDef->GetAtomicMass() == 217){
     flagParticle = 22; //At217
     G4cout << "ParticleName: " <<  partDef->GetParticleName() << G4endl;}
  else if (partDef && partDef->GetAtomicNumber() == 83 && partDef->GetAtomicMass() == 213){
     flagParticle = 23; //Bi213
     G4cout << "ParticleName: " <<  partDef->GetParticleName() << G4endl;}
  else if (partDef && partDef->GetAtomicNumber() == 81 && partDef->GetAtomicMass() == 209){
     flagParticle = 24; //Tl209
     G4cout << "ParticleName: " <<  partDef->GetParticleName() << G4endl;}
  else if (partDef && partDef->GetAtomicNumber() == 82 && partDef->GetAtomicMass() == 209){
     flagParticle = 25; //Pb209
     G4cout << "ParticleName: " <<  partDef->GetParticleName() << G4endl;}
  else if (partDef && partDef->GetAtomicNumber() == 84 && partDef->GetAtomicMass() == 213){
     flagParticle = 26; //Po213
     G4cout << "ParticleName: " <<  partDef->GetParticleName() << G4endl;}
  else if (partDef && partDef->GetAtomicNumber() == 83 && partDef->GetAtomicMass() == 209){
     flagParticle = 27; //Bi209
     G4cout << "ParticleName: " <<  partDef->GetParticleName() << G4endl;}


  else if (partDef && partDef->GetAtomicNumber() == 71 && partDef->GetAtomicMass() == 177) 
     flagParticle = 30; //Lu177
  else if (partDef && partDef->GetAtomicNumber() == 72 && partDef->GetAtomicMass() == 177) 
     flagParticle = 31; //Hf177

  else if (partDef && partDef->GetAtomicNumber() == 82 && partDef->GetAtomicMass() == 212 
                   && partDef->GetParticleName() == "Pb212"){
     flagParticle = 40; //Pb212
     G4cout << "ParticleName: " <<   partDef->GetParticleName() << G4endl;}
  else if (partDef && partDef->GetAtomicNumber() == 83 && partDef->GetAtomicMass() == 212
                   && partDef->GetParticleName() == "Bi212"){
     flagParticle = 41; //Bi212
     G4cout << "ParticleName: " <<   partDef->GetParticleName() << G4endl;}
  else if (partDef && partDef->GetAtomicNumber() == 84 && partDef->GetAtomicMass() == 212
                   && partDef->GetParticleName() == "Po212"){
     flagParticle = 42; //Po212
     G4cout << "ParticleName: " <<   partDef->GetParticleName() << G4endl;}
  else if (partDef && partDef->GetAtomicNumber() == 81 && partDef->GetAtomicMass() == 208
                   && partDef->GetParticleName() == "Tl208"){
     flagParticle = 43; //Tl208
     G4cout << "ParticleName: " <<   partDef->GetParticleName() << G4endl;}
  else if (partDef && partDef->GetAtomicNumber() == 82 && partDef->GetAtomicMass() == 208
                   && partDef->GetParticleName() == "Pb208"){
     flagParticle = 44; //Pb208
     G4cout << "ParticleName: " <<   partDef->GetParticleName() << G4endl;}
  
  

  else if (partDef == instance->GetIon("hydrogen"))
     flagParticle = 3; 

  else if (partDef == instance->GetIon("alpha+"))
     flagParticle = 5; 

  else if (partDef == instance->GetIon("helium"))
     flagParticle = 6; 
  
  else if (partDef->GetParticleName() == "anti_nu_e") flagParticle = 7;
  else G4cout <<"Could not find exact match of: " << partDef->GetParticleName() << G4endl;

  //
/*
  x = aTrack->GetPosition().x() / nanometer;
  y = aTrack->GetPosition().y() / nanometer;
  z = aTrack->GetPosition().z() / nanometer;

  dirx = aTrack->GetMomentumDirection().x();
  diry = aTrack->GetMomentumDirection().y();
  dirz = aTrack->GetMomentumDirection().z();

  // Call analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Fill track information ntuple
  analysisManager->FillNtupleDColumn(1, 0, flagParticle);
  analysisManager->FillNtupleDColumn(1, 1, x);
  analysisManager->FillNtupleDColumn(1, 2, y);
  analysisManager->FillNtupleDColumn(1, 3, z);
  analysisManager->FillNtupleDColumn(1, 4, dirx);
  analysisManager->FillNtupleDColumn(1, 5, diry);
  analysisManager->FillNtupleDColumn(1, 6, dirz);
  analysisManager->FillNtupleDColumn(1, 7, aTrack->GetKineticEnergy() / eV);
  analysisManager->FillNtupleIColumn(1, 8, aTrack->GetTrackID());
  analysisManager->FillNtupleIColumn(1, 9, aTrack->GetParentID());
  analysisManager->AddNtupleRow(1);
*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track*) {}
