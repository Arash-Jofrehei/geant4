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
/// \file electromagnetic/TestEm2/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// $Id: SteppingAction.cc 76259 2013-11-08 11:37:28Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "EEShashDetectorConstruction.hh"
#include "EEShashRunAction.hh"

#include "G4RunManager.hh"
#include "G4SteppingManager.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EEShashDetectorConstruction* det, EEShashRunAction* runAct)
  :G4UserSteppingAction(),fDetector(det), fRunAct(runAct)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{ 
  
  // energy deposit
  //
  G4double dEStep = step->GetTotalEnergyDeposit();

  
  if (dEStep > 0.) {
    G4ThreeVector prePoint  = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector delta = step->GetPostStepPoint()->GetPosition() - prePoint;
    prePoint += G4UniformRand()*delta;
    G4double x = prePoint.x(), y = prePoint.y(), z = prePoint.z();
    G4double radius = std::sqrt(x*x + y*y);
    //   G4double offset = 0.5*fDetector->GetfullLength();
    G4double offset = 0.;
     G4int SlideNb = 1;
    //   G4int SlideNb = G4int((z + offset)/fDetector->GetdLlength());
    G4int RingNb  = 1;        
    //    G4int RingNb  = G4int(radius/fDetector->GetdRlength());        
   
    fRunAct->FillPerStep(dEStep,SlideNb,RingNb);
    }
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
