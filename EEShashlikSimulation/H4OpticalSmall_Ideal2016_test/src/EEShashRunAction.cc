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
// $Id: EEShashRunAction.cc 75215 2013-10-29 16:07:06Z gcosmo $
//
/// \file EEShashRunAction.cc
/// \brief Implementation of the EEShashRunAction class

// For reading environment variables
#include <iostream>
#include <cstdlib>

#include "EEShashRunAction.hh"
#include "EEShashAnalysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "common.h"

// This is too dangerous
// #define nLayers 15
// #define nBGOs 24
// #define nLayers 1
// #define nBGOs 24
// #define nFibers 4

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EEShashRunAction::EEShashRunAction( )
 : G4UserRunAction()
{ 

  // Get geometry definitions from main
  extern int nLayers;
  extern int nBGOs;
  extern int nFibres;
  std::cout << "Printing from EEShashRunAction.cc :" << G4endl;
  std::cout << "  Using nLayers = " << nLayers << G4endl;
  std::cout << "  Using nBGOs   = " << nBGOs << G4endl;
  std::cout << "  Using nFibres = " << nFibres << G4endl;

  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in EEShashAnalysis.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFirstHistoId(1);

  // Book histograms, ntuple
  //
  
  // Creating histograms
  analysisManager->CreateH1("edep_abs","Edep in absorber", 100, 0., 800*MeV);
  analysisManager->CreateH1("edep_act","Edep in active", 100, 0., 100*MeV);
  analysisManager->CreateH1("trkl_abs","trackL in absorber", 100, 0., 1*m);
  analysisManager->CreateH1("trkl_act","trackL in active", 100, 0., 50*cm);


  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EEShashRunAction::~EEShashRunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EEShashRunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  
  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();


  std::cout << "EEShashRunAction::BeginOfRunAction() finished. Event creation starting." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EEShashRunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->GetH1(1) ) {
    G4cout << "\n ----> print histograms statistic ";
    if(isMaster) {
      G4cout << "for the entire run \n" << G4endl; 
    }
    else {
      G4cout << "for the local thread \n" << G4endl; 
    }
    
    G4cout << " EAbs : mean = " 
       << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy") 
       << " rms = " 
       << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;
    
    G4cout << " Eact : mean = " 
       << G4BestUnit(analysisManager->GetH1(2)->mean(), "Energy") 
       << " rms = " 
       << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Energy") << G4endl;
    
    G4cout << " LAbs : mean = " 
      << G4BestUnit(analysisManager->GetH1(3)->mean(), "Length") 
      << " rms = " 
      << G4BestUnit(analysisManager->GetH1(3)->rms(),  "Length") << G4endl;

    G4cout << " Lact : mean = " 
      << G4BestUnit(analysisManager->GetH1(4)->mean(), "Length") 
      << " rms = " 
      << G4BestUnit(analysisManager->GetH1(4)->rms(),  "Length") << G4endl;
  }



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
