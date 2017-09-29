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
// $Id: EEShashEventAction.cc 75604 2013-11-04 13:17:26Z gcosmo $
// 
/// \file EEShashEventAction.cc
/// \brief Implementation of the EEShashEventAction class

#include "common.h"
#include "EEShashEventAction.hh"
#include "EEShashCalorimeterSD.hh"
#include "EEShashCalorHit.hh"
#include "EEShashAnalysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

#include "CreateTree.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EEShashEventAction::EEShashEventAction()
 : G4UserEventAction(),
   fAbsHCID(-1),
   fActHCID(-1),
   fBgoHCID(-1),
   fFibrHCIDCore(-1),
   fFibrHCIDClad(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EEShashEventAction::~EEShashEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EEShashCalorHitsCollection* 
EEShashEventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  EEShashCalorHitsCollection* hitsCollection 
    = static_cast<EEShashCalorHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("EEShashEventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EEShashEventAction::PrintEventStatistics(
                              G4double absEdep, G4double absTrackLength,
                              G4double actEdep, G4double actTrackLength,
			      G4double bgoEdep, G4double bgoTrackLength,
			      G4double fibrEdepCore, G4double fibrTrackLengthCore,
			      G4double fibrEdepClad, G4double fibrTrackLengthClad) const
{
  // print event statistics
  G4cout
     << "   Absorber: total energy: " 
     << std::setw(7) << G4BestUnit(absEdep, "Energy")
     << "       total track length: " 
     << std::setw(7) << G4BestUnit(absTrackLength, "Length")
     << G4endl
     << "        act: total energy: " 
     << std::setw(7) << G4BestUnit(actEdep, "Energy")
     << "       total track length: " 
     << std::setw(7) << G4BestUnit(actTrackLength, "Length")
     << G4endl
     << "        bgo: total energy: " 
     << std::setw(7) << G4BestUnit(bgoEdep, "Energy")
     << "       total track length: " 
     << std::setw(7) << G4BestUnit(bgoTrackLength, "Length")
     << G4endl
     << "        fibr core: total energy: " 
     << std::setw(7) << G4BestUnit(fibrEdepCore, "Energy")
     << "       total track length: " 
     << std::setw(7) << G4BestUnit(fibrTrackLengthCore, "Length")
     << G4endl
     << "        fibr clad: total energy: " 
     << std::setw(7) << G4BestUnit(fibrEdepClad, "Energy")
     << "       total track length: " 
     << std::setw(7) << G4BestUnit(fibrTrackLengthClad, "Length")
     << G4endl;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EEShashEventAction::BeginOfEventAction(const G4Event* /*event*/)
{
  CreateTree::Instance() -> Clear();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EEShashEventAction::EndOfEventAction(const G4Event* event)
{  
  // Use the geometry variables from main
  extern int nLayers;
  extern int nBGOs;
  extern int nFibres;
  /*
  std::cout << "  Printing from EEShashEventAction.cc::EndOfEventAction() :" << G4endl;
  std::cout << "    Using nLayers = " << nLayers << G4endl;
  std::cout << "    Using nBGOs   = " << nBGOs << G4endl;
  std::cout << "    Using nFibres = " << nFibres << G4endl;
  */

  // Get hits collections IDs (only once)
  if ( fAbsHCID == -1 ) {
    fAbsHCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("AbsHitsCollection");
    fActHCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("ActHitsCollection");
    fAbs2HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("AbsHitsCollection2");
    fAct2HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("ActHitsCollection2");
    fAbs3HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("AbsHitsCollection3");
    fAct3HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("ActHitsCollection3");
//    fBgoHCID 
//      = G4SDManager::GetSDMpointer()->GetCollectionID("BgoHitsCollection");
    fFibrHCIDCore 
      = G4SDManager::GetSDMpointer()->GetCollectionID("FibrHitsCollectionCore");
    fFibrHCIDClad 
      = G4SDManager::GetSDMpointer()->GetCollectionID("FibrHitsCollectionClad");
    fAPDHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("APDHitsCollection");

    /*    fScint1HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("Scint1HitsCollection");
    fHodo11HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("Hodo11HitsCollection");
    fHodo12HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("Hodo12HitsCollection");
    */
  }

  // Get hits collections
  EEShashCalorHitsCollection* absHC = GetHitsCollection(fAbsHCID, event);
  EEShashCalorHitsCollection* actHC = GetHitsCollection(fActHCID, event);

  EEShashCalorHitsCollection* Abs2HC = GetHitsCollection(fAbs2HCID, event);
  EEShashCalorHitsCollection* Act2HC = GetHitsCollection(fAct2HCID, event);

  EEShashCalorHitsCollection* Abs3HC = GetHitsCollection(fAbs3HCID, event);
  EEShashCalorHitsCollection* Act3HC = GetHitsCollection(fAct3HCID, event);


  EEShashCalorHitsCollection* bgoHC = GetHitsCollection(fBgoHCID, event);
  EEShashCalorHitsCollection* fibrHCCore = GetHitsCollection(fFibrHCIDCore, event);
  EEShashCalorHitsCollection* fibrHCClad = GetHitsCollection(fFibrHCIDClad, event);
  EEShashCalorHitsCollection* APDHC = GetHitsCollection(fAPDHCID, event);

  // EEShashCalorHitsCollection* scint1HC = GetHitsCollection(fScint1HCID, event);
  // EEShashCalorHitsCollection* hodo11HC = GetHitsCollection(fHodo11HCID, event);
  //  EEShashCalorHitsCollection* hodo12HC = GetHitsCollection(fHodo12HCID, event);

  // Get hit with total values
  EEShashCalorHit* absHit = (*absHC)[absHC->entries()-1];
  EEShashCalorHit* actHit = (*actHC)[actHC->entries()-1];

  EEShashCalorHit* Abs2Hit = (*Abs2HC)[Abs2HC->entries()-1];
  EEShashCalorHit* Act2Hit = (*Act2HC)[Act2HC->entries()-1];

  EEShashCalorHit* Abs3Hit = (*Abs3HC)[Abs3HC->entries()-1];
  EEShashCalorHit* Act3Hit = (*Act3HC)[Act3HC->entries()-1];

  //  EEShashCalorHit* bgoHit = (*bgoHC)[bgoHC->entries()-1];
  EEShashCalorHit* fibrHitCore = (*fibrHCCore)[fibrHCCore->entries()-1];
  EEShashCalorHit* fibrHitClad = (*fibrHCClad)[fibrHCClad->entries()-1];
  // EEShashCalorHit* scint1Hit = (*scint1HC)[scint1HC->entries()-1];
  // EEShashCalorHit* hodo11Hit = (*hodo11HC)[hodo11HC->entries()-1];
  // EEShashCalorHit* hodo12Hit = (*hodo12HC)[hodo11HC->entries()-1];



  // Print per event (modulo n)
  //
  G4int eventID = event->GetEventID();
  G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "------------------------------------" << G4endl;     
    G4cout << "------------------------------------" << G4endl;     
    G4cout << "---> End of event: " << eventID << G4endl;     

    PrintEventStatistics(
      absHit->GetEdep(), absHit->GetTrackLength(),
      actHit->GetEdep(), actHit->GetTrackLength(),
      //      bgoHit->GetEdep(), bgoHit->GetTrackLength(),
      0,0,
      fibrHitCore->GetEdep(), fibrHitCore->GetTrackLength(),
      fibrHitClad->GetEdep(), fibrHitClad->GetTrackLength());
    G4cout << "------------------------------------" << G4endl;     

    G4cout << "  Verifying that PrintEventStatistics changes units when it wants" << G4endl;
    G4cout << "    fibrHitCore->GetEdep() = " << fibrHitCore->GetEdep() << G4endl;
    G4cout << "    fibrHitClad->GetEdep() = " << fibrHitClad->GetEdep() << G4endl;
 
  }  
  
  // Fill histograms, ntuple
  //
  std::cout << "Printing from EEShashEventAction::EndOfEventAction(): Filling branches" << std::endl;

  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 
  // fill histograms
  analysisManager->FillH1(1, absHit->GetEdep());
  analysisManager->FillH1(2, actHit->GetEdep());
  analysisManager->FillH1(3, absHit->GetTrackLength());
  analysisManager->FillH1(4, actHit->GetTrackLength());
  

  //std::cout << "xPosition = " << xBeamPos << std::endl;
  //std::cout << "yPosition = " << yBeamPos << std::endl;


  CreateTree::Instance() -> Event = event->GetEventID();
  CreateTree::Instance() -> nLayers = nLayers;


  CreateTree::Instance() -> Eabs_3x3 = absHit->GetEdep();
  CreateTree::Instance() -> Eact_3x3 = actHit->GetEdep();

  CreateTree::Instance() -> Eabs_CentralXtal=Abs2Hit->GetEdep();
  CreateTree::Instance() -> Eact_CentralXtal=Act2Hit->GetEdep();

  CreateTree::Instance() -> Eabs_1x3 = Abs3Hit->GetEdep();
  CreateTree::Instance() -> Eact_1x3 = Act3Hit->GetEdep();


  CreateTree::Instance() -> xPosition = xBeamPos;
  CreateTree::Instance() -> yPosition = yBeamPos;

  for(int i=0;i<nFibres;++i){
    EEShashCalorHit* APDHit = (*APDHC)[i]; 
    CreateTree::Instance() -> EAPD.push_back(APDHit->GetEdep());
  }

  CreateTree::Instance()->Fill(); 
  
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
