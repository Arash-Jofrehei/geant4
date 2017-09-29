#include "SteppingAction.hh"
#include "common.h"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"

#include "TMath.h"
#include "CreateTree.h"

using namespace std;
using namespace CLHEP;

//G4double fibre0;
//G4double fibre1;
//G4double fibre2;
//G4double fibre3;

int to_int (string name)
{
  int Result ; // int which will contain the result
  stringstream convert (name) ;
  string dummy ;
  convert >> dummy ;
  convert >> Result ;
  return Result ;
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

SteppingAction::SteppingAction (EEShashDetectorConstruction* detectorConstruction,
				const G4int& scint, const G4int& cher):
  fDetectorConstruction(detectorConstruction),
  propagateScintillation(scint),
  propagateCerenkov(cher)
{}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

SteppingAction::~SteppingAction ()
{}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

void SteppingAction::UserSteppingAction (const G4Step * theStep)
{
  G4Track* theTrack = theStep->GetTrack () ;
  G4int trackID = theTrack->GetTrackID();
  TrackInformation* theTrackInfo = (TrackInformation*)(theTrack->GetUserInformation());
  G4ParticleDefinition* particleType = theTrack->GetDefinition () ;

  //Pre point = first step in a volume
  //PostStep point identifies step out of volume
  G4StepPoint * thePrePoint = theStep->GetPreStepPoint () ;
  //  G4StepPoint * thePostPoint = theStep->GetPostStepPoint () ;
  const G4ThreeVector & thePrePosition = thePrePoint->GetPosition () ;
  G4VPhysicalVolume * thePrePV = thePrePoint->GetPhysicalVolume () ;
  //  G4VPhysicalVolume * thePostPV = thePostPoint->GetPhysicalVolume () ;
  G4String thePrePVName = "" ; 
  if ( thePrePV ) thePrePVName = thePrePV -> GetName () ;
  //  G4String thePostPVName = "" ; 
  //  if ( thePostPV ) thePostPVName = thePostPV -> GetName () ;

  
  G4int nStep = theTrack -> GetCurrentStepNumber();
  G4TouchableHandle theTouchable = thePrePoint->GetTouchableHandle();
  G4int copyNo = theTouchable->GetCopyNumber();    



  //-------------
  // get position
  G4double global_x = thePrePosition.x()/mm;
  G4double global_y = thePrePosition.y()/mm;
  G4double global_z = thePrePosition.z()/mm;

  G4TransportationManager* transportMgr ; 
  transportMgr = G4TransportationManager::GetTransportationManager() ;
  G4PropagatorInField * fieldPropagator = transportMgr->GetPropagatorInField() ;
  if(nStep>70000){
    //       std::cout<<"mortacci nstep"<<nStep<<" particle:"<<particleType->GetParticleName()<<" volume:"<<theTrack->GetLogicalVolumeAtVertex()->GetName()<<" id:"<<trackID<<" position:"<<global_x<<" "<<global_y<<" "<<global_z<<" energy:"<<theTrack->GetTotalEnergy()/eV<<std::endl;
    theTrack->SetTrackStatus(fStopAndKill);
  }
  
  //  if(theTrack->GetCreatorProcess() && particleType->GetParticleName()=="gamma")  std::cout<<"particle:"<<particleType->GetParticleName()<<" volume:"<<theTrack->GetLogicalVolumeAtVertex()->GetName()<<" id:"<<trackID<<" energy:"<<theTrack->GetTotalEnergy()/eV<<" process:"<<theTrack->GetCreatorProcess()->GetProcessName()<<std::endl;
  //  std::cout<<"particle: "<<particleType->GetParticleName()<<" step:"<<nStep<<" "<<direction.theta()<<" "<<direction.phi()<<std::endl;
  //  if(theTrack->GetCreatorProcess())  std::cout<<"particle:"<<particleType->GetParticleName()<<" "<<theTrack->GetCreatorProcess()->GetProcessName()<<" step:"<<nStep<<" "<<direction.theta()<<" "<<direction.phi()<<std::endl;


  //count particle in apd and get timing. done for just one fibre
  if(theTrack->GetLogicalVolumeAtVertex()->GetName().contains("APD") && (nStep==1) && copyNo == 0){
    CreateTree::Instance() -> Time_deposit_APD.push_back(theTrack->GetGlobalTime()/nanosecond);
    CreateTree::Instance() ->nParticlesAPD++;
  }



  return ;
}


