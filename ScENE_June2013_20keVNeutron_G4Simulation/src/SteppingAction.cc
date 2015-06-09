#include "SteppingAction.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VProcess.hh"
#include "G4TouchableHandle.hh"
#include "G4TouchableHistory.hh"
#include "G4EventManager.hh"


SteppingAction::SteppingAction()
{;}

SteppingAction::~SteppingAction()
{;}


void SteppingAction::UserSteppingAction(const G4Step * pStep)
{
        G4Track* pTrack = pStep->GetTrack();
    

    G4StepPoint       *thePrePoint = pStep->GetPreStepPoint();
    G4StepPoint      *thePostPoint = pStep->GetPostStepPoint();
    G4VPhysicalVolume    *thePrePV = thePrePoint->GetPhysicalVolume();
    G4String          thePrePVname = thePrePV->GetName();
    G4ThreeVector         Position = pStep->GetPostStepPoint()->GetPosition(); 
    G4ThreeVector   TrackDirection = pTrack->GetMomentum(); 
    G4ThreeVector        Direction = pStep->GetPreStepPoint()->GetMomentum();
    G4ThreeVector   CheckDirection = pStep->GetPreStepPoint()->GetMomentum();
    G4int width = 10, precision = 4;
    
    
    
    
        
    /*
    std::ofstream       outTest("StepOutput.txt",std::ios::app);
    


    outTest<<std::setw(2)
    <<(G4EventManager::GetEventManager())->GetConstCurrentEvent()->GetEventID()<<std::setw(width)<<std::fixed<<std::setprecision(precision)
    <<pTrack->GetParentID()<<std::setw(width)<<std::fixed<<std::setprecision(precision)
    <<pTrack->GetTrackID()<<std::setw(width)<<std::fixed<<std::setprecision(precision)
    <<pStep->GetTrack()->GetCurrentStepNumber()<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
    <<(Position.x()/mm)<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
    <<(Position.y()/mm)<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
    <<(Position.z()/mm)<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
    <<(pTrack->GetKineticEnergy()/keV)<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
    <<(pStep->GetTotalEnergyDeposit()/keV)<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
    <<pTrack->GetGlobalTime()<<"    "<<std::setw(width)
    <<pTrack->GetDefinition()->GetParticleName()<<"    "<<std::setw(width)
    <<pStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()<<"    "<<std::setw(width)
    //<<pTrack->GetCreatorProcess()->GetProcessName()<<"    "<<std::setw(width)
    <<thePrePVname<<"    "<<std::setw(width)
    <<(TrackDirection.x())<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
    <<(TrackDirection.y())<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
    <<(TrackDirection.z())<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
    <<G4endl;

    outTest.close();
   
    */
     
     

    /*
    // check if it is alive
    if(theTrack->GetTrackStatus()!=fAlive) { return; }
    
    // check if it is primary
    if(theTrack->GetParentID()!=0) { return; }
   
    // check if it is NOT muon
    G4ParticleDefinition * particleType = theTrack->GetDefinition();
    if((particleType==G4MuonPlus::MuonPlusDefinition())
       ||(particleType==G4MuonMinus::MuonMinusDefinition()))
    { return; }
    
    // check if it is entering to the calorimeter volume
    G4StepPoint * thePrePoint = theStep->GetPreStepPoint();
    G4VPhysicalVolume * thePrePV = thePrePoint->GetPhysicalVolume();
    G4String thePrePVname = thePrePV->GetName();
    if(thePrePVname(0,4)=="calo") { return; }
    G4StepPoint * thePostPoint = theStep->GetPostStepPoint();
    G4VPhysicalVolume * thePostPV = thePostPoint->GetPhysicalVolume();
    G4String thePostPVname = thePostPV->GetName();
    if(thePostPVname(0,4)!="calo") { return; }
    
    // then suspend the track
    theTrack->SetTrackStatus(fSuspend);
      */
}


