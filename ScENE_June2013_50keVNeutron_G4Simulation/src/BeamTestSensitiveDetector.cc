
#include "BeamTestSensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"  
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"

#include <map>
using std::map;

BeamTestSensitiveDetector::BeamTestSensitiveDetector(const G4String& name)
   :G4VSensitiveDetector(name)
{
   collectionName.insert( "MonitorCollection" );
   fHitsCollectionID = -1;
}

BeamTestSensitiveDetector::~BeamTestSensitiveDetector() {}



void BeamTestSensitiveDetector::Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent)
{
     // Create a new collection
   fHitsCollection = new BeamTestHitsCollection(SensitiveDetectorName, collectionName[0]);
 
  if ( fHitsCollectionID < 0 )
       fHitsCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
 
  // Add collection to the event
  hitsCollectionOfThisEvent->AddHitsCollection(fHitsCollectionID, fHitsCollection);
    
    m_ParentParticleType.clear();
 
}



G4bool BeamTestSensitiveDetector::ProcessHits(G4Step* pStep, G4TouchableHistory*)
{
    G4Track* pTrack = pStep->GetTrack();

    G4StepPoint*      pStepPoint = pStep->GetPostStepPoint();
    G4TouchableHandle      touch = pStepPoint->GetTouchableHandle();

    
    // Create Hit
    BeamTestHit* pHit = new BeamTestHit();
    
    pHit->SetParentID(pTrack->GetParentID());
    pHit->SetTrackID(pTrack->GetTrackID());
    pHit->SetParticleType(pTrack->GetDefinition()->GetParticleName());
    pHit->SetStepNumber(pStep->GetTrack()->GetCurrentStepNumber());
    
    
    //-------- Set Physical Volume Name ----------
    if(touch)
    {
        G4VPhysicalVolume*    volume = touch->GetVolume();
        pHit->SetPhysicalVolumeName(volume->GetName());
    }
    else
        pHit->SetPhysicalVolumeName(G4String("Null"));

    
    
    if( !m_ParentParticleType.count(pTrack->GetTrackID()) )
        m_ParentParticleType[pTrack->GetTrackID()] = pTrack->GetDefinition()->GetParticleName();
    
    
    //------- Set Particle Parent Type ----------
     if(pTrack->GetParentID())
     pHit->SetParentType(m_ParentParticleType[pTrack->GetParentID()]);
     else 
     pHit->SetParentType(G4String("none"));
     
    
	if(pTrack->GetCreatorProcess())
		pHit->SetCreatorProcess(pTrack->GetCreatorProcess()->GetProcessName());
	else
		pHit->SetCreatorProcess(G4String("Null"));
    
    
    
	pHit->SetDepositingProcess(pStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName());    
    pHit->SetDirection(pTrack->GetMomentumDirection());    
	pHit->SetPosition(pStep->GetPostStepPoint()->GetPosition());
	pHit->SetEnergyDeposited(pStep->GetTotalEnergyDeposit());
	pHit->SetKineticEnergy(pTrack->GetKineticEnergy());
	pHit->SetTime(pTrack->GetGlobalTime());
    
    fHitsCollection->insert( pHit ); 
	
    return true;    
}

void BeamTestSensitiveDetector::EndOfEvent(G4HCofThisEvent*) {}
