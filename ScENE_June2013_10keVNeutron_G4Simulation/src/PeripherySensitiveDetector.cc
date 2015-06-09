
#include "PeripherySensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"  
#include "G4StepPoint.hh"
#include "G4TouchableHandle.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
//#include "PeripheryHit.hh"


PeripherySensitiveDetector::PeripherySensitiveDetector(const G4String& name)
   :G4VSensitiveDetector(name)
{
   collectionName.insert( "PeripherySensitiveDetectorCollection" );
   PeripheryHitsCollectionID = -1;
}

PeripherySensitiveDetector::~PeripherySensitiveDetector() {}

void PeripherySensitiveDetector::Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent)
{
     // Create a new collection
   pHitsCollection = new PeripheryHitsCollection(SensitiveDetectorName, collectionName[0]);
 
  if ( PeripheryHitsCollectionID < 0 )
       PeripheryHitsCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID(pHitsCollection);
 
  // Add collection to the event
  hitsCollectionOfThisEvent->AddHitsCollection(PeripheryHitsCollectionID, pHitsCollection);
 
}

G4bool PeripherySensitiveDetector::ProcessHits(G4Step* pStep, G4TouchableHistory*)
{
    G4Track* pTrack = pStep->GetTrack();
    
    G4StepPoint*      pStepPoint = pStep->GetPostStepPoint();
    G4TouchableHandle      touch = pStepPoint->GetTouchableHandle();
     
    
    // Create Hit 
    PeripheryHit* pHit = new PeripheryHit();
    //PeripheryHit* pHit;
    
    
    pHit->SetParentID(pTrack->GetParentID());
    pHit->SetTrackID(pTrack->GetTrackID());
    pHit->SetParticleType(pTrack->GetDefinition()->GetParticleName());
    pHit->SetStepNumber(pStep->GetTrack()->GetCurrentStepNumber());
    
    
    if(touch)
    {
        G4VPhysicalVolume*    volume = touch->GetVolume();
        pHit->SetPhysicalVolumeName(volume->GetName());
    }
    else
    pHit->SetPhysicalVolumeName(G4String("Null"));
    
    
    /*    
     if(pTrack->GetParentID())
     pHit->SetParentType(m_hParticleTypes[pTrack->GetParentID()]);
     else
     pHit->SetParentType(G4String("none"));
     */ 
    
	if(pTrack->GetCreatorProcess())
		pHit->SetCreatorProcess(pTrack->GetCreatorProcess()->GetProcessName());
	else
		pHit->SetCreatorProcess(G4String("Null"));
    
	pHit->SetDepositingProcess(pStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName());
	pHit->SetPosition(pStep->GetPostStepPoint()->GetPosition());
	pHit->SetEnergyDeposited(pStep->GetTotalEnergyDeposit());
	pHit->SetKineticEnergy(pTrack->GetKineticEnergy());
	pHit->SetTime(pTrack->GetGlobalTime());
    
    pHitsCollection->insert( pHit ); 
	
    return true;    
}

void PeripherySensitiveDetector::EndOfEvent(G4HCofThisEvent*) {}
