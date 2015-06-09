#include "TrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4VPhysicalVolume.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4TrackStatus.hh"



void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
    
    
    /*
     
    July 19th 2013
     
     The limits applied on Tracks of ScENE G4 simulation:
     
     1. The propagation of neutron and associated secondary particles in shielding:
     
     artifical cuts with threshold on energy and time
     
     electron energy threshold 20 keV
     protron energy threshold 20 keV
     
     GlobalTime Limit: 1000 ns
     
    */
    
    
    ///*
    G4double      electronKinEneCutThreshold = 40*keV;
    G4double        protonKinEneCutThreshold = 40*keV;
    G4double              GlobalTimeTrackCut = 1000;  //---Unit: ns---
    
    
    G4double     dKineticEnergy = aTrack->GetKineticEnergy();
    G4double        dGlobalTime = aTrack->GetGlobalTime();
    
    G4String           Particle = aTrack->GetDefinition()->GetParticleName();
    G4String         VolumeName = aTrack->GetVolume()->GetName();
    G4String            Process = G4String("Null");
    
    if(aTrack->GetCreatorProcess())
		G4String      Process = aTrack->GetCreatorProcess()->GetProcessName();

		

    //*/
    
    /*=======
     enum G4TrackStatus
     1. fAlive:                         continue tracking
     2. fStopButAlive                   Invoke active rest physics process and kill current track afterward
     3. fStopAndKill                    kill current track and sencondaries are alive
     4. fKillTrackAndSencondaries       Kill the current track and also associated secondaries
     5. fSuspend                        Suspend current track and send it back to stack. associated secondaries are also sent to stack
     6. fPostponeToNextEvent            Postpone tracking of curren track to next event, associated secondaries are sent to the stack.
     =================
     */
    
    
    ///*
    if(dGlobalTime>=GlobalTimeTrackCut){
    ((G4Track*)aTrack)->SetTrackStatus(fKillTrackAndSecondaries);   
    }
    
    if(Process==G4String("RadioactiveDecay"))   //----Ignore the radioactive decay contributions-----
    ((G4Track*)aTrack)->SetTrackStatus(fKillTrackAndSecondaries); 
    
    
    //-----Simplify the physics process in collimator-----
    if( ( (Particle==G4String("e-"))&&((VolumeName==G4String("EJ1_Shielding"))||(VolumeName==G4String("EJ2_Shielding"))||(VolumeName==G4String("EJ3_Shielding")))&&(dKineticEnergy<electronKinEneCutThreshold) ) )
    {
        //G4TrackStatus trackStatus = fKillTrackAndSecondaries;
        ((G4Track*)aTrack)->SetTrackStatus(fKillTrackAndSecondaries);
        
    }
    
    else if( (Particle==G4String("proton"))&&(dKineticEnergy<protonKinEneCutThreshold)&&((VolumeName==G4String("EJ1_Shielding"))||(VolumeName==G4String("EJ2_Shielding"))||(VolumeName==G4String("EJ3_Shielding"))) )
    {
        //((G4Track*)aTrack)->SetTrackStatus(fStopButAlive);
        ((G4Track*)aTrack)->SetTrackStatus(fKillTrackAndSecondaries);


    }
    //*/
    
    /*
    G4int    EventID = (G4EventManager::GetEventManager())->GetConstCurrentEvent()->GetEventID();
    G4int   ParentID = aTrack->GetParentID();
    G4int    TrackID = aTrack->GetTrackID();
    
    G4double     dKineticEnergy = aTrack->GetKineticEnergy();
    G4double        dGlobalTime = aTrack->GetGlobalTime();
    
    G4String           Particle = aTrack->GetDefinition()->GetParticleName();
    G4String         VolumeName = aTrack->GetVolume()->GetName();
    
    G4ThreeVector      Position = aTrack->GetVertexPosition();
    
    
    G4String     CreatorProcess;
    const G4VProcess* creatorProcess = aTrack->GetCreatorProcess();
    
    if(!creatorProcess) 
        CreatorProcess = G4String("null");   
    else
        CreatorProcess = aTrack->GetCreatorProcess()->GetProcessName();

    
    std::ofstream       outTest("TrackOutput.txt",std::ios::app);
    
    G4int width = 10, precision = 4;
     
    outTest<<std::setw(2)
    <<EventID<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
    <<ParentID<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
    <<TrackID<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
    <<Position.x()<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
    <<Position.y()<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
    <<Position.z()<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
    <<(dKineticEnergy/keV)<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
    <<dGlobalTime<<"    "<<std::setw(width)
    <<Particle<<"    "<<std::setw(width)
    <<CreatorProcess<<"    "<<std::setw(width)
    <<VolumeName<<"    "<<std::setw(width)
    <<G4endl;
    
    
    outTest.close();
     */

}

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack){}
