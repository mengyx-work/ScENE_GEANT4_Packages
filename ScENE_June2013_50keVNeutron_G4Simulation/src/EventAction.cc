                                                                                                      

#include "EventAction.hh"
#include "BeamTestHit.hh"
#include  "PeripheryHit.hh"

#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "TMath.h"

#include <iostream>


class AnalysisManager;

EventAction::EventAction(AnalysisManager* AnalysisManagerPointer)
{
    pAnalysisManager = AnalysisManagerPointer;
           //FileName = pAnalysisManager->GetDataFileName();
}



EventAction::~EventAction() {}


void EventAction::BeginOfEventAction(const G4Event*) 
{
    //------HitsCollectionID is used to define the Groups of HitsCollections, which's separated by different SensitiveDetector------
    
    //--- Collect the fHitCollectionID ---
          G4SDManager * SDman = G4SDManager::GetSDMpointer();
            fHitsCollectionID = SDman->GetCollectionID("MonitorCollection");
    
    //--- Collect the PeripheryHitsCollectionID ---
    G4SDManager * PeripherySDman = G4SDManager::GetSDMpointer();
       PeripheryHitsCollectionID = PeripherySDman->GetCollectionID("PeripherySensitiveDetectorCollection");
    
  

}



void EventAction::EndOfEventAction(const G4Event* event)
{
    
  G4HCofThisEvent *hitsCollectionOfThisEvent = event->GetHCofThisEvent();
    
    
 
    
    
    //======================================Section for PeripheryHitsCollection===========================================    
    
    if(PeripheryHitsCollectionID>=0)
    {
        // std::ofstream       outTest("PeripheryHitCollection.txt",std::ios::app);
        
        // std::ofstream outSortedTest("Sorted_PeripheryHitCollection.txt",std::ios::app);
        
        
        PeripheryHitsCollection* hitsCollection = dynamic_cast<PeripheryHitsCollection*>(hitsCollectionOfThisEvent->GetHC(PeripheryHitsCollectionID));
        
        G4int           EventID;
        G4int           ParentID;
        G4int           TrackID;
        G4String        DepositProcess;
        G4String        VolumeName;
        
        G4int numberOfHits = hitsCollection->GetSize();
        
        EventID = event->GetEventID();
        
    
               
               
        /*----Check the String Find Function------
        std::size_t found;
        
        std::string test(G4String("rough_one"));
        
        found = test.find("rough");
        
        if (found==std::string::npos)
        std::cout<<"the string test result is nothing!  "<<std::endl;
        else
            std::cout<<"found! "<<std::endl;
        */
        
        
        G4int            InteractingVolumeInt = 0;
        G4int        PeripheryObject_Physical = 1;
    
        G4bool    Scattered_PeripheryObject_Physical = false;

        
        for(int i=0; i<numberOfHits; i++)
        {  
            PeripheryHit* aHit = (*hitsCollection)[i];
            
                  ParentID = aHit->GetParentID();
                   TrackID = aHit->GetTrackID();
            DepositProcess = aHit->GetDepositingProcess();
                VolumeName = aHit->GetVolumeName();

            std::string     stringVolumeName(VolumeName);

            
            if( (Scattered_PeripheryObject_Physical==false)&&(DepositProcess!=G4String("Transportation"))&&(ParentID==0)&&(TrackID==1) )
            {
                Scattered_PeripheryObject_Physical = true;
                              InteractingVolumeInt = InteractingVolumeInt + PeripheryObject_Physical;
            }
          
            
            if(Scattered_PeripheryObject_Physical==true)
            {
                pAnalysisManager->InteractingVolumeInt = InteractingVolumeInt;
                stringVolumeName.clear();
                break;
            }
            
            if( (Scattered_PeripheryObject_Physical==false)&&(i==(numberOfHits-1)) )
            {
                pAnalysisManager->InteractingVolumeInt = InteractingVolumeInt;
                stringVolumeName.clear();
                break;
            }


        }//---end of the hits loop----
      
        
               
    }
    
    //====================================================================================================================
    
    
    
    
    
    
    
    
    //=================================Section for fHitCollection========================================
    if (fHitsCollectionID>=0)
    {
               
        
        //std::ofstream outFile(FileName.c_str(),std::ios::app);
        
        /*
        std::ofstream outFile("XenonChamber_GEANT4_RawOutput.txt",std::ios::app);

        
        if(!outFile) 
        {
            G4cout<<"//////////////Error opening outFile/////////////////"<<G4endl;
            return;
        }
        */
        
                
        
        bool EJ301_Triggered = false;  //--- to verify the coincidence ---
        bool   TPC_Triggered = false;  //--- to verify the coincidence ---

        
        G4int width = 10, precision = 4;
        
        G4int           EventID;
        G4int           ParentID;
        G4int           TrackID;	
        G4int           StepNumber;
        
        G4ThreeVector   Position;
        G4ThreeVector   Direction;
        G4double        dKinEnergy;
        G4double        dDepositeEnergy;
        G4String        Particle;
        
        G4String        PhysicalVolumeName;
        G4String        CreatorProcess;
        G4String        DepositProcess;
        G4String        PreStepProcess;
        G4String        ParentParticle;
        G4double        globalTime;

        
        
        BeamTestHitsCollection* hitsCollection = dynamic_cast<BeamTestHitsCollection*>(hitsCollectionOfThisEvent->GetHC(fHitsCollectionID));
        
        
        G4int   numberOfHits = hitsCollection->GetSize();
        
                     EventID = event->GetEventID();
        
        
        //--------Root Data structure-----------
        
        pAnalysisManager->EventID = EventID;
        
        pAnalysisManager->ParentID.clear();
        pAnalysisManager->TrackNbr.clear();
        pAnalysisManager->StepNbr.clear();
        pAnalysisManager->PositionX.clear();
        pAnalysisManager->PositionY.clear();
        pAnalysisManager->PositionZ.clear();
        pAnalysisManager->DirectionX.clear();
        pAnalysisManager->DirectionY.clear();
        pAnalysisManager->DirectionZ.clear();
        pAnalysisManager->dDepositeEnergy.clear();
        pAnalysisManager->dKinEnergy.clear();
        pAnalysisManager->globalTime.clear();
        pAnalysisManager->Particle.clear();
        pAnalysisManager->ParentParticle.clear();
        pAnalysisManager->CreatorProcess.clear();
        pAnalysisManager->DepositProcess.clear();
        pAnalysisManager->PhysicalVolume.clear();

        
        
        
        
        std::vector<int>             tempParentID;
        std::vector<int>             tempTrackNbr;
        std::vector<int>             tempStepNbr;

        std::vector<double>           tempPositionX;
        std::vector<double>           tempPositionY;
        std::vector<double>           tempPositionZ;
        std::vector<double>          tempDirectionX;
        std::vector<double>          tempDirectionY;
        std::vector<double>          tempDirectionZ;
        std::vector<double>      tempDepositeEnergy;
        std::vector<double>           tempKinEnergy;
        std::vector<double>          tempGlobalTime;
        
        std::vector<std::string>          tempParticle;
        std::vector<std::string>          tempCreatorProcess;
        std::vector<std::string>          tempDepositProcess;
        std::vector<std::string>          tempParentParticle;
        std::vector<std::string>          tempPhysicalVolume;
        
        
        std::vector<int>             SortedParentID;
        std::vector<int>             SortedTrackNbr;
        std::vector<int>             SortedStepNbr;
      

        std::vector<double>           SortedPositionX;
        std::vector<double>           SortedPositionY;
        std::vector<double>           SortedPositionZ;
        std::vector<double>          SortedDirectionX;
        std::vector<double>          SortedDirectionY;
        std::vector<double>          SortedDirectionZ;
        std::vector<double>      SortedDepositeEnergy;
        std::vector<double>          SortedKinEnergy;
        std::vector<double>          SortedGlobalTime;
        std::vector<std::string>                SortedParticle;
        std::vector<std::string>                SortedCreatorProcess;
        std::vector<std::string>                SortedDepositProcess;
        std::vector<std::string>                SortedParentParticle;
        std::vector<std::string>                SortedPhysicalVolume;

        
        
        
        tempParentID.clear();
        tempTrackNbr.clear();
        tempStepNbr.clear();
        tempPositionX.clear();
        tempPositionY.clear();
        tempPositionZ.clear();
        tempDirectionX.clear();
        tempDirectionY.clear();
        tempDirectionZ.clear();
        tempDepositeEnergy.clear();
        tempKinEnergy.clear();
        tempGlobalTime.clear();
        tempParticle.clear();
        tempCreatorProcess.clear();
        tempDepositProcess.clear();
        tempParentParticle.clear();
        tempPhysicalVolume.clear();  
        
        SortedParentID.clear();
        SortedTrackNbr.clear();
        SortedStepNbr.clear();
        SortedPositionX.clear();
        SortedPositionY.clear();
        SortedPositionZ.clear();
        SortedDirectionX.clear();
        SortedDirectionY.clear();
        SortedDirectionZ.clear();
        SortedDepositeEnergy.clear();
        SortedGlobalTime.clear();
        SortedParticle.clear();
        SortedCreatorProcess.clear();
        SortedDepositProcess.clear();
        SortedParentParticle.clear();
        SortedPhysicalVolume.clear();
        
    //-----------------------------------------------------    


        
           
        //------------- to check the coincidence with EJ301 neutron detector -------------------     
        for(int i=0; i<numberOfHits; i++)
        {
            BeamTestHit* aHit = (*hitsCollection)[i];
            
            Position = aHit->GetPosition()/mm;
            dDepositeEnergy = aHit->GetEnergyDeposited()/keV;
            
            //if( ( (Position.y()<-300)||(Position.y()>300) ) && (dDepositeEnergy>0.1) )
            
                if( ( TMath::Sqrt(Position.y()*Position.y()+Position.x()*Position.x())>200 ) && (dDepositeEnergy>0.1) )
                        EJ301_Triggered = true;   //--- hit in the coincident EJ301 ---

                if( ( TMath::Sqrt(Position.y()*Position.y()+Position.x()*Position.x())<50 ) && (dDepositeEnergy>0.1) )
                        TPC_Triggered = true;
            
                if((EJ301_Triggered==true)&&(TPC_Triggered==true))
                        break;

        }
        
        
        
        //------------ only if EJ301 neutron detector & TPC are both triggered by energy deposit------------------     
        //if(EJ301_Triggered)
        
        if((EJ301_Triggered==true)&&(TPC_Triggered==true))
        {
            
            for(int i=0; i<numberOfHits; i++)
            {  
                BeamTestHit* aHit = (*hitsCollection)[i];
                
                      ParentID = aHit->GetParentID();
                       TrackID = aHit->GetTrackID();
                    StepNumber = aHit->GetStepNumber();
                      Position = aHit->GetPosition()/mm;
                     Direction = aHit->GetDirection();
               dDepositeEnergy = aHit->GetEnergyDeposited()/keV;
                    dKinEnergy = aHit->GetKineticEnergy()/keV;
                      Particle = aHit->GetParticleType();
                CreatorProcess = aHit->GetCreatorProcess();
                DepositProcess = aHit->GetDepositingProcess();
                ParentParticle = aHit->GetParentType();
                    globalTime = aHit->GetTime();
            PhysicalVolumeName = aHit->GetVolumeName();
                
                
                std::string     stringVolumeName(PhysicalVolumeName);
                
                
                    if(dDepositeEnergy>0.01)    //== the 0.01keV is threshold of single step enegy despoit for step recording ====
                    {        
                            /*
                            outFile<<std::setw(2)
                            <<EventID<<std::setw(width)<<std::fixed<<std::setprecision(precision)
                            <<ParentID<<std::setw(width)<<std::fixed<<std::setprecision(precision)
                            <<TrackID<<std::setw(width)<<std::fixed<<std::setprecision(precision)
                            <<StepNumber<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
                            <<Position.x()<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
                            <<Position.y()<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
                            <<Position.z()<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
                            <<dKinEnergy<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
                            <<dDepositeEnergy<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
                            <<globalTime<<"    "<<std::setw(width)
                            <<ParentParticle<<"    "<<std::setw(width)	
                            <<Particle<<"    "<<std::setw(width)	
                            <<DepositProcess<<"    "<<std::setw(width)	
                            <<CreatorProcess<<"    "<<std::setw(width)	
                            <<Direction.x()<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
                            <<Direction.y()<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
                            <<Direction.z()<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
                            <<G4endl;
                             */
                    
                            tempParentID.push_back(ParentID);
                            tempTrackNbr.push_back(TrackID);
                            tempStepNbr.push_back(StepNumber);
                            tempPositionX.push_back(Position.x());
                            tempPositionY.push_back(Position.y());
                            tempPositionZ.push_back(Position.z());
                            tempDirectionX.push_back(Direction.x());
                            tempDirectionY.push_back(Direction.y());
                            tempDirectionZ.push_back(Direction.z());
                            tempDepositeEnergy.push_back(dDepositeEnergy);
                            tempKinEnergy.push_back(dKinEnergy);
                            tempGlobalTime.push_back(globalTime);
                            tempParticle.push_back(Particle);
                            tempCreatorProcess.push_back(CreatorProcess);
                            tempDepositProcess.push_back(DepositProcess);
                            tempParentParticle.push_back(ParentParticle);
                            tempPhysicalVolume.push_back(stringVolumeName);
                             
                    
                    }
                
               
                
            }
            
            
            //-------------------Sort the steps by global time--------------------------------
            
                        double                 TempTime;                /* Temporary value for swapping two array */
                           int                 TempElementOrder = -1;
                    const  int                 TotStepNbr = tempParentID.size();
            
                          bool                   SortFlag = true;
            
                        std::vector<int> ElementOrder;
            
                        ElementOrder.clear();
            
                                for(int q=0; q<TotStepNbr; q++)
                                ElementOrder.push_back(q);
     
            
            
            
                        for(int m=0; m<TotStepNbr-1; m++)
                        {
                                    SortFlag = true;
                
                                        for(int n=0; n<TotStepNbr-m-1; n++)
                                        {
                                                    if (tempGlobalTime.at(n)>tempGlobalTime.at(n+1))
                                                    {
                                                                        TempTime = tempGlobalTime.at(n);
                                                                TempElementOrder = ElementOrder.at(n);
                        
                                                            tempGlobalTime.at(n) = tempGlobalTime.at(n+1);
                                                              ElementOrder.at(n) = ElementOrder.at(n+1);
                        
                                                          tempGlobalTime.at(n+1) = TempTime;
                                                            ElementOrder.at(n+1) = TempElementOrder;
                        
                                                            SortFlag = false;
                                                    } // end if
                    
                                        } // end for n = ...
                
                            
                            if (SortFlag == true)
                                break;
                
                        } // end for m = ...
            
            //------------------------------------------------------------------------------------------  
            
                    for(int v=0; v<TotStepNbr; v++)
                    {
                            int element = ElementOrder.at(v);
                
                
                            SortedParentID.push_back(tempParentID.at(element));
                            SortedTrackNbr.push_back(tempTrackNbr.at(element));
                            SortedStepNbr.push_back(tempStepNbr.at(element));
                
                
                            SortedKinEnergy.push_back(tempKinEnergy.at(element));
                            SortedDepositeEnergy.push_back(tempDepositeEnergy.at(element));
                            SortedGlobalTime.push_back(tempGlobalTime.at(v));
                            SortedPositionX.push_back(tempPositionX.at(element));
                            SortedPositionY.push_back(tempPositionY.at(element));
                            SortedPositionZ.push_back(tempPositionZ.at(element));
                            SortedDirectionX.push_back(tempDirectionX.at(element));
                            SortedDirectionY.push_back(tempDirectionY.at(element));
                            SortedDirectionZ.push_back(tempDirectionZ.at(element));

                
                            SortedParticle.push_back(tempParticle.at(element));
                            SortedCreatorProcess.push_back(tempCreatorProcess.at(element));
                            SortedDepositProcess.push_back(tempDepositProcess.at(element));
                            SortedParentParticle.push_back(tempParentParticle.at(element));
                            SortedPhysicalVolume.push_back(tempPhysicalVolume.at(element));
                    }
            
            
            //------------------------------------------------------------------------------------------     
            
            pAnalysisManager->ParentID = SortedParentID;
            pAnalysisManager->TrackNbr = SortedTrackNbr;
            pAnalysisManager->StepNbr = SortedStepNbr;
            pAnalysisManager->PositionX = SortedPositionX;
            pAnalysisManager->PositionY = SortedPositionY;
            pAnalysisManager->PositionZ = SortedPositionZ;
            pAnalysisManager->DirectionX = SortedDirectionX;
            pAnalysisManager->DirectionY = SortedDirectionY;
            pAnalysisManager->DirectionZ = SortedDirectionZ;

            pAnalysisManager->dDepositeEnergy = SortedDepositeEnergy;
            pAnalysisManager->dKinEnergy = SortedKinEnergy;
            pAnalysisManager->globalTime = SortedGlobalTime;
            pAnalysisManager->Particle = SortedParticle;
            pAnalysisManager->CreatorProcess = SortedCreatorProcess;
            pAnalysisManager->DepositProcess = SortedDepositProcess;
            pAnalysisManager->ParentParticle = SortedParentParticle;
            pAnalysisManager->PhysicalVolume = SortedPhysicalVolume;
            
            pAnalysisManager->_DataTTree->Fill();

            
        }  //--- end of if(EJ301_Triggered) ---
        
        
        //outFile.close();    
        
    }  //--- end of if(fHitsCollectionID>=0) ---
//====================================================================================================================
    
    


    
    
    
} // --- end of the EndOfEvent function


