                                                                                                      

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
            
            /*
                VolumeName = aHit->GetVolumeName();

            std::string     stringVolumeName(VolumeName);
            */
            
            if( (Scattered_PeripheryObject_Physical==false)&&(DepositProcess!=G4String("Transportation"))&&(ParentID==0)&&(TrackID==1) )
            {
                Scattered_PeripheryObject_Physical = true;
                              InteractingVolumeInt = InteractingVolumeInt + PeripheryObject_Physical;
            }
          
            
            if(Scattered_PeripheryObject_Physical==true)
            {
                pAnalysisManager->InteractingVolumeInt = InteractingVolumeInt;
                //stringVolumeName.clear();
                break;
            }
            
            if( (Scattered_PeripheryObject_Physical==false)&&(i==(numberOfHits-1)) )
            {
                pAnalysisManager->InteractingVolumeInt = InteractingVolumeInt;
                //stringVolumeName.clear();
                break;
            }


        }//---end of the hits loop----
      
        
        /*
         G4int            InteractingVolumeInt = 0;
         G4int               AnodeCap_Physical = 1;
         G4int               PTFEWall_Physical = 2;
         G4int              FieldCage_Physical = 4;
         G4int             CathodeCap_Physical = 7;
         
         G4int           InnerSSDewar_Physical = 20;
         G4int           OuterSSDewar_Physical = 40;
         G4int          HeatExchanger_Physical = 100;
         G4int               Coldhead_Physical = 400;
         G4int             CopperRing_Physical = 1000;
         G4int                    PMT_Physical = 4000;   
         G4int             Collimator_Physical = 10000;
         
         G4int width = 10, precision = 4;

         
         G4bool    Scattered_Collimator_Physical = false;
         G4bool      Scattered_AnodeCap_Physical = false;
         G4bool      Scattered_PTFEWall_Physical = false;
         G4bool     Scattered_FieldCage_Physical = false;
         G4bool    Scattered_CathodeCap_Physical = false;
         G4bool    Scattered_CopperRing_Physical = false;
         
         G4bool       Scattered_InnerSSDewar_Physical = false;
         G4bool       Scattered_OuterSSDewar_Physical = false;
         G4bool           Scattered_Coldhead_Physical = false;
         G4bool      Scattered_HeatExchanger_Physical = false;
         G4bool                Scattered_PMT_Physical = false;
         */
        
        
        
        
        /*---reserved parameters------
         G4int           StepNumber;
         G4ThreeVector   Position;
         G4double        dKinEnergy;
         G4double        dDepositeEnergy;
         G4String        Particle;
         G4String        CreatorProcess;
         G4double        globalTime;
         */

        
        /*
        for(int i=0; i<numberOfHits; i++)
        {  
                    PeripheryHit* aHit = (*hitsCollection)[i];
            
                              ParentID = aHit->GetParentID();
                               TrackID = aHit->GetTrackID();
                        DepositProcess = aHit->GetDepositingProcess();
                            VolumeName = aHit->GetVolumeName();
        */    
            
            /*--------backup the output results-------------
             
                            StepNumber = aHit->GetStepNumber();
                              Position = aHit->GetPosition()/mm;
                       dDepositeEnergy = aHit->GetEnergyDeposited()/keV;
                            dKinEnergy = aHit->GetKineticEnergy()/keV;
                              Particle = aHit->GetParticleType();
                        CreatorProcess = aHit->GetCreatorProcess();
                            globalTime = aHit->GetTime();
            
            *///---------------------------------------------
          
        /*
             std::string     stringVolumeName(VolumeName);
                                     
            if( (Scattered_Collimator_Physical==false)&&(VolumeName==G4String("Collimator_tube_Physical"))&&(DepositProcess!=G4String("Transportation"))&&(ParentID==0)&&(TrackID==1) )
            {
                Scattered_Collimator_Physical = true;
                InteractingVolumeInt = InteractingVolumeInt + Collimator_Physical;
            }
        
            if( (Scattered_AnodeCap_Physical==false)&&( (VolumeName==G4String("AnodeCap_Physical")) || (VolumeName==G4String("UpperQuartzWindow_Physical")) )&&(DepositProcess!=G4String("Transportation"))&&(ParentID==0)&&(TrackID==1) )
            {
                Scattered_AnodeCap_Physical = true;
                       InteractingVolumeInt = InteractingVolumeInt + AnodeCap_Physical;
            }
            
            if( (Scattered_PTFEWall_Physical==false)&&(VolumeName==G4String("PTFEWall_Physical"))&&(DepositProcess!=G4String("Transportation"))&&(ParentID==0)&&(TrackID==1) )
            {
                Scattered_PTFEWall_Physical = true;
                       InteractingVolumeInt = InteractingVolumeInt + PTFEWall_Physical;
            }
            
            if( (Scattered_FieldCage_Physical==false)&&(VolumeName==G4String("FieldCage_Physical"))&&(DepositProcess!=G4String("Transportation"))&&(ParentID==0)&&(TrackID==1) )
            {
                Scattered_FieldCage_Physical = true;
                        InteractingVolumeInt = InteractingVolumeInt + FieldCage_Physical;
            }
            
            if( (Scattered_CathodeCap_Physical==false)&&( (VolumeName==G4String("CathodeCap_Physical")) || (VolumeName==G4String("LowerQuartzWindow_Physical")) )&&(DepositProcess!=G4String("Transportation"))&&(ParentID==0)&&(TrackID==1) )
            {
                Scattered_CathodeCap_Physical = true;
                         InteractingVolumeInt = InteractingVolumeInt + CathodeCap_Physical;
            }
            
            
            if( (Scattered_PMT_Physical==false)&&( (VolumeName==G4String("PMT_QuartzWindow_Physical"))||(VolumeName==G4String("PMT_Shell_Physical")) )&&(DepositProcess!=G4String("Transportation"))&&(ParentID==0)&&(TrackID==1) )
            {
                Scattered_PMT_Physical = true;
                  InteractingVolumeInt = InteractingVolumeInt + PMT_Physical;
            }
            
            
            if( (Scattered_CopperRing_Physical==false)&&(stringVolumeName.find("CopperRing_Physical")!=std::string::npos)&&(DepositProcess!=G4String("Transportation"))&&(ParentID==0)&&(TrackID==1) )
            {
                Scattered_CopperRing_Physical = true;
                         InteractingVolumeInt = InteractingVolumeInt + CopperRing_Physical;
            }
            
            
            if( (Scattered_InnerSSDewar_Physical==false)&&(stringVolumeName.find("SSInnerDewar")!=std::string::npos)&&(DepositProcess!=G4String("Transportation"))&&(ParentID==0)&&(TrackID==1) )
            {
                Scattered_InnerSSDewar_Physical = true;
                InteractingVolumeInt = InteractingVolumeInt + InnerSSDewar_Physical;
            }
             
             
            if( (Scattered_OuterSSDewar_Physical==false)&&(stringVolumeName.find("SSOuterDewar")!=std::string::npos)&&(DepositProcess!=G4String("Transportation"))&&(ParentID==0)&&(TrackID==1) )
            {
                Scattered_OuterSSDewar_Physical = true;
                           InteractingVolumeInt = InteractingVolumeInt + OuterSSDewar_Physical;
            }

            

            if( (Scattered_HeatExchanger_Physical==false)&&(stringVolumeName.find("HeatExchanger")!=std::string::npos)&&(DepositProcess!=G4String("Transportation"))&&(ParentID==0)&&(TrackID==1) )
            {
               Scattered_HeatExchanger_Physical = true;
                           InteractingVolumeInt = InteractingVolumeInt + HeatExchanger_Physical;
                
            }
            
            if( (Scattered_Coldhead_Physical==false)&&(stringVolumeName.find("Coldhead")!=std::string::npos)&&(DepositProcess!=G4String("Transportation"))&&(ParentID==0)&&(TrackID==1) )
            {
                    Scattered_Coldhead_Physical = true;
                           InteractingVolumeInt = InteractingVolumeInt + Coldhead_Physical;
                
            }

         
         
         pAnalysisManager->InteractingVolumeInt = InteractingVolumeInt;

            stringVolumeName.clear();
            
           
                        
        } 
        */
        
        

        
                
        
    }
    
    //====================================================================================================================
    
    
    
    
    
    
    
    
    //=================================Section for fHitCollection========================================
    if (fHitsCollectionID>=0)
    {
               
        
        //std::ofstream outFile(FileName.c_str(),std::ios::app);
        
        /*
        std::ofstream outFile("outFileResult.txt",std::ios::app);

        
        if(!outFile) 
        {
            G4cout<<"//////////////Error opening outFile/////////////////"<<G4endl;
            return;
        }
        */
        
                
        
        bool EJ301_Triggered = false;  //--- to verify the coincidence ---
        
        //G4int width = 10, precision = 4;
        
        G4int           EventID;
        G4int           ParentID;
        G4int           TrackID;	
        G4int           StepNumber;
        
        G4ThreeVector   Position;
        G4double        dDepositeEnergy;
        G4String        Particle;
        
        G4String        PhysicalVolumeName;
        G4String        DepositProcess;
        G4double        globalTime;

        /*
        G4ThreeVector   Direction;
        G4double        dKinEnergy;

        G4String        ParentParticle;
        G4String        CreatorProcess;
        G4String        PreStepProcess;
         */
        
        
        BeamTestHitsCollection* hitsCollection = dynamic_cast<BeamTestHitsCollection*>(hitsCollectionOfThisEvent->GetHC(fHitsCollectionID));
        
        //BeamTestHitsCollection* hitsCollection = (BeamTestHitsCollection*)(hitsCollectionOfThisEvent->GetHC(fHitsCollectionID));
        
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
        
        pAnalysisManager->dDepositeEnergy.clear();
        pAnalysisManager->globalTime.clear();
        pAnalysisManager->Particle.clear();
        pAnalysisManager->DepositProcess.clear();
        pAnalysisManager->PhysicalVolume.clear();

        /*
        pAnalysisManager->DirectionX.clear();
        pAnalysisManager->DirectionY.clear();
        pAnalysisManager->DirectionZ.clear();
        pAnalysisManager->dKinEnergy.clear();
        pAnalysisManager->ParentParticle.clear();
        pAnalysisManager->CreatorProcess.clear();
         */
        
        
        std::vector<int>             tempParentID;
        std::vector<int>             tempTrackNbr;
        std::vector<int>             tempStepNbr;

        std::vector<double>           tempPositionX;
        std::vector<double>           tempPositionY;
        std::vector<double>           tempPositionZ;
        
        std::vector<double>      tempDepositeEnergy;
        std::vector<double>          tempGlobalTime;
        
        std::vector<std::string>          tempParticle;
        std::vector<std::string>          tempDepositProcess;
        std::vector<std::string>          tempPhysicalVolume;
        
        /*
        std::vector<std::string>          tempCreatorProcess;
        std::vector<std::string>          tempParentParticle;

        std::vector<double>          tempDirectionX;
        std::vector<double>          tempDirectionY;
        std::vector<double>          tempDirectionZ;
        std::vector<double>           tempKinEnergy;
         */
        
        
        std::vector<int>             SortedParentID;
        std::vector<int>             SortedTrackNbr;
        std::vector<int>             SortedStepNbr;
      

        std::vector<double>           SortedPositionX;
        std::vector<double>           SortedPositionY;
        std::vector<double>           SortedPositionZ;
       
        std::vector<double>      SortedDepositeEnergy;
        std::vector<double>          SortedGlobalTime;
        
        std::vector<std::string>                SortedParticle;
        std::vector<std::string>                SortedDepositProcess;
        std::vector<std::string>                SortedPhysicalVolume;

        /*
        std::vector<std::string>                SortedCreatorProcess;
        std::vector<std::string>                SortedParentParticle;

        std::vector<double>          SortedDirectionX;
        std::vector<double>          SortedDirectionY;
        std::vector<double>          SortedDirectionZ;
        std::vector<double>          SortedKinEnergy;
         */
        
        
        tempParentID.clear();
        tempTrackNbr.clear();
        tempStepNbr.clear();
        tempPositionX.clear();
        tempPositionY.clear();
        tempPositionZ.clear();
      
        tempDepositeEnergy.clear();
        tempGlobalTime.clear();
        tempParticle.clear();
        tempDepositProcess.clear();
        tempPhysicalVolume.clear();  
        
        SortedParentID.clear();
        SortedTrackNbr.clear();
        SortedStepNbr.clear();
        SortedPositionX.clear();
        SortedPositionY.clear();
        SortedPositionZ.clear();
                
        SortedDepositeEnergy.clear();
        SortedGlobalTime.clear();
        SortedParticle.clear();
        SortedDepositProcess.clear();
        SortedPhysicalVolume.clear();
        
        /*
        tempDirectionX.clear();
        tempDirectionY.clear();
        tempDirectionZ.clear();
        
        tempKinEnergy.clear();
        tempCreatorProcess.clear();
        tempParentParticle.clear();

        
        SortedDirectionX.clear();
        SortedDirectionY.clear();
        SortedDirectionZ.clear();

        SortedParentParticle.clear();
         */
        
    //-----------------------------------------------------    


        
           
        //------------- to check the coincidence with EJ301 neutron detector -------------------     
        for(int i=0; i<numberOfHits; i++)
        {
            BeamTestHit* aHit = (*hitsCollection)[i];
            
            Position = aHit->GetPosition()/mm;
            dDepositeEnergy = aHit->GetEnergyDeposited()/keV;
            
            //if( ( (Position.y()<-300)||(Position.y()>300) ) && (dDepositeEnergy>0.1) )
            
                if( ( TMath::Sqrt(Position.y()*Position.y()+Position.x()*Position.x())>500 ) && (dDepositeEnergy>0.1) ){
                        EJ301_Triggered = true;   //--- hit in the coincident EJ301 ---
                        break;
                }
        }
        
        
        
        //------------ only if EJ301 neutron detector is triggered by energy deposit------------------     
        if(EJ301_Triggered)
        {
            
            for(int i=0; i<numberOfHits; i++)
            {  
                BeamTestHit* aHit = (*hitsCollection)[i];
                
                      ParentID = aHit->GetParentID();
                       TrackID = aHit->GetTrackID();
                    StepNumber = aHit->GetStepNumber();
                      Position = aHit->GetPosition()/mm;
               dDepositeEnergy = aHit->GetEnergyDeposited()/keV;
                      Particle = aHit->GetParticleType();
                DepositProcess = aHit->GetDepositingProcess();
                    globalTime = aHit->GetTime();
            PhysicalVolumeName = aHit->GetVolumeName();
               
                /*
                     Direction = aHit->GetDirection();
                    dKinEnergy = aHit->GetKineticEnergy()/keV;
                CreatorProcess = aHit->GetCreatorProcess();
                ParentParticle = aHit->GetParentType();
                */
                
                
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
                           
                            tempDepositeEnergy.push_back(dDepositeEnergy);
                            tempGlobalTime.push_back(globalTime);
                            tempParticle.push_back(Particle);
                            tempDepositProcess.push_back(DepositProcess);
                            tempPhysicalVolume.push_back(stringVolumeName);
                             
                        /*
                        tempKinEnergy.push_back(dKinEnergy);
                        tempDirectionX.push_back(Direction.x());
                        tempDirectionY.push_back(Direction.y());
                        tempDirectionZ.push_back(Direction.z());
                        tempCreatorProcess.push_back(CreatorProcess);
                        tempParentParticle.push_back(ParentParticle);
                        */
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
                
                
                            SortedDepositeEnergy.push_back(tempDepositeEnergy.at(element));
                            SortedGlobalTime.push_back(tempGlobalTime.at(v));
                            SortedPositionX.push_back(tempPositionX.at(element));
                            SortedPositionY.push_back(tempPositionY.at(element));
                            SortedPositionZ.push_back(tempPositionZ.at(element));
                           
                
                            SortedParticle.push_back(tempParticle.at(element));
                            SortedDepositProcess.push_back(tempDepositProcess.at(element));
                            SortedPhysicalVolume.push_back(tempPhysicalVolume.at(element));
                        
                        /*
                        SortedDirectionX.push_back(tempDirectionX.at(element));
                        SortedDirectionY.push_back(tempDirectionY.at(element));
                        SortedDirectionZ.push_back(tempDirectionZ.at(element));
                        SortedKinEnergy.push_back(tempKinEnergy.at(element));
                        SortedCreatorProcess.push_back(tempCreatorProcess.at(element));
                        SortedParentParticle.push_back(tempParentParticle.at(element));
                        */


                        
                    }
            
            
            
            
            //------------------------------------------------------------------------------------------     
            
            pAnalysisManager->ParentID = SortedParentID;
            pAnalysisManager->TrackNbr = SortedTrackNbr;
            pAnalysisManager->StepNbr = SortedStepNbr;
            pAnalysisManager->PositionX = SortedPositionX;
            pAnalysisManager->PositionY = SortedPositionY;
            pAnalysisManager->PositionZ = SortedPositionZ;
            
      
            pAnalysisManager->dDepositeEnergy = SortedDepositeEnergy;
            pAnalysisManager->globalTime = SortedGlobalTime;
            pAnalysisManager->Particle = SortedParticle;
            pAnalysisManager->DepositProcess = SortedDepositProcess;
            pAnalysisManager->PhysicalVolume = SortedPhysicalVolume;
            
            /*
            pAnalysisManager->DirectionX = SortedDirectionX;
            pAnalysisManager->DirectionY = SortedDirectionY;
            pAnalysisManager->DirectionZ = SortedDirectionZ;
            pAnalysisManager->ParentParticle = SortedParentParticle;
            pAnalysisManager->dKinEnergy = SortedKinEnergy;
            pAnalysisManager->CreatorProcess = SortedCreatorProcess;
             */
            
            
            pAnalysisManager->_DataTTree->Fill();

            
        }  //--- end of if(EJ301_Triggered) ---
        
        
        //outFile.close();    
        
    }  //--- end of if(fHitsCollectionID>=0) ---
//====================================================================================================================
    
    


    
    
    
} // --- end of the EndOfEvent function


