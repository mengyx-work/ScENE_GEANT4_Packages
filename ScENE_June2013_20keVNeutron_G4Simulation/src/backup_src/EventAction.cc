                                                                                                      

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
           FileName = pAnalysisManager->GetDataFileName();
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
        
        
        G4int width = 10, precision = 4;
        
        G4int            InteractingVolumeInt = 0;
        G4int               SSFlange_Physical = 20;
        G4int                SSDewar_Physical = 40;
        G4int               AnodeCap_Physical = 1;
        G4int               PTFEWall_Physical = 2;
        G4int              FieldCage_Physical = 4;
        G4int             CathodeCap_Physical = 7;
        G4int             CopperRing_Physical = 100;
        G4int                    PMT_Physical = 400;
        G4int           QuartzWindow_Physical = 200;
        
        
        G4bool      Scattered_SSFlange_Physical = false;
        G4bool       Scattered_SSDewar_Physical = false;
        G4bool      Scattered_AnodeCap_Physical = false;
        G4bool      Scattered_PTFEWall_Physical = false;
        G4bool     Scattered_FieldCage_Physical = false;
        G4bool    Scattered_CathodeCap_Physical = false;
        G4bool    Scattered_CopperRing_Physical = false;
        G4bool  Scattered_QuartzWindow_Physical = false;
        G4bool           Scattered_PMT_Physical = false;
        
        
        
        G4int           EventID;
        G4int           ParentID;
        G4int           TrackID;	
        G4int           StepNumber;
        G4ThreeVector   Position;
        G4double        dKinEnergy;
        G4double        dDepositeEnergy;
        G4String        Particle;
        G4String        CreatorProcess;
        G4String        DepositProcess;
        G4String        VolumeName;
        G4double        globalTime;
        
       
        
        
        G4int numberOfHits = hitsCollection->GetSize();
        
        EventID = event->GetEventID();
        
        
        const double ScatteringAngle = 90;           //--- Angle in degree --
        const double   NeutronEnergy = 500;          //--- Unit: keV ---
        const double     NeutronMass = 939565.38;    //--- Unit: keV ---
        const double  EJ301_Distance = 0.5;          //--- Unit: m ---
        const double    LAr_Distance = 0.5;          //--- Unit: m ---
        const double    Ar_AtomicNbr = 40;
        const double           Theta = ScatteringAngle/180*TMath::Pi();
        
        const double  RecoilEnergy = 2*NeutronEnergy/(1+Ar_AtomicNbr)/(1+Ar_AtomicNbr)*(TMath::Sin(Theta)*TMath::Sin(Theta) + Ar_AtomicNbr - TMath::Cos(Theta)*TMath::Sqrt(Ar_AtomicNbr*Ar_AtomicNbr - TMath::Sin(Theta)*TMath::Sin(Theta)));
        
        
        const double  ScatteredNeutronSpeed = 3*TMath::Sqrt(2*(NeutronEnergy-RecoilEnergy)/NeutronMass)/10;    //--- Unit: m/ns ---
        const double           NeutronSpeed = 3*TMath::Sqrt(2*NeutronEnergy/NeutronMass)/10;                   //--- Unit: m/ns ---
        
        const double               Scattering_TOF = LAr_Distance/NeutronSpeed;
        const double             TOF_Window_Mid   = EJ301_Distance/ScatteredNeutronSpeed;
        
        
        /*
        std::size_t found;
        
        std::string test(G4String("rough_one"));
        
        found = test.find("rough");
        
        if (found==std::string::npos)
        std::cout<<"the string test result is nothing!  "<<std::endl;
        else
            std::cout<<"found! "<<std::endl;
        */
        
        
        for(int i=0; i<numberOfHits; i++)
        {  
                    PeripheryHit* aHit = (*hitsCollection)[i];
            
                              ParentID = aHit->GetParentID();
                               TrackID = aHit->GetTrackID();
                        DepositProcess = aHit->GetDepositingProcess();
                            VolumeName = aHit->GetVolumeName();
            /*
                            StepNumber = aHit->GetStepNumber();
                              Position = aHit->GetPosition()/mm;
                       dDepositeEnergy = aHit->GetEnergyDeposited()/keV;
                            dKinEnergy = aHit->GetKineticEnergy()/keV;
                              Particle = aHit->GetParticleType();
                        CreatorProcess = aHit->GetCreatorProcess();
                            globalTime = aHit->GetTime();
            
            */
            
             std::string     stringVolumeName(VolumeName);
            
            
            if( (Scattered_SSFlange_Physical==false)&&(VolumeName==G4String("SSFlange_Physical"))&&(DepositProcess!=G4String("Transportation"))&&(ParentID==0)&&(TrackID==1) )
            {
                Scattered_SSFlange_Physical = true;
                       InteractingVolumeInt = InteractingVolumeInt + SSFlange_Physical;
            }
            
            if( (Scattered_SSDewar_Physical==false)&&(VolumeName==G4String("SSDewar_Physical"))&&(DepositProcess!=G4String("Transportation"))&&(ParentID==0)&&(TrackID==1) )
            {
                Scattered_SSDewar_Physical = true;
                      InteractingVolumeInt = InteractingVolumeInt + SSDewar_Physical;
            }
            
            if( (Scattered_AnodeCap_Physical==false)&&(VolumeName==G4String("AnodeCap_Physical"))&&(DepositProcess!=G4String("Transportation"))&&(ParentID==0)&&(TrackID==1) )
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
            
            if( (Scattered_CathodeCap_Physical==false)&&(VolumeName==G4String("CathodeCap_Physical"))&&(DepositProcess!=G4String("Transportation"))&&(ParentID==0)&&(TrackID==1) )
            {
                Scattered_CathodeCap_Physical = true;
                         InteractingVolumeInt = InteractingVolumeInt + CathodeCap_Physical;
            }
            
            
            if( (Scattered_QuartzWindow_Physical==false)&&( (VolumeName==G4String("UpperQuartzWindow_Physical"))||(VolumeName==G4String("LowerQuartzWindow_Physical")) )&&(DepositProcess!=G4String("Transportation"))&&(ParentID==0)&&(TrackID==1) )
            {
                Scattered_QuartzWindow_Physical = true;
                           InteractingVolumeInt = InteractingVolumeInt + QuartzWindow_Physical;
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

            
            stringVolumeName.clear();
            
            
            /*
             outTest<<std::setw(2)
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
             <<Particle<<"    "<<std::setw(width)	
             <<CreatorProcess<<"    "<<std::setw(width)
             <<DepositProcess<<"    "<<std::setw(width)
             <<VolumeName<<"    "<<std::setw(width)
             <<InteractingVolumeInt
             <<G4endl;
             */
            
        } 
        
        
        //outTest.close();
        
        
        pAnalysisManager->InteractingVolumeInt = InteractingVolumeInt;
        
        
        /*
         for(int p=0; p<TotStepNbr; p++)
         {
         outSortedTest<<std::setw(2)
         <<EventID<<std::setw(width)<<std::fixed<<std::setprecision(precision)
         <<SortedParentID[p]<<std::setw(width)<<std::fixed<<std::setprecision(precision)
         <<SortedTrackNbr[p]<<std::setw(width)<<std::fixed<<std::setprecision(precision)
         <<SortedStepNbr[p]<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
         <<SortedPositionX[p]<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
         <<SortedPositionY[p]<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
         <<SortedPositionZ[p]<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
         <<SortedKinEnergy[p]<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
         <<SortedDepositeEnergy[p]<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
         <<SortedGlobalTime[p]<<"    "<<std::setw(width)
         <<SortedParticle[p]<<"    "<<std::setw(width)	
         <<SortedCreatorProcess[p]<<"    "<<std::setw(width)
         <<SortedDepositProcess[p]<<"    "<<std::setw(width)
         <<G4endl;
         
         }
         
         
         
         for(int p=0; p<TotStepNbr; p++)
         {
         outSortedTest<<std::setw(2)
         <<EventID<<std::setw(width)<<std::fixed<<std::setprecision(precision)
         <<tempParentID[p]<<std::setw(width)<<std::fixed<<std::setprecision(precision)
         <<tempTrackNbr[p]<<std::setw(width)<<std::fixed<<std::setprecision(precision)
         <<tempStepNbr[p]<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
         <<tempPositionX[p]<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
         <<tempPositionY[p]<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
         <<tempPositionZ[p]<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
         <<tempKinEnergy[p]<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
         <<tempDepositeEnergy[p]<<"    "<<std::setw(width)<<std::fixed<<std::setprecision(precision)
         <<tempGlobalTime[p]<<"    "<<std::setw(width)
         <<tempParticle[p]<<"    "<<std::setw(width)	
         <<tempCreatorProcess[p]<<"    "<<std::setw(width)
         <<tempDepositProcess[p]<<"    "<<std::setw(width)
         <<G4endl;
         
         }
         */
        
        
        
    }
    
    //====================================================================================================================
    
    
    
    
    
    
    //=================================Section for fHitCollection========================================
    if (fHitsCollectionID>=0)
    {
        //std::ofstream outFile("Random_Result.dat",std::ios::app);
        
        /*
        std::ofstream outFile(FileName.c_str(),std::ios::app);
        
        if(!outFile) 
        {
            G4cout<<"//////////////Error opening outFile/////////////////"<<G4endl;
            return;
        }
        */
        
                
        
        bool EJ301_Triggered = false;  //--- to verify the coincidence ---
        
        G4int width = 10, precision = 4;
        
        G4int           EventID;
        G4int           ParentID;
        G4int           TrackID;	
        G4int           StepNumber;
        G4ThreeVector   Position;
        G4double        dKinEnergy;
        G4double        dDepositeEnergy;
        G4String        Particle;
        G4String        CreatorProcess;
        G4String        DepositProcess;
        G4double        globalTime;

        
        
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
        pAnalysisManager->dKinEnergy.clear();
        pAnalysisManager->globalTime.clear();
        pAnalysisManager->Particle.clear();
        pAnalysisManager->CreatorProcess.clear();
        pAnalysisManager->DepositProcess.clear();
        
        
        
        
        std::vector<int>            tempParentID;
        std::vector<int>            tempTrackNbr;
        std::vector<int>             tempStepNbr;
        std::vector<double>           tempPositionX;
        std::vector<double>           tempPositionY;
        std::vector<double>           tempPositionZ;
        std::vector<double>      tempDepositeEnergy;
        std::vector<double>           tempKinEnergy;
        std::vector<double>          tempGlobalTime;
        std::vector<std::string>                tempParticle;
        std::vector<std::string>          tempCreatorProcess;
        std::vector<std::string>          tempDepositProcess;
        
        
        std::vector<int>            SortedParentID;
        std::vector<int>            SortedTrackNbr;
        std::vector<int>             SortedStepNbr;
        std::vector<double>           SortedPositionX;
        std::vector<double>           SortedPositionY;
        std::vector<double>           SortedPositionZ;
        std::vector<double>      SortedDepositeEnergy;
        std::vector<double>          SortedKinEnergy;
        std::vector<double>          SortedGlobalTime;
        std::vector<std::string>                SortedParticle;
        std::vector<std::string>                SortedCreatorProcess;
        std::vector<std::string>                SortedDepositProcess;
        
        
        
        tempParentID.clear();
        tempTrackNbr.clear();
        tempStepNbr.clear();
        tempPositionX.clear();
        tempPositionY.clear();
        tempPositionZ.clear();
        tempDepositeEnergy.clear();
        tempKinEnergy.clear();
        tempGlobalTime.clear();
        tempParticle.clear();
        tempCreatorProcess.clear();
        tempDepositProcess.clear();
        
        
        
        SortedParentID.clear();
        SortedTrackNbr.clear();
        SortedStepNbr.clear();
        SortedPositionX.clear();
        SortedPositionY.clear();
        SortedPositionZ.clear();
        SortedDepositeEnergy.clear();
        SortedGlobalTime.clear();
        SortedParticle.clear();
        SortedCreatorProcess.clear();
        SortedDepositProcess.clear();
        
    //-----------------------------------------------------    


        
           
        //------------- to check the coincidence with EJ301 neutron detector -------------------     
        for(int i=0; i<numberOfHits; i++)
        {
            BeamTestHit* aHit = (*hitsCollection)[i];
            
            Position = aHit->GetPosition()/mm;
            dDepositeEnergy = aHit->GetEnergyDeposited()/keV;
            
            if( ( (Position.y()<-300)||(Position.y()>300) ) && (dDepositeEnergy>0.1) )
                EJ301_Triggered = true;   //--- hit in the coincident EJ301 ---
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
                    dKinEnergy = aHit->GetKineticEnergy()/keV;
                      Particle = aHit->GetParticleType();
                CreatorProcess = aHit->GetCreatorProcess();
                DepositProcess = aHit->GetDepositingProcess();
                    globalTime = aHit->GetTime();
                
                
                
                    if(dDepositeEnergy>0)    //== the 0.01keV is threshold of single step enegy despoit for step recording ====  
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
                            <<Particle<<"    "<<std::setw(width)	
                            <<DepositProcess<<"    "<<std::setw(width)	
                            <<CreatorProcess
                            <<G4endl;
                             */
                    
                            tempParentID.push_back(ParentID);
                            tempTrackNbr.push_back(TrackID);
                            tempStepNbr.push_back(StepNumber);
                            tempPositionX.push_back(Position.x());
                            tempPositionY.push_back(Position.y());
                            tempPositionZ.push_back(Position.z());
                            tempDepositeEnergy.push_back(dDepositeEnergy);
                            tempKinEnergy.push_back(dKinEnergy);
                            tempGlobalTime.push_back(globalTime);
                            tempParticle.push_back(Particle);
                            tempCreatorProcess.push_back(CreatorProcess);
                            tempDepositProcess.push_back(DepositProcess);
                             
                    
                    }
                
               
                
            }
            
            
            //-------------------Sort the steps by global time--------------------------------
            
                        double                        TempTime; /* Temporary value for swapping two array */
                           int           TempElementOrder = -1;
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
                
                
                            SortedParticle.push_back(tempParticle.at(element));
                            SortedCreatorProcess.push_back(tempCreatorProcess.at(element));
                            SortedDepositProcess.push_back(tempDepositProcess.at(element));
                
                            //std::cout<<element<<"    "<<tempParentID.at(v)<<"    "<<SortedParentID.at(v)<<std::endl;
                    }
            
            
            //------------------------------------------------------------------------------------------     
            
            pAnalysisManager->ParentID = SortedParentID;
            pAnalysisManager->TrackNbr = SortedTrackNbr;
            pAnalysisManager->StepNbr = SortedStepNbr;
            pAnalysisManager->PositionX = SortedPositionX;
            pAnalysisManager->PositionY = SortedPositionY;
            pAnalysisManager->PositionZ = SortedPositionZ;
            pAnalysisManager->dDepositeEnergy = SortedDepositeEnergy;
            pAnalysisManager->dKinEnergy = SortedKinEnergy;
            pAnalysisManager->globalTime = SortedGlobalTime;
            pAnalysisManager->Particle = SortedParticle;
            pAnalysisManager->CreatorProcess = SortedCreatorProcess;
            pAnalysisManager->DepositProcess = SortedDepositProcess;
            
            
            pAnalysisManager->_DataTTree->Fill();

            
        }  //--- end of if(EJ301_Triggered) ---
        
        
        //outFile.close();    
        
    }  //--- end of if(fHitsCollectionID>=0) ---
//====================================================================================================================
    
    


    
    
    
} // --- end of the EndOfEvent function


