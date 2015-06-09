#ifndef __ANALYSISMANAGER_H__
#define __ANALYSISMANAGER_H__

#include <vector>
#include <string>
#include <iostream>

#include "globals.hh"
#include "string.h"
#include "TFile.h"
#include "TTree.h"

class TFile;
class TTree;

class AnalysisManager{
    
public:
    void     Initialize();
    void     Finalize();   
    
     ~AnalysisManager();

 void      SetDataFileName(const G4String FileName){DataFileName = FileName;}
 G4String  GetDataFileName(){return DataFileName;}


public:
    
        int EventID;
        int InteractingVolumeInt;
    
    std::vector<int>                     ParentID;
    std::vector<int>                     TrackNbr;
    std::vector<int>                      StepNbr;
    std::vector<double>                 PositionX;
    std::vector<double>                 PositionY;
    std::vector<double>                 PositionZ;
    std::vector<double>                DirectionX;
    std::vector<double>                DirectionY;
    std::vector<double>                DirectionZ;
    


    std::vector<double>     dDepositeEnergy;
    std::vector<double>          dKinEnergy;
    std::vector<double>          globalTime;
    
    std::vector<std::string>          PhysicalVolume;
    std::vector<std::string>          ParentParticle;
    std::vector<std::string>                Particle;
    std::vector<std::string>          CreatorProcess;
    std::vector<std::string>          DepositProcess;
    
    //std::vector<int>
    //std::vector<int>
    
    
    TFile *_DataTFile;
    TTree *_DataTTree;
    
    
    
private:
    G4String DataFileName;
    
   };


#endif __ANALYSISMANAGER_H__
