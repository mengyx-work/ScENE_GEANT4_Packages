#include "RunAction.hh"
#include <iostream>
#include <fstream>
//#include  <Randomize.h>

/*
RunAction::RunAction(){
    
}
*/

RunAction::~RunAction(){
    
}


void    RunAction::BeginOfRunAction(const G4Run *pRun)
{
    std::cout<<"Start of Run!"<<std::endl;
    pAnalysisManager->Initialize();
    
    /*
    std::ofstream       outTrack("TrackOutput.txt", std::ios::trunc);
    std::ofstream       outStep("StepOutput.txt", std::ios::trunc);
    
    outTrack<<"  EventID     ParentID    TrackID     PositionX       PositionY       PositionZ       KineticEnergy       GlobalTime      Particle        CreatorProcess      VolumName"<<G4endl;
    
    outStep<<"  EventID     ParenID      TrackID      StepID       PositionX       PositionY       PositionZ       KineticEnergy      DepositEnergy       GlobalTime      Particle        CreatorProcess      VolumName"<<G4endl;
    
    outTrack.close();
    outStep.close();
    */
    
}


void   RunAction::EndOfRunAction(const G4Run *pRun)
{
    std::cout<<"End of Run!"<<std::endl;
    pAnalysisManager->Finalize();
    
}