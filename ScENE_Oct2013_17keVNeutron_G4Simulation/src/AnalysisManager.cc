#include "AnalysisManager.hh"

void AnalysisManager::Initialize()
{
    _DataTFile = new TFile(DataFileName.c_str(), "RECREATE");
    _DataTTree = new TTree("G4Data", "G4 Raw Step Data");
    
    //gROOT->ProcessLine("#include <vector>");
   /* 
    std::vector<int> ParentID;
    std::vector<int> TrackNbr;
    std::vector<int>  StepNbr;
    */

    
    _DataTTree->Branch("EventID",                                                     &EventID);
    _DataTTree->Branch("InteractingVolumeInt",                           &InteractingVolumeInt);
    _DataTTree->Branch("ParentID",                                                   &ParentID);
    _DataTTree->Branch("TrackNbr",                                                   &TrackNbr);
    _DataTTree->Branch("StepNbr",                                                     &StepNbr);
    _DataTTree->Branch("PositionX",                                                 &PositionX);
    _DataTTree->Branch("PositionY",                                                 &PositionY);
    _DataTTree->Branch("PositionZ",                                                 &PositionZ);
    _DataTTree->Branch("DirectionX",                                               &DirectionX);
    _DataTTree->Branch("DirectionY",                                               &DirectionY);
    _DataTTree->Branch("DirectionZ",                                               &DirectionZ);
    _DataTTree->Branch("dDepositeEnergy",                                     &dDepositeEnergy);
    _DataTTree->Branch("dKinEnergy",                                               &dKinEnergy);
    _DataTTree->Branch("globalTime",                                               &globalTime);
    _DataTTree->Branch("ParentParticle",                                       &ParentParticle);
    _DataTTree->Branch("Particle",                                                   &Particle);
    _DataTTree->Branch("CreatorProcess",                                       &CreatorProcess);
    _DataTTree->Branch("DepositProcess",                                       &DepositProcess);
    _DataTTree->Branch("PhysicalVolume",                                       &PhysicalVolume);
    
    _DataTTree->AutoSave();
                   
}


void AnalysisManager::Finalize()
{
    _DataTFile->Write();
    _DataTFile->Close();
}