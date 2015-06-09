//#include "DetectorConstruction.hh"
#include "SCENE_TwoPhasePrototype_CoincidentGeometry.hh"
//#include "SimpDetectorConstruction.hh"
#include "EventAction.hh"
#include "RunAction.hh"
//#include "StackingAction.hh"
//#include "SteppingAction.hh"
#include "TrackingAction.hh"




#include "AnalysisManager.hh"
//#include "Xenon100PhysicsList.hh"
//#include "Shielding.hh"
//#include "QGSP_BERT_HP.hh"
//#include "QGSP_BIC_HP.hh"
#include "LBE.hh"

#include <sys/time.h>
#include <iostream>

#include "PrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "Randomize.hh"


#ifdef G4UI_USE_TCSH
#include "G4UItcsh.hh"
#endif
#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

//    using namespace std;

int main(int argc,char** argv) {
    
    /*
    if (argc!=3) 
    {
        G4cout<<"Abort! some parameter is missing!"<<G4endl;
        G4cout<<"Standard Command: ./Executalble + Macro + DataFileName "<<G4endl;
        return 0;
    }
    */
    
    
    // Choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    
    unsigned long int randomeSeed= time(NULL);

    //cout<<"the sec output is: "<<sec<<endl;
    
    CLHEP::HepRandom::setTheSeed(randomeSeed);

  
    
    
  // Run manager
  G4RunManager * runManager = new G4RunManager;
  
  // Mandatory initialization classes
  DetectorConstruction* detector = new DetectorConstruction;
  runManager->SetUserInitialization(detector);
  //runManager->SetUserInitialization(new QGSP_BIC_HP); 
  runManager->SetUserInitialization(new LBE);  
  //runManager->SetUserInitialization(new Xenon100PhysicsList);
  //runManager->SetUserInitialization(new Shielding);
  //runManager->SetUserInitialization(new QGSP_BERT_HP);
   
 

    
    
    G4UIsession* session=0;
    
    if (argc==1)   // Define UI session for interactive mode.
    {
        // G4UIterminal is a (dumb) terminal.
#if defined(G4UI_USE_XM)
        session = new G4UIXm(argc,argv);
#elif defined(G4UI_USE_TCSH)
        session = new G4UIterminal(new G4UItcsh);      
#else
        session = new G4UIterminal();
#endif
    }
    
    /*
     #ifdef G4VIS_USE
     // Visualization manager
     G4VisManager* visManager = new G4VisExecutive;
     visManager->Initialize();
     #endif
     */
    
    AnalysisManager *pAnalysisManager = new AnalysisManager;
    
    pAnalysisManager->SetDataFileName(argv[2]);
    
    
    // User action classes
    runManager->SetUserAction(new PrimaryGeneratorAction);
    runManager->SetUserAction(new RunAction(pAnalysisManager));
    runManager->SetUserAction(new EventAction(pAnalysisManager));
    //runManager->SetUserAction(new SteppingAction);
    runManager->SetUserAction(new TrackingAction);
    //runManager->SetUserAction(new TrackingAction(pAnalysisManager));
                              
    
    // Initialize G4 kernel
    runManager->Initialize();
    
    G4UImanager* UI = G4UImanager::GetUIpointer();  
    
    
    
#ifdef G4VIS_USE
    // Visualization manager
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
#endif
    
    
    
    
    if (session)   // Define UI session for interactive mode.
    {
        session->SessionStart();
        delete session;
    }
    else           // Batch mode
    { 
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UI->ApplyCommand(command+fileName);
    }
    
    
    // Job termination
#ifdef G4VIS_USE
    delete visManager;
#endif
    
    delete runManager;
    //wdelete pAnalysisManager;
    
    return 0;
}

