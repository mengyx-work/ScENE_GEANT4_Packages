
#ifndef EVENTACTION_HH
#define EVENTACTION_HH

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "string.h"

#include "AnalysisManager.hh"

class EventAction : public G4UserEventAction {

public:

  // Constructor
  EventAction(AnalysisManager* AnalysisManagerPointer);

  // Destructor
  ~EventAction();

  // Metohds
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);
  
private:

  // Data member
  G4int             fHitsCollectionID;
  G4int     PeripheryHitsCollectionID;
    
  G4String                  FileName;
 AnalysisManager   *pAnalysisManager;    
  
};

#endif

    
