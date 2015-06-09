//#ifndef PERIPHERYSENSITIVEDETECTOR_HH
//#define PERIPHERYSENSITIVEDETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "PeripheryHit.hh"
//#include  "BeamTestHit.hh"

class G4HCofThisEvent;
class G4TouchableHistory;
class G4Step;
class PeripheryHit;

class PeripherySensitiveDetector : public G4VSensitiveDetector {

public:

  // Constructor
  PeripherySensitiveDetector(const G4String& name);

  // Destructor
  virtual ~PeripherySensitiveDetector();
  
  // Methods
  virtual void Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent);
  
  virtual G4bool ProcessHits(G4Step* pStep,G4TouchableHistory* history);

  virtual void EndOfEvent(G4HCofThisEvent* hitsCollectionOfThisEvent);
  
private:
  
  // Data members
  PeripheryHitsCollection*           pHitsCollection;
  G4int                     PeripheryHitsCollectionID;
  //BeamtTestHit *hit;

};



//#endif

