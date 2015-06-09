//#ifndef BEAMTESTSILICONMONITOR_HH
//#define BEAMTESTSILICONMONITOR_HH

#include "G4VSensitiveDetector.hh"
#include "BeamTestHit.hh"

class G4HCofThisEvent;
class G4TouchableHistory;
class G4Step;

class BeamTestSensitiveDetector : public G4VSensitiveDetector {

public:

  // Constructor
  BeamTestSensitiveDetector(const G4String& name);

  // Destructor
  virtual ~BeamTestSensitiveDetector();
  
  // Methods
  virtual void Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent);
  
  virtual G4bool ProcessHits(G4Step* pStep,G4TouchableHistory* history);

  virtual void EndOfEvent(G4HCofThisEvent* hitsCollectionOfThisEvent);
  
private:
  
  // Data members
  BeamTestHitsCollection* fHitsCollection;
  G4int fHitsCollectionID;
    
    std::map<int, G4String> m_ParentParticleType;

};



//#endif

