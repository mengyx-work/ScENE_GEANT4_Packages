
#ifndef BEAMTESTEMCALORIMETER_HH
#define BEAMTESTEMCALORIMETER_HH

#include "G4VSensitiveDetector.hh"
#include "BeamTestEmCalorimeterHit.hh"

class G4HCofThisEvent;
class G4TouchableHistory;
class G4Step;

class BeamTestEmCalorimeter : public G4VSensitiveDetector {

public:

  // Constructor
  BeamTestEmCalorimeter(const G4String& name);

  // Destructor
  virtual ~BeamTestEmCalorimeter();
  
  // Methods
  virtual void Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent);
  
  virtual G4bool ProcessHits(G4Step* aStep,G4TouchableHistory* history);

  virtual void EndOfEvent(G4HCofThisEvent* hitsCollectionOfThisEvent);
  
private:
  
  // Data members
  BeamTestEmCalorimeterHitsCollection* fHitsCollection;
  G4int fHitsCollectionID;

};




#endif

