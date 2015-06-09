#ifndef PRIMARYGENERATORACTION_HH
#define RIMARYGENERATORACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"


class G4ParticleGun;
//class G4GeneralParticleSource;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {

public:

  // Constructor
  PrimaryGeneratorAction();    

  // Destructor
  virtual ~PrimaryGeneratorAction();
  
  // Methods
  void GeneratePrimaries(G4Event*);
  

private:

  //G4GeneralParticleSource* particleSource;
  G4ParticleGun* particleGun;

};

#endif


