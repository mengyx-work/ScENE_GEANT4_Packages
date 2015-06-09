
#ifndef BEAMTESTCELLPARAMETERISATION_HH
#define BEAMTESTCELLPARAMETERISATION_HH

#include "G4VPVParameterisation.hh"
#include "globals.hh"
#include <vector>

class G4VPhysicalVolume;

class BeamTestCellParameterisation : public G4VPVParameterisation {

public:

  // Constructor
  BeamTestCellParameterisation();

  // Destructor
  virtual ~BeamTestCellParameterisation();

  // Method
  virtual void ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const;
  
private:

  // Data members
  std::vector<G4double> xCell;
  std::vector<G4double> yCell;

};

#endif


