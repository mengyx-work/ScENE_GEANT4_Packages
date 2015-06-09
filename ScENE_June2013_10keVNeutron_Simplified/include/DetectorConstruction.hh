
#ifndef DETECTORCONSTRUCTION_HH 
#define DETECTORCONSTRUCTION_HH 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction {

public:
  
  // Constructor
  DetectorConstruction();
  
  // Destructor
  virtual ~DetectorConstruction();
  
  // Method
  virtual G4VPhysicalVolume* Construct();
  
private:

  // Helper methods
  void DefineMaterials();
  void SetupGeometry();
  
  // World logical and physical volumes
  G4LogicalVolume*   fpWorldLogical;
  G4VPhysicalVolume* fpWorldPhysical;
  
};

#endif

