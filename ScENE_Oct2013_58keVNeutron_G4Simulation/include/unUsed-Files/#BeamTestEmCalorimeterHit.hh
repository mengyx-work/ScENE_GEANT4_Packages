
#ifndef BEAMTESTEMCALORIMETERHIT_HH
#define BEAMTESTEMCALORIMETERHIT_HH

#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4THitsCollection.hh"
#include "G4Transform3D.hh"
#include "G4VHit.hh"

class G4AttDef;
class G4AttValue;

class BeamTestEmCalorimeterHit : public G4VHit {

public:
  
  // Constructors
  BeamTestEmCalorimeterHit();
  BeamTestEmCalorimeterHit(G4int id);

  // Destructor
  virtual ~BeamTestEmCalorimeterHit();
  
  inline void *operator new(size_t);
  inline void operator delete(void *aHit);

  // Methods
  //virtual void Draw();

  virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
  virtual std::vector<G4AttValue>* CreateAttValues() const;

  virtual void Print();

  // Deposited energy
  inline void AddDepositedEnergy(G4double energy) {fDepositedEnergy += energy;}
  inline G4double GetDepositedEnergy() const {return fDepositedEnergy;}

  // Position vector
  inline void SetPosition(G4ThreeVector position) {fPosition = position;}
  inline G4ThreeVector GetPosition() const {return fPosition;}

  // Rotation matrix
  inline void SetRotation(G4RotationMatrix rotation) {fRotation = rotation;}
  inline G4RotationMatrix GetRotation() const {return fRotation;}

  // Logical volume
  inline void SetLogicalVolume(G4LogicalVolume* volume) {pLogicalVolume = volume;}
  inline const G4LogicalVolume* GetLogicalVolume() const {return pLogicalVolume;}
  
private:
  
  // Data members
  G4int fCellID;
  G4double fDepositedEnergy;
  G4ThreeVector fPosition;
  G4RotationMatrix fRotation;
  const G4LogicalVolume* pLogicalVolume;
  
};

typedef G4THitsCollection<BeamTestEmCalorimeterHit> BeamTestEmCalorimeterHitsCollection;

extern G4Allocator<BeamTestEmCalorimeterHit> BeamTestEmCalorimeterHitAllocator;

inline void* BeamTestEmCalorimeterHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*)BeamTestEmCalorimeterHitAllocator.MallocSingle();
  return aHit;
}

inline void BeamTestEmCalorimeterHit::operator delete(void* aHit)
{
  BeamTestEmCalorimeterHitAllocator.FreeSingle((BeamTestEmCalorimeterHit*) aHit);
}

#endif


