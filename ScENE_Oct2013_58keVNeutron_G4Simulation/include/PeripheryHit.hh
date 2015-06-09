//#ifndef PERIPHERYHIT_HH
//#define PERIPHERYHIT_HH

#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4THitsCollection.hh"
#include "G4ParticleDefinition.hh"
#include "G4VHit.hh"

//#include "UCLANeutronHit.hh"

class PeripheryHit : public G4VHit {

public:
  
  // Constructors
  PeripheryHit();

  // Destructor
  virtual ~PeripheryHit();
  
  inline void *operator new(size_t);
  inline void operator delete(void *aHit);

  // Methods
  virtual void Draw();
  virtual void Print();

  // Incidence Information 
    inline	void SetTrackID(G4int iTrackID) { m_iTrackID = iTrackID; };
    inline	void SetParentID(G4int iParentID) { m_iParentID = iParentID; };
    inline void SetStepNumber(G4int StepN) {m_iStepNumber = StepN;}
    
    inline	void SetParticleType(const G4String &hParticleType) { m_pParticleType = new G4String(hParticleType); }
    //inline	void SetParentType(const G4String &hParentType) { m_pParentType = new G4String(hParentType); }
    inline	void SetCreatorProcess(const G4String &hProcessName) { m_pCreatorProcess = new G4String(hProcessName); }
    inline	void SetDepositingProcess(const G4String &hProcessName) { m_pDepositingProcess = new G4String(hProcessName); }
    inline	void SetPhysicalVolumeName(const G4String &hVolumeName) {m_pPhysicalVolumeName = new G4String(hVolumeName);}
    
    inline	void SetPosition(G4ThreeVector hPosition) { m_hPosition = hPosition; };
    inline	void SetEnergyDeposited(G4double dEnergyDeposited) { m_dEnergyDeposited = dEnergyDeposited; };
    inline	void SetKineticEnergy(G4double dKineticEnergy) { m_dKineticEnergy = dKineticEnergy; };
    inline	void SetTime(G4double dTime) { m_dTime = dTime; };

    
    
    
    
    
    
    inline	G4int GetTrackID() { return m_iTrackID; };
    inline	G4int GetParentID() { return m_iParentID; };
    inline G4int GetStepNumber()  {return m_iStepNumber;}
    inline	G4String &GetParticleType() { return *m_pParticleType; }
    inline	G4String &GetVolumeName() { return *m_pPhysicalVolumeName; }
    inline	G4String &GetCreatorProcess() { return *m_pCreatorProcess; }
    inline	G4String &GetDepositingProcess() { return *m_pDepositingProcess; }
    inline	G4ThreeVector GetPosition() { return m_hPosition; };
    inline	G4double GetEnergyDeposited() { return m_dEnergyDeposited; };      
    inline	G4double GetKineticEnergy() { return m_dKineticEnergy; };      
    inline	G4double GetTime() { return m_dTime; };      
    
    
private:
    
    G4int m_iTrackID;
    G4int m_iParentID;
    G4int m_iStepNumber;
    
    G4String *m_pParticleType;
    G4String *m_pPhysicalVolumeName;
    G4String *m_pCreatorProcess;
    G4String *m_pDepositingProcess;
    
    
    G4ThreeVector m_hPosition;
    G4double m_dEnergyDeposited;
    G4double m_dKineticEnergy;
    G4double m_dTime;


};

typedef G4THitsCollection<PeripheryHit> PeripheryHitsCollection;

extern G4Allocator<PeripheryHit> PeripheryHitAllocator;

inline void* PeripheryHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*)PeripheryHitAllocator.MallocSingle();
  return aHit;
}

inline void PeripheryHit::operator delete(void* aHit)
{
   PeripheryHitAllocator.FreeSingle((PeripheryHit*) aHit);
}

//#endif
