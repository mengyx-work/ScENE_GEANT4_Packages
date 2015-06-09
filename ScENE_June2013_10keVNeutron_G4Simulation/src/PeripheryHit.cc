#include "PeripheryHit.hh"

G4Allocator<PeripheryHit> PeripheryHitAllocator;

PeripheryHit::PeripheryHit()
{ 
    m_iTrackID = -1;
    m_iParentID = -1;
    m_iStepNumber = -1;
    m_hPosition = G4ThreeVector(0,0,0);
    m_dEnergyDeposited = -1;
    m_dKineticEnergy = -1;
    m_dTime = -1;
}



PeripheryHit::~PeripheryHit()
{
    if(m_pParticleType)       delete m_pParticleType;
    if(m_pCreatorProcess)     delete m_pCreatorProcess;
	if(m_pDepositingProcess)  delete m_pDepositingProcess;
    if(m_pPhysicalVolumeName) delete m_pPhysicalVolumeName;
}



void PeripheryHit::Draw() {}



void PeripheryHit::Print()
{
	/*
   G4cout << "Incidence Particle Name and Kinetic Energy " << fIPD->GetParticleName() 
          << " " << fIKEnergy/MeV << " MeV" << G4endl; 
   G4cout << "Incidence position in Target monitor " << fIPosition/cm << " in cm" << G4endl<<G4endl; 
   //G4cout << "Incidence Direction " << fIMomentumD << G4endl<<G4endl;
     */
}
