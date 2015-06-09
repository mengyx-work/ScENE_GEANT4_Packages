
#include "BeamTestHit.hh"

G4Allocator<BeamTestHit> BeamTestHitAllocator;

BeamTestHit::BeamTestHit()
{ 
    m_iTrackID = -1;
    m_iParentID = -1;
    m_iStepNumber = -1;
    m_hPosition = G4ThreeVector(0,0,0);
    m_hDirection = G4ThreeVector(0,0,0);
    m_dEnergyDeposited = -1;
    m_dKineticEnergy = -1;
    m_dTime = -1;
}



BeamTestHit::~BeamTestHit()
{
    if(m_pParticleType)       delete m_pParticleType;
    if(m_pCreatorProcess)     delete m_pCreatorProcess;
	if(m_pDepositingProcess)  delete m_pDepositingProcess;
    if(m_pParentType)         delete m_pParentType;
    if(m_pPhysicalVolumeName) delete m_pPhysicalVolumeName;
}



void BeamTestHit::Draw() {}



void BeamTestHit::Print()
{
	/*
   G4cout << "Incidence Particle Name and Kinetic Energy " << fIPD->GetParticleName() 
          << " " << fIKEnergy/MeV << " MeV" << G4endl; 
   G4cout << "Incidence position in Target monitor " << fIPosition/cm << " in cm" << G4endl<<G4endl; 
   //G4cout << "Incidence Direction " << fIMomentumD << G4endl<<G4endl;
     */
}
