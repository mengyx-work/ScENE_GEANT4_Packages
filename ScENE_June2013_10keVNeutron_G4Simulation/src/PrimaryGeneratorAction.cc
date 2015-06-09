#include "PrimaryGeneratorAction.hh"
#include "stdlib.h"
#include <stdio.h>
#include <time.h>
#include <fstream>
#include <iostream>
#include "G4Event.hh"


#include "G4GeneralParticleSource.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"

#include "TTree.h"
#include "TMath.h"
#include "TRandom1.h"
#include "TGraph.h"

//------
/*
 Comments on the Primary Generator:
 the various default functions in use.
 
 1. the PrimaryGeneratorAction inheriated function is used to 
 define particle source properies during RunAction. 
 Properies which are the same during the whole run
 
 2. the generatePrimaries. Give a each individual particle the
 corresponding initial values.
 
*/
//------

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    
  //---General Particle Source------
    /*  
  particleSource = new G4GeneralParticleSource();
  G4cout<<"General source used, so the input Neutron Energy is NOT used "<<G4endl;
 
    //particleGun = new G4ParticleGun();
    //G4cout<<"particle used"<<G4endl;
  */
    
    
    
    //---zAxis Position is basedon: const G4double ActiveLAr_Volume_Height     = 3.*25.4*mm -----
    
    
    const      double     Neutron_xAxis_Position = -66*cm;  //---Notre Dame June 2013 setup 66 cm-----
    const      double     Neutron_zAxis_Position = (31.4-1.5*2.54)*cm;
    
 
    
    G4int n_particle = 1;
    particleGun  = new G4ParticleGun(n_particle);
    
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName;
    
    G4ParticleDefinition* particle = particleTable->FindParticle(particleName="neutron");
    
    particleGun->SetParticleDefinition(particle);

    particleGun->SetParticlePosition(G4ThreeVector(Neutron_xAxis_Position, 0., Neutron_zAxis_Position));
    
  }


PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete particleGun;
    //delete particleSource;
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    
    //===============================================================

    //------------Random Seed Generation for Neutron Energy and Direction-------------
    
    const      double                          k = 1.7474;  //--- fitting k from ATOMIC DATA AND NUCLEAR DATA TABLE Ep = 2.35 MeV----
    const      double                  NormConst = k/(TMath::Exp(k) - TMath::Exp(-k));
    const      double                    PDF_Max = NormConst*TMath::Exp(k);
  
    
    TRandom1  *RandomTheta = new TRandom1(0);
    TRandom1    *RandomPhi = new TRandom1(0);
    TRandom1  *RandomLimit = new TRandom1(0);
    
    bool  Accepted = false;
    
    double       randomSeed;
    double      randomLimit;
    double    NeutronEnergy;
    double            Theta;
    double         CosTheta;
    double         SinTheta;
    
    double   Phi = 2*RandomPhi->Uniform()*TMath::Pi();
    
    
    while(!Accepted)
    {
        randomSeed = ( 1.0 - RandomTheta->Uniform() );   //--Uniform excludes point 0 ---
        //cout<<"the random seed in loop: "<<randomSeed<<endl;
        randomLimit = RandomLimit->Uniform();
        
        //--- theta from 90 to 0 deg ------
        Theta = randomSeed*TMath::Pi()/2;
        //---------------------------------
        
        CosTheta = TMath::Cos(Theta);
        SinTheta = TMath::Sin(Theta);
        
        //----default pdf function f(x) = A*exp(k*x) ---
        if( randomLimit*PDF_Max<=NormConst*TMath::Exp(k*CosTheta) )
            Accepted = true;
        
    }
    
    
    double positionX = CosTheta;
    double positionY = SinTheta*TMath::Cos(Phi);
    double positionZ = SinTheta*TMath::Sin(Phi);
    
    

    
    const int ArraySize = 21;
    
    double CosThetaValue[ArraySize] ={1, 0.984, 0.961, 0.94, 0.925, 0.904, 0.8713, 0.8227, 0.7922, 0.7585, 0.7162, 0.6423, 0.5853, 0.5383, 0.4807, 0.4048, 0.3149, 0.2447, 0.156, 0.0659, 0};
    
    double NeutronEnergyValue[ArraySize] = {0.6193, 0.6129, 0.6037, 0.5955, 0.59, 0.5817, 0.568, 0.5515, 0.5406, 0.5287, 0.514, 0.4893, 0.471, 0.4564, 0.439, 0.417, 0.3923, 0.374, 0.352, 0.331, 0.316};

    
    TGraph *graph = new TGraph(ArraySize, CosThetaValue, NeutronEnergyValue);
    
    NeutronEnergy = graph->Eval(CosTheta);
    
    
    
    delete RandomTheta;
    delete RandomPhi;
    delete RandomLimit;
    delete graph;

    /*
    std::ofstream output("direction.txt",std::ios::app);
    output<<positionX<<"  "<<positionY<<"  "<<positionZ<<G4endl;
    */
    
    
    //-----------------Particle Gun Settings-----------------
    ///*  
    particleGun->SetParticleMomentumDirection(G4ThreeVector(positionX, positionY, positionZ));
    particleGun->SetParticleEnergy(NeutronEnergy*MeV);   
        
    //*/
   
    
    //particleSource->GeneratePrimaryVertex(anEvent);
    particleGun->GeneratePrimaryVertex(anEvent);
    
    
    
    //===================================================================================   
  
    
    
    /*
    //------- To Generate Multiple Vertex (track) within a single event ------
     
     // 1. repeatedly add new PrimaryVertex into the particle gun
    for(int k=0; k<4; k++)
    {
     particleGun->SetParticleMomentumDirection(G4ThreeVector(positionX, positionY, positionZ));
     particleGun->SetParticleEnergy(NeutronEnergy*MeV);
     
     particleGun->GeneratePrimaryVertex(anEvent);
        
    }
     
     //G4PrimaryVertex *Vertex = anEvent->GetPrimaryVertex();
     
     //G4int      particleNbr = Vertex->GetNumberOfParticle();
     
     // 2. to look up the specific Vertex#
     G4int VertexNbr = anEvent->GetNumberOfPrimaryVertex();
     
     // 3. for each vertex, to set the correct globaltime
     for(int v=0; v<VertexNbr; v++)
     {
     
     G4PrimaryVertex *Vertex = anEvent->GetPrimaryVertex(v);
     Vertex->SetT0(v*100);
     }
     
     //G4cout<<"The Real Vertex# is: "<<anEvent->GetNumberOfPrimaryVertex()<<G4endl;
     */
    //===================================================================================
    
    


    
    /*
    G4double  Vertex_Time = 10;

    G4ThreeVector Vertex_Position;
    
    
    Vertex_Position.SetX(-5*cm);
    Vertex_Position.SetY(0.*cm);
    Vertex_Position.SetZ(0.*cm);

    G4PrimaryVertex *Vetex = new G4PrimaryVertex(Vertex_Position, Vertex_Time);
    
    
    delete  Vertex;
    */
    
    
    //------------------------------------------------------------------------------------
    
    G4int EventID = anEvent->GetEventID();
    
    const G4int UnitNumber = 50000;
    
    if(EventID%UnitNumber==0)
    {
        time_t now = time(0);
        //G4cout<<"Current time is: "<<ctime(&now)<<G4endl;
        G4cout<<"Event# "<<EventID<<" is simulated at "<<ctime(&now)<<G4endl;
        
    }
   //------------------------------------------------------------------------------------

}
