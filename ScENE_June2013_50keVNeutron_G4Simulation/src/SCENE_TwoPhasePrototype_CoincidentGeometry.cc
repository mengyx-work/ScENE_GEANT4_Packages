    
//#include "DetectorConstruction.hh" 
#include "SCENE_TwoPhasePrototype_CoincidentGeometry.hh"
#include "BeamTestSensitiveDetector.hh"
#include "PeripherySensitiveDetector.hh"

#include "G4Polyhedra.hh"
#include "G4Polycone.hh"

#include "G4Tubs.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"     
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4VisAttributes.hh"

#include "G4SDManager.hh"
#include "stdlib.h"
#include "string.h"
#include <cmath>




DetectorConstruction::DetectorConstruction()
  :fpWorldLogical(0)
  ,fpWorldPhysical(0)
{}

DetectorConstruction::~DetectorConstruction() {}



G4VPhysicalVolume*DetectorConstruction::Construct()
{
  // Material Definition
  DefineMaterials();  

  // Geometry Definition
  SetupGeometry();   

  // Return world volume
  return fpWorldPhysical;  
}







void DetectorConstruction::DefineMaterials()
{
  G4String symbol;             
  G4double a, z, density;     
  G4int ncomponents;
  G4double fractionmass;
  
    //================================== elements ===================================
  new G4Material(     "Silicon",    z=14., a=28.0855*g/mole,  density=2.330*g/cm3);
  new G4Material(        "Iron",    z=26., a=55.845*g/mole,   density=7.87*g/cm3);
  new G4Material(    "Titanium",    z=22., a=47.90*g/mole,    density=4.540*g/cm3);  
  new G4Material( "LiquidArgon",    z=18,  a=39.948*g/mole,   density=1.43*g/cm3);
  new G4Material(    "GasArgon",    z=18,  a=39.948*g/mole,   density=1.449e-3*g/cm3);
  new G4Material(      "Copper",    z=29,  a=63.546*g/mole,   density=8.94*g/cm3);
  new G4Material(   "Aluminium",    z=13,  a=26.98*g/mole,    density=2.70*g/cm3);
  new G4Material(    "Tungsten",    z=74,  a=183.85*g/mole,   density=19.25*g/cm3);
  
    
	G4Element *H  = new G4Element("Hydrogen",  "H",  1.,  1.0079*g/mole);
	G4Element *C  = new G4Element("Carbon",    "C",  6.,  12.011*g/mole);
	G4Element *N  = new G4Element("Nitrogen",  "N",  7.,  14.007*g/mole);
	G4Element *O  = new G4Element("Oxygen",    "O",  8.,  15.999*g/mole);
	G4Element *F  = new G4Element("Fluorine",  "F",  9.,  18.998*g/mole);
	G4Element *Si = new G4Element("Silicon",   "Si", 14., 28.086*g/mole);
	G4Element *Cr = new G4Element("Chromium",  "Cr", 24., 51.996*g/mole);
	G4Element *Mn = new G4Element("Manganese", "Mn", 25., 54.938*g/mole);
	G4Element *Fe = new G4Element("Iron",      "Fe", 26., 55.85*g/mole);   //-- density 7.874*g/cm3 --
	G4Element *Ni = new G4Element("Nickel",    "Ni", 28., 58.693*g/mole);  //-- density 8.9*g/cm3 --
    G4Element *Mo = new G4Element("Nickel",    "Ni", 42., 95.96*g/mole);   //-- density 10.28*g/cm3 --

    //G4Element *Al = new G4Element("Aluminium", "Al", 13., 26.982*g/mole);
	//G4Element *Xe = new G4Element("Xenon",     "Xe", 54., 131.293*g/mole);
    //G4Element *Cu = new G4Element("Copper",    "Cu", 29., 63.546*g/mole);
	//G4Element *Pb = new G4Element("Lead",      "Pb", 82., 207.2*g/mole);
	//G4Element *Gd = new G4Element("Gadolinium","Gd", 64., 157.25*g/mole);

// PTFE
    G4Material* Teflon = new G4Material("Teflon", 2.2*g/cm3, 2, kStateSolid);
	Teflon->AddElement(C, 0.240183);
	Teflon->AddElement(F, 0.759817);
    
// polyethylene 
    G4Material* polyethylene = new G4Material("polyethylene", 0.94*g/cm3, 2, kStateSolid);
    polyethylene->AddElement(C, 1);
    polyethylene->AddElement(H, 2);
    
 
// paraffin
    G4Material* paraffin = new G4Material("paraffin", 0.9*g/cm3, 2, kStateSolid);
    paraffin->AddElement(C, 1);
    paraffin->AddElement(H, 2);
    
    
    //  quartz
    G4Material* quartz = new G4Material("quartz", 2.201*g/cm3, 2, kStateSolid);
    quartz->AddElement(Si, 1);
    quartz->AddElement(O, 2);
    
    
    // Stainless Steel
    G4Material *SS304LSteel = new G4Material("SS304LSteel", 8.00*g/cm3, 5, kStateSolid);
	SS304LSteel->AddElement(Fe, 0.65);
	SS304LSteel->AddElement(Cr, 0.20);
	SS304LSteel->AddElement(Ni, 0.12);
	SS304LSteel->AddElement(Mn, 0.02);
	SS304LSteel->AddElement(Si, 0.01);


    // False Stainless Steel
    G4Material *SS_Pole = new G4Material("SS_Pole", 0.80*g/cm3, 5, kStateSolid);
	SS_Pole->AddElement(Fe, 0.65);
	SS_Pole->AddElement(Cr, 0.20);
	SS_Pole->AddElement(Ni, 0.12);
	SS_Pole->AddElement(Mn, 0.02);
	SS_Pole->AddElement(Si, 0.01);
    
    
    
    // 316L Steel
    G4Material *SS316LSteel = new G4Material("SS316LSteel", 7000/12.7/51.6/5.5*g/cm3, 5, kStateSolid);      //-- Heat Exchanger 127W x 516H x 55D-----
	SS316LSteel->AddElement(Fe, 0.67);
	SS316LSteel->AddElement(Cr, 0.17);
	SS316LSteel->AddElement(Ni, 0.12);
	SS316LSteel->AddElement(Mn, 0.02);
	SS316LSteel->AddElement(Mo, 0.02);
    
    
    
// special Kovar (Co free material), the 3" PMT Shell Material, based on assumption of 8.35*g/cm3 and two components Fe and Ni only
// -- x+y=1--
// -- 7.874*x + 8.9*y = 8.35 --
    
    G4Material *Kovar = new G4Material("Kovar",8.35*g/cm3, 2, kStateSolid);
    Kovar->AddElement(Ni, 0.464);
    Kovar->AddElement(Fe, 0.536);
 
    
  //EJ301 (the liquid scintillator for neutron detector)
    G4Material* EJ301 = new G4Material("EJ301", density=0.874*g/cm3, 2, kStateLiquid);
    EJ301->AddElement(H, 0.548);
    EJ301->AddElement(C, 0.452);
    
  // Air
  G4Material* Air = new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);
  
  // Define vacuum
  G4Material* Vacuum = new G4Material("Vacuum", density= 1.e-5*g/cm3, ncomponents=1, kStateGas, STP_Temperature, 2.e-7*bar);
  Vacuum->AddMaterial(Air, fractionmass=1.);
  
  // Dump material information
  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}








void DetectorConstruction::SetupGeometry()
{
  
    
         G4Material* air = G4Material::GetMaterial("Air");
         G4Material* LAr = G4Material::GetMaterial("LiquidArgon");
         G4Material* GAr = G4Material::GetMaterial("GasArgon");
         G4Material*  Al = G4Material::GetMaterial("Aluminium");
      G4Material* quartz = G4Material::GetMaterial("quartz");
       G4Material* EJ301 = G4Material::GetMaterial("EJ301");
        G4Material* PTFE = G4Material::GetMaterial("Teflon");
       G4Material* steel = G4Material::GetMaterial("SS304LSteel");
  G4Material* pole_steel = G4Material::GetMaterial("SS_Pole");
      G4Material* vacuum = G4Material::GetMaterial("Vacuum");
G4Material* polyethylene = G4Material::GetMaterial("polyethylene");
    G4Material* paraffin = G4Material::GetMaterial("paraffin");
      G4Material* copper = G4Material::GetMaterial("Copper");
       G4Material* Kovar = G4Material::GetMaterial("Kovar");
      G4Material* SS316L = G4Material::GetMaterial("SS316LSteel");  
    G4Material* Tungsten = G4Material::GetMaterial("Tungsten");
    
    
    
    /*  //--- Cylindral world volume-----
    G4double world_Radius = 3*m;
    G4double world_halfHeight = 1.5*m;

 G4Tubs* worldSolid = new G4Tubs("World_Solid", 0*m, world_Radius, world_halfHeight, 0*deg, 360*deg);
     fpWorldLogical = new G4LogicalVolume(worldSolid, air, "World_Logical");
    fpWorldPhysical = new G4PVPlacement(0, G4ThreeVector(), fpWorldLogical, "World_Physical", 0, false, 0);
    */
    
    
    
    G4double world_halfLength = 3*m;
    
    G4Box* worldSolid = new G4Box("World_Solid", world_halfLength, world_halfLength, world_halfLength);
       fpWorldLogical = new G4LogicalVolume(worldSolid, vacuum, "World_Logical");
      fpWorldPhysical = new G4PVPlacement(0, G4ThreeVector(), fpWorldLogical, "World_Physical", 0, false, 0);

    
    
    //---------------------------------------------
    /*
     Geometry Description on LAr TPC and EJ301.
     
     1.  the z=0 plane is set to liquid/gas interface
     
     
     April 11th
     1. PE collimator is added in geometry.
     2. in particle gun setting, neutron is 50cm away from origion.
     3. particle gun offset in zAxis is 38cm. (total length of field cage is 76cm)
     4. PE collimator is set to 3cm away from neutron beam.
     
     
     May 30th
     1. The SS dewar and flange geometry are based on FNAL design.
     2. Coldhead and heat exchanger are also added.
     
     
     
     Nov. 2nd
     
     The Center of InnerSSDewar is the Origion (0, 0, 0)
     OuterSSDewar location is based on the SSHoldingBar length.
     
     Coldhead is moved 15*mm in xAxis direction to avoid overlap with heat exchanger
     
     
     Nov. 28th
     #1. change the distance btw LiF & TPC center to 58 cm
     #2. for 90-deg scattering case, change the distance btw TPC center and EJ301 active volume center to 58 cm
     #3. adjust the LiF right to the PE collimator (no distance in between)
     #4. From upper surface of the bottom quartz window to the bottom flange is 1 foot.
     */
    //---------------------------------------------
    
    
    
    
    // the Stainless Steel Color 
    //G4VisAttributes* SS_Color = new G4VisAttributes(G4Colour::Black());
    G4VisAttributes* SS_Color = new G4VisAttributes(G4Colour(0.2, 0.2, 0.2, 1.0));
    SS_Color->SetForceSolid(true);

    // the Aluminum Frame Color 
    //G4VisAttributes* Al_Color = new G4VisAttributes(G4Colour::Black());
    G4VisAttributes* Al_Color = new G4VisAttributes(G4Colour(0.4, 0.4, 0.4, 0.));
    Al_Color->SetForceSolid(true);
    
    
    // Color for EJs
    //G4VisAttributes* EJColor = new G4VisAttributes(G4Colour(0.1,0.8,0.1,0.1));
    //EJColor->SetForceSolid(false);
    //EJColor->SetForceWireframe(false);
    
    
    //----------- SS Dewar and flange --------------
    
    
    //--- InnerDewar ---
    const G4double      SSInnerDewar_OuterRadius = 4/2*25.4*mm;     //-- 4" OD---
    const G4double      SSInnerDewar_Thickness = 0.065*25.4*mm;     //-- 0.065"---
    const G4double      SSInnerDewar_TotHeight =    30*25.4*mm;     //-- 30" height---
    

    //---CF600 Flange---
    const G4double      SSInnerDewar_CF600Flange_OuterRadius = 5.97/2*25.4*mm; //-- 6" OD --
    const G4double      SSInnerDewar_CF600Flange_TotHeight = 0.78*25.4*mm;     //-- 0.78" thick --
    
    //--- OuterDewar ---
    const G4double      SSOuterDewar_OuterRadius =   14/2*25.4*mm;      //-- 14" OD --
    const G4double      SSOuterDewar_Thickness   =  0.188*25.4*mm;      //-- 0.188" thick --
    const G4double      SSOuterDewar_TotHeight   =   6*12*25.4*mm;      //-- 6' long --
    
    //--- CF1650 Flange ---
    const G4double      SSOuterDewar_CF1650Flange_OuterRadius = 16.5/2*25.4*mm;     //-- 16.5" OD --
    const G4double      SSOuterDewar_CF1650Flange_TotHeight   =   1.12*25.4*mm;
    
    
    //--- SS bar 1” X 0.75” X 29” ---
    const G4double      DewarHoldingBar_Height =   29*25.4*mm;
    const G4double      DewarHoldingBar_Width  =    1*25.4*mm;
    const G4double      DewarHoldingBar_Long   = 0.75*25.4*mm;
    
    
    
    //--- Relative Location parameters ---
    
    const G4double      SSInnerDewar_xCenterPosition = 0*mm;
    const G4double      SSInnerDewar_yCenterPosition = 0*mm;
    const G4double      SSInnerDewar_zCenterPosition = 0*mm;
    
    
    const G4double      SSOuterDewar_xCenterPosition = 0*mm;
    const G4double      SSOuterDewar_yCenterPosition = 0*mm;
    const G4double      SSOuterDewar_zCenterPosition = SSInnerDewar_TotHeight/2 + DewarHoldingBar_Height - SSOuterDewar_TotHeight/2;

    
    const G4double      HoldingBar_xPosition = SSInnerDewar_xCenterPosition + 0.*mm;
    const G4double      HoldingBar_yPosition = SSInnerDewar_yCenterPosition + SSInnerDewar_CF600Flange_OuterRadius + DewarHoldingBar_Long/2;
    const G4double      HoldingBar_zPosition = SSInnerDewar_zCenterPosition + SSInnerDewar_TotHeight/2 + DewarHoldingBar_Height/2;
    

    
    
    
    G4Tubs  *SSInnerDewar_CF600Flange = new G4Tubs("SSInnerDewar_CF600Flange",  0.*mm, SSInnerDewar_CF600Flange_OuterRadius,  SSInnerDewar_CF600Flange_TotHeight/2, 0.*deg, 360*deg);
    G4Tubs *SSOuterDewar_CF1650Flange = new G4Tubs("SSOuterDewar_CF1650Flange", 0.*mm, SSOuterDewar_CF1650Flange_OuterRadius, SSOuterDewar_CF1650Flange_TotHeight/2, 0.*deg, 360*deg);

    
    G4LogicalVolume  *SSInnerDewar_CF600Flange_Logical = new G4LogicalVolume(SSInnerDewar_CF600Flange, steel, "SSInnerDewar_CF600Flange_Logical");
    G4LogicalVolume *SSInnerDewar_CF1650Flange_Logical = new G4LogicalVolume(SSOuterDewar_CF1650Flange, steel, "SSInnerDewar_CF1650Flange_Logical");
    
    
    new G4PVPlacement(0, G4ThreeVector(SSInnerDewar_xCenterPosition, SSInnerDewar_yCenterPosition, (SSInnerDewar_zCenterPosition+SSInnerDewar_TotHeight/2+SSInnerDewar_CF600Flange_TotHeight/2)), SSInnerDewar_CF600Flange_Logical, "SSInnerDewar Top CF600Flange Physical", fpWorldLogical, true, 1, true);
    
    new G4PVPlacement(0, G4ThreeVector(SSInnerDewar_xCenterPosition, SSInnerDewar_yCenterPosition, (SSInnerDewar_zCenterPosition-SSInnerDewar_TotHeight/2-SSInnerDewar_CF600Flange_TotHeight/2)), SSInnerDewar_CF600Flange_Logical, "SSInnerDewar Bottom CF600Flange Physical", fpWorldLogical, true, 1, true);
    
    new G4PVPlacement(0, G4ThreeVector(SSOuterDewar_xCenterPosition, SSOuterDewar_yCenterPosition, (SSOuterDewar_zCenterPosition+SSOuterDewar_TotHeight/2+SSOuterDewar_CF1650Flange_TotHeight/2)), SSInnerDewar_CF1650Flange_Logical, "SSInnerDewar Top CF1650Flange Physical", fpWorldLogical, true, 1, true);

    new G4PVPlacement(0, G4ThreeVector(SSOuterDewar_xCenterPosition, SSOuterDewar_yCenterPosition, (SSOuterDewar_zCenterPosition-SSOuterDewar_TotHeight/2-SSOuterDewar_CF1650Flange_TotHeight/2)), SSInnerDewar_CF1650Flange_Logical, "SSInnerDewar Bottom CF1650Flange Physical", fpWorldLogical, true, 1, true);

  
    
    //SSInnerDewar_CF600Flange_Logical->SetVisAttributes(G4Colour::Grey());
   //SSInnerDewar_CF1650Flange_Logical->SetVisAttributes(G4Colour::Grey());
    
     SSInnerDewar_CF600Flange_Logical->SetVisAttributes(SS_Color);
    SSInnerDewar_CF1650Flange_Logical->SetVisAttributes(SS_Color);
    
    
                     
    G4Tubs   *SSInnerDewar_CF600FlangeRing = new G4Tubs("SSInnerDewar_CF600FlangeRing", SSInnerDewar_OuterRadius, SSInnerDewar_CF600Flange_OuterRadius, SSInnerDewar_CF600Flange_TotHeight/2, 0.*deg, 360*deg);
    G4Tubs  *SSOuterDewar_CF1650FlangeRing = new G4Tubs("SSOuterDewar_CF1650FlangeRing", SSOuterDewar_OuterRadius, SSOuterDewar_CF1650Flange_OuterRadius, SSOuterDewar_CF1650Flange_TotHeight/2, 0.*deg, 360*deg);

                      
    G4LogicalVolume  *SSInnerDewar_CF600FlangeRing_Logical = new G4LogicalVolume(SSInnerDewar_CF600FlangeRing, steel, "SSInnerDewar_CF600FlangeRing_Logical");
    G4LogicalVolume *SSOuterDewar_CF1650FlangeRing_Logical = new G4LogicalVolume(SSOuterDewar_CF1650FlangeRing, steel, "SSOuterDewar_CF1650FlangeRing_Logical");
       
    
    new G4PVPlacement(0, G4ThreeVector(SSInnerDewar_xCenterPosition, SSInnerDewar_yCenterPosition, (SSInnerDewar_zCenterPosition+SSInnerDewar_TotHeight/2-SSInnerDewar_CF600Flange_TotHeight/2)), SSInnerDewar_CF600FlangeRing_Logical, "SSInnerDewar Top CF600FlangeRing Physical", fpWorldLogical, true, 1, true);
                                        
    new G4PVPlacement(0, G4ThreeVector(SSInnerDewar_xCenterPosition, SSInnerDewar_yCenterPosition, (SSInnerDewar_zCenterPosition-SSInnerDewar_TotHeight/2+SSInnerDewar_CF600Flange_TotHeight/2)), SSInnerDewar_CF600FlangeRing_Logical, "SSInnerDewar Bottom CF600FlangeRing Physical", fpWorldLogical, true, 1, true);
                                   
    new G4PVPlacement(0, G4ThreeVector(SSOuterDewar_xCenterPosition, SSOuterDewar_yCenterPosition, (SSOuterDewar_zCenterPosition+SSOuterDewar_TotHeight/2-SSOuterDewar_CF1650Flange_TotHeight/2)), SSOuterDewar_CF1650FlangeRing_Logical, "SSOuterDewar Top CF1650FlangeRing Physical", fpWorldLogical, true, 1, true);
    
    new G4PVPlacement(0, G4ThreeVector(SSOuterDewar_xCenterPosition, SSOuterDewar_yCenterPosition, (SSOuterDewar_zCenterPosition-SSOuterDewar_TotHeight/2+SSOuterDewar_CF1650Flange_TotHeight/2)), SSOuterDewar_CF1650FlangeRing_Logical, "SSOuterDewar Bottom CF1650FlangeRing Physical", fpWorldLogical, true, 1, true);

    
     //SSInnerDewar_CF600FlangeRing_Logical->SetVisAttributes(G4Colour::Grey());
    //SSOuterDewar_CF1650FlangeRing_Logical->SetVisAttributes(G4Colour::Grey());
    
     SSInnerDewar_CF600FlangeRing_Logical->SetVisAttributes(SS_Color);
    SSOuterDewar_CF1650FlangeRing_Logical->SetVisAttributes(SS_Color);

    
    
    
    G4Tubs *SSInnerDewar = new G4Tubs("SSInnerDewar", SSInnerDewar_OuterRadius-SSInnerDewar_Thickness, SSInnerDewar_OuterRadius, SSInnerDewar_TotHeight/2, 0.*deg, 360*deg);
    G4Tubs *SSOuterDewar = new G4Tubs("SSOuterDewar", SSOuterDewar_OuterRadius-SSOuterDewar_Thickness, SSOuterDewar_OuterRadius, SSOuterDewar_TotHeight/2, 0.*deg, 360*deg);

                      
    G4LogicalVolume *SSInnerDewar_Logical = new G4LogicalVolume(SSInnerDewar, steel, "SSInnerDewar_Logical");
    G4LogicalVolume *SSOuterDewar_Logical = new G4LogicalVolume(SSOuterDewar, steel, "SSOuterDewar_Logical");
    
    
    new G4PVPlacement(0, G4ThreeVector(SSInnerDewar_xCenterPosition, SSInnerDewar_yCenterPosition, SSInnerDewar_zCenterPosition), SSInnerDewar_Logical, "SSInnerDewar Physical", fpWorldLogical, true, 1, true);
    new G4PVPlacement(0, G4ThreeVector(SSOuterDewar_xCenterPosition, SSOuterDewar_yCenterPosition, SSOuterDewar_zCenterPosition), SSOuterDewar_Logical, "SSOuterDewar Physical", fpWorldLogical, true, 1, true);
    
    
    G4Box *HoldingBar = new G4Box("HoldingBar", DewarHoldingBar_Width/2, DewarHoldingBar_Long/2, DewarHoldingBar_Height/2);
    
    G4LogicalVolume *HoldingBar_Logical = new G4LogicalVolume(HoldingBar, steel, "HoldingBar_Logical");
                 
    new G4PVPlacement(0, G4ThreeVector(HoldingBar_xPosition, HoldingBar_yPosition, HoldingBar_zPosition), HoldingBar_Logical, "SSInnerDewar Top HoldingBar Physical", fpWorldLogical, true, 1, true);
    new G4PVPlacement(0, G4ThreeVector(HoldingBar_xPosition, -HoldingBar_yPosition, HoldingBar_zPosition), HoldingBar_Logical, "SSInnerDewar Bottom HoldingBar Physical", fpWorldLogical, true, 1, true);
    
    
          //SSInnerDewar_Logical->SetVisAttributes(G4Colour::Grey());
          //SSOuterDewar_Logical->SetVisAttributes(G4Colour::Grey());
            //HoldingBar_Logical->SetVisAttributes(G4Colour::Grey());
    
    SSInnerDewar_Logical->SetVisAttributes(SS_Color);
    SSOuterDewar_Logical->SetVisAttributes(SS_Color);
      HoldingBar_Logical->SetVisAttributes(SS_Color);

    
   
    
    //------------------------heat exchanger and attachments----------------------------------------

    //---** heat exchanger geometry, back and right angle is set to position x=0, y=0 **---
    
    //-- bar dimensions 1.25” X 0.75” X 15”
    //-- plate dimensions 8.25” X 6” X 0.625”
    //-- The bottom of the plate is 3.25” from the bottom of the bar
    //-- Heat Exchanger 127W x 516H x 55D
    
    
    
    G4Box  *HeatExchanger_HoldingBar   = new G4Box("HeatExchanger_HoldingBar", 1.25*25.4*mm/2, 0.75*25.4*mm/2, 15*25.4*mm/2);
    G4Box  *HeatExchanger_HoldingPlate = new G4Box("HeatExchanger_HoldingPlate", 6*25.4*mm/2,  0.625*25.4*mm/2,8.25*25.4*mm/2);
    G4Box  *HeatExchanger              = new G4Box("HeatExchanger", 127*mm/2, 55*mm/2, 516*mm/2);

    
    G4LogicalVolume *HeatExchanger_HoldingBar_Logical   = new G4LogicalVolume(HeatExchanger_HoldingBar, steel, "HeatExchanger_HoldingBar_Logical");
    G4LogicalVolume *HeatExchanger_HoldingPlate_Logical = new G4LogicalVolume(HeatExchanger_HoldingPlate, steel, "HeatExchanger_HoldingPlate_Logical");
    G4LogicalVolume *HeatExchanger_Logical              = new G4LogicalVolume(HeatExchanger, SS316L, "HeatExchanger_Logical");

    
    new G4PVPlacement(0, G4ThreeVector(-127*mm/2-6*25.4*mm/2+1.25*25.4*mm/2, -0.625*25.4*mm-0.75*25.4*mm/2, SSOuterDewar_zCenterPosition+SSOuterDewar_TotHeight/2-15*25.4/2*mm), HeatExchanger_HoldingBar_Logical, "HeatExchanger HoldingBar Physical", fpWorldLogical, true, 0, true);
    
    new G4PVPlacement(0, G4ThreeVector(-127*mm/2, -0.625*25.4*mm/2, SSOuterDewar_zCenterPosition+SSOuterDewar_TotHeight/2-15*25.4*mm-3.25*25.4*mm+8.25*25.4*mm/2), HeatExchanger_HoldingPlate_Logical, "HeatExchanger HoldingPlate Physical", fpWorldLogical, true, 0, true);
    
    new G4PVPlacement(0, G4ThreeVector(-127*mm/2, 55*mm/2, SSOuterDewar_zCenterPosition+SSOuterDewar_TotHeight/2-15*25.4*mm-3.25*25.4*mm+8.25*25.4*mm/2), HeatExchanger_Logical, "HeatExchanger Physical", fpWorldLogical, true, 0, true);

    
      //HeatExchanger_HoldingBar_Logical->SetVisAttributes(G4Colour::Grey());
    //HeatExchanger_HoldingPlate_Logical->SetVisAttributes(G4Colour::Grey());
                 //HeatExchanger_Logical->SetVisAttributes(G4Colour::White());
    
  
    HeatExchanger_HoldingBar_Logical->SetVisAttributes(SS_Color);
    HeatExchanger_HoldingPlate_Logical->SetVisAttributes(SS_Color);
    HeatExchanger_Logical->SetVisAttributes(SS_Color);

    
    
    
       
    
    
    
    
    //---------------------Coldhead, PT60 from Cryomech-------------------------------------------
    
    
    const G4double  Coldehead_xCenterPosition = -SSInnerDewar_xCenterPosition + SSInnerDewar_OuterRadius+15*mm;
    const G4double  Coldehead_yCenterPosition = SSInnerDewar_yCenterPosition;
    const G4double  Coldehead_zCenterPosition = SSOuterDewar_zCenterPosition + SSOuterDewar_TotHeight/2 + SSOuterDewar_CF1650Flange_TotHeight;
    
    const G4double  Coldehead_ColumnHeight = 176.1*mm-7.7*mm-SSOuterDewar_CF1650Flange_TotHeight;
    const G4double  Coldehead_FlangeRadius = 124.5*mm/2;
    
    G4Tubs *Coldhead_Flange  = new G4Tubs("Coldhead_Flange", 0.*mm, Coldehead_FlangeRadius,   12.7*mm/2, 0.*deg, 360*deg);
    G4Tubs *Coldhead_Column  = new G4Tubs("Coldhead_Column", 0.*mm, 70.6*mm/4, Coldehead_ColumnHeight/2, 0.*deg, 360*deg);
    
    G4Box  *Coldehead_TopBox = new G4Box("Coldehead_TopBox", Coldehead_FlangeRadius*1.4142/2, Coldehead_FlangeRadius*1.4142/2, 205.6*mm/2);
    
    
    G4LogicalVolume *Coldhead_Flange_Logical  = new G4LogicalVolume(Coldhead_Flange,  steel, "Coldhead_Flange_Logical");
    G4LogicalVolume *Coldehead_TopBox_Logical = new G4LogicalVolume(Coldehead_TopBox, steel, "Coldehead_TopBox_Logical");
    G4LogicalVolume *Coldhead_Column_Logical  = new G4LogicalVolume(Coldhead_Column,  copper, "Coldhead_Column_Logical");

    
    new G4PVPlacement(0, G4ThreeVector(Coldehead_xCenterPosition, Coldehead_yCenterPosition, Coldehead_zCenterPosition+12.7*mm/2), Coldhead_Flange_Logical, "Coldhead SSFlange Physical", fpWorldLogical, true, 0, true);

    new G4PVPlacement(0, G4ThreeVector(Coldehead_xCenterPosition, Coldehead_yCenterPosition, Coldehead_zCenterPosition+12.7*mm+205.6*mm/2), Coldehead_TopBox_Logical, "Coldehead TopBox Physical", fpWorldLogical, true, 0, true);
    
    new G4PVPlacement(0, G4ThreeVector(Coldehead_xCenterPosition - 70.6*mm/4, Coldehead_yCenterPosition, SSOuterDewar_zCenterPosition + SSOuterDewar_TotHeight/2 - Coldehead_ColumnHeight/2), Coldhead_Column_Logical, "Coldhead Column Physical 1", fpWorldLogical, true, 1, true);
    
    new G4PVPlacement(0, G4ThreeVector(Coldehead_xCenterPosition + 70.6*mm/4, Coldehead_yCenterPosition, SSOuterDewar_zCenterPosition + SSOuterDewar_TotHeight/2 - Coldehead_ColumnHeight/2), Coldhead_Column_Logical, "Coldhead Column Physical 2", fpWorldLogical, true, 1, true);
    
    
    //Coldhead_Flange_Logical->SetVisAttributes(G4Colour::Grey());
    //Coldhead_Column_Logical->SetVisAttributes(G4Colour::Grey());
    //Coldehead_TopBox_Logical->SetVisAttributes(G4Colour::White());
    
    Coldhead_Flange_Logical->SetVisAttributes(SS_Color);
    Coldhead_Column_Logical->SetVisAttributes(SS_Color);
    Coldehead_TopBox_Logical->SetVisAttributes(SS_Color);
    
    

    const G4int     Coldhead_zPlaneNbr = 7;
    const G4double  Coldehead_zPlane[Coldhead_zPlaneNbr] = {0.*mm,   -0.31*25.4*mm,  -0.31*25.4*mm,   (-0.31-0.41)*25.4*mm,  (-0.31-0.41)*25.4*mm,  (-0.31-0.41-0.65)*25.4*mm, (-0.31-1.76)*25.4*mm};
    const G4double  Coldehead_rInner[Coldhead_zPlaneNbr] = {0.*mm,     0.*mm,      0.*mm,       0.*mm,        0.*mm,      0.*mm,        0.*mm};
    const G4double  Coldehead_rOuter[Coldhead_zPlaneNbr] = {2.78/2*25.4*mm,  2.78/2*25.4*mm,   4.348/2*25.4*mm,    4.348/2*25.4*mm,             4/2*25.4*mm,              4/2*25.4*mm,  0.*mm};
    

    G4Polycone              *Coldhead = new G4Polycone("Coldhead", 0.*deg, 360.0*deg, Coldhead_zPlaneNbr, Coldehead_zPlane, Coldehead_rInner, Coldehead_rOuter);
    
    G4LogicalVolume* Coldhead_Logical = new G4LogicalVolume(Coldhead ,	 copper,  "Coldhead_Logical");
    
    new G4PVPlacement(0,  G4ThreeVector(Coldehead_xCenterPosition, Coldehead_yCenterPosition, SSOuterDewar_zCenterPosition + SSOuterDewar_TotHeight/2 - Coldehead_ColumnHeight),  Coldhead_Logical,   "Coldhead Physical",  fpWorldLogical,  true,  1, true);	 
  
    
    
    const G4int     Coldhead_Cup_zPlaneNbr = 4;
    const G4double  Coldhead_Cup_zPlane[Coldhead_Cup_zPlaneNbr] = {0.*mm,                -4.7*25.4*mm,            -4.7*25.4*mm,         (-4.7-0.08)*25.4*mm};
    const G4double  Coldhead_Cup_rInner[Coldhead_Cup_zPlaneNbr] = {2*25.4*mm,             2*25.4*mm,               0.*mm,                0.*mm};
    const G4double  Coldhead_Cup_rOuter[Coldhead_Cup_zPlaneNbr] = {(2+0.08)*25.4*mm,     (2+0.08)*25.4*mm,         (2+0.08)*25.4*mm,     (2+0.08)*25.4*mm};
    
    
    G4Polycone              *Coldhead_Cup = new G4Polycone("Coldhead_Cup", 0.*deg, 360.0*deg, Coldhead_Cup_zPlaneNbr, Coldhead_Cup_zPlane, Coldhead_Cup_rInner, Coldhead_Cup_rOuter);
    
    G4LogicalVolume* Coldhead_Cup_Logical = new G4LogicalVolume(Coldhead_Cup ,	 steel,  "Coldhead_Cup_Logical");   
    
    new G4PVPlacement(0,  G4ThreeVector(Coldehead_xCenterPosition, Coldehead_yCenterPosition, SSOuterDewar_zCenterPosition + SSOuterDewar_TotHeight/2 - Coldehead_ColumnHeight-0.72*25.4*mm),  Coldhead_Cup_Logical,   "Coldhead Cup Physical",  fpWorldLogical,  true,  1, true);
    
    
    
       Coldhead_Cup_Logical->SetVisAttributes(G4Colour::Yellow());
           Coldhead_Logical->SetVisAttributes(G4Colour::Yellow());

    
    
    
    
    
    
    //-----------------------------------------------------------------------------------------
    
    //---  Oct. 11th
    //---  SCENE prototype TPC, the orign is center/middle of Inner SS Dewar. -----
    //---  TPC_zAxis_Offset adjusts the location of TPC with respect to the orign   -----
    //---  TPC zAixis Center is the interface of Liquid/Gas ---
    
    //---  New G4 geometry function is use.
    
    //---- From upper surface of the bottom quartz window to the bottom flange is 1 foot. ---
    
    
    
    const G4double TPC_zAxis_Offset    = 0.*25.4*mm;  // equals to SSInnerDewar_TotHeight/2+12*25.4*mm+3.*25.4*mm;
    
    const G4double LAr_Volume_Height   = SSInnerDewar_TotHeight/2-TPC_zAxis_Offset;
    const G4double LAr_Volume_Radius   = SSInnerDewar_OuterRadius-SSInnerDewar_Thickness;
    
    const G4double GAr_Volume_Height   = SSInnerDewar_TotHeight/2+TPC_zAxis_Offset;
    const G4double GAr_Volume_Radius   = SSInnerDewar_OuterRadius-SSInnerDewar_Thickness;
    
    
    //const G4double LArSubtraction_Volume_Height   = 3.*25.4*mm;
    //const G4double LArSubtraction_Volume_Diameter = 2.7*25.4*mm;
    
    
    
    
    
    
    //--------------TPC Dimensions--------------
    
    const G4double ActiveGAr_Volume_Height     = 0.25*25.4*mm;
    const G4double ActiveGAr_Volume_Diameter   = 2.7*25.4*mm;
    
    const G4double ActiveLAr_Volume_Height     = 3.*25.4*mm;
    const G4double ActiveLAr_Volume_Diameter   = 2.7*25.4*mm;
     
    const G4double QuartzWindow_Volume_Height        = 0.6*25.4*mm;
    const G4double QuartzWindow_Volume_Diameter      = 3.*25.4*mm;

    const G4double PTFERing_Volume_Height            = 2*25.4*mm;
    const G4double PTFERing_Volume_InnerDiameter     = 3*25.4*mm;
    const G4double PTFERing_Volume_OuterDiameter     = 3.8*25.4*mm;

    const G4double PTFEColumn_Volume_xHalfLength     = 0.5/2*25.4*mm;
    const G4double PTFEColumn_Volume_yHalfLength     = 0.5/2*25.4*mm;
    
    const G4double PTFEColumn_LArVolume_zHalfLength  = ActiveLAr_Volume_Height/2;
    const G4double PTFEColumn_GArVolume_zHalfLength  = ActiveGAr_Volume_Height/2;
    
     
    //---------------3" PMT Solid Part ------------
    
    const G4double  TopPMT_zOffset  =  -TPC_zAxis_Offset+ActiveGAr_Volume_Height+QuartzWindow_Volume_Height;
    const G4double  BtmPMT_zOffset  =  -TPC_zAxis_Offset-ActiveLAr_Volume_Height-QuartzWindow_Volume_Height;
    
    
    const G4double  PMT_TotLength   = 123.063*mm;
    
    const G4int     PMTzPlaneNbr = 8;
    const G4double  PMTShellThickness = 0.5*mm;
    //const G4double  PMTzPlane[PMTzPlaneNbr] = {0.*mm+PMT_zOffset, -17.4115*mm+PMT_zOffset, -35.041*mm+PMT_zOffset, -39.9343*mm+PMT_zOffset, -116.*mm+PMT_zOffset, -116.078*mm+PMT_zOffset, -116.713*mm+PMT_zOffset, -PMT_TotLength+PMT_zOffset};
    const G4double  PMTzPlane[PMTzPlaneNbr] = {0.*mm, -17.4115*mm, -35.041*mm, -39.9343*mm, -116.*mm, -116.078*mm, -116.713*mm, -PMT_TotLength};
    
    const G4double  PMTrInner[PMTzPlaneNbr] = {36.576*mm, 36.576*mm, 28.0201*mm, 25.146*mm, 25.146*mm, 0.*mm, 0.*mm, 25.126*mm};
    const G4double  PMTrOuter[PMTzPlaneNbr] = {36.576*mm+PMTShellThickness, 36.576*mm+PMTShellThickness, 28.0201*mm+PMTShellThickness, 25.146*mm+PMTShellThickness, 25.146*mm+PMTShellThickness, 25.146*mm+PMTShellThickness, 25.146*mm+PMTShellThickness, 25.146*mm+PMTShellThickness};
    
    const G4double  PMT_QuartzWindow_Radius = 35.576*mm;
    const G4double  PMT_QuartzWindow_TotLength = 3.4*mm;
    
      //-----------------------------------------------

    
    const G4double LArSubtraction_Volume_Height     = PMT_TotLength+ActiveLAr_Volume_Height+QuartzWindow_Volume_Height;
    const G4double GArSubtraction_Volume_Height     = PMT_TotLength+ActiveGAr_Volume_Height+QuartzWindow_Volume_Height;
    
    const G4double LArSubtraction_Volume_Diameter   = QuartzWindow_Volume_Diameter;
    const G4double GArSubtraction_Volume_Diameter   = QuartzWindow_Volume_Diameter;
    
    
    //--------------------------------------------------------------------------------------------
    
    
       
    //---------------Dead LAr and GAr ------------
    
    G4Tubs             *LAr_Volume  = new G4Tubs("LAr_Volume",       0.*mm,  LAr_Volume_Radius,   LAr_Volume_Height/2, 0.*deg, 360*deg);
    G4Tubs             *GAr_Volume  = new G4Tubs("GAr_Volume",       0.*mm,  GAr_Volume_Radius,   GAr_Volume_Height/2, 0.*deg, 360*deg);
    
    
    G4Tubs  *LArSubtraction_Volume  = new G4Tubs("LArSubtraction_Volume", 0.*mm, LArSubtraction_Volume_Diameter/2,   LArSubtraction_Volume_Height/2, 0.*deg, 360*deg);
    G4Tubs  *GArSubtraction_Volume  = new G4Tubs("GArSubtraction_Volume", 0.*mm, GArSubtraction_Volume_Diameter/2,   GArSubtraction_Volume_Height/2, 0.*deg, 360*deg);

    
     //---------LAr PTFE Part ------------    
    
    
    G4Tubs   *PTFERing_Volume  = new G4Tubs("PTFERing_Volume",    PTFERing_Volume_InnerDiameter/2,  PTFERing_Volume_OuterDiameter/2,   PTFERing_Volume_Height/2, 0.*deg, 360*deg);
    
    G4Box  *PTFEColumn_BtmVolume  = new G4Box("PTFEColumn_BtmVolume", PTFEColumn_Volume_xHalfLength,  PTFEColumn_Volume_yHalfLength,  PTFEColumn_LArVolume_zHalfLength);

    
    G4RotationMatrix  Rotat_PTFEColumn;
    G4ThreeVector    Vector_PTFEColumn;
    
    
    Vector_PTFEColumn.setX(PTFERing_Volume_OuterDiameter/2-PTFEColumn_Volume_yHalfLength); 
    Vector_PTFEColumn.setY(0.*mm); 
    Vector_PTFEColumn.setZ(PTFEColumn_LArVolume_zHalfLength+PTFERing_Volume_Height/2);
    
    
    G4VSolid*  PTFE_Leg1_Volume = new G4UnionSolid("PTFE_Leg1_Volume", PTFERing_Volume, PTFEColumn_BtmVolume,  G4Transform3D(Rotat_PTFEColumn, Vector_PTFEColumn)); 
    
    
    Rotat_PTFEColumn.rotateZ(120.0*deg);
    
    Vector_PTFEColumn.setX((PTFERing_Volume_OuterDiameter/2-PTFEColumn_Volume_yHalfLength)*std::cos(120./180.*3.1415926)); 
    Vector_PTFEColumn.setY((PTFERing_Volume_OuterDiameter/2-PTFEColumn_Volume_yHalfLength)*std::sin(120./180.*3.1415926)); 
    Vector_PTFEColumn.setZ(PTFEColumn_LArVolume_zHalfLength+PTFERing_Volume_Height/2);
    
    
    G4VSolid*  PTFE_Leg2_Volume = new G4UnionSolid("PTFE_Leg2_Volume", PTFE_Leg1_Volume, PTFEColumn_BtmVolume,  G4Transform3D(Rotat_PTFEColumn, Vector_PTFEColumn)); 
    
    
    Rotat_PTFEColumn.rotateZ(120.0*deg);
    
    Vector_PTFEColumn.setX((PTFERing_Volume_OuterDiameter/2-PTFEColumn_Volume_yHalfLength)*std::cos(240./180.*3.1415926)); 
    Vector_PTFEColumn.setY((PTFERing_Volume_OuterDiameter/2-PTFEColumn_Volume_yHalfLength)*std::sin(240./180.*3.1415926)); 
    Vector_PTFEColumn.setZ(PTFEColumn_LArVolume_zHalfLength+PTFERing_Volume_Height/2);
    
    
    G4VSolid*  PTFE_BtmVolume = new G4UnionSolid("PTFE_BtmVolume", PTFE_Leg2_Volume, PTFEColumn_BtmVolume,  G4Transform3D(Rotat_PTFEColumn, Vector_PTFEColumn)); 
    

    
   //--------------GAr PTFE part------------------- 
    
    
    G4Box  *PTFEColumn_TopVolume  = new G4Box("PTFEColumn_TopVolume", PTFEColumn_Volume_xHalfLength,  PTFEColumn_Volume_yHalfLength,  PTFEColumn_GArVolume_zHalfLength);
    
    Rotat_PTFEColumn.rotateZ(120.0*deg);
    
    Vector_PTFEColumn.setX(PTFERing_Volume_OuterDiameter/2-PTFEColumn_Volume_yHalfLength); 
    Vector_PTFEColumn.setY(0.*mm); 
    Vector_PTFEColumn.setZ(-PTFEColumn_GArVolume_zHalfLength-PTFERing_Volume_Height/2);
    
    
    G4VSolid*  PTFE_Leg1_TopVolume = new G4UnionSolid("PTFE_Leg1_TopVolume", PTFERing_Volume, PTFEColumn_TopVolume,  G4Transform3D(Rotat_PTFEColumn, Vector_PTFEColumn)); 
    
    
    Rotat_PTFEColumn.rotateZ(120.0*deg);
    
    Vector_PTFEColumn.setX((PTFERing_Volume_OuterDiameter/2-PTFEColumn_Volume_yHalfLength)*std::cos(120./180.*3.1415926)); 
    Vector_PTFEColumn.setY((PTFERing_Volume_OuterDiameter/2-PTFEColumn_Volume_yHalfLength)*std::sin(120./180.*3.1415926)); 
    Vector_PTFEColumn.setZ(-PTFEColumn_GArVolume_zHalfLength-PTFERing_Volume_Height/2);    
    
    G4VSolid*  PTFE_Leg2_TopVolume = new G4UnionSolid("PTFE_Leg2_TopVolume", PTFE_Leg1_TopVolume, PTFEColumn_TopVolume,  G4Transform3D(Rotat_PTFEColumn, Vector_PTFEColumn)); 
    
    
    Rotat_PTFEColumn.rotateZ(120.0*deg);
    
    Vector_PTFEColumn.setX((PTFERing_Volume_OuterDiameter/2-PTFEColumn_Volume_yHalfLength)*std::cos(240./180.*3.1415926)); 
    Vector_PTFEColumn.setY((PTFERing_Volume_OuterDiameter/2-PTFEColumn_Volume_yHalfLength)*std::sin(240./180.*3.1415926)); 
    Vector_PTFEColumn.setZ(-PTFEColumn_GArVolume_zHalfLength-PTFERing_Volume_Height/2);
    
    
    G4VSolid*  PTFE_TopVolume = new G4UnionSolid("PTFE_Volume", PTFE_Leg2_TopVolume, PTFEColumn_TopVolume,  G4Transform3D(Rotat_PTFEColumn, Vector_PTFEColumn)); 
    
    
    //---------------------------------------------
    
    
    G4VSolid*  LArSubtraction_Part= new G4UnionSolid("LArSubtraction_Part", LArSubtraction_Volume, PTFE_BtmVolume, 0, G4ThreeVector(0., 0., LArSubtraction_Volume_Height/2-PTFERing_Volume_Height/2-PTFEColumn_LArVolume_zHalfLength*2)); 
    
    G4VSolid*  GArSubtraction_Part= new G4UnionSolid("GArSubtraction_Part", GArSubtraction_Volume, PTFE_TopVolume, 0, G4ThreeVector(0., 0., -GArSubtraction_Volume_Height/2+PTFERing_Volume_Height/2+PTFEColumn_GArVolume_zHalfLength*2)); 

    
    //G4LogicalVolume *Test_Logical = new G4LogicalVolume(LArSubtraction_Part, quartz, "Test_Logical");
    
    
    //new G4PVPlacement(0, G4ThreeVector(0., 0., 0.),  Test_Logical, "Test Volume Physical", fpWorldLogical, true, 1, true);
    
    
    G4SubtractionSolid* NewLAr_Volume = new G4SubtractionSolid("NewLAr_Volume", LAr_Volume, LArSubtraction_Part, 0, G4ThreeVector(0.,0.,(LAr_Volume_Height-LArSubtraction_Volume_Height)/2));
    
    G4SubtractionSolid* NewGAr_Volume = new G4SubtractionSolid("NewGAr_Volume", GAr_Volume, GArSubtraction_Part, 0, G4ThreeVector(0.,0.,(-GAr_Volume_Height+GArSubtraction_Volume_Height)/2));

        
    G4LogicalVolume *LAr_Volume_Logical = new G4LogicalVolume(NewLAr_Volume, LAr, "LAr_Volume_Logical");
    G4LogicalVolume *GAr_Volume_Logical = new G4LogicalVolume(NewGAr_Volume, GAr, "GAr_Volume_Logical");

    	LAr_Volume_Logical->SetVisAttributes(G4VisAttributes::Invisible);
        GAr_Volume_Logical->SetVisAttributes(G4VisAttributes::Invisible);
    
    
    new G4PVPlacement(0, G4ThreeVector(0., 0., -LAr_Volume_Height/2-TPC_zAxis_Offset),  LAr_Volume_Logical, "LAr Volume Physical", fpWorldLogical, false, 1, true);
    new G4PVPlacement(0, G4ThreeVector(0., 0.,  GAr_Volume_Height/2-TPC_zAxis_Offset),  GAr_Volume_Logical, "GAr Volume Physical", fpWorldLogical, false, 1, true);
    
    
    G4LogicalVolume *PTFE_BtmVolume_Logical = new G4LogicalVolume(PTFE_BtmVolume, PTFE, "PTFE_BtmVolume_Logical");
    G4LogicalVolume *PTFE_TopVolume_Logical = new G4LogicalVolume(PTFE_TopVolume, PTFE, "GAr_Volume_Logical");
    
    
    //PTFE_BtmVolume_Logical->SetVisAttributes(G4Colour::Grey());
    //PTFE_TopVolume_Logical->SetVisAttributes(G4Colour::Grey());
    
    PTFE_BtmVolume_Logical->SetVisAttributes(G4Colour::Magenta());
    PTFE_TopVolume_Logical->SetVisAttributes(G4Colour::Magenta());
    
    
    new G4PVPlacement(0, G4ThreeVector(0., 0., -PTFERing_Volume_Height/2-TPC_zAxis_Offset-PTFEColumn_LArVolume_zHalfLength*2),  PTFE_BtmVolume_Logical, "Bottom  PTFE Volume Physical", fpWorldLogical, true, 1, true);
    
    new G4PVPlacement(0, G4ThreeVector(0., 0.,  PTFERing_Volume_Height/2+PTFEColumn_GArVolume_zHalfLength*2-TPC_zAxis_Offset),  PTFE_TopVolume_Logical, "GAr Volume Physical", fpWorldLogical, true, 1, true);
    
    
        //---------------Active LAr and GAr ------------
    
      G4Tubs  *ActiveLAr_Volume  = new G4Tubs("ActiveLAr_Volume",    0.*mm,  ActiveLAr_Volume_Diameter/2,   ActiveLAr_Volume_Height/2, 0.*deg, 360*deg);
      G4Tubs  *ActiveGAr_Volume  = new G4Tubs("ActiveGAr_Volume",    0.*mm,  ActiveGAr_Volume_Diameter/2,   ActiveGAr_Volume_Height/2, 0.*deg, 360*deg);
    
        G4LogicalVolume *ActiveLAr_Volume_Logical = new G4LogicalVolume(ActiveLAr_Volume, LAr, "ActiveLAr_Volume_Logical");
        G4LogicalVolume *ActiveGAr_Volume_Logical = new G4LogicalVolume(ActiveGAr_Volume, GAr, "ActiveGAr_Volume_Logical");
    
    ActiveGAr_Volume_Logical->SetVisAttributes(G4Colour::Green());
    ActiveLAr_Volume_Logical->SetVisAttributes(G4Colour::Blue());

    new G4PVPlacement(0, G4ThreeVector(0., 0., -ActiveLAr_Volume_Height/2-TPC_zAxis_Offset),  ActiveLAr_Volume_Logical, "ActiveLAr Volume Physical", fpWorldLogical, true, 1, true);
    new G4PVPlacement(0, G4ThreeVector(0., 0.,  ActiveGAr_Volume_Height/2-TPC_zAxis_Offset),  ActiveGAr_Volume_Logical, "ActiveGAr Volume Physical", fpWorldLogical, true, 1, true);
    
    
        //---------------Quartz Window ------------------
    
       G4Tubs  *QuartzWindow_Volume  = new G4Tubs("QuartzWindow_Volume",    0.*mm,  QuartzWindow_Volume_Diameter/2,   QuartzWindow_Volume_Height/2, 0.*deg, 360*deg);
    
       G4LogicalVolume *QuartzWindow_Volume_Logical = new G4LogicalVolume(QuartzWindow_Volume, quartz, "QuartzWindow_Volume_Logical");
        
    QuartzWindow_Volume_Logical->SetVisAttributes(G4Colour::Red());
    
    G4RotationMatrix  Rotat_QuartzWindow;
    G4ThreeVector    Vector_QuartzWindow;
    
     Rotat_QuartzWindow.rotateY(180.0*deg);
    
    Vector_QuartzWindow.setX(0.*mm); 
    Vector_QuartzWindow.setY(0.*mm); 
    Vector_QuartzWindow.setZ(QuartzWindow_Volume_Height/2-TPC_zAxis_Offset+ActiveGAr_Volume_Height);

    
    new G4PVPlacement(G4Transform3D(Rotat_QuartzWindow, Vector_QuartzWindow),  QuartzWindow_Volume_Logical, "Top Quartz Window Physical",  fpWorldLogical, true, 1,  true);
    
    new G4PVPlacement(0, G4ThreeVector(0,0,-ActiveLAr_Volume_Height-TPC_zAxis_Offset-QuartzWindow_Volume_Height/2),   QuartzWindow_Volume_Logical, "Bottom Quartz Window Physical",  fpWorldLogical,  true,  1,  true);	
    
   

    
    //-----------------3" PMT Solid Part ---------------
    
    G4Polycone              *PMT_Shell = new G4Polycone("PMT_Shell", 0.*deg, 360.0*deg, PMTzPlaneNbr, PMTzPlane, PMTrInner, PMTrOuter);
    
    G4LogicalVolume* PMT_Shell_Logical = new G4LogicalVolume(PMT_Shell,	 Kovar,  "PMT_Shell_Logical");   
    
    
    G4RotationMatrix  Rotat_TopPMT;
    G4ThreeVector    Vector_TopPMT;
    
    Rotat_TopPMT.rotateY(180.0*deg);
    Vector_TopPMT.setX(0.*mm); Vector_TopPMT.setY(0.*mm); Vector_TopPMT.setZ(TopPMT_zOffset);
    
    
    new G4PVPlacement(G4Transform3D(Rotat_TopPMT, Vector_TopPMT),  PMT_Shell_Logical, "Top PMT Shell Physical",  fpWorldLogical, true, 1,  true);
    new G4PVPlacement(0, G4ThreeVector(0,0,BtmPMT_zOffset),        PMT_Shell_Logical, "Bottom PMT Shell Physical",  fpWorldLogical,  true,  1,  true);	                           
    
    
    G4Tubs  *PMT_QuartzWindow = new G4Tubs("PMT_QuartzWindow", 0.*mm, PMT_QuartzWindow_Radius, PMT_QuartzWindow_TotLength/2, 0.*deg, 360*deg);
    
    G4LogicalVolume *PMT_QuartzWindow_Logical = new G4LogicalVolume(PMT_QuartzWindow, quartz, "PMT_QuartzWindow_Logical");
    
    PMT_QuartzWindow_Logical->SetVisAttributes(G4Colour::Red());
           //PMT_Shell_Logical->SetVisAttributes(G4Colour::Yellow());
           PMT_Shell_Logical->SetVisAttributes(G4Colour::Cyan());
    
    new G4PVPlacement(0, G4ThreeVector(0,0, TopPMT_zOffset+PMT_QuartzWindow_TotLength/2), PMT_QuartzWindow_Logical, "Top PMT QuartzWindow Physical", fpWorldLogical, true, 1,  true);
    
    new G4PVPlacement(0, G4ThreeVector(0,0, BtmPMT_zOffset-PMT_QuartzWindow_TotLength/2), PMT_QuartzWindow_Logical, "Bottom PMT QuartzWindow Physical", fpWorldLogical, true, 1,  true);
    
    
    //----------------------------------------------- 


    

    
    
    
    
    //--------------------------------------------------------------------------------------------   
     
    
     //--- for 2"DIA x 2" Chamber of Scintillator, neutron detector ---
    
     const G4double                EJ_Radius = 2.5*2.54*cm;
     const G4double            EJ_halfHeight = 2.5*2.54*cm;
    const G4double   TPC_Center_zAxis_Offset = -1.5*2.54*cm; //----- TPC center is the LAr center -------
    
    
    
    /*
     
     const double            TopEJ_Angle = 24.72/180*TMath::Pi();
     
     
     double LiF_PositionX = -64.77;
     double LiF_PositionY = 0.;
     double LiF_PositionZ = 31.36;
     
     
     double TopEJ301Center_PositionX = 71.12*TMath::Cos(TopEJ_Angle);
     double TopEJ301Center_PositionY = 0.;
     double TopEJ301Center_PositionZ = 71.12*TMath::Sin(TopEJ_Angle);
     
     double SideEJ301Center_PositionX =  40.663;
     double SideEJ301Center_PositionY =  54.9253;   //--- (+/-) ---
     double SideEJ301Center_PositionZ = -19.6923;
     
     */
    
    
    
    G4Tubs*                         EJ_tube = new G4Tubs("EJ_tube", 0.*cm, EJ_Radius, EJ_halfHeight, 0.*deg, 360.*deg);
    
    G4LogicalVolume*        EJ_tube_Logical = new G4LogicalVolume(EJ_tube, EJ301, "EJ_tube_Logical");
    
    
    EJ_tube_Logical->SetVisAttributes(G4Colour::Red());

    
    G4RotationMatrix R_EJ1, R_EJ2, R_EJ3;
    G4ThreeVector    T_EJ1, T_EJ2, T_EJ3;
    
    
    const double            PiValue = 3.1415926;
    
    
    const G4double       zAxis_TPC_Offset = -1.5*2.54*cm;

    // --- only  EJ#4 exists in this 3D geometry (uses name EJ#4 )---
    // --- this is specifically used for the simulation of EJ#4 at ~90 deg ( 50keVnr recoil ) ---
    // --- the specific postion is based on the measurement from elog#4111 ---
    
    const G4double       EJ4_xAxis_Pos = 14.7*cm;
    const G4double       EJ4_yAxis_Pos = -77.3*cm;
    const G4double       EJ4_zAxis_Pos =  5.9*cm + TPC_Center_zAxis_Offset;

    
    new G4PVPlacement(0, G4ThreeVector(EJ4_xAxis_Pos, EJ4_yAxis_Pos, EJ4_zAxis_Pos), EJ_tube_Logical, "EJ301_4", fpWorldLogical, true, 1,  true);

    //----------------------------------------------------------------------------------------------------------------------
   
    
    
    
    
    
   
    
    //--------------------------- SCENE PE & Paraffin Shielding Dimensions -------------------------------------
    
    
    const G4double           TPC_Center_Ground_Distance = 142.91*cm;
    const G4double   zAxis_TPC_Center_BeamLine_Distance = 31.4*cm;

    const double    EJ2_Shielding_Angle = 20;
    
    //----- For EJ#1 & EJ#2, PE Shielding  Dimension Dia. 8.75" x 8.75" --------
    
    const G4double     PE_Polythene_Shielding_Radius       = 8.75*2.54/2*cm;
    const G4double     PE_Polythene_Shielding_Half_Length  = 8.75*2.54/2*cm;
    
    
    const G4double    EJ1_Shielding_Btm_Ground_Distance = 166*cm;
    const G4double    EJ2_Shielding_Btm_Ground_Distance = 141*cm;

    
    const G4double     EJ1_Polythene_Shielding_xAxisPos = -SSOuterDewar_OuterRadius-PE_Polythene_Shielding_Radius;
    const G4double     EJ1_Polythene_Shielding_zAxisPos =  EJ1_Shielding_Btm_Ground_Distance-TPC_Center_Ground_Distance+1.5*2.54*cm+PE_Polythene_Shielding_Half_Length;

    const G4double     EJ2_Polythene_Shielding_xAxisPos = -(PE_Polythene_Shielding_Radius+SSOuterDewar_OuterRadius)*std::sin(EJ2_Shielding_Angle/180.*PiValue);
    const G4double     EJ2_Polythene_Shielding_yAxisPos = -(PE_Polythene_Shielding_Radius+SSOuterDewar_OuterRadius)*std::cos(EJ2_Shielding_Angle/180.*PiValue);
    const G4double     EJ2_Polythene_Shielding_zAxisPos =  EJ2_Shielding_Btm_Ground_Distance-TPC_Center_Ground_Distance+1.5*2.54*cm+PE_Polythene_Shielding_Half_Length;

    
    
    G4Tubs* PE_Shielding_tube = new G4Tubs("PE_Shielding_tube", 0.*cm,  PE_Polythene_Shielding_Radius,  PE_Polythene_Shielding_Half_Length,  0.*deg,  360.*deg);
    
    G4LogicalVolume*  PE_Shielding_Logical = new G4LogicalVolume(PE_Shielding_tube, polyethylene, "PE_Shielding_Logical");
    
    
    new G4PVPlacement(0, G4ThreeVector(EJ1_Polythene_Shielding_xAxisPos, 0, EJ1_Polythene_Shielding_zAxisPos), PE_Shielding_Logical, "EJ1_Shielding", fpWorldLogical, true, 1,  true);
    
    
    new G4PVPlacement(0, G4ThreeVector(EJ2_Polythene_Shielding_xAxisPos, EJ2_Polythene_Shielding_yAxisPos, EJ2_Polythene_Shielding_zAxisPos), PE_Shielding_Logical, "EJ2_Shielding", fpWorldLogical, false, 1,  true);

    
    
    
    
    //---------- For EJ#3 Paraffin Shielding 7.5" x 8.25" x 5.25" -------------
    
    const G4double       ProtonBeam_Line_zAxis_Pos = (31.4-1.5*2.54)*cm;

    const G4double     EJ3_Paraffin_Shielding_yAxis_Length = 7.5*2.54*cm;
    const G4double     EJ3_Paraffin_Shielding_zAxis_Length = 8.25*2.54*cm;
    const G4double     EJ3_Paraffin_Shielding_xAxis_Length = 5.25*2.54*cm;
    
    const G4double     EJ3_Paraffin_Shielding_yAxisPos = SSOuterDewar_OuterRadius + 11*2.54*cm - EJ3_Paraffin_Shielding_yAxis_Length/2;
    const G4double     EJ3_Paraffin_Shielding_xAxisPos = 1*2.54*cm + EJ3_Paraffin_Shielding_xAxis_Length/2;
    const G4double     EJ3_Paraffin_Shielding_zAxisPos = ProtonBeam_Line_zAxis_Pos - 18.5*2.54*cm + EJ3_Paraffin_Shielding_zAxis_Length;


    G4Box  *EJ3_Shielding_Box = new G4Box("EJ3_Shielding_Box", EJ3_Paraffin_Shielding_xAxis_Length/2, EJ3_Paraffin_Shielding_yAxis_Length/2, EJ3_Paraffin_Shielding_zAxis_Length/2);
    
    G4LogicalVolume*  EJ3_Shielding_Logical = new G4LogicalVolume(EJ3_Shielding_Box, paraffin, "EJ3_Shielding_Logical");
    
    new G4PVPlacement(0, G4ThreeVector(EJ3_Paraffin_Shielding_xAxisPos, EJ3_Paraffin_Shielding_yAxisPos, EJ3_Paraffin_Shielding_zAxisPos), EJ3_Shielding_Logical, "EJ3_Shielding", fpWorldLogical, false, 1,  true);
    
    

     PE_Shielding_Logical->SetVisAttributes(G4Colour::Grey());
    EJ3_Shielding_Logical->SetVisAttributes(G4Colour::Grey());

    

    //--------------------------------------------------------------------------------------------------------------
    
    
    
    
    
    
    //----------------------Aluminum Frame at the Bottom-----------------------------------------
    
    const G4double Al_Frame_Btm_Ground_Distance =   6*2.54*cm;
    
    const G4double Al_Frame_SideLength          = 3.5*2.54*cm;
    const G4double Al_Frame_xAxis_Length        =  32*2.54*cm;
    const G4double Al_Frame_yAxis_Length        =  31*2.54*cm;

    const G4double Al_Frame_LongBeam_yAxis_Length = 60*2.54*cm;

    
    const G4double Al_Frame_xPosition = (Al_Frame_xAxis_Length+Al_Frame_SideLength)/2;
    const G4double Al_Frame_zPosition = -TPC_Center_Ground_Distance-1.5*2.54*cm+Al_Frame_Btm_Ground_Distance+Al_Frame_SideLength/2;
    
    
    G4Box  *Al_Frame_Inner_Hollow = new G4Box("Al_Frame_Inner_Hollow", (Al_Frame_xAxis_Length-2*Al_Frame_SideLength)/2, (Al_Frame_yAxis_Length-2*Al_Frame_SideLength)/2, Al_Frame_SideLength);

    G4Box  *Al_Frame_Full_Box = new G4Box("Al_Frame_Full_Box", Al_Frame_xAxis_Length/2, Al_Frame_yAxis_Length/2, Al_Frame_SideLength/2);
    
       G4SubtractionSolid* Al_Frame_Square = new G4SubtractionSolid("Al_Frame_Square", Al_Frame_Full_Box, Al_Frame_Inner_Hollow, 0, G4ThreeVector(0.,0.,0));
    
                 G4Box  *Al_Frame_LongBeam = new G4Box("Al_Frame_LongBeam", Al_Frame_SideLength/2, Al_Frame_LongBeam_yAxis_Length/2, Al_Frame_SideLength/2);

    
    G4LogicalVolume   *Al_Frame_Square_Logic = new G4LogicalVolume(Al_Frame_Square, Al, "Al_Frame_Square_Logic");

    G4LogicalVolume *Al_Frame_LongBeam_Logic = new G4LogicalVolume(Al_Frame_LongBeam, Al, "Al_Frame_LongBeam_Logic");

    
    Al_Frame_Square_Logic->SetVisAttributes(Al_Color);
  Al_Frame_LongBeam_Logic->SetVisAttributes(Al_Color);
    
    new G4PVPlacement(0, G4ThreeVector(0., 0., -TPC_Center_Ground_Distance-1.5*2.54*cm+Al_Frame_Btm_Ground_Distance+Al_Frame_SideLength/2), Al_Frame_Square_Logic, "Al Main Frame Physical", fpWorldLogical, true, 1, true);

    new G4PVPlacement(0, G4ThreeVector(-Al_Frame_xPosition, 0., Al_Frame_zPosition), Al_Frame_LongBeam_Logic, "Al LongBeam Physical 1", fpWorldLogical, true, 1, true);
    
    new G4PVPlacement(0, G4ThreeVector(Al_Frame_xPosition, 0., Al_Frame_zPosition), Al_Frame_LongBeam_Logic, "Al LongBeam Physical 2", fpWorldLogical, true, 1, true);

    
    //-------------------------------------------------------------------------------------------

    
    
    
    
    
    
    
    //---------------------SS Poles to Hold the Dewars-------------------------------------------
    
    
    const G4double Front_Pole_Radius = 1.8/2*2.54*cm;
    const G4double  Side_Pole_Radius =   2/2*2.54*cm;

    const G4double Front_Pole_Proton_Beam_Distance = 8*2.54*cm;
    const G4double  Side_Pole_Proton_Beam_Distance = 8*2.54*cm + 6.*cm;
    
    const G4double Front_Pole_Height =  TPC_Center_Ground_Distance+zAxis_TPC_Center_BeamLine_Distance-Front_Pole_Proton_Beam_Distance-Al_Frame_Btm_Ground_Distance-Al_Frame_SideLength;
    const G4double  Side_Pole_Height =  TPC_Center_Ground_Distance+zAxis_TPC_Center_BeamLine_Distance-Side_Pole_Proton_Beam_Distance-Al_Frame_Btm_Ground_Distance-Al_Frame_SideLength;

    

        
    const G4double Front_Pole_xPosition = 6*2.54*cm + Front_Pole_Radius + SSOuterDewar_OuterRadius;
    const G4double Front_Pole_yPosition = 0.;
    const G4double Front_Pole_zPosition = Al_Frame_zPosition+Al_Frame_SideLength/2+Front_Pole_Height/2;
    
    const G4double Side_Pole_xPosition = 0.;
    const G4double Side_Pole_yPosition = (Al_Frame_yAxis_Length - Al_Frame_SideLength)/2;
    const G4double Side_Pole_zPosition = Al_Frame_zPosition+Al_Frame_SideLength/2+Side_Pole_Height/2;
    
    
    
    G4Tubs *SS_Front_Pole = new G4Tubs("SS_Front_Pole", 0, Front_Pole_Radius, Front_Pole_Height/2, 0.*deg, 360*deg);
    G4Tubs  *SS_Side_Pole = new G4Tubs("SS_Side_Pole",  0, Side_Pole_Radius, Side_Pole_Height/2, 0.*deg, 360*deg);

    
    G4LogicalVolume *SS_Front_Pole_Logic = new G4LogicalVolume(SS_Front_Pole, steel, "SS_Front_Pole_Logic");
    G4LogicalVolume *SS_Side_Pole_Logic =  new G4LogicalVolume(SS_Side_Pole, steel, "SS_Side_Pole_Logic");

     SS_Side_Pole_Logic->SetVisAttributes(Al_Color);    
    SS_Front_Pole_Logic->SetVisAttributes(Al_Color);
    
    
    new G4PVPlacement(0, G4ThreeVector(Front_Pole_xPosition, Front_Pole_yPosition, Front_Pole_zPosition), SS_Front_Pole_Logic, "Front Pole Physical", fpWorldLogical, true, 1, true);
    

    new G4PVPlacement(0, G4ThreeVector(-Front_Pole_xPosition, Front_Pole_yPosition, Front_Pole_zPosition), SS_Front_Pole_Logic, "Back Pole Physical", fpWorldLogical, true, 1, true);

    
    new G4PVPlacement(0, G4ThreeVector(Side_Pole_xPosition, Side_Pole_yPosition, Side_Pole_zPosition), SS_Side_Pole_Logic, "Side Pole Physical 1", fpWorldLogical, true, 1, true);

    new G4PVPlacement(0, G4ThreeVector(Side_Pole_xPosition, -Side_Pole_yPosition, Side_Pole_zPosition), SS_Side_Pole_Logic, "Side Pole Physical 2", fpWorldLogical, true, 1, true);

    
    //-------------------------------------------------------------------------------------------
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    //===============================================

	// Invisible world volume.
	fpWorldLogical->SetVisAttributes(G4VisAttributes::Invisible);
	
	
	//G4VisAttributes*  PMT_Attributes = new G4VisAttributes(G4Colour::Blue());  
	//PMT_Attributes->SetForceSolid(true);
    
    //===============================================
    
    // Create the TPC & EJs sensitive detector
    G4VSensitiveDetector* monitor = new BeamTestSensitiveDetector("Detector"); 
    
    // Create Periphery Volumes sensitive detector
    G4VSensitiveDetector* PeripheryDetector = new PeripherySensitiveDetector("PeripheryDetector"); 

    
    // Register detector with manager
    G4SDManager::GetSDMpointer()->AddNewDetector(monitor);
    G4SDManager::GetSDMpointer()->AddNewDetector(PeripheryDetector);
    

    
    
    ///*
    // Attach detector to volume defining calorimeter cells
    ActiveGAr_Volume_Logical->SetSensitiveDetector(monitor);
    ActiveLAr_Volume_Logical->SetSensitiveDetector(monitor);
             EJ_tube_Logical->SetSensitiveDetector(monitor);


    
      LAr_Volume_Logical->SetSensitiveDetector(PeripheryDetector);
      GAr_Volume_Logical->SetSensitiveDetector(PeripheryDetector);

    
    PTFE_BtmVolume_Logical->SetSensitiveDetector(PeripheryDetector);
    PTFE_TopVolume_Logical->SetSensitiveDetector(PeripheryDetector);
    
    
QuartzWindow_Volume_Logical->SetSensitiveDetector(PeripheryDetector);
          PMT_Shell_Logical->SetSensitiveDetector(PeripheryDetector); 
   PMT_QuartzWindow_Logical->SetSensitiveDetector(PeripheryDetector);
     //*/ 
    
    
    //---------------External Parts------------------
   ///* 
    
    
                          // fpWorldLogical->SetSensitiveDetector(PeripheryDetector);
    
         SSInnerDewar_CF600Flange_Logical->SetSensitiveDetector(PeripheryDetector);
        SSInnerDewar_CF1650Flange_Logical->SetSensitiveDetector(PeripheryDetector);
     SSInnerDewar_CF600FlangeRing_Logical->SetSensitiveDetector(PeripheryDetector);
    SSOuterDewar_CF1650FlangeRing_Logical->SetSensitiveDetector(PeripheryDetector);
                     SSInnerDewar_Logical->SetSensitiveDetector(PeripheryDetector);
                     SSOuterDewar_Logical->SetSensitiveDetector(PeripheryDetector);
                       HoldingBar_Logical->SetSensitiveDetector(PeripheryDetector);
    
    
         HeatExchanger_HoldingBar_Logical->SetSensitiveDetector(PeripheryDetector);
       HeatExchanger_HoldingPlate_Logical->SetSensitiveDetector(PeripheryDetector);
                    HeatExchanger_Logical->SetSensitiveDetector(PeripheryDetector);
    
                  Coldhead_Flange_Logical->SetSensitiveDetector(PeripheryDetector);
                 Coldehead_TopBox_Logical->SetSensitiveDetector(PeripheryDetector);
                  Coldhead_Column_Logical->SetSensitiveDetector(PeripheryDetector);
                         Coldhead_Logical->SetSensitiveDetector(PeripheryDetector);
                     Coldhead_Cup_Logical->SetSensitiveDetector(PeripheryDetector);
    
    
    EJ3_Shielding_Logical->SetSensitiveDetector(PeripheryDetector);
     PE_Shielding_Logical->SetSensitiveDetector(PeripheryDetector);

    
         SS_Side_Pole_Logic->SetSensitiveDetector(PeripheryDetector);
        SS_Front_Pole_Logic->SetSensitiveDetector(PeripheryDetector);
      Al_Frame_Square_Logic->SetSensitiveDetector(PeripheryDetector);
    Al_Frame_LongBeam_Logic->SetSensitiveDetector(PeripheryDetector);
    
    //delete SS_Color;
    
//*/
    
    //===============================================
    
  /*
  ////////////////////////////////////////////////////////////////////////
  // Visualisation attributes
 
  // Invisible world volume.
  fpWorldLogical->SetVisAttributes(G4VisAttributes::Invisible);

  // HandsOn4: Calorimeter attributes 
  // Invisible calorimeter mother volume
  //calorimeterLogical->SetVisAttributes(G4VisAttributes::Invisible);
 
  // Calorimeter cells - green with transparency
  G4VisAttributes* calorimeterAttributes =
     new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 0.1));
  calorimeterAttributes->SetForceSolid(true);
  //cellLogical->SetVisAttributes(calorimeterAttributes);
 cellLogical->SetVisAttributes(G4VisAttributes::Invisible);

  // Target Volume - light blue
  G4VisAttributes* targetAttributes = new G4VisAttributes(G4Colour(0.0,0.5,0.5,1.0));
  targetAttributes->SetForceSolid(true);
  //targetLogical->SetVisAttributes(targetAttributes);
 targetLogical->SetVisAttributes(G4VisAttributes::Invisible);


  // Evacuation chamber  - magenta
  G4VisAttributes* evacChamberAttributes = new G4VisAttributes(G4Colour::Magenta());
  evacChamberAttributes->SetForceSolid(true);
  //evacChamberLogical->SetVisAttributes(evacChamberAttributes);
 evacChamberLogical->SetVisAttributes(G4VisAttributes::Invisible);

  // Silicon Monitor  - cyan
  G4VisAttributes* siliconMonitorAttributes = new G4VisAttributes(G4Colour::Cyan());
  siliconMonitorAttributes->SetForceSolid(true);
  //siliconMonitorLogical->SetVisAttributes(siliconMonitorAttributes);
siliconMonitorLogical->SetVisAttributes(G4VisAttributes::Invisible);



  // Beam Pipe Vacuum - yellow
  G4VisAttributes* beamPipeAttributes = new G4VisAttributes(G4Colour::Yellow());
  beamPipeAttributes->SetForceSolid(true);
  //beamPipeLogical->SetVisAttributes(beamPipeAttributes);
  beamPipeLogical->SetVisAttributes(G4VisAttributes::Invisible);  
*/  
 

}
