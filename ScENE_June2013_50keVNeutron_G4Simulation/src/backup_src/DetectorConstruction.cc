
#include "DetectorConstruction.hh" 

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
#include "G4VisAttributes.hh"

#include "G4SDManager.hh"
#include "stdlib.h"
#include "string.h"

//#include "G4GDMLParser.hh"


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
  new G4Material(    "GasArgon",    z=18,  a=39.948*g/mole,    density=1.449e-3*g/cm3);
  new G4Material(      "Copper",    z=29,  a=63.546*g/mole,   density=8.94*g/cm3);
    
    
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
      G4Material* quartz = G4Material::GetMaterial("quartz");
       G4Material* EJ301 = G4Material::GetMaterial("EJ301");
        G4Material* PTFE = G4Material::GetMaterial("Teflon");
       G4Material* steel = G4Material::GetMaterial("SS304LSteel");
      G4Material* vacuum = G4Material::GetMaterial("Vacuum");
G4Material* polyethylene = G4Material::GetMaterial("polyethylene");
      G4Material* copper = G4Material::GetMaterial("Copper");
       G4Material* Kovar = G4Material::GetMaterial("Kovar");
      G4Material* SS316L = G4Material::GetMaterial("SS316LSteel");  

    
    
    
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
     2.  the total length of gas Ar is 8mm and liquid Ar is 68mm
     
     
     April 11th
     1. PE collimator is added in geometry.
     2. in particle gun setting, neutron is 50cm away from origion.
     3. particle gun offset in zAxis is 38cm. (total length of field cage is 76cm)
     4. PE collimator is set to 3cm away from neutron beam.
     
     
     May 30th
     1. The SS dewar and flange geometry are based on FNAL design.
     2. Coldhead and heat exchanger are also added.
     */
    //---------------------------------------------
    
    
    //----------- SS Dewar and flange --------------

    
    ///*
    //--- InnerDewar ---
    const G4double      SSInnerDewar_OuterRadius = 4*25.4*mm;    //-- 4" OD---
    const G4double      SSInnerDewar_Thickness = 0.065*25.4*mm;  //-- 0.065"---
    const G4double      SSInnerDewar_TotHeight = 30*25.4*mm;     //-- 30" height---
    

    //---CF600 Flange---
    const G4double      SSInnerDewar_CF600Flange_OuterRadius = 5.97*25.4*mm;
    const G4double      SSInnerDewar_CF600Flange_TotHeight = 0.78*25.4*mm;
    
    //--- OuterDewar ---
    const G4double      SSOuterDewar_OuterRadius = 14*25.4*mm;      //-- 14" OD --
    const G4double      SSOuterDewar_Thickness = 0.188*25.4*mm;
    const G4double      SSOuterDewar_TotHeight = 6*12*25.4*mm;      //-- 6' long --
    
    //--- CF1650 Flange ---
    const G4double      SSOuterDewar_CF1650Flange_OuterRadius = 16.5*25.4*mm;
    const G4double      SSOuterDewar_CF1650Flange_TotHeight = 1.12*25.4*mm;
    
    //--- SS bar 1” X 0.75” X 29” ---
    const G4double      DewarHoldingBar_Height = 29*25.4*mm;
    const G4double      DewarHoldingBar_Width  = 1*25.4*mm;
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
    

    
    
    
    G4Tubs  *SSInnerDewar_CF600Flange = new G4Tubs("SSInnerDewar_CF600Flange", 0.*mm, SSInnerDewar_CF600Flange_OuterRadius, SSInnerDewar_CF600Flange_TotHeight/2, 0.*deg, 360*deg);
    G4Tubs *SSOuterDewar_CF1650Flange = new G4Tubs("SSOuterDewar_CF1650Flange", 0.*mm, SSOuterDewar_CF1650Flange_OuterRadius, SSOuterDewar_CF1650Flange_TotHeight/2, 0.*deg, 360*deg);

    
    G4LogicalVolume  *SSInnerDewar_CF600Flange_Logical = new G4LogicalVolume(SSInnerDewar_CF600Flange, steel, "SSInnerDewar_CF600Flange_Logical");
    G4LogicalVolume *SSInnerDewar_CF1650Flange_Logical = new G4LogicalVolume(SSOuterDewar_CF1650Flange, steel, "SSInnerDewar_CF1650Flange_Logical");
    
    
    new G4PVPlacement(0, G4ThreeVector(SSInnerDewar_xCenterPosition, SSInnerDewar_yCenterPosition, (SSInnerDewar_zCenterPosition+SSInnerDewar_TotHeight/2+SSInnerDewar_CF600Flange_TotHeight/2)), SSInnerDewar_CF600Flange_Logical, "SSInnerDewar Top CF600Flange Physical", fpWorldLogical, true, 0);
    
    new G4PVPlacement(0, G4ThreeVector(SSInnerDewar_xCenterPosition, SSInnerDewar_yCenterPosition, (SSInnerDewar_zCenterPosition-SSInnerDewar_TotHeight/2-SSInnerDewar_CF600Flange_TotHeight/2)), SSInnerDewar_CF600Flange_Logical, "SSInnerDewar Bottom CF600Flange Physical", fpWorldLogical, true, 0);
    
    new G4PVPlacement(0, G4ThreeVector(SSOuterDewar_xCenterPosition, SSOuterDewar_yCenterPosition, (SSOuterDewar_zCenterPosition+SSOuterDewar_TotHeight/2+SSOuterDewar_CF1650Flange_TotHeight/2)), SSInnerDewar_CF1650Flange_Logical, "SSInnerDewar Top CF1650Flange Physical", fpWorldLogical, true, 0);

    new G4PVPlacement(0, G4ThreeVector(SSOuterDewar_xCenterPosition, SSOuterDewar_yCenterPosition, (SSOuterDewar_zCenterPosition-SSOuterDewar_TotHeight/2-SSOuterDewar_CF1650Flange_TotHeight/2)), SSInnerDewar_CF1650Flange_Logical, "SSInnerDewar Bottom CF1650Flange Physical", fpWorldLogical, true, 0);

    
    
                     
    G4Tubs   *SSInnerDewar_CF600FlangeRing = new G4Tubs("SSInnerDewar_CF600FlangeRing", SSInnerDewar_OuterRadius, SSInnerDewar_CF600Flange_OuterRadius, SSInnerDewar_CF600Flange_TotHeight/2, 0.*deg, 360*deg);
    G4Tubs  *SSOuterDewar_CF1650FlangeRing = new G4Tubs("SSOuterDewar_CF1650FlangeRing", SSOuterDewar_OuterRadius, SSOuterDewar_CF1650Flange_OuterRadius, SSOuterDewar_CF1650Flange_TotHeight/2, 0.*deg, 360*deg);

                      
    G4LogicalVolume  *SSInnerDewar_CF600FlangeRing_Logical = new G4LogicalVolume(SSInnerDewar_CF600FlangeRing, steel, "SSInnerDewar_CF600FlangeRing_Logical");
    G4LogicalVolume *SSOuterDewar_CF1650FlangeRing_Logical = new G4LogicalVolume(SSOuterDewar_CF1650FlangeRing, steel, "SSOuterDewar_CF1650FlangeRing_Logical");
       
    
    new G4PVPlacement(0, G4ThreeVector(SSInnerDewar_xCenterPosition, SSInnerDewar_yCenterPosition, (SSInnerDewar_zCenterPosition+SSInnerDewar_TotHeight/2-SSInnerDewar_CF600Flange_TotHeight/2)), SSInnerDewar_CF600FlangeRing_Logical, "SSInnerDewar Top CF600FlangeRing Physical", fpWorldLogical, true, 0);
                                        
    new G4PVPlacement(0, G4ThreeVector(SSInnerDewar_xCenterPosition, SSInnerDewar_yCenterPosition, (SSInnerDewar_zCenterPosition-SSInnerDewar_TotHeight/2+SSInnerDewar_CF600Flange_TotHeight/2)), SSInnerDewar_CF600FlangeRing_Logical, "SSInnerDewar Bottom CF600FlangeRing Physical", fpWorldLogical, true, 0);
                                   
    new G4PVPlacement(0, G4ThreeVector(SSOuterDewar_xCenterPosition, SSOuterDewar_yCenterPosition, (SSOuterDewar_zCenterPosition+SSOuterDewar_TotHeight/2-SSOuterDewar_CF1650Flange_TotHeight/2)), SSOuterDewar_CF1650FlangeRing_Logical, "SSOuterDewar Top CF1650FlangeRing Physical", fpWorldLogical, true, 0);
    
    new G4PVPlacement(0, G4ThreeVector(SSOuterDewar_xCenterPosition, SSOuterDewar_yCenterPosition, (SSOuterDewar_zCenterPosition-SSOuterDewar_TotHeight/2+SSOuterDewar_CF1650Flange_TotHeight/2)), SSOuterDewar_CF1650FlangeRing_Logical, "SSOuterDewar Bottom CF1650FlangeRing Physical", fpWorldLogical, true, 0);

    
    
    
    G4Tubs *SSInnerDewar = new G4Tubs("SSInnerDewar", SSInnerDewar_OuterRadius-SSInnerDewar_Thickness, SSInnerDewar_OuterRadius, SSInnerDewar_TotHeight/2, 0.*deg, 360*deg);
    G4Tubs *SSOuterDewar = new G4Tubs("SSOuterDewar", SSOuterDewar_OuterRadius-SSOuterDewar_Thickness, SSOuterDewar_OuterRadius, SSOuterDewar_TotHeight/2, 0.*deg, 360*deg);

                      
    G4LogicalVolume *SSInnerDewar_Logical = new G4LogicalVolume(SSInnerDewar, steel, "SSInnerDewar_Logical");
    G4LogicalVolume *SSOuterDewar_Logical = new G4LogicalVolume(SSOuterDewar, steel, "SSOuterDewar_Logical");
    
    
    new G4PVPlacement(0, G4ThreeVector(SSInnerDewar_xCenterPosition, SSInnerDewar_yCenterPosition, SSInnerDewar_zCenterPosition), SSInnerDewar_Logical, "SSInnerDewar Physical", fpWorldLogical, true, 0);
    new G4PVPlacement(0, G4ThreeVector(SSOuterDewar_xCenterPosition, SSOuterDewar_yCenterPosition, SSOuterDewar_zCenterPosition), SSOuterDewar_Logical, "SSOuterDewar Physical", fpWorldLogical, true, 0);
    
    
    G4Box *HoldingBar = new G4Box("HoldingBar", DewarHoldingBar_Width/2, DewarHoldingBar_Long/2, DewarHoldingBar_Height/2);
    
    G4LogicalVolume *HoldingBar_Logical = new G4LogicalVolume(HoldingBar, steel, "HoldingBar_Logical");
                 
    new G4PVPlacement(0, G4ThreeVector(HoldingBar_xPosition, HoldingBar_yPosition, HoldingBar_zPosition), HoldingBar_Logical, "SSInnerDewar Top HoldingBar Physical", fpWorldLogical, true, 0);
    new G4PVPlacement(0, G4ThreeVector(HoldingBar_xPosition, -HoldingBar_yPosition, HoldingBar_zPosition), HoldingBar_Logical, "SSInnerDewar Bottom HoldingBar Physical", fpWorldLogical, true, 0);
    
    
    
    
    
    //------------------------heat exchanger and attachments----------------------------------------

    //-- heat exchanger geometry, back and right angle is set to position x=0, y=0---
    
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

    
    new G4PVPlacement(0, G4ThreeVector(-127*mm/2-6*25.4*mm/2+1.25*25.4*mm/2, -0.625*25.4*mm-0.75*25.4*mm/2, SSOuterDewar_zCenterPosition+SSOuterDewar_TotHeight/2-15*25.4/2*mm), HeatExchanger_HoldingBar_Logical, "HeatExchanger HoldingBar Physical", fpWorldLogical, true, 0);
    
    new G4PVPlacement(0, G4ThreeVector(-127*mm/2, -0.75*25.4*mm/2, SSOuterDewar_zCenterPosition+SSOuterDewar_TotHeight/2-15*25.4*mm-3.25*25.4*mm+8.25*25.4*mm/2), HeatExchanger_HoldingPlate_Logical, "HeatExchanger HoldingPlate Physical", fpWorldLogical, true, 0);
    
    new G4PVPlacement(0, G4ThreeVector(-127*mm/2, 55*mm/2, SSOuterDewar_zCenterPosition+SSOuterDewar_TotHeight/2-15*25.4*mm-3.25*25.4*mm+8.25*25.4*mm/2), HeatExchanger_Logical, "HeatExchanger Physical", fpWorldLogical, true, 0);


    
    //---------------------Coldhead, PT60 from Cryomech-------------------------------------------
    
    
    const G4double  Coldehead_xCenterPosition = -SSInnerDewar_xCenterPosition + SSInnerDewar_OuterRadius;
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

    
    new G4PVPlacement(0, G4ThreeVector(Coldehead_xCenterPosition, Coldehead_yCenterPosition, Coldehead_zCenterPosition+12.7*mm/2), Coldhead_Flange_Logical, "Coldhead SSFlange Physical", fpWorldLogical, true, 0);

    new G4PVPlacement(0, G4ThreeVector(Coldehead_xCenterPosition, Coldehead_yCenterPosition, Coldehead_zCenterPosition+12.7*mm+205.6*mm/2), Coldehead_TopBox_Logical, "Coldehead TopBox Physical", fpWorldLogical, true, 0);
    
    new G4PVPlacement(0, G4ThreeVector(Coldehead_xCenterPosition - 70.6*mm/4, Coldehead_yCenterPosition, SSOuterDewar_zCenterPosition + SSOuterDewar_TotHeight/2 - Coldehead_ColumnHeight/2), Coldhead_Column_Logical, "Coldhead Column Physical 1", fpWorldLogical, true, 0);
    
    new G4PVPlacement(0, G4ThreeVector(Coldehead_xCenterPosition + 70.6*mm/4, Coldehead_yCenterPosition, SSOuterDewar_zCenterPosition + SSOuterDewar_TotHeight/2 - Coldehead_ColumnHeight/2), Coldhead_Column_Logical, "Coldhead Column Physical 2", fpWorldLogical, true, 0);
    
    
    

    const G4int     Coldhead_zPlaneNbr = 7;
    const G4double  Coldehead_zPlane[Coldhead_zPlaneNbr] = {0.*mm,           -0.31*25.4*mm,    -0.31*25.4*mm,       (-0.31-0.41)*25.4*mm,     (-0.31-0.41)*25.4*mm,     (-0.31-0.41-0.65)*25.4*mm, (-0.31-1.76)*25.4*mm};
    const G4double  Coldehead_rInner[Coldhead_zPlaneNbr] = {0.*mm,           0.*mm,            0.*mm,               0.*mm,                       0.*mm,                   0.*mm,        0.*mm};
    const G4double  Coldehead_rOuter[Coldhead_zPlaneNbr] = {2.78/2*25.4*mm,  2.78/2*25.4*mm,   4.348/2*25.4*mm,    4.348/2*25.4*mm,             4/2*25.4*mm,              4/2*25.4*mm,  0.*mm};
    

    G4Polycone              *Coldhead = new G4Polycone("Coldhead", 0.*deg, 360.0*deg, Coldhead_zPlaneNbr, Coldehead_zPlane, Coldehead_rInner, Coldehead_rOuter);
    
    G4LogicalVolume* Coldhead_Logical = new G4LogicalVolume(Coldhead ,	 copper,  "Coldhead_Logical");   
    
    new G4PVPlacement(0,  G4ThreeVector(Coldehead_xCenterPosition, Coldehead_yCenterPosition, SSOuterDewar_zCenterPosition + SSOuterDewar_TotHeight/2 - Coldehead_ColumnHeight),  Coldhead_Logical,   "Coldhead Physical",  fpWorldLogical,  true,  0);	 
  
    
    
    const G4int     Coldhead_Cup_zPlaneNbr = 4;
    const G4double  Coldhead_Cup_zPlane[Coldhead_Cup_zPlaneNbr] = {0.*mm,                -4.7*25.4*mm,            -4.7*25.4*mm,         (-4.7-0.08)*25.4*mm};
    const G4double  Coldhead_Cup_rInner[Coldhead_Cup_zPlaneNbr] = {2*25.4*mm,             2*25.4*mm,               0.*mm,                0.*mm};
    const G4double  Coldhead_Cup_rOuter[Coldhead_Cup_zPlaneNbr] = {(2+0.08)*25.4*mm,     (2+0.08)*25.4*mm,         (2+0.08)*25.4*mm,     (2+0.08)*25.4*mm};
    
    
    G4Polycone              *Coldhead_Cup = new G4Polycone("Coldhead_Cup", 0.*deg, 360.0*deg, Coldhead_Cup_zPlaneNbr, Coldhead_Cup_zPlane, Coldhead_Cup_rInner, Coldhead_Cup_rOuter);
    
    G4LogicalVolume* Coldhead_Cup_Logical = new G4LogicalVolume(Coldhead_Cup ,	 steel,  "Coldhead_Cup_Logical");   
    
    new G4PVPlacement(0,  G4ThreeVector(Coldehead_xCenterPosition, Coldehead_yCenterPosition, SSOuterDewar_zCenterPosition + SSOuterDewar_TotHeight/2 - Coldehead_ColumnHeight-0.72*25.4*mm),  Coldhead_Cup_Logical,   "Coldhead Cup Physical",  fpWorldLogical,  true,  0);
    
    
    //-------------------------------------------------------------------------------------
    //*/
    
    
    
    
    //---------------------------TPC components and PMT parts-------------------------------------
    
    
    
    const G4double      CathodeCap_TotLength = 15.4*mm;
    const G4double    CathodeCap_OuterRadius = 43*mm;
    
          G4double   QuartzWindow_OuterRadius = 40*mm;
    const G4double     QuartzWindow_TotLength = 0.5*mm;
    
    const G4double    AnodeCap_TotLength = 14*mm;
    const G4double  AnodeCap_OuterRadius = CathodeCap_OuterRadius;
    const G4double  AnodeCap_InnerRadius = QuartzWindow_OuterRadius;
    const G4double   AnodeCap_PTFELength = 5*mm;

    const G4double    PTFEWall_InnerRadius = 38.2*mm;
    const G4double    PTFEWall_OuterRadius = 43*mm;
    const G4double      PTFEWall_TotLength = 68*mm;
    
    const G4double     FieldCage_InnerRadius = 32*mm;
    const G4double     FieldCage_OuterRadius = 74.2/2*mm;
    const G4double       FieldCage_TotLength = 76.*mm;
    const G4double  zFieldCageTop_LArSurface = 8.*mm;
    
    const G4double  CopperRing_TotLength = 3.*mm;
    const G4double  CopperRing_Thickness = 0.5*mm;
    const G4double   CopperRing_Distance = 4.*mm;
    
    const G4double  CopperRing_InnerRadius = FieldCage_OuterRadius;
    const G4double  CopperRing_OuterRadius = CopperRing_InnerRadius + CopperRing_Thickness;
    
    const G4double        ArgonRadius = FieldCage_InnerRadius;
    const G4double      GAr_TotLength = 8.*mm;
    const G4double      LAr_TotLength = FieldCage_TotLength - GAr_TotLength;  //-- also equal to PTFEWall_TotLength---
    
    
    
    //--- For PE Collimator 22cm x OD-22cm  ID-2.54cm 
    const G4double     Polythene_Cylinder_Inner_Radius = 2.54/2*cm;
    const G4double     Polythene_Cylinder_Outer_Radius = 11*cm;
    const G4double     Polythene_Cylinder_Half_Length  = 11*cm;
    const G4double     Polythene_Cylinder_xAxis_Distance = -70*cm + Polythene_Cylinder_Half_Length + 3*cm;
    const G4double     Polythene_Cylinder_zAxis_offset = -38.0*mm;
    
    
    //--- for 2"DIA x 2" Chamber of Scintillator, neutron detector ---  
    const G4double        EJ301_Radius = 2.54*cm;     
    const G4double        EJ301_halfHeight = 2.54*cm;
    
    const G4double       Distance_EJ301_Away = 50*cm;
    const G4double       EJ301_zOffset = -(FieldCage_TotLength/2 - GAr_TotLength);
    

    
    const G4double  UpperQuartzWindow_zOffset = QuartzWindow_TotLength/2 + GAr_TotLength;
    const G4double  LowerQuartzWindow_zOffset = QuartzWindow_TotLength/2 + LAr_TotLength;
    
    const G4double     FieldCage_zOffset = FieldCage_TotLength/2 - zFieldCageTop_LArSurface;
    
    const G4double	  CathodeCap_zOffset = -PTFEWall_TotLength - CathodeCap_TotLength - QuartzWindow_TotLength;
    
    const G4double           PMT_zOffset = -LAr_TotLength - 2*mm;

   
    
    
    
    //----PE Collimator ---------
    
    G4Tubs* Collimator_tube = new G4Tubs("Collimator_tube", Polythene_Cylinder_Inner_Radius,  Polythene_Cylinder_Outer_Radius,  Polythene_Cylinder_Half_Length,  0.*deg,  360.*deg);

    G4LogicalVolume*  Collimator_tube_Logical = new G4LogicalVolume(Collimator_tube, polyethylene, "Collimator_tube_Logical");

    G4RotationMatrix R_collimator;
    G4ThreeVector    T_collimator;
    
    R_collimator.rotateY(90.0*deg);
    T_collimator.setX(Polythene_Cylinder_xAxis_Distance); T_collimator.setY(0); T_collimator.setZ(Polythene_Cylinder_zAxis_offset);
    
    new G4PVPlacement(G4Transform3D(R_collimator, T_collimator),  Collimator_tube_Logical, "Collimator_tube_Physical", fpWorldLogical, false, 0);

   
    
    
        
    
    
    
    //----LAr and GAr ------
    
    
    G4Tubs *LiquidArgon = new G4Tubs("LiquidArgon", 0.*mm, ArgonRadius, LAr_TotLength/2, 0.*deg, 360*deg);
    
    G4LogicalVolume *LiquidArgon_Logical = new G4LogicalVolume(LiquidArgon, LAr, "LiquidArgon_Logical");
    
    //new G4PVPlacement(0, G4ThreeVector(0,0,-LAr_TotLength/2), LiquidArgon_Logical, "LiquidArgon_Physical", fpWorldLogical, true, 0);

    
    
    G4Tubs *GasArgon = new G4Tubs("GasArgon", 0.*mm, ArgonRadius, GAr_TotLength/2, 0.*deg, 360*deg);
    
    G4LogicalVolume *GasArgon_Logical = new G4LogicalVolume(GasArgon, GAr, "GasArgon_Logical");
    
   // new G4PVPlacement(0, G4ThreeVector(0,0, GAr_TotLength/2), GasArgon_Logical, "GasArgon_Physical", fpWorldLogical, true, 0);
    
    
    
    G4Tubs *FullLiquidArgon = new G4Tubs("FullLiquidArgon", 0.*mm, ArgonRadius, FieldCage_TotLength/2, 0.*deg, 360*deg);
    
    G4LogicalVolume *FullLiquidArgon_Logical = new G4LogicalVolume(FullLiquidArgon, LAr, "FullLiquidArgon_Logical");
    
    new G4PVPlacement(0, G4ThreeVector(0,0,-(FieldCage_TotLength/2-GAr_TotLength)), FullLiquidArgon_Logical, "FullLiquidArgon_Physical", fpWorldLogical, true, 0);
            
    //----------------------
    
    
    
    
    //----Anode Cap -------
    
    const G4int     AnodeCapPlaneNbr = 4;
    const G4double   AnodeCapPlane[AnodeCapPlaneNbr] = {0.*mm, (AnodeCap_TotLength-AnodeCap_PTFELength), (AnodeCap_TotLength-AnodeCap_PTFELength), AnodeCap_TotLength};
    const G4double  AnodeCaprInner[AnodeCapPlaneNbr] = {AnodeCap_InnerRadius, AnodeCap_InnerRadius, 0.0*mm, 0.0*mm};
    const G4double  AnodeCaprOuter[AnodeCapPlaneNbr] = {AnodeCap_OuterRadius, AnodeCap_OuterRadius, AnodeCap_OuterRadius, AnodeCap_OuterRadius};
    
    
    G4Polycone *AnodeCap = new G4Polycone("AnodeCap", 0.*deg, 360.0*deg, AnodeCapPlaneNbr, AnodeCapPlane, AnodeCaprInner, AnodeCaprOuter);
    
    G4LogicalVolume *AnodeCap_Logical = new G4LogicalVolume(AnodeCap, PTFE, "AnodeCap_Logical");
    
    new G4PVPlacement(0, G4ThreeVector(0,0,0), AnodeCap_Logical, "AnodeCap_Physical", fpWorldLogical, true, 0);
  
    //----------------------
    
    
    
    //----PTFE Wall ------
    
    
    G4Tubs *PTFEWall = new G4Tubs("PTFEWall", PTFEWall_InnerRadius, PTFEWall_OuterRadius, PTFEWall_TotLength/2, 0.*deg, 360*deg);

    G4LogicalVolume *PTFEWall_Logical = new G4LogicalVolume(PTFEWall, polyethylene, "PTFEWall_Logical");
    
    new G4PVPlacement(0, G4ThreeVector(0,0,-PTFEWall_TotLength/2), PTFEWall_Logical, "PTFEWall_Physical", fpWorldLogical, true, 0);

    //--------------------
    
    
  
	
    //----- Field Cage -----
    
    G4Tubs *FieldCage = new G4Tubs("FieldCage", FieldCage_InnerRadius, FieldCage_OuterRadius, FieldCage_TotLength/2, 0.*deg, 360*deg);
    
    G4LogicalVolume *FieldCage_Logical = new G4LogicalVolume(FieldCage, polyethylene, "FieldCage_Logical");
    
    new G4PVPlacement(0, G4ThreeVector(0,0, -FieldCage_zOffset), FieldCage_Logical, "FieldCage_Physical", fpWorldLogical, true, 0);
   //----------------------
    
    
    
    
    
    //---- Copper Ring (18 in total)----
    
    G4Tubs *CopperRing = new G4Tubs("CopperRing", CopperRing_InnerRadius, CopperRing_OuterRadius, CopperRing_TotLength/2, 0.*deg, 360*deg);
    
    G4LogicalVolume *CopperRing_Logical = new G4LogicalVolume(CopperRing, copper, "CopperRing_Logical");
    
    G4double InitialPosition = 1*CopperRing_Distance;
    
    
    for(int i=0; i<18; i++)
    {
        char *VolumeName = new char[30];
        sprintf(VolumeName, "CopperRing_Physical %i", i);
        
        G4double Distance = InitialPosition-i*CopperRing_Distance;
        
    new G4PVPlacement(0, G4ThreeVector(0,0,Distance), CopperRing_Logical, VolumeName, fpWorldLogical, true, 0);
    }
    //----------------------
    
    
    
    //----Cathode Cap------
    
	const G4int		CathodeCapPlaneNbr = 4;
	const G4double	CathodeCapPlane[CathodeCapPlaneNbr] = {0.*mm,  CathodeCap_TotLength-1*mm, CathodeCap_TotLength-1*mm, CathodeCap_TotLength};
	const G4double       CathodeCaprInner[CathodeCapPlaneNbr] = {39.37*mm, 39.37*mm, 31.7*mm, 31.7*mm};
	const G4double       CathodeCaprOuter[CathodeCapPlaneNbr] = {CathodeCap_OuterRadius, CathodeCap_OuterRadius, CathodeCap_OuterRadius, CathodeCap_OuterRadius};
	
	
	G4Polycone *CathodeCap = new G4Polycone("CathodeCap", 0.*deg, 360.0*deg, CathodeCapPlaneNbr, CathodeCapPlane, CathodeCaprInner, CathodeCaprOuter);
    
    G4LogicalVolume *CathodeCap_Logical = new G4LogicalVolume(CathodeCap, PTFE, "CathodeCap_Logical");
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, CathodeCap_zOffset), CathodeCap_Logical, "CathodeCap_Physical", fpWorldLogical, true, 0);
	
   //--------------------
	
    
    
    //-----quartz window -------
    
    QuartzWindow_OuterRadius = CopperRing_OuterRadius; //--- to prevent the overlap, set quartzWindowRadius the same as CopperRingRadius ----
    
    G4Tubs *QuartzWindow = new G4Tubs("QuartzWindow", 0.*mm, QuartzWindow_OuterRadius, QuartzWindow_TotLength/2, 0.*deg, 360*deg);
    
    G4LogicalVolume *QuartzWindow_Logical = new G4LogicalVolume(QuartzWindow, quartz, "QuartzWindow_Logical");
    
    new G4PVPlacement(0, G4ThreeVector(0,0, UpperQuartzWindow_zOffset), QuartzWindow_Logical, "UpperQuartzWindow_Physical", fpWorldLogical, true, 0);
    
    new G4PVPlacement(0, G4ThreeVector(0,0, -LowerQuartzWindow_zOffset), QuartzWindow_Logical, "LowerQuartzWindow_Physical", fpWorldLogical, true, 0);
    
    //--------------------
	

    
    
    //----3" PMT Solid Part ----
    
    const G4int     PMTzPlaneNbr = 8;
    const G4double  PMTShellThickness = 0.5*mm;
    const G4double  PMTzPlane[PMTzPlaneNbr] = {0.*mm+PMT_zOffset, -17.4115*mm+PMT_zOffset, -35.041*mm+PMT_zOffset, -39.9343*mm+PMT_zOffset, -116.*mm+PMT_zOffset, -116.078*mm+PMT_zOffset, -116.713*mm+PMT_zOffset, -123.063*mm+PMT_zOffset};
    const G4double  PMTrInner[PMTzPlaneNbr] = {36.576*mm, 36.576*mm, 28.0201*mm, 25.146*mm, 25.146*mm, 0.*mm, 0.*mm, 25.126*mm};
    const G4double  PMTrOuter[PMTzPlaneNbr] = {36.576*mm+PMTShellThickness, 36.576*mm+PMTShellThickness, 28.0201*mm+PMTShellThickness, 25.146*mm+PMTShellThickness, 25.146*mm+PMTShellThickness, 25.146*mm+PMTShellThickness, 25.146*mm+PMTShellThickness, 25.146*mm+PMTShellThickness};
    
    const G4double  PMT_QuartzWindow_Radius = 35.576*mm;
    const G4double  PMT_QuartzWindow_TotLength = 3.4*mm;
    
    
    G4Polycone              *PMT_Shell = new G4Polycone("PMT_Shell", 0.*deg, 360.0*deg, PMTzPlaneNbr, PMTzPlane, PMTrInner, PMTrOuter);
    
    G4LogicalVolume* PMT_Shell_Logical = new G4LogicalVolume(PMT_Shell,	 Kovar,  "PMT_Shell_Logical");   
    
                                         new G4PVPlacement(0,  G4ThreeVector(0,0,0),  PMT_Shell_Logical,   "PMT_Shell_Physical",  fpWorldLogical,  true,  0);	                           
    
    
    G4Tubs                 * PMT_QuartzWindow = new G4Tubs("PMT_QuartzWindow", 0.*mm, PMT_QuartzWindow_Radius, PMT_QuartzWindow_TotLength/2, 0.*deg, 360*deg);
    
    G4LogicalVolume *PMT_QuartzWindow_Logical = new G4LogicalVolume(PMT_QuartzWindow, quartz, "PMT_QuartzWindow_Logical");
    
                                                new G4PVPlacement(0, G4ThreeVector(0,0, PMT_zOffset-PMT_QuartzWindow_TotLength/2), PMT_QuartzWindow_Logical, "PMT_ QuartzWindow_Physical", fpWorldLogical, true, 0);

    
    //----------------------

    
    
    
    //-----------EJ301, neutron detector-----
    
    
    G4Tubs*                        EJ301_tube = new G4Tubs("EJ301_tube", 0.*cm, EJ301_Radius, EJ301_halfHeight, 0.*deg, 360.*deg);
    
    G4LogicalVolume*       EJ301_tube_Logical = new G4LogicalVolume(EJ301_tube, EJ301, "EJ301_tube_Logical");
    
    
    G4RotationMatrix Ra, Rb;
    Ra.rotateX(90.0*deg);
    Rb.rotateX(90.0*deg);
    
    G4ThreeVector Ta, Tb, Tc;
    Ta.setX(0.*cm); Ta.setY(Distance_EJ301_Away); Ta.setZ(EJ301_zOffset);
    Tb.setX(0.*cm); Tb.setY(-Distance_EJ301_Away); Tb.setZ(EJ301_zOffset);
   
    
    ///*
    new G4PVPlacement(G4Transform3D(Rb,Ta),
                      EJ301_tube_Logical, "Top_EJ301_tube_Physical",
                      fpWorldLogical,false,0);
    
    new G4PVPlacement(G4Transform3D(Ra,Tb),
                      EJ301_tube_Logical, "Bottom_EJ301_tube_Physical",
                      fpWorldLogical,false,0);
    //*/
    
    
    
    //---------------------------------------
    
        
    //===============================================

	// Invisible world volume.
	fpWorldLogical->SetVisAttributes(G4VisAttributes::Invisible);
	
	
	G4VisAttributes*  PMT_Attributes = new G4VisAttributes(G4Colour::Blue());  
	PMT_Attributes->SetForceSolid(true);
    
	G4VisAttributes*  PTFE_Attributes = new G4VisAttributes(G4Colour::Green()); 
	//polyethylene_Attributes->Ser
    
	PMT_Shell_Logical->SetVisAttributes(PMT_Attributes);
	
	AnodeCap_Logical->SetVisAttributes(PTFE_Attributes);
	CathodeCap_Logical->SetVisAttributes(PTFE_Attributes);
	PTFEWall_Logical->SetVisAttributes(PTFE_Attributes);
	
	
    FieldCage_Logical->SetVisAttributes(G4Colour::White());
    

    EJ301_tube_Logical->SetVisAttributes(G4Colour::Red());
    
/*
    GasArgon_Logical->SetVisAttributes(G4Colour::Yellow());
	FullLiquidArgon_Logical->SetVisAttributes(G4Colour::Yellow());
    LiquidArgon_Logical->SetVisAttributes(G4Colour::Yellow());
*/
    
  /*  
           Coldhead_Logical->SetVisAttributes(G4Colour::Yellow());
    Coldhead_Column_Logical->SetVisAttributes(G4Colour::Yellow());
         CopperRing_Logical->SetVisAttributes(G4Colour::Yellow());
    
    
    HeatExchanger_Logical->SetVisAttributes(G4Colour::Grey());
   */
    
    //===============================================
    
    // Create a new BetamTestSiliconMonitor sensitive detector
    G4VSensitiveDetector* monitor = new BeamTestSensitiveDetector("Detector"); 
    
    // Create a new BetamTestSiliconMonitor sensitive detector
    G4VSensitiveDetector* PeripheryDetector = new PeripherySensitiveDetector("PeripheryDetector"); 

    // Get pointer to detector manager
    
    // Register detector with manager
    G4SDManager::GetSDMpointer()->AddNewDetector(monitor);
    G4SDManager::GetSDMpointer()->AddNewDetector(PeripheryDetector);
    
    // Attach detector to volume defining calorimeter cells
    FullLiquidArgon_Logical->SetSensitiveDetector(monitor);
         EJ301_tube_Logical->SetSensitiveDetector(monitor);


  Collimator_tube_Logical->SetSensitiveDetector(PeripheryDetector);
    
      AnodeCap_Logical->SetSensitiveDetector(PeripheryDetector);
    CathodeCap_Logical->SetSensitiveDetector(PeripheryDetector);
      PTFEWall_Logical->SetSensitiveDetector(PeripheryDetector);
     FieldCage_Logical->SetSensitiveDetector(PeripheryDetector);
    CopperRing_Logical->SetSensitiveDetector(PeripheryDetector);

    QuartzWindow_Logical->SetSensitiveDetector(PeripheryDetector);
       PMT_Shell_Logical->SetSensitiveDetector(PeripheryDetector);
PMT_QuartzWindow_Logical->SetSensitiveDetector(PeripheryDetector);
 
    
    
    //----External Parts--------
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



  // Beam Window - brown
  G4VisAttributes* beamWindowAttributes = new G4VisAttributes(G4Colour(0.5,0.0,0.0,1.0));
  beamWindowAttributes->SetForceSolid(true);
  //beamWindowLogical->SetVisAttributes(beamWindowAttributes);
  beamWindowLogical->SetVisAttributes(G4VisAttributes::Invisible);


  // Beam Pipe Vacuum - yellow
  G4VisAttributes* beamPipeAttributes = new G4VisAttributes(G4Colour::Yellow());
  beamPipeAttributes->SetForceSolid(true);
  //beamPipeLogical->SetVisAttributes(beamPipeAttributes);
  beamPipeLogical->SetVisAttributes(G4VisAttributes::Invisible);  
*/  
 

}
