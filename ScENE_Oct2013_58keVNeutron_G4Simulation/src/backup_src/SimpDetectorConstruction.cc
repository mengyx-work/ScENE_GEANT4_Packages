
#include "SimpDetectorConstruction.hh" 
//#include "BeamTestCellParameterisation.hh"
//#include "BeamTestEmCalorimeter.hh"

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

    
    
    
    
    G4double world_Radius = 3*m;
    G4double world_halfHeight = 1.5*m;

 G4Tubs* worldSolid = new G4Tubs("World_Solid", 0*m, world_Radius, world_halfHeight, 0*deg, 360*deg);
     fpWorldLogical = new G4LogicalVolume(worldSolid, air, "World_Logical");
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
     */
    //---------------------------------------------
    
    const G4double     FieldCage_InnerRadius = 32*mm;
    const G4double     FieldCage_OuterRadius = 74.2/2*mm;
    const G4double       FieldCage_TotLength = 76.*mm;
    const G4double  zFieldCageTop_LArSurface = 8.*mm;
    const G4double               ArgonRadius = FieldCage_InnerRadius;
    const G4double             GAr_TotLength = 8.*mm;
    const G4double             LAr_TotLength = FieldCage_TotLength - GAr_TotLength;  //-- also equal to PTFEWall_TotLength---
    
    
    
    //--- For PE Collimator 22cm x OD-22cm  ID-2.54cm 
    const G4double     Polythene_Cylinder_Inner_Radius = 2.54/2*cm;
    const G4double     Polythene_Cylinder_Outer_Radius = 11*cm;
    const G4double     Polythene_Cylinder_Half_Length  = 11*cm;
    const G4double     Polythene_Cylinder_xAxis_Distance = -50*cm + Polythene_Cylinder_Half_Length + 3*cm;
    const G4double     Polythene_Cylinder_zAxis_offset = -38.0*mm;
    


    //----PE Collimator ---------
    
    G4Tubs* Collimator_tube = new G4Tubs("Collimator_tube", Polythene_Cylinder_Inner_Radius,  Polythene_Cylinder_Outer_Radius,  Polythene_Cylinder_Half_Length,  0.*deg,  360.*deg);

    G4LogicalVolume*  Collimator_tube_Logical = new G4LogicalVolume(Collimator_tube, polyethylene, "Collimator_tube_Logical");

    G4RotationMatrix R_collimator;
    G4ThreeVector    T_collimator;
    
    R_collimator.rotateY(90.0*deg);
    T_collimator.setX(Polythene_Cylinder_xAxis_Distance); T_collimator.setY(0); T_collimator.setZ(Polythene_Cylinder_zAxis_offset);
    
    //new G4PVPlacement(G4Transform3D(R_collimator, T_collimator),  Collimator_tube_Logical, "Collimator_tube_Physical", fpWorldLogical, false, 0);

    
    
    
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
    
    
    
    
        //===============================================

	// Invisible world volume.
	fpWorldLogical->SetVisAttributes(G4VisAttributes::Invisible);
	

    GasArgon_Logical->SetVisAttributes(G4Colour::Yellow());
	FullLiquidArgon_Logical->SetVisAttributes(G4Colour::Yellow());
    LiquidArgon_Logical->SetVisAttributes(G4Colour::Yellow());

   
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
