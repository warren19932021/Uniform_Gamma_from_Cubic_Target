//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file MuStoppedDetectorConstruction.cc
/// \brief Implementation of the MuStoppedDetectorConstruction class

#include "MuStoppedDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "MuStoppedRunAction.hh"
#include "G4VisAttributes.hh"


#include "MuStoppedTrackingAction.hh" // Include the full definition


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MuStoppedDetectorConstruction::MuStoppedDetectorConstruction()
: G4VUserDetectorConstruction(),
  fTheLogicalVolumeWithEnergyRecord(0),
  BlockLength(70),
  BlockWidth(70),
  BlockHeight(70),

  TargetMaterial("G4_Ti")
{


   fMesesenger = new G4GenericMessenger(this,"/detector/","Detector Construction");   
   
   fMesesenger->DeclareProperty("BlockLength",BlockLength,"Side Thcikness of CupTarget");
   fMesesenger->DeclareProperty("BlockWidth",BlockWidth,"Front Thcikness of CupTarget");
   fMesesenger->DeclareProperty("BlockHeight",BlockHeight,"Height of CupTarget");


   fMesesenger->DeclareProperty("TargetMaterial",TargetMaterial,"Target Material");

    DefineMaterials();
    DefineMaterialsFromAndi();
    isAndyGeDet = true; 
    if(isAndyGeDet) 
    {
      andiDet = new GeometryGeDet(GeometryGeDet::Iktp_nType, false);
      //andiDet->SetDefaults();
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MuStoppedDetectorConstruction::~MuStoppedDetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MuStoppedDetectorConstruction::SetRunAction(MuStoppedRunAction* pRunAction)
{
    if (!pRunAction) G4cout << "pRunAction is nullptr" << G4endl;
    runAction = pRunAction;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MuStoppedDetectorConstruction::SetTrackingAction(MuStoppedTrackingAction* pTrackingAction)
{
    if (!pTrackingAction) G4cout << "pTrackingAction is nullptr" << G4endl;
    trackingAction = pTrackingAction;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MuStoppedDetectorConstruction::DefineMaterials()
{

    G4NistManager * nist = G4NistManager::Instance();
    worldMat = nist->FindOrBuildMaterial("G4_AIR");
/*    
    G4Isotope* isotopeTe124 = new G4Isotope("Te124",52,124,123.904*g / mole);   
    G4Element* elementTe124 = new G4Element("Tellurium124","Te124",1);
    elementTe124->AddIsotope(isotopeTe124,1);   
    Target_material = new G4Material("Target_material",6.24*g/cm3,1);
    Target_material->AddElement(elementTe124,1.0);
*/

    Target_material = nist->FindOrBuildMaterial(TargetMaterial);
    
    G4cout<<"TargetMaterial: ============================    "<<Target_material->GetName()<<G4endl;

}    


void MuStoppedDetectorConstruction::DefineMaterialsFromAndi(){
  // get NIST-Material-Manager
  G4NistManager* nist = G4NistManager::Instance();


  // defining Materials
  // ------------------------------  
  //here one should define the material he/she wants to use 
  //in the geometry, if not read from the DB
  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density, fractionmass, abundance;
  G4String name, symbol;
  G4int ncomponents, natoms;


  // >>>>>>>> Elements
  G4Element* elH = nist->FindOrBuildElement("H");
  G4Element* elO = nist->FindOrBuildElement("O");
  G4Element* elC = nist->FindOrBuildElement("C");
  G4Element* elB = nist->FindOrBuildElement("B");
  G4Element* elN = nist->FindOrBuildElement("N");
  G4Element* elAr = nist->FindOrBuildElement("Ar");
  G4Element* elSi = nist->FindOrBuildElement("Si");
  G4Element* elCr = nist->FindOrBuildElement("Cr");
  G4Element* elMn = nist->FindOrBuildElement("Mn");
  G4Element* elFe = nist->FindOrBuildElement("Fe");
  G4Element* elNi = nist->FindOrBuildElement("Ni");
  G4Element* elLi = nist->FindOrBuildElement("Li");
  G4Element* elPb = nist->FindOrBuildElement("Pb");
  G4Element* elAl = nist->FindOrBuildElement("Al");
  G4Element* elMg = nist->FindOrBuildElement("Mg");
  G4Element* elCu = nist->FindOrBuildElement("Cu");

  G4Isotope* Ge70 = new G4Isotope(name = "Ge70", 32, 70, 69.92*g / mole);
  G4Isotope* Ge72 = new G4Isotope(name = "Ge72", 32, 72, 71.92*g / mole);
  G4Isotope* Ge73 = new G4Isotope(name = "Ge73", 32, 73, 73.0*g / mole);
  G4Isotope* Ge74 = new G4Isotope(name = "Ge74", 32, 74, 74.0*g / mole);
  G4Isotope* Ge76 = new G4Isotope(name = "Ge76", 32, 76, 76.0*g / mole);

  a = 75.71*g / mole;
  G4int nIsotopes = 5;
  G4Element* elGeNat = new G4Element(name = "naturalGermanium", symbol = "GeNat", nIsotopes);
  elGeNat->AddIsotope(Ge70, abundance = 20.9*perCent);
  elGeNat->AddIsotope(Ge72, abundance = 27.5*perCent);
  elGeNat->AddIsotope(Ge73, abundance = 7.7*perCent);
  elGeNat->AddIsotope(Ge74, abundance = 36.3*perCent);
  elGeNat->AddIsotope(Ge76, abundance = 7.6*perCent);

  G4Element* elGeEnr = new G4Element(name = "enrichedGermanium", symbol = "GeEnr", nIsotopes);
  elGeEnr->AddIsotope(Ge70, abundance = 0.0*perCent);
  elGeEnr->AddIsotope(Ge72, abundance = 0.1*perCent);
  elGeEnr->AddIsotope(Ge73, abundance = 0.2*perCent);
  elGeEnr->AddIsotope(Ge74, abundance = 13.1*perCent);
  elGeEnr->AddIsotope(Ge76, abundance = 86.6*perCent);


  // >>>>>>>> Materials
  // StainlessSteel
  G4Material* StainlessSteel = new G4Material("StainlessSteel", density = 8.06*g / cm3, ncomponents = 6);
  StainlessSteel->AddElement(elC, fractionmass = 0.001);
  StainlessSteel->AddElement(elSi, fractionmass = 0.007);
  StainlessSteel->AddElement(elCr, fractionmass = 0.18);
  StainlessSteel->AddElement(elMn, fractionmass = 0.01);
  StainlessSteel->AddElement(elFe, fractionmass = 0.712);
  StainlessSteel->AddElement(elNi, fractionmass = 0.09);

  // Vacuum
  a = 14.000674*g / mole;
  density = 1.E-10*g / cm3;
  G4Material* Vacuum = new G4Material(name = "Vacuum", z = 7., a, density);

  // Natural Germanium
  density = 5.32*g / cm3;
  G4Material* NaturalGe = new G4Material(name = "NaturalGe", density, 1);
  NaturalGe->AddElement(elGeNat, natoms = 1);

  // Enriched Germanium
  density = 5.54*g / cm3;
  G4Material* EnrichedGe = new G4Material(name = "EnrichedGe", density, 1);
  EnrichedGe->AddElement(elGeEnr, natoms = 1);

  // Germanium/Lithium mixture (1%, det-dead layer)
  density = 5.323*g / cm3;
  G4Material* GeLi = new G4Material(name = "GeLi", density, ncomponents = 2);
  GeLi->AddElement(elGeNat, fractionmass = 0.99);
  GeLi->AddElement(elLi, fractionmass = 0.01);

  // Germanium/Bor mixture (1%, det-dead layer)
  density = 5.323*g / cm3;
  G4Material* GeB = new G4Material(name = "GeB", density, ncomponents = 2);
  GeB->AddElement(elGeNat, fractionmass = 0.99);
  GeB->AddElement(elB, fractionmass = 0.01);

  // PolyEthylen Terephthalat  (PET, Mylar)
  G4Material* PET = new G4Material(name = "PET", density = 1370 * kg / m3, ncomponents = 3);
  PET->AddElement(elC, natoms = 10);
  PET->AddElement(elH, natoms = 8);
  PET->AddElement(elO, natoms = 4);

  // Polylactide  ---- TODO: density!!!
  G4Material* PLA = new G4Material(name = "PLA", density = 1400 * kg / m3, ncomponents = 3);
  PLA->AddElement(elC, natoms = 3);
  PLA->AddElement(elH, natoms = 4);
  PLA->AddElement(elO, natoms = 2);

  // Nylon  ---- TODO: Density, better source (wikipedia right now)
  // PA 6 : [NH−(CH2)5−CO]n || AUCH MOEGLICH: G4_NYLON-6/6   (ist das Selbe T_T)
  G4Material* Nylon = new G4Material(name = "Nylon", density = 1150 * kg / m3, ncomponents = 4);
  Nylon->AddElement(elC, natoms = 6);
  Nylon->AddElement(elH, natoms = 11);
  Nylon->AddElement(elO, natoms = 1);
  Nylon->AddElement(elN, natoms = 1);

  // Styropor
  G4Material* Styropor = new G4Material(name = "Styropor", density = 20 * kg / m3, ncomponents = 2);
  Styropor->AddElement(elH, fractionmass = 12.8*perCent);
  Styropor->AddElement(elC, fractionmass = 87.2*perCent);

  // Polycarbonate
  G4Material* Polycarbonate = new G4Material(name = "Polycarbonate", density = 1200 * kg / m3, ncomponents = 3);
  Polycarbonate->AddElement(elH, fractionmass = 5.5491*perCent);
  Polycarbonate->AddElement(elC, fractionmass = 75.5751*perCent);
  Polycarbonate->AddElement(elO, fractionmass = 18.8758*perCent);
  
  // CRFP (Carbon Fiber Reinforced Polymer): M55 Quasiisotropic Layup  ---------------------------- TODO?!
  G4Material* CFRP = new G4Material("CFRP", 1.66*g / cm3, 1);
  CFRP->AddElement(elC, 1);

  // Impure Aluminium (test with 0,5% Pb as impurity)
  density = 2.7*g / cm3;
  G4Material* ImpureAl = new G4Material(name = "ImpureAl", density, ncomponents = 2);
  // FractionMass_A = Abundance_A * AtomicMass_A / ( Abundance_A * AtomicMass_A + Abundance_B * AtomicMass_B)
  // std::cout << "Atomasse Al.: " << elAl->GetAtomicMassAmu() << " | Atommasse Pb: " << elPb->GetAtomicMassAmu() << std::endl;
  ImpureAl->AddElement(elAl, fractionmass = 0.96284);
  ImpureAl->AddElement(elPb, fractionmass = 0.03716);

  // Impure Aluminium2 - See Norm (test with 0,5% Pb as impurity)
  density = 2.7*g / cm3;
  G4Material* ImpureAl2 = new G4Material(name = "ImpureAl2", density, ncomponents = 7);
  // FractionMass_A = Abundance_A * AtomicMass_A / ( Abundance_A * AtomicMass_A + Abundance_B * AtomicMass_B)
  // std::cout << "Atomasse Al.: " << elAl->GetAtomicMassAmu() << " | Atommasse Pb: " << elPb->GetAtomicMassAmu() << std::endl;
  ImpureAl2->AddElement(elAl, fractionmass = 0.99);
  ImpureAl2->AddElement(elSi, fractionmass = 0.0027);
  ImpureAl2->AddElement(elFe, fractionmass = 0.0027);
  ImpureAl2->AddElement(elCu, fractionmass = 0.0006);
  ImpureAl2->AddElement(elCr, fractionmass = 0.0006);
  ImpureAl2->AddElement(elMn, fractionmass = 0.0017);
  ImpureAl2->AddElement(elMg, fractionmass = 0.0017); 

  // Customized Air
  // Option 1: Calculate Density manually with general Gas-Equation; (air_temp and air_pressure are going to be iggnored in constructor)
  // use customized Air-Components etc!
  G4double air_temp = 293.0*kelvin;
  G4double air_pressure = 1.0205 * bar;
  density = (air_pressure / pascal) * 1.0 / ((air_temp / kelvin) * 287.058) * kg / m3;
  G4Material* Custom_Air_01 = new G4Material("Custom_Air_01", density, ncomponents = 4, kStateGas, air_temp, air_pressure);
  Custom_Air_01->AddElement(elN, fractionmass = 0.781);
  Custom_Air_01->AddElement(elO, fractionmass = 0.209);
  Custom_Air_01->AddElement(elAr, fractionmass = 0.0093);
  Custom_Air_01->AddElement(elC, fractionmass = 0.0007);

  // Option 2: G4-Intern calcualtion of density; only give pressure and temperature!!
  // use Standard-G4-Air Components
  G4Material* Custom_Air_02 = nist->ConstructNewGasMaterial("Custom_Air_02", "G4_AIR", air_temp, air_pressure);

  // Nist-Basic Materials
  G4Material* Copper = nist->FindOrBuildMaterial("G4_Cu");
  G4Material* Beryllium = nist->FindOrBuildMaterial("G4_Be");
  G4Material* Air = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* Helium = nist->FindOrBuildMaterial("G4_He");
  G4Material* Aluminium = nist->FindOrBuildMaterial("G4_Al");
  G4Material* Tungsten = nist->FindOrBuildMaterial("G4_W");
  G4Material* Tin = nist->FindOrBuildMaterial("G4_Sn");

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;  
}


void MuStoppedDetectorConstruction:: ConstructAndyGeDet()
{

//--------------------------------------  construct Andy detector
// the surface is 1.23 mm away from the durface of block target position (0,0,0)
   andiDet->CreateGeometry(logicalWorld);
   logicalAndyDetector = andiDet->GetActiveLogical();
   G4cout<<"logicalAndyDetector->GetName(): "<<logicalAndyDetector->GetName()<<G4endl;
}


G4VPhysicalVolume* MuStoppedDetectorConstruction::Construct()
{  
    
    //---- World------
    solidWorld =  new G4Box("SolidWorld",50.0*m,50.0*m,50.0*m); 
    logicalWorld = new G4LogicalVolume(solidWorld, worldMat,"logicWorld");
    physWorld = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),logicalWorld,"physWolrd",0,false,0,true); 

    G4VisAttributes* worldVisAttributes = new G4VisAttributes(false); // tranparent
    logicalWorld->SetVisAttributes(worldVisAttributes);


   //---- solidTarget------
    solidTarget =  new G4Box("SolidTarget",BlockLength/2,BlockWidth/2,BlockHeight/2); 
    logical_block_target = new G4LogicalVolume(solidTarget,Target_material,"logical_block_target");
    phys_block_Target = new G4PVPlacement(0,G4ThreeVector(0.,0.,BlockHeight/2),logical_block_target,"physTarget",logicalWorld,false,0,true); 

    G4cout<<"oooooo0000000000------DEBUG::Detector Setup Output------oooooo0000000000"<<G4endl;
    G4cout<<"BlockLength: "<<BlockLength<<G4endl;
    G4cout<<"BlockWidth: "<<BlockWidth<<G4endl;
    G4cout<<"BlockHeight: "<<BlockHeight<<G4endl;
    G4cout<<"TargetMaterial: "<<TargetMaterial<<G4endl;
    G4cout<<"oooooo0000000000---------------------------------------ooooooo0000000000"<<G4endl;


   //---GeDet----
   if(isAndyGeDet)
      ConstructAndyGeDet();



    // Set Block shape Target as TheLogicalVolumeWithEnergyRecord volume
    //fTheLogicalVolumeWithEnergyRecord = logical_block_target;
    fTheLogicalVolumeWithEnergyRecord = logicalAndyDetector;
   
    G4cout<<"fTheLogicalVolumeWithEnergyRecord->GetName() in MuStoppedDetectorConstruction: "<<fTheLogicalVolumeWithEnergyRecord->GetName()<<G4endl;

    //always return the physical World
    return physWorld;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MuStoppedDetectorConstruction::ConstructSDandField()
{
    
   // senDet = new MySensitiveDetector("Test",runAction,trackingAction);
   // logicalAndyDetector->SetSensitiveDetector(senDet);


    // Create the sensitive detector
    senDet = new MySensitiveDetector("Test", runAction, trackingAction );
    logicalAndyDetector->SetSensitiveDetector(senDet);

    // Set the sensitive detector in the tracking action
    trackingAction->SetSensitiveDetector(senDet);

}    