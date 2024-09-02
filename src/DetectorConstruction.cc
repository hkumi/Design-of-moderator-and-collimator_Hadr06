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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fMaterial(0), fLAbsor(0),  fWorldMat(0), fPWorld(0), fDetectorMessenger(0)
{
  fRadius = 3*mm;
  fWorldSize = 3*m;
  DefineMaterials();
  SetMaterial("CF4");  
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // specific element name for thermal neutronHP
  // (see G4ParticleHPThermalScatteringNames.cc)
  G4NistManager *nist = G4NistManager::Instance();
  G4int ncomponents, natoms;
   // boron
  G4Isotope* B10 = new G4Isotope("B10", 5, 10);
  G4Isotope* B11 = new G4Isotope("B11", 5, 11);
  G4Element* B = new G4Element("Boron", "B", ncomponents=2);
  B->AddIsotope(B10, 19.9*perCent);
  B->AddIsotope(B11, 80.1*perCent);
  G4Material* boron = new G4Material("boron", 2.46*g/cm3, ncomponents=1, kStateSolid,293*kelvin, 1*atmosphere);
  boron->AddElement(B, natoms=1);

  // pressurized water
  G4Element* H  = new G4Element("TS_H_of_Water" ,"H" , 1., 1.0079*g/mole);
  G4Element* O  = new G4Element("Oxygen"        ,"O" , 8., 16.00*g/mole);
  G4Material* H2O = 
  new G4Material("Water_ts", 1.000*g/cm3, ncomponents=2,
                            kStateLiquid, 593*kelvin, 150*bar);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);

  // heavy water
  G4Isotope* H2 = new G4Isotope("H2",1,2);
  G4Element* D  = new G4Element("TS_D_of_Heavy_Water", "D", 1);
  D->AddIsotope(H2, 100*perCent);  
  G4Material* D2O = new G4Material("HeavyWater", 1.11*g/cm3, ncomponents=2,
                        kStateLiquid, 293.15*kelvin, 1*atmosphere);
  D2O->AddElement(D, natoms=2);
  D2O->AddElement(O, natoms=1);

  // graphite
  G4Isotope* C12 = new G4Isotope("C12", 6, 12);  
  G4Element* C   = new G4Element("TS_C_of_Graphite","C", ncomponents=1);
  C->AddIsotope(C12, 100.*perCent);
  G4Material* graphite = 
  new G4Material("graphite", 2.27*g/cm3, ncomponents=1,
                         kStateSolid, 293*kelvin, 1*atmosphere);
  graphite->AddElement(C, natoms=1);
  
  //NE213
  G4Material* ne213 = 
  new G4Material("NE213", 0.874*g/cm3, ncomponents=2);
  ne213->AddElement(H,  9.2*perCent);
  ne213->AddElement(C, 90.8*perCent);  
      
  // example of vacuum
  fWorldMat = new G4Material("Galactic", 1, 1.01*g/mole, universe_mean_density,
                 kStateGas, 2.73*kelvin, 3.e-18*pascal);

 ///G4cout << *(G4Material::GetMaterialTable()) << G4endl;

 // polyethilene
  G4Element* Hpe = new G4Element("TS_H_of_Polyethylene", "H", 1, 1.0079*g/mole);
  G4Element* Cpe = new G4Element("Carbon", "C", 6, 12.01*g/mole);
  polyethylene = new G4Material("polyethylene", 0.93*g/cm3, ncomponents=2, kStateSolid, 293*kelvin, 1*atmosphere);
  polyethylene->AddElement(Hpe, natoms=4);
  polyethylene->AddElement(Cpe, natoms=2);

  G4double Vdens = 1.e-25*g/cm3;
  G4double Vpres = 1.e-19*pascal;
  G4double Vtemp = 0.1*kelvin;
  G4double a, z;

   // Define the lead material
  leadMaterial = new G4Material("Lead", 82, 207.2 * g/mole, 11.35 * g/cm3);

  // vacuum
  Vacc = new G4Material("Galactic", z=1, a=1.01*g/mole, Vdens, kStateGas, Vtemp, Vpres);

  // borated polyethilene
  b_polyethylene = new G4Material("b_polyethylene",0.94*g/cm3,ncomponents=4,kStateSolid,293*kelvin,1*atmosphere);
  b_polyethylene->AddElement(Hpe, 11.6*perCent);
  b_polyethylene->AddElement(Cpe, 61.2*perCent);
  b_polyethylene->AddElement(B, 5*perCent);
  b_polyethylene->AddElement(O, 22.2*perCent);

   // Define the lead material
  leadMaterial = new G4Material("Lead", 82, 207.2 * g/mole, 11.35 * g/cm3);

  //...............creating the materials for the scintillator..............................
  //----------------------------------- CarbonTetrafluoride ------------------------
  G4Element *F  = new G4Element("Fluorine","F",9.,18.998*g/mole);
  G4double pressure = 0.0328947*atmosphere; //25Torr
  G4double temperature = 293.15*kelvin; // 
  CF4 = new G4Material("CF4", 0.1223*mg/cm3,2,kStateGas,temperature,pressure);
  CF4->AddElement(C,1);
  CF4->AddElement(F,4);

  // ------------------------------------ Polypropilene ------------------------------------

  nist->FindOrBuildMaterial("G4_POLYPROPYLENE");
  PP = G4Material::GetMaterial("G4_POLYPROPYLENE");

  
//----------------------------------- Aluminium ------------------------------------
 
   G4double density;
   Aluminium = new G4Material("Al", z = 13., a = 26.98 * g / mole,
                        density = 2.7 * g / cm3);
 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::MaterialWithSingleIsotope( G4String name,
                           G4String symbol, G4double density, G4int Z, G4int A)
{
 // define a material from an isotope
 //
 G4int ncomponents;
 G4double abundance, massfraction;

 G4Isotope* isotope = new G4Isotope(symbol, Z, A);
 
 G4Element* element  = new G4Element(name, symbol, ncomponents=1);
 element->AddIsotope(isotope, abundance= 100.*perCent);
 
 G4Material* material = new G4Material(name, density, ncomponents=1);
 material->AddElement(element, massfraction=100.*perCent);

 return material;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  // World
  //
  G4Box*
  sWorld = new G4Box("World",                           //name
                   3*m/2,3*m/2,3*m/2);   //dimensions
                   
  G4LogicalVolume*
  lWorld = new G4LogicalVolume(sWorld,                  //shape
                             fWorldMat,                 //material
                             "World");                  //name

  fPWorld = new G4PVPlacement(0,                        //no rotation
                            G4ThreeVector(),            //at (0,0,0)
                            lWorld,                     //logical volume
                            "World",                    //name
                            0,                          //mother volume
                            false,                      //no boolean operation
                            0);                         //copy number
 

  //The HDPE_block1

  auto fblockSize = 10*cm;


  auto HDPE_Box1 = new G4Box("HDPE1",                             //its name
                   10*cm/2,fblockSize/2,10*cm/2);   //its dimensions

  auto HDPE_LV1 = new G4LogicalVolume(HDPE_Box1,                     //its shape
                              polyethylene,                      //its material
                             "HDPE1");                  //its name

  auto HDPE_PV1 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0*cm,0*cm,5*cm),            //at (0,0,0)
                             HDPE_LV1,                      //its logical volume
                            "HDPE1",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number


  //The HDPE_block2


  auto HDPE_Box2 = new G4Box("HDPE2",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  auto HDPE_LV2 = new G4LogicalVolume(HDPE_Box2,                     //its shape
                              polyethylene,                      //its material
                             "HDPE2");                  //its name

  auto HDPE_PV2 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,10*cm,5*cm),            //at (0,0,0)
                             HDPE_LV2,                      //its logical volume
                            "HDPE2",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

 //The HDPE_block3


  auto HDPE_Box3 = new G4Box("HDPE3",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  auto HDPE_LV3 = new G4LogicalVolume(HDPE_Box3,                     //its shape
                              polyethylene,                      //its material
                             "HDPE3");                  //its name

  auto HDPE_PV3 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,-10*cm,5*cm),            //at (0,0,0)
                             HDPE_LV3,                      //its logical volume
                            "HDPE3",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number
          
  
  //The HDPE_block4


  auto HDPE_Box4 = new G4Box("HDPE4",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  auto HDPE_LV4 = new G4LogicalVolume(HDPE_Box4,                     //its shape
                              polyethylene,                      //its material
                             "HDPE4");                  //its name

  auto HDPE_PV4 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,10*cm,15*cm),            //at (0,0,0)
                             HDPE_LV4,                      //its logical volume
                            "HDPE4",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

//The HDPE_block5


  auto HDPE_Box5 = new G4Box("HDPE5",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  auto HDPE_LV5 = new G4LogicalVolume(HDPE_Box5,                     //its shape
                              polyethylene,                      //its material
                             "HDPE5");                  //its name

  auto HDPE_PV5 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,-10*cm,15*cm),            //at (0,0,0)
                             HDPE_LV5,                      //its logical volume
                            "HDPE5",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number
  //The HDPE_block6


  auto HDPE_Box6 = new G4Box("HDPE6",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  auto HDPE_LV6 = new G4LogicalVolume(HDPE_Box6,                     //its shape
                              polyethylene,                      //its material
                             "HDPE6");                  //its name

  auto HDPE_PV6 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,0,5*cm),            //at (0,0,0)
                             HDPE_LV6,                      //its logical volume
                            "HDPE6",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

  //The HDPE_block7


  auto HDPE_Box7 = new G4Box("HDPE7",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  auto HDPE_LV7 = new G4LogicalVolume(HDPE_Box7,                     //its shape
                              polyethylene,                      //its material
                             "HDPE7");                  //its name

  auto HDPE_PV7 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,0,5*cm),            //at (0,0,0)
                             HDPE_LV7,                      //its logical volume
                            "HDPE7",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number
  //The HDPE_block8


  auto HDPE_Box8 = new G4Box("HDPE8",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  auto HDPE_LV8 = new G4LogicalVolume(HDPE_Box8,                     //its shape
                              polyethylene,                      //its material
                             "HDPE8");                  //its name

  auto HDPE_PV8 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,10*cm,5*cm),            //at (0,0,0)
                             HDPE_LV8,                      //its logical volume
                            "HDPE8",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number


     //The HDPE_block9


  auto HDPE_Box9 = new G4Box("HDPE9",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  auto HDPE_LV9 = new G4LogicalVolume(HDPE_Box9,                     //its shape
                              polyethylene,                      //its material
                             "HDPE9");                  //its name

  auto HDPE_PV9 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,-10*cm,5*cm),            //at (0,0,0)
                             HDPE_LV9,                      //its logical volume
                            "HDPE9",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

  //The HDPE_block10

  auto HDPE_Box10 = new G4Box("HDPE10",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  auto HDPE_LV10 = new G4LogicalVolume(HDPE_Box10,                     //its shape
                              polyethylene,                      //its material
                             "HDPE10");                  //its name

  auto HDPE_PV10 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,10*cm,5*cm),            //at (0,0,0)
                             HDPE_LV10,                      //its logical volume
                            "HDPE10",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number 

  //The HDPE_block11

  auto HDPE_Box11 = new G4Box("HDPE11",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  auto HDPE_LV11 = new G4LogicalVolume(HDPE_Box11,                     //its shape
                              polyethylene,                      //its material
                             "HDPE11");                  //its name

  auto HDPE_PV11 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,-10*cm,5*cm),            //at (0,0,0)
                             HDPE_LV11,                      //its logical volume
                            "HDPE11",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number


//The HDPE_block12


  auto HDPE_Box12 = new G4Box("HDPE12",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  auto HDPE_LV12 = new G4LogicalVolume(HDPE_Box12,                     //its shape
                              polyethylene,                      //its material
                             "HDPE12");                  //its name

  auto HDPE_PV12 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,-10*cm,15*cm),            //at (0,0,0)
                             HDPE_LV12,                      //its logical volume
                            "HDPE12",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number
 //The HDPE_block13


  auto HDPE_Box13 = new G4Box("HDPE13",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  auto HDPE_LV13 = new G4LogicalVolume(HDPE_Box13,                     //its shape
                              polyethylene,                      //its material
                             "HDPE13");                  //its name

  auto HDPE_PV13 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,10*cm,15*cm),            //at (0,0,0)
                             HDPE_LV13,                      //its logical volume
                            "HDPE13",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

  //The HDPE_block14


  auto HDPE_Box14 = new G4Box("HDPE14",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  auto HDPE_LV14 = new G4LogicalVolume(HDPE_Box14,                     //its shape
                              polyethylene,                      //its material
                             "HDPE14");                  //its name

  auto HDPE_PV14 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,-10*cm,15*cm),            //at (0,0,0)
                             HDPE_LV14,                      //its logical volume
                            "HDPE14",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number



  //The HDPE_block15

  auto HDPE_Box15 = new G4Box("HDPE15",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  auto HDPE_LV15 = new G4LogicalVolume(HDPE_Box15,                     //its shape
                              polyethylene,                      //its material
                             "HDPE15");                  //its name

  auto HDPE_PV15 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,10*cm,15*cm),            //at (0,0,0)
                             HDPE_LV15,                      //its logical volume
                            "HDPE15",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

//The HDPE_block16


  auto HDPE_Box16 = new G4Box("HDPE16",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  auto HDPE_LV16 = new G4LogicalVolume(HDPE_Box16,                     //its shape
                              polyethylene,                      //its material
                             "HDPE16");                  //its name

  auto HDPE_PV16 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,0*cm,15*cm),            //at (0,0,0)
                             HDPE_LV16,                      //its logical volume
                            "HDPE16",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number
   //The HDPE_block17

  auto HDPE_Box17 = new G4Box("HDPE17",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  auto HDPE_LV17 = new G4LogicalVolume(HDPE_Box17,                     //its shape
                              polyethylene,                      //its material
                             "HDPE17");                  //its name

  auto HDPE_PV17 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,0*cm,15*cm),            //at (0,0,0)
                             HDPE_LV17,                      //its logical volume
                            "HDPE17",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

    //The HDPE_block18


  auto HDPE_Box18 = new G4Box("HDPE18",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  auto HDPE_LV18 = new G4LogicalVolume(HDPE_Box18,                     //its shape
                              polyethylene,                      //its material
                             "HDPE18");                  //its name

  auto HDPE_PV18 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,-10*cm,25*cm),            //at (0,0,0)
                             HDPE_LV18,                      //its logical volume
                            "HDPE18",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

    //The HDPE_block19


  auto HDPE_Box19 = new G4Box("HDPE19",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  auto HDPE_LV19 = new G4LogicalVolume(HDPE_Box19,                     //its shape
                              polyethylene,                      //its material
                             "HDPE19");                  //its name

  auto HDPE_PV19 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,-10*cm,25*cm),            //at (0,0,0)
                             HDPE_LV19,                      //its logical volume
                            "HDPE19",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);
    //The HDPE_block20


  auto HDPE_Box20 = new G4Box("HDPE20",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  auto HDPE_LV20 = new G4LogicalVolume(HDPE_Box20,                     //its shape
                              polyethylene,                      //its material
                             "HDPE20");                  //its name

  auto HDPE_PV20 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,10*cm,25*cm),            //at (0,0,0)
                             HDPE_LV20,                      //its logical volume
                            "HDPE20",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

 //The HDPE_block21

  auto HDPE_Box21 = new G4Box("HDPE21",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  auto HDPE_LV21 = new G4LogicalVolume(HDPE_Box21,                     //its shape
                              polyethylene,                      //its material
                             "HDPE21");                  //its name

  auto HDPE_PV21 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,10*cm,25*cm),            //at (0,0,0)
                             HDPE_LV21,                      //its logical volume
                            "HDPE21",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

//The HDPE_block22


  auto HDPE_Box22 = new G4Box("HDPE22",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  auto HDPE_LV22 = new G4LogicalVolume(HDPE_Box22,                     //its shape
                              polyethylene,                      //its material
                             "HDPE22");                  //its name

  auto HDPE_PV22 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,0*cm,25*cm),            //at (0,0,0)
                             HDPE_LV22,                      //its logical volume
                            "HDPE22",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

  //The HDPE_block23

  auto HDPE_Box23 = new G4Box("HDPE23",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  auto HDPE_LV23 = new G4LogicalVolume(HDPE_Box23,                     //its shape
                              polyethylene,                      //its material
                             "HDPE23");                  //its name

  auto HDPE_PV23 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,0*cm,25*cm),            //at (0,0,0)
                             HDPE_LV23,                      //its logical volume
                            "HDPE23",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number


 //The HDPE_block24
  

  auto HDPE_Box24 = new G4Box("HDPE24",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  auto HDPE_LV24 = new G4LogicalVolume(HDPE_Box24,                     //its shape
                              polyethylene,                      //its material
                             "HDPE24");                  //its name

  auto HDPE_PV24 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,10*cm,25*cm),            //at (0,0,0)
                             HDPE_LV24,                      //its logical volume
                            "HDPE24",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

//The HDPE_block25


  auto HDPE_Box25 = new G4Box("HDPE25",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  auto HDPE_LV25 = new G4LogicalVolume(HDPE_Box25,                     //its shape
                              polyethylene,                      //its material
                             "HDPE25");                  //its name

  auto HDPE_PV25 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,-10*cm,25*cm),            //at (0,0,0)
                             HDPE_LV25,                      //its logical volume
                            "HDPE25",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number 
  //The lead1
  auto fLeadSize = 10*cm;


  auto Lead_Box = new G4Box("Lead",                             //its name
                   fLeadSize/2,fLeadSize/2, 5*cm/2);   //its dimensions

  auto Lead_LV = new G4LogicalVolume(Lead_Box,                     //its shape
                              leadMaterial,                      //its material
                             "Lead");                  //its name

  auto Lead_PV = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,0,12.5*cm),            //at (0,0,0)
                             Lead_LV,                      //its logical volume
                            "Lead",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number


   G4VisAttributes* red = new G4VisAttributes(G4Colour::Red());

   red->SetVisibility(true);
   red->SetForceAuxEdgeVisible(true);

   Lead_LV->SetVisAttributes(red);
 
    //The Borated polythylene_block1 with pinhole

  auto BoratedSize = 30*cm;
  auto Borated_thickness = 3*cm;
  auto Borated_Box1 = new G4Box("Borated1",                             //its name
                   BoratedSize/2,  BoratedSize/2,Borated_thickness/2);   //its dimensions



  auto Hole = new G4Tubs("BoxHole", 0.0*cm, 1.5*cm, 1.5*cm, 0*deg, 360*deg);  // the diameter of the exit(pinhole) is 3cm. In G4 we use halfsize of the radius.
  auto Hole_1 = new G4Tubs("BoxHole", 0.0*cm, 1.5*cm, 1.5*cm, 0*deg, 360*deg);  // the diameter of the exit(pinhole) is 3cm. In G4 we use halfsize of the radius. 
  auto Hole_LV1 = new G4LogicalVolume(Hole_1,                     //its shape
                              Vacc,                      //its material
                             "H1");                  //its name

 

  auto Hole_LV = new G4LogicalVolume(Hole,                     //its shape
                              Vacc,                      //its material
                             "H1");                  //its name

  auto Borated_LV1 = new G4LogicalVolume(Borated_Box1,                     //its shape
                              b_polyethylene,                      //its material
                             "Borated1", 0,0,0);                  //its name



  auto Borated_PV1 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0*cm,0*cm,31.5*cm),            //at (0,0,0)
                             Borated_LV1,                      //its logical volume
                            "Borated1",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

  auto Hole_PV = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0*cm,0*cm,0*cm),            //at (0,0,0)
                             Hole_LV,                      //its logical volume
                            "H1",                    //its name

                            Borated_LV1,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number
   G4VisAttributes* green = new G4VisAttributes(G4Colour::Green());

   green->SetVisibility(true);
   green->SetForceAuxEdgeVisible(true);

   Borated_LV1->SetVisAttributes(green);

   //The lead2
  

  auto LeadSize = 3*cm;
  auto Lead_Box2 = new G4Box("Lead2",                             //its name
                   30*cm/2,LeadSize/2,33*cm/2);   //its dimensions

  auto Lead_LV2 = new G4LogicalVolume(Lead_Box2,                     //its shape
                              leadMaterial,                      //its material
                             "Lead2");                  //its name

  auto Lead_PV2 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0*cm,16.5*cm,16.5*cm),            //at (0,0,0)
                             Lead_LV2,                      //its logical volume
                            "Lead2",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number


  G4VisAttributes* yellow= new G4VisAttributes(G4Colour::Yellow());

  yellow->SetVisibility(true);
  yellow->SetForceAuxEdgeVisible(true);

  Lead_LV2->SetVisAttributes(red);
 

   //The lead3
 
  auto Lead_Box3 = new G4Box("Lead3",                             //its name
                   3*cm/2,30*cm/2, 33*cm/2);   //its dimensions

  auto Lead_LV3 = new G4LogicalVolume(Lead_Box3,                     //its shape
                              leadMaterial,                      //its material
                             "Lead3");                  //its name

  auto Lead_PV3 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(16.5*cm,0*cm,16.5*cm),            //at (0,0,0)
                             Lead_LV3,                      //its logical volume
                            "Lead3",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number


   
  Lead_LV3->SetVisAttributes(red);


   //The lead4

  auto Lead_Box4 = new G4Box("Lead4",                             //its name
                   3*cm/2,30*cm/2, 33*cm/2);   //its dimensions

  auto Lead_LV4 = new G4LogicalVolume(Lead_Box4,                     //its shape
                              leadMaterial,                      //its material
                             "Lead4");                  //its name

  auto Lead_PV4 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-16.5*cm,0*cm,16.5*cm),            //at (0,0,0)
                             Lead_LV4,                      //its logical volume
                            "Lead4",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number


   
  Lead_LV4->SetVisAttributes(red);
  

  //The lead5

  auto Lead_Box5 = new G4Box("Lead5",                             //its name
                   30*cm/2,30*cm/2, 3*cm/2);   //its dimensions

  auto Lead_LV5 = new G4LogicalVolume(Lead_Box5,                     //its shape
                              leadMaterial,                      //its material
                             "Lead5");                  //its name

  auto Lead_PV5 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0*cm,0*cm,34.5*cm),            //at (0,0,0)
                             Lead_LV5,                      //its logical volume
                            "Lead5",                    //its name
                            lWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

  auto Hole_PV2 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0*cm,0*cm,0*cm),            //at (0,0,0)
                             Hole_LV1,                      //its logical volume
                            "H2",                    //its name
                             Lead_LV5,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number


 Lead_LV5->SetVisAttributes(red);

 //place a convertor for neutron-proton conversion. 
  auto shield = new G4Box("shield", 70/2*cm, 70/2*cm, 3.0/2*mm);
  auto lShield = new G4LogicalVolume(shield, polyethylene, "Shield");    
  auto pShield = new G4PVPlacement(0,
                                               G4ThreeVector(0.*cm, 0.*cm, 112*cm),
                                               lShield,
                                               "Shield",
                                               lWorld,
                                               false,
                                               0,true);


//placing the aluminium coating.
  // The Coating
  auto coatingThickness = 100*nm; //100nm
  auto size = 70*cm;
  auto coatingSolid = new G4Box("coatingSolid", size/2, size/2, coatingThickness/2); // Adjust the size to cover the polyethylene shield completely
  auto lCoating = new G4LogicalVolume(coatingSolid, Aluminium, "AluminiumCoating");

  // Place aluminum coating around the polyethylene shield
  auto pCoating = new G4PVPlacement(0,
                                               G4ThreeVector(0.*cm, 0.*cm, 112.30001*cm),
                                               lCoating,
                                               "AluminiumCoating",
                                               lWorld,
                                               false,
                                               0,true);
  
  lCoating->SetVisAttributes(red);            

  // Absorber
  //

  auto 
  sAbsor = new G4Box("Absorber",                     //name
                     70/2*cm, 70/2*cm,fRadius/2);   //dimensions

  fLAbsor = new G4LogicalVolume(sAbsor,                  //shape
                             fMaterial,                 //material
                             fMaterial->GetName());     //name
                               
     new G4PVPlacement(0,                         //no rotation
                           G4ThreeVector(0*cm,0*cm,112.60002*cm),             //at (0,0,0)
                           fLAbsor,                     //logical volume
                           fMaterial->GetName(),        //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0,true);                          //copy number
  fScoringVolume = fLAbsor;


//place the polyprophlene foil here
  // The Foil.
  G4double foilThickness = 100*nm; //100nm
  G4double size_pp = 70*cm;
  G4Box* foilSolid = new G4Box("foilSolid", size_pp/2, size_pp/2, foilThickness/2); // Adjust the size to >
  G4LogicalVolume* PPCoating = new G4LogicalVolume(foilSolid, PP, "PPCoating");
  G4PVPlacement* foilCoating = new G4PVPlacement(0,
                                               G4ThreeVector(0.*cm, 0.*cm, 112.90003*cm),
                                               PPCoating,
                                               "PPCoating",
                                               lWorld,
                                               false,
                                               0,true);


  PPCoating->SetVisAttributes(yellow); 




  PrintParameters();
  
  //always return the root volume
  //
  return fPWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n The Absorber is " << G4BestUnit(fRadius,"Length")
         << " of " << fMaterial->GetName() 
         << "\n \n" << fMaterial << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    fMaterial = pttoMaterial;
    if(fLAbsor) { fLAbsor->SetMaterial(fMaterial); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetRadius(G4double value)
{
  fRadius = value;
  fWorldSize = 1.1*fRadius;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

