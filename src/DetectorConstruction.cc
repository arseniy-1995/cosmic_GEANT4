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
/// \file DetectorConstruction.cc
/// \brief Implementation of the Cosmic_sim::DetectorConstruction class

#include "DetectorConstruction.hh"
#include "PlasticSD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>
#include <TString.h>

namespace Cosmic_sim
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");

    // Liquid argon material
    G4double a;  // mass of a mole;
    G4double z;  // z=mean number of protons;
    G4double density;

    ////////////////////////////////////////////////////////////////////////////////
    // ELEMENTS :
    ////////////////////////////////////////////////////////////////////////////////

    G4Element *elH = new G4Element("Hydrogen", "H", 1, 1.01 * g / mole);
    G4Element *elC = new G4Element("Carbon", "C", 6, 12.01 * g / mole);
    G4Element *elN = new G4Element("Nitrogen", "N", 7, 14.01 * g / mole);
    G4Element *elO = new G4Element("Oxygen", "O", 8, 16.00 * g / mole);
    G4Element *elNa = new G4Element("Sodium", "Na", 11, 22.99*g/mole);
    G4Element *elCr = new G4Element("Chromium", "Cr", 24, 52.0 * g / mole);
    G4Element *elFe = new G4Element("Iron", "Fe", 26, 55.85 * g / mole);
    G4Element *elNi = new G4Element("Nickel", "Ni", 28, 58.69 * g / mole);
    G4Element *elCu = new G4Element("Copper", "Cu", 29, 63.55 * g / mole);
    G4Element *elZn = new G4Element("Zinc", "Zn", 30, 65.39 * g / mole);
    G4Element *elSi = new G4Element("Silicon", "Si", 14, 28.1 * g / mole);
    G4Element *elAl = new G4Element("Aluminium", "Al", 13, 26.9 * g / mole);
    G4Element* elW = new G4Element("Tungsten", "W",  z=74., a = 183.84*g/mole);
    G4Element* elPb = new G4Element("Lead","Pb", z=82., a= 207.20*g/mole);
    G4Element *elGe = new G4Element("German", "Ge", 32, 72.61*g/mole);
    G4Element *elI =  new G4Element("Iodine", "I", 53, 126.9*g/mole);
    G4Element *elCs = new G4Element("Cesium", "Cs", 55, 132.91*g/mole);
    G4Element *elBi = new G4Element("Bismuth", "Bi", 83, 208.98*g/mole);
    G4Element *elMg = new G4Element("Magnesium", "Mg", 12, 24.3*g/mole);


  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

  new G4Material("Vacuum", 1, a = 1.0 * g / mole, density = 1.0e-08 * g / cm3, kStateGas);


    G4Material* Air = new G4Material("Air", density = 1.29e-03 * g / cm3, 2);
    Air->AddElement(elN, .7);
    Air->AddElement(elO, .3);

    G4Material* Scintil = new G4Material("Scintillator", density = 1.032 * g / cm3, 2);

    //Scintil->AddElement(elH, 1.102);
    //Scintil->AddElement(elC, 1.0);
    Scintil->AddElement(elH, 10);
    Scintil->AddElement(elC, 9);


    G4Material* Scintil_BC_422 = new G4Material("Scintillator_BC_422", density = 1.032 * g / cm3, 2);

    //Scintil->AddElement(elH, 1.102);
    //Scintil->AddElement(elC, 1.0);
    Scintil_BC_422->AddElement(elH, 10);
    Scintil_BC_422->AddElement(elC, 9);

    //  G4Material* Air  = man->FindOrBuildMaterial("G4_AIR");

    G4Material* CONCRETE = nistManager->FindOrBuildMaterial("G4_CONCRETE"); //бетон, поиск в стандартной таблице


  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  fNofLayers = 10;
  G4double absoThickness = 10.*mm;
  G4double gapThickness =  5.*mm;
  G4double calorSizeXY  = 10.*cm;

  auto layerThickness = absoThickness + gapThickness;
  auto calorThickness = fNofLayers * layerThickness;
 // auto worldSizeXY = 1.2 * calorSizeXY;
 // auto worldSizeZ  = 1.2 * calorThickness;

    auto worldSizeX = 7.0 * m;
    auto worldSizeY = 8.0 * m;
    auto worldSizeZ  = 7.0 * m;

  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto absorberMaterial = G4Material::GetMaterial("G4_Pb");
  auto gapMaterial = G4Material::GetMaterial("liquidArgon");

  auto plasticMaterial = G4Material::GetMaterial("Scintillator");
  auto betonMaterial = G4Material::GetMaterial("G4_CONCRETE");


  if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial || !plasticMaterial || !betonMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }

  //
  // World
  //
  auto worldS
    = new G4Box("World",           // its name
                 worldSizeX/2, worldSizeY/2, worldSizeZ/2); // its size

  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name

  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps


   // Plastics

    G4double plasticFatNsys2SizeX = 1000. *mm, plasticFatNsys2SizeY = 120.* mm, plasticFatNsys2SizeZ = 120.* mm;
    G4double plasticFatNsys2PositionX = 0, plasticFatNsys2PositionY = 1.* m, plasticFatNsys2PositionZ = 0;

    G4double plasticFatNsys1SizeX = 1006. *mm, plasticFatNsys1SizeY = 200.* mm, plasticFatNsys1SizeZ = 200.* mm;
    G4double plasticFatNsys1PositionX = 0, plasticFatNsys1PositionY = -1.* m, plasticFatNsys1PositionZ = 0;

    G4double plasticThinNsys2SizeX[3] = {0,60.* cm,60.* cm}, plasticThinNsys2SizeY[3] = {0,1.* cm,1.* cm}, plasticThinNsys2SizeZ[3] = {0,80.* cm,80.* cm};
    G4double plasticThinNsys2PositionX[3] = {0,0,0}, plasticThinNsys2PositionY[3] = {0,0.5*m,0.5*m +0.5 *cm + plasticThinNsys2SizeY[1]}, plasticThinNsys2PositionZ[3] = {0,0,0};

    G4double plasticThinNsys1SizeX[3] = {0,30.* cm,30.* cm}, plasticThinNsys1SizeY[3] = {0,1.* cm,1.* cm}, plasticThinNsys1SizeZ[3] = {0,50.* cm,50.* cm};
    G4double plasticThinNsys1PositionX[3] = {0,15.0*cm +0.5* cm,-15.0*cm -0.5* cm}, plasticThinNsys1PositionY[3] = {0,-0.5*m,-0.5*m}, plasticThinNsys1PositionZ[3] = {0,0,0};

    G4double betonSizeX = 4. *m, betonSizeY = 1.* m, betonSizeZ = 4.* m;
    G4double betonPositionX = 0, betonPositionY = 3.5* m, betonPositionZ = 0;

    // ТОЛСТЫЕ ПЛАСТИКИ
    auto plastic_fat_nsys2S
            = new G4Box("plastic_fat_nsys2",     // its name
                        plasticFatNsys2SizeX/2, plasticFatNsys2SizeY/2, plasticFatNsys2SizeZ/2); // its size

    auto plastic_fat_nsys2LV
            = new G4LogicalVolume(
                    plastic_fat_nsys2S,     // its solid
                    plasticMaterial,  // its material
                    "plastic_fat_nsys2LV");   // its name

    new G4PVPlacement(
            0,                // no rotation
            G4ThreeVector(plasticFatNsys2PositionX,plasticFatNsys2PositionY, plasticFatNsys2PositionZ),  // at (0,0,0)
            plastic_fat_nsys2LV,          // its logical volume
            "plastic_fat_nsys2",    // its name
            worldLV,          // its mother  volume
            false,            // no boolean operation
            0,                // copy number
            fCheckOverlaps);  // checking overlaps



    auto plastic_fat_nsys1S
            = new G4Box("plastic_fat_nsys1",     // its name
                        plasticFatNsys1SizeX/2, plasticFatNsys1SizeY/2, plasticFatNsys1SizeZ/2); // its size

    auto plastic_fat_nsys1LV
            = new G4LogicalVolume(
                    plastic_fat_nsys1S,     // its solid
                    plasticMaterial,  // its material
                    "plastic_fat_nsys1LV");   // its name

    new G4PVPlacement(
            0,                // no rotation
            G4ThreeVector(plasticFatNsys1PositionX,plasticFatNsys1PositionY, plasticFatNsys1PositionZ),  // at (0,0,0)
            plastic_fat_nsys1LV,          // its logical volume
            "plastic_fat_nsys1",    // its name
            worldLV,          // its mother  volume
            false,            // no boolean operation
            0,                // copy number
            fCheckOverlaps);  // checking overlaps

// ТОНКИЕ ПЛАСТИКИ

    G4Box  *plastic_thin_nsys2S[3], *plastic_thin_nsys1S[3];
    G4LogicalVolume *plastic_thin_nsys2LV[3], *plastic_thin_nsys1LV[3];

for (int i=1;i<=2;i++) {
    plastic_thin_nsys2S[i]
            = new G4Box(Form("plastic_thin%i_nsys2", i),     // its name
                        plasticThinNsys2SizeX[i] / 2, plasticThinNsys2SizeY[i] / 2,
                        plasticThinNsys2SizeZ[i] / 2); // its size

    plastic_thin_nsys2LV[i]
            = new G4LogicalVolume(
                    plastic_thin_nsys2S[i],     // its solid
                    plasticMaterial,  // its material
                    Form("plastic_thin%i_nsys2LV",i));   // its name

    new G4PVPlacement(
            0,                // no rotation
            G4ThreeVector(plasticThinNsys2PositionX[i], plasticThinNsys2PositionY[i],
                          plasticThinNsys2PositionZ[i]),  // at (0,0,0)
            plastic_thin_nsys2LV[i],          // its logical volume
            Form("plastic_thin%i_nsys2",i),    // its name
            worldLV,          // its mother  volume
            false,            // no boolean operation
            0,                // copy number
            fCheckOverlaps);  // checking overlaps

  plastic_thin_nsys1S[i]
            = new G4Box(Form("plastic_thin%i_nsys1",i),     // its name
                        plasticThinNsys1SizeX[i] / 2, plasticThinNsys1SizeY[i] / 2,
                        plasticThinNsys1SizeZ[i] / 2); // its size

   plastic_thin_nsys1LV[i]
            = new G4LogicalVolume(
                    plastic_thin_nsys1S[i],     // its solid
                    plasticMaterial,  // its material
                    Form("plastic_thin%i_nsys1LV",i));   // its name

    new G4PVPlacement(
            0,                // no rotation
            G4ThreeVector(plasticThinNsys1PositionX[i],plasticThinNsys1PositionY[i], plasticThinNsys1PositionZ[i]),  // at (0,0,0)
            plastic_thin_nsys1LV[i],          // its logical volume
            Form("plastic_thin%i_nsys1",i),    // its name
            worldLV,          // its mother  volume
            false,            // no boolean operation
            0,                // copy number
            fCheckOverlaps);  // checking overlaps


}

// БЕТОН


    auto betonS
            = new G4Box("beton",     // its name
                        betonSizeX/2, betonSizeY/2, betonSizeZ/2); // its size

    auto betonLV
            = new G4LogicalVolume(
                    betonS,     // its solid
                    betonMaterial,  // its material
                    "betonLV");   // its name

    new G4PVPlacement(
            0,                // no rotation
            G4ThreeVector(betonPositionX,betonPositionY, betonPositionZ),  // at (0,0,0)
            betonLV,          // its logical volume
            "beton",    // its name
            worldLV,          // its mother  volume
            false,            // no boolean operation
            0,                // copy number
            fCheckOverlaps);  // checking overlaps


  /*

  //
  // Calorimeter
  //
  auto calorimeterS
    = new G4Box("Calorimeter",     // its name
                 calorSizeXY/2, calorSizeXY/2, calorThickness/2); // its size

  auto calorLV
    = new G4LogicalVolume(
                 calorimeterS,     // its solid
                 defaultMaterial,  // its material
                 "Calorimeter");   // its name

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 calorLV,          // its logical volume
                 "Calorimeter",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  //
  // Layer
  //
  auto layerS
    = new G4Box("Layer",           // its name
                 calorSizeXY/2, calorSizeXY/2, layerThickness/2); //its size

  auto layerLV
    = new G4LogicalVolume(
                 layerS,           // its solid
                 defaultMaterial,  // its material
                 "Layer");         // its name

  new G4PVReplica(
                 "Layer",          // its name
                 layerLV,          // its logical volume
                 calorLV,          // its mother
                 kZAxis,           // axis of replication
                 fNofLayers,        // number of replica
                 layerThickness);  // witdth of replica

  //
  // Absorber
  //
  auto absorberS
    = new G4Box("Abso",            // its name
                 calorSizeXY/2, calorSizeXY/2, absoThickness/2); // its size

  auto absorberLV
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "AbsoLV");        // its name

   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -gapThickness/2), // its position
                 absorberLV,       // its logical volume
                 "Abso",           // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  //
  // Gap
  //
  auto gapS
    = new G4Box("Gap",             // its name
                 calorSizeXY/2, calorSizeXY/2, gapThickness/2); // its size

  auto gapLV
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "GapLV");         // its name

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., absoThickness/2), // its position
                 gapLV,            // its logical volume
                 "Gap",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  //
  // print parameters
  //
  G4cout
    << G4endl
    << "------------------------------------------------------------" << G4endl
    << "---> The calorimeter is " << fNofLayers << " layers of: [ "
    << absoThickness/mm << "mm of " << absorberMaterial->GetName()
    << " + "
    << gapThickness/mm << "mm of " << gapMaterial->GetName() << " ] " << G4endl
    << "------------------------------------------------------------" << G4endl;

   */

  //
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
 // calorLV->SetVisAttributes(simpleBoxVisAtt);

    plastic_fat_nsys1LV->SetVisAttributes(simpleBoxVisAtt);
    plastic_fat_nsys2LV->SetVisAttributes(simpleBoxVisAtt);

    plastic_thin_nsys1LV[1]->SetVisAttributes(simpleBoxVisAtt);
    plastic_thin_nsys2LV[1]->SetVisAttributes(simpleBoxVisAtt);

    plastic_thin_nsys1LV[2]->SetVisAttributes(simpleBoxVisAtt);
    plastic_thin_nsys2LV[2]->SetVisAttributes(simpleBoxVisAtt);

    betonLV->SetVisAttributes(simpleBoxVisAtt);

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  //
  // Sensitive detectors
  //

    auto plastic_fat_nsys1SD
            = new PlasticSD("plastic_fat_nsys1SD", "plastic_fat_nsys1HitsCollection", fNofLayers);
    G4SDManager::GetSDMpointer()->AddNewDetector(plastic_fat_nsys1SD);
    SetSensitiveDetector("plastic_fat_nsys1LV",plastic_fat_nsys1SD);

/*
  auto absoSD
    = new CalorimeterSD("AbsorberSD", "AbsorberHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(absoSD);
  SetSensitiveDetector("AbsoLV",absoSD);

  auto gapSD
    = new CalorimeterSD("GapSD", "GapHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(gapSD);
  SetSensitiveDetector("GapLV",gapSD);

 */

  //
  // Magnetic field
  //
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
