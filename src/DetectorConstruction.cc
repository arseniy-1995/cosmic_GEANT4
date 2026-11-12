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
/// \brief Implementation of the Cosmic::DetectorConstruction class

#include "DetectorConstruction.hh"
#include "PlasticSD.hh"
#include "HadronCalorimeterSD.hh"
#include "ChamberSD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include <G4EllipticalTube.hh>
#include <G4Trap.hh>
#include <G4Trd.hh>
#include <TString.h>
#include "G4AutoDelete.hh"
#include "G4Box.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"

#include "G4MultiUnion.hh"
#include "G4Transform3D.hh"


#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>
//#include <TString.h>

namespace Cosmic {

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    G4ThreadLocal
    G4GlobalMagFieldMessenger
            *
            DetectorConstruction::fMagFieldMessenger = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    DetectorConstruction::DetectorConstruction() {
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    DetectorConstruction::~DetectorConstruction() {
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    G4VPhysicalVolume *DetectorConstruction::Construct() {
        // Define materials
        DefineMaterials();

        // Define rotation
        DefineRotationMatrices();

        // Define vis attributes
        DefineVisAttributes();


        // Define volumes
        return DefineVolumes();
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void DetectorConstruction::DefineMaterials() {
        // Lead material defined using NIST Manager
        auto nistManager = G4NistManager::Instance();
        nistManager->FindOrBuildMaterial("G4_Pb");
        nistManager->FindOrBuildMaterial("G4_W");

        // Liquid argon material
        G4double a;  // mass of a mole;
        G4double z;  // z=mean number of protons;
        G4double density;

        G4int ncomponents, natoms;
        G4double fractionmass;
        G4double temperature, pressure;

        ////////////////////////////////////////////////////////////////////////////////
        // ELEMENTS :
        ////////////////////////////////////////////////////////////////////////////////

        G4Element *elH = new G4Element("Hydrogen", "H", 1, 1.01 * g / mole);
        G4Element *elD = new G4Element("Deuterium", "D", 1);
        G4Element *elC = new G4Element("Carbon", "C", 6, 12.01 * g / mole);
        G4Element *elBe = new G4Element("Beryllium", "Be", 4, 9.1218 * g / mole);
        G4Element *elN = new G4Element("Nitrogen", "N", 7, 14.01 * g / mole);
        G4Element *elO = new G4Element("Oxygen", "O", 8, 16.00 * g / mole);
        G4Element *elNa = new G4Element("Sodium", "Na", 11, 22.99 * g / mole);
        G4Element *elCr = new G4Element("Chromium", "Cr", 24, 52.0 * g / mole);
        G4Element *elFe = new G4Element("Iron", "Fe", 26, 55.85 * g / mole);
        G4Element *elNi = new G4Element("Nickel", "Ni", 28, 58.69 * g / mole);
        G4Element *elCu = new G4Element("Copper", "Cu", 29, 63.55 * g / mole);
        G4Element *elZn = new G4Element("Zinc", "Zn", 30, 65.39 * g / mole);
        G4Element *elSi = new G4Element("Silicon", "Si", 14, 28.1 * g / mole);
        G4Element *elAl = new G4Element("Aluminium", "Al", 13, 26.9 * g / mole);
        G4Element *elW = new G4Element("Tungsten", "W", z = 74., a = 183.84 * g / mole);
        G4Element *elPb = new G4Element("Lead", "Pb", z = 82., a = 207.20 * g / mole);
        G4Element *elGe = new G4Element("German", "Ge", 32, 72.61 * g / mole);
        G4Element *elI = new G4Element("Iodine", "I", 53, 126.9 * g / mole);
        G4Element *elCs = new G4Element("Cesium", "Cs", 55, 132.91 * g / mole);
        G4Element *elBi = new G4Element("Bismuth", "Bi", 83, 208.98 * g / mole);
        G4Element *elMg = new G4Element("Magnesium", "Mg", 12, 24.3 * g / mole);


        new G4Material("liquidArgon", z = 18., a = 39.95 * g / mole, density = 1.390 * g / cm3);
        // The argon by NIST Manager is a gas with a different density

        // Vacuum
        Galactic = new G4Material("Galactic", z = 1., a = 1.01 * g / mole, density = universe_mean_density, kStateGas, 2.73 * kelvin, 3.e-18 * pascal);

        Vacuum = new G4Material("Vacuum", 1, a = 1.0 * g / mole, density = 1.0e-08 * g / cm3, kStateGas);
        // Vacuum
        auto vacuum_mt = new G4MaterialPropertiesTable();
        vacuum_mt->AddProperty("RINDEX", "Air");
        Vacuum->SetMaterialPropertiesTable(vacuum_mt);
        new G4Material("Copper", 29, 63.55 * g / mole, 8.96 * g / cm3);

        G4Material *Air = new G4Material("Air", density = 1.29e-03 * g / cm3, 2);
        Air->AddElement(elN, .7);
        Air->AddElement(elO, .3);

        G4Material *Scintil = new G4Material("Scintillator", density = 1.032 * g / cm3, 2);

        //Scintil->AddElement(elH, 1.102);
        //Scintil->AddElement(elC, 1.0);
        Scintil->AddElement(elH, 10);
        Scintil->AddElement(elC, 9);


        G4Material *Scintil_BC_422 = new G4Material("Scintillator_BC_422", density = 1.032 * g / cm3, 2);

        //Scintil->AddElement(elH, 1.102);
        //Scintil->AddElement(elC, 1.0);
        Scintil_BC_422->AddElement(elH, 10);
        Scintil_BC_422->AddElement(elC, 9);

        //  G4Material* Air  = man->FindOrBuildMaterial("G4_AIR");

        G4Material *Concrete = nistManager->FindOrBuildMaterial("G4_CONCRETE"); //бетон, поиск в стандартной таблице


        G4Material *Wood = nistManager->FindOrBuildMaterial("G4_CELLULOSE_CELLOPHANE");

        Alumin = new G4Material("Aluminum", 13, 26.98 * g / mole, 2.7 * g / cm3);
        Titan = new G4Material("Titanium", 22, 47.87 * g / mole, 4.54 * g / cm3);

        // G4Material *Beryllium = new G4Material("Beryllium", 4, 9.1218 * g / mole, 4.54 * g / cm3);
        Beryllium = new G4Material("Beryllium", 1.85 * g / cm3, 1);
        Beryllium->AddElement(elBe, 1.0);

        // отражающая обертка
        Mylar = new G4Material("Mylar", density = 1.39 * g / cm3, 3);
        Mylar->AddElement(elC, 5);
        Mylar->AddElement(elO, 2);
        Mylar->AddElement(elH, 4);


        Steel = new G4Material("Steel", density = 8.89 * g / cm3, 3);
        Steel->AddElement(elCr, .197);
        Steel->AddElement(elFe, .704);
        Steel->AddElement(elNi, .099);

        Iron = new G4Material("Iron", density = 8.89 * g / cm3, 1);
        Iron->AddElement(elFe, 1.0);

        // для камер

        SiO2 = new G4Material("SiO2", 2.2 * g / cm3, 2);
        SiO2->AddElement(elSi, 1);
        SiO2->AddElement(elO, 2);

        GasAr = new G4Material("ArgonGas", 18, 39.95 * g / mole, 1.7839e-03 * g / cm3);

        GasCO2 = new G4Material("CarbonicGas", 1.977e-03 * g / cm3, 2);
        GasCO2->AddElement(elC, 1);
        GasCO2->AddElement(elO, 2);

        G4Material *GasIsobutane = new G4Material("IsobutaneGas", 2.486e-03 * g / cm3, 2);
        GasIsobutane->AddElement(elC, 4);
        GasIsobutane->AddElement(elH, 10);

        GasArCO2 = new G4Material("ArgonCarbonicGas", 1.98e-03 * g / cm3, 2);
        // GasArCO2->AddMaterial(GasAr, .2); // .8 // .82
      //  GasArCO2->AddMaterial(GasCO2, .8); // .2 // .18
        GasArCO2->AddMaterial(GasAr, .82);
        GasArCO2->AddMaterial(GasCO2, .18);


        G4Material *GasArIsobutane = new G4Material("ArgonIsobutaneGas", 1.910278e-03 * g / cm3, 2);
        GasArIsobutane->AddMaterial(GasAr, .82);
        GasArIsobutane->AddMaterial(GasIsobutane, .18);

        G4Material *Convertor_LQ = nistManager->FindOrBuildMaterial("G4_Pb");


        // газ дейтерий для мишени

        G4Isotope *isoD = new G4Isotope("Deuteron", 1, 2, 2.01410178 * g / mole);
        elD->AddIsotope(isoD, 100. * perCent);
        static const G4double torr = atmosphere / 760.0;
        pressure = 1e-6 * torr;
        G4Material *deuteriumGas = new G4Material("DeuteriumGas", (pressure / (1.0 * atmosphere)) * 0.167 * kg / m3, 1,
                                                  kStateGas, 77.0 * kelvin, pressure);
        // deuteriumGas->AddElement(elD, 2); // В молекуле D2 два атома дейтерия
        deuteriumGas->AddElement(elD, 1); // атомарный дейтерий

        // газ водород для мишени


        pressure = 1e-6 * torr;
        G4Material *customH2 = new G4Material("HydrogenGas", (pressure / (1.0 * atmosphere)) * 0.083748 * mg / cm3, 1,
                                              kStateGas, 77.0 * kelvin, pressure);
        // customH2->AddElement(elH, 2); // Молекула H2 состоит из 2 атомов
        customH2->AddElement(elH, 1);


        // материал каптон для GEM

        G4Material *Kapton = new G4Material("Kapton", 1.42 * g / cm3, 4); // 1.032 * g / cm3 1.42 * g/cm3;
        //  Kapton->AddElement(elH, 0.0273);
        //  Kapton->AddElement(elC, 0.7213);
        //  Kapton->AddElement(elN,  0.0765);
        //  Kapton->AddElement(elO, 0.1749);

        Kapton->AddElement(elH, 0.026362);
        Kapton->AddElement(elC, 0.691133);
        Kapton->AddElement(elN, 0.073270);
        Kapton->AddElement(elO, 0.209235);

        // Kapton->AddElement(elC, 22);
        // Kapton->AddElement(elH, 10);
        // Kapton->AddElement(elO, 5);
        // Kapton->AddElement(elN, 2);

        // Print materials
        G4cout << *(G4Material::GetMaterialTable()) << G4endl;
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


    void DetectorConstruction::DefineRotationMatrices() {
        ////////////////////////////////////////////////////////////////////////////////
        // ROTATION MATRICES :
        ////////////////////////////////////////////////////////////////////////////////

        //   Rotate30X.rotateX(-30*deg);
        Rotate180X.rotateX(180 * deg);
        Rotate180Y.rotateY(180 * deg);
        Rotate180Z.rotateZ(180 * deg);
        Rotate270X.rotateX(-90 * deg);
        Rotate90X.rotateX(90 * deg);
        Rotate9X.rotateX(9 * deg);
        Rotate99X.rotateX(99 * deg);
        Rotate45X.rotateX(45 * deg);
        Rotate45X180Z.rotateX(45 * deg);
        Rotate45X180Z.rotateZ(180 * deg);
        Rotate10X.rotateX(10 * deg);
        Rotate25X.rotateX(25 * deg);
        Rotate65X.rotateX(65 * deg);
        Rotate90Y.rotateY(90 * deg);
        Rotate270Y.rotateY(-90 * deg);
        Rotate90X180Z.rotateX(90 * deg);
        Rotate90X180Z.rotateZ(180 * deg);
        Rotate90Y180Z.rotateY(90 * deg);
        Rotate90Y180Z.rotateZ(180 * deg);
        Rotate90Z.rotateZ(90 * deg);
        Rotate270Z.rotateZ(-90 * deg);
        Rotate270Y180X.rotateY(-90 * deg);
        Rotate270Y180X.rotateX(180 * deg);
        Rotate90Y180X.rotateY(90 * deg);
        Rotate90Y180X.rotateX(180 * deg);
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void DetectorConstruction::DefineVisAttributes() {


        G4Colour steel_col(0., 0.0, 1.0);
        G4Colour iron_col(0.5, 0.5, 0.8);
        G4Colour alum_col(0.0, 0.8, 0.8);
        G4Colour titanfoil_col(0.8, 0.8, 0.8);
        G4Colour berylliumfoil_col(0.8, 0.8, 0.8);
        G4Colour copper_col(1.0, 0.0, 0.0);
        G4Colour convertor_col(0.0, 0.5, 0.5);
        G4Colour brass_col(0.75, 0.75, 0.1);
        G4Colour stef_col(0.742, 0.699, 0.121);
        G4Colour plastic_col(1.0, 0.0, 1.0);
        G4Colour gas_col(0.93, 0.77, 0.57);
        G4Colour mylar_col(0.5, 0.2, 0.2);
        G4Colour alfoil_col(0.9, 0.9, 0.9);
        G4Colour blackpaper_col(0.2, 0.2, 0.2);
        G4Colour concrete_col(0., 0.1, 0.8);
        G4Colour WOOD_col(0., 0.9, 0.1);
        G4Colour CONVERTOR_LQ_col(0.9, 0.1, 0.1);

        Gas_VisAtt = new G4VisAttributes(gas_col);
        Foil_VisAtt = new G4VisAttributes(alfoil_col);
        Mylar_VisAtt = new G4VisAttributes(mylar_col);
        Stef_VisAtt = new G4VisAttributes(stef_col);
        Steel_VisAtt = new G4VisAttributes(steel_col);
        Mag_VisAtt = new G4VisAttributes(copper_col);
        Shield_VisAtt = new G4VisAttributes(copper_col);
        Iron_VisAtt = new G4VisAttributes(iron_col);
        Alum_VisAtt = new G4VisAttributes(alum_col);
        Plastic_VisAtt = new G4VisAttributes(plastic_col);
        Convertor_VisAtt = new G4VisAttributes(convertor_col);
        Concrete_VisAtt = new G4VisAttributes(concrete_col);
        WOOD_VisAtt = new G4VisAttributes(WOOD_col);
        Convertor_LQ_VisAtt = new G4VisAttributes(CONVERTOR_LQ_col);
        ProCover_VisAtt = new G4VisAttributes(blackpaper_col); // обертка майлар
        TitanFoil_VisAtt = new G4VisAttributes(titanfoil_col);
        BerylliumFoil_VisAtt = new G4VisAttributes(berylliumfoil_col);


    }


    G4VPhysicalVolume *DetectorConstruction::DefineVolumes() {
        // Geometry parameters
        // fNofLayers =10;
        // G4double absoThickness = 10.*mm;
        // G4double gapThickness =  5.*mm;
        // G4double calorSizeXY  = 10.*cm;

//  auto layerThickness = absoThickness + gapThickness;
//  auto calorThickness = fNofLayers * layerThickness;
        // auto worldSizeXY = 1.2 * calorSizeXY;
        // auto worldSizeZ  = 1.2 * calorThickness;


        auto worldSizeX = 10.0 * m;
        auto worldSizeY = 35.0 * m;
        auto worldSizeZ = 10.0 * m;

        // Get materials
        auto defaultMaterial = G4Material::GetMaterial("Galactic");
        auto absorberMaterial = G4Material::GetMaterial("G4_Pb");
        auto gapMaterial = G4Material::GetMaterial("liquidArgon");

        auto plasticMaterial = G4Material::GetMaterial("Scintillator");
        auto betonMaterial = G4Material::GetMaterial("G4_CONCRETE");
        auto AirMaterial = G4Material::GetMaterial("Air");
        auto MylarMaterial = G4Material::GetMaterial("Mylar");

        if (!defaultMaterial || !absorberMaterial || !gapMaterial || !plasticMaterial || !betonMaterial ||
            !AirMaterial || !MylarMaterial) {
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
                            worldSizeX / 2, worldSizeY / 2, worldSizeZ / 2); // its size

        fWorldLV = new G4LogicalVolume(
                worldS,           // its solid
                defaultMaterial,  // its material
                "World");         // its name

        worldPV
                = new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(),  // at (0,0,0)
                                    fWorldLV, // its logical volume
                                    "World",          // its name
                0,                // its mother  volume
                false,            // no boolean operation
                0,                // copy number
                fCheckOverlaps);  // checking overlaps


        // поворот мишени

        auto rot_for_world_target = new G4RotationMatrix();
        rot_for_world_target->rotateY(alpha0_GEM);

        G4ThreeVector size_world_rot(0.15 * m, 0.1 * m, 0.3 * m);

        G4ThreeVector world_rot_pos = G4ThreeVector(0., 0., 0.);


#ifndef isGenElectronGEM
        size_world_rot = G4ThreeVector(0.01 * cm, 0.01 * cm, 0.01 * cm);
#endif

#ifndef isGenElectronGEM
        world_rot_pos = G4ThreeVector(-0.2 * m, -0.2 * m, 0.0 * m);
#endif

        G4Box *world_box_rot = new G4Box("WorldRot", size_world_rot.x(), size_world_rot.y(), size_world_rot.z());

        auto VacuumMaterial = G4Material::GetMaterial("Vacuum");
        // auto fWorldRotMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
        // auto fWorldRotMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
        auto fWorldRotMaterial = VacuumMaterial;
        fWorldLV_rot = new G4LogicalVolume(world_box_rot, fWorldRotMaterial, "World_rot_log", 0, 0, 0);
        fWorldLV_rot->SetVisAttributes(Gas_VisAtt);

        auto worldPV_rot = new G4PVPlacement(rot_for_world_target, world_rot_pos, fWorldLV_rot, "World_rot", fWorldLV,
                                             false, fCheckOverlaps);


        // Plastics

        // ТОЛСТЫЕ ПЛАСТИКИ

        // Neutron/Proton E/dE counters #2   (above the beam)
#ifdef PF2_FAT
        ConstructPlasticFat2();
#endif //PF2_FAT


        ///////////////////////////////////
// Neutron/Proton E/dE counters #1 (below the beam)

#ifdef PF1_FAT
        ConstructPlasticFat1();
#endif //PF1_FAT


        // ТОНКИЕ ПЛАСТИКИ
        // ВЕРХНИЕ ПЛЕЧО
#ifdef PF2_THIN
        ConstructPlasticThin2();
#endif //PF2_THIN

        // нижнее ПЛЕЧО
#ifdef PF1_THIN
        ConstructPlasticThin1();
#endif //PF1_THIN

        /// КОЛОРИМЕТР ВЫНЕСЕН В ОТДЕЛЬНУЮ ФУНКЦИЮ

#if defined(HADCAL1) || defined(HADCAL1)
        G4LogicalVolume *HadronCalorimeter_nsys1LV = ConstructHadronCalorimeter(1);
        G4LogicalVolume *HadronCalorimeter_nsys2LV = ConstructHadronCalorimeter(2);

        // HadronCalorimeter DATA

        // in this version sandwich is square with all strips (X and Z) are the same
        // and number of X bars equals to number of Z bars

        G4int NbOfXBars;    // number of bars for phi
        G4int NbOfZBars;    // number of bars for theta
        G4double ScintSizeX;// 10.*cm;
        G4double ScintSizeZ;// 10.*cm;

        G4double GapX = 0.5 * mm;
        G4double GapY = 0.25 * mm;

        G4double ScintThickness = 7. * mm;
        ScintSizeX = 79.5 * mm;
        ScintSizeZ = NX_BARS * (ScintSizeX + GapX);

        G4double IronThickness = 16. * mm;
        G4double IronSizeX = 16. * cm;
        G4double IronSizeZ = ScintSizeZ + 2. * 10. * cm;

        G4int NbOfLayers = N_LAYERS;// Hadron Calo layers */
        NbOfXBars = NX_BARS;  // number of bars for phi
        NbOfZBars = NZ_BARS;  // number of bars for theta

        G4double LayerStep = 35. * mm;

        // Compute sizes
        G4double StripSizeX = IronSizeX;//(ScintSizeX+GapX)*2;
        G4double StripSizeZ = IronSizeZ + 2 * GapX;
        G4double StripSizeY = (ScintThickness + GapY) * 2 + IronThickness + 4 * GapY;

        G4double LayerSizeX = StripSizeX * NbOfXBars / 2;
        G4double LayerSizeY = StripSizeY;
        G4double LayerSizeZ = StripSizeZ;

        G4double SandSizeX = IronSizeZ;
        G4double SandSizeY = LayerStep * (NbOfLayers + 0.5);// no ACC layers
        G4double SandSizeZ = IronSizeZ;

        G4double VertPos = 150.0 * cm;// 150.*cm; !! now to fron face of SANDW
        G4double HorPos = 73.8 * cm;

        G4double y_pos = VertPos + SandSizeY / 2.;
        //  new G4PVPlacement(G4Transform3D(RotateNull, G4ThreeVector(0.0 * cm, y_pos, HorPos)),HadronCalorimeterLV,
        //  "HadronCalorimeterNsys2PV", fWorldLV, false, ARM2_IND); // RIA new G4PVPlacement(G4Transform3D(Rotate180Z,
        //  G4ThreeVector(0.0 * cm, -y_pos, HorPos)),HadronCalorimeterLV, "HadronCalorimeterNsys1PV", fWorldLV, false,
        //  ARM1_IND); // RIA

#endif

#ifdef HADCAL1

        G4int pCopyNo_HADCAL_nsys1 = ARM1_IND;
        // G4int pCopyNo_HADCAL_nsys1 = 0;

        G4ThreeVector pos_HAD;
#ifdef RUN21
            pos_HAD=  G4ThreeVector(0.0 * cm, -171. * cm, 84. *cm); // y-position and z-position modyfied by Gauzshtein
#endif
#ifdef RUN23
        pos_HAD = G4ThreeVector(0.0 * cm, -171. * cm, 78.4 * cm);
#endif

        new G4PVPlacement(G4Transform3D(Rotate180Z, pos_HAD), HadronCalorimeter_nsys1LV, "Sand_phys_nsys1", fWorldLV,
                          false, pCopyNo_HADCAL_nsys1, fCheckOverlaps);
#endif
#ifdef HADCAL2
        G4int pCopyNo_HADCAL_nsys2 = ARM2_IND;
        //G4int pCopyNo_HADCAL_nsys2 = 0;
#ifdef RUN21
            pos_HAD =  G4ThreeVector(0.0 * cm, 171. * cm, 78. * cm);
#endif
#ifdef RUN23
        pos_HAD = G4ThreeVector(0.0 * cm, 171. * cm, 78.4 * cm);
#endif

        new G4PVPlacement(G4Transform3D(RotateNull, pos_HAD), HadronCalorimeter_nsys2LV, "Sand_phys_nsys2", fWorldLV,
                          false, pCopyNo_HADCAL_nsys2, fCheckOverlaps);
#endif


        // Vertex Chambers
#if defined(VCARM1) && defined(RUN21)
        auto VCBox_log_nsys1 = ConstructVC(VCGas_log_nsys1LV);
        //G4int pCopyNo_VC_nsys1 = ARM1_IND;
        G4int pCopyNo_VC_nsys1 = 0;
        new G4PVPlacement(G4Transform3D(Rotate180Z, G4ThreeVector(0.0 * cm, -8.2 * cm, 0. * cm)), VCBox_log_nsys1,
                          "VCBox_phys_nsys1", fWorldLV, false, pCopyNo_VC_nsys1, fCheckOverlaps);
#endif

#if defined(VCARM2) && defined(RUN21)
        auto VCBox_log_nsys2 = ConstructVC(VCGas_log_nsys2LV);
        //G4int pCopyNo_VC_nsys2 = ARM2_IND;
        G4int pCopyNo_VC_nsys2 = 0;
        new G4PVPlacement(G4Transform3D(RotateNull, G4ThreeVector(0.0 * cm, 8.2 * cm, 0. * cm)), VCBox_log_nsys2,
                          "VCBox_phys_nsys2", fWorldLV, false, pCopyNo_VC_nsys2, fCheckOverlaps);
#endif


        // КАМЕРЫ

#ifdef DCARM1 // Конструкторы
        WCTheta1_gas_nsys1LV = NULL;
        G4LogicalVolume *WCTheta1_nsys1LV = ConstructWC(NW1_WRS * 2. * cm, 20. * cm, WC1_IND, WCTheta1_gas_nsys1LV);

        WCPhi1_gas_nsys1LV = NULL;
        G4LogicalVolume *WCPhi1_nsys1LV = ConstructWC(NW2_WRS * 2. * cm, 45. * cm, WC2_IND, WCPhi1_gas_nsys1LV);

        WCTheta2_gas_nsys1LV = NULL;
        G4LogicalVolume *WCTheta2_nsys1LV = ConstructWC(NW3_WRS * 2. * cm, 43. * cm, WC3_IND, WCTheta2_gas_nsys1LV);

#endif // DCARM1

#ifdef DCARM2 // Конструкторы

        WCTheta1_gas_nsys2LV = NULL;
        G4LogicalVolume *WCTheta1_nsys2LV = ConstructWC(NW1_WRS * 2. * cm, 20. * cm, WC1_IND, WCTheta1_gas_nsys2LV);

        WCPhi1_gas_nsys2LV = NULL;
        G4LogicalVolume *WCPhi1_nsys2LV = ConstructWC(NW2_WRS * 2. * cm, 45. * cm, WC2_IND, WCPhi1_gas_nsys2LV);

        WCTheta2_gas_nsys2LV = NULL;
        G4LogicalVolume *WCTheta2_nsys2LV = ConstructWC(NW3_WRS * 2. * cm, 43. * cm, WC3_IND, WCTheta2_gas_nsys2LV);

#endif // DCARM1

//////////////////////////////////
// Place chambers !! /////////////
//////////////////////////////////

// ARM #1
#ifdef DCARM1

        G4int pCopyNo_WC_nsys1 = ARM1_IND;
        // G4int pCopyNo_WC_nsys1 = 0;

#ifdef RUN21
        new G4PVPlacement(G4Transform3D(Rotate180Z,
                                        G4ThreeVector(0.0 * cm, -14.8 * cm, 9.4 * cm)), WCTheta1_nsys1LV,
                          "WCTheta1a_nsys1", fWorldLV, false, pCopyNo_WC_nsys1, fCheckOverlaps);
        new G4PVPlacement(G4Transform3D(Rotate180Z,
                                        G4ThreeVector(0.0 * cm, -34.2 * cm, 19.4 * cm)), WCTheta2_nsys1LV,
                          "WCTheta2a_nsys1", fWorldLV, false, pCopyNo_WC_nsys1, fCheckOverlaps);
        new G4PVPlacement(G4Transform3D(Rotate90Y180Z,
                                        G4ThreeVector(0.0 * cm, -26.5 * cm, 15.2 * cm)), WCPhi1_nsys1LV,
                          "WCPhi1a_nsys1", fWorldLV, false, pCopyNo_WC_nsys1, fCheckOverlaps);

#endif

 #ifdef RUN23
        new G4PVPlacement(G4Transform3D(Rotate180Z,
                                        G4ThreeVector(0.0 * cm, -14.8 * cm - 2. * mm, 9.4 * cm)),
                          WCTheta1_nsys1LV, "WCTheta1a_nsys1", fWorldLV, false, pCopyNo_WC_nsys1, fCheckOverlaps);
        new G4PVPlacement(G4Transform3D(Rotate180Z,
                                        G4ThreeVector(0.0 * cm, -34.2 * cm - 2. * mm, 18.8 * cm)),
                          WCTheta2_nsys1LV, "WCTheta2a_nsys1", fWorldLV, false, pCopyNo_WC_nsys1, fCheckOverlaps);
        new G4PVPlacement(G4Transform3D(Rotate90Y180Z,
                                        G4ThreeVector(0.0 * cm, -26.5 * cm - 2. * mm, 15.2 * cm)),
                          WCPhi1_nsys1LV, "WCPhi1a_nsys1", fWorldLV, false, pCopyNo_WC_nsys1, fCheckOverlaps);

#endif

#endif
// ARM #2
#ifdef DCARM2

        G4int pCopyNo_WC_nsys2 = ARM2_IND;
        //G4int pCopyNo_WC_nsys2 = 0;

#ifdef RUN21
        new G4PVPlacement(G4Transform3D(RotateNull,
                                        G4ThreeVector(0.0 * cm, 14.8 * cm, 9.4 * cm)), WCTheta1_nsys2LV,
                          "WCTheta1b_nsys2", fWorldLV, false, pCopyNo_WC_nsys2, fCheckOverlaps);
        new G4PVPlacement(G4Transform3D(RotateNull,
                                        G4ThreeVector(0.0 * cm, 34.2 * cm , 19.4 * cm)), WCTheta2_nsys2LV,
                          "WCTheta2b_nsys2", fWorldLV, false, pCopyNo_WC_nsys2, fCheckOverlaps);
        new G4PVPlacement(G4Transform3D(Rotate90Y,
                                        G4ThreeVector(0.0 * cm, 26.5 * cm, 15.2 * cm)), WCPhi1_nsys2LV,
                          "WCPhi1b_nsys2", fWorldLV, false, pCopyNo_WC_nsys2, fCheckOverlaps);

#endif


#ifdef RUN23
        new G4PVPlacement(G4Transform3D(RotateNull,
                                        G4ThreeVector(0.0 * cm, 14.8 * cm + 7. * mm, 9.4 * cm)),
                          WCTheta1_nsys2LV, "WCTheta1b_nsys2", fWorldLV, false, pCopyNo_WC_nsys2, fCheckOverlaps);
        new G4PVPlacement(G4Transform3D(RotateNull,
                                        G4ThreeVector(0.0 * cm, 34.2 * cm + 7. * mm, 18.8 * cm)),
                          WCTheta2_nsys2LV, "WCTheta2b_nsys2", fWorldLV, false, pCopyNo_WC_nsys2, fCheckOverlaps);
        new G4PVPlacement(G4Transform3D(Rotate90Y,
                                        G4ThreeVector(0.0 * cm, 26.5 * cm  + 7. * mm, 15.2 * cm)),
                          WCPhi1_nsys2LV, "WCPhi1b_nsys2", fWorldLV, false, pCopyNo_WC_nsys2, fCheckOverlaps);

#endif

#endif


        // Электронные счетчики LO-поляриметра
#if defined(LOWQ1) || defined(LOWQ2)

        auto LQBox_log_nsys1 = ConstructLOWQ(1); // нижний
        auto LQBox_log_nsys2 = ConstructLOWQ(2); // верхний

        // PLACE LQ


        auto rmx = new G4RotationMatrix();
        rmx->rotateZ(180. * deg);
#if defined(RUN21)
            G4double x_pos_LQ1 = 0.0; // нижний
            G4double y_pos_LQ1 = 32.2 * cm;
            G4double z_pos_LQ1 = 63.6 * cm;

            G4double x_pos_LQ2 = 0.0; // верхний
            G4double y_pos_LQ2 = 32.2 * cm;
            G4double z_pos_LQ2 = 63.6 * cm;

            rmx->rotateX(-60. * deg);
#endif

#if defined(RUN23)

        G4double x_pos_LQ1 = 0.0; // нижний
        G4double y_pos_LQ1 = 35.6 * cm;

        G4double converter_th = 0.0 * cm;
#ifdef LOWQ_CONVERTOR
        converter_th = 1.8 * cm;
#endif


        G4double z_pos_LQ1 = 70.7 * cm - converter_th / 2; // учесть сдвиг/толщину конвертера

        G4double x_pos_LQ2 = 0.0; // верхний
        G4double y_pos_LQ2 = 31.9 * cm;
        G4double z_pos_LQ2 = 65.7 * cm - converter_th / 2; // учесть сдвиг/толщину конвертера
        // rmx->rotateX(-57. * deg); // 90 - 33
        rmx->rotateX(-63.2 * deg); // 90 - (180-153.2)
#endif


#ifdef LOWQ1
        auto vol_phys_LQnsys1 = new G4PVPlacement(rmx,
                                     G4ThreeVector(x_pos_LQ1, -y_pos_LQ1, z_pos_LQ1), LQBox_log_nsys1,
                                                  "LQ_phys", fWorldLV, false, 0, fCheckOverlaps);
#endif
#ifdef LOWQ2
        rmx = new G4RotationMatrix();
#if defined(RUN21)
            rmx->rotateX(-60. * deg); // 90-30
#endif

#if defined(RUN23)
        //rmx->rotateX(-57. * deg); //90-33
        rmx->rotateX(-64.1 * deg); // 90 - (180-154.1)
#endif
        auto vol_phys_LQnsys2 = new G4PVPlacement(rmx,
                                                  G4ThreeVector(x_pos_LQ2, y_pos_LQ2, z_pos_LQ2), LQBox_log_nsys2,
                                                  "LQ_phys", fWorldLV, false, 0, fCheckOverlaps);
#endif

#endif // LOWQ1 or LOWQ2


#ifdef TARGET
        ConstructTarget();
#endif // TARGET

        //МАГНИТЫ

#ifdef TARGET
        //////////////////////////////////////////////////////////////////////
        // 		Magnet
        //////////////////////////////////////////////////////////////////////
#ifdef MAGNET_CELL
        ConstructMagnetCell();
#endif

#endif // TARGET
#ifdef MAGNET_D2

        ConstructMagnetD2();

#endif // MAGNET_D2

#ifdef GEM
        ConstructGEM();
#endif

#ifdef MAGNET_D1
        ConstructMagnetD1();
#endif

#ifdef MAGNET_D3
        ConstructMagnetD3();
#endif

#ifdef MAGNET_QF
        ConstructMagnetQF();
#endif

#ifdef MAGNET_QD
        ConstructMagnetQD();
#endif


        // MRPC

#ifdef MRPC1
        ConstructMRPC();
#endif

// БЕТОН

#ifdef CONCRETE


///// БАЛКИ
        G4double betonSizeX = 5. * m, betonSizeY = 1. * m, betonSizeZ = 5. * m;
        G4double betonPositionX = 0, betonPositionY = 3.5 * m, betonPositionZ = 0;

        auto betonS
                = new G4Box("beton",     // its name
                            betonSizeX / 2, betonSizeY / 2, betonSizeZ / 2); // its size

        auto betonLV
                = new G4LogicalVolume(
                        betonS,     // its solid
                        betonMaterial,  // its material
                        "betonLV");   // its name

        new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(betonPositionX, betonPositionY, betonPositionZ),  // at (0,0,0)
                betonLV,          // its logical volume
                "beton",    // its name
                          fWorldLV, // its mother  volume
                          false,            // no boolean operation
                0,                // copy number
                fCheckOverlaps);  // checking overlaps




        //////// КРЫША

        G4double beton_roofSizeX = 5. * m, beton_roofSizeY = 20. * cm, beton_roofSizeZ = 5. * m;
        G4double beton_roofPositionX = 0, beton_roofPositionY = 12.0 * m + 3.5 * m, beton_roofPositionZ = 0;

        auto beton_roofS
                = new G4Box("beton_roof",     // its name
                            beton_roofSizeX / 2, beton_roofSizeY / 2, beton_roofSizeZ / 2); // its size

        auto beton_roofLV
                = new G4LogicalVolume(
                        beton_roofS,     // its solid
                        betonMaterial,  // its material
                        "beton_roofLV");   // its name

        new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(beton_roofPositionX, beton_roofPositionY, beton_roofPositionZ),  // at (0,0,0)
                beton_roofLV,          // its logical volume
                "beton_roof",    // its name
                          fWorldLV, // its mother  volume
                          false,            // no boolean operation
                0,                // copy number
                fCheckOverlaps);  // checking overlaps

#endif //CONCRETE
/*
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
        fWorldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
        fWorldLV_rot->SetVisAttributes(G4VisAttributes::GetInvisible());

        auto simpleBoxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
        simpleBoxVisAtt->SetVisibility(true);
        // calorLV->SetVisAttributes(simpleBoxVisAtt);

#ifdef PF1_FAT
        plastic_fat_nsys1LV->SetVisAttributes(Plastic_VisAtt);
#endif
#ifdef PF2_FAT
        plastic_fat_nsys2LV_120->SetVisAttributes(Plastic_VisAtt);
        plastic_fat_nsys2LV_125->SetVisAttributes(Plastic_VisAtt);
#endif
        //plastic_thin_nsys1LV[0]->SetVisAttributes(simpleBoxVisAtt);
        // plastic_thin_nsys2LV[0]->SetVisAttributes(simpleBoxVisAtt);

        //plastic_thin_nsys1LV[0]->SetVisAttributes(simpleBoxVisAtt);
        // plastic_thin_nsys2LV[1]->SetVisAttributes(simpleBoxVisAtt);
#ifdef CONCRETE
        betonLV->SetVisAttributes(simpleBoxVisAtt);
#endif //CONCRETE
        //
        // Always return the physical World
        //


        return worldPV;
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void DetectorConstruction::ConstructPlasticFat1() {

        ///////////////////////////////////
// Neutron/Proton E/dE counters #1 (below the beam)

#ifdef PF1_FAT


        auto plasticMaterial = G4Material::GetMaterial("Scintillator");
        auto betonMaterial = G4Material::GetMaterial("G4_CONCRETE");
        auto AirMaterial = G4Material::GetMaterial("Air");
        auto MylarMaterial = G4Material::GetMaterial("Mylar");

        //  fNofLayers_plastic_fat_nsys1 = 6;

        G4double plasticFatNsys1SizeX = 1060. * mm, plasticFatNsys1SizeY = 200. * mm, plasticFatNsys1SizeZ = 200. * mm;


        G4double plasticFatNsys1PositionX = 0.;
        G4double plasticFatNsys1PositionY = -102.4 * cm;
        G4double plasticFatNsys1dz = plasticFatNsys1SizeZ + 0.4 * cm;
        // G4double plasticFatNsys1PositionZ_initial = 43.8 * cm - 3.5 * plasticFatNsys1dz; // RIA
        // G4double plasticFatNsys1PositionZ_initial = 6.6 * cm; //Gauzshtein
        G4double plasticFatNsys1PositionZ_initial = -4.5 * cm; //Yurchenko
//        G4double plasticFatNsys1PositionZ_final =
//                plasticFatNsys1PositionZ_initial + fNofLayers_plastic_fat_nsys1 * (plasticFatNsys1dz + 0.01 * cm);

        G4double plasticFatNsys1PositionZ_final = 115.5 * cm; // Yurchenko
        //G4double plasticFatNsys1PositionZ = plasticFatNsys1PositionZ_initial +
        //                                    (plasticFatNsys1PositionZ_final - plasticFatNsys1PositionZ_initial) / 2.;

        G4double plasticFatNsys1PositionZ = plasticFatNsys1PositionZ_initial +
                                            (plasticFatNsys1PositionZ_final - plasticFatNsys1PositionZ_initial) / 2.;

        // Это объем всех шести счетчиков с пленкой
        auto plastic_fat_nsys1_boxallS = new G4Box("plastic_fat_nsys1_boxS",
                                                   1.0 * (plasticFatNsys1SizeX + 0.5 * cm) / 2.,
                                                   1.0 * (plasticFatNsys1SizeY + 0.5 * cm) / 2.,
                                                   fNofLayers_plastic_fat_nsys1 * (plasticFatNsys1SizeZ + 0.4 * cm) /
                                                   2.); // точно равно plasticFatNsys1dz, иначе программа зависает при большом счете
        auto plastic_fat_nsys1_boxallLV = new G4LogicalVolume(plastic_fat_nsys1_boxallS, AirMaterial,
                                                              "plastic_fat_nsys1_boxallLV");
        plastic_fat_nsys1_boxallLV->SetVisAttributes(G4VisAttributes::GetInvisible());

        auto rmx = new G4RotationMatrix();
        rmx->rotateX(180. * deg);
        auto plastic_fat_nsys1_boxallPV
                = new G4PVPlacement(
                        rmx,                // no rotation
                        G4ThreeVector(plasticFatNsys1PositionX, plasticFatNsys1PositionY,
                                      plasticFatNsys1PositionZ),  // at (0,0,0)
                        plastic_fat_nsys1_boxallLV,          // its logical volume
                        "plastic_fat_nsys1_boxallPV",    // its name
                              fWorldLV, // its mother  volume
                              false,            // no boolean operation
                        0,                // copy number
                        fCheckOverlaps);  // checking overlaps

        // Это объем счетчик + пленка
        auto plastic_fat_nsys1_boxS = new G4Box("plastic_fat_nsys1_boxS", (plasticFatNsys1SizeX + 0.2 * cm) / 2.,
                                                (plasticFatNsys1SizeY + 0.2 * cm) / 2.,
                                                (plasticFatNsys1SizeZ + 0.2 * cm) / 2.);
        auto plastic_fat_nsys1_boxLV = new G4LogicalVolume(plastic_fat_nsys1_boxS, AirMaterial,
                                                           "plastic_fat_nsys1_boxLV");
        plastic_fat_nsys1_boxLV->SetVisAttributes(G4VisAttributes::GetInvisible());

        //  Placement РАЗМЕЩЕНИЕ

        auto plastic_fat_nsys1_boxPV =
                new G4PVReplica(
                        "plastic_fat_nsys1_boxPV",          // its name
                        plastic_fat_nsys1_boxLV,          // its logical volume
                        plastic_fat_nsys1_boxallLV,          // its mother
                        kZAxis,           // axis of replication
                        fNofLayers_plastic_fat_nsys1,        // number of replica
                        plasticFatNsys1dz);  // witdth of replic

        // Это объем самого счетчика
        auto plastic_fat_nsys1S
                = new G4Box("plastic_fat_nsys1S",     // its name
                            plasticFatNsys1SizeX / 2., plasticFatNsys1SizeY / 2.,
                            plasticFatNsys1SizeZ / 2.); // its size

        plastic_fat_nsys1LV
                = new G4LogicalVolume(
                plastic_fat_nsys1S,     // its solid
                plasticMaterial,  // its material
                "plastic_fat_nsys1LV");   // its name

        auto plastic_fat_nsys1PV
                = new G4PVPlacement(
                        0,                // no rotation
                        G4ThreeVector(0.0, 0.0, 0.0),  // at (0,0,0)
                        plastic_fat_nsys1LV,          // its logical volume
                        "plastic_fat_nsys1PV",    // its name
                        plastic_fat_nsys1_boxLV,          // its mother  volume
                        false,            // no boolean operation
                        0,                // copy number
                        fCheckOverlaps);  // checking overlaps


        // Это объем пленки

        //  G4VisAttributes *ProCover_VisAtt = new G4VisAttributes(blackpaper_col);
        auto plastic_fat_nsys1_coverS = new G4Box("plastic_fat_nsys1_coverS", plasticFatNsys1SizeX / 2., 0.15 / 2. * mm,
                                                  plasticFatNsys1SizeZ / 2.);
        auto plastic_fat_nsys1_coverLV = new G4LogicalVolume(plastic_fat_nsys1_coverS, MylarMaterial,
                                                             "plastic_fat_nsys1_coverLV");
        plastic_fat_nsys1_coverLV->SetVisAttributes(ProCover_VisAtt);
        // plastic_fat_nsys1_coverLV->SetVisAttributes(new G4VisAttributes(G4Color(0.1, 0.1, 0.5)));

        new G4PVPlacement(0, G4ThreeVector(0.0, +plasticFatNsys1SizeY / 2. + 0.1 * mm, 0.0), plastic_fat_nsys1_coverLV,
                          "plastic_fat_nsys1_coverPV", plastic_fat_nsys1_boxLV, false, 0, fCheckOverlaps);
        new G4PVPlacement(0, G4ThreeVector(0.0, -plasticFatNsys1SizeY / 2. - 0.1 * mm, 0.0), plastic_fat_nsys1_coverLV,
                          "plastic_fat_nsys1_coverPV", plastic_fat_nsys1_boxLV, false, 0, fCheckOverlaps);

        auto rmx1 = new G4RotationMatrix();
        rmx1->rotateX(90. * deg);
        new G4PVPlacement(rmx1, G4ThreeVector(0.0, 0.0, +plasticFatNsys1SizeZ / 2. + 0.1 * mm),
                          plastic_fat_nsys1_coverLV, "plastic_fat_nsys1_coverPV", plastic_fat_nsys1_boxLV, false, 0,
                          fCheckOverlaps);
        new G4PVPlacement(rmx1, G4ThreeVector(0.0, 0.0, -plasticFatNsys1SizeZ / 2. - 0.1 * mm),
                          plastic_fat_nsys1_coverLV, "plastic_fat_nsys1_coverPV", plastic_fat_nsys1_boxLV, false, 0,
                          fCheckOverlaps);



        /*
        for (int i = 1; i <= 6; i++) {
            G4String Phname;
           // Phname = Form("plastic_fat_nsys1_%i_PV", i);
            Phname = "plastic_fat_nsys1PV";
            plastic_fat_nsys1PV = new G4PVPlacement(0, G4ThreeVector( plasticFatNsys1PositionX,  plasticFatNsys1PositionY,  plasticFatNsys1PositionZ),plastic_fat_nsys1_boxLV, Phname, fWorldLV, false, ARM1_IND + DE_IND + i -1);
            plasticFatNsys1PositionZ +=   plasticFatNsys1dz;
        }

    */

#endif // #ifdef PF1
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void DetectorConstruction::ConstructPlasticFat2() {

        // ТОЛСТЫЕ ПЛАСТИКИ

        // Neutron/Proton E/dE counters #2   (above the beam)
#ifdef PF2_FAT


        auto plasticMaterial = G4Material::GetMaterial("Scintillator");
        auto betonMaterial = G4Material::GetMaterial("G4_CONCRETE");
        auto AirMaterial = G4Material::GetMaterial("Air");
        auto MylarMaterial = G4Material::GetMaterial("Mylar");

        // fNofLayers_plastic_fat_nsys2 = 8;

        G4double plasticFatNsys2SizeX_120 = 1000. * mm, plasticFatNsys2SizeY_120 = 120. * mm, plasticFatNsys2SizeZ_120 =
                120. * mm;
        G4double plasticFatNsys2SizeX_125 = 1000. * mm, plasticFatNsys2SizeY_125 = 120. * mm, plasticFatNsys2SizeZ_125 =
                125. * mm;

        G4double plasticFatNsys2PositionX = 0.;
        G4double plasticFatNsys2PositionY = +92. * cm;
        G4double plasticFatNsys2dz = plasticFatNsys2SizeZ_125 + 1.0 * cm;
        //  G4double plasticFatNsys2PositionZ_initial = 36.3 * cm - 3.5 * plasticFatNsys2dz; // RIA
        // G4double plasticFatNsys2PositionZ_initial = 6.2 * cm; //Gauzshtein
        G4double plasticFatNsys2PositionZ_initial =
                0.1 * cm; //Yurchenko, последняя геометрия (некоторые счетчики 125 и 120)
        // G4double plasticFatNsys2PositionZ_final =
        //         plasticFatNsys2PositionZ_initial + fNofLayers_plastic_fat_nsys2 * plasticFatNsys2dz;
        G4double plasticFatNsys2PositionZ_final = 101.5 * cm;

        G4double plasticFatNsys2PositionZ = plasticFatNsys2PositionZ_initial +
                                            (plasticFatNsys2PositionZ_final - plasticFatNsys2PositionZ_initial) / 2.;

        // Это объем всех восьми счетчиков с пленкой
        auto plastic_fat_nsys2_boxallS = new G4Box("plastic_fat_nsys2_boxallS",
                                                   1.0 * (plasticFatNsys2SizeX_125 + 1.0 * cm) / 2.,
                                                   1.0 * (plasticFatNsys2SizeY_125 + 1.0 * cm) / 2.,
                                                   fNofLayers_plastic_fat_nsys2 *
                                                   (plasticFatNsys2SizeZ_125 + 1.0 * cm) /
                                                   2.);
        auto plastic_fat_nsys2_boxallLV = new G4LogicalVolume(plastic_fat_nsys2_boxallS, AirMaterial,
                                                              "plastic_fat_nsys2_boxallLV");
        plastic_fat_nsys2_boxallLV->SetVisAttributes(G4VisAttributes::GetInvisible());

        auto rmx = new G4RotationMatrix();
        rmx->rotateX(180. * deg);

        auto plastic_fat_nsys2_boxallPV
                = new G4PVPlacement(
                        0,                // no rotation
                        G4ThreeVector(plasticFatNsys2PositionX, plasticFatNsys2PositionY,
                                      plasticFatNsys2PositionZ),  // at (0,0,0)
                        plastic_fat_nsys2_boxallLV,          // its logical volume
                        "plastic_fat_nsys2_boxallPV",    // its name
                              fWorldLV, // its mother  volume
                              false,            // no boolean operation
                        0,                // copy number
                        fCheckOverlaps);  // checking overlaps



        // Это объем счетчик + пленка

        // для счетчиков 120
        auto plastic_fat_nsys2_boxS_120 = new G4Box("plastic_fat_nsys2_boxS",
                                                    (plasticFatNsys2SizeX_120 + 0.3 * cm) / 2.,
                                                    (plasticFatNsys2SizeY_120 + 0.3 * cm) / 2.,
                                                    (plasticFatNsys2SizeZ_120 + 0.3 * cm) / 2.);
        auto plastic_fat_nsys2_boxLV_120 = new G4LogicalVolume(plastic_fat_nsys2_boxS_120, AirMaterial,
                                                               "plastic_fat_nsys2_boxLV");
        plastic_fat_nsys2_boxLV_120->SetVisAttributes(G4VisAttributes::GetInvisible());


        // Это объем счетчик + пленка

        // для счетчиков 125
        auto plastic_fat_nsys2_boxS_125 = new G4Box("plastic_fat_nsys2_boxS",
                                                    (plasticFatNsys2SizeX_125 + 0.3 * cm) / 2.,
                                                    (plasticFatNsys2SizeY_125 + 0.3 * cm) / 2.,
                                                    (plasticFatNsys2SizeZ_125 + 0.3 * cm) / 2.);
        auto plastic_fat_nsys2_boxLV_125 = new G4LogicalVolume(plastic_fat_nsys2_boxS_125, AirMaterial,
                                                               "plastic_fat_nsys2_boxLV");
        plastic_fat_nsys2_boxLV_125->SetVisAttributes(G4VisAttributes::GetInvisible());



//        auto plastic_fat_nsys2_boxPV = new G4PVReplica(
//                "plastic_fat_nsys2_boxPV",          // its name
//                plastic_fat_nsys2_boxLV_120,          // its logical volume
//                plastic_fat_nsys2_boxallLV,          // its mother
//                kZAxis,           // axis of replication
//                fNofLayers_plastic_fat_nsys2,        // number of replica
//                plasticFatNsys2dz);  // witdth of replica


        // размещение

        // auto plastic_fat_nsys2_boxPV = new G4PVPlacement(0, G4ThreeVector( plasticFatNsys2PositionX,
        // plasticFatNsys2PositionY,  plasticFatNsys2PositionZ),plastic_fat_nsys2_boxLV, Phname, fWorldLV, false,
        // ARM2_IND + DE_IND + i -1, fCheckOverlaps);

        G4double dz_temp = 0.0 * mm;

        #ifdef RUN21
        G4double dz_temp2 = 0.0 * cm;
        #endif

        #ifdef RUN23
        G4double dz_temp2 = -8.0 * cm;
        #endif


        G4double plasticFatNsys2PositionZ_temp =
                -(plasticFatNsys2PositionZ_final - plasticFatNsys2PositionZ_initial) / 2. + 125.0 * mm / 2. + 1.0 * mm +
                dz_temp + dz_temp2;


//////// РАЗМЕЩЕНИЕ

// счетчик №8
        auto plastic_fat_nsys2_boxPV
                = new G4PVPlacement(
                        0,                // no rotation
                        G4ThreeVector(0, 0,
                                      plasticFatNsys2PositionZ_temp),  // at (0,0,0)
                        plastic_fat_nsys2_boxLV_125,          // its logical volume
                        "plastic_fat_nsys2_boxPV",    // its name
                        plastic_fat_nsys2_boxallLV,          // its mother  volume
                        false,            // no boolean operation
                        8 - 1,                // copy number
                        fCheckOverlaps);  // checking overlaps


// счетчик №7
        plasticFatNsys2PositionZ_temp += 125.0 * mm / 2 + 125.0 * mm / 2 + 3.5 * mm + dz_temp;

        plastic_fat_nsys2_boxPV
                = new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(0, 0,
                              plasticFatNsys2PositionZ_temp),  // at (0,0,0)
                plastic_fat_nsys2_boxLV_125,          // its logical volume
                "plastic_fat_nsys2_boxPV",    // its name
                plastic_fat_nsys2_boxallLV,          // its mother  volume
                false,            // no boolean operation
                7 - 1,                // copy number
                fCheckOverlaps);  // checking overlaps


// счетчик №6

        plasticFatNsys2PositionZ_temp += 120.0 * mm / 2 + 125.0 * mm / 2 + 4.0 * mm + dz_temp;
        plastic_fat_nsys2_boxPV
                = new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(0, 0,
                              plasticFatNsys2PositionZ_temp),  // at (0,0,0)
                plastic_fat_nsys2_boxLV_120,          // its logical volume
                "plastic_fat_nsys2_boxPV",    // its name
                plastic_fat_nsys2_boxallLV,          // its mother  volume
                false,            // no boolean operation
                6 - 1,                // copy number
                fCheckOverlaps);  // checking overlaps


// счетчик №5

        plasticFatNsys2PositionZ_temp += 120.0 * mm / 2 + 120.0 * mm / 2 + 6.0 * mm + dz_temp;
        plastic_fat_nsys2_boxPV
                = new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(0, 0,
                              plasticFatNsys2PositionZ_temp),  // at (0,0,0)
                plastic_fat_nsys2_boxLV_120,          // its logical volume
                "plastic_fat_nsys2_boxPV",    // its name
                plastic_fat_nsys2_boxallLV,          // its mother  volume
                false,            // no boolean operation
                5 - 1,                // copy number
                fCheckOverlaps);  // checking overlaps

// счетчик №4

        plasticFatNsys2PositionZ_temp += 120.0 * mm / 2 + 120.0 * mm / 2 + 7.0 * mm + dz_temp;
        plastic_fat_nsys2_boxPV
                = new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(0, 0,
                              plasticFatNsys2PositionZ_temp),  // at (0,0,0)
                plastic_fat_nsys2_boxLV_120,          // its logical volume
                "plastic_fat_nsys2_boxPV",    // its name
                plastic_fat_nsys2_boxallLV,          // its mother  volume
                false,            // no boolean operation
                4 - 1,                // copy number
                fCheckOverlaps);  // checking overlaps


// счетчик №3

        plasticFatNsys2PositionZ_temp += 125.0 * mm / 2 + 120.0 * mm / 2 + 4.0 * mm + dz_temp;
        plastic_fat_nsys2_boxPV
                = new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(0, 0,
                              plasticFatNsys2PositionZ_temp),  // at (0,0,0)
                plastic_fat_nsys2_boxLV_125,          // its logical volume
                "plastic_fat_nsys2_boxPV",    // its name
                plastic_fat_nsys2_boxallLV,          // its mother  volume
                false,            // no boolean operation
                3 - 1,                // copy number
                fCheckOverlaps);  // checking overlaps

// счетчик №2
        plasticFatNsys2PositionZ_temp += 125.0 * mm / 2 + 125.0 * mm / 2 + 3.5 * mm + dz_temp;
        plastic_fat_nsys2_boxPV
                = new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(0, 0,
                              plasticFatNsys2PositionZ_temp),  // at (0,0,0)
                plastic_fat_nsys2_boxLV_125,          // its logical volume
                "plastic_fat_nsys2_boxPV",    // its name
                plastic_fat_nsys2_boxallLV,          // its mother  volume
                false,            // no boolean operation
                2 - 1,                // copy number
                fCheckOverlaps);  // checking overlaps

// счетчик №1
        plasticFatNsys2PositionZ_temp += 125.0 * mm / 2 + 125.0 * mm / 2 + 3.0 * mm + dz_temp;
        plastic_fat_nsys2_boxPV
                = new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(0, 0,
                              plasticFatNsys2PositionZ_temp),  // at (0,0,0)
                plastic_fat_nsys2_boxLV_125,          // its logical volume
                "plastic_fat_nsys2_boxPV",    // its name
                plastic_fat_nsys2_boxallLV,          // its mother  volume
                false,            // no boolean operation
                1 - 1,                // copy number
                fCheckOverlaps);  // checking overlaps

////////// КОНЕЦ РАЗМЕЩЕНИЯ


// для счетчиков 120
        // Это объем с самим пластиком
        auto plastic_fat_nsys2S_120
                = new G4Box("plastic_fat_nsys2S",     // its name
                            plasticFatNsys2SizeX_120 / 2., plasticFatNsys2SizeY_120 / 2.,
                            plasticFatNsys2SizeZ_120 / 2.); // its size

        plastic_fat_nsys2LV_120
                = new G4LogicalVolume(
                plastic_fat_nsys2S_120,     // its solid
                plasticMaterial,  // its material
                "plastic_fat_nsys2LV");   // its name

        auto plastic_fat_nsys2PV_120
                = new G4PVPlacement(
                        0,                // no rotation
                        G4ThreeVector(0.0, 0.0, 0.0),  // at (0,0,0)
                        plastic_fat_nsys2LV_120,          // its logical volume
                        "plastic_fat_nsys2PV",    // its name
                        plastic_fat_nsys2_boxLV_120,          // its mother  volume
                        false,            // no boolean operation
                        0,                // copy number
                        fCheckOverlaps);  // checking overlaps


// для счетчиков 125
// Это объем с самим пластиком
        auto plastic_fat_nsys2S_125
                = new G4Box("plastic_fat_nsys2S",     // its name
                            plasticFatNsys2SizeX_125 / 2., plasticFatNsys2SizeY_125 / 2.,
                            plasticFatNsys2SizeZ_125 / 2.); // its size

        plastic_fat_nsys2LV_125
                = new G4LogicalVolume(
                plastic_fat_nsys2S_125,     // its solid
                plasticMaterial,  // its material
                "plastic_fat_nsys2LV");   // its name

        auto plastic_fat_nsys2PV_125
                = new G4PVPlacement(
                        0,                // no rotation
                        G4ThreeVector(0.0, 0.0, 0.0),  // at (0,0,0)
                        plastic_fat_nsys2LV_125,          // its logical volume
                        "plastic_fat_nsys2PV",    // its name
                        plastic_fat_nsys2_boxLV_125,          // its mother  volume
                        false,            // no boolean operation
                        0,                // copy number
                        fCheckOverlaps);  // checking overlaps




// пленка для счетчиков 120

        //  G4VisAttributes *ProCover_VisAtt = new G4VisAttributes(blackpaper_col);
        auto plastic_fat_nsys2_coverS_120 = new G4Box("plastic_fat_nsys2_coverS", plasticFatNsys2SizeX_120 / 2.,
                                                      0.15 / 2. * mm,
                                                      plasticFatNsys2SizeZ_120 / 2.);
        auto plastic_fat_nsys2_coverLV_120 = new G4LogicalVolume(plastic_fat_nsys2_coverS_120, MylarMaterial,
                                                                 "plastic_fat_nsys2_coverLV");
        plastic_fat_nsys2_coverLV_120->SetVisAttributes(ProCover_VisAtt);
        // plastic_fat_nsys2_coverLV->SetVisAttributes(new G4VisAttributes(G4Color(0.1, 0.1, 0.5)));

        new G4PVPlacement(0, G4ThreeVector(0.0, +plasticFatNsys2SizeY_120 / 2. + 0.1 * mm, 0.0),
                          plastic_fat_nsys2_coverLV_120,
                          "plastic_fat_nsys2_coverPV", plastic_fat_nsys2_boxLV_120, false, 0, fCheckOverlaps);
        new G4PVPlacement(0, G4ThreeVector(0.0, -plasticFatNsys2SizeY_120 / 2. - 0.1 * mm, 0.0),
                          plastic_fat_nsys2_coverLV_120,
                          "plastic_fat_nsys2_coverPV", plastic_fat_nsys2_boxLV_120, false, 0, fCheckOverlaps);

        auto rmx2 = new G4RotationMatrix();
        rmx2->rotateX(90. * deg);
        new G4PVPlacement(rmx2, G4ThreeVector(0.0, 0.0, +plasticFatNsys2SizeZ_120 / 2. + 0.1 * mm),
                          plastic_fat_nsys2_coverLV_120, "plastic_fat_nsys2_coverPV", plastic_fat_nsys2_boxLV_120,
                          false, 0,
                          fCheckOverlaps);
        new G4PVPlacement(rmx2, G4ThreeVector(0.0, 0.0, -plasticFatNsys2SizeZ_120 / 2. - 0.1 * mm),
                          plastic_fat_nsys2_coverLV_120, "plastic_fat_nsys2_coverPV", plastic_fat_nsys2_boxLV_120,
                          false, 0,
                          fCheckOverlaps);



        // пленка для счетчиков 125

        //  G4VisAttributes *ProCover_VisAtt = new G4VisAttributes(blackpaper_col);
        auto plastic_fat_nsys2_coverS_125 = new G4Box("plastic_fat_nsys2_coverS", plasticFatNsys2SizeX_125 / 2.,
                                                      0.15 / 2. * mm,
                                                      plasticFatNsys2SizeZ_125 / 2.);
        auto plastic_fat_nsys2_coverLV_125 = new G4LogicalVolume(plastic_fat_nsys2_coverS_125, MylarMaterial,
                                                                 "plastic_fat_nsys2_coverLV");
        plastic_fat_nsys2_coverLV_125->SetVisAttributes(ProCover_VisAtt);
        // plastic_fat_nsys2_coverLV->SetVisAttributes(new G4VisAttributes(G4Color(0.1, 0.1, 0.5)));

        new G4PVPlacement(0, G4ThreeVector(0.0, +plasticFatNsys2SizeY_125 / 2. + 0.1 * mm, 0.0),
                          plastic_fat_nsys2_coverLV_125,
                          "plastic_fat_nsys2_coverPV", plastic_fat_nsys2_boxLV_125, false, 0, fCheckOverlaps);
        new G4PVPlacement(0, G4ThreeVector(0.0, -plasticFatNsys2SizeY_125 / 2. - 0.1 * mm, 0.0),
                          plastic_fat_nsys2_coverLV_125,
                          "plastic_fat_nsys2_coverPV", plastic_fat_nsys2_boxLV_125, false, 0, fCheckOverlaps);

        auto rmx3 = new G4RotationMatrix();
        rmx3->rotateX(90. * deg);

        new G4PVPlacement(rmx3, G4ThreeVector(0.0, 0.0, +plasticFatNsys2SizeZ_125 / 2. + 0.1 * mm),
                          plastic_fat_nsys2_coverLV_125, "plastic_fat_nsys2_coverPV", plastic_fat_nsys2_boxLV_125,
                          false, 0,
                          fCheckOverlaps);
        new G4PVPlacement(rmx3, G4ThreeVector(0.0, 0.0, -plasticFatNsys2SizeZ_125 / 2. - 0.1 * mm),
                          plastic_fat_nsys2_coverLV_125, "plastic_fat_nsys2_coverPV", plastic_fat_nsys2_boxLV_125,
                          false, 0,
                          fCheckOverlaps);


        /////////////

/*
    for (int i = 1; i <= 8; i++) {

        G4String Phname;
       // Phname = Form("plastic_fat_nsys2_%i_PV", i);
        Phname = "plastic_fat_nsys2PV";
        plastic_fat_nsys2PV = new G4PVPlacement(0, G4ThreeVector( plasticFatNsys2PositionX,  plasticFatNsys2PositionY,  plasticFatNsys2PositionZ),plastic_fat_nsys2_boxLV, Phname, fWorldLV, false, ARM2_IND + DE_IND + i -1, fCheckOverlaps);
        plasticFatNsys2PositionZ +=   plasticFatNsys2dz;

    }

*/

#endif // #ifdef PF2


    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void Cosmic::DetectorConstruction::ConstructPlasticThin1() {

        // нижнее ПЛЕЧО
#ifdef PF1_THIN


        auto plasticMaterial = G4Material::GetMaterial("Scintillator");
        auto betonMaterial = G4Material::GetMaterial("G4_CONCRETE");
        auto AirMaterial = G4Material::GetMaterial("Air");
        auto MylarMaterial = G4Material::GetMaterial("Mylar");


#ifdef RUN21
        G4double plasticThinNsys1SizeX = 30. * cm, plasticThinNsys1SizeY = 1. * cm, plasticThinNsys1SizeZ = 50. * cm;
        G4double plasticThinNsys1PositionX = 0, plasticThinNsys1PositionY = -45.4 * cm, plasticThinNsys1PositionZ =
                27.0 * cm;

#endif

#ifdef RUN23
        G4double plasticThinNsys1SizeX = 56. * cm, plasticThinNsys1SizeY = 1. * cm, plasticThinNsys1SizeZ = 82. * cm;
        G4double plasticThinNsys1PositionX = 0, plasticThinNsys1PositionY = -49.8 * cm, plasticThinNsys1PositionZ = 21.1 * cm; // 29.5 * cm;

#endif

#ifdef RUN21
        auto plastic_thin_nsys1_boxallS = new G4Box("plastic_thin_nsys1_boxallS",
                                                    fNofLayers_plastic_thin_nsys1 * (plasticThinNsys1SizeX + 4.0 * cm) /
                                                    2., (plasticThinNsys1SizeY + 4.0 * cm) / 2.,
                                                    (plasticThinNsys1SizeZ + 4.0 * cm) / 2.);

#endif

#ifdef RUN23

            // Это объем 2 счетчиков с пленокй
         auto plastic_thin_nsys1_boxallS = new G4Box("plastic_thin_nsys1_boxallS",
                                                     (plasticThinNsys1SizeX + 4.0 * cm) / 2.,
                                                     fNofLayers_plastic_thin_nsys2 * (plasticThinNsys1SizeY + 4.0 * cm) /
                                                     2., (plasticThinNsys1SizeZ + 4.0 * cm) / 2.);

#endif

        auto plastic_thin_nsys1_boxallLV = new G4LogicalVolume(plastic_thin_nsys1_boxallS, AirMaterial,
                                                               "plastic_thin_nsys1_boxallLV");
        plastic_thin_nsys1_boxallLV->SetVisAttributes(G4VisAttributes::GetInvisible());
        //   plastic_thin_nsys1_boxallLV->SetVisAttributes(ProCover_VisAtt);
        // PLACE AC

        //   auto plastic_thin_nsys1_boxallPV =
        //           new G4PVPlacement(0, G4ThreeVector(plasticThinNsys1PositionX, plasticThinNsys1PositionY,
        //           plasticThinNsys1PositionZ),
        //                            plastic_thin_nsys1_boxallLV, "plastic_thin_nsys1_boxPV", fWorldLV, false,
        //                            ARM1_IND);


        auto plastic_thin_nsys1_boxallPV
                = new G4PVPlacement(
                        0,                // no rotation
                        G4ThreeVector(plasticThinNsys1PositionX, plasticThinNsys1PositionY,
                                      plasticThinNsys1PositionZ),  // at (0,0,0)
                        plastic_thin_nsys1_boxallLV,          // its logical volume
                        "plastic_thin_nsys1_boxallPV",    // its name
                              fWorldLV, // its mother  volume
                              false,            // no boolean operation
                        0,                // copy number
                        fCheckOverlaps);  // checking overlaps

// Это объем счетчик + пленка

        auto plastic_thin_nsys1_boxS = new G4Box("plastic_thin_nsys1_boxS", (plasticThinNsys1SizeX + 0.2 * cm) / 2.,
                                                 (plasticThinNsys1SizeY + 0.2 * cm) / 2.,
                                                 (plasticThinNsys1SizeZ + 0.2 * cm) / 2.);
        auto plastic_thin_nsys1_boxLV = new G4LogicalVolume(plastic_thin_nsys1_boxS, AirMaterial,
                                                            "plastic_thin_nsys1_boxLV");
        plastic_thin_nsys1_boxLV->SetVisAttributes(G4VisAttributes::GetInvisible());
        plastic_thin_nsys1_boxLV->SetVisAttributes(ProCover_VisAtt);

#ifdef RUN21
        auto rmx_thin1 = new G4RotationMatrix();
        rmx_thin1->rotateZ(180. * deg);
        rmx_thin1->rotateX(3. * deg);
        new G4PVPlacement(rmx_thin1, G4ThreeVector(plasticThinNsys1SizeX / 2. + 0.1 * cm, 0.0, 0.0),
                          plastic_thin_nsys1_boxLV, "plastic_thin_nsys1_boxPV", plastic_thin_nsys1_boxallLV, false, 0,
                          fCheckOverlaps);
        new G4PVPlacement(rmx_thin1, G4ThreeVector(-plasticThinNsys1SizeX / 2. - 0.1 * cm, 0.0, 0.0),
                          plastic_thin_nsys1_boxLV, "plastic_thin_nsys1_boxPV", plastic_thin_nsys1_boxallLV, false, 1,
                          fCheckOverlaps);

#endif

#ifdef RUN23

        auto rmx_thin1 = new G4RotationMatrix();
        rmx_thin1->rotateZ(180. * deg);
        rmx_thin1->rotateX(3. * deg);

        //    rmx_thin1->rotateX(0. * deg);

        new G4PVPlacement(rmx_thin1, G4ThreeVector(0.0, -plasticThinNsys1SizeY / 2 - 0.1 * cm, 0.0),
                          plastic_thin_nsys1_boxLV, "plastic_thin_nsys1_boxPV", plastic_thin_nsys1_boxallLV,
                          false, 1, fCheckOverlaps); // здесь нумерация слоев наоборот

        new G4PVPlacement(rmx_thin1, G4ThreeVector(0.0, +plasticThinNsys1SizeY / 2 + 0.1 * cm, 0.0), plastic_thin_nsys1_boxLV,
                          "plastic_thin_nsys1_boxPV", plastic_thin_nsys1_boxallLV, false, 0, fCheckOverlaps);

#endif


        //  G4Box *plastic_thin_nsys1S[2];

               // Это объем самого счетчика
        auto plastic_thin_nsys1S = new G4Box("plastic_thin_nsys1S", plasticThinNsys1SizeX / 2.,
                                             plasticThinNsys1SizeY / 2., plasticThinNsys1SizeZ / 2.);
        //  plastic_thin_nsys1S[1] = new G4Box("plastic_thin_nsys1S_2", plasticThinNsys1SizeX / 2., plasticThinNsys1SizeY / 2., plasticThinNsys1SizeZ / 2.);
        plastic_thin_nsys1LV = new G4LogicalVolume(plastic_thin_nsys1S, plasticMaterial, "plastic_thin_nsys1LV");
        plastic_thin_nsys1LV->SetVisAttributes(Plastic_VisAtt);

        new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), plastic_thin_nsys1LV, "plastic_thin_nsys1PV",
                          plastic_thin_nsys1_boxLV, false, -1, fCheckOverlaps);



        //   G4VisAttributes *ProCover_VisAtt = new G4VisAttributes(blackpaper_col);

        auto plastic_thin_nsys1_coverS = new G4Box("plastic_thin_nsys1_coverS", plasticThinNsys1SizeX / 2.,
                                                   0.15 / 2. * mm, plasticThinNsys1SizeZ / 2.);
        auto plastic_thin_nsys1_coverLV = new G4LogicalVolume(plastic_thin_nsys1_coverS, MylarMaterial,
                                                              "plastic_thin_nsys1_coverLV");
        plastic_thin_nsys1_coverLV->SetVisAttributes(ProCover_VisAtt);

        // plastic_thin_nsys1_coverLV->SetVisAttributes(new G4VisAttributes(G4Color(0.1, 0.1, 0.5)));
        //  new G4PVPlacement(0, G4ThreeVector(0.0, +plasticThinNsys1SizeY / 2. + 0.1 * mm, 0.0),plastic_thin_nsys1_coverLV, "plastic_thin_nsys1_coverPV_1", plastic_thin_nsys1LV[0], false, -1, fCheckOverlaps);
        //  new G4PVPlacement(0, G4ThreeVector(0.0, -plasticThinNsys1SizeZ / 2. - 0.1 * mm, 0.0),plastic_thin_nsys1_coverLV, "plastic_thin_nsys1_coverPV_1", plastic_thin_nsys1LV[0], false, -1, fCheckOverlaps);

        //  new G4PVPlacement(0, G4ThreeVector(0.0, +plasticThinNsys1SizeY / 2. + 0.1 * mm, 0.0),plastic_thin_nsys1_coverLV, "plastic_thin_nsys1_coverPV_2", plastic_thin_nsys1LV[1], false, -1, fCheckOverlaps);
        //  new G4PVPlacement(0, G4ThreeVector(0.0, -plasticThinNsys1SizeZ / 2. - 0.1 * mm, 0.0),plastic_thin_nsys1_coverLV, "plastic_thin_nsys1_coverPV_2", plastic_thin_nsys1LV[1], false, -1, fCheckOverlaps);

        new G4PVPlacement(0, G4ThreeVector(0.0, +plasticThinNsys1SizeY / 2. + 0.1 * mm, 0.0),
                          plastic_thin_nsys1_coverLV, "plastic_thin_nsys1_coverPV", plastic_thin_nsys1_boxLV, false, -1,
                          fCheckOverlaps);
        new G4PVPlacement(0, G4ThreeVector(0.0, -plasticThinNsys1SizeY / 2. - 0.1 * mm, 0.0),
                          plastic_thin_nsys1_coverLV, "plastic_thin_nsys1_coverPV", plastic_thin_nsys1_boxLV, false, -1,
                          fCheckOverlaps);




        //  auto rmx_thin1 = new G4RotationMatrix();
        //   rmx_thin1->rotateZ(180. * deg);
        //  rmx_thin1->rotateX(3. * deg);

        //  new G4PVPlacement(rmx_thin1, G4ThreeVector(plasticThinNsys1SizeX /2. + 0.1* cm, 0.0, 0.0),plastic_thin_nsys1LV, "plastic_thin_nsys1PV_all1", plastic_thin_nsys1_boxallLV, false, AC_IND);
        // new G4PVPlacement(rmx_thin1, G4ThreeVector(-plasticThinNsys1SizeX /2. - 0.1* cm, 0.0, 0.0),plastic_thin_nsys1LV, "plastic_thin_nsys1PV_all2", plastic_thin_nsys1_boxallLV, false, AC_IND + 1);



#endif //ifdef PF_THIN
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void DetectorConstruction::ConstructPlasticThin2() {

        // ТОНКИЕ ПЛАСТИКИ
        // ВЕРХНЕЕ ПЛЕЧО
#ifdef PF2_THIN


        auto plasticMaterial = G4Material::GetMaterial("Scintillator");
        auto betonMaterial = G4Material::GetMaterial("G4_CONCRETE");
        auto AirMaterial = G4Material::GetMaterial("Air");
        auto MylarMaterial = G4Material::GetMaterial("Mylar");

#ifdef RUN22
        G4double plasticThinNsys2SizeX = 56. * cm, plasticThinNsys2SizeY = 1. * cm, plasticThinNsys2SizeZ = 82. * cm;
        G4double plasticThinNsys2PositionX = 0, plasticThinNsys2PositionY = 49.8 * cm, plasticThinNsys2PositionZ =
                29.5 * cm;

#endif

#ifdef RUN23
        G4double plasticThinNsys2SizeX = 56. * cm, plasticThinNsys2SizeY = 1. * cm, plasticThinNsys2SizeZ = 82. * cm;
        G4double plasticThinNsys2PositionX = 0, plasticThinNsys2PositionY = 46.7 * cm, plasticThinNsys2PositionZ = 21.1 * cm; //29.5 * cm;

#endif

        G4double plasticThinNsys2dy = plasticThinNsys2SizeY + 1.0 * cm;

// Это объем 2 счетчиков с пленокй
        auto plastic_thin_nsys2_boxallS = new G4Box("plastic_thin_nsys2_boxallS",
                                                    (plasticThinNsys2SizeX + 4.0 * cm) / 2.,
                                                    fNofLayers_plastic_thin_nsys2 * (plasticThinNsys2SizeY + 4.0 * cm) /
                                                    2., (plasticThinNsys2SizeZ + 4.0 * cm) / 2.);
        auto plastic_thin_nsys2_boxallLV = new G4LogicalVolume(plastic_thin_nsys2_boxallS, AirMaterial,
                                                               "plastic_thin_nsys2_boxallLV");
        plastic_thin_nsys2_boxallLV->SetVisAttributes(G4VisAttributes::GetInvisible());
        //   plastic_thin_nsys2_boxallLV->SetVisAttributes(ProCover_VisAtt);
        // PLACE AC
        //  auto plastic_thin_nsys2_boxallPV =
        //          new G4PVPlacement(0, G4ThreeVector(plasticThinNsys2PositionX, plasticThinNsys2PositionY,
        //          plasticThinNsys2PositionZ),
        //                            plastic_thin_nsys2_boxallLV, "plastic_thin_nsys2_boxPV", fWorldLV, false,
        //                            ARM2_IND);

        auto plastic_thin_nsys2_boxallPV
                = new G4PVPlacement(
                        0,                // no rotation
                        G4ThreeVector(plasticThinNsys2PositionX, plasticThinNsys2PositionY,
                                      plasticThinNsys2PositionZ),  // at (0,0,0)
                        plastic_thin_nsys2_boxallLV,          // its logical volume
                        "plastic_thin_nsys2_boxallPV",    // its name
                              fWorldLV, // its mother  volume
                              false,            // no boolean operation
                        0,                // copy number
                        fCheckOverlaps);  // checking overlaps



        //  new G4PVPlacement(0, G4ThreeVector(0.0, plasticThinNsys2SizeY /2. + 0.1* cm, 0.0),plastic_thin_nsys2LV[0], "plastic_thin_nsys2PV_all1", plastic_thin_nsys2_boxallLV, false, AC_IND);
        // new G4PVPlacement(0, G4ThreeVector(0.0, -plasticThinNsys2SizeY /2 - 0.1* cm, 0.0),plastic_thin_nsys2LV[1], "plastic_thin_nsys2PV_all2", plastic_thin_nsys2_boxallLV, false, AC_IND + 1);

// Это объем счетчик + пленка
        auto plastic_thin_nsys2_boxS = new G4Box("plastic_thin_nsys2_boxS", (plasticThinNsys2SizeX + 0.2 * cm) / 2.,
                                                 (plasticThinNsys2SizeY + 0.2 * cm) / 2.,
                                                 (plasticThinNsys2SizeZ + 0.2 * cm) / 2.);
        auto plastic_thin_nsys2_boxLV = new G4LogicalVolume(plastic_thin_nsys2_boxS, AirMaterial,
                                                            "plastic_thin_nsys2_boxLV");
        plastic_thin_nsys2_boxLV->SetVisAttributes(G4VisAttributes::GetInvisible());
        plastic_thin_nsys2_boxLV->SetVisAttributes(ProCover_VisAtt);


        //  Placement РАЗМЕЩЕНИЕ

        auto Yvector = G4ThreeVector(0.0, -1.0, 0.0);
        auto plastic_thin_nsys2_boxPV =
                //   new G4PVReplica(
                //           "plastic_thin_nsys2_boxPV",          // its name
                //          plastic_thin_nsys2_boxLV,          // its logical volume
                //           plastic_thin_nsys2_boxallLV,          // its mother
                //           kYAxis,           // axis of replication
                //           fNofLayers_plastic_thin_nsys2,        // number of replica
                //           plasticThinNsys2dy);  // witdth of replic

        new G4PVPlacement(0, G4ThreeVector(0.0, -plasticThinNsys2SizeY / 2 - 0.1 * cm, 0.0),
                                  plastic_thin_nsys2_boxLV, "plastic_thin_nsys2_boxPV", plastic_thin_nsys2_boxallLV,
                                  false, 0, fCheckOverlaps);

        new G4PVPlacement(0, G4ThreeVector(0.0, +plasticThinNsys2SizeY / 2 + 0.1 * cm, 0.0), plastic_thin_nsys2_boxLV,
                          "plastic_thin_nsys2_boxPV", plastic_thin_nsys2_boxallLV, false, 1, fCheckOverlaps);


        //  new G4PVPlacement(0, G4ThreeVector(0.0, plasticThinNsys2SizeY /2. + 0.1* cm, 0.0),plastic_thin_nsys2LV, "plastic_thin_nsys2PV_all1", plastic_thin_nsys2_boxallLV, false, 1);
        //   new G4PVPlacement(0, G4ThreeVector(0.0, -plasticThinNsys2SizeY /2 - 0.1* cm, 0.0),plastic_thin_nsys2LV, "plastic_thin_nsys2PV_all2", plastic_thin_nsys2_boxallLV, false, 2);

        // Это объем самого счетчика
        auto plastic_thin_nsys2S = new G4Box("plastic_thin_nsys2S", plasticThinNsys2SizeX / 2.,
                                             plasticThinNsys2SizeY / 2., plasticThinNsys2SizeZ / 2.);
        plastic_thin_nsys2LV = new G4LogicalVolume(plastic_thin_nsys2S, plasticMaterial, "plastic_thin_nsys2LV");
        plastic_thin_nsys2LV->SetVisAttributes(Plastic_VisAtt);
        new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), plastic_thin_nsys2LV, "plastic_thin_nsys2PV",
                          plastic_thin_nsys2_boxLV, false, 0, fCheckOverlaps);


        //   G4VisAttributes *ProCover_VisAtt = new G4VisAttributes(blackpaper_col);

        // это объем пленки
        auto plastic_thin_nsys2_coverS = new G4Box("plastic_thin_nsys2_coverS", plasticThinNsys2SizeX / 2.,
                                                   0.15 / 2. * mm, plasticThinNsys2SizeZ / 2.);
        auto plastic_thin_nsys2_coverLV = new G4LogicalVolume(plastic_thin_nsys2_coverS, MylarMaterial,
                                                              "plastic_thin_nsys2_coverLV");
        plastic_thin_nsys2_coverLV->SetVisAttributes(ProCover_VisAtt);

        // plastic_thin_nsys2_coverLV->SetVisAttributes(new G4VisAttributes(G4Color(0.1, 0.1, 0.5)));
        //  new G4PVPlacement(0, G4ThreeVector(0.0, +plasticThinNsys2SizeY / 2. + 0.1 * mm, 0.0),plastic_thin_nsys2_coverLV, "plastic_thin_nsys2_coverPV_1", plastic_thin_nsys2LV[0], false, -1, fCheckOverlaps);
        //  new G4PVPlacement(0, G4ThreeVector(0.0, -plasticThinNsys2SizeZ / 2. - 0.1 * mm, 0.0),plastic_thin_nsys2_coverLV, "plastic_thin_nsys2_coverPV_1", plastic_thin_nsys2LV[0], false, -1, fCheckOverlaps);

        //  new G4PVPlacement(0, G4ThreeVector(0.0, +plasticThinNsys2SizeY / 2. + 0.1 * mm, 0.0),plastic_thin_nsys2_coverLV, "plastic_thin_nsys2_coverPV_2", plastic_thin_nsys2LV[1], false, -1, fCheckOverlaps);
        //  new G4PVPlacement(0, G4ThreeVector(0.0, -plasticThinNsys2SizeZ / 2. - 0.1 * mm, 0.0),plastic_thin_nsys2_coverLV, "plastic_thin_nsys2_coverPV_2", plastic_thin_nsys2LV[1], false, -1, fCheckOverlaps);

        new G4PVPlacement(0, G4ThreeVector(0.0, +plasticThinNsys2SizeY / 2. + 0.1 * mm, 0.0),
                          plastic_thin_nsys2_coverLV, "plastic_thin_nsys2_coverPV", plastic_thin_nsys2_boxLV, false, -1,
                          fCheckOverlaps);
        new G4PVPlacement(0, G4ThreeVector(0.0, -plasticThinNsys2SizeY / 2. - 0.1 * mm, 0.0),
                          plastic_thin_nsys2_coverLV, "plastic_thin_nsys2_coverPV", plastic_thin_nsys2_boxLV, false, -1,
                          fCheckOverlaps);

        auto rmx3 = new G4RotationMatrix();
        rmx3->rotateX(90. * deg);
        //  new G4PVPlacement(rmx3, G4ThreeVector(0.0, 0.0,+plasticThinNsys2SizeZ / 2. + 0.1 * mm),plastic_thin_nsys2_coverLV, "plastic_thin_nsys2_coverPV", plastic_thin_nsys2_boxLV, false, -1, fCheckOverlaps);
        //  new G4PVPlacement(rmx3, G4ThreeVector(0.0, 0.0, -plasticThinNsys2SizeY / 2. - 0.1 * mm),plastic_thin_nsys2_coverLV, "plastic_thin_nsys2_coverPV", plastic_thin_nsys2_boxLV, false, -1, fCheckOverlaps);


        // new G4PVPlacement(0, G4ThreeVector(0.0, plasticThinNsys2SizeY /2. + 0.1* cm, 0.0),plastic_thin_nsys2LV, "plastic_thin_nsys2PV_all1", plastic_thin_nsys2_boxallLV, false, 1);
        //  new G4PVPlacement(0, G4ThreeVector(0.0, -plasticThinNsys2SizeY /2 - 0.1* cm, 0.0),plastic_thin_nsys2LV, "plastic_thin_nsys2PV_all2", plastic_thin_nsys2_boxallLV, false, 2);


        // G4PVPlacement *plastic_thin_nsys2PV[2];
        // plastic_thin_nsys2PV[0] = new G4PVPlacement(0, G4ThreeVector(0.0, plasticThinNsys2SizeY /2. + 0.1* cm, 0.0),plastic_thin_nsys2LV, "plastic_thin_nsys2PV_all1", plastic_thin_nsys2_boxallLV, false, AC_IND);
        // plastic_thin_nsys2PV[1] = new G4PVPlacement(0, G4ThreeVector(0.0, -plasticThinNsys2SizeY /2 - 0.1* cm, 0.0),plastic_thin_nsys2LV, "plastic_thin_nsys2PV_all2", plastic_thin_nsys2_boxallLV, false, AC_IND + 1);


#endif //ifdef PF_THIN


    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Адронный калориметр
    G4LogicalVolume *DetectorConstruction::ConstructHadronCalorimeter(G4int nsys = 1) {

        auto ScintilMaterial = G4Material::GetMaterial("Scintillator");
        auto betonMaterial = G4Material::GetMaterial("G4_CONCRETE");
        auto AirMaterial = G4Material::GetMaterial("Air");
        auto MylarMaterial = G4Material::GetMaterial("Mylar");
        auto SteelMaterial = G4Material::GetMaterial("Steel");

        G4LogicalVolume *Volume_log;
        G4VPhysicalVolume *vol_phys;

        G4LogicalVolume *scint_HadCalLV;

        G4double x_pos, y_pos, z_pos;
        G4double dx, dy, dz;

        G4int NbOfXBars;    // number of bars for phi
        G4int NbOfZBars;    // number of bars for theta
        G4double ScintSizeX;// 10.*cm;
        G4double ScintSizeZ;// 10.*cm;

        //------------------------------------
        //         HARDRON SANDWICH        //////
        //------------------------------------

        /*
        // FEU
        G4LogicalVolume *lFEU;

        G4Tubs *FEU30_box = new G4Tubs("FEU30", 0.0 * cm, 4.5 * cm, 28.0 / 2.0 * cm, 0., 2.0 * M_PI);
        G4LogicalVolume *FEU30_log = new G4LogicalVolume(FEU30_box, AirMaterial, "FEU30_log", 0, 0, 0);
        FEU30_log->SetVisAttributes(G4VisAttributes::GetInvisible());

        G4Tubs *FEU30_tub = new G4Tubs("FEU30tub", 3.0 * cm, 4.0 * cm, 25.0 / 2.0 * cm, 0., 2.0 * M_PI);
        lFEU = new G4LogicalVolume(FEU30_tub, SteelMaterial, "FEU30t_log", 0, 0, 0);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.1 * cm),lFEU, "t", FEU30_log, false, -1,fCheckOverlaps);
        G4Tubs *FEU30_lid = new G4Tubs("FEU30lid", 0.0, 4.0 * cm, 0.1 / 2.0 * cm, 0., 2.0 * M_PI);
        lFEU = new G4LogicalVolume(FEU30_lid, SteelMaterial, "FEU30l_log", 0, 0, 0);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., +(25.0 / 2 + 0.1) * cm),lFEU, "l", FEU30_log, false, -1,fCheckOverlaps);


        G4Tubs *FEU63_box = new G4Tubs("FEU63", 0.0 * cm, 8.5 * cm, 38.0 / 2.0 * cm, 0., 2.0 * M_PI);
        G4LogicalVolume *FEU63_log = new G4LogicalVolume(FEU63_box, AirMaterial, "FEU63_log", 0, 0, 0);
        FEU63_log->SetVisAttributes(G4VisAttributes::GetInvisible());

        G4Tubs *FEU63_tub = new G4Tubs("FEU63tub", 7.0 * cm, 8.0 * cm, 35.0 / 2.0 * cm, 0., 2.0 * M_PI);
        lFEU = new G4LogicalVolume(FEU63_tub, SteelMaterial, "FEU63t_log", 0, 0, 0);
        lFEU->SetVisAttributes(Iron_VisAtt);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.1 * cm),lFEU, "t", FEU63_log, false, -1,fCheckOverlaps);

        G4Tubs *FEU63_lid = new G4Tubs("FEU63lid", 0.0, 8.0 * cm, 0.1 / 2.0 * cm, 0., 2.0 * M_PI);
        lFEU = new G4LogicalVolume(FEU63_lid, SteelMaterial, "FEU63l_log", 0, 0, 0);
        lFEU->SetVisAttributes(Iron_VisAtt);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., +(35.0 / 2 + 0.1) * cm),lFEU, "l", FEU63_log, false, -1,fCheckOverlaps);

        //   G4double FEUlength=35.*cm;

    */

        ////////////////////////////////////////////////////////////////////////////////
        //------------------------------------------------------------------------------
        /// SANDWICH
        //------------------------------------------------------------------------------
        ////////////////////////////////////////////////////////////////////////////////

        //
        //        Y^
        //         |
        //         +---->Z
        //       /
        //    X|_

        // in this version sandwich is square with all strips (X and Z) are the same
        // and number of X bars equals to number of Z bars
        G4double GapX = 0.5 * mm;
        G4double GapY = 0.25 * mm;

        G4double ScintThickness = 7. * mm;
        ScintSizeX = 79.5 * mm;
        ScintSizeZ = NX_BARS * (ScintSizeX + GapX);

        G4double IronThickness = 16. * mm;
        G4double IronSizeX = 16. * cm;
        G4double IronSizeZ = ScintSizeZ + 2. * 10. * cm;

        G4int NbOfLayers = N_LAYERS;// Hadron Calo layers */
        NbOfXBars = NX_BARS;  // number of bars for phi
        NbOfZBars = NZ_BARS;  // number of bars for theta

        G4double LayerStep = 35. * mm;

        // Compute sizes
        G4double StripSizeX = IronSizeX;//(ScintSizeX+GapX)*2;
        G4double StripSizeZ = IronSizeZ + 2. * GapX;
        G4double StripSizeY = (ScintThickness + GapY) * 2. + IronThickness + 4. * GapY;

        G4double LayerSizeX = StripSizeX * NbOfXBars / 2.;
        G4double LayerSizeY = StripSizeY;
        G4double LayerSizeZ = StripSizeZ;

        G4double SandSizeX = IronSizeZ + 1. * mm;
        G4double SandSizeY = LayerStep * (NbOfLayers + 1.0) + 1. * mm;// no ACC layers
        G4double SandSizeZ = IronSizeZ + 1. * mm;

        G4double VertPos = 150.0 * cm;// 150.*cm; !! now to fron face of SANDW
        G4double HorPos = 73.8 * cm;

        G4Box *ubox = new G4Box("Sand", SandSizeX / 2., SandSizeY / 2., SandSizeZ / 2.);
        G4LogicalVolume *sand_vol = new G4LogicalVolume(ubox, AirMaterial, "Sand", 0, 0, 0);
        sand_vol->SetVisAttributes(G4VisAttributes::GetInvisible());

        ubox = new G4Box("Absor", StripSizeX / 2., IronThickness / 2., IronSizeZ / 2.);
        G4LogicalVolume *absor_vol = new G4LogicalVolume(ubox, SteelMaterial, "Absor", 0, 0, 0);
        absor_vol->SetVisAttributes(Convertor_VisAtt);

        ubox = new G4Box("Scint", ScintSizeX / 2., ScintThickness / 2., ScintSizeZ / 2.);
        scint_HadCalLV = new G4LogicalVolume(ubox, ScintilMaterial, "Scint", 0, 0, 0);
        scint_HadCalLV->SetVisAttributes(Plastic_VisAtt);

        ubox = new G4Box("Strip", StripSizeX / 2., StripSizeY / 2., StripSizeZ / 2.);
        G4LogicalVolume *strip_log = new G4LogicalVolume(ubox, AirMaterial, "Strip", 0, 0, 0);
        strip_log->SetVisAttributes(G4VisAttributes::GetInvisible());

        // Standard layer : iron bar and 2 pairs of scints above and below
        // iron bar
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), absor_vol, "BarI", strip_log, false, 0,
                                     fCheckOverlaps);
        y_pos = (IronThickness + ScintThickness) / 2. + GapY;
        x_pos = (ScintSizeX + GapX) / 2.;
        // 2 scint strips on bottom of  iron bar внизу
        vol_phys = new G4PVPlacement(0, G4ThreeVector(-x_pos, -y_pos, 0.), scint_HadCalLV, "BarS", strip_log, false,
                                     0 * N_UNITS, fCheckOverlaps);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(x_pos, -y_pos, 0.), scint_HadCalLV, "BarS", strip_log, false,
                                     1 * N_UNITS, fCheckOverlaps);
        // 2 scint strips on top of  iron bar вверху
        vol_phys = new G4PVPlacement(0, G4ThreeVector(-x_pos, y_pos, 0.), scint_HadCalLV, "BarS", strip_log, false,
                                     2 * N_UNITS, fCheckOverlaps);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(x_pos, y_pos, 0.), scint_HadCalLV, "BarS", strip_log, false,
                                     3 * N_UNITS, fCheckOverlaps);


        ubox = new G4Box("Layer", LayerSizeX / 2., LayerSizeY / 2., LayerSizeZ / 2.);
        G4LogicalVolume *layer_log = new G4LogicalVolume(ubox, AirMaterial, "Layer", 0, 0, 0);
        layer_log->SetVisAttributes(G4VisAttributes::GetInvisible());

        new G4PVReplica("Bars", strip_log, layer_log, kXAxis, N_UNITS, StripSizeX);

        // half-layer
        dy = IronThickness + GapY + ScintThickness;
        ubox = new G4Box("StripH", StripSizeX / 2., dy / 2., StripSizeZ / 2.);
        G4LogicalVolume *stripH_log = new G4LogicalVolume(ubox, AirMaterial, "StripH", 0, 0, 0);
        stripH_log->SetVisAttributes(G4VisAttributes::GetInvisible());
        //  stripH_log->SetVisAttributes(new G4VisAttributes(G4Color(1.0, 0.0, 0.0)));

        dy = (dy - IronThickness) / 2.;
        y_pos = dy;
        vol_phys =
            new G4PVPlacement(0, G4ThreeVector(0., y_pos, 0), absor_vol, "BarI", stripH_log, false, 0,
                                     fCheckOverlaps);
        y_pos -= IronThickness / 2.;
        y_pos -= GapY;
        y_pos -= ScintThickness / 2.;
        x_pos = (ScintSizeX + GapX) / 2.;
        vol_phys = new G4PVPlacement(0, G4ThreeVector(-x_pos, y_pos, 0.), scint_HadCalLV, "BarS", stripH_log, false, 0,
                                     fCheckOverlaps);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(x_pos, y_pos, 0.), scint_HadCalLV, "BarS", stripH_log, false,
                                     N_UNITS, fCheckOverlaps);

        ubox = new G4Box("LayerH", LayerSizeX / 2., LayerSizeY / 2., LayerSizeZ / 2.);
        G4LogicalVolume *layerH_log = new G4LogicalVolume(ubox, AirMaterial, "LayerH", 0, 0, 0);
        layerH_log->SetVisAttributes(G4VisAttributes::GetInvisible());
        //layerH_log->SetVisAttributes(new G4VisAttributes(G4Color(1.0, 0.0, 0.0)));

        new G4PVReplica("BarsH", stripH_log, layerH_log, kXAxis, N_UNITS, StripSizeX);

        // outer (iron-less) half-layer // первый слой
        ubox = new G4Box("StripO", StripSizeX / 2., StripSizeY / 2., StripSizeZ / 2.);
        G4LogicalVolume *stripO_log = new G4LogicalVolume(ubox, AirMaterial, "StripO", 0, 0, 0);
        stripO_log->SetVisAttributes(G4VisAttributes::GetInvisible());
        //stripO_log->SetVisAttributes(new G4VisAttributes(G4Color(1.0, 0.0, 0.0)));

        y_pos = -StripSizeY / 2.;
        y_pos += ScintThickness / 2.;
        x_pos = (ScintSizeX + GapX) / 2.;
        vol_phys = new G4PVPlacement(0, G4ThreeVector(-x_pos, y_pos, 0.), scint_HadCalLV, "BarS", stripO_log, false, 0,
                                     fCheckOverlaps);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(x_pos, y_pos, 0.), scint_HadCalLV, "BarS", stripO_log, false,
                                     N_UNITS, fCheckOverlaps);
        y_pos += ScintThickness / 2.;
        y_pos += IronThickness / 2.;
        // NO IRON HERE !!
        //  vol_phys = new G4PVPlacement(0, G4ThreeVector(0.,y_pos,0),absor_vol, "BarI", stripH_log, false, 0);

        ubox = new G4Box("LayerO", LayerSizeX / 2., LayerSizeY / 2., LayerSizeZ / 2.);
        G4LogicalVolume *layerO_log = new G4LogicalVolume(ubox, AirMaterial, "LayerO", 0, 0, 0);
        layerO_log->SetVisAttributes(G4VisAttributes::GetInvisible());

        new G4PVReplica("BarsO", stripO_log, layerO_log, kXAxis, N_UNITS, StripSizeX);

        ///
        G4int nx = HCX_IND, nz = HCZ_IND;

        // Размножение по Y

        y_pos = -SandSizeY / 2. + LayerStep / 2.;
        for (G4int i = 0; i < 1 * NbOfLayers + 1; y_pos += LayerStep, i++)
        {
            // нечетные 1,3,5,7,9
                if ((i & 1) == 1)
            {

                //  G4cout << "!!!" << i << G4endl;
                // X-layer
                        if (i != NbOfLayers)
                        {
                                // not last layer?
                                vol_phys = new G4PVPlacement(G4Transform3D(RotateNull, G4ThreeVector(0., y_pos, 0.)),
                                                             layer_log,
                                                             "BarsX", sand_vol, false, nx, fCheckOverlaps);
                                nx += NbOfXBars * 2;
                        }
                        else
                        {
                                // if last layer - use half-layer
                                vol_phys = new G4PVPlacement(G4Transform3D(RotateNull, G4ThreeVector(0., y_pos, 0.)),
                                                             layerO_log,
                                                             "BarsHX", sand_vol, false, nx, fCheckOverlaps);
                                nx += NbOfXBars;
                        }
                }
                else
                    // четные 0,2,4,6,8,10
                {
                        // Z-layer
                        if (i == 0)
                        {
                                // first layer ?
                                vol_phys = new G4PVPlacement(
                                        G4Transform3D(Rotate90Y180Z, G4ThreeVector(0., y_pos + dy - 0.15 * mm, 0.)),
                                        layerH_log,
                                        "BarsHZ", sand_vol, false, nz, fCheckOverlaps);
                                nz += NbOfZBars;
                    continue;
                }
                if (i != NbOfLayers) {// not last layer?
                    vol_phys = new G4PVPlacement(G4Transform3D(Rotate90Y180Z, G4ThreeVector(0., y_pos, 0.)), layer_log,
                                                 "BarsZ", sand_vol, false, nz, fCheckOverlaps);
                    nz += NbOfZBars * 2;
                } else {// if last layer - use half-layer
                    vol_phys = new G4PVPlacement(G4Transform3D(Rotate90Y, G4ThreeVector(0., y_pos, 0.)), layerO_log,
                                                 "BarsOZ", sand_vol, false, nz, fCheckOverlaps);
                    nz += NbOfZBars;
                }
            }

        }

        G4cerr << "***** nx=" << nx - HCX_IND << " nz=" << nz - HCZ_IND;
        G4cerr << " HCX_IND=" << HCX_IND << " HCZ_IND=" << HCZ_IND << " LQ_IND=" << LQ_IND;
        G4cerr << " ARM1_IND=" << ARM1_IND << " ARM2_IND=" << ARM2_IND << G4endl;


        if (nsys == 1) scint_HadCal_nsys1LV = scint_HadCalLV;
        if (nsys == 2) scint_HadCal_nsys2LV = scint_HadCalLV;


        return sand_vol;
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// MRP-детекторы
    void DetectorConstruction::ConstructMRPC_old() {

        // RPC
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
        //BEGIN RPC// for simulated mRPC
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//

        auto ScintilMaterial = G4Material::GetMaterial("Scintillator");

        G4int NbOfXBars;    // number of bars for phi
        G4int NbOfZBars;    // number of bars for theta
        G4double ScintSizeX;// 10.*cm;
        G4double ScintSizeZ;// 10.*cm;

        G4double GapX = 0.5 * mm;
        G4double GapY = 0.25 * mm;

        G4double ScintThickness = 7. * mm;
        ScintSizeX = 79.5 * mm;
        ScintSizeZ = NX_BARS * (ScintSizeX + GapX);

        G4double IronThickness = 16. * mm;
        G4double IronSizeX = 16. * cm;
        G4double IronSizeZ = ScintSizeZ + 2. * 10. * cm;

        G4int NbOfLayers = N_LAYERS;// Hadron Calo layers */
        NbOfXBars = NX_BARS;  // number of bars for phi
        NbOfZBars = NZ_BARS;  // number of bars for theta

        G4double LayerStep = 35. * mm;

        // Compute sizes
        G4double StripSizeX = IronSizeX;//(ScintSizeX+GapX)*2;
        G4double StripSizeZ = IronSizeZ + 2 * GapX;
        G4double StripSizeY = (ScintThickness + GapY) * 2 + IronThickness + 4 * GapY;

        G4double LayerSizeX = StripSizeX * NbOfXBars / 2;
        G4double LayerSizeY = StripSizeY;
        G4double LayerSizeZ = StripSizeZ;

        G4double SandSizeX = IronSizeZ + 1 * mm;
        G4double SandSizeY = LayerStep * (NbOfLayers + 1.0) + 1 * mm;// no ACC layers
        G4double SandSizeZ = IronSizeZ + 1 * mm;

        G4double pl_thick = 1.0 * cm;    // y-axis
        G4double pl_width = SandSizeZ;    // z-axis
        G4double pl_length = SandSizeX;    // x-axis

        G4double VertPos = 150.0 * cm;// 150.*cm; !! now to fron face of SANDW
        G4double HorPos = 73.8 * cm;

//G4Material* rpc_material=Scintil;

        auto ubox = new G4Box("RPC", pl_length / 2., pl_thick / 2., pl_width / 2.);
        G4LogicalVolume *RPC_log = new G4LogicalVolume(ubox, ScintilMaterial, "RPC_log", 0, 0, 0);
        RPC_log->SetVisAttributes(Plastic_VisAtt);

// PLACE RPC
// Position -- close to HS
        G4double x_pos = 0.0;
        G4double y_pos = VertPos - 0.5 * LayerSizeY;
//G4cout<<"*** RPC_Y = "<<y_pos<<G4endl;
        G4double z_pos = y_pos * (tan((90. - 45.) * deg) + tan((90. - 85.) * deg)) / 2.;
//G4cerr << "========> AC : y_pos="<<y_pos/cm<<" cm  z_pos="<<z_pos/cm<<" cm"<<G4endl;

//#ifdef RPC
/*
i=TOF_IND  + ARM1_IND;
new G4PVPlacement(G4Transform3D(RotateNull,
   		  G4ThreeVector(x_pos, -y_pos, z_pos)),
		  RPC_log,
                  "RPC1_phys",
                  logicWorld,
                  false,
                  i);
*/

        G4int i = TOF_IND + ARM2_IND;
        new G4PVPlacement(G4Transform3D(Rotate180Z,
                                        G4ThreeVector(x_pos, y_pos, z_pos)),
                          RPC_log,
                          "RPC2_phys", fWorldLV,
                          false,
                          -1, fCheckOverlaps);
//                  i);

//#endif
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
        //END RPC// for simulated mRPC
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//

    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void DetectorConstruction::ConstructMRPC() {

        // RPC
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
        //BEGIN RPC// for simulated mRPC
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//

        auto ScintilMaterial = G4Material::GetMaterial("Scintillator");

        G4int NbOfXBars;    // number of bars for phi
        G4int NbOfZBars;    // number of bars for theta
        G4double ScintSizeX;// 10.*cm;
        G4double ScintSizeZ;// 10.*cm;

        G4double GapX = 0.5 * mm;
        G4double GapY = 0.25 * mm;

        G4double ScintThickness = 7. * mm;
        ScintSizeX = 79.5 * mm;
        ScintSizeZ = NX_BARS * (ScintSizeX + GapX);

        G4double IronThickness = 16. * mm;
        G4double IronSizeX = 16. * cm;
        G4double IronSizeZ = ScintSizeZ + 2. * 10. * cm;

        G4int NbOfLayers = N_LAYERS;// Hadron Calo layers */
        NbOfXBars = NX_BARS;  // number of bars for phi
        NbOfZBars = NZ_BARS;  // number of bars for theta

        G4double LayerStep = 35. * mm;

        // Compute sizes
        G4double StripSizeX = IronSizeX;//(ScintSizeX+GapX)*2;
        G4double StripSizeZ = IronSizeZ + 2 * GapX;
        G4double StripSizeY = (ScintThickness + GapY) * 2 + IronThickness + 4 * GapY;

        G4double LayerSizeX = StripSizeX * NbOfXBars / 2;
        G4double LayerSizeY = StripSizeY;
        G4double LayerSizeZ = StripSizeZ;

        G4double SandSizeX = IronSizeZ + 1 * mm;
        G4double SandSizeY = LayerStep * (NbOfLayers + 1.0) + 1 * mm;// no ACC layers
        G4double SandSizeZ = IronSizeZ + 1 * mm;

        G4double pl_thick = 1.0 * cm;    // y-axis
        G4double pl_width = SandSizeZ;    // z-axis
        G4double pl_length = SandSizeX;    // x-axis

//G4Material* rpc_material=Scintil;

        auto ubox = new G4Box("RPC", pl_length / 2., pl_thick / 2., pl_width / 2.);
        G4LogicalVolume *RPC_log = new G4LogicalVolume(ubox, ScintilMaterial, "RPC_log", 0, 0, 0);
        RPC_log->SetVisAttributes(Plastic_VisAtt);

        G4double VertPos = 150.0 * cm;// 150.*cm; !! now to fron face of SANDW
        G4double HorPos = 73.8 * cm;

// PLACE RPC
// Position -- close to HS
        G4double x_pos = 0.0;
        G4double y_pos = VertPos - 0.5 * LayerSizeY;
//G4cout<<"*** RPC_Y = "<<y_pos<<G4endl;
        G4double z_pos = y_pos * (tan((90. - 45.) * deg) + tan((90. - 85.) * deg)) / 2.;
//G4cerr << "========> AC : y_pos="<<y_pos/cm<<" cm  z_pos="<<z_pos/cm<<" cm"<<G4endl;



        G4int i = TOF_IND + ARM1_IND;
        new G4PVPlacement(G4Transform3D(RotateNull,
                                        G4ThreeVector(x_pos, -y_pos, z_pos)),
                          RPC_log,
                          "RPC1_phys",
                          fWorldLV, false,
                          i, fCheckOverlaps);


        i = TOF_IND + ARM2_IND;
        new G4PVPlacement(G4Transform3D(Rotate180Z,
                                        G4ThreeVector(x_pos, y_pos, z_pos)),
                          RPC_log,
                          "RPC2_phys", fWorldLV,
                          false,
                          -1, fCheckOverlaps);
        //                  i);


        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
        //END RPC// for simulated mRPC
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Вершинные камеры
    G4LogicalVolume *DetectorConstruction::ConstructVC(G4LogicalVolume *&forSD) {

        auto AirMaterial = G4Material::GetMaterial("Air");
        auto MylarMaterial = G4Material::GetMaterial("Mylar");
        auto SteelMaterial = G4Material::GetMaterial("Steel");
        auto SiO2Material = G4Material::GetMaterial("SiO2");
        auto GasArCO2Material = G4Material::GetMaterial("ArgonCarbonicGas");
        auto GasArIsobutaMaterial = G4Material::GetMaterial("ArgonIsobutaneGas");

       // auto GasVCMaterial = G4Material::GetMaterial("ArgonCarbonicGas");
        auto GasVCMaterial = G4Material::GetMaterial("ArgonIsobutaneGas");



        ////////////////////////////////////////////////////////////////////////////////
// VOLUME : Vertex chamber
////////////////////////////////////////////////////////////////////////////////


        G4Box *vbox =
                new G4Box("VCBox", 12.0 / 2. * cm, 3.0 / 2. * cm, 53.0 / 2. * cm);

        G4LogicalVolume *
                VCBox_log = new G4LogicalVolume(vbox, AirMaterial, "VCBox_log", 0, 0, 0);
        VCBox_log->SetVisAttributes(G4VisAttributes::GetInvisible());

//------------------------------------------------------------------------------

        vbox = new G4Box("frame", 11.4 / 2. * cm, 0.4 / 2. * cm, 52.0 / 2. * cm);

        G4Box *vcwin = new G4Box("framew", 9.4 / 2. * cm, 0.4025 / 2. * cm, 43.6 / 2. * cm);

        G4SubtractionSolid *VCframe =
                new G4SubtractionSolid("VCframe_0", vbox, vcwin,
                                       G4Transform3D(RotateNull, G4ThreeVector(0., 0., 0.)));

        vbox = new G4Box("stef", 114. / 2. * mm, 20.0 / 2. * mm, 48.0 / 2. * cm);

        vcwin = new G4Box("stefw", 94. / 2. * mm, 20.02 / 2. * mm, 43.6 / 2. * cm);

        G4SubtractionSolid *VCStef =
                new G4SubtractionSolid("VCframe_0", vbox, vcwin,
                                       G4Transform3D(RotateNull, G4ThreeVector(0., 0., 0.)));

        G4LogicalVolume *VC_log = new G4LogicalVolume(VCframe, SteelMaterial, "VCFrame_log", 0, 0, 0);
        new G4PVPlacement(0, G4ThreeVector(0., 12.0 * mm, 0.),
                          VC_log, "VCFrame_phys", VCBox_log, false, -1, fCheckOverlaps);

//   G4VisAttributes *VCVisAtt = new G4VisAttributes(steel_col);
        VC_log->SetVisAttributes(Steel_VisAtt);

        VC_log = new G4LogicalVolume(VCStef, SiO2Material, "VCStef_log", 0, 0, 0);
        new G4PVPlacement(0, G4ThreeVector(0., 0.0 * cm, 0.),
                          VC_log, "VCStef_phys", VCBox_log, false, -1, fCheckOverlaps);

//   VCVisAtt = new G4VisAttributes(stef_col);
        VC_log->SetVisAttributes(Stef_VisAtt);

        vbox = new G4Box("VCGas_trap", 2. * NVC_WRS / 2. * mm, 10.0 / 2. * mm, 43.6 / 2. * cm);

        VC_log = new G4LogicalVolume(vbox, GasVCMaterial, "VCGas_log", 0, 0, 0);
        new G4PVPlacement(0, G4ThreeVector(0.0 * cm, 0.0 * cm, 0.0 * cm),
                          VC_log, "VCtrap_phys", VCBox_log, false, VC_IND, fCheckOverlaps);

        vbox = new G4Box("VCGas_cell", 2. / 2. * mm, 10.0 / 2. * mm, 43.6 / 2. * cm);

        G4LogicalVolume *
                VCGas_log = new G4LogicalVolume(vbox, GasVCMaterial, "VCGas_log", 0, 0, 0);
        new G4PVReplica("VCcells", VCGas_log, VC_log, kXAxis, NVC_WRS, 2. * mm);


//   VCVisAtt = new G4VisAttributes(gas_col);
        VCGas_log->SetVisAttributes(Gas_VisAtt);

        vbox = new G4Box("VCMylar_box1", 94.0 / 2. * mm, 0.02 / 2. * mm, 43.6 / 2. * cm);

        VC_log = new G4LogicalVolume(vbox, MylarMaterial, "VCMylar1_log", 0, 0, 0);
        new G4PVPlacement(0, G4ThreeVector(0.0, -5.01 * mm, 0.0),
                          VC_log, "VCMylar1_phys", VCBox_log, false, -1, fCheckOverlaps);
        new G4PVPlacement(0, G4ThreeVector(0.0, 5.01 * mm, 0.0),
                          VC_log, "VCMylar1_phys", VCBox_log, false, -1, fCheckOverlaps);

//   VCVisAtt = new G4VisAttributes(mylar_col);
//   VCVisAtt->SetForceWireframe(true);
        VC_log->SetVisAttributes(Mylar_VisAtt);

        vbox = new G4Box("VCMylar_box1", 94.0 / 2. * mm, 0.07 / 2. * mm, 43.6 / 2. * cm);

        VC_log = new G4LogicalVolume(vbox, MylarMaterial, "VCMylar1_log", 0, 0, 0);
        new G4PVPlacement(0, G4ThreeVector(0.0, -10.0 * mm, 0.0),
                          VC_log, "VCMylar1_phys", VCBox_log, false, -1, fCheckOverlaps);
        new G4PVPlacement(0, G4ThreeVector(0.0, 10.0 * mm, 0.0),
                          VC_log, "VCMylar1_phys", VCBox_log, false, -1, fCheckOverlaps);

        VC_log->SetVisAttributes(Mylar_VisAtt);

        forSD = VCGas_log;
        return VCBox_log;

    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Дрейфовые камеры
    G4LogicalVolume *DetectorConstruction::ConstructWC(G4double Lwin,
                                                       G4double Wwin, G4int ind, G4LogicalVolume *&forSD) {


        auto AirMaterial = G4Material::GetMaterial("Air");
        auto MylarMaterial = G4Material::GetMaterial("Mylar");
        auto SteelMaterial = G4Material::GetMaterial("Steel");
        auto SiO2Material = G4Material::GetMaterial("SiO2");
        auto GasArCO2Material = G4Material::GetMaterial("ArgonCarbonicGas");
        auto GasArIsobutaMaterial = G4Material::GetMaterial("ArgonIsobutaneGas");

        // auto GasWCMaterial = G4Material::GetMaterial("ArgonCarbonicGas");
        auto GasWCMaterial = G4Material::GetMaterial("ArgonIsobutaneGas");



////////////////////////////////////////////////////////////////////////////////
// VOLUME : Tracking chamber
////////////////////////////////////////////////////////////////////////////////
//
//        Y^
//         |
//         +---->Z
//       /
//    X|_


        G4double Lframe, Wframe, Gap, Swid;
        G4double x_pos, y_pos, z_pos;

//Lwin=40.*cm; Wwin=15.*cm;
        Swid = 4. * cm;
        Gap = 6. * mm;
        Lframe = Lwin + Swid * 2;
        Wframe = Wwin + Swid * 2;

        G4double thk = Gap * 8.;
        G4double foil_thk = 0.03 * mm;
        G4double mylar_thk = 0.05 * mm;

        G4Box *WCBox_box =
                new G4Box("WCBox_box", Wframe / 2., thk / 2., Lframe / 2.);

        G4LogicalVolume *
                WCBox_log = new G4LogicalVolume(WCBox_box, GasWCMaterial, "WCBox_log", 0, 0, 0);
        WCBox_log->SetVisAttributes(G4VisAttributes::GetInvisible());

//------------------------------------------------------------------------------

        G4Box *
                WCBox = new G4Box("WC_S1", Swid / 2., thk / 2., Lwin / 2.);
        G4LogicalVolume *
                WC_log = new G4LogicalVolume(WCBox, SiO2Material, "S1", 0, 0, 0);
        WC_log->SetVisAttributes(Stef_VisAtt);

        x_pos = Wwin / 2. + Swid / 2.;
        G4VPhysicalVolume *WCelem;
        WCelem = new G4PVPlacement(0, G4ThreeVector(x_pos, 0., 0.0), WC_log, "WCS1a",
                                   WCBox_log, false, -1, fCheckOverlaps);
        WCelem = new G4PVPlacement(0, G4ThreeVector(-x_pos, 0., 0.0), WC_log, "WCS1b",
                                   WCBox_log, false, -1, fCheckOverlaps);

        x_pos = Wwin + 2 * Swid;
        WCBox = new G4Box("WC_S2", x_pos / 2., thk / 2., Swid / 2.);
        WC_log = new G4LogicalVolume(WCBox, SiO2Material, "S2", 0, 0, 0);
        WC_log->SetVisAttributes(Stef_VisAtt);
        z_pos = Lwin / 2. + Swid / 2.;
        WCelem = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, z_pos), WC_log, "WCS2a",
                                   WCBox_log, false, -1, fCheckOverlaps);
        WCelem = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -z_pos), WC_log, "WCS2b",
                                   WCBox_log, false, -1, fCheckOverlaps);

//----------------------------------------
        G4Box *
                WCFoil_box = new G4Box("WCFoil_box", Wwin / 2., foil_thk, Lwin / 2.);

        G4Box *
                WCMylar_box = new G4Box("WCMylar_box", Wwin / 2., mylar_thk, Lwin / 2.);

        G4LogicalVolume *
                WCFoil = new G4LogicalVolume(WCFoil_box, MylarMaterial, "WCfoil", 0, 0, 0);
        WCFoil->SetVisAttributes(Foil_VisAtt);
        G4LogicalVolume *
                WCMylar = new G4LogicalVolume(WCMylar_box, MylarMaterial, "WCmylar", 0, 0, 0);
        WCMylar->SetVisAttributes(Mylar_VisAtt);

        y_pos = Gap * 4 - 0.05 * mm;
        WCelem = new G4PVPlacement(0, G4ThreeVector(0.0, y_pos, 0.0), WCMylar, "WCM1a",
                                   WCBox_log, false, -1, fCheckOverlaps);
        WCelem = new G4PVPlacement(0, G4ThreeVector(0.0, -y_pos, 0.0), WCMylar, "WCM1b",
                                   WCBox_log, false, -1, fCheckOverlaps);
        y_pos = Gap * 1 - 0.05 * mm;
        WCelem = new G4PVPlacement(0, G4ThreeVector(0.0, y_pos, 0.0), WCFoil, "WCF1a",
                                   WCBox_log, false, -1, fCheckOverlaps);
        WCelem = new G4PVPlacement(0, G4ThreeVector(0.0, -y_pos, 0.0), WCFoil, "WCF1b",
                                   WCBox_log, false, -1, fCheckOverlaps);
        y_pos = Gap * 3 - 0.05 * mm;
        WCelem = new G4PVPlacement(0, G4ThreeVector(0.0, y_pos, 0.0), WCFoil, "WCF2a",
                                   WCBox_log, false, -1, fCheckOverlaps);
        WCelem = new G4PVPlacement(0, G4ThreeVector(0.0, -y_pos, 0.0), WCFoil, "WCF2b",
                                   WCBox_log, false, -1, fCheckOverlaps);
//------------------------------------------------------------------------------
// for sensitive detector -- wire planes
        x_pos = Wwin;//-4.*mm;
        y_pos = 2. * Gap * 0.9;
        z_pos = Lwin;//-4.*mm;
        G4Box *
                WCplane_box = new G4Box("WCplane_box", x_pos / 2., y_pos / 2., z_pos / 2.);
        G4Box *
                WCcell_box = new G4Box("WCcell_box", x_pos / 2., y_pos / 2., 2. * cm / 2.);

        G4LogicalVolume *
                WCplane_log = new G4LogicalVolume(WCplane_box, GasWCMaterial, "WCplane", 0, 0, 0);
        WCplane_log->SetVisAttributes(G4VisAttributes::GetInvisible());

        G4LogicalVolume *
                WCcell_log = new G4LogicalVolume(WCcell_box, GasWCMaterial, "WCcell", 0, 0, 0);
//   WCcell_log->SetVisAttributes(Gas_VisAtt);
        WCcell_log->SetVisAttributes(Steel_VisAtt);

        G4int nr = Lwin / (2. * cm) + 0.1;

        new G4PVReplica("WCcells", WCcell_log, WCplane_log, kZAxis, nr, 2. * cm);

        y_pos = -2. * Gap;
//   WCelem = new G4PVPlacement(0,G4ThreeVector(0.0,y_pos,0.0),WCplane_log,"WCpl1",WCBox_log,false,ind);
        y_pos = 0.;
        WCelem = new G4PVPlacement(0, G4ThreeVector(0.0, y_pos, 0.0), WCplane_log, "WCpl2", WCBox_log, false, ind,
                                   fCheckOverlaps);

        y_pos = 2. * Gap;
//   WCelem = new G4PVPlacement(0,G4ThreeVector(0.0,y_pos,0.0),WCplane_log,"WCpl3",WCBox_log,false,ind);


        forSD = WCcell_log;
        return WCBox_log;

    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructTarget()
    {

#ifdef TARGET
        ////////////////////////////////////////////////////////////////////////////////
        // VOLUME : Main Block - target place
        ////////////////////////////////////////////////////////////////////////////////

        auto SteelMaterial = G4Material::GetMaterial("Steel");
        auto AluminMaterial = G4Material::GetMaterial("Aluminum");
        auto VacuumMaterial = G4Material::GetMaterial("Vacuum");
        auto TitanMaterial = G4Material::GetMaterial("Titanium");
        auto BerylliumMaterial = G4Material::GetMaterial("Beryllium");
        auto CopperMaterial = G4Material::GetMaterial("Copper"); // материал для медной фольги

        G4LogicalVolume *Volume_log;
        G4VPhysicalVolume *vol_phys;

        G4double w_lng = 24.0 * cm; //21.6*cm;		// titanium window length
        G4double w_shft = -10. * cm;

        G4double dx, dz, d_y1, d_y2;

        G4double w_pos = w_lng / 2 + w_shft;


#ifdef CUBIK

        G4Box *VolumeOut_box =
                new G4Box("VolumeOut_box1", 6.3 / 2 * cm, 5.0 / 2 * cm, 52.8 / 2 * cm);

        G4Box *VolumeIn_box =
                new G4Box("VolumeIn_box", (6.3 - 0.6) / 2 * cm, (5.0 - 0.6) / 2 * cm, 52.81 / 2 * cm);

        dx = 52.8 / 2. * cm - (w_pos + (w_lng + 1.4 * cm) / 2.) - 1.0 * cm;
        dz = 52.8 / 2. * cm - 0.5 * cm - dx / 2.;

        d_y1 = 1.0 * cm - 0.4 * mm; // inlet tube outer radius
        d_y2 = d_y1 - 0.2 * cm; // inlet tube inner radius
        G4Tubs *InletHole_tub =
                new G4Tubs("InletHole_tub", 0.0 * cm, d_y1, 0.31 / 2 * cm, 0 * deg, 360 * deg);
        G4Tubs *Inlet_tub =
                new G4Tubs("Inlet_tub", d_y2, d_y1, 20. / 2. * cm, 0 * deg, 360 * deg);

        G4SubtractionSolid *Volume_0 =
                new G4SubtractionSolid("Volume_0", VolumeOut_box, VolumeIn_box);

        G4SubtractionSolid *Volume_1 =
                new G4SubtractionSolid("Volume_1", Volume_0, InletHole_tub,
                                       G4Transform3D(Rotate90Y, G4ThreeVector(3.0 * cm, 0.0 * cm, 0.0 * cm)));

        G4SubtractionSolid *Volume_3;
        Volume_3 = Volume_1;

        G4Box *Cavity_box =
                new G4Box("Cavity_box", 4.0 / 2 * cm, 5.05 / 2 * cm, w_lng / 2);

        G4SubtractionSolid *Volume =
                new G4SubtractionSolid("Volume", Volume_3, Cavity_box,
                                       G4Transform3D(RotateNull, G4ThreeVector(0.0, 0.0, w_pos)));

        Volume_log = new G4LogicalVolume(Volume, SteelMaterial, "Volume_log", 0, 0, 0);

        vol_phys = new G4PVPlacement(0, G4ThreeVector(), Volume_log, "Volume_phys", fWorldLV_rot, false, -1, fCheckOverlaps);

        //  G4VisAttributes *Volume_logVisAtt = new G4VisAttributes(steel_col);
        Volume_log->SetVisAttributes(Steel_VisAtt);

// Размещение трубки для инжекции атомов

        G4LogicalVolume *InTub_log = new G4LogicalVolume(Inlet_tub, AluminMaterial,
                                                         "InTub_log", 0, 0, 0);
#ifdef INJECT_TUBE

        vol_phys = new G4PVPlacement(G4Transform3D(Rotate90Y180X, G4ThreeVector((1.2 + 10.) * cm, 0.0 * cm, 0.0 * cm)),
                                     InTub_log, "InTub_phys", fWorldLV_rot, false, -1, fCheckOverlaps);

#endif

        G4VisAttributes * Alum_VisAtt_ = new G4VisAttributes(*Alum_VisAtt);
      Alum_VisAtt_->SetForceSolid();
        InTub_log->SetVisAttributes(Alum_VisAtt_);
      // InTub_log->SetVisAttributes(Alum_VisAtt);

      // Windows
        G4Box *Window_box = new G4Box("Window_box", 5.2 / 2 * cm, 3.0 / 2 * mm, (w_lng + 1.4 * cm) / 2);
        G4LogicalVolume *Window_log = new G4LogicalVolume(Window_box, VacuumMaterial, "Window_log", 0, 0, 0);
        Window_log->SetVisAttributes(G4VisAttributes::GetInvisible());


        G4Box *W1b = new G4Box("W1b", 2.0 / 2 * mm, 2.0 / 2 * mm, w_lng / 2);
        G4LogicalVolume *W1b_log = new G4LogicalVolume(W1b, SteelMaterial, "W1b_log", 0, 0, 0);
        W1b_log->SetVisAttributes(Steel_VisAtt);

        vol_phys = new G4PVPlacement(0, G4ThreeVector(2.1 * cm, -0.5 * mm, 0.0),
                                     W1b_log, "W1b", Window_log, false, -1, fCheckOverlaps);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(-2.1 * cm, -0.5 * mm, 0.0),
                                     W1b_log, "W1b", Window_log, false, -1, fCheckOverlaps);

        G4Box *W2b = new G4Box("W2b", 4.4 / 2 * cm, 2.0 / 2 * mm, 3.0 / 2 * mm);
        G4LogicalVolume *W2b_log = new G4LogicalVolume(W2b, SteelMaterial, "W2b_log", 0, 0, 0);
        W2b_log->SetVisAttributes(Steel_VisAtt);

        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0, -0.5 * mm, (w_lng + 3. * mm) / 2.),
                                     W2b_log, "W2b", Window_log, false, -1, fCheckOverlaps);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0, -0.5 * mm, -(w_lng + 3. * mm) / 2.),
                                     W2b_log, "W2b", Window_log, false, -1, fCheckOverlaps);

        G4Box *W3b = new G4Box("W3b", 6.0 / 2 * mm, 1.0 / 2 * mm, w_lng / 2);
        G4LogicalVolume *W3b_log = new G4LogicalVolume(W3b, SteelMaterial, "W3b_log", 0, 0, 0);
        W3b_log->SetVisAttributes(Steel_VisAtt);

        vol_phys = new G4PVPlacement(0, G4ThreeVector(2.3 * cm, 1.0 * mm, 0.0),
                                     W3b_log, "W3b", Window_log, false, -1, fCheckOverlaps);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(-2.3 * cm, 1.0 * mm, 0.0),
                                     W3b_log, "W3b", Window_log, false, -1, fCheckOverlaps);

        G4Box *W4b = new G4Box("W4b", 5.2 / 2 * cm, 1.0 / 2 * mm, 7.0 / 2 * mm);
        G4LogicalVolume *W4b_log = new G4LogicalVolume(W4b, SteelMaterial, "W4b_log", 0, 0, 0);
        W4b_log->SetVisAttributes(Steel_VisAtt);

        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0, 1.0 * mm, (w_lng + 7. * mm) / 2.),
                                     W4b_log, "W4b", Window_log, false, -1, fCheckOverlaps);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0, 1.0 * mm, -(w_lng + 7. * mm) / 2.),
                                     W4b_log, "W4b", Window_log, false, -1, fCheckOverlaps);


        vol_phys = new G4PVPlacement(0,
                                     G4ThreeVector(0.0, 2.65 * cm, w_pos), Window_log, "Volume_phys", fWorldLV_rot,
                                     false, -1, fCheckOverlaps);


        vol_phys = new G4PVPlacement(G4Transform3D(Rotate180Z, G4ThreeVector(0.0, -2.65 * cm, w_pos)), Window_log,
                                     "Volume_phys", fWorldLV_rot, false, -1, fCheckOverlaps);


        //   G4UnionSolid *Volume_9 =
//   new G4UnionSolid("Volume_9", Volume_8, Flange,
// 			G4Transform3D(RotateNull, G4ThreeVector(0.0*cm, 0.0*cm, -26.7*cm)));

//   G4UnionSolid *Volume_10 =
//   new G4UnionSolid("Volume_10", Volume_9, Flange,
//  			G4Transform3D(Rotate180X, G4ThreeVector(0.0*cm, 0.0*cm, 26.7*cm)));


#endif

        ////////////////////////////////////////////////////////////////////////////////
        // COUPPER FOILS SHIELD
        ////////////////////////////////////////////////////////////////////////////////

#if defined(SHIELD)

        G4double thick_shield = 30. * um;

        G4Box *CoipperShield_box1 =
                new G4Box("CoipperShield_box1", 6.3 / 2 * cm, thick_shield / 2,
                          (52.8 * cm - (w_lng + 1.4 * cm) - w_pos - 2.0 * cm) / 2. / 2.);

        G4Box *CoipperShield_box2 =
                new G4Box("CoipperShield_box2", 6.3 / 2 * cm, thick_shield / 2,
                          (52.8 * cm - (w_lng + 1.4 * cm) + w_pos + 1.5 * cm) / 2. / 2.);

        G4Box *CoipperShield_box3 =
                new G4Box("CoipperShield_box3", 6.3 / 2 * cm, 6.3 / 2 * cm, thick_shield / 2);




        G4Box *CoipperShield_box4 =
                new G4Box("CoipperShield_box4", thick_shield / 2, 6.3 / 2 * cm, 52.8 / 2 * cm + 0.5 * mm);

        G4Box *CoipperShield_box5 =
                new G4Box("CoipperShield_box5", thick_shield / 2, 6.3 / 2 * cm, 52.8 / 4 * cm - d_y1 / 2 - 1.0 * mm);

//G4cout<< w_pos <<"!!!!" << G4endl;

        G4LogicalVolume *
                CoipperShield_log1 = new G4LogicalVolume(CoipperShield_box1, CopperMaterial, "CoipperShield_log", 0, 0,
                                                         0);

        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0 * cm, 2.81 * cm + 0.2 * mm, w_pos + w_lng / 4 - 4 * cm / 4 +
                                                                                      (52.8 * cm - (w_lng + 1.4 * cm)) /
                                                                                      2.),
                              CoipperShield_log1, "CoipperShield_phys", fWorldLV_rot, false, -1, fCheckOverlaps);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0 * cm, -2.81 * cm - 0.2 * mm, w_pos + w_lng / 4 - 4 * cm / 4 +
                                                                                       (52.8 * cm -
                                                                                        (w_lng + 1.4 * cm)) / 2.),
                              CoipperShield_log1, "CoipperShield_phys", fWorldLV_rot, false, -1, fCheckOverlaps);

        G4LogicalVolume *
                CoipperShield_log2 = new G4LogicalVolume(CoipperShield_box2, CopperMaterial, "CoipperShield_log", 0, 0,
                                                         0);

        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0 * cm, 2.81 * cm + 0.2 * mm, w_pos - w_lng / 4 - 4 * cm / 4 -
                                                                                      (52.8 * cm - (w_lng + 1.4 * cm)) /
                                                                                      2.),
                              CoipperShield_log2, "CoipperShield_phys", fWorldLV_rot, false, -1, fCheckOverlaps);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0 * cm, -2.81 * cm - 0.2 * mm, w_pos - w_lng / 4 - 4 * cm / 4 -
                                                                                       (52.8 * cm -
                                                                                        (w_lng + 1.4 * cm)) / 2.),
                              CoipperShield_log2, "CoipperShield_phys", fWorldLV_rot, false, -1, fCheckOverlaps);

        G4LogicalVolume *
                CoipperShield_log3 = new G4LogicalVolume(CoipperShield_box3, CopperMaterial, "CoipperShield_log", 0, 0,
                                                         0);

        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0 * cm, 0.0 * cm, (52.8 + 0.2) / 2 * cm), CoipperShield_log3,
                                     "CoipperShield_phys", fWorldLV_rot, false, -1, fCheckOverlaps);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0 * cm, 0.0 * cm, -(52.8 + 0.2) / 2 * cm), CoipperShield_log3,
                                     "CoipperShield_phys", fWorldLV_rot, false, -1, fCheckOverlaps);


        G4LogicalVolume *
                CoipperShield_log4 = new G4LogicalVolume(CoipperShield_box4, CopperMaterial, "CoipperShield_log", 0, 0,
                                                         0);

        vol_phys = new G4PVPlacement(0, G4ThreeVector(-3.15 * cm - 0.2 * mm, 0.0 * cm, 0.0 * cm), CoipperShield_log4,
                                     "CoipperShield_phys", fWorldLV_rot, false, -1, fCheckOverlaps);

        G4LogicalVolume *
                CoipperShield_log5 = new G4LogicalVolume(CoipperShield_box5, CopperMaterial, "CoipperShield_log", 0, 0,
                                                         0);

        vol_phys = new G4PVPlacement(0, G4ThreeVector(3.15 * cm + 0.2 * mm, 0.0 * cm, -(52.8 /2  * cm +
                                                                                       1.5* d_y1) / 2.),
                              CoipperShield_log5, "CoipperShield_phys", fWorldLV_rot, false, -1, fCheckOverlaps);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(3.15 * cm + 0.2 * mm, 0.0 * cm, (52.8 / 2 * cm +
                                             1.5* d_y1) / 2.),
                              CoipperShield_log5, "CoipperShield_phys", fWorldLV_rot, false, -1, fCheckOverlaps);

        CoipperShield_log1->SetVisAttributes(Shield_VisAtt);
        CoipperShield_log2->SetVisAttributes(Shield_VisAtt);
        CoipperShield_log3->SetVisAttributes(Shield_VisAtt);
        CoipperShield_log4->SetVisAttributes(Shield_VisAtt);
        CoipperShield_log5->SetVisAttributes(Shield_VisAtt);


#endif //SHIELD

        ////////////////////////////////////////////////////////////////////////////////
        // TITAN FOILS
        ////////////////////////////////////////////////////////////////////////////////

#if defined(FOIL)

        // титановая фольга 70 мкм закрывает вылет из окна
#if defined(RUN21_for_TARGET)

        G4Box *TitanFoil_box =
                new G4Box("TitanFoil_box", 5.2 / 2 * cm, 0.007 / 2 * cm, (w_lng + 1.4 * cm) / 2.);

        G4LogicalVolume *
                TitanFoil_log = new G4LogicalVolume(TitanFoil_box, TitanMaterial, "TitanFoil_log", 0, 0, 0);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0 * cm, 2.81 * cm, w_pos), TitanFoil_log, "TitanFoil_phys",
                                     fWorldLV, false, -1, fCheckOverlaps);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0 * cm, -2.81 * cm, w_pos), TitanFoil_log, "TitanFoil_phys",
                                     fWorldLV, false, -1, fCheckOverlaps);


        TitanFoil_log->SetVisAttributes(TitanFoil_VisAtt);
#endif

        // бериллиевая фольга 300 мкм закрывает вылет из окна

#if defined(RUN23_for_TARGET)

        G4Box *BerylliumFoil_box =
                new G4Box("BerylliumFoil_box", 5.2 / 2 * cm, 300. / 2 * um, (w_lng + 1.4 * cm) / 2.);

        G4LogicalVolume *
                BerylliumFoil_log = new G4LogicalVolume(BerylliumFoil_box, BerylliumMaterial, "BerylliumFoil_log", 0, 0,
                                                        0);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0 * cm, 2.82 * cm, w_pos), BerylliumFoil_log,
                                     "BerylliumFoil_phys", fWorldLV_rot, false, -1, fCheckOverlaps);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0 * cm, -2.82 * cm, w_pos), BerylliumFoil_log,
                                     "BerylliumFoil_phys", fWorldLV_rot, false, -1, fCheckOverlaps);


        BerylliumFoil_log->SetVisAttributes(BerylliumFoil_VisAtt);
#endif

#endif //FOIL


        ////////////////////////////////////////////////////////////////////////////////
        // VOLUME : Storage cell
        ////////////////////////////////////////////////////////////////////////////////

#define CELL_THICK 0.05

#ifdef CELL

        auto CellMaterial = VacuumMaterial;
        CellMaterial = G4NistManager::Instance()->FindOrBuildMaterial("DeuteriumGas");

        G4Tubs *Cell_tube =
                new G4Tubs("Cell", 13. / 2. * mm, (13. + CELL_THICK) / 2. * mm, 40. / 2. * cm, 90. * deg, 180. * deg);
        G4Tubs *Cell_tube_in = new G4Tubs("Cell_in", 0.0, 13. / 2. * mm, 40. / 2. * cm, 90. * deg, 180. * deg);

        G4LogicalVolume *Cell_log = new G4LogicalVolume(Cell_tube, AluminMaterial, "Cell_log", 0, 0, 0);
        Cell_log->SetVisAttributes(Alum_VisAtt);

        G4LogicalVolume *Cell_in_log = new G4LogicalVolume(Cell_tube_in, CellMaterial, "Cell_in_log", 0, 0, 0);

        vol_phys = new G4PVPlacement(0, G4ThreeVector(-(24. - 13.) / 2. * mm, 0., 0.), Cell_log, "Cell", fWorldLV_rot,
                                     false, -1, fCheckOverlaps);
        vol_phys = new G4PVPlacement(&Rotate180Z, G4ThreeVector((24. - 13.) / 2. * mm, 0., 0.), Cell_log, "Cell",
                                     fWorldLV_rot, false, -1, fCheckOverlaps);
        // вакуум внутри ячейки

        vol_phys = new G4PVPlacement(0, G4ThreeVector(-(24. - 13.) / 2. * mm, 0., 0.), Cell_in_log, "Cell_in",
                                     fWorldLV_rot, false, -1, fCheckOverlaps);
        vol_phys = new G4PVPlacement(&Rotate180Z, G4ThreeVector((24. - 13.) / 2. * mm, 0., 0.), Cell_in_log, "Cell_in",
                                     fWorldLV_rot, false, -1, fCheckOverlaps);


        G4Box *Cell_plate = new G4Box("Cplate", (24. - 13.) / 2 * mm, CELL_THICK / 2 * mm, 40. / 2. * cm);
        G4Box *Cell_plate_in = new G4Box("Cplate_in", (24. - 13.) / 2 * mm, (13.) / 2 * mm, 40. / 2. * cm);

        Cell_log = new G4LogicalVolume(Cell_plate, AluminMaterial, "Cell2_log", 0, 0, 0);
        Cell_log->SetVisAttributes(Alum_VisAtt);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0., (13. + CELL_THICK) / 2., 0.), Cell_log, "Cell", fWorldLV_rot,
                                     false, -1, fCheckOverlaps);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0., -(13. + CELL_THICK) / 2., 0.), Cell_log, "Cell", fWorldLV_rot,
                                     false, -1, fCheckOverlaps);

        // вакуум внутри ячейки

        Cell_log = new G4LogicalVolume(Cell_plate_in, CellMaterial, "Cell2_in_log", 0, 0, 0);
        // Cell_log->SetVisAttributes();
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), Cell_log, "Cell_in", fWorldLV_rot, false, -1,
                                     fCheckOverlaps);


#endif // CELL

        // вакуумная камера ускорителя

#ifdef VACUUM_TUBE

        // G4double Xvac_tube_vepp3 = 50. * mm; // 55
        // G4double Yvac_tube_vepp3 = 27 * mm;  // 27
        //  G4double Xvac_tube_vepp3 = 38. * mm; // 55
        //  G4double Yvac_tube_vepp3 = 20 * mm;  // 27
        G4double Xvac_tube_in_vepp3 = 46. * mm; // сечение внутри промежутка
        G4double Yvac_tube_in_vepp3 = 30 * mm;
        G4double Xvac_tube_out_vepp3 = 36. * mm; // сечение на выход/вход в ВЭПП-3
        G4double Yvac_tube_out_vepp3 = 20 * mm;

        G4double Zvac_tube1_vepp3 = 0.019 * m;
        G4double Zvac_tube12_vepp3 = 0.068 * m;
        G4double Zvac_tube13_vepp3 = 0.1 * m;

        G4double Zvac_tube2_vepp3 = 0.378 * m;
        G4double Zvac_tube22_vepp3 = 0.065 * m;
        G4double Zvac_tube23_vepp3 = 0.1 * m;

        G4double deltaZvac_tube_vepp3 = 25.0 * cm;

        G4double VAC_TUBE_THICK = 2.0 * mm;

        G4EllipticalTube *vac_tube1_vepp3_out = new G4EllipticalTube(
            "vac_tube1_vepp3_out", Xvac_tube_in_vepp3 / 2. * mm, Yvac_tube_in_vepp3 / 2. * mm, Zvac_tube1_vepp3 / 2);
        G4EllipticalTube *vac_tube1_vepp3_in = new G4EllipticalTube(
            "vac_tube1_vepp3_in", vac_tube1_vepp3_out->GetDx() - VAC_TUBE_THICK,
            vac_tube1_vepp3_out->GetDy() - VAC_TUBE_THICK, vac_tube1_vepp3_out->GetDz() + 0.2 * mm);
        G4SubtractionSolid *vac_tube1_vepp3 =
            new G4SubtractionSolid("vac_tube1_vepp3", vac_tube1_vepp3_out, vac_tube1_vepp3_in);

        G4EllipticalTube *vac_tube12_vepp3_out =
            new G4EllipticalTube("vac_tube12_vepp3_out", Xvac_tube_out_vepp3 / 2. * mm, Yvac_tube_out_vepp3 / 2. * mm,
                                 Zvac_tube12_vepp3 / 2);
        G4EllipticalTube *vac_tube12_vepp3_in = new G4EllipticalTube(
            "vac_tube12_vepp3_in", vac_tube12_vepp3_out->GetDx() - VAC_TUBE_THICK,
            vac_tube12_vepp3_out->GetDy() - VAC_TUBE_THICK, vac_tube12_vepp3_out->GetDz() + 0.2 * mm);
        G4SubtractionSolid *vac_tube12_vepp3 =
            new G4SubtractionSolid("vac_tube12_vepp3", vac_tube12_vepp3_out, vac_tube12_vepp3_in);


        G4EllipticalTube *vac_tube13_vepp3_out =
            new G4EllipticalTube("vac_tube13_vepp3_out", Xvac_tube_out_vepp3 / 2. * mm, Yvac_tube_out_vepp3 / 2. * mm,
                                 Zvac_tube13_vepp3 / 2);
        G4EllipticalTube *vac_tube13_vepp3_in = new G4EllipticalTube(
            "vac_tube13_vepp3_in", vac_tube13_vepp3_out->GetDx() - VAC_TUBE_THICK,
            vac_tube13_vepp3_out->GetDy() - VAC_TUBE_THICK, vac_tube13_vepp3_out->GetDz() + 0.2 * mm);
        G4SubtractionSolid *vac_tube13_vepp3 =
            new G4SubtractionSolid("vac_tube13_vepp3", vac_tube13_vepp3_out, vac_tube13_vepp3_in);

        /////

        G4EllipticalTube *vac_tube2_vepp3_out = new G4EllipticalTube(
            "vac_tube2_vepp3_out", Xvac_tube_out_vepp3 / 2. * mm, Yvac_tube_out_vepp3 / 2. * mm, Zvac_tube2_vepp3 / 2);
        G4EllipticalTube *vac_tube2_vepp3_in = new G4EllipticalTube(
            "vac_tube2_vepp3_in", vac_tube2_vepp3_out->GetDx() - VAC_TUBE_THICK,
            vac_tube2_vepp3_out->GetDy() - VAC_TUBE_THICK, vac_tube2_vepp3_out->GetDz() + 0.2 * mm);
        G4SubtractionSolid *vac_tube2_vepp3 =
            new G4SubtractionSolid("vac_tube2_vepp3", vac_tube2_vepp3_out, vac_tube2_vepp3_in);

        G4EllipticalTube *vac_tube22_vepp3_out =
            new G4EllipticalTube("vac_tube22_vepp3_out", Xvac_tube_out_vepp3 / 2. * mm, Yvac_tube_out_vepp3 / 2. * mm,
                                 Zvac_tube22_vepp3 / 2);
        G4EllipticalTube *vac_tube22_vepp3_in = new G4EllipticalTube(
            "vac_tube22_vepp3_in", vac_tube22_vepp3_out->GetDx() - VAC_TUBE_THICK,
            vac_tube22_vepp3_out->GetDy() - VAC_TUBE_THICK, vac_tube22_vepp3_out->GetDz() + 0.2 * mm);
        G4SubtractionSolid *vac_tube22_vepp3 =
            new G4SubtractionSolid("vac_tube22_vepp3", vac_tube22_vepp3_out, vac_tube22_vepp3_in);


        G4EllipticalTube *vac_tube23_vepp3_out =
            new G4EllipticalTube("vac_tube23_vepp3_out", Xvac_tube_out_vepp3 / 2. * mm, Yvac_tube_out_vepp3 / 2. * mm,
                                 Zvac_tube23_vepp3 / 2);
        G4EllipticalTube *vac_tube23_vepp3_in = new G4EllipticalTube(
            "vac_tube23_vepp3_in", vac_tube23_vepp3_out->GetDx() - VAC_TUBE_THICK,
            vac_tube23_vepp3_out->GetDy() - VAC_TUBE_THICK, vac_tube23_vepp3_out->GetDz() + 0.2 * mm);
        G4SubtractionSolid *vac_tube23_vepp3 =
            new G4SubtractionSolid("vac_tube23_vepp3", vac_tube23_vepp3_out, vac_tube23_vepp3_in);


        G4LogicalVolume *vac_tube1_in_vepp3_log =
            new G4LogicalVolume(vac_tube1_vepp3_in, VacuumMaterial, "vac_tube1_vepp3_in_log", 0, 0, 0);
        G4LogicalVolume *vac_tube2_in_vepp3_log =
            new G4LogicalVolume(vac_tube2_vepp3_in, VacuumMaterial, "vac_tube2_vepp3_in_log", 0, 0, 0);

        G4LogicalVolume *vac_tube12_in_vepp3_log =
            new G4LogicalVolume(vac_tube12_vepp3_in, VacuumMaterial, "vac_tube12_vepp3_in_log", 0, 0, 0);
        G4LogicalVolume *vac_tube13_in_vepp3_log =
            new G4LogicalVolume(vac_tube13_vepp3_in, VacuumMaterial, "vac_tube13_vepp3_in_log", 0, 0, 0);

        G4LogicalVolume *vac_tube22_in_vepp3_log =
            new G4LogicalVolume(vac_tube22_vepp3_in, VacuumMaterial, "vac_tube22_vepp3_in_log", 0, 0, 0);
        G4LogicalVolume *vac_tube23_in_vepp3_log =
            new G4LogicalVolume(vac_tube23_vepp3_in, VacuumMaterial, "vac_tube23_vepp3_in_log", 0, 0, 0);


        G4LogicalVolume *vac_tube1_out_vepp3_log =
            new G4LogicalVolume(vac_tube1_vepp3_out, VacuumMaterial, "vac_tube1_vepp3_out_log", 0, 0, 0);
        G4LogicalVolume *vac_tube2_out_vepp3_log =
            new G4LogicalVolume(vac_tube2_vepp3_out, VacuumMaterial, "vac_tube2_vepp3_out_log", 0, 0, 0);

        G4LogicalVolume *vac_tube12_out_vepp3_log =
            new G4LogicalVolume(vac_tube12_vepp3_out, VacuumMaterial, "vac_tube12_vepp3_out_log", 0, 0, 0);
        G4LogicalVolume *vac_tube13_out_vepp3_log =
            new G4LogicalVolume(vac_tube13_vepp3_out, VacuumMaterial, "vac_tube13_vepp3_out_log", 0, 0, 0);

        G4LogicalVolume *vac_tube22_out_vepp3_log =
            new G4LogicalVolume(vac_tube22_vepp3_out, VacuumMaterial, "vac_tube22_vepp3_out_log", 0, 0, 0);
        G4LogicalVolume *vac_tube23_out_vepp3_log =
            new G4LogicalVolume(vac_tube23_vepp3_out, VacuumMaterial, "vac_tube23_vepp3_out_log", 0, 0, 0);


        G4LogicalVolume *vac_tube1_vepp3_log =
            new G4LogicalVolume(vac_tube1_vepp3, AluminMaterial, "vac_tube1_vepp3_log", 0, 0, 0);
        G4LogicalVolume *vac_tube2_vepp3_log =
            new G4LogicalVolume(vac_tube2_vepp3, AluminMaterial, "vac_tube2_vepp3_log", 0, 0, 0);

        G4LogicalVolume *vac_tube12_vepp3_log =
            new G4LogicalVolume(vac_tube12_vepp3, AluminMaterial, "vac_tube12_vepp3_log", 0, 0, 0);
        G4LogicalVolume *vac_tube13_vepp3_log =
            new G4LogicalVolume(vac_tube13_vepp3, AluminMaterial, "vac_tube13_vepp3_log", 0, 0, 0);

        G4LogicalVolume *vac_tube22_vepp3_log =
            new G4LogicalVolume(vac_tube22_vepp3, AluminMaterial, "vac_tube22_vepp3_log", 0, 0, 0);
        G4LogicalVolume *vac_tube23_vepp3_log =
            new G4LogicalVolume(vac_tube23_vepp3, AluminMaterial, "vac_tube23_vepp3_log", 0, 0, 0);


        auto rot1_for_world_target = new G4RotationMatrix(Rotate180Z);
        rot1_for_world_target->rotateY(-alpha0_GEM);
        auto rot2_for_world_target = new G4RotationMatrix(Rotate180Z);
        rot2_for_world_target->rotateY(theta0_GEM + 0.0 * alpha0_GEM);

        //  new G4PVPlacement(rot1_for_world_target, G4ThreeVector(2.2 * cm, 0., -(Zvac_tube1_vepp3 /2  +
        //  deltaZvac_tube_vepp3)),vac_tube1_in_vepp3_log, "vac_tube1_vepp3_in", fWorldLV, false, -1, fCheckOverlaps);
        //  new G4PVPlacement(rot2_for_world_target, G4ThreeVector(0.15 * cm, 0., Zvac_tube2_vepp3 / 2 +
        //  deltaZvac_tube_vepp3 + 41.0 * cm),vac_tube2_in_vepp3_log, "vac_tube2_vepp3_in", fWorldLV, false, -1,
        //  fCheckOverlaps);

        //  new G4PVPlacement(rot1_for_world_target, G4ThreeVector(1.3 * cm, 0., -(Zvac_tube1_vepp3 /2  +
        //  deltaZvac_tube_vepp3)),vac_tube1_out_vepp3_log, "vac_tube1_vepp3_out", fWorldLV, false, -1, fCheckOverlaps);
        //  new G4PVPlacement(rot2_for_world_target, G4ThreeVector(0.8 * cm, 0., Zvac_tube2_vepp3 / 2 +
        //  deltaZvac_tube_vepp3 + 41.0 * cm + 8.2 * cm),vac_tube2_out_vepp3_log, "vac_tube2_vepp3_out", fWorldLV,
        //  false, -1, fCheckOverlaps);

        new G4PVPlacement(rot1_for_world_target,
                          G4ThreeVector(0.95 * cm, 0., -(Zvac_tube1_vepp3 / 2 + deltaZvac_tube_vepp3 + 6.0 * cm)),
                          vac_tube1_vepp3_log, "vac_tube1_vepp3", fWorldLV, false, -1, fCheckOverlaps);
        new G4PVPlacement(
            rot1_for_world_target,
            G4ThreeVector(1.3 * cm + 0.8 * cm, 0.,
                          -(Zvac_tube12_vepp3 / 2 + Zvac_tube1_vepp3 + deltaZvac_tube_vepp3 + 6.0 * cm + 16.2 * cm)),
            vac_tube12_out_vepp3_log, "vac_tube12_vepp3_out", fWorldLV, false, -1, fCheckOverlaps);
        new G4PVPlacement(0,
                          G4ThreeVector(1.3 * cm + 1.2 * cm, 0.,
                                        -(Zvac_tube13_vepp3 / 2 + Zvac_tube1_vepp3 + Zvac_tube12_vepp3 +
                                          deltaZvac_tube_vepp3 + 6.0 * cm + 16.2 + 26.8 * cm)),
                          vac_tube13_out_vepp3_log, "vac_tube13_vepp3_out", fWorldLV, false, -1, fCheckOverlaps);


        new G4PVPlacement(rot2_for_world_target,
                          G4ThreeVector(0.06 * cm, 0., Zvac_tube2_vepp3 / 2 + deltaZvac_tube_vepp3 + 41.0 * cm),
                          vac_tube2_vepp3_log, "vac_tube2_vepp3", fWorldLV, false, -1, fCheckOverlaps);
        new G4PVPlacement(rot2_for_world_target,
                          G4ThreeVector(0.8 * cm + 0.6 * cm, 0.,
                                        Zvac_tube22_vepp3 / 2 + Zvac_tube2_vepp3 + deltaZvac_tube_vepp3 + 14.7 * cm +
                                            41.0 * cm + 1.6 * cm),
                          vac_tube22_out_vepp3_log, "vac_tube22_vepp3_out", fWorldLV, false, -1, fCheckOverlaps);
        new G4PVPlacement(0,
                          G4ThreeVector(0.8 * cm + 1.9 * cm, 0.,
                                        Zvac_tube23_vepp3 / 2 + Zvac_tube2_vepp3 + Zvac_tube22_vepp3 +
                                            deltaZvac_tube_vepp3 + 14.7 * cm + 41.0 * cm + 1.6 * cm + 14.8 * cm),
                          vac_tube23_out_vepp3_log, "vac_tube23_vepp3_out", fWorldLV, false, -1, fCheckOverlaps);

        ///

        new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), vac_tube1_vepp3_log, "vac_tube1_vepp3", vac_tube1_out_vepp3_log,
                          false, -1, fCheckOverlaps);
        new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), vac_tube2_vepp3_log, "vac_tube2_vepp3", vac_tube2_out_vepp3_log,
                          false, -1, fCheckOverlaps);

        new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), vac_tube12_vepp3_log, "vac_tube12_vepp3",
                          vac_tube12_out_vepp3_log, false, -1, fCheckOverlaps);
        new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), vac_tube13_vepp3_log, "vac_tube13_vepp3",
                          vac_tube13_out_vepp3_log, false, -1, fCheckOverlaps);

        new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), vac_tube22_vepp3_log, "vac_tube22_vepp3",
                          vac_tube22_out_vepp3_log, false, -1, fCheckOverlaps);
        new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), vac_tube23_vepp3_log, "vac_tube23_vepp3",
                          vac_tube23_out_vepp3_log, false, -1, fCheckOverlaps);


        auto VacTubeVepp3_VissAtt = new G4VisAttributes(true);
        VacTubeVepp3_VissAtt->SetColour(Alum_VisAtt->GetColor().GetRed(), Alum_VisAtt->GetColor().GetGreen(),
                                        Alum_VisAtt->GetColor().GetBlue(), 0.4);
        VacTubeVepp3_VissAtt->SetForceSolid(true);
        // VacTubeVepp3_VissAtt->SetForceWireframe(true);


        vac_tube1_vepp3_log->SetVisAttributes(VacTubeVepp3_VissAtt);
        vac_tube2_vepp3_log->SetVisAttributes(VacTubeVepp3_VissAtt);

        vac_tube12_vepp3_log->SetVisAttributes(VacTubeVepp3_VissAtt);
        vac_tube13_vepp3_log->SetVisAttributes(VacTubeVepp3_VissAtt);

        vac_tube22_vepp3_log->SetVisAttributes(VacTubeVepp3_VissAtt);
        vac_tube23_vepp3_log->SetVisAttributes(VacTubeVepp3_VissAtt);

        auto VacTubeInVepp3_VissAtt = new G4VisAttributes(true);
        // VacTubeInVepp3_VissAtt->SetColour(Gas_VisAtt->GetColor().GetRed(),Gas_VisAtt->GetColor().GetGreen(),Gas_VisAtt->GetColor().GetBlue(),
        // 0.2);
        VacTubeInVepp3_VissAtt->SetColour(0., 0., 1.0, 0.2);
        VacTubeInVepp3_VissAtt->SetForceSolid(true);
        // VacTubeInVepp3_VissAtt->SetForceWireframe(true);


        vac_tube1_in_vepp3_log->SetVisAttributes(VacTubeInVepp3_VissAtt);
        vac_tube2_in_vepp3_log->SetVisAttributes(VacTubeInVepp3_VissAtt);

        vac_tube12_in_vepp3_log->SetVisAttributes(VacTubeInVepp3_VissAtt);
        vac_tube13_in_vepp3_log->SetVisAttributes(VacTubeInVepp3_VissAtt);

        vac_tube22_in_vepp3_log->SetVisAttributes(VacTubeInVepp3_VissAtt);
        vac_tube23_in_vepp3_log->SetVisAttributes(VacTubeInVepp3_VissAtt);

        vac_tube1_out_vepp3_log->SetVisAttributes(VacTubeInVepp3_VissAtt);
        vac_tube2_out_vepp3_log->SetVisAttributes(VacTubeInVepp3_VissAtt);

        vac_tube12_out_vepp3_log->SetVisAttributes(VacTubeInVepp3_VissAtt);
        vac_tube13_out_vepp3_log->SetVisAttributes(VacTubeInVepp3_VissAtt);

        vac_tube22_out_vepp3_log->SetVisAttributes(VacTubeInVepp3_VissAtt);
        vac_tube23_out_vepp3_log->SetVisAttributes(VacTubeInVepp3_VissAtt);

        //////////

#endif

#endif // TARGET

    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructMagnetCell()
    {

        // МАГНИТЫ

#ifdef TARGET
        //////////////////////////////////////////////////////////////////////
        // 		Magnet
        //////////////////////////////////////////////////////////////////////
#ifdef MAGNET_CELL

        auto CopperMaterial = G4Material::GetMaterial("Copper");

        G4double z_size_mag = 50. * cm;
        G4double y_size_mag = 10. * cm;
        G4double x_size_mag = 6. * cm;

        G4double d_z_mag = 2. * cm;
        G4double x_pos_mag = 8.5 * cm;    // was 10.5 cm
        G4double d_y1_mag = 1. * cm;
        G4double y_pos_mag = (d_y1_mag + y_size_mag / 2.) / 2.;
        G4double z_pos_mag = 0. * cm;    // was 5 cm


        G4Box *solidMag1 = new G4Box("Mag1", x_size_mag / 2., y_size_mag / 2., z_size_mag / 2.);
        G4LogicalVolume *logicMag1 = new G4LogicalVolume(solidMag1, CopperMaterial, "Mag1");
        new G4PVPlacement(0, G4ThreeVector(x_pos_mag, 0., z_pos_mag), logicMag1, "Mag1", fWorldLV_rot, false, -1,
                          fCheckOverlaps);
        new G4PVPlacement(0, G4ThreeVector(-x_pos_mag, 0., z_pos_mag), logicMag1, "Mag1", fWorldLV_rot, false, -1,
                          fCheckOverlaps);

        G4Box *solidMag1a = new G4Box("Mag1a", x_size_mag / 2. + d_z_mag, (y_size_mag / 2. - d_y1_mag) / 2.,
                                      z_size_mag / 2. + d_z_mag);
        G4LogicalVolume *logicMag1a = new G4LogicalVolume(solidMag1a, CopperMaterial, "Mag1a");

        new G4PVPlacement(0, G4ThreeVector(x_pos_mag, y_pos_mag, z_pos_mag), logicMag1a, "Mag1a", fWorldLV_rot, false,
                          -1, fCheckOverlaps);
        new G4PVPlacement(0, G4ThreeVector(x_pos_mag, -y_pos_mag, z_pos_mag), logicMag1a, "Mag1a", fWorldLV_rot, false,
                          -1, fCheckOverlaps);
        new G4PVPlacement(0, G4ThreeVector(-x_pos_mag, y_pos_mag, z_pos_mag), logicMag1a, "Mag1a", fWorldLV_rot, false,
                          -1, fCheckOverlaps);
        new G4PVPlacement(0, G4ThreeVector(-x_pos_mag, -y_pos_mag, z_pos_mag), logicMag1a, "Mag1a", fWorldLV_rot, false,
                          -1, fCheckOverlaps);

        G4Tubs *solidMag2 = new G4Tubs("Mag2", 0., x_size_mag / 2., z_size_mag / 2., 0., M_PI);
        G4LogicalVolume *logicMag2 = new G4LogicalVolume(solidMag2, CopperMaterial, "Mag2");
        new G4PVPlacement(0, G4ThreeVector(x_pos_mag, y_size_mag / 2., z_pos_mag), logicMag2, "Mag2", fWorldLV_rot,
                          false, -1, fCheckOverlaps);
        new G4PVPlacement(0, G4ThreeVector(-x_pos_mag, y_size_mag / 2., z_pos_mag), logicMag2, "Mag2", fWorldLV_rot,
                          false, -1, fCheckOverlaps);
        new G4PVPlacement(G4Transform3D(Rotate180Z, G4ThreeVector(x_pos_mag, -y_size_mag / 2., z_pos_mag)), logicMag2,
                          "Mag2", fWorldLV_rot, false, -1, fCheckOverlaps);
        new G4PVPlacement(G4Transform3D(Rotate180Z, G4ThreeVector(-x_pos_mag, -y_size_mag / 2., z_pos_mag)), logicMag2,
                          "Mag2", fWorldLV_rot, false, -1, fCheckOverlaps);


        logicMag1->SetVisAttributes(Steel_VisAtt);
        logicMag1a->SetVisAttributes(Mag_VisAtt);
        logicMag2->SetVisAttributes(Steel_VisAtt);


#endif // MAGNET_CELL

#endif // TARGET
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void DetectorConstruction::ConstructMagnetD2()
    {

#ifdef MAGNET_D2
        auto Hall_log = fWorldLV;

        G4Colour stst_col(0.53, 0.53, 0.84, 0.5);
        G4VisAttributes *Volume_logVisAtt = new G4VisAttributes(stst_col);

        ////

        G4LogicalVolume *Volume_log = nullptr;
        G4PVPlacement *Volume_phys = nullptr;
        G4Box *b_box = nullptr;
        G4Box *s_box = nullptr;
        G4LogicalVolume *VolVac = nullptr;
        G4SubtractionSolid *Vols1 = nullptr;
        G4SubtractionSolid *Vols2 = nullptr;
        G4EllipticalTube *eltub = nullptr;
        G4EllipticalTube *eltub1 = nullptr;
        G4double x_pos = 0.;
        G4double y_pos = 0.;
        G4double z_pos = 0.;


        /////////////

        /// Upstream pumping volume

        b_box = new G4Box("b_box", 32. * cm / 2., 25. * cm / 2, 6. * cm / 2.);
        VolVac = new G4LogicalVolume(b_box, Vacuum, "Vacuum_log", 0, 0, 0);
        VolVac->SetVisAttributes(G4VisAttributes::GetInvisible());
        s_box = new G4Box("s_box", (32. * cm - 3.05 * mm) / 2., (25. * cm - 6. * mm) / 2, (6. * cm - 6. * mm) / 2.);
        Vols2 =
            new G4SubtractionSolid("Vols2", b_box, s_box, G4Transform3D(RotateNull, G4ThreeVector(3. * mm, 0., 0.)));
        eltub = new G4EllipticalTube("el1", 46. * mm / 2., 30. * mm / 2., 3.02 * mm / 2);

        auto rot_for_world_target = new G4RotationMatrix();
        // rot_for_world_target->rotateY(alpha0_GEM);
        //  rot_for_world_target->rotateY(-alpha0_GEM);
        rot_for_world_target = &RotateNull;

        Vols1 = new G4SubtractionSolid(
            "Vols1", Vols2, eltub,
            G4Transform3D(*rot_for_world_target, G4ThreeVector(1.45 * cm / 2, 0., 6. * cm / 2. - 3. * mm / 2.)));
        Vols2 = new G4SubtractionSolid(
            "Vols2", Vols1, eltub,
            G4Transform3D(*rot_for_world_target, G4ThreeVector(1.85 * cm / 2, 0., -(6. * cm / 2. - 3. * mm / 2.))));

        Volume_log = new G4LogicalVolume(Vols2, Steel, "Volume_log", 0, 0, 0);
        new G4PVPlacement(0, G4ThreeVector(), Volume_log, "dnstream_vac", VolVac, false, 0, fCheckOverlaps);
        z_pos = -(50. * cm / 2. + b_box->GetZHalfLength());
#ifdef UPPUMP
        Volume_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., z_pos), VolVac, "upstream_vac1", Hall_log, false, 0,
                                        fCheckOverlaps);
#endif
        Volume_log->SetVisAttributes(Volume_logVisAtt);
        z_pos -= b_box->GetZHalfLength();

        // entry elliptical vacuum tube
        eltub = new G4EllipticalTube("el0o", 46. * mm / 2., 30. * mm / 2., 160 * mm / 2);
        eltub1 =
            new G4EllipticalTube("el0i", eltub->GetDx() - 2. * mm, eltub->GetDy() - 2. * mm, eltub->GetDz() + 0.2 * mm);
        VolVac = new G4LogicalVolume(eltub, Vacuum, "Vacuum_log", 0, 0, 0);
        VolVac->SetVisAttributes(G4VisAttributes::GetInvisible());
        Vols1 = new G4SubtractionSolid("Vols1", eltub, eltub1);
        Volume_log = new G4LogicalVolume(Vols1, Steel, "Volume_log", 0, 0, 0);
        // #ifdef UPPUMP
        auto rot_for_world_target_ = new G4RotationMatrix();
        rot_for_world_target_->rotateY(alpha0_GEM);
        // rot_for_world_target_ = &RotateNull;
        new G4PVPlacement(G4Transform3D(RotateNull, G4ThreeVector(0., 0., 0.)), Volume_log, "dddstream_vac", VolVac,
                          false, 0, fCheckOverlaps);
        // #endif
        z_pos -= eltub->GetDz();

#ifdef UPPUMP
        // я задаю трубку в другом месте
        //  Volume_phys = new G4PVPlacement(rot_for_world_target_, G4ThreeVector(2.5 * cm / 2, 0., z_pos), VolVac,
        //  "dndnstream_vac", Hall_log, false, 0, fCheckOverlaps);
#endif
        Volume_log->SetVisAttributes(Volume_logVisAtt);

        b_box = new G4Box("b_box", 30. * cm / 2., 25. * cm / 2, 20. * cm / 2.);
        s_box = new G4Box("s_box", (30. * cm - 6. * mm) / 2., (25. * cm - 6. * mm) / 2, (20. * cm - 6. * mm) / 2.);
        Vols2 = new G4SubtractionSolid("Vols2", b_box, s_box);
        Volume_log = new G4LogicalVolume(Vols2, Steel, "Volume_log", 0, 0, 0);
#ifdef UPPUMP
        Volume_phys =
            new G4PVPlacement(0, G4ThreeVector(32. * cm / 2. + 30. * cm / 2., 0., -(50. * cm / 2 + 20. * cm / 2.)),
                              Volume_log, "upstream_vac2", Hall_log, false, 0, fCheckOverlaps);
#endif
        Volume_log->SetVisAttributes(Volume_logVisAtt);


        ////

        G4Trd *trd1 = nullptr;
        G4Trd *trd2 = nullptr;

        //////////////

        /// Downstream pumping volume
        G4double r_ang = 25. * deg;
        trd1 = new G4Trd("trd1", 8. * cm, 8. * cm, (11. - 7. * tan(r_ang)) * cm, 11. * cm, 7. * cm / 2);
        VolVac = new G4LogicalVolume(trd1, Vacuum, "Vacuum_log", 0, 0, 0);
        VolVac->SetVisAttributes(G4VisAttributes::GetInvisible());

        trd2 = new G4Trd("trd2", 8. * cm - 3. * mm / 2., 8. * cm - 3. * mm / 2.,
                         (11. - 7. * tan(r_ang)) * cm - 3. * mm / tan(45. * deg + r_ang / 2.),
                         11. * cm - 3. * mm / tan(45. * deg - r_ang / 2.), 7. * cm / 2 - 3. * mm);
        Vols1 = new G4SubtractionSolid("Vols1", trd1, trd2,
                                       G4Transform3D(RotateNull, G4ThreeVector(3.02 * mm / 2., 0., 0.)));
        eltub = new G4EllipticalTube("el1", 46. * mm / 2., 30. * mm / 2., 3.01 * mm / 2);
        Vols2 = new G4SubtractionSolid(
            "Vols2", Vols1, eltub,
            G4Transform3D(RotateNull, G4ThreeVector(-0.8 * cm, 0., -7. * cm / 2. + 3. * mm / 2.)));
        b_box = new G4Box("bb2", 36. * mm / 2., 20. * mm / 2., 3.01 * mm / 2);

        Vols1 = new G4SubtractionSolid(
            "Vols1", Vols2, b_box, G4Transform3D(RotateNull, G4ThreeVector(-5. * mm, 0., 7. * cm / 2. - 3. * mm / 2.)));

        Volume_log = new G4LogicalVolume(Vols1, Steel, "Volume_log", 0, 0, 0);
#ifdef DOWNPUMP
        new G4PVPlacement(0, G4ThreeVector(), Volume_log, "dnstream_vac", VolVac, false, 0);
#endif
        z_pos = 50. * cm / 2. + 7. * cm / 2.;
        Volume_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., z_pos), VolVac, "dnstream_vac", Hall_log, false, 0);
        Volume_log->SetVisAttributes(Volume_logVisAtt);
        z_pos += 7. * cm / 2.;

        b_box = new G4Box("b_box", 18. * cm / 2., 22. * cm / 2, 7. * cm / 2.);
        VolVac = new G4LogicalVolume(b_box, Vacuum, "Vacuum_log", 0, 0, 0);
        VolVac->SetVisAttributes(G4VisAttributes::GetInvisible());
        s_box = new G4Box("s_box", (18. * cm - 3. * mm) / 2., 22. * cm / 2 - 3. * mm, 7. * cm / 2. - 3. * mm);
        Vols2 = new G4SubtractionSolid("Vols2", b_box, s_box,
                                       G4Transform3D(RotateNull, G4ThreeVector(3.02 * mm / 2., 0., 0.)));
        trd2 = new G4Trd("trd3", 3.05 * mm / 2., 3.05 * mm / 2.,
                         (11. - 7. * tan(r_ang)) * cm - 3. * mm / tan(45. * deg + r_ang / 2.),
                         11. * cm - 3. * mm / tan(45. * deg - r_ang / 2.), 7. * cm / 2 - 3. * mm);
        Vols1 = new G4SubtractionSolid("Vols1", Vols2, trd2,
                                       G4Transform3D(RotateNull, G4ThreeVector(-(18. * cm / 2 - 3. * mm / 2), 0., 0.)));

        Volume_log = new G4LogicalVolume(Vols1, Steel, "Volume_log", 0, 0, 0);
#ifdef DOWNPUMP
        new G4PVPlacement(0, G4ThreeVector(), Volume_log, "dnstream_vac", VolVac, false, 0);
#endif
        Volume_phys = new G4PVPlacement(0, G4ThreeVector(8. * cm + 18. * cm / 2, 0., 50. * cm / 2. + 7. * cm / 2.),
                                        VolVac, "dnstream_vac1", Hall_log, false, 0);
        Volume_log->SetVisAttributes(Volume_logVisAtt);

        b_box = new G4Box("b_box", 12. * cm / 2., 22. * cm / 2, 40. * cm / 2.);
        s_box = new G4Box("s_box", (12. * cm - 6. * mm) / 2., (22. * cm - 6. * mm) / 2, (40. * cm - 6. * mm) / 2.);
        Vols2 = new G4SubtractionSolid("Vols2", b_box, s_box);
        Volume_log = new G4LogicalVolume(Vols2, Steel, "Volume_log", 0, 0, 0);
#ifdef DOWNPUMP
        Volume_phys =
            new G4PVPlacement(0, G4ThreeVector(8. * cm + 18. * cm + 12. * cm / 2., 0., 50. * cm / 2. + 40. * cm / 2.),
                              Volume_log, "dnstream_vac2", Hall_log, false, 0);
#endif
        Volume_log->SetVisAttributes(Volume_logVisAtt);

        // flange for vacuum gauge
        auto Vin_tube = new G4Tubs("Vfl1_box", 50. * mm / 2., 55. * mm / 2., 30. * mm / 2., 0, 360. * deg);
        Volume_log = new G4LogicalVolume(Vin_tube, Steel, "Vfl1_log", 0, 0, 0);
        Volume_log->SetVisAttributes(Volume_logVisAtt);
#ifdef DOWNPUMP
        Volume_phys = new G4PVPlacement(
            G4Transform3D(Rotate90X,
                          G4ThreeVector(8. * cm + 18. * cm + 12. * cm / 2., 22. * cm / 2 + 30. * mm / 2.,
                                        50. * cm / 2. + 40. * cm - 125. * mm)),
            Volume_log, "Vfl1", Hall_log, false, 0);
#endif
        Vin_tube = new G4Tubs("Vfl2_box", 0., 90. * mm / 2., 25. * mm / 2., 0, 360. * deg);
        Volume_log = new G4LogicalVolume(Vin_tube, Steel, "Vfl2_log", 0, 0, 0);
        Volume_log->SetVisAttributes(Volume_logVisAtt);
#ifdef DOWNPUMP
        Volume_phys = new G4PVPlacement(
            G4Transform3D(Rotate90X,
                          G4ThreeVector(8. * cm + 18. * cm + 12. * cm / 2., 22. * cm / 2 + 30. * mm + 25. * mm / 2.,
                                        50. * cm / 2. + 40. * cm - 125. * mm)),
            Volume_log, "Vfl1", Hall_log, false, 0);
#endif

        /////////////

        /// D2 MAGNET
        G4double d2_zpos = 50. * cm,
                 d2_dz = 32. *
            cm; // реальный  размер магнита 302 мм, 320 мм это эффективный размер с учетом вываливания поля
        G4double dx = 49. * cm;
        G4double dy = 30. * cm;
        G4double dz = d2_dz;
        G4double dx1 = 34. * cm, dy1 = 22. * cm;

        b_box = new G4Box("baba", 6. * cm / 2, 4 * cm / 2., (d2_zpos - dz / 2. - z_pos) / 2.);
        VolVac = new G4LogicalVolume(b_box, Vacuum, "Vacuum_log", 0, 0, 0);
        VolVac->SetVisAttributes(G4VisAttributes::GetInvisible());
        s_box = new G4Box("bsbs", b_box->GetXHalfLength() - 3. * mm, b_box->GetYHalfLength() - 3. * mm,
                          b_box->GetZHalfLength() + 0.2 * mm);
        Vols1 = new G4SubtractionSolid("Vols1", b_box, s_box);
        Volume_log = new G4LogicalVolume(Vols1, Steel, "Volume_log", 0, 0, 0);

        G4Colour D2_stst_col(0.53, 0.53, 0.84, 0.6);
        G4VisAttributes *D2_Volume_logVisAtt = new G4VisAttributes(D2_stst_col);
        // D2_Volume_logVisAtt->SetForceSolid(true);
        Volume_log->SetVisAttributes(D2_Volume_logVisAtt);
        new G4PVPlacement(0, G4ThreeVector(), Volume_log, "d_vac", VolVac, false, 0);
        z_pos += b_box->GetZHalfLength();
        Volume_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., z_pos), VolVac, "dnstream_vac1", Hall_log, false, 0);


        //  ubox=new G4Box("D2box",dx/2.,dy/2.,dz/2.);
        trd1 = new G4Trd("D2box", dx1 / 2., dx / 2., dy1 / 2., dy / 2., (dz - 1.8 * cm) / 2.);
        auto trd1_ = new G4Trd("D2box", dx1 / 2., dx / 2., dy1 / 2., dy / 2., (dz - 0.0 * cm) / 2.);

        VolVac = new G4LogicalVolume(trd1, Vacuum, "Vacuum_log", 0, 0, 0);
        auto VolVac_ = new G4LogicalVolume(trd1_, Vacuum, "Vacuum_log", 0, 0, 0);
        // VolVac->SetVisAttributes(G4VisAttributes::GetInvisible());
        VolVac_->SetVisAttributes(G4VisAttributes::GetInvisible());
        dx = 18. * cm;
        dy = 5. * cm;
        dx1 = 10. * cm;
        trd2 = new G4Trd("D2box", dx1 / 2., dx / 2., dy / 2., dy / 2., (dz - 1.8 * cm) / 2. + 0.2 * mm);
        Vols1 = new G4SubtractionSolid("Vols1", trd1, trd2);

        G4LogicalVolume *fLogicDipole = new G4LogicalVolume(Vols1, Iron, "DipoleD2_log", 0, 0, 0);
        //  fLogicDipole->SetVisAttributes(new G4VisAttributes(G4Colour(0.5, 0.7, 1.0)));
        auto D2dipole_col = G4Colour(0.5, 0.7, 1.0, 0.4);
        G4VisAttributes *D2dipole_Volume_logVisAtt = new G4VisAttributes(D2dipole_col);
        D2dipole_Volume_logVisAtt->SetForceSolid(true);
        fLogicDipole->SetVisAttributes(D2dipole_Volume_logVisAtt);
        new G4PVPlacement(0, G4ThreeVector(), fLogicDipole, "MagD2", VolVac, false, -1);

        fLogicD2gap = new G4LogicalVolume(trd2, Vacuum, "D2gap_log", 0, 0, 0);
        new G4PVPlacement(0, G4ThreeVector(), fLogicD2gap, "D2gap", VolVac, false, 0);
        // fLogicD2gap->SetVisAttributes(G4VisAttributes::GetInvisible());

        // дополнительный объем после магнита для вываливания поля

        auto box1 = new G4Box("D2box_add1", dx1 / 2., dy / 2., 18.0 * mm / 4 - 0.001 * mm);
        auto box2 = new G4Box("D2box_add2", dx / 2., dy / 2., 18.0 * mm / 4 - 0.001 * mm);

        fLogicD2gap_add1 = new G4LogicalVolume(box1, Vacuum, "D2gap_add1_log", 0, 0, 0);
        fLogicD2gap_add2 = new G4LogicalVolume(box2, Vacuum, "D2gap_add2_log", 0, 0, 0);

        // fLogicD2gap_add1->SetVisAttributes(new G4VisAttributes(G4Colour(1.0, 1.0, 1.0)));
        //   fLogicD2gap_add2->SetVisAttributes(new G4VisAttributes(G4Colour(1.5, 0.7, 1.0)));

        new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), VolVac, "VolVac_VolVac_", VolVac_, false, 0);


        new G4PVPlacement(0, G4ThreeVector(0., 0., -(d2_dz - 18.0 * mm / 2) / 2), fLogicD2gap_add1, "D2gap_add1",
                          VolVac_, false, 0);
        new G4PVPlacement(0, G4ThreeVector(0., 0., (d2_dz - 18.0 * mm / 2) / 2), fLogicD2gap_add2, "D2gap_add2",
                          VolVac_, false, 0);


        x_pos = 0.0;
        y_pos = 0. * cm;
        z_pos = d2_zpos;

        Volume_phys = new G4PVPlacement(0, G4ThreeVector(x_pos, y_pos, z_pos), VolVac_, "MagD2", Hall_log, false, -1);

        //

        // магнитное поле в SD


        /// Pumping volume behind D2 magnet
        dy = 30. * cm;
        dz = 23. * cm; // 23. cm
        dx1 = 34. * cm, dy1 = 22. * cm;
        G4double alf = 6. * deg;

        z_pos = 97 * cm - dz / 2.;

        /////

        G4Trap *b_trap = nullptr;
        G4Trap *s_trap = nullptr;

#ifdef WINDOW_E
        // part of vacuum chamber
        b_trap = new G4Trap("b_trap", 2. * cm, 0. * mm + 1. * (z_pos - dz / 2. - (d2_zpos + d2_dz / 2.)),
                            18. * cm + 0.1 * cm, 12. * cm);
        dx = 2 * cm;
        VolVac = new G4LogicalVolume(b_trap, Vacuum, "Vacuum_log", 0, 0, 0);
        VolVac->SetVisAttributes(G4VisAttributes::GetInvisible());
        s_trap =
            new G4Trap("s_trap", (b_trap->GetZHalfLength() - 2.0 * mm) * 2, (b_trap->GetYHalfLength1() + 0.2 * mm) * 2,
                       (b_trap->GetXHalfLength1() - 2. * mm) * 2, (b_trap->GetXHalfLength2() - 2. * mm) * 2);
        Vols2 = new G4SubtractionSolid("Vols2", b_trap, s_trap);


        Volume_log = new G4LogicalVolume(Vols2, Steel, "Volume_log", 0, 0, 0);
        Volume_log->SetVisAttributes(Volume_logVisAtt);
        new G4PVPlacement(0, G4ThreeVector(), Volume_log, "dndnstream_vac", VolVac, false, 0);
        Volume_phys =
            new G4PVPlacement(G4Transform3D(Rotate270X,
                                            G4ThreeVector((b_trap->GetXHalfLength1() + b_trap->GetXHalfLength2()) / 2. -
                                                              3. * cm - 0.2 * cm,
                                                          0., z_pos - dz / 2. - b_trap->GetYHalfLength1())),
                              VolVac, "dndnstream_vac", Hall_log, false, 0);

        // window for e'

        G4double w_th = 0.3 * mm; // 300 микрон
        b_box = new G4Box("bee", 11. * cm / 2., 18. * mm / 2., (10 * mm + w_th) / 2.);
        VolVac = new G4LogicalVolume(b_trap, Vacuum, "Vacuum_log", 0, 0, 0);
        VolVac->SetVisAttributes(G4VisAttributes::GetInvisible());
        s_box = new G4Box("see", b_box->GetXHalfLength() - 1. * mm, b_box->GetYHalfLength() - 1. * mm,
                          b_box->GetXHalfLength() + 0.2 * mm);
        Vols2 = new G4SubtractionSolid("Vols2", b_box, s_box);
        Volume_log = new G4LogicalVolume(Vols2, Steel, "Volume_log", 0, 0, 0);
        Volume_log->SetVisAttributes(Volume_logVisAtt);
        new G4PVPlacement(0, G4ThreeVector(), Volume_log, "dw_vac", VolVac, false, 0);
        s_box = new G4Box("swe", b_box->GetXHalfLength() + 1. * mm, b_box->GetYHalfLength() + 1. * mm, w_th / 2.);
        Volume_log = new G4LogicalVolume(s_box, Beryllium, "Volume_log", 0, 0, 0);
        new G4PVPlacement(0, G4ThreeVector(0., 0., b_box->GetZHalfLength() - w_th / 2.), Volume_log, "dw", VolVac,
                          false, 0);
        Volume_log->SetVisAttributes(Foil_VisAtt);
        Volume_phys = new G4PVPlacement(G4Transform3D(RotateNull,
                                                      G4ThreeVector(b_trap->GetXHalfLength1() * 2 - 3.0 * cm -
                                                                        1.8 * cm - b_box->GetXHalfLength(), // -3.0 * cm
                                                                    0., z_pos - dz / 2.)),
                                        VolVac, "dndnstream_vac", Hall_log, false, 0);

#endif

#ifdef PUMP

        // pumping volume itself
        b_trap = new G4Trap("b_trap", dy, dz, 50. * cm, 50. * cm - dz * tan(alf));
        dx = 2 * cm;
        VolVac = new G4LogicalVolume(b_trap, Vacuum, "Vacuum_log", 0, 0, 0);
        VolVac->SetVisAttributes(G4VisAttributes::GetInvisible());
        s_trap =
            new G4Trap("s_trap", (b_trap->GetZHalfLength() - 3.0 * mm) * 2, (b_trap->GetYHalfLength1() - 3. * mm) * 2,
                       (b_trap->GetXHalfLength1() - 3. * mm) * 2, (b_trap->GetXHalfLength2() - 3. * mm) * 2);
        Vols2 = new G4SubtractionSolid("Vols2", b_trap, s_trap, G4Transform3D(RotateNull, G4ThreeVector(0., 0., 0.)));
        b_box = new G4Box("bbb", 60. * mm / 2., 20. * mm / 2., 3.1 * mm / 2.);
        x_pos = b_trap->GetXHalfLength2() - (b_trap->GetXHalfLength1() - b_trap->GetXHalfLength2()) / 2. -
            b_box->GetXHalfLength() - 3 * mm;
        Vols1 = new G4SubtractionSolid(
            "Vols1", Vols2, b_box,
            G4Transform3D(Rotate270X, G4ThreeVector(x_pos, b_trap->GetYHalfLength1() - 3. * mm / 2., 0.)));

        eltub = new G4EllipticalTube("el1", 36. * mm / 2., 20. * mm / 2., 3.1 * mm / 2);
        x_pos = b_trap->GetXHalfLength1() - dx;
        Vols2 = new G4SubtractionSolid(
            "Vols2", Vols1, eltub,
            G4Transform3D(Rotate270X,
                          G4ThreeVector(x_pos - 1.15 * cm, -(b_trap->GetYHalfLength1() - 3. * mm / 2.), 0.)));

        Volume_log = new G4LogicalVolume(Vols2, Steel, "Volume_log", 0, 0, 0);
        new G4PVPlacement(0, G4ThreeVector(), Volume_log, "dndnstream_vac", VolVac, false, 0);
        //   x_pos *=cos(alf);//z_pos+=x_pos*sin(alf);
        // dx = 3.5 * cm;
        dx = 1.75 * cm;
        Volume_phys = new G4PVPlacement(G4Transform3D(Rotate270X, G4ThreeVector(-x_pos + dx, 0., z_pos)), VolVac,
                                        "dndnstream_vac", Hall_log, false, 0);
        Volume_log->SetVisAttributes(Volume_logVisAtt);


        auto rot_for_world_target2 = new G4RotationMatrix();
        rot_for_world_target2->rotateY(-theta0_GEM);
        // rot_for_world_target2 = &RotateNull;

        // exit elliptical vacuum tube
        eltub = new G4EllipticalTube("el2o", 36. * mm / 2., 20. * mm / 2., 160 * mm / 2);
        eltub1 =
            new G4EllipticalTube("el2i", eltub->GetDx() - 2. * mm, eltub->GetDy() - 2. * mm, eltub->GetDz() + 0.2 * mm);
        VolVac = new G4LogicalVolume(eltub, Vacuum, "Vacuum_log", 0, 0, 0);
        VolVac->SetVisAttributes(G4VisAttributes::GetInvisible());
        Vols1 = new G4SubtractionSolid("Vols1", eltub, eltub1);
        //  			G4Transform3D(RotateNull, G4ThreeVector(x_pos, b_trap->GetYHalfLength1()-3.*mm/2., 0.)));
        Volume_log = new G4LogicalVolume(Vols1, Steel, "Volume_log", 0, 0, 0);
        new G4PVPlacement(G4Transform3D(RotateNull, G4ThreeVector()), Volume_log, "dddstream_vac", VolVac, false, 0);
        z_pos += b_trap->GetYHalfLength1() + eltub->GetDz();

        // я задаю трубку в другом месте
        // Volume_phys = new G4PVPlacement(rot_for_world_target2, G4ThreeVector(dx - 0.75 * cm, 0., z_pos),
        //                                VolVac, "dndnstream_vac", Hall_log, false, 0);
        Volume_log->SetVisAttributes(Volume_logVisAtt);


#endif

#ifdef E_PLASTIC

        /// E-arm sandwich
        ///

        G4Colour steel_col(0., 0.0, 1.0);
        G4Colour convertor_col(0.0, 0.5, 0.5);
        G4Colour plastic_col(1.0, 0.0, 1.0);

#define X_SHIFT 99. * mm

        G4double Lg = -3 * cm;
        dz = 7 * cm;
        x_pos = 70. * mm + 127. * mm / 2.;
        z_pos = 95. * cm;

        G4Box *LQbox = new G4Box("LQbox", 430 * mm / 2., 130. * mm / 2., dz / 2.);
        G4LogicalVolume *LQbox_log = new G4LogicalVolume(LQbox, Air, "LQbox_log", 0, 0, 0);
        LQbox_log->SetVisAttributes(G4VisAttributes::GetInvisible());
        new G4PVPlacement(G4Transform3D(Rotate270Z, G4ThreeVector(x_pos, 0, z_pos)), LQbox_log, "LQbox", Hall_log,
                          false, 0);


        G4double a = 311. - 99.;
        G4Trap *LQiron1 = new G4Trap("LQiron1", 3. * mm, 127. * mm, (99. + a) * mm, (41. + a) * mm);
        Volume_log = new G4LogicalVolume(LQiron1, Steel, "LQiron_log", 0, 0, 0);
        Volume_log->SetVisAttributes(new G4VisAttributes(steel_col));

        x_pos = a * mm / 2. + 99. * mm - (99. + 41.) * mm / 4.;
        x_pos = X_SHIFT - x_pos;
        new G4PVPlacement(0, G4ThreeVector(x_pos, 0, Lg), Volume_log, "LQiron1", LQbox_log, false, 0);

        G4double x1 = 75. * mm / 2., x2 = 169. * mm / 2., y1 = 110. * mm / 2., y2 = y1 - 10 * mm;

        auto ubox = new G4Box("SBox", x2 + 2 * mm, y1 + 2 * mm, dz);

        std::vector<G4TwoVector> polygon(6);
        polygon.clear();
        polygon.push_back(G4TwoVector(-x2, -y1));
        polygon.push_back(G4TwoVector(-x2, -y2));
        polygon.push_back(G4TwoVector(-x1, y1));
        polygon.push_back(G4TwoVector(x1, y1));
        polygon.push_back(G4TwoVector(x2, -y2));
        polygon.push_back(G4TwoVector(x2, -y1));

        G4ExtrudedSolid *LQtun =
            new G4ExtrudedSolid("LQ_tun1", polygon, 10. * mm / 2., G4TwoVector(0., 0.), 1., G4TwoVector(0., 0.), 1.);
        G4LogicalVolume *LQtun_log = new G4LogicalVolume(LQtun, Tungsten, "LQ_tun1log", 0, 0, 0);
        LQtun_log->SetVisAttributes(new G4VisAttributes(convertor_col));

        Lg += 5. * mm + 10. * mm / 2.;
        x_pos = X_SHIFT - (9. + 169.5 / 2.) * mm;
        new G4PVPlacement(0, G4ThreeVector(x_pos, 0, Lg), LQtun_log, "LQtun1", LQbox_log, false, 0);

        LQtun = new G4ExtrudedSolid("LQ_tun2", polygon, 5. * mm / 2., G4TwoVector(0., 0.), 1., G4TwoVector(0., 0.), 1.);
        LQtun_log = new G4LogicalVolume(LQtun, Tungsten, "LQ_tun2log", 0, 0, 0);
        LQtun_log->SetVisAttributes(new G4VisAttributes(convertor_col));

        G4Box *LQlayer = new G4Box("LQlayer1", 200. * mm / 2., 110. * mm / 2., 5. * mm / 2.);
        G4LogicalVolume *LQlayer_log = new G4LogicalVolume(LQlayer, Air, "LQlayer1_log", 0, 0, 0);
        LQlayer_log->SetVisAttributes(G4VisAttributes::GetInvisible());

        G4LogicalVolume *LQscints[7];
        G4int n = 0; //(ARM_LQ|PL_LQ_1)<<DATA_SHIFT;
        for (G4int i = 0; i < 7; i++)
        {
            G4double lng = (90. + 13. * i) * mm;
            ubox = new G4Box("LQplst", lng / 2., 15. * mm / 2., 4.5 * mm / 2.);
            LQscints[i] = new G4LogicalVolume(ubox, Scintil, "LQscint_log", 0, 0, 0);
            LQscints[i]->SetVisAttributes(new G4VisAttributes(plastic_col));
        }
        for (G4int i = 0; i < 7; i++)
        {
            G4double xx[7] = {52.6, 43.1, 40.8, 35.0, 28.0, 20.2, 13.7};
            G4double lng = 90. + 13. * i;
            y_pos = (3 - i) * 15. * mm;
            x_pos = X_SHIFT - (xx[i] + lng / 2.) * mm;
            new G4PVPlacement(0, G4ThreeVector(x_pos, y_pos, 0), LQscints[i], "LQscint", LQlayer_log, false, 1);
        }

        Lg += 10. * mm / 2. + 5. * mm / 2.;
        new G4PVPlacement(0, G4ThreeVector(0, 0, Lg), LQlayer_log, "LQlayer1", LQbox_log, false, n);

        LQtun = new G4ExtrudedSolid("LQ_tun2", polygon, 5. * mm / 2., G4TwoVector(0., 0.), 1., G4TwoVector(0., 0.), 1.);
        LQtun_log = new G4LogicalVolume(LQtun, Tungsten, "LQ_tun2log", 0, 0, 0);
        LQtun_log->SetVisAttributes(new G4VisAttributes(convertor_col));

        Lg += 5. * mm / 2. + 5. * mm / 2.;
        x_pos = X_SHIFT - (8.5 + 169.5 / 2.) * mm;
        new G4PVPlacement(0, G4ThreeVector(x_pos, 0, Lg), LQtun_log, "LQtun2", LQbox_log, false, 0);
        //---------
        LQlayer_log = new G4LogicalVolume(LQlayer, Air, "LQlayer2_log", 0, 0, 0);
        LQlayer_log->SetVisAttributes(G4VisAttributes::GetInvisible());
        for (G4int i = 0; i < 7; i++)
        {
            G4double xx[7] = {52.1, 43.6, 41.1, 35.0, 28.2, 21.0, 13.0};
            G4double lng = 90. + 13. * i;
            y_pos = (3 - i) * 15. * mm;
            x_pos = X_SHIFT - (xx[i] + lng / 2.) * mm;
            new G4PVPlacement(0, G4ThreeVector(x_pos, y_pos, 0), LQscints[i], "LQscint", LQlayer_log, false, 1);
        }
        Lg += 5. * mm / 2. + 5. * mm / 2.;
        new G4PVPlacement(0, G4ThreeVector(0, 0, Lg), LQlayer_log, "LQlayer2", LQbox_log, false, n);
        //--------
        Lg += 5. * mm / 2. + 5. * mm / 2.;
        x_pos = X_SHIFT - (9. + 169.5 / 2.) * mm;
        new G4PVPlacement(0, G4ThreeVector(x_pos, 0, Lg), LQtun_log, "LQtun3", LQbox_log, false, 0);
        //---------
        LQlayer_log = new G4LogicalVolume(LQlayer, Air, "LQlayer3_log", 0, 0, 0);
        LQlayer_log->SetVisAttributes(G4VisAttributes::GetInvisible());
        for (G4int i = 0; i < 7; i++)
        {
            G4double xx[7] = {51.6, 44.6, 40.8, 35.2, 28.1, 21.5, 12.0};
            G4double lng = 90. + 13. * i;
            y_pos = (3 - i) * 15. * mm;
            x_pos = X_SHIFT - (xx[i] + lng / 2.) * mm;
            new G4PVPlacement(0, G4ThreeVector(x_pos, y_pos, 0), LQscints[i], "LQscint", LQlayer_log, false, 1);
        }
        Lg += 5. * mm / 2. + 5. * mm / 2.;
        new G4PVPlacement(0, G4ThreeVector(0, 0, Lg), LQlayer_log, "LQlayer3", LQbox_log, false, n);
        //----------

        Lg += 5. * mm / 2. + 5. * mm / 2.;
        x_pos = X_SHIFT - (9. + 169.5 / 2.) * mm;
        new G4PVPlacement(0, G4ThreeVector(x_pos, 0, Lg), LQtun_log, "LQtun4", LQbox_log, false, 0);
        //------------
        LQlayer_log = new G4LogicalVolume(LQlayer, Air, "LQlayer4_log", 0, 0, 0);
        LQlayer_log->SetVisAttributes(G4VisAttributes::GetInvisible());
        for (G4int i = 0; i < 7; i++)
        {
            double xx[7] = {51.5, 45.8, 40.8, 35.8, 28.5, 20.8, 11.9};
            double lng = 90. + 13. * i;
            y_pos = (3 - i) * 15. * mm;
            x_pos = X_SHIFT - (xx[i] + lng / 2.) * mm;
            new G4PVPlacement(0, G4ThreeVector(x_pos, y_pos, 0), LQscints[i], "LQscint", LQlayer_log, false, 1);
        }
        Lg += 5. * mm / 2. + 5. * mm / 2.;
        new G4PVPlacement(0, G4ThreeVector(0, 0, Lg), LQlayer_log, "LQlayer4", LQbox_log, false, n);


#endif

#ifdef GAMMA_PLASTIC

#define BGO_IND 16
        G4Colour convertor_col(0.0, 0.5, 0.5);

        auto Convertor_VisAtt = new G4VisAttributes(convertor_col);
        G4VisAttributes *BGO_VisAtt = new G4VisAttributes(G4Colour(0.75, 0.5, 0.5));

        ///-------------- 0deg - gamma detector -----------------------------
        // steel wall of vacuum chamber
        G4double z_wall = 97. * cm;
        dz = 2 * mm; // wall thickness
        auto ubox = new G4Box("GWall", 6 * cm / 2., 10 * cm / 2., dz / 2.);
        Volume_log = new G4LogicalVolume(ubox, Steel, "GWall_log", 0, 0, 0);
        Volume_log->SetVisAttributes(Steel_VisAtt);
        z_pos = z_wall;
        //   vol_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,z_pos), Volume_log, "GWall_phys", Hall_log, false, -1);
        z_pos += dz / 2;

        // tungsten layer
        dz = 4.0 * mm; // tungsten thickness
        dx = 3. * cm;
        dy = 3. * cm;
        ubox = new G4Box("GW", dx / 2., dy / 2., dz / 2.);
        Volume_log = new G4LogicalVolume(ubox, Steel, "GWall_log", 0, 0, 0);
        Volume_log->SetVisAttributes(Convertor_VisAtt);
        z_pos += dz / 2;
        auto vol_phys =
            new G4PVPlacement(0, G4ThreeVector(0., 0., z_pos + 0.2 * mm), Volume_log, "GW_phys", Hall_log, false, -1);
        z_pos += dz / 2;

        // BGO
        dz = 6.0 * cm;
        dx = 2. * cm;
        dy = 2. * cm;
        ubox = new G4Box("GBGO", dx / 2., dy / 2., dz / 2.);
        G4LogicalVolume *BGO_log = new G4LogicalVolume(ubox, Steel, "GBGO_log", 0, 0, 0);
        BGO_log->SetVisAttributes(BGO_VisAtt);
        z_pos += dz / 2 + 3 * mm;
        G4int n_ind = BGO_IND;
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., z_pos), BGO_log, "GBGO_phys", Hall_log, false, n_ind);


        ///------------------------------------------------------------------
#endif // GAMMA


#ifdef GAMMA_PLASTIC2

        G4double z_wall = 97. * cm;

        // BGO

        G4VisAttributes *BGO_VisAtt = new G4VisAttributes(G4Colour(0.75, 0.5, 0.5));
        // dz=6.0*cm;
        dz = 4.0 * cm;
        dx = 2. * cm;
        dy = 2. * cm;
        auto ubox = new G4Box("GBGO", dx / 2., dy / 2., dz / 2.);
        G4LogicalVolume *BGO_log = new G4LogicalVolume(ubox, Steel, "GBGO_log", 0, 0, 0);
        BGO_log->SetVisAttributes(BGO_VisAtt);

        z_pos = z_wall;
        z_pos += dz / 2 + 3 * mm - 0.6 * cm;

        x_pos = -(ubox->GetXHalfLength() + 3.0 * cm);

        new G4PVPlacement(&Rotate90Y, G4ThreeVector(x_pos, 0., z_pos), BGO_log, "GBGO_phys", Hall_log, false, 0);


        ///------------------------------------------------------------------
#endif //

#endif // MAGNET_D2
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void DetectorConstruction::ConstructMagnetD1()
    {
#ifdef MAGNET_D1
        auto Hall_log = fWorldLV;

        G4Colour stst_col(0.53, 0.53, 0.84, 0.5);
        G4VisAttributes *Volume_logVisAtt = new G4VisAttributes(stst_col);

        ////

        G4LogicalVolume *Volume_log = nullptr;
        G4PVPlacement *Volume_phys = nullptr;

        G4double d1_zpos = -62. * cm, d1_dz = 12. * cm;
        G4double dx = 38.0 * cm, dy = 44.0 * cm;
        G4double dx1 = 12.0 * cm, dy1 = 2.2 * cm; // зазор

        G4double x_pos = 0., y_pos = 0., z_pos = 0.;

        auto b_box = new G4Box("bD1_box", dx / 2, dy / 2, d1_dz / 2);
        auto VolVac = new G4LogicalVolume(b_box, Vacuum, "Vacuum_log", 0, 0, 0);

        auto s_box = new G4Box("sD1_box", dx1 / 2, dy1 / 2, d1_dz / 2 + 0.2 * mm);

        auto Vols = new G4SubtractionSolid(
            "Vols_D1", b_box, s_box, G4Transform3D(RotateNull, G4ThreeVector(dx / 2 - dx1 / 2 - 3.0 * cm, 0.0, 0.)));
        G4LogicalVolume *fLogicDipole = new G4LogicalVolume(Vols, Iron, "DipoleD1_log", 0, 0, 0);
        auto D1dipole_col = G4Colour(0.5, 0.7, 1.0, 0.4);
        G4VisAttributes *D1dipole_Volume_logVisAtt = new G4VisAttributes(D1dipole_col);
        D1dipole_Volume_logVisAtt->SetForceSolid(true);
        fLogicDipole->SetVisAttributes(D1dipole_Volume_logVisAtt);
        new G4PVPlacement(0, G4ThreeVector(), fLogicDipole, "MagD1", VolVac, false, -1);

        // VolVac = new G4LogicalVolume(s_box, Vacuum, "Vacuum_log", 0, 0, 0);
        //   VolVac->SetVisAttributes(G4VisAttributes::GetInvisible());


        fLogicD1gap = new G4LogicalVolume(s_box, Vacuum, "D1gap_log", 0, 0, 0);
        new G4PVPlacement(0, G4ThreeVector(dx / 2 - dx1 / 2 - 3.0 * cm, 0.0, 0.), fLogicD1gap, "D1gap", VolVac, false,
                          0);
        // fLogicD1gap->SetVisAttributes(G4VisAttributes::GetInvisible());

        x_pos = -(dx / 2 - dx1 / 2 - 3.0 * cm) + 3.0 * cm;
        y_pos = 0.0 * cm;
        z_pos = d1_zpos;

        auto rot_for_world_target = new G4RotationMatrix();
        rot_for_world_target->rotateY(alpha0_GEM);
        rot_for_world_target = &RotateNull;
        Volume_phys = new G4PVPlacement(rot_for_world_target, G4ThreeVector(x_pos, y_pos, z_pos), VolVac, "MagD1",
                                        Hall_log, false, -1);

        //

        // магнитное поле в SD

#endif
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void DetectorConstruction::ConstructMagnetD3()
    {
#ifdef MAGNET_D3
        auto Hall_log = fWorldLV;

        G4Colour stst_col(0.53, 0.53, 0.84, 0.5);
        G4VisAttributes *Volume_logVisAtt = new G4VisAttributes(stst_col);

        ////

        G4LogicalVolume *Volume_log = nullptr;
        G4PVPlacement *Volume_phys = nullptr;

        G4double d3_zpos = 134. * cm, d3_dz = 14.5 * cm;
        G4double dx = 38.0 * cm, dy = 44.0 * cm;
        G4double dx1 = 12.0 * cm, dy1 = 2.2 * cm; // зазор

        G4double x_pos = 0., y_pos = 0., z_pos = 0.;

        auto b_box = new G4Box("bD3_box", dx / 2, dy / 2, d3_dz / 2);
        auto VolVac = new G4LogicalVolume(b_box, Vacuum, "Vacuum_log", 0, 0, 0);

        auto s_box = new G4Box("sD3_box", dx1 / 2, dy1 / 2, d3_dz / 2 + 0.2 * mm);

        auto Vols = new G4SubtractionSolid(
            "Vols_D3", b_box, s_box, G4Transform3D(RotateNull, G4ThreeVector(dx / 2 - dx1 / 2 - 3.0 * cm, 0.0, 0.)));
        G4LogicalVolume *fLogicDipole = new G4LogicalVolume(Vols, Iron, "DipoleD3_log", 0, 0, 0);
        auto D3dipole_col = G4Colour(0.5, 0.7, 1.0, 0.4);
        G4VisAttributes *D3dipole_Volume_logVisAtt = new G4VisAttributes(D3dipole_col);
        D3dipole_Volume_logVisAtt->SetForceSolid(true);
        fLogicDipole->SetVisAttributes(D3dipole_Volume_logVisAtt);
        new G4PVPlacement(0, G4ThreeVector(), fLogicDipole, "MagD3", VolVac, false, -1);

        // VolVac = new G4LogicalVolume(s_box, Vacuum, "Vacuum_log", 0, 0, 0);
        //   VolVac->SetVisAttributes(G4VisAttributes::GetInvisible());


        fLogicD3gap = new G4LogicalVolume(s_box, Vacuum, "D3gap_log", 0, 0, 0);
        new G4PVPlacement(0, G4ThreeVector(dx / 2 - dx1 / 2 - 3.0 * cm, 0.0, 0.), fLogicD3gap, "D3gap", VolVac, false,
                          0);
        // fLogicD3gap->SetVisAttributes(G4VisAttributes::GetInvisible());

        x_pos = -(dx / 2 - dx1 / 2 - 3.0 * cm) + 3.4 * cm;
        y_pos = 0.0 * cm;
        z_pos = d3_zpos;

        auto rot_for_world_target = new G4RotationMatrix();
        rot_for_world_target->rotateY(-theta0_GEM);
        rot_for_world_target = &RotateNull;
        Volume_phys = new G4PVPlacement(rot_for_world_target, G4ThreeVector(x_pos, y_pos, z_pos), VolVac, "MagD3",
                                        Hall_log, false, -1);

        //

        // магнитное поле в SD

#endif
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void DetectorConstruction::ConstructMagnetQF()
    {
#ifdef MAGNET_QF

        auto Hall_log = fWorldLV;

        G4Colour stst_col(0.53, 0.53, 0.84, 0.5);
        G4VisAttributes *Volume_logVisAtt = new G4VisAttributes(stst_col);

        ////

        G4LogicalVolume *Volume_log = nullptr;
        G4PVPlacement *Volume_phys = nullptr;

        G4double qf_zpos = -41. * cm, qf_dz = 16.0 * cm;
        G4double dx = 20.4 * cm, dy = 32.2 * cm;
        G4double dx1 = 7.4 * cm, dy1 = 4. * cm; // зазор


        auto b_box = new G4Box("bQF_box", dx / 2, dy / 2, qf_dz / 2);
        auto VolVac = new G4LogicalVolume(b_box, Vacuum, "Vacuum_log", 0, 0, 0);

        auto s_box = new G4Box("sQF_box", dx1 / 2, dy1 / 2, qf_dz / 2 + 0.2 * mm);

        auto Vols =
            new G4SubtractionSolid("Vols_QF", b_box, s_box, G4Transform3D(RotateNull, G4ThreeVector(0.0, 0.0, 0.0)));
        G4LogicalVolume *fLogicDipole = new G4LogicalVolume(Vols, Iron, "QudrupoleQF_log", 0, 0, 0);
        auto QFqudrupole_col = G4Colour(0.7, 0.5, 1.0, 0.4);
        G4VisAttributes *QFqudrupole_Volume_logVisAtt = new G4VisAttributes(QFqudrupole_col);
        QFqudrupole_Volume_logVisAtt->SetForceSolid(true);
        fLogicDipole->SetVisAttributes(QFqudrupole_Volume_logVisAtt);
        new G4PVPlacement(0, G4ThreeVector(), fLogicDipole, "MagQF", VolVac, false, -1);

        // VolVac = new G4LogicalVolume(s_box, Vacuum, "Vacuum_log", 0, 0, 0);
        //   VolVac->SetVisAttributes(G4VisAttributes::GetInvisible());


        fLogicQFgap = new G4LogicalVolume(s_box, Vacuum, "QFgap_log", 0, 0, 0);
        new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.), fLogicQFgap, "QFgap", VolVac, false, 0);
        // fLogicQFgap->SetVisAttributes(G4VisAttributes::GetInvisible());

        x_pos_qf = dx1 / 2 - 2.5 * cm;
        y_pos_qf = 0.0 * cm;
        z_pos_qf = qf_zpos;

        auto rot_for_world_target = new G4RotationMatrix();
        rot_for_world_target->rotateY(alpha0_GEM);
        // rot_for_world_target = &RotateNull;
        Volume_phys = new G4PVPlacement(rot_for_world_target, G4ThreeVector(x_pos_qf, y_pos_qf, z_pos_qf), VolVac,
                                        "MagQF", Hall_log, false, -1);

        //

        // магнитное поле в SD


#endif
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void DetectorConstruction::ConstructMagnetQD()
    {
#ifdef MAGNET_QD

        auto Hall_log = fWorldLV;

        G4Colour stst_col(0.53, 0.53, 0.84, 0.5);
        G4VisAttributes *Volume_logVisAtt = new G4VisAttributes(stst_col);

        ////

        G4LogicalVolume *Volume_log = nullptr;
        G4PVPlacement *Volume_phys = nullptr;

        G4double qd_zpos = 112. * cm, qd_dz = 16.0 * cm;
        G4double dx = 20.4 * cm, dy = 32.2 * cm;
        G4double dx1 = 7.4 * cm, dy1 = 4. * cm; // зазор


        auto b_box = new G4Box("bQD_box", dx / 2, dy / 2, qd_dz / 2);
        auto VolVac = new G4LogicalVolume(b_box, Vacuum, "Vacuum_log", 0, 0, 0);

        auto s_box = new G4Box("sQD_box", dx1 / 2, dy1 / 2, qd_dz / 2 + 0.2 * mm);

        auto Vols =
            new G4SubtractionSolid("Vols_QD", b_box, s_box, G4Transform3D(RotateNull, G4ThreeVector(0.0, 0.0, 0.0)));
        G4LogicalVolume *fLogicDipole = new G4LogicalVolume(Vols, Iron, "QudrupoleQD_log", 0, 0, 0);
        auto QDqudrupole_col = G4Colour(0.7, 0.5, 1.0, 0.4);
        G4VisAttributes *QDqudrupole_Volume_logVisAtt = new G4VisAttributes(QDqudrupole_col);
        QDqudrupole_Volume_logVisAtt->SetForceSolid(true);
        fLogicDipole->SetVisAttributes(QDqudrupole_Volume_logVisAtt);
        new G4PVPlacement(0, G4ThreeVector(), fLogicDipole, "MagQD", VolVac, false, -1);

        // VolVac = new G4LogicalVolume(s_box, Vacuum, "Vacuum_log", 0, 0, 0);
        //   VolVac->SetVisAttributes(G4VisAttributes::GetInvisible());


        fLogicQDgap = new G4LogicalVolume(s_box, Vacuum, "QDgap_log", 0, 0, 0);
        new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.), fLogicQDgap, "QDgap", VolVac, false, 0);
        // fLogicQDgap->SetVisAttributes(G4VisAttributes::GetInvisible());

        x_pos_qd = dx1 / 2 - 2.5 * cm;
        y_pos_qd = 0.0 * cm;
        z_pos_qd = qd_zpos;

        auto rot_for_world_target = new G4RotationMatrix();
        rot_for_world_target->rotateY(-theta0_GEM);
        // rot_for_world_target = &RotateNull;
        Volume_phys = new G4PVPlacement(rot_for_world_target, G4ThreeVector(x_pos_qd, y_pos_qd, z_pos_qd), VolVac,
                                        "MagQD", Hall_log, false, -1);

        //

        // магнитное поле в SD


#endif
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void DetectorConstruction::ConstructGEM()
    {

        auto AirMaterial = G4Material::GetMaterial("Air");
        auto MylarMaterial = G4Material::GetMaterial("Mylar");
        auto SteelMaterial = G4Material::GetMaterial("Steel");
        auto SiO2Material = G4Material::GetMaterial("SiO2");
        auto GasArCO2Material = G4Material::GetMaterial("ArgonCarbonicGas");
        auto GasArIsobutaMaterial = G4Material::GetMaterial("ArgonIsobutaneGas");

        G4NistManager *nist = G4NistManager::Instance();
        G4Material *KaptonMaterial = nist->FindOrBuildMaterial("G4_KAPTON");

        auto CopperMaterial = G4Material::GetMaterial("Copper"); // материал для медной фольги

        // auto GasVCMaterial = G4Material::GetMaterial("ArgonCarbonicGas");
        auto GasVCMaterial = G4Material::GetMaterial("ArgonIsobutaneGas");

        G4Box *vbox_all = new G4Box("box_all_GEM", 60.0 / 2. * mm, 30.0 / 2. * mm, 180.0 / 2. * mm);
        G4LogicalVolume *box_all_log = new G4LogicalVolume(vbox_all, AirMaterial, "box_all_GEM_log", 0, 0, 0);

        G4Box *vbox = new G4Box("box_GEM", 40.0 / 2. * mm, 8.0 / 2. * mm, 160.0 / 2. * mm);
        box_GEM_log = new G4LogicalVolume(vbox, GasArCO2Material, "box_GEM_log", 0, 0, 0);
        // box_log->SetVisAttributes(G4VisAttributes::GetInvisible());

        G4Box *vbox_gap = new G4Box("box_gap_GEM", vbox->GetXHalfLength() - 0.02 * mm, 0.05 / 2. * mm,
                                    vbox->GetZHalfLength() - 0.02 * mm);

        // медь + каптон

        G4Box *vbox_gap_all = new G4Box("box_gap_all_GEM", vbox->GetXHalfLength() - 0.01 * mm,
                                        0.05 / 2. * mm + 0.01 * mm, vbox->GetZHalfLength() - 0.01 * mm);


        // медный слой
        G4Box *vbox_cuprum_gap = new G4Box("box_cuprum_GEM", vbox->GetXHalfLength() - 0.01 * mm, 0.003 / 2. * mm,
                                           vbox->GetZHalfLength() - 0.01 * mm);

        G4LogicalVolume *box_gap_log = new G4LogicalVolume(vbox_gap, KaptonMaterial, "box_gap_GEM_log", 0, 0, 0);

        G4LogicalVolume *box_gap_all_log =
            new G4LogicalVolume(vbox_gap_all, GasArCO2Material, "box_gap_all_GEM_log", 0, 0, 0);


        G4LogicalVolume *box_cuprum_log =
            new G4LogicalVolume(vbox_cuprum_gap, CopperMaterial, "box_cuprum_GEM_log", 0, 0, 0);


        G4Colour cuprum_col(0., 0.0, 1.0, 0.2);
        auto Cuprum_VisAtt = new G4VisAttributes(cuprum_col);
        Cuprum_VisAtt->SetForceSolid(true);
        box_cuprum_log->SetVisAttributes(Cuprum_VisAtt);
        box_gap_log->SetVisAttributes(Alum_VisAtt);


        new G4PVPlacement(
            0, G4ThreeVector(0., -(vbox_gap->GetYHalfLength() + 2.0 * (vbox_cuprum_gap->GetYHalfLength())), 0.),
            box_cuprum_log, "GEMcuprum", box_gap_all_log, false, 0, fCheckOverlaps);

        new G4PVPlacement(
            0, G4ThreeVector(0., (vbox_gap->GetYHalfLength() + 2.0 * (vbox_cuprum_gap->GetYHalfLength())), 0.),
            box_cuprum_log, "GEMcuprum", box_gap_all_log, false, 0, fCheckOverlaps);

        new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), box_gap_log, "GEM_gap_allFrame_phys", box_gap_all_log, false, 1,
                          fCheckOverlaps);

        ///

        // вычитание дырок

        auto *vbox_hole = new G4Tubs("box_hole_GEM", 0., 0.06 * mm, vbox_gap->GetYHalfLength(), 0. * deg, 360. * deg);

        G4LogicalVolume *vbox_hole_log = new G4LogicalVolume(vbox_hole, GasArCO2Material, "box_hole_GEM_log", 0, 0, 0);


        // auto vbox_gap_all_minus_hole = new G4SubtractionSolid("box_gap_all_minus_hole_GEM", vbox_gap_all ,vbox_hole);
        G4VSolid *vbox_gap_all_minus_hole = vbox_gap_all;

        G4String solidName = "";
        G4double xPos = 0.;
        G4double zPos = 0.;

        G4ThreeVector translation_hole(0., 0., 0.);

        G4RotationMatrix *rot_hole = new G4RotationMatrix();
        rot_hole->rotateX(90. * deg);

        // for (G4int i_x = 1; i_x <= 6; i_x++) // 320
        // {
        //      xPos = -(vbox_gap_all->GetXHalfLength())+ i_x * 0.5 * mm;
        //     for (G4int i_z = 1; i_z <= 10; i_z++) // 320
        //     {
        //        // solidName = "box_gap_all_minus_hole_GEM" + std::to_string(i_x) + std::to_string(i_z);
        //
        //             solidName = Form("box_gap_all_minus_hole_GEM_%i_%i", i_x, i_z);
        //          zPos = -(vbox_gap_all->GetZHalfLength())+ i_z * 0.5 * mm;
        //
        //         translation_hole = G4ThreeVector(xPos, 0. , zPos);
        //
        //
        //       //  vbox_gap_all_minus_hole = new G4SubtractionSolid( solidName, vbox_gap_all_minus_hole ,vbox_hole,
        //       rot_hole, translation_hole);
        //     }
        // }

        G4LogicalVolume *box_gap_all_minus_hole_log =
            new G4LogicalVolume(vbox_gap_all_minus_hole, GasArCO2Material, "box_gap_all_minus_hole_GEM_log", 0, 0, 0);


        // другой спосбо, добавление дырок с воздухом

        G4MultiUnion *multiUnion = new G4MultiUnion("MyMultiUnion");

        // G4Transform3D transform_hole= G4Transform3D::Identity;

        // G4RotationMatrix rot_hole_;
        //  rot_hole_.rotateX(90. * deg);
        //  G4Transform3D transform_hole(rot_hole_, translation_hole);

        G4Transform3D baseTr = G4Transform3D::Identity; // без смещения и вращения
        // multiUnion->AddNode(*vbox_gap, baseTr);


        for (G4int i_x = 1; i_x < 10; i_x++) // 320
        {
            xPos = -(vbox_gap_all->GetXHalfLength()) + i_x * 0.5 * mm;
            for (G4int i_z = 1; i_z < 50; i_z++) // 320
            {
                // solidName = "box_gap_all_minus_hole_GEM" + std::to_string(i_x) + std::to_string(i_z);

                solidName = Form("box_gap_all_minus_hole_GEM_%i_%i", i_x, i_z);
                zPos = -(vbox_gap_all->GetZHalfLength()) + i_z * 0.5 * mm;

                translation_hole = G4ThreeVector(xPos, 0., zPos);
                // transform_hole = G4Transform3D (rot_hole_, translation_hole);

                G4RotationMatrix rot_hole_;
                rot_hole_.rotateX(90. * deg);
                G4Transform3D transform_hole(rot_hole_, translation_hole);


                multiUnion->AddNode(*vbox_hole, transform_hole);

                // new G4PVPlacement(rot_hole, translation_hole,
                //                vbox_hole_log,  solidName, box_gap_log, false, i_x + i_z, false);
            }
        }

        multiUnion->Voxelize();


        G4LogicalVolume *vbox_hole_mu_log =
            new G4LogicalVolume(multiUnion, GasArCO2Material, "box_hole_mu_GEM_log", 0, 0, 0);

        // new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),vbox_hole_mu_log,
        //"GEM_gap_allFrame_phys", box_gap_all_log, false, 1, fCheckOverlaps);


        //

        ///

        new G4PVPlacement(0, G4ThreeVector(0., -(vbox->GetYHalfLength()) + vbox_gap->GetYHalfLength() + 2. * mm, 0.),
                          box_gap_all_log, "GEMgap1", box_GEM_log, false, 0, fCheckOverlaps);

        new G4PVPlacement(
            0,
            G4ThreeVector(0., -(vbox->GetYHalfLength()) + 3. * (vbox_gap->GetYHalfLength()) + 2. * mm + 1.5 * mm, 0.),
            box_gap_all_log, "GEMgap2", box_GEM_log, false, 0, fCheckOverlaps);

        new G4PVPlacement(
            0,
            G4ThreeVector(
                0., -(vbox->GetYHalfLength()) + 5. * (vbox_gap->GetYHalfLength()) + 2. * mm + 1.5 * mm + 1.5 * mm, 0.),
            box_gap_all_log, "GEMgap3", box_GEM_log, false, 0, fCheckOverlaps);


        new G4PVPlacement(&Rotate180Y, G4ThreeVector(0., 0., 0.), box_GEM_log, "GEM_allFrame_phys", box_all_log, false,
                          0, fCheckOverlaps);


        auto rot = new G4RotationMatrix();
        rot->rotateX(90 * deg);
        rot->rotateY(90 * deg);

        auto pos1 = G4ThreeVector((X10_GEM + vbox->GetZHalfLength() - 0.0 * cm), 0,
                                  L0_GEM + D_GEM + Z1_GEM + vbox->GetYHalfLength() + 0 * mm);
        auto pos2 = G4ThreeVector((X20_GEM + vbox->GetZHalfLength() - 0.0 * cm), 0,
                                  L0_GEM + D_GEM + Z1_GEM + L1_GEM + vbox->GetYHalfLength() + 0. * mm);
        auto pos3 = G4ThreeVector((X30_GEM + vbox->GetZHalfLength() - 0.0 * cm), 0,
                                  L0_GEM + D_GEM + Z1_GEM + L1_GEM + L2_GEM + vbox->GetYHalfLength() + 0. * mm);

        new G4PVPlacement(rot, pos1, box_all_log, "GEM1Frame_phys", fWorldLV, false, 1, fCheckOverlaps);

        new G4PVPlacement(rot, pos2, box_all_log, "GEM2Frame_phys", fWorldLV, false, 2, fCheckOverlaps);

        new G4PVPlacement(rot, pos3, box_all_log, "GEM3Frame_phys", fWorldLV, false, 3, fCheckOverlaps);
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.......

    G4LogicalVolume *DetectorConstruction::ConstructLOWQ(G4int nsys = 1) {

#if defined(LOWQ1) || defined(LOWQ2)
        ////////////////////////////////////////////////////////////////////////////////

        auto ScintilMaterial = G4Material::GetMaterial("Scintillator");
        auto AirMaterial = G4Material::GetMaterial("Air");
        auto MylarMaterial = G4Material::GetMaterial("Mylar");
        auto WoodMaterial = G4Material::GetMaterial("G4_CELLULOSE_CELLOPHANE");
        auto ConvertorLQMaterial = G4Material::GetMaterial("G4_Pb");


// Размер Пластика
        G4double pl1_thick = 1.0 * cm;//2.0*cm;	// y-axis
        G4double pl1_width = 50.0 * cm;    // x-axis
        G4double pl1_length = 23.4 * cm; // z-axis

        G4double pl2_thick = 2.0 * cm;//2.0*cm;	// y-axis
        G4double pl2_width = 50.0 * cm;    // x-axis
        G4double pl2_length = 23.5 * cm; // z-axis

#ifdef RUN23
        pl2_thick = 2.0 * cm; //2.0*cm;	// y-axis
        pl1_thick = 2.0 * cm; //2.0*cm;	// y-axis
#endif


        G4double LQ_box_thick = 0.0;	// y-axis
        G4double LQ_box_width = 0.0;    // x-axis
        G4double LQ_box_length = 0.0;    // z-axis

        // Размер Конвертера
        G4double convertor_thick = 1.8 * cm; // y-axis =3*6 mm
        G4double convertor_width = 50.0 * cm;    // x-axis
        G4double convertor_length = 23.5 * cm; // z-axis

#ifdef RUN21
        convertor_thick = 0. * cm;	// y-axis
#endif


        G4Box *ubox;
        G4VPhysicalVolume *vol_phys;

        G4LogicalVolume* LQ_log_nsys1_layer1 = nullptr;
        G4LogicalVolume* LQ_log_nsys1_layer2 = nullptr;
        G4LogicalVolume *LQ_log_nsys2 = nullptr;

        G4LogicalVolume *LQBox_log = nullptr;

        if (nsys==1){

            LQ_box_thick =  (pl1_thick + 1.2 * cm) / 2. +  (pl2_thick + 1.2 * cm) / 2. + convertor_thick;	// y-axis
            LQ_box_width = (pl1_width + 1.0 * cm) / 2.;    // x-axis
            LQ_box_length = (pl1_length + 1.0 * cm) / 2.;    // z-axis

            auto placement1 = G4ThreeVector(0.0, 0.6 * cm, 0.0);
            auto placement2 = G4ThreeVector(0.0, -1.2 * cm, 0.0);


            // Это объем двух сцинтилляторов с пленокой
            ubox = new G4Box("LQ_box", LQ_box_width, LQ_box_thick,LQ_box_length);
            LQBox_log = new G4LogicalVolume(ubox, AirMaterial, "LQBox_log", 0, 0, 0);
            LQBox_log->SetVisAttributes(G4VisAttributes::GetInvisible());

            // Это объем одного сцинтияллтора + пленка
            auto ubox1_cover = new G4Box("LQ_box_cover", pl1_width / 2. + 0.5 * mm, pl1_thick / 2. + 0.5 * mm, pl1_length / 2. + 0.5 * mm);
            auto ubox2_cover = new G4Box("LQ_box_cover", pl2_width / 2. + 0.5 * mm, pl2_thick / 2. + 0.5 * mm, pl2_length / 2. + 0.5 * mm);
            auto LQBoxCover1_log = new G4LogicalVolume(ubox1_cover, AirMaterial, "LQBoxCover1_log", 0, 0, 0);
            auto LQBoxCover2_log = new G4LogicalVolume(ubox2_cover, AirMaterial, "LQBoxCover2_log", 0, 0, 0);
            // LQBoxCover_log->SetVisAttributes(Plastic_VisAtt);

            G4int copy_number = 0;

#ifdef RUN21
                copy_number = 1;
#endif

            vol_phys = new G4PVPlacement(0, placement1,
                                         LQBoxCover1_log, "LQ1", LQBox_log, false, copy_number, fCheckOverlaps);

#if defined(RUN21)
            vol_phys = new G4PVPlacement(0, placement2,
                                         LQBoxCover2_log, "LQ2", LQBox_log, false, copy_number - 1, fCheckOverlaps);
#endif //RUN21

            // Это объем одного сцинтиллятора
            auto ubox1 = new G4Box("LQ1", pl1_width / 2., pl1_thick / 2., pl1_length / 2.);
            auto ubox2 = new G4Box("LQ2", pl2_width / 2., pl2_thick / 2., pl2_length / 2.);
            LQ_log_nsys1_layer1 = new G4LogicalVolume(ubox1, ScintilMaterial, "LQ1_log", 0, 0, 0);
            LQ_log_nsys1_layer1->SetVisAttributes(Plastic_VisAtt);


            vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0),
                                         LQ_log_nsys1_layer1, "LQ1", LQBoxCover1_log, false, 0, fCheckOverlaps);

            LQ_log_nsys1_layer2 = new G4LogicalVolume(ubox2, ScintilMaterial, "LQ2_log", 0, 0, 0);
            LQ_log_nsys1_layer2->SetVisAttributes(Plastic_VisAtt);

#if defined(RUN21)
            vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0),
                                         LQ_log_nsys1_layer2 , "LQ2", LQBoxCover2_log, false, copy_number - 1, fCheckOverlaps);
#endif //RUN21


#if defined(RUN23) && defined(LOWQ_CONVERTOR)

            // Это объем конвертера перед сцинтиллятором
            auto convertor_box = new G4Box("CONVERTOR_LQ", convertor_width / 2., convertor_thick / 2., convertor_length / 2.);
            auto convertorLV = new G4LogicalVolume(convertor_box, ConvertorLQMaterial, "CONVERTOR_LQ_log", 0, 0, 0);
            convertorLV->SetVisAttributes(Convertor_LQ_VisAtt);

            auto converterPV = new G4PVPlacement(0,
                                                 G4ThreeVector(0.0, -(placement1.y() + pl1_thick / 2. + convertor_thick / 2. + 0.5 * cm - 1.*cm), 0.0),
                                                 convertorLV, "CONVERTOR_LQ_phys", LQBox_log, false, 0, fCheckOverlaps);
#endif //RUN23


            // Это объем пленки
//   G4VisAttributes *ProCover_VisAtt = new G4VisAttributes(blackpaper_col);
            ubox = new G4Box("LQCover1", pl1_width / 2., 0.15 / 2. * mm, pl1_length / 2.);
            G4LogicalVolume *LQCover1_log = new G4LogicalVolume(ubox, MylarMaterial, "LQCover1_log", 0, 0, 0);
            //   LQCover_log->SetVisAttributes(ProCover_VisAtt);
            LQCover1_log->SetVisAttributes(new G4VisAttributes(G4Color(0.1, 0.1, 0.5)));
            new G4PVPlacement(0, G4ThreeVector(0.0, +pl1_thick / 2. + 0.1 * mm, 0.0),
                              LQCover1_log, "LQCover1_phys", LQBoxCover1_log, false, -1, fCheckOverlaps);
            new G4PVPlacement(0, G4ThreeVector(0.0, -pl1_thick / 2. - 0.1 * mm, 0.0),
                              LQCover1_log, "LQCover1_phys", LQBoxCover1_log, false, -1, fCheckOverlaps);


            ubox = new G4Box("LQCover2", pl2_width / 2., 0.15 / 2. * mm, pl2_length / 2.);
            G4LogicalVolume *LQCover2_log = new G4LogicalVolume(ubox, MylarMaterial, "LQCover2_log", 0, 0, 0);
            //   LQCover_log->SetVisAttributes(ProCover_VisAtt);
            LQCover2_log->SetVisAttributes(new G4VisAttributes(G4Color(0.1, 0.1, 0.5)));

#if defined(RUN21)
            new G4PVPlacement(0, G4ThreeVector(0.0, +pl2_thick / 2. + 0.1 * mm, 0.0),
                              LQCover2_log, "LQCover2_phys", LQBoxCover2_log, false, -1, fCheckOverlaps);
            new G4PVPlacement(0, G4ThreeVector(0.0, -pl2_thick / 2. - 0.1 * mm, 0.0),
                              LQCover2_log, "LQCover2_phys", LQBoxCover2_log, false, -1, fCheckOverlaps);
#endif //RUN21

        }

/////////////
/////////////


        if (nsys == 2)
        {
            // Размер Дерева
            G4double wood_thick = 4.5 * cm; // y-axis
            G4double wood_width = 50.0 * cm; // x-axis
                G4double wood_length = 23.5 * cm; // z-axis

#if defined(RUN23)
                wood_thick = 0.5 * cm; // y-axis
#endif

            LQ_box_thick = 2 * (pl1_thick + 1.2 * cm) / 2. + wood_thick + convertor_thick;	// y-axis
            LQ_box_width = (pl1_width + 1.0 * cm) / 2.;    // x-axis
            LQ_box_length = (pl1_length + 1.0 * cm) / 2.;    // z-axis



            auto placement1 = G4ThreeVector(0.0, -0.6 * cm, 0.0);
            auto placement2 = G4ThreeVector(0.0, 0.6 * cm + wood_thick + 1.0*cm + convertor_thick, 0.0);


            // Это объем двух сцинтилляторов с пленокой
            ubox = new G4Box("LQ_box", LQ_box_width, LQ_box_thick, LQ_box_length);
            LQBox_log = new G4LogicalVolume(ubox, AirMaterial, "LQBox_log", 0, 0, 0);
            LQBox_log->SetVisAttributes(G4VisAttributes::GetInvisible());

            // Это объем одного сцинтияллтора + пленка
            ubox = new G4Box("LQ_box_cover", pl1_width / 2. + 0.5 * mm, pl1_thick / 2. + 0.5 * mm, pl1_length / 2. + 0.5 * mm);
            auto LQBoxCover_log = new G4LogicalVolume(ubox, AirMaterial, "LQBoxCover_log", 0, 0, 0);
            // LQBoxCover_log->SetVisAttributes(Plastic_VisAtt);

                G4int copy_number = 0;

#ifdef RUN21
                copy_number = 1;
#endif


                vol_phys = new G4PVPlacement(0, placement1, LQBoxCover_log, "LQ1", LQBox_log, false, copy_number,
                                             fCheckOverlaps);

#if defined(RUN21) // два сцинтиллятора только в заходе 21 года
            vol_phys = new G4PVPlacement(0, placement2,
                                         LQBoxCover_log, "LQ2", LQBox_log, false, copy_number - 1, fCheckOverlaps);
#endif //RUN21
            // Это объем одного сцинтияллтора
            ubox = new G4Box("LQ", pl1_width / 2., pl1_thick / 2., pl1_length / 2.);
            LQ_log_nsys2 = new G4LogicalVolume(ubox, ScintilMaterial, "LQ_log", 0, 0, 0);
            LQ_log_nsys2 ->SetVisAttributes(Plastic_VisAtt);
            vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), LQ_log_nsys2, "LQ1", LQBoxCover_log, false, 1,
                                         fCheckOverlaps);


            // Это объем дерева между сцинтилляторами
                auto wood_box = new G4Box("WOOD_LQ", wood_width / 2., wood_thick / 2., wood_length / 2.);
            auto woodLV = new G4LogicalVolume(wood_box, WoodMaterial, "WOOD_LQ_log", 0, 0, 0);
                woodLV->SetVisAttributes(WOOD_VisAtt);

                // дерево только для верхней системы в заходе 2021

#if defined(RUN21)
            auto woodPV = new G4PVPlacement(0,
                                                G4ThreeVector(0.0, placement1.y() + pl1_thick / 2. + wood_thick / 2. + 0.5 * cm, 0.0),
                                                woodLV, "WOOD_LQ_phys", LQBox_log, false, 0, fCheckOverlaps);
#endif // RUN21

#if defined(RUN23) && defined(LOWQ_CONVERTOR)

                // Это объем конвертера перед сцинтиллятором
            auto convertor_box = new G4Box("CONVERTOR_LQ", convertor_width / 2., convertor_thick / 2., convertor_length / 2.);
            auto convertorLV = new G4LogicalVolume(convertor_box, ConvertorLQMaterial, "CONVERTOR_LQ_log", 0, 0, 0);
            convertorLV->SetVisAttributes(Convertor_LQ_VisAtt);

            auto converterPV = new G4PVPlacement(0,
                G4ThreeVector(0.0, -(placement1.y() + pl1_thick / 2. + convertor_thick / 2. + 0.5 * cm + 1.*cm), 0.0),
                                                 convertorLV, "CONVERTOR_LQ_phys", LQBox_log, false, 0, fCheckOverlaps);
#endif //RUN23




            // Это объем пленки
            //   G4VisAttributes *ProCover_VisAtt = new G4VisAttributes(blackpaper_col);
            ubox = new G4Box("LQCover", pl1_width / 2., 0.15 / 2. * mm, pl1_length / 2.);
            G4LogicalVolume *LQCover_log = new G4LogicalVolume(ubox, MylarMaterial, "LQCover_log", 0, 0, 0);
            // LQCover_log->SetVisAttributes(ProCover_VisAtt);
            LQCover_log->SetVisAttributes(new G4VisAttributes(G4Color(0.1, 0.1, 0.5)));
            new G4PVPlacement(0, G4ThreeVector(0.0, +pl1_thick / 2. + 0.1 * mm, 0.0),
                              LQCover_log, "LQCover_phys", LQBoxCover_log, false, -1, fCheckOverlaps);
            new G4PVPlacement(0, G4ThreeVector(0.0, -pl1_thick / 2. - 0.1 * mm, 0.0),
                              LQCover_log, "LQCover_phys", LQBoxCover_log, false, -1, fCheckOverlaps);
        }


        if (nsys == 1)
        {
            scint_LQ_nsys1LV_layer1 = LQ_log_nsys1_layer1;
                scint_LQ_nsys1LV_layer2 = LQ_log_nsys1_layer2;
        }
        if (nsys == 2) scint_LQ_nsys2LV = LQ_log_nsys2;


       return LQBox_log;


#endif // LOWQ

    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void DetectorConstruction::ConstructSDandField() {
        auto sdManager = G4SDManager::GetSDMpointer();
        G4String SDname;
        sdManager->SetVerboseLevel(1);
        //
        // Sensitive detectors
        //

        // BEGIN sensitive detector for proton plastic


        // if(!aplasticSD)
        // {
        //   auto aplasticSD = new PlasticSD("plasticSD", "plasticHitsCollection", fNofLayers_plastic_fat_nsys1);
        //   auto aplasticSD = new PlasticSD("plasticSD", "plasticHitsCollection", this);


        //  }

 auto det_constr = this;

#ifndef isGenLQ
#ifdef PF1_FAT
        //  SetSensitiveDetector("plastic_fat_nsys1LV", aplasticSD);
        auto aplastic_fat_nsys1SD =
     new PlasticSD(SDname = "/plastic_fat_nsys1SD", "plastic_fat_nsys1HitsCollection", fNofLayers_plastic_fat_nsys1);
        aplastic_fat_nsys1SD->SetDiscrThres(plastic_fat_threshold_TDC);
        sdManager->AddNewDetector(aplastic_fat_nsys1SD);
        plastic_fat_nsys1LV->SetSensitiveDetector(aplastic_fat_nsys1SD);
#endif

#ifdef PF2_FAT
        // SetSensitiveDetector("plastic_fat_nsys2LV", aplasticSD);
        auto aplastic_fat_nsys2SD = new PlasticSD(SDname = "/plastic_fat_nsys2SD", "plastic_fat_nsys2HitsCollection",
                                                  fNofLayers_plastic_fat_nsys2);
        aplastic_fat_nsys2SD->SetDiscrThres(plastic_fat_threshold_TDC);
        sdManager->AddNewDetector(aplastic_fat_nsys2SD);
        plastic_fat_nsys2LV_120->SetSensitiveDetector(aplastic_fat_nsys2SD);
        plastic_fat_nsys2LV_125->SetSensitiveDetector(aplastic_fat_nsys2SD);
#endif
#endif // isGenLQ


#ifdef PF1_THIN
        // SetSensitiveDetector("plastic_fat_nsys2LV", aplasticSD);
        auto aplastic_thin_nsys1SD = new PlasticSD(SDname = "/plastic_thin_nsys1SD", "plastic_thin_nsys1HitsCollection",
                                                   fNofLayers_plastic_thin_nsys1);
        aplastic_thin_nsys1SD->SetDiscrThres(plastic_thin_threshold_TDC);
        sdManager->AddNewDetector(aplastic_thin_nsys1SD);
        plastic_thin_nsys1LV->SetSensitiveDetector(aplastic_thin_nsys1SD);
#endif

#ifdef PF2_THIN
        // SetSensitiveDetector("plastic_fat_nsys2LV", aplasticSD);
        auto aplastic_thin_nsys2SD = new PlasticSD(SDname = "/plastic_thin_nsys2SD", "plastic_thin_nsys2HitsCollection",
                                                   fNofLayers_plastic_thin_nsys2);
        aplastic_thin_nsys2SD->SetDiscrThres(plastic_thin_threshold_TDC);
        sdManager->AddNewDetector(aplastic_thin_nsys2SD);
        plastic_thin_nsys2LV->SetSensitiveDetector(aplastic_thin_nsys2SD);
#endif


#ifdef LOWQ1

        auto aplastic_LQ_nsys1SD = new PlasticSD(SDname = "/plastic_LQ_nsys1SD", "plastic_LQ_nsys1HitsCollection",
                                                   fNofLayers_plastic_LQ_nsys1);
        aplastic_LQ_nsys1SD->SetDiscrThres(plastic_LQ_threshold_TDC);
        sdManager->AddNewDetector(aplastic_LQ_nsys1SD);
        scint_LQ_nsys1LV_layer1->SetSensitiveDetector(aplastic_LQ_nsys1SD);

#ifdef RUN21
        scint_LQ_nsys1LV_layer2->SetSensitiveDetector(aplastic_LQ_nsys1SD);
#endif



#endif // LOWQ1


#ifdef LOWQ2

        auto aplastic_LQ_nsys2SD = new PlasticSD(SDname = "/plastic_LQ_nsys2SD", "plastic_LQ_nsys2HitsCollection",
                                                 fNofLayers_plastic_LQ_nsys2);
        aplastic_LQ_nsys2SD->SetDiscrThres(plastic_LQ_threshold_TDC);
        sdManager->AddNewDetector(aplastic_LQ_nsys2SD);
        scint_LQ_nsys2LV->SetSensitiveDetector(aplastic_LQ_nsys2SD);

#endif // LOWQ2


#ifndef isGenLQ
#ifdef HADCAL1
        auto ahadron_calorimeter_nsys1SD = new HadronCalorimeterSD(SDname = "/hadron_calorimeter_nsys1SD",
                                                                   "hadron_calorimeter_nsys1HitsCollection", 1);
        ahadron_calorimeter_nsys1SD->SetDiscrThres(HadronCalorimeter_threshold_TDC);
        sdManager->AddNewDetector(ahadron_calorimeter_nsys1SD);
        scint_HadCal_nsys1LV->SetSensitiveDetector(ahadron_calorimeter_nsys1SD);
#endif


#ifdef HADCAL2
        auto ahadron_calorimeter_nsys2SD = new HadronCalorimeterSD(SDname = "/hadron_calorimeter_nsys2SD",
                                                                   "hadron_calorimeter_nsys2HitsCollection", 2);
        ahadron_calorimeter_nsys2SD->SetDiscrThres(HadronCalorimeter_threshold_TDC);
        sdManager->AddNewDetector(ahadron_calorimeter_nsys2SD);
        scint_HadCal_nsys2LV->SetSensitiveDetector(ahadron_calorimeter_nsys2SD);
#endif
#endif // isGenLQ


        // wire chambers
#ifdef DCARM1
        auto aW_chamber_nsys1SD = new ChamberSD(SDname = "/W_chamber_nsys1SD", "W_Chamber_nsys1HitsCollection", 1);
        aW_chamber_nsys1SD->SetDiscrThres(WChamber_threshold_TDC);
        sdManager->AddNewDetector(aW_chamber_nsys1SD);

        WCTheta1_gas_nsys1LV->SetSensitiveDetector(aW_chamber_nsys1SD);
        WCPhi1_gas_nsys1LV->SetSensitiveDetector(aW_chamber_nsys1SD);
        WCTheta2_gas_nsys1LV->SetSensitiveDetector(aW_chamber_nsys1SD);
#endif

#ifdef DCARM2
        auto aW_chamber_nsys2SD = new ChamberSD(SDname = "/W_chamber_nsys2SD", "W_Chamber_nsys2HitsCollection", 2);
        aW_chamber_nsys2SD->SetDiscrThres(WChamber_threshold_TDC);
        sdManager->AddNewDetector(aW_chamber_nsys2SD);

        WCTheta1_gas_nsys2LV->SetSensitiveDetector(aW_chamber_nsys2SD);
        WCPhi1_gas_nsys2LV->SetSensitiveDetector(aW_chamber_nsys2SD);
        WCTheta2_gas_nsys2LV->SetSensitiveDetector(aW_chamber_nsys2SD);
#endif


#if defined(VCARM1) && defined(RUN21)
        auto aV_chamber_nsys1SD = new ChamberSD(SDname = "/V_chamber_nsys1SD", "V_Chamber_nsys1HitsCollection", 1);
        aV_chamber_nsys1SD->SetDiscrThres(VChamber_threshold_TDC);
        sdManager->AddNewDetector(aV_chamber_nsys1SD);

        VCGas_log_nsys1LV->SetSensitiveDetector(aV_chamber_nsys1SD);
#endif
#if defined(VCARM2) && defined(RUN21)
        auto aV_chamber_nsys2SD = new ChamberSD(SDname = "/V_chamber_nsys2SD", "V_Chamber_nsys2HitsCollection", 2);
        aV_chamber_nsys2SD->SetDiscrThres(VChamber_threshold_TDC);
        sdManager->AddNewDetector(aV_chamber_nsys2SD);

        VCGas_log_nsys2LV->SetSensitiveDetector(aV_chamber_nsys2SD);
#endif


        // END sensitive detector for proton plastic


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


    G4int GetType(G4int n) {
        n = n % ARM2_IND;
        if (n < VC_IND) return DT_SCINT;
        if (n < LQ_IND) return DT_WIRE;
        return DT_EMCAL;
    }

    G4int GetType2(G4int n) {

        if (n == 11) return DT_SCINT_THIN;
        if (n == 22) return DT_SCINT_FAT;
        if (n == 22) return DT_SCINT_LQ;

        return DT_SCINT_THIN;
    }

    G4double GetAttenuL(G4int n) {
        G4double r = 1000. * mm;
        if (GetType(n) == DT_SCINT) r = 2000. * mm;

        return r;
    }

    G4double GetDiscrThr(G4int n) {
        G4double r = 0.25 * MeV;

        if (GetType(n) == DT_SCINT) r = 1.0 * MeV;
        if (GetType(n) == DT_WIRE) r = 0.1 * keV;
        if (GetType(n) == DT_EMCAL) r = 0.5 * MeV;


        if (GetType2(n) == DT_SCINT_THIN) r = 0.5 * MeV;
        if (GetType2(n) == DT_SCINT_FAT) r = 15. * MeV;
        if (GetType2(n) == DT_SCINT_LQ) r = 0.5 * MeV;

        return r;
    }


    G4double GetDiscrThr2(G4int n) {
        G4double r = 0.25 * MeV;


        if (GetType2(n) == DT_SCINT_THIN) r = 0.5 * MeV;
        if (GetType2(n) == DT_SCINT_FAT) r = 15. * MeV;
        if (GetType2(n) == DT_SCINT_LQ) r = 0.5 * MeV;

        return r;
    }
}



