//
// Created by Арсений Юрченко on 12.11.2026.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TaggingSystem/forTaggingSystem.hh"

#include "Constants.hh"
#include "DetectorConstruction.hh"

//#include "ChamberSD.hh"
#include "DetectorConstruction.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
//#include "HadronCalorimeterSD.hh"
//#include "PlasticSD.hh"

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

#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>
// #include <TString.h>


#include "PrimaryGeneratorAction.hh"

#include <TMath.h>

#include "PrimaryGeneratorMessenger.hh"

#include "Constants.hh"
#include "DetectorConstruction.hh"
#include "EventInfo.hh"
#include "G4BetheHeitler5DModel.hh"
#include "G4BetheHeitlerModel.hh"
#include "G4Event.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4eBremsstrahlungRelModel.hh"
#include "Randomize.hh"

#include "G4GenericMessenger.hh"
#include "G4ParticleMomentum.hh"
#include "G4ParticleTypes.hh"
#include "G4ThreeVector.hh"


#include "../EcoMug/EcoMug.h"

// #pragma once


namespace Cosmic
{

    // class BetheHeitlerGenerator;
    // class DetectorConstruction;


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

        G4double w_lng = 24.0 * cm; // 21.6*cm;		// titanium window length
        G4double w_shft = -10. * cm;

        G4double dx, dz, d_y1, d_y2;

        G4double w_pos = w_lng / 2 + w_shft;


#ifdef CUBIK

        G4Box *VolumeOut_box = new G4Box("VolumeOut_box1", 6.3 / 2 * cm, 5.0 / 2 * cm, 52.8 / 2 * cm);

        G4Box *VolumeIn_box = new G4Box("VolumeIn_box", (6.3 - 0.6) / 2 * cm, (5.0 - 0.6) / 2 * cm, 52.81 / 2 * cm);

        dx = 52.8 / 2. * cm - (w_pos + (w_lng + 1.4 * cm) / 2.) - 1.0 * cm;
        dz = 52.8 / 2. * cm - 0.5 * cm - dx / 2.;

        d_y1 = 1.0 * cm - 0.4 * mm; // inlet tube outer radius
        d_y2 = d_y1 - 0.2 * cm; // inlet tube inner radius
        G4Tubs *InletHole_tub = new G4Tubs("InletHole_tub", 0.0 * cm, d_y1, 0.31 / 2 * cm, 0 * deg, 360 * deg);
        G4Tubs *Inlet_tub = new G4Tubs("Inlet_tub", d_y2, d_y1, 20. / 2. * cm, 0 * deg, 360 * deg);

        G4SubtractionSolid *Volume_0 = new G4SubtractionSolid("Volume_0", VolumeOut_box, VolumeIn_box);

        G4SubtractionSolid *Volume_1 = new G4SubtractionSolid(
            "Volume_1", Volume_0, InletHole_tub, G4Transform3D(Rotate90Y, G4ThreeVector(3.0 * cm, 0.0 * cm, 0.0 * cm)));

        G4SubtractionSolid *Volume_3;
        Volume_3 = Volume_1;

        G4Box *Cavity_box = new G4Box("Cavity_box", 4.0 / 2 * cm, 5.05 / 2 * cm, w_lng / 2);

        G4SubtractionSolid *Volume = new G4SubtractionSolid("Volume", Volume_3, Cavity_box,
                                                            G4Transform3D(RotateNull, G4ThreeVector(0.0, 0.0, w_pos)));

        Volume_log = new G4LogicalVolume(Volume, SteelMaterial, "Volume_log", 0, 0, 0);

        vol_phys =
            new G4PVPlacement(0, G4ThreeVector(), Volume_log, "Volume_phys", fWorldLV_rot, false, -1, fCheckOverlaps);

        //  G4VisAttributes *Volume_logVisAtt = new G4VisAttributes(steel_col);
        Volume_log->SetVisAttributes(Steel_VisAtt);

        // Размещение трубки для инжекции атомов

        G4LogicalVolume *InTub_log = new G4LogicalVolume(Inlet_tub, AluminMaterial, "InTub_log", 0, 0, 0);
#ifdef INJECT_TUBE

        vol_phys = new G4PVPlacement(G4Transform3D(Rotate90Y180X, G4ThreeVector((1.2 + 10.) * cm, 0.0 * cm, 0.0 * cm)),
                                     InTub_log, "InTub_phys", fWorldLV_rot, false, -1, fCheckOverlaps);

#endif

        G4VisAttributes *Alum_VisAtt_ = new G4VisAttributes(*Alum_VisAtt);
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

        vol_phys = new G4PVPlacement(0, G4ThreeVector(2.1 * cm, -0.5 * mm, 0.0), W1b_log, "W1b", Window_log, false, -1,
                                     fCheckOverlaps);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(-2.1 * cm, -0.5 * mm, 0.0), W1b_log, "W1b", Window_log, false, -1,
                                     fCheckOverlaps);

        G4Box *W2b = new G4Box("W2b", 4.4 / 2 * cm, 2.0 / 2 * mm, 3.0 / 2 * mm);
        G4LogicalVolume *W2b_log = new G4LogicalVolume(W2b, SteelMaterial, "W2b_log", 0, 0, 0);
        W2b_log->SetVisAttributes(Steel_VisAtt);

        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0, -0.5 * mm, (w_lng + 3. * mm) / 2.), W2b_log, "W2b",
                                     Window_log, false, -1, fCheckOverlaps);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0, -0.5 * mm, -(w_lng + 3. * mm) / 2.), W2b_log, "W2b",
                                     Window_log, false, -1, fCheckOverlaps);

        G4Box *W3b = new G4Box("W3b", 6.0 / 2 * mm, 1.0 / 2 * mm, w_lng / 2);
        G4LogicalVolume *W3b_log = new G4LogicalVolume(W3b, SteelMaterial, "W3b_log", 0, 0, 0);
        W3b_log->SetVisAttributes(Steel_VisAtt);

        vol_phys = new G4PVPlacement(0, G4ThreeVector(2.3 * cm, 1.0 * mm, 0.0), W3b_log, "W3b", Window_log, false, -1,
                                     fCheckOverlaps);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(-2.3 * cm, 1.0 * mm, 0.0), W3b_log, "W3b", Window_log, false, -1,
                                     fCheckOverlaps);

        G4Box *W4b = new G4Box("W4b", 5.2 / 2 * cm, 1.0 / 2 * mm, 7.0 / 2 * mm);
        G4LogicalVolume *W4b_log = new G4LogicalVolume(W4b, SteelMaterial, "W4b_log", 0, 0, 0);
        W4b_log->SetVisAttributes(Steel_VisAtt);

        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0, 1.0 * mm, (w_lng + 7. * mm) / 2.), W4b_log, "W4b",
                                     Window_log, false, -1, fCheckOverlaps);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0, 1.0 * mm, -(w_lng + 7. * mm) / 2.), W4b_log, "W4b",
                                     Window_log, false, -1, fCheckOverlaps);


        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0, 2.65 * cm, w_pos), Window_log, "Volume_phys", fWorldLV_rot,
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

        G4Box *CoipperShield_box1 = new G4Box("CoipperShield_box1", 6.3 / 2 * cm, thick_shield / 2,
                                              (52.8 * cm - (w_lng + 1.4 * cm) - w_pos - 2.0 * cm) / 2. / 2.);

        G4Box *CoipperShield_box2 = new G4Box("CoipperShield_box2", 6.3 / 2 * cm, thick_shield / 2,
                                              (52.8 * cm - (w_lng + 1.4 * cm) + w_pos + 1.5 * cm) / 2. / 2.);

        G4Box *CoipperShield_box3 = new G4Box("CoipperShield_box3", 6.3 / 2 * cm, 6.3 / 2 * cm, thick_shield / 2);


        G4Box *CoipperShield_box4 =
            new G4Box("CoipperShield_box4", thick_shield / 2, 6.3 / 2 * cm, 52.8 / 2 * cm + 0.5 * mm);

        G4Box *CoipperShield_box5 =
            new G4Box("CoipperShield_box5", thick_shield / 2, 6.3 / 2 * cm, 52.8 / 4 * cm - d_y1 / 2 - 1.0 * mm);

        // G4cout<< w_pos <<"!!!!" << G4endl;

        G4LogicalVolume *CoipperShield_log1 =
            new G4LogicalVolume(CoipperShield_box1, CopperMaterial, "CoipperShield_log", 0, 0, 0);

        vol_phys =
            new G4PVPlacement(0,
                              G4ThreeVector(0.0 * cm, 2.81 * cm + 0.2 * mm,
                                            w_pos + w_lng / 4 - 4 * cm / 4 + (52.8 * cm - (w_lng + 1.4 * cm)) / 2.),
                              CoipperShield_log1, "CoipperShield_phys", fWorldLV_rot, false, -1, fCheckOverlaps);
        vol_phys =
            new G4PVPlacement(0,
                              G4ThreeVector(0.0 * cm, -2.81 * cm - 0.2 * mm,
                                            w_pos + w_lng / 4 - 4 * cm / 4 + (52.8 * cm - (w_lng + 1.4 * cm)) / 2.),
                              CoipperShield_log1, "CoipperShield_phys", fWorldLV_rot, false, -1, fCheckOverlaps);

        G4LogicalVolume *CoipperShield_log2 =
            new G4LogicalVolume(CoipperShield_box2, CopperMaterial, "CoipperShield_log", 0, 0, 0);

        vol_phys =
            new G4PVPlacement(0,
                              G4ThreeVector(0.0 * cm, 2.81 * cm + 0.2 * mm,
                                            w_pos - w_lng / 4 - 4 * cm / 4 - (52.8 * cm - (w_lng + 1.4 * cm)) / 2.),
                              CoipperShield_log2, "CoipperShield_phys", fWorldLV_rot, false, -1, fCheckOverlaps);
        vol_phys =
            new G4PVPlacement(0,
                              G4ThreeVector(0.0 * cm, -2.81 * cm - 0.2 * mm,
                                            w_pos - w_lng / 4 - 4 * cm / 4 - (52.8 * cm - (w_lng + 1.4 * cm)) / 2.),
                              CoipperShield_log2, "CoipperShield_phys", fWorldLV_rot, false, -1, fCheckOverlaps);

        G4LogicalVolume *CoipperShield_log3 =
            new G4LogicalVolume(CoipperShield_box3, CopperMaterial, "CoipperShield_log", 0, 0, 0);

        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0 * cm, 0.0 * cm, (52.8 + 0.2) / 2 * cm), CoipperShield_log3,
                                     "CoipperShield_phys", fWorldLV_rot, false, -1, fCheckOverlaps);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0 * cm, 0.0 * cm, -(52.8 + 0.2) / 2 * cm), CoipperShield_log3,
                                     "CoipperShield_phys", fWorldLV_rot, false, -1, fCheckOverlaps);


        G4LogicalVolume *CoipperShield_log4 =
            new G4LogicalVolume(CoipperShield_box4, CopperMaterial, "CoipperShield_log", 0, 0, 0);

        vol_phys = new G4PVPlacement(0, G4ThreeVector(-3.15 * cm - 0.2 * mm, 0.0 * cm, 0.0 * cm), CoipperShield_log4,
                                     "CoipperShield_phys", fWorldLV_rot, false, -1, fCheckOverlaps);

        G4LogicalVolume *CoipperShield_log5 =
            new G4LogicalVolume(CoipperShield_box5, CopperMaterial, "CoipperShield_log", 0, 0, 0);

        vol_phys =
            new G4PVPlacement(0, G4ThreeVector(3.15 * cm + 0.2 * mm, 0.0 * cm, -(52.8 / 2 * cm + 1.5 * d_y1) / 2.),
                              CoipperShield_log5, "CoipperShield_phys", fWorldLV_rot, false, -1, fCheckOverlaps);
        vol_phys =
            new G4PVPlacement(0, G4ThreeVector(3.15 * cm + 0.2 * mm, 0.0 * cm, (52.8 / 2 * cm + 1.5 * d_y1) / 2.),
                              CoipperShield_log5, "CoipperShield_phys", fWorldLV_rot, false, -1, fCheckOverlaps);

        CoipperShield_log1->SetVisAttributes(Shield_VisAtt);
        CoipperShield_log2->SetVisAttributes(Shield_VisAtt);
        CoipperShield_log3->SetVisAttributes(Shield_VisAtt);
        CoipperShield_log4->SetVisAttributes(Shield_VisAtt);
        CoipperShield_log5->SetVisAttributes(Shield_VisAtt);


#endif // SHIELD

        ////////////////////////////////////////////////////////////////////////////////
        // TITAN FOILS
        ////////////////////////////////////////////////////////////////////////////////

#if defined(FOIL)

        // титановая фольга 70 мкм закрывает вылет из окна
#if defined(RUN21_for_TARGET)

        G4Box *TitanFoil_box = new G4Box("TitanFoil_box", 5.2 / 2 * cm, 0.007 / 2 * cm, (w_lng + 1.4 * cm) / 2.);

        G4LogicalVolume *TitanFoil_log = new G4LogicalVolume(TitanFoil_box, TitanMaterial, "TitanFoil_log", 0, 0, 0);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0 * cm, 2.81 * cm, w_pos), TitanFoil_log, "TitanFoil_phys",
                                     fWorldLV, false, -1, fCheckOverlaps);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0 * cm, -2.81 * cm, w_pos), TitanFoil_log, "TitanFoil_phys",
                                     fWorldLV, false, -1, fCheckOverlaps);


        TitanFoil_log->SetVisAttributes(TitanFoil_VisAtt);
#endif

        // бериллиевая фольга 300 мкм закрывает вылет из окна

#if defined(RUN23_for_TARGET)

        G4Box *BerylliumFoil_box = new G4Box("BerylliumFoil_box", 5.2 / 2 * cm, 300. / 2 * um, (w_lng + 1.4 * cm) / 2.);

        G4LogicalVolume *BerylliumFoil_log =
            new G4LogicalVolume(BerylliumFoil_box, BerylliumMaterial, "BerylliumFoil_log", 0, 0, 0);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0 * cm, 2.82 * cm, w_pos), BerylliumFoil_log,
                                     "BerylliumFoil_phys", fWorldLV_rot, false, -1, fCheckOverlaps);
        vol_phys = new G4PVPlacement(0, G4ThreeVector(0.0 * cm, -2.82 * cm, w_pos), BerylliumFoil_log,
                                     "BerylliumFoil_phys", fWorldLV_rot, false, -1, fCheckOverlaps);


        BerylliumFoil_log->SetVisAttributes(BerylliumFoil_VisAtt);
#endif

#endif // FOIL


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
        G4double x_pos_mag = 8.5 * cm; // was 10.5 cm
        G4double d_y1_mag = 1. * cm;
        G4double y_pos_mag = (d_y1_mag + y_size_mag / 2.) / 2.;
        G4double z_pos_mag = 0. * cm; // was 5 cm


        G4Box *solidMag1 = new G4Box("Mag1", x_size_mag / 2., y_size_mag / 2., z_size_mag / 2.);
        G4LogicalVolume *logicMag1 = new G4LogicalVolume(solidMag1, CopperMaterial, "Mag1");
        new G4PVPlacement(0, G4ThreeVector(x_pos_mag, 0., z_pos_mag), logicMag1, "Mag1", fWorldLV_rot, false, -1,
                          fCheckOverlaps);
        new G4PVPlacement(0, G4ThreeVector(-x_pos_mag, 0., z_pos_mag), logicMag1, "Mag1", fWorldLV_rot, false, -1,
                          fCheckOverlaps);

        G4Box *solidMag1a =
            new G4Box("Mag1a", x_size_mag / 2. + d_z_mag, (y_size_mag / 2. - d_y1_mag) / 2., z_size_mag / 2. + d_z_mag);
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


    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void PrimaryGeneratorAction::GenerateElectronGEM(G4Event *event)
    {

        G4double xx_cell = 0.0, yy_cell = 0.0, zz_cell = 0.0;

        G4ThreeVector initial_pos(0.0, 0.0, 0.0);
        G4ThreeVector initial_momentum_direction(0.0, 0.0, 0.0);

        // initial_pos = G4ThreeVector(0.0, 20.0 * cm, 0.0);
        initial_momentum_direction = G4ThreeVector(0.0, 0.0, 1.0);


        G4double alpha_GEM = alpha0_GEM;
        G4double theta_GEM = 0.0;


        //

        G4double kinetic_energy = 100. * MeV;


        G4double Egamma = 400. * MeV;


        // Egamma = virtual_photon_random(360., 650.) * MeV;
        Egamma = virtual_photon_random(400., 640.) * MeV;
        //  Egamma = virtual_photon_random(150., 780.) * MeV;

        //  Egamma = 0.;
        //  Egamma = virtual_photon_random(350., 680.) * MeV;

        kinetic_energy = Ebeam * MeV - Egamma;

        if (1)
        {
            if (G4UniformRand() > 0.5)
            {
                //  kinetic_energy = (G4UniformRand() > 0.5) ? 160.0 * MeV : 400.0 * MeV;
                kinetic_energy = Ebeam * MeV - Egamma;
            }
            else
            {
                kinetic_energy = Ebeam * MeV;
            }
        }

        //   kinetic_energy = Ebeam * MeV;

        // kinetic_energy = 350.0 * MeV;
        //  kinetic_energy = 200.0 * MeV;

        // небольшая размазка, которая зависит от энергии пучка
        // alpha_GEM = G4RandGauss::shoot(alpha0_GEM, 2.0* (1. - kinetic_energy/ (Ebeam * MeV)) * mrad);
        // если E=Ebeam, alpha_GEM=alpha0_GEM
        // если E=0.5*Ebeam, alpha_GEM=alpha0_GEM +- 1.0 mrad
        // если E=0.2*Ebeam, alpha_GEM=alpha0_GEM +- 0.8 mrad

        // или по формуле для тормозного излучения
        // G4double sigma_alpha_GEM = ((Me * MeV) / kinetic_energy) * rad;
        //  if (kinetic_energy == Ebeam * MeV) sigma_alpha_GEM = 0.;
        //  alpha_GEM = G4RandGauss::shoot(alpha0_GEM,  sigma_alpha_GEM);

        // генерируем вектор тормозного фотона

        G4double energy_Br = 0.0;
        G4double theta_BR = 0.0;
        G4double phi_BR = 0.0;

        G4double e_min_cut = 10 * keV;
        // e_min_cut = 150.0 * MeV;
        G4double e_max_cut = 100 * keV;
        // e_max_cut = -1; // без ката
        // e_max_cut = 700.0 * MeV;
        e_max_cut = 0.2 * kinetic_energy;


        //  BetheHeitlerGenerator event_br;

        event_br->sampleBremsstrahlungGenerator(kinetic_energy, e_min_cut, e_max_cut);
        // BremsstrahlungGenerator(kinetic_energy,e_min_cut ,e_max_cut, energy_Br, theta_BR, phi_BR);
        energy_Br = event_br->photon_energy;
        theta_BR = event_br->photon_theta;
        phi_BR = event_br->photon_phi;


        G4ThreeVector dir_BR(0.0, 0.0, 1.0);
        dir_BR.setTheta(theta_BR);
        // dir_BR.setTheta(alpha0_GEM+theta_BR);
        dir_BR.setPhi(phi_BR);
        // dir_BR = G4ThreeVector(sin(theta_BR)*cos(phi_BR), sin(theta_BR)*sin(phi_BR), cos(alpha_GEM));

        // для генерации за магнитом D2
        alpha_GEM = alpha0_GEM;
        // theta_GEM = asin(((Ebeam*MeV /kinetic_energy))*(sin(alpha0_GEM)+sin(theta0_GEM)) - sin(alpha_GEM)) * rad;

        initial_momentum_direction = G4ThreeVector(-sin(alpha_GEM), 0.0, cos(alpha_GEM));
        // initial_momentum_direction +=dir_BR;  // размазка тормозным излучением
        // initial_momentum_direction.rotate(initial_momentum_direction, )
        auto initial_momentum_direction_old = G4ThreeVector(initial_momentum_direction);
        initial_momentum_direction.rotateY(theta_BR);
        initial_momentum_direction.rotate(phi_BR, initial_momentum_direction_old);
        alpha_GEM = initial_momentum_direction.getTheta();

        theta_GEM = asin(((Ebeam * MeV / kinetic_energy)) * (sin(alpha0_GEM) + sin(theta0_GEM)) - sin(alpha_GEM)) * rad;
        initial_momentum_direction = G4ThreeVector(sin(theta_GEM), 0.0, cos(theta_GEM));


        G4double R_GEM = D_GEM / (sin(alpha_GEM) + sin(theta_GEM));
        G4double R0_GEM = D_GEM / (sin(alpha0_GEM) + sin(theta0_GEM));

        G4double X0_GEM =
            R0_GEM * (cos(alpha0_GEM) - cos(theta0_GEM)); // X на вызоде из магнита D2 для невозмущенной траектории
        G4double X_GEM = R_GEM * (cos(alpha_GEM) - cos(theta_GEM)); // X на вызоде из магнита D2

        //  X_GEM = R_GEM - R_GEM * cos(asin(D_GEM/R_GEM));


        initial_pos = G4ThreeVector((1.0 * mm + (X0_GEM + X_GEM)), 0.0,
                                    L0_GEM + D_GEM + 0 * (Z1_GEM + L_GEM + 20. / 2. * mm + 2. * mm));


        // генерация из мишени для расчета в поле D2, перед магнитом D2

        //  xx_cell = G4RandGauss::shoot(meanX_beam, Xsigma_beam) * mm; // в см
        // yy_cell = G4RandGauss::shoot(meanY_beam, Ysigma_beam) * mm;

        //  G4double r2 = G4UniformRand();
        // zz_cell = -cell_z_size_2 + r2 * l_zz_cell; // равномерно

        // kinetic_energy = 800.0 * MeV;
        // kinetic_energy = 300.0 * MeV;
        // initial_pos = G4ThreeVector(0.0, 0.0, 0.0);

        G4double Xbeam = 0., Ybeam = 0.;
        // initial_pos = G4ThreeVector(0.0, 0.0, 0. * cm);
        //  initial_pos = G4ThreeVector(xx_cell, yy_cell, zz_cell);
        auto initial_pos_in_cell = GenVertex(Xbeam, Ybeam, Xsigma_beam, Ysigma_beam);
        initial_pos_in_cell.rotateY(-alpha0_GEM); // поворот ячейки

        initial_pos = initial_pos_in_cell;


        alpha_GEM = alpha0_GEM;
        initial_momentum_direction = G4ThreeVector(-sin(alpha_GEM), 0.0, cos(alpha_GEM));

        G4LorentzVector Electron_4vec_initial(0, 0, 0, 0); // электрон до тормозного изулучения
        G4ThreeVector Electron_3vec_initial(0.0, 0.0, 1.0); // импульс

        G4LorentzVector Electron_4vec_final(0, 0, 0, 0); // электрон после тормозного изулучения
        G4ThreeVector Electron_3vec_final(0.0, 0.0, 1.0); // импульс

        G4LorentzVector PhotonBr_4vec(0, 0, 0, 0); // фотон тормозного излучения
        G4ThreeVector PhotonBr_3vec(0.0, 0.0, 1.0); // импульс


        // электрон до тормозного излучения
        G4double Pe_initial = sqrt(kinetic_energy * (kinetic_energy + 2. * Me * MeV)); // модуль импульс электрона
        Electron_3vec_initial.setTheta(initial_momentum_direction.getTheta());
        Electron_3vec_initial.setPhi(initial_momentum_direction.getPhi());
        Electron_3vec_initial.setMag(Pe_initial);
        // Electron_4vec_initial.setVect(Pe_vec_initial);
        // Electron_4vec_initial.setE(Eelectron_initial + Me); // полная энергия
        Electron_4vec_initial.setVectMag(Electron_3vec_initial, Me * MeV);

        // фотон тормозного излучения
        G4double Pgamma_br = energy_Br;

        // PhotonBr_3vec.setTheta(theta_BR);
        // PhotonBr_3vec.setPhi(phi_BR);
        PhotonBr_3vec = dir_BR;
        PhotonBr_3vec.setMag(Pgamma_br);
        // PhotonBr_3vec = G4ThreeVector(0.,0.,0.); // зануление
        PhotonBr_4vec.setVectMag(PhotonBr_3vec, 0.);

        // PhotonBr_4vec.setVectMag(G4ThreeVector(0.,0.,0.),0.); // зануление

        // поворачиваем угол фотона по ходу движения начального электрона

        // G4ThreeVector v1(0.0, 0.0, 0.0); // исходный вектор
        // G4ThreeVector v2(0.0, 0.0, 0.0); // целевой вектор
        //
        // v1 = G4ThreeVector(dir_BR);
        // v2 = G4ThreeVector(initial_momentum_direction);
        //
        // // 1. Нормализуем оба вектора
        // G4ThreeVector u1 = v1.unit();
        // G4ThreeVector u2 = v2.unit();
        //
        // // 2. Находим ось вращения через векторное произведение (cross product)
        // G4ThreeVector axis = u1.cross(u2);
        // G4double axisMag = axis.mag();
        //
        // G4RotationMatrix rotMat;
        //
        // if (axisMag > 1e-6) {
        //     axis /= axisMag; // нормализуем ось
        //     // 3. Угол между векторами через скалярное произведение (dot product)
        //     G4double cosAngle = u1.dot(u2);
        //     // Защита от выхода за пределы [-1, 1] из-за погрешностей округления
        //     if (cosAngle > 1.0) cosAngle = 1.0;
        //     if (cosAngle < -1.0) cosAngle = -1.0;
        //     G4double angle = std::acos(cosAngle);
        //
        //     // 4. Заполняем матрицу поворота
        //     rotMat.rotate(angle, axis);
        // } else if (u1.dot(u2) < 0) {
        //     // Если векторы направлены строго в противоположные стороны,
        //     // выбираем любую перпендикулярную ось
        //     G4ThreeVector perpAxis = u1.orthogonal();
        //     rotMat.rotate(M_PI, perpAxis);
        // }
        //
        // // Теперь v1, повернутый с помощью rotMat, будет сонаправлен с v2:
        // G4ThreeVector rotatedV1 = rotMat * v1;

        // обратно перенормируем
        // rotatedV1.setMag(PhotonBr_3vec.mag());

        // PhotonBr_3vec = G4ThreeVector(rotatedV1);
        // PhotonBr_3vec.setMag(Pgamma_br);

        // PhotonBr_4vec.setVect(Pe_vec_initial);
        // PhotonBr_4vec.setE(Eelectron_initial + Me); // полная энергия
        // PhotonBr_4vec.setVectMag(PhotonBr_3vec, 0.0);

        // поворачиваем угол фотона по ходу движения начального электрона

        // PhotonBr_4vec.rotateUz(initial_momentum_direction.unit());
        PhotonBr_4vec = G4LorentzVector(PhotonBr_4vec.rotateUz(initial_momentum_direction.unit()));

        initial_momentum_direction_old = G4ThreeVector(initial_momentum_direction);
        Electron_4vec_final = Electron_4vec_initial - PhotonBr_4vec; // вычетание четырех векторов
        Electron_3vec_final = Electron_4vec_final.getV();
        initial_momentum_direction = Electron_3vec_final.unit(); // нормализуем
        G4double kinetic_energy_initial = kinetic_energy;
        kinetic_energy = Electron_4vec_final.getT() - Me * MeV;
        // kinetic_energy = sqrt(pow(Electron_4vec_final.getV().mag(),2.) + pow(Me * MeV,2.))  - Me * MeV;

        if (kinetic_energy_initial == (Ebeam * MeV)) // для энергии пучка не делаем размазку тормозным излучением
        {
            initial_momentum_direction = initial_momentum_direction_old;
            kinetic_energy = kinetic_energy_initial;
        }


        if (true or (kinetic_energy_initial != (Ebeam * MeV)))
        {
            G4cerr << "initial_kinetic_energy = " << kinetic_energy_initial
                   << " final_kinetic_energy = " << kinetic_energy
                   << " initial_dir = " << initial_momentum_direction_old.getTheta() / mrad << " "
                   << initial_momentum_direction_old.getPhi() / mrad
                   << " final_dir = " << initial_momentum_direction.getTheta() / mrad << " "
                   << initial_momentum_direction.getPhi() / mrad << " BR_dir = " << dir_BR.getTheta() / mrad << " "
                   << dir_BR.getPhi() / mrad << " BR_dir_rot = " << PhotonBr_4vec.getV().getTheta() / mrad << " "
                   << PhotonBr_4vec.getV().getPhi() / mrad << G4endl

                   << G4endl;
        }

        // G4cerr << "initial_kinetic_energy = " << kinetic_energy_initial / MeV << " final_kinetic_energy = " <<
        // kinetic_energy / MeV << " " <<  (Electron_4vec_final.getT() - Me * MeV) / MeV<< G4endl;

        // initial_momentum_direction +=dir_BR;  // размазка тормозным излучением
        // initial_momentum_direction.rotateY(theta_BR);
        initial_momentum_direction_old = G4ThreeVector(initial_momentum_direction);
        // initial_momentum_direction.rotateY(theta_BR);
        //   initial_momentum_direction.rotate(phi_BR, initial_momentum_direction_old);
        // initial_momentum_direction.rotateUz(dir_BR);
        alpha_GEM = initial_momentum_direction.getTheta();

        //  G4cerr << "!!! = " << kinetic_energy / MeV << " " << 1000. * theta_BR  << " " << 1000.*(dir_BR.getTheta())
        //  << " "  << 1000.* alpha0_GEM << " "  <<  1000.*alpha_GEM<< G4endl;
        // G4cerr << "!!! = "  << 1000.*theta0_BR << " " << 1000.* sigma_theta_BR << " " << 1000.*theta_BR << " "
        // <<landau_value << G4endl;

        // для генерации за магнитом D1
        //  if (0)
        if (kinetic_energy == Ebeam * MeV)
        {
            //  kinetic_energy = 150 * MeV;
            kinetic_energy = 800 * MeV;
            initial_pos = G4ThreeVector(2.6 * cm + G4RandGauss::shoot(0., Xsigma_beam),
                                        G4RandGauss::shoot(0., Ysigma_beam), -80. * cm);
            // initial_pos = G4ThreeVector(2.5 * cm, 0.0 * cm, -80. * cm);
            initial_momentum_direction = G4ThreeVector(0.0, 0.0, 1.0);
            // initial_momentum_direction +=dir_BR;  // размазка тормозным излучением
            //  auto initial_momentum_direction_old =  G4ThreeVector(initial_momentum_direction);
            // initial_momentum_direction.rotateY(theta_BR);
            // initial_momentum_direction.rotate(phi_BR, initial_momentum_direction_old);
        }

        // для генерации за магнитом D3
        if (0)
        {
            // kinetic_energy = 150 * MeV;
            kinetic_energy = 800 * MeV;
            // initial_pos = G4ThreeVector(1.6 * cm + G4RandGauss::shoot(0.,Xsigma_beam ),
            // G4RandGauss::shoot(0.,Ysigma_beam ), 120. * cm);
            initial_pos = G4ThreeVector(1.6 * cm, 0.0 * cm, 120. * cm);

            // initial_momentum_direction = G4ThreeVector(0.0, 0.0, 1.0);
            //   G4double sigma_alpha_GEM = ((Me * MeV) / kinetic_energy) * rad;
            //  if (kinetic_energy == Ebeam * MeV) sigma_alpha_GEM = 0.;
            //  alpha_GEM = G4RandGauss::shoot(alpha0_GEM,  sigma_alpha_GEM);
            alpha_GEM = alpha0_GEM;
            theta_GEM =
                asin(((Ebeam * MeV / kinetic_energy)) * (sin(alpha0_GEM) + sin(theta0_GEM)) - sin(alpha_GEM)) * rad;

            initial_momentum_direction = G4ThreeVector(-sin(alpha_GEM), 0.0, cos(alpha_GEM));
            // initial_momentum_direction +=dir_BR;  // размазка тормозным излучение
            auto initial_momentum_direction_old = G4ThreeVector(initial_momentum_direction);
            initial_momentum_direction.rotateY(theta_BR);
            initial_momentum_direction.rotate(phi_BR, initial_momentum_direction_old);
            alpha_GEM = initial_momentum_direction.getTheta();
            theta_GEM =
                asin(((Ebeam * MeV / kinetic_energy)) * (sin(alpha0_GEM) + sin(theta0_GEM)) - sin(alpha_GEM)) * rad;
            initial_momentum_direction = G4ThreeVector(-sin(theta_GEM), 0.0, cos(theta_GEM));
        }
        //  G4cerr << "!!! = " << (X0_GEM + X_GEM ) / mm << G4endl;

        // G4cerr << "!!! = " << kinetic_energy / MeV << " " << X_GEM / mm << " " << theta_GEM << G4endl;

        fParticleGun->SetParticleEnergy(kinetic_energy);

        G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
        G4ParticleDefinition *particle = particleTable->FindParticle("e-");
        fParticleGun->SetParticleDefinition(particle);

        fParticleGun->SetParticleMomentumDirection(initial_momentum_direction);
        fParticleGun->SetParticlePosition(initial_pos);

        fParticleGun->GeneratePrimaryVertex(event);

        // отдельно гененируем тормозной фотон


        G4ParticleDefinition *gamma_particle = particleTable->FindParticle("gamma");
        fParticleGun->SetParticleDefinition(gamma_particle);
        fParticleGun->SetParticleEnergy(energy_Br);
        fParticleGun->SetParticlePosition(initial_pos);
        fParticleGun->SetParticleMomentumDirection((PhotonBr_4vec.getV()).unit());

        if (kinetic_energy_initial == (Ebeam * MeV)) // рисуем тормозные фотоны только от пучка
        {
            fParticleGun->SetParticlePosition(initial_pos_in_cell);
            auto rot_vec = G4ThreeVector(0.0, 0.0, 1.0);
            rot_vec.rotateY(-alpha0_GEM);
            //  rot_vec.rotateY(-5 * degree);
            auto PhotonBr_4vec_from_Ebeam = G4LorentzVector();
            PhotonBr_4vec_from_Ebeam.setVectMag(PhotonBr_3vec, 0.);
            PhotonBr_4vec_from_Ebeam.rotateUz(rot_vec.unit());
            fParticleGun->SetParticleMomentumDirection((PhotonBr_4vec_from_Ebeam.getV()).unit());
            //  fParticleGun->SetParticleMomentumDirection(rot_vec);
            fParticleGun->GeneratePrimaryVertex(event);
        }
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


    // void PrimaryGeneratorAction::BremsstrahlungGenerator(G4double electron_energy0, G4double e_min_cut, G4double
    // e_max_cut, G4double &photon_energy, G4double &photon_theta, G4double &photon_phi)
    // void PrimaryGeneratorAction::BremsstrahlungGenerator()
    void BetheHeitlerGenerator::sampleBremsstrahlungGenerator(G4double electron_energy0, G4double e_min_cut,
                                                              G4double e_max_cut)

    {

        // формула Бете — Гайтлера для тормозного излучения релятивистких эектронов

        // electron_energy Полная энергия налетающего электрона (МэВ)
        // e_cut           Инфракрасный порог обрезания энергии фотона

        if (electron_energy0 <= ((Me * MeV) + e_min_cut))
        {
            G4cerr << "Energy Electron is small" << G4endl;
        }

        // Кинетическая энергия и границы доли энергии y = E_gamma / E_e
        G4double E0_norm = electron_energy0 / (Me * MeV);
        y_min = e_min_cut / electron_energy0;
        y_max = e_max_cut / electron_energy0;

        if (e_max_cut == -1 or y_max >= (1.0 - 1.0 / E0_norm))
        {
            //  y_max = 1.0 - ((Me * MeV) / electron_energy);
            y_max = 1.0 - 1.0 / E0_norm;
            //  y_max = 1.0 - y_min;
        }


        G4double u_min = TMath::Log(y_min);
        G4double u_max = TMath::Log(y_max);

        //  G4cerr << "u_min = " << u_min << " u_max" << u_max << G4endl;

        // Оценка мажоранты (максимума функции g(y) на интервале [y_min, y_max])
        // Максимум достигается на левой границе y = y_min
        G4double f_max = evaluate_pdf_br(y_min, E0_norm) * 1.05;
        // f_max = 4.0 / 3.0;


        // G4cerr <<"y_min ... y_max = " << y_min << "..." << y_max << G4endl;

        if (y_min >= y_max)
        {
            G4cerr << "Incorrect e_cut" << G4endl;
        }

        // 1. Розыгрыш энергии фотона (метод отбора с мажорантой 1/y)
        G4double y = 0.0;
        G4double u = 0.0;
        G4double ln_ratio = TMath::Log(y_max / y_min);

        while (true)
        {
            G4double r1 = G4UniformRand();
            G4double r2 = G4UniformRand();

            // Генерация y согласно плотности p(y) ~ 1/y
            // y = y_min * TMath::Exp(r1 * ln_ratio);

            u = u_min + r1 * (u_max - u_min);
            y = TMath::Exp(u);

            // Функция отклонения g(y) = 1 - y + 0.75 * y^2 (максимум <= 1 на)
            G4double g_y = 1.0 - y + 0.75 * y * y;

            // if (r2 <= g_y) {
            //     break;
            // }

            if (r2 * f_max <= evaluate_pdf_br(y, E0_norm))
            {
                break;
            }

            // G4cerr << "!!!" << G4endl;
        }

        photon_energy = y * electron_energy0;
        Eelectron_final = electron_energy0 - photon_energy;

        // G4cerr <<"photon_energy (MeV) = " << photon_energy / MeV << " electron_energy (MeV) = " << electron_energy /
        // MeV << G4endl;

        // 2. Розыгрыш угла вылета фотона (приближение Цая)
        // Лоренц-фактор электрона
        G4double gamma = electron_energy0 / (Me * MeV);
        gamma = E0_norm;

        // Обратная функция распределения для f(u) ~ u / (1 + u^2)^2:
        // CDF(u) = u^2 / (1 + u^2) = r3  =>  u = sqrt(r3 / (1 - r3))
        G4double r3 = G4UniformRand();
        G4double r4 = G4UniformRand();
        // Защита от деления на ноль при r3 -> 1.0
        if (r3 > 0.999999999)
            r3 = 0.999999999;

        // Упрощенное угловое распределение dN/d(theta) ~ theta / (theta^2 + 1/gamma^2)^2
        G4double u2 = TMath::Sqrt(r3 / (1.0 - r3));
        photon_theta = (1.0 / gamma) * u2;

        // Ограничение угла физическим пределом [0, PI]
        if (photon_theta > M_PI)
            photon_theta = M_PI;

        // Азимутальный угол изотропен
        photon_phi = 2.0 * M_PI * r4;

        photon_theta *= rad;
        photon_phi *= rad;

        /// мой изначальный способ

        // G4double theta0_BR = ((Me * MeV) / kinetic_energy) * rad; // пик распределения
        // theta0_BR /= sqrt(7.0); // пишут, что нужно делить на корень из 7
        // // G4double sigma_theta_BR = TMath::Sqrt(theta0_BR / rad) * rad; // ширину оцениваем как корень
        // //G4double sigma_theta_BR = 1.0 * mrad;
        // G4double sigma_theta_BR = 0.1 *theta0_BR;
        // // Стандартное распределение Ландау в CLHEP генерируется вокруг наиболее вероятного значения (MPV)
        // // ≈ -0.222782 и со стандартным параметром масштаба (ширины) ξ = 1.
        // G4double landau_value = CLHEP::RandLandau::shoot(G4Random::getTheEngine()) + 0.222782;
        // landau_value = fabs(landau_value);
        // //G4double landau_value = G4RandLandau::shoot();
        // // G4double theta_BR = theta0_BR + sigma_theta_BR * landau_value;  // ассиметричное распределение, подгоняем
        // распределением Ландау G4double theta_BR = fabs(G4RandGauss::shoot(theta0_BR / rad,  sigma_theta_BR)) * rad;
        // // theta_BR = 0.;
        //
        // // theta_BR = -1* deg;
        // // G4cerr << "!!! = " << 1000.* (alpha_GEM - alpha0_GEM) << G4endl;
        // G4double phi_BR = 2.0*M_PI * G4UniformRand();

        // через встроенный в GEANT класс

        auto br_model = new G4eBremsstrahlungRelModel();
        // G4BetheHeitlerModel() G4BetheHeitler5DModel
        auto photon_dp = new G4DynamicParticle();
        auto photon_cc = new G4MaterialCutsCouple();
        //  br_model->SampleSecondaries(photon_dp, photon_cc);
    }


} // namespace Cosmic
