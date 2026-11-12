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
/// \file DetectorConstruction.hh
/// \brief Definition of the Cosmic::DetectorConstruction class

#ifndef CosmicDetectorConstruction_h
#define CosmicDetectorConstruction_h 1

#include "TaggingSystem/forTaggingSystem.hh"

#include "Constants.hh"
#include "HadronCalorimeterSD.hh"
#include "PlasticSD.hh"

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"


class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;
class G4UserLimits;

#pragma once

namespace Cosmic
{

/// Detector construction class to define materials and geometry.
/// The calorimeter is a box made of a given number of layers. A layer consists
/// of an absorber plate and of a detection gap. The layer is replicated.
///
/// Four parameters define the geometry of the calorimeter :
///
/// - the thickness of an absorber plate,
/// - the thickness of a gap,
/// - the number of layers,
/// - the transverse size of the calorimeter (the input face is a square).
///
/// In ConstructSDandField() sensitive detectors of CalorimeterSD type
/// are created and associated with the Absorber and Gap volumes.
/// In addition a transverse uniform magnetic field is defined
/// via G4GlobalMagFieldMessenger class.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction() override;

  public:
    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    G4int GetType(G4int n);
    G4int GetType2(G4int n);
    G4double GetAttenuL(G4int n);
    G4double GetDiscrThr(G4int n);
    G4double GetDiscrThr2(G4int n);


  private:
    G4LogicalVolume *fWorldLV;
    G4LogicalVolume *fWorldLV_rot = nullptr; // для повернутой мишени
    G4VPhysicalVolume *worldPV;


    // methods
    //
    void DefineMaterials();

    void DefineRotationMatrices();

    void DefineVisAttributes();

    G4VPhysicalVolume *DefineVolumes();

    G4LogicalVolume *ConstructHadronCalorimeter(G4int nsys);

    void ConstructMRPC_old();

    void ConstructMRPC();

    G4LogicalVolume *ConstructVC(G4LogicalVolume *&forSD);

    G4LogicalVolume *ConstructWC(G4double Lwin, G4double Wwin, G4int ind, G4LogicalVolume *&forSD);

    void ConstructPlasticFat1();

    void ConstructPlasticFat2();

    void ConstructPlasticThin1();

    void ConstructPlasticThin2();

    void ConstructTarget();

    void ConstructMagnetCell();

    void ConstructMagnetD2();

    void ConstructMagnetD1();
    void ConstructMagnetD3();

    void ConstructMagnetQF();
    void ConstructMagnetQD();

    void ConstructGEM();


    G4LogicalVolume * ConstructLOWQ(G4int nsys);

    G4Material *GasAr = nullptr;
    G4Material *liqAr = nullptr;
    G4Material *Tungsten = nullptr;
    G4Material *Alumin = nullptr;
    G4Material *Beryllium = nullptr;
    G4Material *Glass = nullptr;
    G4Material *DurAlumin = nullptr;
    G4Material *Titan = nullptr;
    G4Material *Copper = nullptr;
    G4Material *Air = nullptr;
    G4Material *Scintil = nullptr;
    G4Material *Scintil_simple_plastic = nullptr;
    G4MaterialPropertiesTable *MPT_BC_422 = nullptr;
    G4Material *Scintil_BC_422 = nullptr;
    G4MaterialPropertiesTable *MPT_BC_420 = nullptr;
    G4Material *Scintil_BC_420 = nullptr;
    G4Material *LightGuide_material = nullptr;
    G4Material *OpticalGlue_material = nullptr;
    G4Material *OpticalGrease_material = nullptr;
    G4Material *Silicon = nullptr;
    G4Material *Mylar = nullptr;
    G4Material *AlMylar = nullptr;
    G4Material *BlackPaper = nullptr;
    G4Material *aluminum_mirror = nullptr;
    G4Material *BGO = nullptr;
    G4Material *NaI = nullptr;
    G4Material *CsI = nullptr;
    G4Material *Brass = nullptr;
    G4Material *Steel = nullptr;
    G4Material *Iron = nullptr;
    G4Material *GasCO2 = nullptr;
    G4Material *GasArCO2 = nullptr;
    G4Material *MgO = nullptr;
    G4Material *SiO2 = nullptr;
    G4Material *Vacuum = nullptr;
    G4Material *Galactic = nullptr;
    G4Material *H2O = nullptr;
    G4Material *PbWO4 = nullptr;
    G4Material *Aerog = nullptr;


    //      G4RotationMatrix Rotate30X;
    G4RotationMatrix Rotate180X;
    G4RotationMatrix Rotate180Y;
    G4RotationMatrix Rotate180Z;
    G4RotationMatrix Rotate270X;
    G4RotationMatrix Rotate90X;
    G4RotationMatrix Rotate9X;
    G4RotationMatrix Rotate99X;
    G4RotationMatrix Rotate45X;
    G4RotationMatrix Rotate45X180Z;
    G4RotationMatrix Rotate10X;
    G4RotationMatrix Rotate25X;
    G4RotationMatrix Rotate65X;
    G4RotationMatrix Rotate90Y;
    G4RotationMatrix Rotate270Y;
    G4RotationMatrix RotateNull;
    G4RotationMatrix Rotate90X180Z;
    G4RotationMatrix Rotate90Y180Z;
    G4RotationMatrix Rotate90Z;
    G4RotationMatrix Rotate90Y180X;
    G4RotationMatrix Rotate270Z;
    G4RotationMatrix Rotate270Y180X;

// VIS Atribyte

        G4VisAttributes *Steel_VisAtt = nullptr;
    G4VisAttributes *Mag_VisAtt = nullptr;
    G4VisAttributes *Shield_VisAtt = nullptr;
    G4VisAttributes *Iron_VisAtt = nullptr;
    G4VisAttributes *Alum_VisAtt = nullptr;
    G4VisAttributes *Gas_VisAtt = nullptr;
    G4VisAttributes *Foil_VisAtt = nullptr;
    G4VisAttributes *Mylar_VisAtt = nullptr;
    G4VisAttributes *Stef_VisAtt = nullptr;
    G4VisAttributes *Plastic_VisAtt = nullptr;
    G4VisAttributes *Convertor_VisAtt = nullptr;
    G4VisAttributes *Concrete_VisAtt = nullptr;
    G4VisAttributes *WOOD_VisAtt = nullptr;
    G4VisAttributes *Convertor_LQ_VisAtt = nullptr;
    G4VisAttributes *ProCover_VisAtt = nullptr;
    G4VisAttributes *TitanFoil_VisAtt = nullptr;
    G4VisAttributes *BerylliumFoil_VisAtt = nullptr;


    // PlasticSD* aplasticSD;


        // data members
        //
        static G4ThreadLocal G4GlobalMagFieldMessenger *fMagFieldMessenger;
        // magnetic field messenger

        G4bool fCheckOverlaps = true; // option to activate checking of volumes overlaps

        // Логические Объемы for SD

        G4LogicalVolume *plastic_fat_nsys1LV = nullptr;
        G4LogicalVolume *plastic_fat_nsys2LV_120 = nullptr;
        G4LogicalVolume *plastic_fat_nsys2LV_125 = nullptr;
        G4LogicalVolume *plastic_thin_nsys1LV = nullptr;
        G4LogicalVolume *plastic_thin_nsys2LV = nullptr;

        G4LogicalVolume *box_GEM_log = nullptr;

        G4LogicalVolume *scint_HadCal_nsys1LV = nullptr;
        G4LogicalVolume *scint_HadCal_nsys2LV = nullptr;

    G4LogicalVolume* scint_LQ_nsys1LV_layer1 = nullptr;
    G4LogicalVolume* scint_LQ_nsys1LV_layer2 = nullptr;
    G4LogicalVolume* scint_LQ_nsys2LV = nullptr;


    G4LogicalVolume* WCTheta1_gas_nsys1LV = nullptr;
    G4LogicalVolume* WCPhi1_gas_nsys1LV = nullptr;
    G4LogicalVolume* WCTheta2_gas_nsys1LV = nullptr;
    G4LogicalVolume* VCGas_log_nsys1LV = nullptr;

    G4LogicalVolume* WCTheta1_gas_nsys2LV = nullptr;
    G4LogicalVolume *WCPhi1_gas_nsys2LV = nullptr;
        G4LogicalVolume *WCTheta2_gas_nsys2LV = nullptr;
        G4LogicalVolume *VCGas_log_nsys2LV = nullptr;


        // для магнита D2 GEM

        G4LogicalVolume *fLogicD2gap = nullptr;
        G4LogicalVolume *fLogicD2gap_add1 = nullptr; // дополнительный объем после магнита для вываливания поля
        G4LogicalVolume *fLogicD2gap_add2 = nullptr;

        G4LogicalVolume *fLogicD1gap = nullptr;
        G4LogicalVolume *fLogicD3gap = nullptr;

        G4LogicalVolume *fLogicQFgap = nullptr;
        G4LogicalVolume *fLogicQDgap = nullptr;

        G4double x_pos_qf = 0., y_pos_qf = 0., z_pos_qf = 0.;
        G4double x_pos_qd = 0., y_pos_qd = 0., z_pos_qd = 0.;
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

