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


#include "Constants.hh"
#include "PlasticSD.hh"
#include "HadronCalorimeterSD.hh"

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"


class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;
class G4UserLimits;

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
    G4double GetAttenuL(G4int n);
    G4double GetDiscrThr(G4int n);


  private:

   G4LogicalVolume *worldLV;
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

    void ConstructMagnet();
    G4LogicalVolume * ConstructLOWQ(G4int nsys);

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

        G4VisAttributes *Steel_VisAtt;
        G4VisAttributes *Mag_VisAtt;
        G4VisAttributes *Shield_VisAtt;
        G4VisAttributes *Iron_VisAtt;
        G4VisAttributes *Alum_VisAtt;
        G4VisAttributes *Gas_VisAtt;
        G4VisAttributes *Foil_VisAtt;
        G4VisAttributes *Mylar_VisAtt;
        G4VisAttributes *Stef_VisAtt;
        G4VisAttributes *Plastic_VisAtt;
        G4VisAttributes *Convertor_VisAtt;
        G4VisAttributes *Concrete_VisAtt;
        G4VisAttributes *WOOD_VisAtt;
        G4VisAttributes *Convertor_LQ_VisAtt;
        G4VisAttributes *ProCover_VisAtt;
        G4VisAttributes *TitanFoil_VisAtt;
        G4VisAttributes *BerylliumFoil_VisAtt;


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


        G4LogicalVolume *scint_HadCal_nsys1LV = nullptr;
        G4LogicalVolume *scint_HadCal_nsys2LV = nullptr;

    G4LogicalVolume* scint_LQ_nsys1LV_layer1 = nullptr;
    G4LogicalVolume* scint_LQ_nsys1LV_layer2 = nullptr;
    G4LogicalVolume *scint_LQ_nsys2LV = nullptr;


        G4LogicalVolume *WCTheta1_gas_nsys1LV = nullptr ;
        G4LogicalVolume *WCPhi1_gas_nsys1LV = nullptr;
        G4LogicalVolume *WCTheta2_gas_nsys1LV = nullptr;
        G4LogicalVolume *VCGas_log_nsys1LV = nullptr;

        G4LogicalVolume *WCTheta1_gas_nsys2LV = nullptr ;
        G4LogicalVolume *WCPhi1_gas_nsys2LV = nullptr;
        G4LogicalVolume *WCTheta2_gas_nsys2LV = nullptr;
        G4LogicalVolume *VCGas_log_nsys2LV = nullptr;



};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

