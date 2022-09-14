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

#include "PlasticSD.hh"

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#define PF2_FAT
#define PF1_FAT

//#define PF1_THIN
//#define PF2_THIN

//#define HADCAL1 // Адронный калориметр
//#define HADCAL2

//#define DCARM1 // Дрейфовая камера
//#define DCARM2
//#define VCARM1 // Вершинная камера
//#define VCARM2

#define NX_BARS 22
#define NZ_BARS 22
#define N_UNITS (NX_BARS / 2) /* assuming square calorimetr */
#define N_LAYERS 10
#define NVC_WRS 42
#define NW1_WRS 18
#define NW2_WRS 18
#define NW3_WRS 27
#define NCX_BRS 7
#define N_LQ 2

#define AC_IND 0
#define N_AC 2
#define DE_IND (AC_IND + N_AC)
#define N_DE 8
#define HCX_IND (DE_IND + N_DE) /* X bars */
#define N_HCX NX_BARS *N_LAYERS
#define HCZ_IND (HCX_IND + N_HCX) /* Z bars */
#define N_HCZ NZ_BARS *N_LAYERS
#define VC_IND (HCZ_IND + N_HCZ)    /* Vertex Ch */
#define WC1_IND (VC_IND + NVC_WRS)  /* Drift Chambers */
#define WC2_IND (WC1_IND + NW1_WRS) /* Drift Chambers */
#define WC3_IND (WC2_IND + NW2_WRS) /* Drift Chambers */
#define LQ_IND (WC3_IND + NW3_WRS)  /* LQ layers */
#define ARM1_IND 0
#define ARM2_IND (LQ_IND + N_LQ)

#define NSENS ARM2_IND * 2

// detector element type
#define DT_WIRE 0
#define DT_SCINT 1
#define DT_EMCAL 2

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
    // methods
    //
    void DefineMaterials();
    void DefineRotationMatrices();
    void DefineVisAttributes();

    G4VPhysicalVolume* DefineVolumes();

    G4LogicalVolume *ConstructHadronCalorimeter();
    G4LogicalVolume *ConstructVC();
    G4LogicalVolume *ConstructWC(G4double Lwin, G4double Wwin, G4int ind, G4LogicalVolume*& forSD);


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
        G4VisAttributes *Iron_VisAtt;
        G4VisAttributes *Alum_VisAtt;
        G4VisAttributes *Gas_VisAtt;
        G4VisAttributes *Foil_VisAtt;
        G4VisAttributes *Mylar_VisAtt;
        G4VisAttributes *Stef_VisAtt;
        G4VisAttributes *Plastic_VisAtt;
        G4VisAttributes *Convertor_VisAtt;
        G4VisAttributes *CONCRETE_VisAtt;
        G4VisAttributes *ProCover_VisAtt;

       // PlasticSD* aplasticSD;



    // data members
    //
    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger;
                                      // magnetic field messenger

    G4bool fCheckOverlaps = true; // option to activate checking of volumes overlaps
   // G4int  fNofLayers = -1;     // number of layers

    G4int  fNofLayers_plastic_fat_nsys1 = -1; // число толстых пластиков
    G4int  fNofLayers_plastic_fat_nsys2 = -1;


        // for SD


        G4LogicalVolume* plastic_fat_nsys1LV = nullptr;
        G4LogicalVolume* plastic_fat_nsys2LV = nullptr;

        G4LogicalVolume *WCTheta1_gas;
        G4LogicalVolume *WCPhi1_gas;
        G4LogicalVolume *WCTheta2_gas;
        G4LogicalVolume *VCGas_log;




};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

