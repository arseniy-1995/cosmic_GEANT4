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
/// \file Constants.hh
/// \brief Definition of B5 example constants.

#ifndef CosmicConstants_h
#define CosmicConstants_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"

namespace Cosmic {


#define PF2_FAT
#define PF1_FAT

#define PF1_THIN
#define PF2_THIN

#define HADCAL1 // Адронный калориметр
#define HADCAL2

#define DCARM1 // Дрейфовая камера
#define DCARM2
#define VCARM1 // Вершинная камера
#define VCARM2

#define TARGET
#define MAGNET

#define LOWQ1 // электронные плечи ЛО-КУ поляриметра
#define LOWQ2


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



    constexpr G4int fNofLayers_plastic_fat_nsys1 = 6; // число толстых пластиков
    constexpr G4int fNofLayers_plastic_fat_nsys2 = 8;

    constexpr G4int fNofLayers_plastic_thin_nsys1 = 2; // число тонких пластиков
    constexpr G4int fNofLayers_plastic_thin_nsys2 = 2;

    constexpr G4int fNofLayers_plastic_LQ_nsys1 = 2; // число пластиков LQ электронного плеча
    constexpr G4int fNofLayers_plastic_LQ_nsys2 = 2;

    constexpr G4int fNofLayers_HadrtonCalorimeter_nsys1 = N_HCX + N_HCZ +100; // число чувствительных элементов в калориметре
    constexpr G4int fNofLayers_HadrtonCalorimeter_nsys2 = N_HCX + N_HCZ + 100;

    constexpr G4double plastic_fat_threshold = 1.0 * MeV; // порог записи в файл
    constexpr G4double plastic_thin_threshold = 0.1 * MeV;
    constexpr G4double plastic_LQ_threshold = 0.1 * MeV;
    constexpr G4double HadronCalorimeter_threshold = 0.5 * MeV;

    constexpr G4double DENSITY_LO = 1.032;	// scintill. density - for light output calculation





}

#endif
