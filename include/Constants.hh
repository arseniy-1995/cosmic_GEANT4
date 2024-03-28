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
#include "G4PhysicalConstants.hh"

namespace Cosmic {

//#define isGenCosmic

//#define isGenLQ
//#define isTrigLQ

#define isGenPNpair
#define isTrigPN

    //#define isTrig_accept_all

    constexpr G4int T2M_zanulenie[3] = {1, 1, 1}; // если 1, то зануления нет, 0 - зануление есть T20, T21, T22

#define SW(N)   " " << std::setw(N)

    //#define RUN21
#define RUN23
#define LOWQ_CONVERTOR

#ifdef isGenCosmic
#define CONCRETE // Бетон
#endif

#define PF1_FAT // Толстые пластики
#define PF2_FAT

#define PF1_THIN // Тонкие пластики
#define PF2_THIN

#define HADCAL1 // Адронный калориметр
#define HADCAL2

    //#define MRPC1 // MRPC-детекторы
    //#define MRPC2

#define DCARM1 // Дрейфовая камера
#define DCARM2
#define VCARM1 // Вершинная камера
#define VCARM2

#define TARGET // Мишень+фольга
#define MAGNET // Магнит мишени
#define CELL
#define SHIELD // Фольга (экран) вокруг кубика
#define FOIL

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

#define N_AC    1
#define TOF_IND (AC_IND+N_AC)        /* MRPC */
#define N_TOF    1

//#define VC_IND (HCZ_IND + N_HCZ)    /* Vertex Ch */
//#define WC1_IND (VC_IND + NVC_WRS)  /* Drift Chambers */
//#define WC2_IND (WC1_IND + NW1_WRS) /* Drift Chambers */
//#define WC3_IND (WC2_IND + NW2_WRS) /* Drift Chambers */

#define LQ_IND (WC3_IND + NW3_WRS)  /* LQ layers */
#define ARM1_IND 0
#define ARM2_IND (LQ_IND + N_LQ)


#define VC_IND 0    /* Vertex Ch */
#define WC1_IND 0 /* Drift Chambers */
#define WC2_IND (WC1_IND + NW1_WRS) /* Drift Chambers */
#define WC3_IND (WC2_IND + NW2_WRS) /* Drift Chambers */



#define NSENS ARM2_IND * 2

// detector element type
#define DT_WIRE 0
#define DT_SCINT 1
#define DT_EMCAL 2

#define DC_RES (0.3*mm)
#define VC_RES (0.6*mm)


    constexpr G4int fNofLayers_plastic_fat_nsys1 = 6; // число толстых пластиков
    constexpr G4int fNofLayers_plastic_fat_nsys2 = 8;

    constexpr G4int fNofLayers_plastic_thin_nsys1 = 2; // число тонких пластиков
    constexpr G4int fNofLayers_plastic_thin_nsys2 = 2;

#ifdef RUN21
    constexpr G4int fNofLayers_plastic_LQ_nsys1 = 2; // число пластиков LQ электронного плеча
    constexpr G4int fNofLayers_plastic_LQ_nsys2 = 2;
#endif // RUN21

#ifdef RUN23
    constexpr G4int fNofLayers_plastic_LQ_nsys1 = 1; // число пластиков LQ электронного плеча
    constexpr G4int fNofLayers_plastic_LQ_nsys2 = 1;
#endif // RUN23

    constexpr G4int fNofLayers_HadrtonCalorimeter_nsys1 =
        N_HCX + N_HCZ + 0; // число чувствительных элементов в калориметре
    constexpr G4int fNofLayers_HadrtonCalorimeter_nsys2 = N_HCX + N_HCZ + 0;

    constexpr G4int fNofLayers_W_Chamber_nsys1 = NW1_WRS + NW2_WRS + NW3_WRS + NVC_WRS + 100;
    constexpr G4int fNofLayers_W_Chamber_nsys2 = NW1_WRS + NW2_WRS + NW3_WRS + NVC_WRS + 100;
    constexpr G4int fNofLayers_V_Chamber_nsys1 = NVC_WRS + 100;
    constexpr G4int fNofLayers_V_Chamber_nsys2 = NVC_WRS + 100;

    // пороги записи в файл
    constexpr G4double plastic_fat_threshold = 0.5 * MeV;
    constexpr G4double plastic_thin_threshold = 0.1 * MeV;
    constexpr G4double plastic_LQ_threshold = 0.1 * MeV;
    constexpr G4double HadronCalorimeter_threshold = 1.0 * MeV;

    //   constexpr G4double WChamber_threshold = 0.5 * keV;
    //   constexpr G4double VChamber_threshold = 0.5 * keV;

    constexpr G4double WChamber_threshold = 0.01 * keV;
    constexpr G4double VChamber_threshold = 0.01 * keV;

    constexpr G4double DENSITY_LO = 1.032;    // scintill. density - for light output calculation


    // FOR GENBOS

#define GENBOS

#define PHIC    90.*0.0174532925199
#define DPHI    35.*0.0174532925199
#define THDN    40.*0.0174532925199
#define XY_RAND
#define VMAX  8
    constexpr G4double cell_z_size_2 = 20.0 * cm;
    constexpr G4double meanX_beam = 0., meanY_beam = 0.;
    //constexpr G4double Xsigma_beam = 0.7 * mm, Ysigma_beam = 0.3 * mm; //for 2GeV
    constexpr G4double Xsigma_beam = 0.47 * mm, Ysigma_beam = 0.18 * mm; // for 800Mev из базы данных захода

// FOR LQ-GENERATOR

#define Md 1875.63// ????? ????????
#define Mp 938.28 // ????? ???????
#define Mn 939.52 // ????? ????????
#define Me 0.511 // ????? ????????

#define Mup 2.7928 // ????????? ?????? ???????
#define Mun -1.9130 // ????????? ?????? ????????
#define Mud 0.8574 // ????????? ?????? ????????

#define Degree .01745329252      /* degree->radians*/
#define Fm 197.32705359 // hc = ???*?? ?????????? ??????????
#define alpha_em 1./137.035989561 // ?????????? ?????? ?????????
#define charge_e 1.6021773349 // ????? ?????????, ????????? e-19 ??????
    constexpr G4double l_zz_cell = 40.0; // длина (полная) накопительной ячейки, см
    constexpr G4double Ebeam = 800.0;  // Энергия электронов, МэВ

#ifdef RUN21
    constexpr G4double delta_add_counter = 5.0; //добавляем несколько см, чтобы захватить большую площадь счетчика
#endif

#ifdef RUN23
    constexpr G4double delta_add_counter = 10.0; //добавляем несколько см, чтобы захватить большую площадь счетчика
#endif

    constexpr G4double l_theta_counter = 12.5 + delta_add_counter; // половина размера счетчика в тета направлении, см 15
    constexpr G4double l_phi_counter = 25.0 + delta_add_counter; // половина размера счетчика в тета направлении, см 25

#ifdef RUN21
    constexpr G4double R0_counter = 83.0; // радиус от центра мишени до середины счетчика, см
    constexpr G4double theta0_counter = 30.0 * Degree; // угол центра счетчика 30 градусов -> радианы
#endif
#ifdef RUN23
    constexpr G4double R0_counter = 80.0; // радиус от центра мишени до середины счетчика, см
    constexpr G4double theta0_counter = 33.0 * Degree; // угол центра счетчика 30 градусов -> радианы
#endif


    // constexpr G4double Pzz1 = 1.0, r = -2.0;
    // constexpr G4double Pzz2 = r * Pzz1;


        constexpr G4double Pzz1 = 0.39, r = -1.7;
    constexpr G4double Pzz2 = r * Pzz1;
}

#endif
