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
/// \file EventAction.hh
/// \brief Definition of the Cosmic::EventAction class

#ifndef CosmicEventAction_h
#define CosmicEventAction_h 1

#include "Constants.hh"
#include "G4UserEventAction.hh"

#include "PlasticHit.hh"
#include "HadronCalorimeterHit.hh"
#include "ChamberHit.hh"

#include "globals.hh"

#include <vector>
#include <array>

const G4int NofLayers_plastic_fat_nsys1 = 6;
const G4int NofLayers_plastic_fat_nsys2 = 8;
const G4int NofLayers_plastic_thin_nsys1 = 2;
const G4int NofLayers_plastic_thin_nsys2 = 2;
const G4int NofLayers_plastic_LQ_nsys1 = 2;
const G4int NofLayers_plastic_LQ_nsys2 = 2;

namespace Cosmic {

/// Event action class
///
/// In EndOfEventAction(), it prints the accumulated quantities of the energy
/// deposit and track lengths of charged particles in Absober and Gap layers
/// stored in the hits collections.

    class EventAction : public G4UserEventAction {
    public:
        EventAction();

        ~EventAction() override;

        void BeginOfEventAction(const G4Event *event) override;

        void EndOfEventAction(const G4Event *event) override;

        ///////// ТОЛСТЫЕ ПЛАСТИКИ
        std::vector<G4double> &GetPlasticFatEdep(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_fatEdep[1];
            return fPlastic_fatEdep[0];
        }
        std::vector<G4double> &GetPlasticFatLO(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_fatLO[1];
            return fPlastic_fatLO[0];
        }
        std::vector<G4double> &GetPlasticFatA1(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_fatA1[1];
            return fPlastic_fatA1[0];
        }
        std::vector<G4double> &GetPlasticFatA2(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_fatA2[1];
            return fPlastic_fatA2[0];
        }
        std::vector<G4double> &GetPlasticFatT1(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_fatT1[1];
            return fPlastic_fatT1[0];
        }
        std::vector<G4double> &GetPlasticFatT2(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_fatT2[1];
            return fPlastic_fatT2[0];
        }
        std::vector<G4double> &GetPlasticFatTrackLength(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_fatTrackLength[1];
            return fPlastic_fatTrackLength[1];
        }
        std::vector<G4double> &GetPlasticFatToF(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_fatToF[1];
            return fPlastic_fatToF[0];
        }
        std::vector<G4double> &GetPlasticFatPos(G4int nsys = 1, G4int index = 1) {
            if (nsys == 1) {
                if (index == 1) return fPlastic_fatXPos[0];
                if (index == 2) return fPlastic_fatYPos[0];
                if (index == 3) return fPlastic_fatZPos[0];
            }
            if (nsys == 2) {
                if (index == 1) return fPlastic_fatXPos[1];
                if (index == 2) return fPlastic_fatYPos[1];
                if (index == 3) return fPlastic_fatZPos[1];
            }
            return fPlastic_fatXPos[0];
        }

        /////////////// ТОНКИЕ ПЛАСТИКИ

        std::vector<G4double> &GetPlasticThinEdep(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_thinEdep[1];
            return fPlastic_thinEdep[0];
        }
        std::vector<G4double> &GetPlasticThinLO(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_thinLO[1];
            return fPlastic_thinLO[0];
        }
        std::vector<G4double> &GetPlasticThinA1(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_thinA1[1];
            return fPlastic_thinA1[0];
        }
        std::vector<G4double> &GetPlasticThinA2(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_thinA2[1];
            return fPlastic_thinA2[0];
        }
        std::vector<G4double> &GetPlasticThinT1(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_thinT1[1];
            return fPlastic_thinT1[0];
        }
        std::vector<G4double> &GetPlasticThinT2(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_thinT2[1];
            return fPlastic_thinT2[0];
        }
        std::vector<G4double> &GetPlasticThinTrackLength(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_thinTrackLength[1];
            return fPlastic_thinTrackLength[1];
        }
        std::vector<G4double> &GetPlasticThinToF(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_thinToF[1];
            return fPlastic_thinToF[0];
        }
        std::vector<G4double> &GetPlasticThinPos(G4int nsys = 1, G4int index = 1) {
            if (nsys == 1) {
                if (index == 1) return fPlastic_thinXPos[0];
                if (index == 2) return fPlastic_thinYPos[0];
                if (index == 3) return fPlastic_thinZPos[0];
            }
            if (nsys == 2) {
                if (index == 1) return fPlastic_thinXPos[1];
                if (index == 2) return fPlastic_thinYPos[1];
                if (index == 3) return fPlastic_thinZPos[1];
            }
            return fPlastic_thinXPos[0];
        }
   

        ///////////////

        std::vector<G4double> &GetPlasticLQEdep(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_LQ_Edep[1];
            return fPlastic_LQ_Edep[0];
        }

        std::vector<G4double> &GetPlasticLQLO(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_LQ_LO[1];
            return fPlastic_LQ_LO[0];
        }


        std::vector<G4double> &GetPlasticLQA1(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_LQ_A1[1];
            return fPlastic_LQ_A1[0];
        }

        std::vector<G4double> &GetPlasticLQA2(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_LQ_A2[1];
            return fPlastic_LQ_A2[0];
        }

        std::vector<G4double> &GetPlasticLQT1(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_LQ_T1[1];
            return fPlastic_LQ_T1[0];
        }

        std::vector<G4double> &GetPlasticLQT2(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_LQ_T2[1];
            return fPlastic_LQ_T2[0];
        }

        std::vector<G4double> &GetPlasticLQTrackLength(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_LQ_TrackLength[1];
            return fPlastic_LQ_TrackLength[0];
        }

        std::vector<G4double> &GetPlasticLQToF(G4int nsys = 1) {
            if (nsys == 2) return fPlastic_LQ_ToF[1];
            return fPlastic_LQ_ToF[0];
        }

        std::vector<G4double> &GetPlasticLQPos(G4int nsys = 1, G4int index = 1) {
            if (nsys == 1) {
                if (index == 1) return fPlastic_LQ_XPos[0];
                if (index == 2) return fPlastic_LQ_YPos[0];
                if (index == 3) return fPlastic_LQ_ZPos[0];
            }
            if (nsys == 2) {
                if (index == 1) return fPlastic_LQ_XPos[1];
                if (index == 2) return fPlastic_LQ_YPos[1];
                if (index == 3) return fPlastic_LQ_ZPos[1];
            }

            return fPlastic_LQ_XPos[0];
        }
        
        ////////////////////

        // Методы для Адронного Калориметра

        std::vector<G4int> &Get_HCX_N(G4int nsys = 1) {
            if (nsys == 2) return fHCX_N[1];
            return fHCX_N[0];
        }

        std::vector<G4int> &Get_HCX_AL(G4int nsys = 1) {
            if (nsys == 2) return fHCX_AL[1];
            return fHCX_AL[0];
        }

        std::vector<G4double> &Get_HCX_Edep(G4int nsys = 1) {
            if (nsys == 2) return fHCX_Edep[1];
            return fHCX_Edep[0];
        }

        std::vector<G4double> &Get_HCX_LO(G4int nsys = 1) {
            if (nsys == 2) return fHCX_LO[1];
            return fHCX_LO[0];
        }

        std::vector<G4double> &Get_HCX_A1(G4int nsys = 1) {
            if (nsys == 2) return fHCX_A1[1];
            return fHCX_A1[0];
        }

        std::vector<G4double> &Get_HCX_A2(G4int nsys = 1) {
            if (nsys == 2) return fHCX_A2[1];
            return fHCX_A2[0];
        }

        std::vector<G4double> &Get_HCX_T1(G4int nsys = 1) {
            if (nsys == 2) return fHCX_T1[1];
            return fHCX_T1[0];
        }

        std::vector<G4double> &Get_HCX_T2(G4int nsys = 1) {
            if (nsys == 2) return fHCX_T2[1];
            return fHCX_T2[0];
        }

        std::vector<G4double> &Get_HCX_TrackLength(G4int nsys = 1) {
            if (nsys == 2) return fHCX_TrackLength[1];
            return fHCX_TrackLength[0];
        }

        std::vector<G4double> &Get_HCX_ToF(G4int nsys = 1) {
            if (nsys == 2) return fHCX_ToF[1];
            return fHCX_ToF[0];
        }

        std::vector<G4double> &Get_HCX_Pos(G4int nsys = 1, G4int index = 1) {

            if (nsys == 1) {
                if (index == 1) return fHCX_XPos[0];
                if (index == 2) return fHCX_YPos[0];
                if (index == 3) return fHCX_ZPos[0];
            }
            if (nsys == 2) {
                if (index == 1) return fHCX_XPos[1];
                if (index == 2) return fHCX_YPos[1];
                if (index == 3) return fHCX_ZPos[1];
            }
            return fHCX_XPos[0];

        }

        ///////


        std::vector<G4int> &Get_HCZ_N(G4int nsys = 1) {
            if (nsys == 2) return fHCZ_N[1];
            return fHCZ_N[0];
        }

        std::vector<G4int> &Get_HCZ_AL(G4int nsys = 1) {
            if (nsys == 2) return fHCZ_AL[1];
            return fHCZ_AL[0];
        }

        std::vector<G4double> &Get_HCZ_Edep(G4int nsys = 1) {
            if (nsys == 2) return fHCZ_Edep[1];
            return fHCZ_Edep[0];
        }

        std::vector<G4double> &Get_HCZ_LO(G4int nsys = 1) {
            if (nsys == 2) return fHCZ_LO[1];
            return fHCZ_LO[0];
        }

        std::vector<G4double> &Get_HCZ_A1(G4int nsys = 1) {
            if (nsys == 2) return fHCZ_A1[1];
            return fHCZ_A1[0];
        }

        std::vector<G4double> &Get_HCZ_A2(G4int nsys = 1) {
            if (nsys == 2) return fHCZ_A2[1];
            return fHCZ_A2[0];
        }

        std::vector<G4double> &Get_HCZ_T1(G4int nsys = 1) {
            if (nsys == 2) return fHCZ_T1[1];
            return fHCZ_T1[0];
        }

        std::vector<G4double> &Get_HCZ_T2(G4int nsys = 1) {
            if (nsys == 2) return fHCZ_T2[1];
            return fHCZ_T2[0];
        }

        std::vector<G4double> &Get_HCZ_TrackLength(G4int nsys = 1) {
            if (nsys == 2) return fHCZ_TrackLength[1];
            return fHCZ_TrackLength[0];
        }

        std::vector<G4double> &Get_HCZ_ToF(G4int nsys = 1) {
            if (nsys == 2) return fHCZ_ToF[1];
            return fHCZ_ToF[0];
        }

        std::vector<G4double> &Get_HCZ_Pos(G4int nsys = 1, G4int index = 1) {

            if (nsys == 1) {
                if (index == 1) return fHCZ_XPos[0];
                if (index == 2) return fHCZ_YPos[0];
                if (index == 3) return fHCZ_ZPos[0];
            }
            if (nsys == 2) {
                if (index == 1) return fHCZ_XPos[1];
                if (index == 2) return fHCZ_YPos[1];
                if (index == 3) return fHCZ_ZPos[1];
            }
            return fHCZ_XPos[0];

        }


        ///////// Методы для камер
        std::vector<G4int> &Get_W_N(G4int nsys = 1, G4int type_chamber = 1) {

            if (type_chamber == 1) // a
            {
                if (nsys == 1) return fWa_N[0];
                if (nsys == 2) return fWa_N[1];
            }
            if (type_chamber == 2) // b
            {
                if (nsys == 1) return fWb_N[0];
                if (nsys == 2) return fWb_N[1];
            }
            if (type_chamber == 3) // c
            {
                if (nsys == 1) return fWc_N[0];
                if (nsys == 2) return fWc_N[1];
            }
            return fWa_N[0];
        }
        std::vector<G4double> &Get_W_Pos(G4int nsys = 1, G4int index = 1, G4int type_chamber = 1) {

            if (type_chamber == 1) // a
            {
                if (nsys == 1) {
                    if (index == 1) return fWa_XPos[0];
                    if (index == 2) return fWa_YPos[0];
                    if (index == 3) return fWa_ZPos[0];
                }
                if (nsys == 2) {
                    if (index == 1) return fWa_XPos[1];
                    if (index == 2) return fWa_YPos[1];
                    if (index == 3) return fWa_ZPos[1];
                }
            }
            if (type_chamber == 2) // b
            {
                if (nsys == 1) {
                    if (index == 1) return fWb_XPos[0];
                    if (index == 2) return fWb_YPos[0];
                    if (index == 3) return fWb_ZPos[0];
                }
                if (nsys == 2) {
                    if (index == 1) return fWb_XPos[1];
                    if (index == 2) return fWb_YPos[1];
                    if (index == 3) return fWb_ZPos[1];
                }
            }
            if (type_chamber == 3) // c
            {
                if (nsys == 1) {
                    if (index == 1) return fWc_XPos[0];
                    if (index == 2) return fWc_YPos[0];
                    if (index == 3) return fWc_ZPos[0];
                }
                if (nsys == 2) {
                    if (index == 1) return fWc_XPos[1];
                    if (index == 2) return fWc_YPos[1];
                    if (index == 3) return fWc_ZPos[1];
                }
            }
            return fWa_XPos[0];

        }

        std::vector<G4int> &Get_VC_N(G4int nsys = 1) {

                if (nsys == 2) return fVC_N[1];
            return fVC_N[0];
        }

        std::vector<G4double> &Get_VC_Pos(G4int nsys = 1, G4int index = 1) {
                if (nsys == 1) {
                    if (index == 1) return fVC_XPos[0];
                    if (index == 2) return fVC_YPos[0];
                    if (index == 3) return fVC_ZPos[0];
                    if (index == 4) return fVC_RPos[0];
                }
                if (nsys == 2) {
                    if (index == 1) return fVC_XPos[1];
                    if (index == 2) return fVC_YPos[1];
                    if (index == 3) return fVC_ZPos[1];
                    if (index == 4) return fVC_RPos[1];
                }
            return fVC_XPos[0];
        }


    private:

        G4double vertex_x, vertex_y, vertex_z;    // vertex position
        G4int number_vertex;
        G4int vertex_index;// number of vertexes (particles)
        G4double vertex_energy, vertex_momentum, vertex_mass;
        G4double vertex_theta, vertex_phi;


        // methods
        PlasticHitsCollection *GetHitsCollection(G4int hcID, const G4Event *event) const;

        void PrintEventStatistics(G4double aplasticEdep, G4double plasticTrackLength) const;

        // energy deposit in calorimeters cells

        // здесь нулевой инедекс [0]  - первая система, первый индекс [1] - вторая система

        // Раньше было так
        //     std::vector<G4double> fPlastic_fat_nsys1Edep{std::vector<G4double>(fNofLayers_plastic_thin_nsys1 + 1, 0.0)};
        //        std::vector<G4double> fPlastic_fat_nsys2Edep{std::vector<G4double>(fNofLayers_plastic_thin_nsys2 + 1, 0.0)};

        std::array<std::vector<G4double>, 2> fPlastic_fatEdep{{std::vector<G4double>(fNofLayers_plastic_fat_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_fat_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_fatLO{{std::vector<G4double>(fNofLayers_plastic_fat_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_fat_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_fatA1{{std::vector<G4double>(fNofLayers_plastic_fat_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_fat_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_fatA2{{std::vector<G4double>(fNofLayers_plastic_fat_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_fat_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_fatT1{{std::vector<G4double>(fNofLayers_plastic_fat_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_fat_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_fatT2{{std::vector<G4double>(fNofLayers_plastic_fat_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_fat_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_fatTrackLength{{std::vector<G4double>(fNofLayers_plastic_fat_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_fat_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_fatToF{{std::vector<G4double>(fNofLayers_plastic_fat_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_fat_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_fatXPos{{std::vector<G4double>(fNofLayers_plastic_fat_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_fat_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_fatYPos{{std::vector<G4double>(fNofLayers_plastic_fat_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_fat_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_fatZPos{{std::vector<G4double>(fNofLayers_plastic_fat_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_fat_nsys2 + 1, 0.0) }};


//////////////

        std::array<std::vector<G4double>, 2> fPlastic_thinEdep{{std::vector<G4double>(fNofLayers_plastic_thin_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_thin_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_thinLO{{std::vector<G4double>(fNofLayers_plastic_thin_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_thin_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_thinA1{{std::vector<G4double>(fNofLayers_plastic_thin_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_thin_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_thinA2{{std::vector<G4double>(fNofLayers_plastic_thin_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_thin_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_thinT1{{std::vector<G4double>(fNofLayers_plastic_thin_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_thin_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_thinT2{{std::vector<G4double>(fNofLayers_plastic_thin_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_thin_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_thinTrackLength{{std::vector<G4double>(fNofLayers_plastic_thin_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_thin_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_thinToF{{std::vector<G4double>(fNofLayers_plastic_thin_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_thin_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_thinXPos{{std::vector<G4double>(fNofLayers_plastic_thin_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_thin_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_thinYPos{{std::vector<G4double>(fNofLayers_plastic_thin_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_thin_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_thinZPos{{std::vector<G4double>(fNofLayers_plastic_thin_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_thin_nsys2 + 1, 0.0) }};


        //////////////

        std::array<std::vector<G4double>, 2> fPlastic_LQ_Edep{{std::vector<G4double>(fNofLayers_plastic_LQ_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_LQ_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_LQ_LO{{std::vector<G4double>(fNofLayers_plastic_LQ_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_LQ_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_LQ_A1{{std::vector<G4double>(fNofLayers_plastic_LQ_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_LQ_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_LQ_A2{{std::vector<G4double>(fNofLayers_plastic_LQ_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_LQ_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_LQ_T1{{std::vector<G4double>(fNofLayers_plastic_LQ_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_LQ_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_LQ_T2{{std::vector<G4double>(fNofLayers_plastic_LQ_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_LQ_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_LQ_TrackLength{{std::vector<G4double>(fNofLayers_plastic_LQ_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_LQ_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_LQ_ToF{{std::vector<G4double>(fNofLayers_plastic_LQ_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_LQ_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_LQ_XPos{{std::vector<G4double>(fNofLayers_plastic_LQ_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_LQ_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_LQ_YPos{{std::vector<G4double>(fNofLayers_plastic_LQ_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_LQ_nsys2 + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fPlastic_LQ_ZPos{{std::vector<G4double>(fNofLayers_plastic_LQ_nsys1 + 1, 0.0), std::vector<G4double>(fNofLayers_plastic_LQ_nsys2 + 1, 0.0) }};



        ////// Для Адронного Калориметра
        

        std::array<std::vector<G4int>, 2> fHCX_N{{std::vector<G4int>(N_HCX + 1, 0.0), std::vector<G4int>(N_HCX + 1, 0.0) }};
        std::array<std::vector<G4int>, 2> fHCX_AL{{std::vector<G4int>(N_HCX + 1, 0.0), std::vector<G4int>(N_HCX + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fHCX_Edep{{std::vector<G4double>(N_HCX + 1, 0.0), std::vector<G4double>(N_HCX + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fHCX_LO{{std::vector<G4double>(N_HCX + 1, 0.0), std::vector<G4double>(N_HCX + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fHCX_A1{{std::vector<G4double>(N_HCX + 1, 0.0), std::vector<G4double>(N_HCX + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fHCX_A2{{std::vector<G4double>(N_HCX + 1, 0.0), std::vector<G4double>(N_HCX + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fHCX_T1{{std::vector<G4double>(N_HCX + 1, 0.0), std::vector<G4double>(N_HCX + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fHCX_T2{{std::vector<G4double>(N_HCX + 1, 0.0), std::vector<G4double>(N_HCX + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fHCX_TrackLength{{std::vector<G4double>(N_HCX + 1, 0.0), std::vector<G4double>(N_HCX + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fHCX_ToF{{std::vector<G4double>(N_HCX + 1, 0.0), std::vector<G4double>(N_HCX + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fHCX_XPos{{std::vector<G4double>(N_HCX + 1, 0.0), std::vector<G4double>(N_HCX + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fHCX_YPos{{std::vector<G4double>(N_HCX + 1, 0.0), std::vector<G4double>(N_HCX + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fHCX_ZPos{{std::vector<G4double>(N_HCX + 1, 0.0), std::vector<G4double>(N_HCX + 1, 0.0) }};

        std::array<std::vector<G4int>, 2> fHCZ_N{{std::vector<G4int>(N_HCZ + 1, 0.0), std::vector<G4int>(N_HCZ + 1, 0.0) }};
        std::array<std::vector<G4int>, 2> fHCZ_AL{{std::vector<G4int>(N_HCZ + 1, 0.0), std::vector<G4int>(N_HCZ + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fHCZ_Edep{{std::vector<G4double>(N_HCZ + 1, 0.0), std::vector<G4double>(N_HCZ + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fHCZ_LO{{std::vector<G4double>(N_HCZ + 1, 0.0), std::vector<G4double>(N_HCZ + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fHCZ_A1{{std::vector<G4double>(N_HCZ + 1, 0.0), std::vector<G4double>(N_HCZ + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fHCZ_A2{{std::vector<G4double>(N_HCZ + 1, 0.0), std::vector<G4double>(N_HCZ + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fHCZ_T1{{std::vector<G4double>(N_HCZ + 1, 0.0), std::vector<G4double>(N_HCZ + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fHCZ_T2{{std::vector<G4double>(N_HCZ + 1, 0.0), std::vector<G4double>(N_HCZ + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fHCZ_TrackLength{{std::vector<G4double>(N_HCZ + 1, 0.0), std::vector<G4double>(N_HCZ + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fHCZ_ToF{{std::vector<G4double>(N_HCZ + 1, 0.0), std::vector<G4double>(N_HCZ + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fHCZ_XPos{{std::vector<G4double>(N_HCZ + 1, 0.0), std::vector<G4double>(N_HCZ + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fHCZ_YPos{{std::vector<G4double>(N_HCZ + 1, 0.0), std::vector<G4double>(N_HCZ + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fHCZ_ZPos{{std::vector<G4double>(N_HCZ + 1, 0.0), std::vector<G4double>(N_HCZ + 1, 0.0) }};

        
///////



        // Для камер


        std::array<std::vector<G4int>, 2> fWa_N{{std::vector<G4int>(NW1_WRS + 1, 0.0), std::vector<G4int>(NW1_WRS + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fWa_XPos{{std::vector<G4double>(NW1_WRS + 1, 0.0), std::vector<G4double>(NW1_WRS + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fWa_YPos{{std::vector<G4double>(NW1_WRS + 1, 0.0), std::vector<G4double>(NW1_WRS + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fWa_ZPos{{std::vector<G4double>(NW1_WRS + 1, 0.0), std::vector<G4double>(NW1_WRS + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fWa_RPos{{std::vector<G4double>(NW1_WRS + 1, 0.0), std::vector<G4double>(NW1_WRS + 1, 0.0) }};


//////

        std::array<std::vector<G4int>, 2> fWb_N{{std::vector<G4int>(NW2_WRS + 1, 0.0), std::vector<G4int>(NW2_WRS + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fWb_XPos{{std::vector<G4double>(NW2_WRS + 1, 0.0), std::vector<G4double>(NW2_WRS + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fWb_YPos{{std::vector<G4double>(NW2_WRS + 1, 0.0), std::vector<G4double>(NW2_WRS + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fWb_ZPos{{std::vector<G4double>(NW2_WRS + 1, 0.0), std::vector<G4double>(NW2_WRS + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fWb_RPos{{std::vector<G4double>(NW2_WRS + 1, 0.0), std::vector<G4double>(NW2_WRS + 1, 0.0) }};

//////

        std::array<std::vector<G4int>, 2> fWc_N{{std::vector<G4int>(NW3_WRS + 1, 0.0), std::vector<G4int>(NW3_WRS + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fWc_XPos{{std::vector<G4double>(NW3_WRS + 1, 0.0), std::vector<G4double>(NW3_WRS + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fWc_YPos{{std::vector<G4double>(NW3_WRS + 1, 0.0), std::vector<G4double>(NW3_WRS + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fWc_ZPos{{std::vector<G4double>(NW3_WRS + 1, 0.0), std::vector<G4double>(NW3_WRS + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fWc_RPos{{std::vector<G4double>(NW3_WRS + 1, 0.0), std::vector<G4double>(NW3_WRS + 1, 0.0) }};


        // Для вершинных камер


        std::array<std::vector<G4int>, 2> fVC_N{{std::vector<G4int>(NVC_WRS + 1, 0.0), std::vector<G4int>(NVC_WRS + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fVC_XPos{{std::vector<G4double>(NVC_WRS + 1, 0.0), std::vector<G4double>(NVC_WRS + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fVC_YPos{{std::vector<G4double>(NVC_WRS + 1, 0.0), std::vector<G4double>(NVC_WRS + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fVC_ZPos{{std::vector<G4double>(NVC_WRS + 1, 0.0), std::vector<G4double>(NVC_WRS + 1, 0.0) }};
        std::array<std::vector<G4double>, 2> fVC_RPos{{std::vector<G4double>(NVC_WRS + 1, 0.0), std::vector<G4double>(NVC_WRS + 1, 0.0) }};


        // data members

        G4int fplastic_fat_nsys1HCID = -1;
        G4int fplastic_fat_nsys2HCID = -1;
        G4int fplastic_thin_nsys1HCID = -1;
        G4int fplastic_thin_nsys2HCID = -1;
        G4int fplastic_LQ_nsys1HCID = -1;
        G4int fplastic_LQ_nsys2HCID = -1;


        G4int fHadronCalorimeter_nsys1HCID = -1;
        G4int fHadronCalorimeter_nsys2HCID = -1;

        G4int fW_Chamber_nsys1HCID = -1; // проволочные (трековые) камеры
        G4int fW_Chamber_nsys2HCID = -1;
        G4int fV_Chamber_nsys1HCID = -1; // вершинные камеры
        G4int fV_Chamber_nsys2HCID = -1;

    };

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


