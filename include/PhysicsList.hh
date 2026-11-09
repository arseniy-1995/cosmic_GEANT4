//
// Created by yurchenko on 06.11.2026.
//

#ifndef PHYSICSLIST_HH
#define PHYSICSLIST_HH

#include "G4VModularPhysicsList.hh"
#include "G4VUserPhysicsList.hh"

#include "FTFP_BERT.hh"
#include "G4OpticalPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4SystemOfUnits.hh"

#include "G4VDiscreteProcess.hh"

class MyPhysicsList : public G4VModularPhysicsList
{
    // class MyPhysicsList: public FTFP_BERT {
    // class MyPhysicsList: public G4VUserPhysicsList {
public:
    // MyPhysicsList() : G4VUserPhysicsList() {};
    // MyPhysicsList();
    MyPhysicsList();
    virtual ~MyPhysicsList() {};
    // ~MyUserPhysicsList();


protected:
    void SetCuts() override;

    void ConstructProcess() override;
    void ConstructParticle() override;
};


class MyCustomInteraction : public G4VDiscreteProcess
{
public:
    MyCustomInteraction(const G4String &name = "MyCustomProcess") : G4VDiscreteProcess(name) {}

    // Возвращает длину свободного пробега частицы
    virtual G4double GetMeanFreePath(const G4Track &, G4double, G4ForceCondition *) override
    {
        return 10.0 * cm; // Или расчет на основе сечения взаимодействия
    }

    // Описывает, что происходит в точке взаимодействия (финальное состояние)
    virtual G4VParticleChange *PostStepDoIt(const G4Track &track, const G4Step &step) override
    {
        aParticleChange.Initialize(track);

        // Пример: уничтожаем первичную частицу
        aParticleChange.ProposeTrackStatus(fStopAndKill);

        return &aParticleChange;
    }
};

#endif // START_TIME_DETECTOR_PHYSICSLIST_HH
