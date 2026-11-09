//
// Created by yurchenko on 06.11.2026.
//

#include "PhysicsList.hh"
#include "G4SystemOfUnits.hh"

#include "G4ProcessManager.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserPhysicsList.hh"

#include "G4EmStandardPhysics.hh"

// Процессы для электронов/позитронов
#include <G4MuonMinus.hh>
#include <G4MuonPlus.hh>
#include <G4OpticalPhoton.hh>

#include "G4DecayPhysics.hh"
#include "G4eBremsstrahlung.hh" // Класс тормозного излучения
#include "G4eIonisation.hh"
#include "G4eMultipleScattering.hh"
#include "G4eplusAnnihilation.hh"

#include "G4ComptonScattering.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4Positron.hh"

#include "FTFP_BERT.hh"
#include "G4OpticalPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

MyPhysicsList::MyPhysicsList()
{

    SetVerboseLevel(1);

    // Добавляем стандартную электромагнитную физику
    RegisterPhysics(new G4EmStandardPhysics());
    // Добавляем процессы распада частиц
    RegisterPhysics(new G4DecayPhysics());


    // Добавляем регистрацию оптических процессов сверху базового списка
    RegisterPhysics(new G4OpticalPhysics());
    // Добавляем радиоактивный распад
    RegisterPhysics(new G4RadioactiveDecayPhysics());

    // G4ProcessManager* pManager = G4Electron::ElectronDefinition()->GetProcessManager();
    // Добавление процесса тормозного излучения
    // G4eBremsstrahlung* brems = new G4eBremsstrahlung();
    // pManager->AddDiscreteProcess(brems);
}
// 1. Регистрируем частицы
void MyPhysicsList::ConstructParticle()
{

    G4Gamma::GammaDefinition();
    G4Electron::ElectronDefinition();
    G4Positron::PositronDefinition();
    G4MuonMinus::MuonMinusDefinition();
    G4MuonPlus::MuonPlusDefinition();
    G4OpticalPhoton::OpticalPhotonDefinition();
}
// 2. Назначаем процессы
void MyPhysicsList::ConstructProcess()
{
    // Обязательный процесс для перемещения всех частиц в геометрии
    // Добавление базового транспорта (обязательно)
    AddTransportation();

    // Пример ручного добавления EM-процессов для Гамма-квантов

    // Получаем указатель на частицы через итератор
    auto particleIterator = GetParticleIterator();
    particleIterator->reset();


    while ((*particleIterator)())
    {
        G4ParticleDefinition *particle = particleIterator->value();
        G4ProcessManager *pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();

        if (particleName == "gamma")
        {
            pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
            pmanager->AddDiscreteProcess(new G4ComptonScattering);

            // Добавляем наш собственный процесс
            // pmanager->AddDiscreteProcess(new MyCustomInteraction());
        }

        // Примечание: Числа -1, 3, 3 задают порядок вызова процессов на стадиях AtRestDoIt, AlongStepDoIt и
        // PostStepDoIt соответственно.

        if (particleName == "e-")
        {
            // Электрон: множественное рассеяние, ионизация и тормозное излучение
            pmanager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
            pmanager->AddProcess(new G4eIonisation, -1, 2, 2);
            pmanager->AddProcess(new G4eBremsstrahlung, -1, 3, 3); // Подключение тормозного излучения
        }
        else if (particleName == "e+")
        {
            // Позитрон: то же самое + аннигиляция
            pmanager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
            pmanager->AddProcess(new G4eIonisation, -1, 2, 2);
            pmanager->AddProcess(new G4eBremsstrahlung, -1, 3, 3); // Подключение тормозного излучения
            pmanager->AddProcess(new G4eplusAnnihilation, 0, -1, 4);
        }
    }
}

// 3. Задаем пространственные пороги (Cuts)
// Установка порогов рождения вторичных частиц
void MyPhysicsList::SetCuts()
{
    // Порог в 1 мм означает: фотоны с пробегом < 1 мм не рождаются как отдельные треки,
    // а их энергия учитывается локально как потери электрона (Continuous energy loss).
    // Запретить рождение вторичных частиц, если их пробег меньше 1 мм
    SetCutsWithDefault();
    // Или вручную: SetCutValue(1.0*mm, "gamma");

    // Задается в единицах длины, Geant4 сам пересчитает это в энергию для каждого материала
    SetCutsWithDefault();
    SetCutValue(1.0 * mm,
                "gamma"); // Фотоны тормозного излучения с пробегом < 1 мм не будут рождаться как отдельные частицы
    SetCutValue(1.0 * mm, "e-");
}
