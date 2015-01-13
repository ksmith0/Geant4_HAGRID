#include "PhysicsList.h"

//Particle Types
#include "G4ParticleDefinition.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4Gamma.hh"
#include "G4OpticalPhoton.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

//Helper for process registration
#include "G4PhysicsListHelper.hh"

//Gamma Interactions
#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"

//Electron/Positron Interactions
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

//Optical Interactions
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpBoundaryProcess.hh"
//#include "G4OpMieHG.hh"


PhysicsList::PhysicsList() :
	G4VUserPhysicsList()
{
 
}

void PhysicsList::ConstructParticle() {
	G4BosonConstructor bConstructor;
	bConstructor.ConstructParticle();

	G4LeptonConstructor lConstructor;
	lConstructor.ConstructParticle();

	G4MesonConstructor mConstructor;
	mConstructor.ConstructParticle();

	G4BaryonConstructor rConstructor;
	rConstructor.ConstructParticle();

	G4IonConstructor iConstructor;
	iConstructor.ConstructParticle(); 

}
void PhysicsList::ConstructProcess() {
	AddTransportation();
	ConstructEM();
	ConstructOp();
}

void PhysicsList::ConstructOp() {
	
	// Get pointer to G4PhysicsListHelper
	G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

	// Construct and register processes for gamma
	G4ParticleDefinition* particle = G4Gamma::GammaDefinition();
	ph->RegisterProcess(new G4Scintillation, particle);

 	particle = G4Electron::ElectronDefinition();
	ph->RegisterProcess(new G4Scintillation, particle);
/*
	particle = G4Positron::PositronDefinition();
	ph->RegisterProcess(new G4Scintillation, particle);
*/
	
	// Construct and register processes for optical photons
	particle = G4OpticalPhoton::OpticalPhotonDefinition();
	ph->RegisterProcess(new G4OpAbsorption, particle);
	ph->RegisterProcess(new G4OpRayleigh, particle);
	ph->RegisterProcess(new G4OpBoundaryProcess, particle);
	//ph->RegisterProcess(new G4OpMieHG, particle);

}

void PhysicsList::ConstructEM() {
	// Get pointer to G4PhysicsListHelper
	G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

	// Construct and register processes for gamma
	G4ParticleDefinition* particle = G4Gamma::GammaDefinition();
	ph->RegisterProcess(new G4PhotoElectricEffect(), particle);
	ph->RegisterProcess(new G4ComptonScattering(), particle);
	ph->RegisterProcess(new G4GammaConversion(), particle);

	//Construct and register electron processes
	particle = G4Electron::ElectronDefinition();
	ph->RegisterProcess(new G4eMultipleScattering(),particle);
	ph->RegisterProcess(new G4eIonisation(),particle);
	ph->RegisterProcess(new G4eBremsstrahlung(),particle);

	//Construct and register positron processes
	particle = G4Positron::PositronDefinition();
	ph->RegisterProcess(new G4eMultipleScattering(),particle);
	ph->RegisterProcess(new G4eIonisation(),particle);
	ph->RegisterProcess(new G4eBremsstrahlung(),particle);
	ph->RegisterProcess(new G4eplusAnnihilation(),particle);
}
