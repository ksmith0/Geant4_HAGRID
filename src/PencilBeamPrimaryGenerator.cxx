#include "PencilBeamPrimaryGenerator.h"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

PencilBeamPrimaryGenerator::PencilBeamPrimaryGenerator() :
	G4VUserPrimaryGeneratorAction(),
	fParticleGun(0)
{
	fParticleGun = new G4ParticleGun(1);

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = particleTable->FindParticle("gamma");
	//G4ParticleDefinition* particle = particleTable->FindParticle("e-");
	fParticleGun->SetParticleDefinition(particle);
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
	fParticleGun->SetParticleEnergy(1.*MeV);
}

PencilBeamPrimaryGenerator::~PencilBeamPrimaryGenerator()
{
	delete fParticleGun;
}

void PencilBeamPrimaryGenerator::GeneratePrimaries(G4Event* anEvent)
{
	
	G4double envSizeXY = 0;
	G4double envSizeZ = 0;
	
	G4double size = 0.8; 
	G4double x0 = size * envSizeXY * (G4UniformRand()-0.5);
	G4double y0 = size * envSizeXY * (G4UniformRand()-0.5);
	G4double z0 = -0.5 * envSizeZ;
	
	fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
	//fParticleGun->SetParticleEnergy(pow(10.,(3*G4UniformRand()-2))*MeV);
	
	fParticleGun->GeneratePrimaryVertex(anEvent);

}
