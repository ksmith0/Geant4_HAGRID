#include "PencilBeamPrimaryGenerator.h"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4RotationMatrix.hh"
#include "Randomize.hh"

#include "G4RootAnalysisManager.hh"

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
	fParticleGun->SetParticleEnergy(1.333*MeV);
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

	G4double phi = 0; //polar
	G4double theta = 0; //azimuthal

	//Isotropic emission
 	phi = acos(2 * G4UniformRand() - 1); //polar
	theta = 2 * 3.14159 * G4UniformRand();

	//Isotropic limited to detector.
	G4double detRadius = 55.8*mm / 2.;
	G4double detDistance = 25*cm;
	G4double	detPhi = 180.*deg;
	G4double	detTheta = 0.*deg;
	phi =  acos(cos((2 * G4UniformRand() -1) * (atan(detRadius / detDistance))));
	theta = 2 * 3.14159 * G4UniformRand();
	
	G4ThreeVector direction(sin(phi)*cos(theta),sin(phi)*sin(theta),cos(phi));
	direction.transform(G4RotationMatrix(0,detPhi,90.*deg + detTheta));
	fParticleGun->SetParticleMomentumDirection(direction);
	//fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(phi)*cos(theta),sin(phi)*sin(theta),cos(phi)));
	
	//Random energy
	G4double energy = pow(10.,(2*G4UniformRand()-1))*MeV;
	energy = 1.333*MeV;
	//Relativistic Doppler Broadening.
	G4double beta = 0.7;
	energy *= sqrt(1 - beta) / (1 - beta * cos(phi));
	fParticleGun->SetParticleEnergy(energy);
//
	fParticleGun->GeneratePrimaryVertex(anEvent);

	G4RootAnalysisManager* analysisManager = G4RootAnalysisManager::Instance();
	analysisManager->FillNtupleDColumn(7,beta);

}
