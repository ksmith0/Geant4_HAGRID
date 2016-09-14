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
	fParticleGun->SetParticleEnergy(0.66167*MeV);
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

	G4ThreeVector direction;

/*
	//Isotropic emission
 	phi = acos(2 * G4UniformRand() - 1); //polar
	theta = 2 * 3.14159 * G4UniformRand();
	direction = G4ThreeVector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
*/

	//Isotropic limited to detector.
	// First we make a cone pointed at polar angle = 0 with constraints on detector size.
	//  Polar angle varies from 0 to size of detector. Azimuthal angle varies across 2 pi.
	//  theta = Polar Angle 
	//  phi = Azimuthal Angle
	G4double detRadius = 55.8*mm / 2.;
	G4double detDistance = 25*cm;
	G4double	detTheta = 0.*deg;
	G4double	detPhi = 0.*deg;
	G4double cosThetaMax = cos(atan(detRadius / detDistance));
	theta = acos((1 - cosThetaMax) * G4UniformRand() + cosThetaMax);
	phi = 2 * 3.14159 * G4UniformRand();
	direction = G4ThreeVector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
	//The cone is then transformed to the correct angle.
	direction.transform(G4RotationMatrix(0,detTheta,90.*deg + detPhi));

	//Then we can set the direction of momentum.
	fParticleGun->SetParticleMomentumDirection(direction);
	
	//Random energy
	G4double energy = pow(10.,(2*G4UniformRand()-1))*MeV;
	energy = 0.66167*MeV;
	//Relativistic Doppler Broadening.
	G4double beta = 0.3;
	//doppler shift based on polar angle theta.
	energy *= sqrt(1 - pow(beta,2)) / (1 - beta * cos(direction.theta()));
	fParticleGun->SetParticleEnergy(energy);

	fParticleGun->GeneratePrimaryVertex(anEvent);

	G4RootAnalysisManager* analysisManager = G4RootAnalysisManager::Instance();
	analysisManager->FillNtupleDColumn(7,beta);

}
