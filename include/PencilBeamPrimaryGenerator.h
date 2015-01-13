#ifndef PENCILBEAMPRIMARYGENERATOR_H
#define PENCILBEAMPRIMARYGENERATOR_H

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"

class PencilBeamPrimaryGenerator : public G4VUserPrimaryGeneratorAction
{
	private:
		G4ParticleGun* fParticleGun;
		
	public:
		PencilBeamPrimaryGenerator();
		virtual ~PencilBeamPrimaryGenerator();

		virtual void GeneratePrimaries(G4Event*);

		const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

};

#endif
