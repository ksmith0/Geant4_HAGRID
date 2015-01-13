#ifndef PENCILBEAMSTEPPINGACTION_H
#define PENCILBEAMSTEPPINGACTION_H

#include "G4UserSteppingAction.hh"
#include "G4LogicalVolume.hh"

class PencilBeamSteppingAction : public G4UserSteppingAction
{
	private:
		std::map< G4String, G4double > fEnDep;
		G4LogicalVolume* fScoringVolume;

	public:
		PencilBeamSteppingAction();
		virtual ~PencilBeamSteppingAction() {};

		virtual void UserSteppingAction(const G4Step*);
		void SetScoringVolume(G4LogicalVolume* scoringVolume)
			{fScoringVolume = scoringVolume;}

};

#endif

