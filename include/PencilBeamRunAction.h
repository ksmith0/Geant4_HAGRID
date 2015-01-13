#ifndef PENCILBEAMRUNACTION_H
#define PENCILBEAMRUNACTION_H

#include "G4UserRunAction.hh"
#include "G4Event.hh"
#include "G4RootAnalysisManager.hh"

class PencilBeamRunAction : public G4UserRunAction
{
	public:
		PencilBeamRunAction();
		virtual ~PencilBeamRunAction();

		void BeginOfRunAction(const G4Run *run);
		void EndOfRunAction(const G4Run *run);
		G4Run* GenerateRun();
	private:
		///Pointer to analysis manager.
		G4RootAnalysisManager *fManager;
		///The number of ntuple columns proceeding the sensitive detectors.
		G4int fSDoffset;
};

#endif

