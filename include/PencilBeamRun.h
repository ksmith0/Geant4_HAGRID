#ifndef PENCILBEAMRUN_H
#define PENCILBEAMRUN_H

#include "G4Run.hh"
#include "G4Event.hh"

class PencilBeamRun : public G4Run {

	public:
		PencilBeamRun(G4int sdOffset);
		virtual ~PencilBeamRun() {};
		///Implementation of event data accumulation.
		void RecordEvent(const G4Event* event);

	private:
		///The number of ntuple columns proceeding the sensitive detectors.
		G4int fSDoffset;
};

#endif
