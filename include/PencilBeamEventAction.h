#ifndef PENCILBEAMEVENTACTION_H
#define PENCILBEAMEVENTACTION_H

#include "G4UserEventAction.hh"
#include "G4Event.hh"

class PencilBeamEventAction : public G4UserEventAction
{
	public:
		PencilBeamEventAction();
		virtual ~PencilBeamEventAction() {};
		
		virtual void EndOfEventAction(const G4Event *evt);
};

#endif
