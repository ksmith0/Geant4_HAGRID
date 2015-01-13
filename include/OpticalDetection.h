#ifndef OPTICALDETECTION_H
#define OPTICALDETECTION_H

#include "G4VPrimitiveScorer.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"

class OpticalDetection : public G4VPrimitiveScorer
{
	private:
      G4int  fHCID;
      G4THitsMap<G4double>* fEvtMap;

	public:
		OpticalDetection(G4String name, G4int depth=0);
		virtual ~OpticalDetection() {}; 
      virtual void Initialize(G4HCofThisEvent*);
      virtual void EndOfEvent(G4HCofThisEvent*);

	protected:
		virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

};

#endif
