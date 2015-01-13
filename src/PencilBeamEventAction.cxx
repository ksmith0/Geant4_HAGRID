#include "PencilBeamEventAction.h"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"

PencilBeamEventAction::PencilBeamEventAction() :
	G4UserEventAction()
{

}
void PencilBeamEventAction::EndOfEventAction(const G4Event *evt) {


	G4HCofThisEvent *hitCollection = evt->GetHCofThisEvent();
	if (!hitCollection) return;

	G4int numCollections = hitCollection->GetNumberOfCollections();
	for (int i=0;i<numCollections;i++) {
		G4THitsMap<G4double>* evtMap = 
			static_cast<G4THitsMap<G4double>*>(hitCollection->GetHC(i));

		std::map<G4int,G4double*>::iterator itr;
		for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
			G4double edep = *(itr->second);
			printf("Collection %d edep %f\n",i,edep);
		}

	}
}
