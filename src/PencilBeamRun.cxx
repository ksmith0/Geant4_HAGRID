#include "PencilBeamRun.h"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4RootAnalysisManager.hh"

PencilBeamRun::PencilBeamRun(G4int sdOffset) :
	G4Run(),
	fSDoffset(sdOffset)
{

}

void PencilBeamRun::RecordEvent(const G4Event *event)
{
	//Get G4 Analysis Manager
	G4RootAnalysisManager* analysisManager = G4RootAnalysisManager::Instance();
	//Get information about primary vertex
	G4PrimaryVertex *primary = event->GetPrimaryVertex();
	G4PrimaryParticle *primaryParticle = primary->GetPrimary();
	analysisManager->FillNtupleDColumn(0,primary->GetX0());
	analysisManager->FillNtupleDColumn(1,primary->GetY0());
	analysisManager->FillNtupleDColumn(2,primary->GetZ0());
	analysisManager->FillNtupleDColumn(3,primaryParticle->GetPx());
	analysisManager->FillNtupleDColumn(4,primaryParticle->GetPy());
	analysisManager->FillNtupleDColumn(5,primaryParticle->GetPz());
	analysisManager->FillNtupleDColumn(6,primaryParticle->GetTotalEnergy());

	//Get the hit collection for this event
	G4HCofThisEvent *hitCollection = event->GetHCofThisEvent();
	if (!hitCollection) return;

	//Get information about hits in sensitive detectors
	G4int numCollections = hitCollection->GetNumberOfCollections();
	for (int i=0;i<numCollections;i++) {
		//Get event map for this collection (sensitive detector)
		G4THitsMap<G4double>* evtMap = 
			static_cast<G4THitsMap<G4double>*>(hitCollection->GetHC(i));

		std::map<G4int,G4double*>::iterator itr;
		//Loop over each event in the map fo this collection (sensitive detector)
		for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
			
			G4double value = *(itr->second);
			//fill column offsetting for primary particle info
			analysisManager->FillNtupleDColumn(i + fSDoffset, value);
		}
	}

	analysisManager->AddNtupleRow();  

	
}
