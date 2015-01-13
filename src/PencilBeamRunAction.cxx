#include "PencilBeamRunAction.h"
#include "PencilBeamRun.h"
#include "G4SDManager.hh"

PencilBeamRunAction::PencilBeamRunAction() :
	G4UserRunAction(),
	fManager(0),
	fSDoffset(7)
{
}

void PencilBeamRunAction::BeginOfRunAction(const G4Run*)
{
	//Construct analysis manager and save pointer
	fManager = G4RootAnalysisManager::Instance();
	fManager->CreateNtuple("test","Testing");
	fManager->CreateNtupleDColumn("primaryX");
	fManager->CreateNtupleDColumn("primaryY");
	fManager->CreateNtupleDColumn("primaryZ");
	fManager->CreateNtupleDColumn("primaryPx");
	fManager->CreateNtupleDColumn("primaryPy");
	fManager->CreateNtupleDColumn("primaryPz");
	fManager->CreateNtupleDColumn("primaryTotEn");

	//Get the number of Sensitive Detectors
	G4HCtable *hitCollectionTable = G4SDManager::GetSDMpointer()->GetHCtable();
	G4int numSDentries = hitCollectionTable->entries();

	if (numSDentries==0) 
		printf("WARNING: No sensitive detectors found.\n");
	else {

		printf("Adding the following sensitive detectors for ROOT output:\n");
		//Loop over each sensitive detector creating a column in the ntuple
		for (int i=0;i<numSDentries;i++) {
			G4String sdName = hitCollectionTable->GetSDname(i);
			G4String hcName = hitCollectionTable->GetHCname(i);
			G4String ntupleName = sdName + "_" + hcName;
			printf("\t%s/%s\n",sdName.c_str(),hcName.c_str());
			fManager->CreateNtupleDColumn(ntupleName);
		}
	}
	fManager->FinishNtuple();

	fManager->OpenFile("test.root");
}

void PencilBeamRunAction::EndOfRunAction(const G4Run*) 
{
	if (fManager) {
		fManager->Write();
		fManager->CloseFile();
		delete fManager;
	}
}
PencilBeamRunAction::~PencilBeamRunAction() 
{

}


G4Run* PencilBeamRunAction::GenerateRun() {
	return new PencilBeamRun(fSDoffset);
}



