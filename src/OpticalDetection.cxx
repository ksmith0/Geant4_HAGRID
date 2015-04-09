#include "OpticalDetection.h"

#include "G4ProcessManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"

OpticalDetection::OpticalDetection(G4String name, G4int depth) :
	G4VPrimitiveScorer(name, depth),
	fHCID(-1)
{

}

G4bool OpticalDetection::ProcessHits(G4Step* step, G4TouchableHistory*)
{
	//A pointer to the Optical Boundary process.
	//	This variable should not change between function calls
	static G4ThreadLocal G4OpBoundaryProcess* opBoundaryProcess=NULL;

	//find the boundaryProcess process only once
	if(!opBoundaryProcess){
		G4ProcessManager* pm
			= step->GetTrack()->GetDefinition()->GetProcessManager();
		G4int nprocesses = pm->GetProcessListLength();
		G4ProcessVector* pv = pm->GetProcessList();
		for(int i=0;i<nprocesses;i++){
			if((*pv)[i]->GetProcessName()=="OpBoundary"){
				opBoundaryProcess = (G4OpBoundaryProcess*)(*pv)[i];
				break;
			}
		}
	}

	//Check that the particle is a optical phone
	G4ParticleDefinition* particleType = step->GetTrack()->GetDefinition();
	if(particleType == G4OpticalPhoton::OpticalPhotonDefinition()){
		if(step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary){
			//Check that the status of the OpBoundary process is Detection
			if (opBoundaryProcess->GetStatus() == Detection) {
				G4int index = GetIndex(step);
				G4double value = 1;
				fEvtMap->add(index,value);
			}
		}
	}

	return true;
}

void OpticalDetection::Initialize(G4HCofThisEvent* HCE)
{
  fEvtMap = new G4THitsMap<G4double>(detector->GetName(), GetName());
  if ( fHCID < 0 ) fHCID = GetCollectionID(0);
  HCE->AddHitsCollection(fHCID, (G4VHitsCollection*)fEvtMap);
}

void OpticalDetection::EndOfEvent(G4HCofThisEvent*)
{

}

