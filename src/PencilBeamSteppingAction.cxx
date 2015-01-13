#include "PencilBeamSteppingAction.h"
#include "G4Step.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4ProcessManager.hh"

PencilBeamSteppingAction::PencilBeamSteppingAction() : 
	G4UserSteppingAction()
{

}
/**When an optical photon is detected its track is terminated without processing the hit.
 * The following UserAction checks if an event is from an optical photon occuring at a 
 * geometrical boundary. If so and the status was \c Detection the ProcessHits of the 
 * corresponding sensitive detector is envoked.
 */
void PencilBeamSteppingAction::UserSteppingAction(const G4Step* step)
{

	//A pointer to the Optical Boundary process.
	//	This variable should not change between function calls
	static G4ThreadLocal G4OpBoundaryProcess* boundaryProcess=NULL;

	//find the boundaryProcess process only once
	if(!boundaryProcess){
		G4ProcessManager* pm
			= step->GetTrack()->GetDefinition()->GetProcessManager();
		G4int nprocesses = pm->GetProcessListLength();
		G4ProcessVector* pv = pm->GetProcessList();
		for(int i=0;i<nprocesses;i++){
			if((*pv)[i]->GetProcessName()=="OpBoundary"){
				boundaryProcess = (G4OpBoundaryProcess*)(*pv)[i];
				break;
			}
		}
	}

	G4ParticleDefinition* particleType = step->GetTrack()->GetDefinition();
	//Check that the particle is a optical phone
	if(particleType==G4OpticalPhoton::OpticalPhotonDefinition()){
		//Check to see if the partcile was actually at a boundary
		//Otherwise the boundaryProcess status may not be valid
		if(step->GetPostStepPoint()->GetStepStatus()==fGeomBoundary){
			G4OpBoundaryProcessStatus boundaryStatus = boundaryProcess->GetStatus();
			//Check that the status of the OpBoundary process is Detection
			if (boundaryStatus == Detection) {
				//Get volume of the current step
				G4LogicalVolume* volume 
					= step->GetPreStepPoint()->GetTouchableHandle()
					->GetVolume()->GetLogicalVolume();
				G4String volName = volume->GetName();
				//Get the post step volume.
				volume 
					= step->GetPostStepPoint()->GetTouchableHandle()
					->GetVolume()->GetLogicalVolume();
				G4String volNamePost = volume->GetName();

				G4double edepStep = step->GetTotalEnergyDeposit();
				fEnDep[volName] += edepStep;

				printf("Edep %s->%s %f %f\n",volName.data(),volNamePost.data(),edepStep,fEnDep[volName]);

			}
		}
	}

}
