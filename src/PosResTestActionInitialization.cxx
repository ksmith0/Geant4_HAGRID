#include "PosResTestActionInitialization.h"
#include "PencilBeamPrimaryGenerator.h"
#include "PencilBeamSteppingAction.h"
#include "PencilBeamRunAction.h"
//#include "PencilBeamEventAction.h"

PosResTestActionInitialization::PosResTestActionInitialization() :
	G4VUserActionInitialization()
{ }

PosResTestActionInitialization::~PosResTestActionInitialization()
{ }

void PosResTestActionInitialization::BuildForMaster() const
{

}

void PosResTestActionInitialization::Build() const
{
	SetUserAction(new PencilBeamPrimaryGenerator);
	//SetUserAction(new PencilBeamSteppingAction);
	SetUserAction(new PencilBeamRunAction);	
	//SetUserAction(new PencilBeamEventAction);	

}
