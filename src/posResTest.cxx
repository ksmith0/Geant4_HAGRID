/**A GEANT4 simulation of a single BrilLanCe 380 cyrstal.
 *
 *
 */

#ifdef G4MULTITHREADED
	#include "G4MTRunManager.hh"
#else
	#include "G4RunManager.hh"
#endif

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "G4UImanager.hh"
#include "PhysicsList.h"

#include "SingleCrystalConstruction.h"
#include "PosResTestActionInitialization.h"

int main(int argc, char** argv)
{

#ifdef G4MULTITHREADED
	G4MTRunManager* runManager = new G4MTRunManager;
#else
	G4RunManager* runManager = new G4RunManager;
#endif

	//Initialize detector construction
	runManager->SetUserInitialization(new SingleCrystalConstruction(false));
	//runManager->SetUserInitialization(new ArgonneBetaDetConstruction(true));

	//Initialize physics list
	runManager->SetUserInitialization(new PhysicsList());
	
	//Initialize user action
	runManager->SetUserInitialization(new PosResTestActionInitialization());
	
	//Initialize kernel
	runManager->Initialize();
  
#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (argc!=1) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode : define UI session
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute share/init_vis.mac"); 
#else
    UImanager->ApplyCommand("/control/execute share/init.mac"); 
#endif
    ui->SessionStart();
    delete ui;
#endif
  }

#ifdef G4VIS_USE
  delete visManager;
#endif
	delete runManager;
	
	return 0;
}
