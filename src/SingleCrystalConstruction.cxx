#include "SingleCrystalConstruction.h"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

//layer visible attributes
//#include "G4VisAttributes.hh"
//#include "G5Colour.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSFlatSurfaceCurrent.hh"
#include "G4PSNofSecondary.hh"
#include "G4PSNofStep.hh"
//#include "G4PSDoseDeposit.hh"
#include "G4SDParticleFilter.hh"
#include "OpticalDetection.h"

#include "HagridCrystal.h"

//#define POLYCONE
//#define DOUBLEWINDOW

SingleCrystalConstruction::SingleCrystalConstruction(G4bool useOptical) :
	G4VUserDetectorConstruction(),
	useOptical_(useOptical)
{ }

SingleCrystalConstruction::~SingleCrystalConstruction()
{ }

G4VPhysicalVolume *SingleCrystalConstruction::Construct()
{
	// Option to switch on/off checking of volumes overlaps
	G4bool checkOverlaps = false;

	// Get nist material manager
	G4NistManager* nist = G4NistManager::Instance();

	// World
	G4double world_sizeXY = 10*cm;
	G4double world_sizeZ  = 100*cm;

	G4double detOffset = 25*cm;

	HagridCrystal *hagrid = new HagridCrystal();

	//-------------------------------------------------------
	//------------ Material Definitions ---------------------
	//-------------------------------------------------------

	//------------ Air --------------------------------------
	G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
	//G4Material* vacuum = nist->BuildMaterialWithNewDensity("Vacuum","G4_AIR",8E-7*g/cm3);

	if (useOptical_) {
		//Create Material Properties Table to store optical properties
		G4MaterialPropertiesTable* air_MPT = new G4MaterialPropertiesTable();

		//Assume constant index of refraction
		const G4int nEntries = 2;
		G4double air_Energy[nEntries] = {3.82495*eV, 2.91425*eV};
		G4double air_rIndex[nEntries] = {1.0002853,1.0002785};
		air_MPT->AddProperty("RINDEX", air_Energy, air_rIndex, nEntries);

		air->SetMaterialPropertiesTable(air_MPT);
	}	


	//-------------------------------------------------------
	//------------ Construct Volumes ------------------------
	//-------------------------------------------------------

	//------------ World Volume -----------------------------
	G4Box* solidWorld =    
		new G4Box("World",                       //its name
				0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);    //its size

	G4LogicalVolume* logicWorld =                         
		new G4LogicalVolume(solidWorld,          //its solid
				air,           //its material
				"World");            //its name

	G4VPhysicalVolume* physWorld = new G4PVPlacement(
			0, 					//no rotation
			G4ThreeVector(),	//at (0,0,0)
			logicWorld,			//contataing logical volume
			"World",				//name
			0,						//no mother volume
			false,				//no boolean operation
			0,						//copy number
			checkOverlaps);	//overlap checking


	hagrid->Construct(logicWorld, G4ThreeVector(0,0,detOffset));//, G4RotationMatrix(0,90.*deg,0));

	//Return the world
	return physWorld;
}


void SingleCrystalConstruction::ConstructSDandField() 
{
	// declare crystal as a MultiFunctionalDetector scorer
	G4MultiFunctionalDetector* cryst = new G4MultiFunctionalDetector("LaBr3Crystal");
	cryst->RegisterPrimitive(new G4PSEnergyDeposit("edep"));
	cryst->RegisterPrimitive(new G4PSNofSecondary("numSecondary"));
	cryst->RegisterPrimitive(new G4PSNofStep("numStep"));
	SetSensitiveDetector("LaBr3Crystal",cryst);

	G4MultiFunctionalDetector* glass = new G4MultiFunctionalDetector("glassFace");
	glass->RegisterPrimitive(new G4PSEnergyDeposit("edep"));
	SetSensitiveDetector("glassFace",glass);

	G4MultiFunctionalDetector* alum = new G4MultiFunctionalDetector("alumExterior");
	alum->RegisterPrimitive(new G4PSEnergyDeposit("edep"));
	SetSensitiveDetector("alumExterior",alum);

#ifdef DOUBLEWINDOW

	G4MultiFunctionalDetector* glass2 = new G4MultiFunctionalDetector("glassFace2");
	glass2->RegisterPrimitive(new G4PSEnergyDeposit("edep"));
	SetSensitiveDetector("glassFace2",glass2);
#endif

	if (useOptical_) {
		G4MultiFunctionalDetector* pmt = new G4MultiFunctionalDetector("PMT");
		pmt->SetFilter(new G4SDParticleFilter("opPhotonFilter","opticalphoton"));
		pmt->RegisterPrimitive(new G4PSEnergyDeposit("edep"));
		pmt->RegisterPrimitive(new OpticalDetection("opDet"));
		G4PSFlatSurfaceCurrent *currentScore = new G4PSFlatSurfaceCurrent("current",fCurrent_In);
		currentScore->DivideByArea(false);
		pmt->RegisterPrimitive(currentScore);
		SetSensitiveDetector("PMT",pmt);
	}
}
