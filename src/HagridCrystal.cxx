#include "HagridCrystal.h"

#include "G4NistManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

#include "G4PVPlacement.hh"
#include "G4Tubs.hh"

#include "G4LogicalBorderSurface.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

HagridCrystal::HagridCrystal():
	useOptical_(false),
	checkOverlaps_(false)
{

	BuildMaterials();
}

HagridCrystal::~HagridCrystal() {

}

void HagridCrystal::BuildMaterials() {

	// Get nist material manager
	G4NistManager* nist = G4NistManager::Instance();

	//------------ Elements ---------------------------------
	G4Element* La = nist->FindOrBuildElement("La");
	G4Element* Br = nist->FindOrBuildElement("Br");
	G4Element* Ce = nist->FindOrBuildElement("Ce"); 

	G4Element* Sb = nist->FindOrBuildElement("Sb"); 
	G4Element* Rb = nist->FindOrBuildElement("Rb"); 
	G4Element* Cs = nist->FindOrBuildElement("Cs"); 

	//------------ Air --------------------------------------
	air_ = nist->FindOrBuildMaterial("G4_AIR");

	if (useOptical_) {
		//Create Material Properties Table to store optical properties
		G4MaterialPropertiesTable* airMPT = new G4MaterialPropertiesTable();

		//Assume constant index of refraction
		const G4int nEntries = 2;
		G4double airEnergy[nEntries] = {3.82495*eV, 2.91425*eV};
		G4double airRindex[nEntries] = {1.0002853,1.0002785};
		airMPT->AddProperty("RINDEX", airEnergy, airRindex, nEntries);

		air_->SetMaterialPropertiesTable(airMPT);
	}	

	
	//--------------- Aluminum -------------------------------
	aluminum_ = nist->FindOrBuildMaterial("G4_Al");

	//--------------- Teflon ---------------------------------
	teflon_ = nist->FindOrBuildMaterial("G4_TEFLON");

	teflonSurface_ = NULL;
	if (useOptical_) {
		G4MaterialPropertiesTable * teflon_MPT = new G4MaterialPropertiesTable();

		//\note Need to find a good source of teflon optical properties.
		//Assume constant index of refraction
		const G4int nEntries = 2;
		G4double teflon_Energy[nEntries] = {3.82495*eV, 2.91425*eV};

		G4double teflon_rIndex[nEntries] = {1.35,1.35};
		G4double teflon_absorption[nEntries] = {2.*mm,2.*mm};
		teflon_MPT->AddProperty("RINDEX", teflon_Energy, teflon_rIndex, nEntries);
		teflon_MPT->AddProperty("ABSLENGTH", teflon_Energy, teflon_absorption, nEntries);

		teflon_->SetMaterialPropertiesTable(teflon_MPT);

		//Create surface properties
		G4MaterialPropertiesTable *teflonSurface_MPT = new G4MaterialPropertiesTable();
		G4double teflon_reflectivity[nEntries] = {0.98, 0.98};
		teflonSurface_MPT->AddProperty("REFLECTIVITY", teflon_Energy, teflon_reflectivity, nEntries);

		//Create the optical surface
		teflonSurface_ = new G4OpticalSurface("teflonSurface");
		teflonSurface_->SetType(dielectric_dielectric);
		//GroundFrontPainted finish specifies only reflection or absorption. 100% Lambertian
		teflonSurface_->SetFinish(groundfrontpainted);
		teflonSurface_->SetMaterialPropertiesTable(teflonSurface_MPT);
	}

	//--------------- LaBr3(Ce) ------------------------------
	G4Material* LaBr3 = new G4Material("LaBr3",5.07*g/cm3,2);
	LaBr3->AddElement(La,1);
	LaBr3->AddElement(Br,3);

	//Create LaBr3(Ce) material by doping LaBr3
	laBr3_Ce_ = new G4Material("laBr3_Ce_",5.08*g/cm3, 2);
	laBr3_Ce_->AddMaterial(LaBr3, 95*perCent);
	laBr3_Ce_->AddElement(Ce, 5*perCent);

	if (useOptical_) {
		//Create Material Properties Table to store optical properties
		G4MaterialPropertiesTable* LaBrCe_MPT = new G4MaterialPropertiesTable();

		// Data from 10.1109/TNS.2012.2193597
		const G4int nEntries = 21;
		G4double LaBr_Energy[nEntries] =  
		{3.82495*eV, 3.76611*eV, 3.70904*eV, 3.65368*eV, 3.59995*eV, 
			3.54778*eV, 3.4971*eV, 3.44784*eV, 3.39996*eV, 3.35338*eV, 
			3.30807*eV, 3.26396*eV, 3.22101*eV, 3.17918*eV, 3.13842*eV, 
			3.09869*eV, 3.05996*eV, 3.02218*eV, 2.98533*eV, 2.94936*eV, 
			2.91425*eV};
		//Fast compoennt
		G4double LaBr_Fast[nEntries]  =  
		{0.00171744, 0.0169987, 0.104353, 0.413358, 1.09586, 2.01313, 2.65712, 
			2.64391, 2.16732, 1.73403, 1.60549, 1.65549, 1.63166, 1.41242, 
			1.05117, 0.674215, 0.375993, 0.184151, 0.0799886, 0.0310985, 
			0.0109151};     
		LaBrCe_MPT->AddProperty("FASTCOMPONENT",LaBr_Energy, LaBr_Fast, nEntries);
		//Index of reflection
		G4double LaBr_RIND[nEntries] =
			{2.53548, 2.50933, 2.48529, 2.4631, 2.44257, 2.42352, 2.40581, 
			2.38929, 2.37386, 2.35942, 2.34587, 2.33314, 2.32116, 2.30986, 
			2.2992, 2.28912, 2.27958, 2.27054, 2.26196, 2.25381, 2.24605} ;
		LaBrCe_MPT->AddProperty("RINDEX",LaBr_Energy, LaBr_RIND,    nEntries);
		//Absorbstion length
		G4double LaBr_ABSL[nEntries] = 
			{0.00214455*cm, 0.00485183*cm, 0.0130231*cm, 0.0408024*cm, 
			0.147022*cm, 0.601092*cm, 2.75427*cm, 13.9861*cm, 77.9037*cm, 
			471.533*cm, 3074.95*cm, 21435.2*cm, 158585.*cm, 1.23706E6*cm, 
			1.01134E7*cm, 8.61758E7*cm, 7.61483E8*cm, 6.94559E9*cm, 
			6.51156E10*cm, 6.25021E11*cm, 6.12049E12*cm};
		/*	
			{35.*cm, 35.*cm, 35.*cm, 35.*cm, 35.*cm, 
			35.*cm, 35.*cm, 35.*cm, 35.*cm, 35.*cm, 
			35.*cm, 35.*cm, 35.*cm, 35.*cm, 
			35.*cm, 35.*cm, 35.*cm, 35.*cm,
			35.*cm, 35.*cm, 35.*cm};
		*/
		/*
			{1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 
			1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 
			1.*mm, 1.*mm, 1.*mm, 1.*mm, 
			1.*mm, 1.*mm, 1.*mm, 1.*mm,
			1.*mm, 1.*mm, 1.*mm};
		*/	

		LaBrCe_MPT->AddProperty("ABSLENGTH", LaBr_Energy, LaBr_ABSL, nEntries);

		G4double LaBr_RayleighAbs[nEntries] = 
			{8.736*cm, 9.035*cm, 9.334*cm, 9.633*cm, 9.932*cm, 10.231*cm, 
			10.53*cm, 10.829*cm, 11.128*cm, 11.427*cm, 11.726*cm, 12.025*cm, 
			12.324*cm, 12.623*cm, 12.922*cm, 13.221*cm, 13.52*cm, 13.819*cm, 
			14.118*cm, 14.417*cm, 14.716*cm};
		LaBrCe_MPT->AddProperty("RAYLEIGH", LaBr_Energy, LaBr_RayleighAbs, nEntries);

		//Following is a calculated quantity that most likely underestiamtes the abs length.
		//G4double LaBr_WLSABS[nEntries] = 
		//	{0.0123851*cm, 0.00283096*cm, 0.0012378*cm, 0.000979046*cm, 
		//	 0.00133068*cm, 0.00296152*cm, 0.0102811*cm, 0.0524681*cm, 
		//	 0.356516*cm, 2.69712*cm, 18.9965*cm, 128.423*cm, 963.998*cm, 
		//	 8687.01*cm, 95426.1*cm, 1.26774E6*cm, 2.00875E7*cm, 
		//	 3.74093E8*cm, 8.07423E9*cm, 1.99343E11*cm, 5.56163E12*cm};
		//LaBrCe_MPT->AddProperty("WLSABSLENGTH", LaBr_Energy, LaBr_ABSL, nEntries);

		//63 photons / keV
		LaBrCe_MPT->AddConstProperty("SCINTILLATIONYIELD",6.3/keV);
		LaBrCe_MPT->AddConstProperty("RESOLUTIONSCALE", 0.71);
		LaBrCe_MPT->AddConstProperty("FASTTIMECONSTANT", 16.*ns);

		laBr3_Ce_->SetMaterialPropertiesTable(LaBrCe_MPT);
	}

	//--------------- Glass Window ---------------------------
	glass_ = nist->FindOrBuildMaterial("G4_Pyrex_Glass");

	if (useOptical_) {
		//Create Material Properties Table to store optical properties
		G4MaterialPropertiesTable *glassMTP = new G4MaterialPropertiesTable();

		const G4int nEntries = 2;
		G4double LaBr_Energy[nEntries] = {3.82495*eV, 2.91425*eV};
		//Index of reflection
		G4double rIndex_Glass[nEntries]=
			{1.49,1.49};
		glassMTP->AddProperty("RINDEX",LaBr_Energy,rIndex_Glass,nEntries);
		//Absorbstion length
		G4double absLength_Glass[nEntries]={420.*cm,420.*cm};
		glassMTP->AddProperty("ABSLENGTH",LaBr_Energy,absLength_Glass,nEntries);

		glass_->SetMaterialPropertiesTable(glassMTP);
	}

	//--------------- PMT ------------------------------------
	//PMT cathode is made of Super Bialkali from Hamatsu R6231-100
	//Just use bialkali for now
	pmtMaterial_ = new G4Material("SuperBialkali",0,3);
	pmtMaterial_->AddElement(Sb,1);
	pmtMaterial_->AddElement(Rb,1);
	pmtMaterial_->AddElement(Cs,1);

	if (useOptical_) {
		//--------------- PMT Optical Properties -----------------
		pmtSurface_ = new G4OpticalSurface("pmtSurface_");
		pmtSurface_->SetModel(unified);
		//Must be dielectric metal for detection
		pmtSurface_->SetType(dielectric_metal);
		pmtSurface_->SetFinish(polished);

		//Define optical surface optical properties 
		G4MaterialPropertiesTable* pmtSurface__MPT = new G4MaterialPropertiesTable();
		const G4int numEntries = 22;
		G4double pmtEnergies[numEntries] = 
			{4.531417*eV, 4.214736*eV, 4.059972*eV, 3.870475*eV, 3.677380*eV, 3.502637*eV, 3.326978*eV,
			 3.168097*eV, 3.065617*eV, 2.891889*eV, 2.794437*eV, 2.692368*eV, 2.587361*eV, 2.518601*eV,
			 2.417636*eV, 2.365904*eV, 2.308279*eV, 2.193781*eV, 2.083541*eV, 1.977935*eV, 1.882518*eV,
			 1.748549*eV};

		G4double pmtReflectivity[numEntries] = {0}; 
		pmtSurface__MPT->AddProperty("REFLECTIVITY",pmtEnergies,pmtReflectivity,numEntries);

		//Efficiencncy from http://www.hamamatsu.com/us/en/technology/innovation/photocathode/index.html
		//Quantum Efficiency 
		G4double pmtEfficiency[numEntries] = 
			{0.198864, 0.284091, 0.313920, 0.339489, 0.348011, 0.349432, 0.348011, 
			0.342330, 0.335227, 0.315341, 0.299716, 0.276989, 0.237216, 0.207386, 
			0.163352, 0.140625, 0.117898, 0.075284, 0.042614, 0.019886, 0.007102,
			0.000000};
			//{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
		//We can use the efficiecny as the reflectivity is 0
		pmtSurface__MPT->AddProperty("EFFICIENCY",pmtEnergies,pmtEfficiency,numEntries);
		
		pmtSurface_->SetMaterialPropertiesTable(pmtSurface__MPT);

	}
}

void HagridCrystal::Construct(G4LogicalVolume* world, const G4ThreeVector &pos, const G4RotationMatrix &rotation/*=0*/) {

	G4double alumThickness = 0.5*mm;
	G4double alumHeight = 58.*mm;
	G4double alumRadius = 55.8/2.*mm;

	G4double crystalHeight = 50.8*mm;
	G4double crystalRadius = 50.8/2.*mm;

	G4double glassThickness = 5.*mm;
	G4double glassRadius = crystalRadius;

#ifdef DOUBLEWINDOW
	alumHeight = crystalHeight + 2. * glassThickness;
#endif

	G4double teflonHeight = alumHeight - alumThickness; 
	G4double teflonThickness = teflonHeight - crystalHeight - glassThickness;
	G4double teflonOuterRadius = alumRadius - alumThickness;

	G4double pmtThickness = 1.*cm;
	G4double pmtRadius = alumRadius;

	if (!useOptical_) pmtThickness = 0;

	//------------ Detector Volume --------------------------
	G4double detectorRadius = alumRadius;
	if (detectorRadius < pmtRadius) detectorRadius = pmtRadius;
	G4double detectorHeight = alumHeight + pmtThickness;

	G4ThreeVector detectorPos = pos + G4ThreeVector(0,0,detectorHeight/2.);

	//We compute the 3D transformation to put the object in the correct place.
	G4Transform3D transform = G4Transform3D(rotation,pos);
	//We add a translation such that the point specified is the face of the detector.
	transform = transform * G4Translate3D(G4ThreeVector(0,0,detectorHeight/2.)); 

	G4Tubs *detector = new G4Tubs(
		"detector",
		0,
		detectorRadius,
		detectorHeight/2.,
		0. * deg,
		360. * deg);
		
	G4LogicalVolume* detectorLogic = new G4LogicalVolume(
		detector,
		air_,
		"detector");					

	new G4PVPlacement(
			transform,
			detectorLogic,
			"detector",
			world,
			false,
			0,
			checkOverlaps_);

	G4VisAttributes *detectorVis = new G4VisAttributes(G4Colour(1,0,0,255/2));
	detectorVis->SetForceWireframe(true);
	detectorLogic->SetVisAttributes(detectorVis);



	//------------ Aluminum Exterior ------------------------
#ifdef POLYCONE
	const G4int numAlumIntersections = 4;
	G4double zAlumExt[numAlumIntersections] = {0, alumThickness, alumThickness, alumHeight};
	G4double rInnerAlumExt[numAlumIntersections] = {0, 0, alumRadius - alumThickness, alumRadius - alumThickness};
	G4double rOuterAlumExt[numAlumIntersections] = {alumRadius,alumRadius,alumRadius,alumRadius};

	G4ThreeVector alumExtPos=G4ThreeVector(0,0,alumThickness/2.);
	G4Polycone *alumExterior = 
		new G4Polycone("alumExterior",
				0.*deg,
				360.*deg,
				numAlumIntersections,
				zAlumExt,
				rInnerAlumExt,
				rOuterAlumExt);
#else
	G4ThreeVector alumExtPos=G4ThreeVector(0,0,(alumHeight - detectorHeight) / 2.);
	G4Tubs *alumExterior = 
		new G4Tubs("alumExterior",
				0,
				alumRadius,
				alumHeight/2.,
				0.*deg,
				360.*deg);		
#endif

	G4LogicalVolume* logicAlumExterior =
		new G4LogicalVolume(alumExterior,
				aluminum_,
				"alumExterior");						
	logicAlumExterior->SetVisAttributes(new G4VisAttributes(G4Colour(0.75,0.75,0.75)));

	new G4PVPlacement(0,
			alumExtPos,
			logicAlumExterior,
			"alumExterior",
			detectorLogic,
			false,
			0,
			checkOverlaps_);					

	//--------------- Teflon Liner ---------------------------

#ifdef POLYCONE
	G4double teflonInnerRadius = crystalRadius;
	const G4int numTeflonIntersections = 4;
	G4double zTeflonExt[numTeflonIntersections] = {0, teflonThickness, teflonThickness, teflonHeight};
	G4double rInnerTeflonExt[numTeflonIntersections] = {0, 0, teflonInnerRadius, teflonInnerRadius};
	G4double rOuterTeflonExt[numTeflonIntersections] = {teflonOuterRadius, teflonOuterRadius, teflonOuterRadius, teflonOuterRadius}; 

	//If mother volume is the world use the following
	//G4ThreeVector teflonLinerPos = alumExtPos + G4ThreeVector(0,0,alumThickness);
	G4ThreeVector teflonLinerPos = G4ThreeVector(0,0,alumThickness);
	G4Polycone *teflonLiner = 
		new G4Polycone("teflonLiner",
				0.*deg,
				360.*deg,
				numTeflonIntersections,
				zTeflonExt,
				rInnerTeflonExt,
				rOuterTeflonExt);
#else
	G4ThreeVector teflonLinerPos = G4ThreeVector(0,0,alumThickness/2.);

	G4Tubs *teflonLiner = 
		new G4Tubs("teflonLiner",
				0,
				teflonOuterRadius,
				teflonHeight/2.,
				0.*deg,
				360.*deg);
#endif

	G4LogicalVolume* logicTeflonLiner =
		new G4LogicalVolume(teflonLiner,
				teflon_,
				"teflonLiner");						

	G4PVPlacement *physTeflon = 
	new G4PVPlacement(0,
			teflonLinerPos,
			logicTeflonLiner,
			"teflonLiner",
			logicAlumExterior, //world,
			false,
			0,
			checkOverlaps_);		

	//--------------- LaBr3(Ce) Cyrstal ----------------------

	//Placement of LaBr3 Crystal
	// If mother volume is the world use the following
	//  G4ThreeVector laBr3CrystalPos= teflonLinerPos + G4ThreeVector(0,0,teflonThickness + crystalHeight/2.);
#ifdef POLYCONE
	G4ThreeVector laBr3CrystalPos=G4ThreeVector(0,0,teflonThickness + crystalHeight/2.);
#else
	G4ThreeVector laBr3CrystalPos=G4ThreeVector(0,0,teflonThickness/2.-glassThickness/2.);
#endif
	G4Tubs *laBr3Crystal = new G4Tubs(
		"LaBr3Crystal", 	//name
		0.,					//inner radius
		crystalRadius, 	//outer radius
		crystalHeight/2.,	//thickness
		0.*deg,				//starting angle
		360.*deg);			//ending angle

	G4LogicalVolume* logicLaBr3Crystal = new G4LogicalVolume(
		laBr3Crystal,
		laBr3_Ce_,
		"LaBr3Crystal");						
	logicLaBr3Crystal->SetVisAttributes(new G4VisAttributes(G4Colour(1,0.33,1)));

	G4PVPlacement *physLaBr3Crystal = 
	new G4PVPlacement(
		0,
		laBr3CrystalPos,
		logicLaBr3Crystal,
		"LaBr3Crystal",
		logicTeflonLiner,//world,
		false,
		0,
		checkOverlaps_);

	if (useOptical_) {
		//Add optical border surface between LaBr3 and teflon.
		new G4LogicalBorderSurface(
				"teflonCrystalSurface",
				physLaBr3Crystal,
				physTeflon,
				teflonSurface_);
	}

//--------------- Glass Window ---------------------------
	//If mother volume is the world use the following
	//G4ThreeVector glassFacePos = teflonLinerPos + G4ThreeVector(0,0,teflonHeight - glassThickness/2.);
#ifdef POLYCONE
	G4ThreeVector glassFacePos = G4ThreeVector(0,0,teflonHeight - glassThickness/2.);
#else
	G4ThreeVector glassFacePos = G4ThreeVector(0,0,teflonHeight/2. - glassThickness/2.);
#endif

	G4Tubs *glassFace = new G4Tubs(
		"glassFace",		//name
		0.,					//inner radius
		glassRadius,		//outer radius
		glassThickness/2.,//thickness
		0.*deg,				//starting angle
		360.*deg);			//ending angle

	G4LogicalVolume* logicGlassFace = new G4LogicalVolume(
		glassFace,
		glass_,
		"glassFace");						
	logicGlassFace->SetVisAttributes(new G4VisAttributes(G4Colour(0.33,1,1)));

	G4PVPlacement *physGlass = 
	new G4PVPlacement(
		0,
		glassFacePos,
		logicGlassFace,
		"glassFace",
		logicTeflonLiner,
		false,
		0,
		checkOverlaps_);

#ifdef DOUBLEWINDOW
	G4ThreeVector glassFacePos2 = G4ThreeVector(0,0,-teflonHeight/2. + glassThickness/2. - alumThickness);
	G4Tubs *glassFace2 = new G4Tubs(
		"glassFace2",		//name
		0.,					//inner radius
		glassRadius,		//outer radius
		glassThickness/2.,//thickness
		0.*deg,				//starting angle
		360.*deg);			//ending angle

	G4LogicalVolume* logicGlassFace2 = new G4LogicalVolume(
		glassFace2,
		glass_,
		"glassFace2");						
	logicGlassFace2->SetVisAttributes(new G4VisAttributes(G4Colour(0.33,1,1)));

	//G4PVPlacement *physGlass2 = 
	new G4PVPlacement(
		0,
		glassFacePos2,
		logicGlassFace2,
		"glassFace2",
		logicTeflonLiner,
		false,
		0,
		checkOverlaps_);


#endif
	
	//--------------- PMT ------------------------------------
	G4PVPlacement *physPMT; 
	if (useOptical_) {
			G4ThreeVector pmtPos = G4ThreeVector(0,0,(detectorHeight - pmtThickness)/ 2.);
			//G4ThreeVector pmtPos = G4ThreeVector(0,0,pmtThickness/2.);

			G4Tubs *pmt = new G4Tubs(
			"PMT",				//name
			0.,					//inner radius
			pmtRadius,			//outer radius
			pmtThickness/2.,	//thickness
			0.*deg,				//starting angle
			360.*deg);			//ending angle

		G4LogicalVolume* logicPMT = new G4LogicalVolume(
			pmt,
			pmtMaterial_,
			"PMT");						

		physPMT = 
			new G4PVPlacement(
				0,
				pmtPos,
				logicPMT,
				"PMT",
				detectorLogic,
				false,
				0,
				checkOverlaps_);

	}

	if (useOptical_) {
		//Create the optical surface as the directional border from first
		//	physical volume to the second physical volume
		new G4LogicalBorderSurface("opSurfGlassPMT", physGlass, physPMT, pmtSurface_);
	}

}

