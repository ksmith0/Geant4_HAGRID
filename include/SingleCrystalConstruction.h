/**
 *
 */

#ifndef SINGLECRYSTALCONSTRUCTION_H
#define SINGLECRYSTALCONSTRUCTION_H

#include "G4VUserDetectorConstruction.hh"

class SingleCrystalConstruction : public G4VUserDetectorConstruction
{
	private:
		G4bool fUseOptical;
		
	public:
		SingleCrystalConstruction(G4bool useOptical=true);
		virtual ~SingleCrystalConstruction();
		
		virtual G4VPhysicalVolume* Construct();
		virtual void ConstructSDandField();
	
};

#endif
