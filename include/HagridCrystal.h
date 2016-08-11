#ifndef HAGRIDCRYSTAL_H
#define HAGRIDCRYSTAL_H

#include "G4Material.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#include "G4LogicalVolume.hh"

#include "G4OpticalSurface.hh"

class HagridCrystal {
	public:
		HagridCrystal();
		virtual ~HagridCrystal();
		void BuildMaterials();
		void Construct(G4LogicalVolume* world, const G4ThreeVector &pos, const G4RotationMatrix &rotation = G4RotationMatrix(0,0,0));

	private:
		bool useOptical_;
		bool checkOverlaps_;
		G4Material *air_;
		G4Material *teflon_;
		G4Material *aluminum_;
		G4Material *glass_;
		G4Material *laBr3_Ce_;
		G4Material *pmtMaterial_;
		G4OpticalSurface *teflonSurface_;
		G4OpticalSurface* pmtSurface_;
};

#endif //HAGRIDCRYSTAL_H
