#ifndef FLATSURFACECURRENT_H
#define FLATSURFACECURRENT_H

#include "G4VPrimitiveScorer.hh"
#include "G4THitsMap.hh"

#include "G4Box.hh"
#include "G4PSDirectionFlag.hh"

class FlatSurfaceCurrent : public G4VPrimitiveScorer
{

	public: // with description
		enum FaceFlag {
			fXFace,
			fXMinusFace,
			fYFace,
			fYMinusFace,
			fZFace,
			fZMinusFace};

		FlatSurfaceCurrent(G4String name, G4int direction, G4int face, G4int depth=0);
      FlatSurfaceCurrent(G4String name, G4int direction, G4int face,
			     const G4String& unit, G4int depth=0);
      virtual ~FlatSurfaceCurrent();

      // Scoring options
      inline void Weighted(G4bool flg=true) { weighted = flg; }
      // Multiply track weight
      inline void DivideByArea(G4bool flg=true) { divideByArea = flg; }
      // Divided By Area

  protected: // with description
      virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      G4int IsSelectedSurface(G4Step*,G4Box*);

  public: 
      virtual void Initialize(G4HCofThisEvent*);
      virtual void EndOfEvent(G4HCofThisEvent*);
      virtual void clear();
      virtual void DrawAll();
      virtual void PrintAll();

      virtual void SetUnit(const G4String& unit);    

  protected:
      virtual void DefineUnitAndCategory();

  private:
      G4int  HCID;
      G4int  fDirection;
		G4int  fFace;
      G4THitsMap<G4double>* EvtMap;
      G4bool weighted;
      G4bool divideByArea;
};

#endif
