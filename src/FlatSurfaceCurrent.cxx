
// FlatSurfaceCurrent
#include "FlatSurfaceCurrent.h"

#include "G4SystemOfUnits.hh"
#include "G4StepStatus.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4UnitsTable.hh"
#include "G4GeometryTolerance.hh"
////////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring only Surface Flux.
//  Current version assumes only for G4Box shape. 
//
// Surface is defined at the surface specified.
// Assuming -Z face we have the following.
// Direction                  -Z   +Z
//   0  IN || OUT            ->|<-  |
//   1  IN                   ->|    |
//   2  OUT                    |<-  |
//
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
// 17-Nov-2005 T.Aso, Bug fix for area definition.
// 31-Mar-2007 T.Aso, Add option for normalizing by the area.
// 2010-07-22   Introduce Unit specification.
// 2014-11-25  K.Smith, Introduce faces
// 
///////////////////////////////////////////////////////////////////////////////


FlatSurfaceCurrent::FlatSurfaceCurrent(G4String name, 
					 G4int direction, G4int face, G4int depth)
    :G4VPrimitiveScorer(name,depth),HCID(-1),fDirection(direction),fFace(face),
     weighted(true),divideByArea(true)
{
    DefineUnitAndCategory();
    SetUnit("percm2");
}

FlatSurfaceCurrent::FlatSurfaceCurrent(G4String name, 
					       G4int direction, G4int face,
					       const G4String& unit, 
					       G4int depth)
    :G4VPrimitiveScorer(name,depth),HCID(-1),fDirection(direction),fFace(face),
     weighted(true),divideByArea(true)
{
    DefineUnitAndCategory();
    SetUnit(unit);
}

FlatSurfaceCurrent::~FlatSurfaceCurrent()
{;}

G4bool FlatSurfaceCurrent::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4StepPoint* preStep = aStep->GetPreStepPoint();
  G4VPhysicalVolume* physVol = preStep->GetPhysicalVolume();
  G4VPVParameterisation* physParam = physVol->GetParameterisation();
  G4VSolid * solid = 0;
  if(physParam)
  { // for parameterized volume
    G4int idx = ((G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable()))
                ->GetReplicaNumber(indexDepth);
    solid = physParam->ComputeSolid(idx, physVol);
    solid->ComputeDimensions(physParam,idx,physVol);
  }
  else
  { // for ordinary volume
    solid = physVol->GetLogicalVolume()->GetSolid();
  }

  G4Box* boxSolid = (G4Box*)(solid);

  G4int dirFlag =IsSelectedSurface(aStep,boxSolid);
  if ( dirFlag > 0 ) {
   if ( fDirection == fCurrent_InOut || fDirection == dirFlag ){
       G4int index = GetIndex(aStep);
       G4TouchableHandle theTouchable = preStep->GetTouchableHandle();
       G4double current = 1.0;
       if ( weighted ) current=preStep->GetWeight(); // Current (Particle Weight)
       if ( divideByArea ){
	   G4double square = 1;
		if (fFace == fXFace || fFace == fXMinusFace) square = 4.*boxSolid->GetYHalfLength()*boxSolid->GetZHalfLength();
		else if (fFace == fYFace || fFace == fYMinusFace) square = 4.*boxSolid->GetXHalfLength()*boxSolid->GetZHalfLength();
		else if (fFace == fZFace || fFace == fZMinusFace) square = 4.*boxSolid->GetXHalfLength()*boxSolid->GetYHalfLength();
	   current = current/square;  // Normalized by Area
       }
       EvtMap->add(index,current);
   }
  }

  return TRUE;
}

G4int FlatSurfaceCurrent::IsSelectedSurface(G4Step* aStep, G4Box* boxSolid){

  G4TouchableHandle theTouchable = 
    aStep->GetPreStepPoint()->GetTouchableHandle();
  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
	G4double halfLength = 0;
	if (fFace == fXFace) halfLength = boxSolid->GetXHalfLength();
	else if (fFace == fXMinusFace) halfLength = -boxSolid->GetXHalfLength();
	else if (fFace == fYFace) halfLength = boxSolid->GetYHalfLength();
	else if (fFace == fYMinusFace) halfLength = -boxSolid->GetYHalfLength();
	else if (fFace == fZFace) halfLength = boxSolid->GetZHalfLength();
	else if (fFace == fZMinusFace) halfLength = -boxSolid->GetZHalfLength();

  if (aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary ){
    // Entering Geometry
    G4ThreeVector stppos1= aStep->GetPreStepPoint()->GetPosition();
    G4ThreeVector localpos1 = 
      theTouchable->GetHistory()->GetTopTransform().TransformPoint(stppos1);
    if((fFace == fXFace || fFace == fXMinusFace) && std::fabs(localpos1.x() + halfLength)<kCarTolerance )
      return fCurrent_In;
    else if((fFace == fYFace || fFace == fYMinusFace) && std::fabs(localpos1.y() + halfLength)<kCarTolerance )
      return fCurrent_In;
    else if((fFace == fZFace || fFace == fZMinusFace) && std::fabs(localpos1.z() + halfLength)<kCarTolerance )
      return fCurrent_In;
  }

  if (aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary ){
    // Exiting Geometry
    G4ThreeVector stppos2= aStep->GetPostStepPoint()->GetPosition();
    G4ThreeVector localpos2 = 
      theTouchable->GetHistory()->GetTopTransform().TransformPoint(stppos2);
	 if((fFace == fXFace || fFace == fXMinusFace) && std::fabs(localpos2.x() + halfLength)<kCarTolerance )
		 return fCurrent_Out;
	 else if((fFace == fYFace || fFace == fYMinusFace) && std::fabs(localpos2.y() + halfLength)<kCarTolerance )
		 return fCurrent_Out;
	 else if((fFace == fZFace || fFace == fZMinusFace) && std::fabs(localpos2.z() + halfLength)<kCarTolerance )
		 return fCurrent_Out;
  }

  return -1;
}

void FlatSurfaceCurrent::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(detector->GetName(), GetName());
  if ( HCID < 0 ) HCID = GetCollectionID(0);
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void FlatSurfaceCurrent::EndOfEvent(G4HCofThisEvent*)
{;}

void FlatSurfaceCurrent::clear(){
  EvtMap->clear();
}

void FlatSurfaceCurrent::DrawAll()
{;}

void FlatSurfaceCurrent::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() <<G4endl; 
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for(; itr != EvtMap->GetMap()->end(); itr++) {
      G4cout << "  copy no.: " << itr->first << " current : " ;
      if ( divideByArea ) {
	  G4cout << *(itr->second)/GetUnitValue() 
		 << " ["<<GetUnit()<<"]";
      }else {
	  G4cout <<  *(itr->second)/GetUnitValue() << " [tracks]";
      }
      G4cout << G4endl;
  }
}

void FlatSurfaceCurrent::SetUnit(const G4String& unit)
{
    if ( divideByArea ) {
	CheckAndSetUnit(unit,"Per Unit Surface");
    } else {
	if (unit == "" ){
	    unitName = unit;
	    unitValue = 1.0;
	}else{
	    G4String msg = "Invalid unit ["+unit+"] (Current  unit is [" +GetUnit()+"] ) for " + GetName();
	    G4Exception("FlatSurfaceCurrent::SetUnit","DetPS0007",JustWarning,msg);
	}
    }
}

void FlatSurfaceCurrent::DefineUnitAndCategory(){
   // Per Unit Surface
   new G4UnitDefinition("percentimeter2","percm2","Per Unit Surface",(1./cm2));
   new G4UnitDefinition("permillimeter2","permm2","Per Unit Surface",(1./mm2));
   new G4UnitDefinition("permeter2","perm2","Per Unit Surface",(1./m2));
}



