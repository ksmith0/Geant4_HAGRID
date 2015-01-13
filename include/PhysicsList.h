#ifndef PHYSICSLIST_H
#define PHYSICSLIST_H

#include "G4VUserPhysicsList.hh"

class PhysicsList : public G4VUserPhysicsList {
	private:

	public:
		PhysicsList();
		virtual ~PhysicsList() {};
		void ConstructProcess();
		void ConstructParticle();
		void ConstructEM();
		void ConstructOp();

};

#endif

