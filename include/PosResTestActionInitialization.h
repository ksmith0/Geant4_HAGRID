#ifndef POSRESTESTACTIONINITIALIZATION
#define POSRESTESTACTIONINITIALIZATION

#include "G4VUserActionInitialization.hh"

class PosResTestActionInitialization : public G4VUserActionInitialization
{
	public:
		PosResTestActionInitialization();
		virtual ~PosResTestActionInitialization();

		virtual void BuildForMaster() const;
		virtual void Build() const;
};

#endif
