#include "../hierarchy.h"

void MultiResolutionHierarchy::convert2Poly()
{
	mPoly = new HE_Polyhedron(mV_tag, F_tag);
}