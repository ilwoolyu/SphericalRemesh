#include <string>

#include "SRemeshCLP.h"
#include "SphericalRemeshing.h"

int main(int argc, char* argv[])
{
	PARSE_ARGS;
	
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " --help" << std::endl;
		return -1;
	}
	
	char *subject = NULL;
	if (!input.empty()) subject = (char *) input.c_str();
	char *dfield = NULL;
	if (!coeff.empty()) dfield = (char *) coeff.c_str();
	char *sphere_ = (char *) sphere.c_str();
	char *reference = NULL;
	if (!ref.empty()) reference = (char *) ref.c_str();
	char *dsphere = NULL;
	if (!deform.empty()) dsphere = (char *) deform.c_str();
	char *colormap = NULL;
	if (!color.empty()) colormap = (char *) color.c_str();
	SphericalRemeshing *SR;

	//if (!input.empty()) cout << "subject: " << subject << endl;
	//cout << "sphere: " << sphere << endl;
	//if (!coeff.empty())  cout << "deformation: " << dfield << endl;
	//if (!ref.empty()) cout << "reference: " << reference << endl;
	
	SR = new SphericalRemeshing(subject, sphere_, dfield, keepC, reference, colormap, property, !nneighbor, deg);

	cout << "Write output surface model..\n";
	if (!output.empty())
	{
		cout << output.c_str() << endl;
		SR->saveDeformedSurface(output.c_str());
	}
	if (dsphere != NULL)
	{
		cout << dsphere << endl;
		SR->saveDeformedSphere(dsphere);
	}
	if (!outputProp.empty())
	{
		cout << outputProp.c_str() << endl;
		SR->saveDeformedProperty(outputProp.c_str());
	}

	delete SR;
	
	return 0;
}

