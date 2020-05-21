#include "SphericalRemeshing.h"
#include "PARSE_ARGS.h"

int main(int argc, char **argv)
{
    PARSE_ARGS(argc, argv);

	char *subject = NULL;
	if (!input.empty()) subject = (char *) input.c_str();
	char *dfield = NULL;
	if (!coeff.empty()) dfield = (char *) coeff.c_str();
	char *sphere_ = (char *) sphere.c_str();
	char *reference = NULL;
	if (!refSphere.empty()) reference = (char *) refSphere.c_str();
	char *dsphere = NULL;
	if (!deform.empty()) dsphere = (char *) deform.c_str();
	char *colormap = NULL;
	if (!color.empty()) colormap = (char *) color.c_str();
	SphericalRemeshing *SR;

	//if (!input.empty()) cout << "subject: " << subject << endl;
	//cout << "sphere: " << sphere << endl;
	//if (!coeff.empty())  cout << "deformation: " << dfield << endl;
	//if (!ref.empty()) cout << "reference: " << reference << endl;
	
	SR = new SphericalRemeshing(subject, sphere_, dfield, keepC, reference, colormap, property, !nneighbor, deg, !quite);

	if (!quite) cout << "Write output surface model..\n";
	if (!output.empty())
	{
		if (!quite) cout << output.c_str() << endl;
		SR->saveDeformedSurface(output.c_str());
	}
	if (dsphere != NULL)
	{
		if (!quite) cout << dsphere << endl;
		SR->saveDeformedSphere(dsphere);
	}
	if (!outputProp.empty())
	{
		if (!quite) cout << outputProp.c_str() << endl;
		SR->saveDeformedProperty(outputProp.c_str());
	}
	if (!bary.empty())
	{
		if (!quite) cout << bary.c_str() << endl;
		SR->saveBary(bary.c_str());
	}

	delete SR;
	
	return 0;
}

