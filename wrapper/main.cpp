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
	char *btype = NULL;
	if (!ftype.empty()) btype = (char *) ftype.c_str();
	SphericalRemeshing *SR;

	//if (!input.empty()) cout << "subject: " << subject << endl;
	//cout << "sphere: " << sphere << endl;
	//if (!coeff.empty())  cout << "deformation: " << dfield << endl;
	//if (!ref.empty()) cout << "reference: " << reference << endl;
	
	SR = new SphericalRemeshing(subject, sphere_, dfield, keepC, reference, colormap, property, !nneighbor, deg, !quite, btype);
	deg = SR->degree();
	if (deg0 == -1 || deg0 > deg) deg0 = deg;
	for (int i = deg0; i <= deg; i++)
	{
		SR->deform(i);

		if (!quite) cout << "Write output surface model..\n";
		if (!output.empty())
		{
			std::string outname = output;
			if (deg0 < deg) outname += to_string(i);
			if (!quite) cout << outname.c_str() << endl;
			SR->saveDeformedSurface(outname.c_str());
		}
		if (dsphere != NULL)
		{
			std::string outname = dsphere;
			if (deg0 < deg) outname += to_string(i);
			if (!quite) cout << outname << endl;
			SR->saveDeformedSphere(outname.c_str());
		}
		if (!outputProp.empty())
		{
			std::string outname = outputProp;
			if (deg0 < deg) outname += to_string(i);
			if (!quite) cout << outname.c_str() << endl;
			SR->saveDeformedProperty(outname.c_str());
		}
		if (!bary.empty())
		{
			std::string outname = bary;
			if (deg0 < deg) outname += to_string(i);
			if (!quite) cout << outname.c_str() << endl;
			SR->saveBary(outname.c_str());
		}
	}
	delete SR;
	
	return 0;
}

