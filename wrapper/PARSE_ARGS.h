#include <string>
#include <vector>
#include "CLI11.hpp"

bool keepC = false;
std::string sphere;
std::string input;
std::string coeff;
int deg = -1;
std::string output;
std::string outputProp;
std::string refSphere;
std::string deform;
std::string color;
bool nneighbor = false;
std::string bary;
std::vector<std::string> property;

void PARSE_ARGS(int argc, char **argv)
{
    
	std::string desc("Spherical Remeshing Tool v1.1\n"
					 "Author: Ilwoo Lyu\n"
					 );

	CLI::App app(desc);

	app.add_flag("--keepColor", keepC, "Keep color of the template model");
	app.add_option("-s,--sphere", sphere, "Provide spherical model (unit sphere) used for the original spherical parameterization (vtk)")->required()->check(CLI::ExistingFile);
	app.add_option("-i,--input", input, "Provide input subject surface model (vtk)")->check(CLI::ExistingFile);
	app.add_option("-c,--coeff", coeff, "Provide spherical harmonic coefficients (txt)")->check(CLI::ExistingFile);
	app.add_option("--deg", deg, "Provide Reconstruction level of spherical harmonicss")->check(CLI::NonNegativeNumber);
	app.add_option("-o,--output", output, "Remeshe subject surface model (vtk)");
	app.add_option("--outputProperty", outputProp, "Write output of properties after remeshing (no extension needed)");
	app.add_option("-r,--ref", refSphere, "Provide reference model (unit sphere) for the final surface remeshing (vtk)")->check(CLI::ExistingFile);
	app.add_option("--deform", deform, "Deform unit sphere by the given deformation field (vtk)");
	app.add_option("--color", color, "Provide color map of the reference model that will assign the same color to the generated surface (txt)");
	app.add_flag("--nneighbor", nneighbor, "Enable nearest neighbor interpolation");
	app.add_option("--bary", bary, "Write barycentric coordinates (txt)");
	app.add_option("-p,--property", property, "Provide properties that are included in the resulting mesh file (txt)")->check(CLI::ExistingFile);

	try
	{
		app.parse(argc, argv);
	}
	catch (const CLI::ParseError &e)
	{
		exit(app.exit(e));
	}
}
