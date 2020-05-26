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
bool quite = false;
std::string bary;
std::vector<std::string> property;

std::string ftype = "";

void PARSE_ARGS(int argc, char **argv)
{
    
	std::string desc("Spherical Remeshing Tool "
					 SREMESH_VERSION "\n"
					 "Author: Ilwoo Lyu\n"
					 );

	CLI::App app(desc);
	app.add_flag("--quite", quite, "Do not print outputs");

	app.add_option("-s,--sphere", sphere, "Specify spherical model (unit sphere) used for the original spherical parameterization (vtk)")->required()->check(CLI::ExistingFile)->group("Inputs");
	app.add_option("-i,--input", input, "Specify input subject surface model (vtk)")->check(CLI::ExistingFile)->group("Inputs");
	app.add_option("-c,--coeff", coeff, "Specify HSD spherical harmonic coefficients (txt)")->check(CLI::ExistingFile)->group("Inputs");
	app.add_option("-r,--ref", refSphere, "Specify reference model (unit sphere) for the final surface remeshing (vtk)")->check(CLI::ExistingFile)->group("Inputs");
	app.add_option("-p,--property", property, "Specify properties that are included in the resulting mesh files (txt)")->check(CLI::ExistingFile)->group("Inputs");
	app.add_option("--color", color, "Specify color map of the reference model that will assign the same color to the generated surface (txt)")->check(CLI::ExistingFile)->group("Inputs");

	app.add_option("-o,--output", output, "Specify output surface model (vtk)")->group("Outputs");
	app.add_option("--outputProperty", outputProp, "Specify output of properties after remeshing (no extension needed)")->group("Outputs");
	app.add_option("--deform", deform, "Specify deformed unit sphere by HSD deformation (vtk)")->group("Outputs");
	app.add_option("--bary", bary, "Specify barycentric coordinates (txt)")->group("Outputs");
	app.add_option("--ftype", ftype, "Specify data type - output will be saved as a binary format (dat)")->check(CLI::IsMember({"int", "float","double"}, CLI::ignore_case, CLI::ignore_underscore))->group("Outputs");

	app.add_flag("--nneighbor", nneighbor, "Enable nearest neighbor interpolation")->group("Remeshing parameters");
	app.add_option("--deg", deg, "Specify level of HSD spherical harmonics")->check(CLI::NonNegativeNumber)->group("Remeshing parameters");
	app.add_flag("--keepColor", keepC, "Keep color of the template model")->group("Remeshing parameters");

	try
	{
		app.parse(argc, argv);
	}
	catch (const CLI::ParseError &e)
	{
		exit(app.exit(e));
	}
}
