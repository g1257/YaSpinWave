#include "../../spf/v7/Engine/Engine.h"
#include "InputCheckMonteCarlo.h"
#include "../../spf/v7/Engine/ParametersEngine.h"
#include "InputNg.h"
#include "MersenneTwister.h"
#include "Heisenberg.h"
#include "SpinModulus.h"

template<typename InputType>
class BogusGeometry {

public:

	BogusGeometry(InputType& io)
	{
		io.readline(jfile,"Jinput=");
		std::ifstream fin(jfile.c_str());
		if (!fin || fin.bad() || !fin.good()) {
			PsimagLite::String msg("BogusGeometry: Cannot open file ");
			msg += jfile + "\n";
			throw PsimagLite::RuntimeError(msg);
		}

		int tmp = 0;
		fin>>tmp;
		if (fin.eof()) {
			PsimagLite::String msg("BogusGeometry: Syntax error in ");
			msg += jfile + "\n";
			throw PsimagLite::RuntimeError(msg);
		}

		fin>>unitCellSize;
		if (fin.eof()) {
			PsimagLite::String msg("BogusGeometry: Syntax error in ");
			msg += jfile + "\n";
			throw PsimagLite::RuntimeError(msg);
		}

		io.readline(tmp,"Verbose=");
		verbose = (tmp > 0);

		io.read(qvector,"Qvector ");
	}

	int unitCellSize;
	PsimagLite::String jfile;
	bool verbose;
	PsimagLite::Vector<double>::Type qvector;
};

typedef double RealType;
typedef PsimagLite::InputNg<yasw::InputCheckMonteCarlo> InputNgType;
typedef BogusGeometry<InputNgType::Readable> GeometryType;
typedef Spf::ParametersEngine<RealType,InputNgType::Readable> ParametersEngineType;
typedef yasw::Heisenberg<ParametersEngineType,GeometryType> ModelType;
typedef ModelType::SpaceConnectorsType SpaceConnectorsType;
typedef PsimagLite::MersenneTwister RngType;
typedef Spf::Engine<ParametersEngineType,
                    ModelType,
                    InputNgType::Readable,
                    RngType> EngineType;
typedef yasw::SpinModulus<PsimagLite::Vector<RealType>::Type> SpinModulusType;

int main(int argc, char* argv[])
{
	PsimagLite::String filename="";
	SizeType pixelSize = 1;

	int opt = 0;
	PsimagLite::String strUsage(argv[0]);
	strUsage += " -f filename [-P pixelSize] [-M modulusFile]";
	yasw::InputCheckMonteCarlo inputCheck;
	PsimagLite::String spinModulusFile;

	while ((opt = getopt(argc, argv,"f:P:M:")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		case 'P':
			pixelSize = atoi(optarg);
			break;
		case 'M':
			spinModulusFile = optarg;
			break;
		default:
			inputCheck.usageMain(strUsage);
			return 1;
		}
	}

	// sanity checks here
	if (filename=="") {
		inputCheck.usageMain(strUsage);
		return 1;
	}

	InputNgType::Writeable ioWriteable(filename,inputCheck);
	InputNgType::Readable io(ioWriteable);

	ParametersEngineType engineParams(io);
	std::cerr<<engineParams;

	GeometryType geometry(io);
	SpaceConnectorsType spaceConnectors(geometry.jfile, pixelSize, geometry.verbose);
	SpinModulusType spinModulus(spinModulusFile, spaceConnectors.rows());
	ModelType model(engineParams, geometry, spaceConnectors, spinModulus());

	EngineType engine(engineParams, model, io);

	engine.main();
}

