#include "../../spf/v7/Engine/Engine.h"
#include "InputCheckMonteCarlo.h"
#include "../../spf/v7/Engine/ParametersEngine.h"
#include "InputNg.h"
#include "MersenneTwister.h"
#include "Heisenberg.h"
#include "Qvectors.h"

template<typename InputType, typename RealType>
class BogusGeometry {

public:

	typedef PsimagLite::Matrix<RealType> MatrixRealType;
	typedef yasw::Qvectors<MatrixRealType> QvectorsType;

	BogusGeometry(InputType& io) : qvectors(io)
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
	}

	int unitCellSize;
	PsimagLite::String jfile;
	bool verbose;
	QvectorsType qvectors;
};

typedef double RealType;
typedef PsimagLite::InputNg<yasw::InputCheckMonteCarlo> InputNgType;
typedef BogusGeometry<InputNgType::Readable, RealType> GeometryType;
typedef Spf::ParametersEngine<RealType,InputNgType::Readable> ParametersEngineType;
typedef yasw::Heisenberg<ParametersEngineType,GeometryType> ModelType;
typedef PsimagLite::MersenneTwister RngType;
typedef Spf::Engine<ParametersEngineType,
                    ModelType,
                    InputNgType::Readable,
                    RngType> EngineType;

int main(int argc, char* argv[])
{
	PsimagLite::String filename="";
	SizeType pixelSize = 1;

	int opt = 0;
	PsimagLite::String strUsage(argv[0]);
	strUsage += " -f filename [-P pixelSize]";
	yasw::InputCheckMonteCarlo inputCheck;

	while ((opt = getopt(argc, argv,"f:P:")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		case 'P':
			pixelSize = atoi(optarg);
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
	ModelType model(engineParams, geometry, pixelSize, io);

	EngineType engine(engineParams, model, io);

	engine.main();
}

