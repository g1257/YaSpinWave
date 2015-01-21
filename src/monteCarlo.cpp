#include "../../spf/v7/Engine/Engine.h"
#include "InputCheckMonteCarlo.h"
#include "../../spf/v7/Engine/ParametersEngine.h"
#include "InputNg.h"
#include "MersenneTwister.h"
#include "Heisenberg.h"

template<typename InputType>
class BogusGeometry {

public:

	BogusGeometry(InputType& io)
	{
		io.readline(unitCellSize,"UnitCellSize=");
		io.readline(jfile,"Jinput=");
		int tmp = 0;
		io.readline(tmp,"Verbose=");
		verbose = (tmp > 0);
	}

	int unitCellSize;
	PsimagLite::String jfile;
	bool verbose;
};

typedef double RealType;
typedef PsimagLite::InputNg<yasw::InputCheckMonteCarlo> InputNgType;
typedef BogusGeometry<InputNgType::Readable> GeometryType;
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
	int opt = 0;
	PsimagLite::String strUsage(argv[0]);
	strUsage += " -f filename";
	yasw::InputCheckMonteCarlo inputCheck;

	while ((opt = getopt(argc, argv,"f:")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
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
	ModelType model(engineParams,geometry,io);

	EngineType engine(engineParams,model,io);

	engine.main();
}

