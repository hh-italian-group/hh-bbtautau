// BTagEfficiency.cxx
#include <boost/format.hpp>
#include <vector>

#include "AnalysisTools/Run/include/program_main.h" // definition of wrappers for the program main and program arguments.
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"


struct Arguments { // list of all program arguments
    REQ_ARG(std::vector<std::string>, input_file); // required argument "input_file"
    REQ_ARG(std::string, output_file); // required argument "output_file"
    OPT_ARG(bool, flag, false); // optional argument "flag" with the default value = false
};

class BTagData : public root_ext::AnalyzerData {
    public:
	explicit BTagData(std::shared_ptr<TFile> _outputFile, const std::string& directoryName = "") :
	    AnalyzerData(_outputFile, directoryName)
	    {
              
	    }
            TH1D_ENTRY(h1,201,-0.5,200.5)
};

class BTagEfficiency { // simple analyzer definition
    public:
	using Event = ntuple::Event;
	using EventTuple = ntuple::EventTuple;


	BTagEfficiency(const Arguments& _args) : args(_args), outfile(root_ext::CreateRootFile(args.output_file())), anaData(outfile)
    {
	// Analyzer initialization (e.g. open input/output files, parse configs...)
    }
	void Run()
	{
	    // analyzer code
	    //std::cout << boost::format("Processing input file '%1%' into output file '%2%' with flag = %3%.\n")
	    //	% args.input_file() % args.output_file() % args.flag();
	    for (const auto& name : args.input_file()){
		std::shared_ptr<TFile> in_file(root_ext::OpenRootFile(name));
		//ntuple::EventTuple myTree("muTau", in_file.get(), true);
		std::shared_ptr<EventTuple> tuple;
		try {
		    tuple = ntuple::CreateEventTuple("muTau", in_file.get(), true, ntuple::TreeState::Full,true);
		} catch(std::exception&) {
		    std::cerr << "WARNING: tree 'mutau not found in file '"<< std::endl;
		}

		for(const Event& event : *tuple){

		    auto event_ptr = std::make_shared<Event>(event);
		    
                    for (auto jet : event_ptr->jets_p4){
                        //std::cout<<" Jet Pt = "<<jet.Pt()<<std::endl;
                        anaData.h1("Jet Pt").Fill(jet.Pt());

                    }
                    
             
		}
	    }
	}
    private:
	Arguments args;
	std::shared_ptr<TFile> outfile;
	BTagData anaData;


};


PROGRAM_MAIN(BTagEfficiency, Arguments) // definition of the main program function
