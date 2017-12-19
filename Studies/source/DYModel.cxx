// simple_analyzer.cxx
#include <boost/format.hpp>
#include "AnalysisTools/Run/include/program_main.h" // definition of wrappers for the program main and program arguments.
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
namespace analysis{
class DYModel { // simple analyzer definition
public:
    std::shared_ptr<TH1F> scale_factor_histo;
    std::map<std::string,double> scale_factor_maps;
    DYModel()
    {
        // Analyzer initialization (e.g. open input/output files, parse configs...)
        std::shared_ptr<TFile> input_file = root_ext::OpenRootFile("hh-bbtautau/McCorrections/data/DY_Scale_factors.root");
        scale_factor_histo = std::make_shared<TH1F>(root_ext::ReadObject<TH1F>(*input_file,"scale_factors"));
        int nbins = scale_factor_histo->GetNbinsX();
        for(int i=1; i<=nbins;i++){
            std::string scale_factor_name = scale_factor_histo->GetXaxis()->GetBinLabel(i);
            double value = scale_factor_histo->GetBinContent(i);
            scale_factor_maps[scale_factor_name] = value;
        }
    }

private:

};
}
PROGRAM_MAIN(analysis::DYModel,analysis::AnalyzerArguments) // definition of the main program function
