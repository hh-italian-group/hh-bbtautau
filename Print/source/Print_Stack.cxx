/*! Print stack with specified name superimposing several files.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "hh-bbtautau/Analysis/include/SampleDescriptorConfigEntryReader.h"
#include "hh-bbtautau/Analysis/include/StackedPlotsProducer.h"

struct Arguments {
    REQ_ARG(std::string, setup);
    REQ_ARG(std::string, sources);
    REQ_ARG(std::string, input);
    REQ_ARG(std::string, output);
    REQ_ARG(analysis::Channel, channel);
    REQ_ARG(std::string, hists);
};

namespace analysis {

class PrintStackBase {
public:
    virtual ~PrintStackBase() {}
    virtual void Run() = 0;
};

template<typename _FirstLeg, typename _SecondLeg>
class PrintStackT : public PrintStackBase {
public:
    using FirstLeg = _FirstLeg;
    using SecondLeg = _SecondLeg;
    using AnaData = ::analysis::EventAnalyzerData<FirstLeg, SecondLeg>;
    using AnaDataCollection = ::analysis::EventAnalyzerDataCollection<AnaData>;
    using PlotsProducer = ::analysis::StackedPlotsProducer<FirstLeg, SecondLeg>;

    static const EventCategorySet& EventCategoriesToProcess()
    {
        static const EventCategorySet categories = {
            EventCategory::TwoJets_Inclusive(), EventCategory::TwoJets_ZeroBtag_Resolved(),
            EventCategory::TwoJets_OneBtag_Resolved(), EventCategory::TwoJets_TwoBtag_Resolved(),
            EventCategory::TwoJets_TwoLooseBtag_Boosted()
        };
        return categories;
    }

    PrintStackT(const Arguments& _args) :
        args(_args), anaDataCollection(args.input())
    {
        ConfigReader config_reader;

        AnalyzerSetupCollection ana_setup_collection;
        AnalyzerConfigEntryReader ana_entry_reader(ana_setup_collection);
        config_reader.AddEntryReader("ANA_DESC", ana_entry_reader, false);

        MvaReaderSetupCollection mva_setup_collection;
        MvaReaderSetupEntryReader mva_entry_reader(mva_setup_collection);
        config_reader.AddEntryReader("MVA", mva_entry_reader, false);

        SampleDescriptorConfigEntryReader sample_entry_reader(sample_descriptors);
        config_reader.AddEntryReader("SAMPLE", sample_entry_reader, true);

        CombinedSampleDescriptorConfigEntryReader combined_entry_reader(cmb_sample_descriptors, sample_descriptors);
        config_reader.AddEntryReader("SAMPLE_CMB", combined_entry_reader, false);

        config_reader.ReadConfig(args.sources());

        if(!ana_setup_collection.count(args.setup()))
            throw exception("Setup '%1%' not found in the configuration file '%2%'.") % args.setup() % args.sources();
        ana_setup = ana_setup_collection.at(args.setup());

        if(ana_setup.mva_setup.size()) {
            if(!mva_setup_collection.count(ana_setup.mva_setup))
                throw exception("MVA setup '%1%' not found.") % ana_setup.mva_setup;
            mva_setup = mva_setup_collection.at(ana_setup.mva_setup);
        }

        CreateEventSubCategoriesToProcess();
        const std::vector<std::string> hist_list = SplitValueList(args.hists(), false, ", ", true);
        hist_names.insert(hist_list.begin(), hist_list.end());
    }

    virtual void Run() override
    {
        std::cout << "Creating plots..." << std::endl;
        const auto samplesToDraw = PlotsProducer::CreateOrderedSampleCollection(
                    ana_setup.draw_sequence, sample_descriptors, cmb_sample_descriptors, ana_setup.signals,
                    ana_setup.data);
        PlotsProducer plotsProducer(anaDataCollection, samplesToDraw, false, true, hist_names);
        const std::set<std::string> signal_names(ana_setup.signals.begin(), ana_setup.signals.end());
        plotsProducer.PrintStackedPlots(args.output(), EventRegion::SignalRegion(), EventCategoriesToProcess(),
                                        sub_categories_to_process, signal_names);
    }

private:
    void CreateEventSubCategoriesToProcess()
    {
        const auto base = EventSubCategory().SetCutResult(SelectionCut::mh, true);
        sub_categories_to_process.insert(base);
        if(mva_setup.is_initialized()) {
            for(const auto& mva_sel : mva_setup->selections) {
                auto sub_category = EventSubCategory(base).SetCutResult(mva_sel.first, true);
                sub_categories_to_process.insert(sub_category);
            }
        }
    }

private:
    Arguments args;
    AnaDataCollection anaDataCollection;
    AnalyzerSetup ana_setup;
    boost::optional<MvaReaderSetup> mva_setup;
    SampleDescriptorCollection sample_descriptors;
    CombinedSampleDescriptorCollection cmb_sample_descriptors;
    EventSubCategorySet sub_categories_to_process;
    std::set<std::string> hist_names;
};

class PrintStack {
public:
    PrintStack(const Arguments& _args) : print_base(MakePrintStack(_args)) {}

    void Run()
    {
        print_base->Run();
    }

private:
    template<int channel_id>
    static std::shared_ptr<PrintStackBase> MakePrintStackForChannel(const Arguments& args)
    {
        using LegInfo = ChannelLegInfo<channel_id>;
        using FirstLeg = typename LegInfo::FirstLeg;
        using SecondLeg = typename LegInfo::SecondLeg;
        using PrintStackClass = PrintStackT<FirstLeg, SecondLeg>;
        return std::make_shared<PrintStackClass>(args);
    }

    static std::shared_ptr<PrintStackBase> MakePrintStack(const Arguments& args)
    {
        using MakeFn = std::shared_ptr<PrintStackBase> (*)(const Arguments&);
        static std::map<Channel, MakeFn> makers = {
            { Channel::ETau, &MakePrintStackForChannel<static_cast<int>(Channel::ETau)> },
            { Channel::MuTau, &MakePrintStackForChannel<static_cast<int>(Channel::MuTau)> },
            { Channel::TauTau, &MakePrintStackForChannel<static_cast<int>(Channel::TauTau)> },
            { Channel::MuMu, &MakePrintStackForChannel<static_cast<int>(Channel::MuMu)> }
        };
        const Channel channel = args.channel();
        if(!makers.count(channel))
            throw exception("Unsupported channel '%1%'.") % channel;
        return makers.at(channel)(args);
    }

private:
    std::shared_ptr<PrintStackBase> print_base;
};

} // namespace analysis

PROGRAM_MAIN(analysis::PrintStack, Arguments)
