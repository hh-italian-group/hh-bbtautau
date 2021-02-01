/*! Post-processing of the analysis results.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Core/include/ProgressReporter.h"
#include "hh-bbtautau/Analysis/include/AnaTupleReader.h"
#include "hh-bbtautau/Analysis/include/EventAnalyzerCore.h"
#include "hh-bbtautau/Analysis/include/EventAnalyzerDataCollection.h"
#include "hh-bbtautau/Analysis/include/LimitsInputProducer.h"
#include "hh-bbtautau/Analysis/include/SampleDescriptorConfigEntryReader.h"
#include "hh-bbtautau/Analysis/include/StackedPlotsProducer.h"

namespace analysis {

struct AnalyzerArguments : CoreAnalyzerArguments {
    REQ_ARG(Channel, channel);
    REQ_ARG(Period, period);
    REQ_ARG(std::string, input);
    REQ_ARG(std::string, output);
    OPT_ARG(bool, shapes, true);
    OPT_ARG(bool, draw, true);
    OPT_ARG(std::string, vars, "");
    OPT_ARG(std::string, input_friends, "");
};

class ProcessAnaTuple : public EventAnalyzerCore {
public:
    using AnaData = ::analysis::EventAnalyzerData;
    using AnaDataCollection = ::analysis::EventAnalyzerDataCollection;
    using PlotsProducer = ::analysis::StackedPlotsProducer;

    ProcessAnaTuple(const AnalyzerArguments& _args) :
        EventAnalyzerCore(_args, _args.channel(), false), args(_args), activeVariables(ParseVarSet(args.vars())),
        outputFile(root_ext::CreateRootFile(args.output() + "_full.root", ROOT::kLZMA, 9))

    {
        for (auto& hist_config : ana_setup.hist_cfg)
            histConfig.Parse(FullPath(hist_config));
        if(!ana_setup.unc_cfg.empty()) {
            ConfigReader config_reader;
            unc_collection = std::make_shared<ModellingUncertaintyCollection>();
            ModellingUncertaintyEntryReader unc_reader(*unc_collection);
            config_reader.AddEntryReader("UNC", unc_reader, true);
            config_reader.ReadConfig(FullPath(ana_setup.unc_cfg));
        }
        std::set_union(ana_setup.norm_unc_sources.begin(), ana_setup.norm_unc_sources.end(),
                        unc_sources_group.begin(), unc_sources_group.end(),
                        std::inserter(unc_sources_total, unc_sources_total.begin()));


    }

    void Run()
    {
        const std::set<std::string> signal_names(ana_setup.signals.begin(), ana_setup.signals.end());
        const std::set<std::string> bkg_names(ana_setup.backgrounds.begin(), ana_setup.backgrounds.end());
        const auto samplesToDraw = PlotsProducer::CreateOrderedSampleCollection(
                    ana_setup.draw_sequence, sample_descriptors, cmb_sample_descriptors, ana_setup.signals,
                    ana_setup.data, args.channel());

        std::ofstream qcd_out(args.output() +"_QCD.txt");

        AnaDataCollection anaDataCollection(outputFile, channelId, activeVariables, histConfig.GetItems(),
                                            false, bkg_names, unc_collection);
        const EventSubCategorySet& subCategories = sub_categories_to_process;

        std::vector<std::string> unc_sources = SplitValueList(args.unc_sources_groups(),false, ",", true);
        for (auto& unc_source : unc_sources){
            std::string input_file = args.input();
            boost::replace_all(input_file, "{UNC_GROUP}", unc_source);
            std::cout << "\n input file \t " ;
            std::cout << input_file << std::endl;
            std::vector<std::string> input_friends = SplitValueList(args.input_friends(),false, ",", true);
            for (auto& f: input_friends) {
                boost::replace_all(f, "{UNC_GROUP}", unc_source);
            }
            if(!input_friends.empty()) {
                std::cout << "Input TTree friends:\n";
                for(const std::string& file_name : input_friends)
                    std::cout << '\t' << file_name << '\n';
                std::cout << std::endl;
            }
            bbtautau::EventTagCreator eventTagger(ana_setup.categories, sub_categories_to_process, unc_sources_total, ana_setup.norm_unc_sources,ana_setup.use_IterativeFit, ana_setup.r_factors_file.at(args.channel()));
            int hastune = 0; // 1 for 2016 etau and tautau, 2 for 2016 muTau, 0 for the rest
            if(args.channel()!= Channel::MuTau && args.period()==Period::Run2016) hastune=1;
            else if (args.channel()== Channel::MuTau && args.period()==Period::Run2016) hastune=2;
            bbtautau::AnaTupleReader tupleReader(input_file, args.channel(), activeVariables, input_friends, eventTagger,
                                                    ana_setup.mdnn_version, ana_setup.norm_unc_sources, hastune);

            if(!tupleReader.GetParametricVariables().empty()) {
                std::cout << "Parametric variables:\n";
                for(const auto& [name, columns] : tupleReader.GetParametricVariables()) {
                    std::cout << '\t' << name << ": ";
                    for(const auto& column : columns) {
                        std::cout << column << " ";
                    }
                    std::cout << '\n';
                }
                std::cout << std::endl;
            }
            if(args.shapes() && limitVariables.empty()) {
                for(const auto& limit_setup : ana_setup.limit_setup){
                    for(const auto& [cat, var] : limit_setup.second) {
                        if(limitVariables.count(var)) continue;
                        bool var_found = false;
                        if(activeVariables.count(var)) {
                            var_found = true;
                        } else if(var.back() == '+') {
                            const std::string var_name = var.substr(0, var.size() - 1);
                            var_found = tupleReader.GetParametricVariables().count(var_name);
                        }

                        if(!var_found) {
                            throw exception("Variable '%1%' is used in the limit setup '%2%', but it is not enabled."
                                            " Consider to enable it or run with option --shapes 0.")
                                            % var % limit_setup.first;
                        }
                        limitVariables.insert(var);
                    }
                }
            }
            std::cout << "\tCreating histograms..." << std::endl;
            ProduceHistograms(anaDataCollection, tupleReader);
        }
        std::cout << "\tProcessing combined samples and QCD... " << std::endl;
        for(const auto& subCategory : subCategories) {
            std::cout << "\t\tsub-category " << subCategory << "\n";
            ProcessCombinedSamples(anaDataCollection, subCategory, ana_setup.cmb_samples);
            for(const auto& sample : sample_descriptors) {
                if(sample.second.sampleType == SampleType::QCD) {
                    EstimateQCD(anaDataCollection, subCategory, sample.second, qcd_out);
                    break;
                }
            }
        }
        if(args.shapes()) {
            std::cout << "\tProducing inputs for limits..." << std::endl;
            LimitsInputProducer limitsInputProducer(anaDataCollection, sample_descriptors,
                                                    cmb_sample_descriptors);
            for(const auto& limit_setup : ana_setup.limit_setup) {
                std::cout << "\t\tsetup_name: " << limit_setup.first <<  std::endl;
                for(const auto& subCategory : subCategories){
                    EventRegionSet regions_forlimits={EventRegion::OS_Isolated()};
                    //regions_forlimits instead of ana_setup.regions
                    limitsInputProducer.Produce(args.output(), limit_setup.first, limit_setup.second, subCategory,
                                                unc_sources_total, regions_forlimits, mva_sel_aliases,
                                                args.period());
                }
            }
        }
        if(args.draw()) {
            std::cout << "\tCreating plots..." << std::endl;
            PlotsProducer plotsProducer(anaDataCollection, samplesToDraw, FullPath(ana_setup.plot_cfg),
                                        ana_setup.plot_page_opt);
            std::string pdf_prefix = args.output();
            plotsProducer.PrintStackedPlots(pdf_prefix, EventRegion::SignalRegion(), ana_setup.categories,
                                            subCategories, signal_names);
        }
        std::cout << "Saving output file..." << std::endl;
    }

private:

    template<typename _Key>
    class MultiHist {
    public:
        using Key = _Key;
        using Real = float;
        using BinV = std::vector<Real>;
        using Map = std::map<Key, BinV>;
        using const_iterator = typename Map::const_iterator;

        MultiHist(const TAxis& _axis) :
            axis(_axis), n_bins_total(static_cast<size_t>(axis.GetNbins()) + 2)
        {
        }

        void Fill(const Key& key, double value, Real weight)
        {
            auto iter = bins.find(key);
            if(iter == bins.end())
                iter = bins.insert({key, BinV(2 * n_bins_total, 0.f)}).first;
            const int bin = axis.FindFixBin(value);
            const size_t n = static_cast<size_t>(bin);
            iter->second.at(n) += weight;
            iter->second.at(n_bins_total + n) += static_cast<float>(std::pow(weight, 2));
        }

        void AddTo(const Key& key, TH1& hist) const
        {
            auto iter = bins.find(key);
            if(iter == bins.end())
                return;
            for(size_t n = 0; n < n_bins_total; ++n) {
                const int i = static_cast<int>(n);
                const double bin_value = hist.GetBinContent(i) + iter->second.at(n);
                const double bin_error2 = std::pow(hist.GetBinError(i), 2) + iter->second.at(n_bins_total + n);
                hist.SetBinContent(i, bin_value);
                hist.SetBinError(i, std::sqrt(bin_error2));
            }
        }

        const_iterator begin() const { return bins.begin(); }
        const_iterator end() const { return bins.end(); }

    private:
        const TAxis axis;
        const size_t n_bins_total;
        Map bins;
    };

    //struct AnaDataFiller : public TObject {
    struct AnaDataFiller : ROOT::Detail::RDF::RActionImpl<AnaDataFiller> {
        using Result_t = bool;
        using Hist = EventAnalyzerData::Entry::Hist;
        using Mutex = Hist::Mutex;
        using DataId = bbtautau::AnaTupleReader::DataId;
        using HKey = std::tuple<EventRegion, UncertaintySource, UncertaintyScale, std::string>;
        using MultiH = MultiHist<HKey>;
        using MapKey = std::tuple<EventCategory, EventSubCategory>;
        using HistMap = std::map<MapKey, MultiH>;

        AnaDataCollection& anaDataCollection;
        const std::string hist_name;
        std::vector<std::shared_ptr<HistMap>> histograms;
        std::shared_ptr<Mutex> mutex;
        std::shared_ptr<Result_t> result;

        AnaDataFiller(AnaDataCollection& _anaDataCollection, const std::string& _hist_name, size_t n_slots) :
                      anaDataCollection(_anaDataCollection), hist_name(_hist_name), histograms(n_slots),
                      mutex(std::make_shared<Mutex>()), result(std::make_shared<bool>(false))
        {
            for(size_t n = 0; n < n_slots; ++n) {
                histograms.at(n) = std::make_shared<HistMap>();
            }
        }

        AnaDataFiller(AnaDataFiller&) = delete;
        AnaDataFiller(AnaDataFiller&&) = default;

        void Initialize() {}
        void InitTask(TTreeReader *, unsigned int) {}
        std::shared_ptr<bool> GetResultPtr() const { return result; }
        std::string GetActionName() const { return "AnaDataFiller"; }
        Result_t& PartialUpdate(unsigned int) { return *result; }


        template<typename T>
        void Exec(unsigned int slot, const bbtautau::EventTags& event_tags, T&& value)
        {
            for(size_t i = 0; i < event_tags.dataIds.size(); ++i) {
                const auto& dataId = event_tags.dataIds.at(i);
                const double weight = event_tags.weights.at(i);
                const HKey h_key(dataId.Get<EventRegion>(), dataId.Get<UncertaintySource>(),
                                 dataId.Get<UncertaintyScale>(), dataId.Get<std::string>());
                MultiH& hist = GetHistogram(slot, dataId);
                //std::lock_guard<Hist::Mutex> lock(hist.GetMutex());
                hist.Fill(h_key, value, static_cast<float>(weight));
            }
        }

        void Finalize()
        {
            std::lock_guard<Mutex> lock(*mutex);

            std::set<DataId> all_dataIds;
            for(size_t slot = 0; slot < histograms.size(); ++slot) {
                for(const auto& [m_key, hist] : *histograms.at(slot)) {
                    for(const auto& [h_key, bins] : hist) {
                        const DataId dataId(std::get<EventCategory>(m_key), std::get<EventSubCategory>(m_key),
                                            std::get<EventRegion>(h_key), std::get<UncertaintySource>(h_key),
                                            std::get<UncertaintyScale>(h_key), std::get<std::string>(h_key));
                        all_dataIds.insert(dataId);
                    }
                }

            }
            for(const auto& dataId : all_dataIds) {
                Hist& base_hist = anaDataCollection.Get(dataId).GetHistogram(hist_name)();
                for(size_t slot = 0; slot < histograms.size(); ++slot) {
                    const MapKey m_key(dataId.Get<EventCategory>(), dataId.Get<EventSubCategory>());
                    const HKey h_key(dataId.Get<EventRegion>(), dataId.Get<UncertaintySource>(),
                                     dataId.Get<UncertaintyScale>(), dataId.Get<std::string>());
                    const auto iter = histograms.at(slot)->find(m_key);
                    if(iter != histograms.at(slot)->end()) {
                        iter->second.AddTo(h_key, base_hist);
                    }
                }
            }

            histograms.clear();
            *result = true;
        }

    private:
        MultiH& GetHistogram(unsigned int slot, const DataId& dataId)
        {
            if(slot >= histograms.size())
                throw exception("Slot is out of range");

            const MapKey m_key(dataId.Get<EventCategory>(), dataId.Get<EventSubCategory>());
            const auto iter = histograms.at(slot)->find(m_key);
            if(iter != histograms.at(slot)->end())
                return iter->second;

            std::lock_guard<Mutex> lock(*mutex);
            Hist& hist = anaDataCollection.Get(dataId).GetHistogram(hist_name)();
            histograms.at(slot)->insert({m_key, *hist.GetXaxis()});
            return histograms.at(slot)->at(m_key);
        }
    };

    template <typename T> using VecType = ROOT::VecOps::RVec<T>;

    void ProduceHistograms(AnaDataCollection& anaDataCollection, bbtautau::AnaTupleReader& tupleReader)
    {
        const auto has_column = [](bbtautau::AnaTupleReader::RDF& df, const std::string& column_name) {
            auto columns = df.GetColumnNames();
            auto iter = std::find(columns.begin(), columns.end(), column_name);
            return iter != columns.end();
        };

        const auto get_df = [&](const std::string& hist_name) {
            auto main_df = tupleReader.GetDataFrame();
            if(has_column(main_df, hist_name)) return main_df;
            for(auto df : tupleReader.GetSkimmedDataFrames())
                if(has_column(df, hist_name)) return df;
            throw exception("ProduceHistograms: Column with name '%1%' not found.") % hist_name;
        };

        std::vector<ROOT::RDF::RResultPtr<bool>> results;
        std::cout << "\t\tAdding: ";
        for(const auto& hist_name : activeVariables) {

            auto branch_type_ptr = tupleReader.TryGetBranchType(hist_name);
            const std::string branch_type = branch_type_ptr ? *branch_type_ptr : "Float_t";

            std::cout << hist_name << " ";
            const std::vector<std::string> branches = { "event_tags", hist_name };
            AnaDataFiller filter(anaDataCollection, hist_name, args.n_threads());
            auto df = get_df(hist_name);
            ROOT::RDF::RResultPtr<bool> result;
            if(branch_type == "Bool_t")
                result = df.Book<bbtautau::EventTags, bool>(std::move(filter), branches);
            else if(branch_type == "Int_t")
                result = df.Book<bbtautau::EventTags, int>(std::move(filter), branches);
            else if(branch_type == "Float_t")
                result = df.Book<bbtautau::EventTags, float>(std::move(filter), branches);
            else
                throw exception("Branch type = '%1%' does not supported.") % branch_type;
            results.push_back(result);
        }
        std::cout << std::endl;

        analysis::tools::ProgressReporter progressReporter(10, std::cout, "Filling histograms...");
        const size_t n_total = tupleReader.GetNumberOfEntries();
        progressReporter.SetTotalNumberOfEvents(n_total);
        static constexpr size_t counts_per_step = 1000;
        size_t n_processed = 0;
        AnaDataFiller::Mutex mutex;
        results.front().OnPartialResultSlot(counts_per_step, [&](unsigned int, bool&) {
            std::lock_guard<AnaDataFiller::Mutex> lock(mutex);
            n_processed += counts_per_step;
            progressReporter.Report(n_processed, false);
        });

        for(auto& result : results)
            result.GetValue();
        progressReporter.Report(n_total, true);
        // std::cout << "number of event loops: " << df.GetNRuns() << std::endl;
    }

    void EstimateQCD(AnaDataCollection& anaDataCollection, const EventSubCategory& subCategory,
                     const SampleDescriptor& qcd_sample, std::ostream& log)
    {
        static const EventRegionSet sidebandRegions = {
            EventRegion::OS_AntiIsolated(), EventRegion::SS_Isolated(), EventRegion::SS_AntiIsolated(),
            // EventRegion::SS_LooseIsolated()
        };
        static const std::set<UncertaintySource> qcdUncSources = { UncertaintySource::None };
        static const std::set<UncertaintyScale> qcdUncScales = { UncertaintyScale::Central };
        const EventSubCategorySet subCategories = { subCategory };

        for(const EventAnalyzerDataId& metaDataId : EventAnalyzerDataId::MetaLoop(ana_setup.categories,
                subCategories, sidebandRegions, qcdUncSources, qcdUncScales)) {
            const auto qcdAnaDataId = metaDataId.Set(qcd_sample.name);
            auto& qcdAnaData = anaDataCollection.Get(qcdAnaDataId);
            for(const auto& sample_name : sample_descriptors) {
                const SampleDescriptor& sample =  sample_name.second;
                if(sample.sampleType == SampleType::QCD) continue;
                if(ana_setup.IsSignal(sample.name)) continue;
                double factor = sample.sampleType == SampleType::Data ? +1 : -1;
                for(const auto& sample_wp : sample.working_points) {
                    const auto anaDataId = metaDataId.Set(sample_wp.full_name);
                    auto& anaData = anaDataCollection.Get(anaDataId);
                    for(const auto& sub_entry : anaData.template GetEntriesEx<TH1D>()) {
                        auto& entry = qcdAnaData.template GetEntryEx<TH1D>(sub_entry.first);
                        for(const auto& hist : sub_entry.second->GetHistograms()) {
                            entry(hist.first).Add(hist.second.get(), factor);
                        }
                    }
                }
            }
        }

        for(const EventAnalyzerDataId& metaDataId : EventAnalyzerDataId::MetaLoop(ana_setup.categories,
                subCategories, qcdUncSources, qcdUncScales)) {
            const auto anaDataId = metaDataId.Set(qcd_sample.name);
            auto& osIsoData = anaDataCollection.Get(anaDataId.Set(EventRegion::OS_Isolated()));
            auto& ssIsoData = anaDataCollection.Get(anaDataId.Set(EventRegion::SS_Isolated()));
            auto& osAntiIsoData = anaDataCollection.Get(anaDataId.Set(EventRegion::OS_AntiIsolated()));
            auto& ssAntiIsoData = anaDataCollection.Get(anaDataId.Set(EventRegion::SS_AntiIsolated()));
            auto& shapeData = anaDataCollection.Get(anaDataId.Set(ana_setup.qcd_shape));

            for(const auto& sub_entry : ssIsoData.template GetEntriesEx<TH1D>()) {
                auto& entry_osIso = osIsoData.template GetEntryEx<TH1D>(sub_entry.first);
                auto& entry_osAntiIso = osAntiIsoData.template GetEntryEx<TH1D>(sub_entry.first);
                auto& entry_ssAntiIso = ssAntiIsoData.template GetEntryEx<TH1D>(sub_entry.first);
                auto& entry_shape = shapeData.template GetEntryEx<TH1D>(sub_entry.first);

                for(const auto& hist : sub_entry.second->GetHistograms()) {
                    log << anaDataId << ": " << sub_entry.first << " " << hist.first << "\n";
                    const auto osAntiIso_integral = Integral(entry_osAntiIso(hist.first), true);
                    const auto ssAntiIso_integral = Integral(entry_ssAntiIso(hist.first), true);
                    if (osAntiIso_integral.GetValue() <= 0 || osAntiIso_integral.IsCompatible(PhysicalValue::Zero)){
                        log << "Warning: OS Anti Iso integral is too small " << hist.first << std::endl;
                        if(ana_setup.qcd_ss_os_sf <= 0 ) continue;
                    }

                    if (ssAntiIso_integral.GetValue() <= 0 || ssAntiIso_integral.IsCompatible(PhysicalValue::Zero)){
                        log << "Warning: SS Anti Iso integral is too small " << hist.first << std::endl;
                        if(ana_setup.qcd_ss_os_sf <= 0) continue;
                    }
                    PhysicalValue k_factor(ana_setup.qcd_ss_os_sf,ana_setup.qcd_ss_os_err);

                    const auto ssIso_integral = analysis::Integral(*hist.second, true);
                    if(ana_setup.qcd_ss_os_sf <=0 ) k_factor = ssIso_integral / ssAntiIso_integral;
                    if (osAntiIso_integral.GetValue() <= 0){
                        log << "Warning: SS Iso integral less or equal 0 for " << hist.first << std::endl;
                        continue;
                    }
                    const auto total_yield = osAntiIso_integral * k_factor;

                    log << anaDataId << ": osAntiIso integral = " << osAntiIso_integral
                        << ", ssAntiIso integral = " << ssAntiIso_integral << ", os/ss sf = " << k_factor
                        << ", ssIso integral = " << ssIso_integral << ", total yield = " << total_yield << std::endl;

                    TH1D shape_hist(entry_shape(hist.first));
                    std::string debug_info, negative_bins_info;
                    if(!FixNegativeContributions(shape_hist, debug_info, negative_bins_info)) {
                        log << debug_info << "\n" << negative_bins_info << "\n";
                        continue;
                    }
                    entry_osIso(hist.first).CopyContent(shape_hist);
                    analysis::RenormalizeHistogram(entry_osIso(hist.first), total_yield.GetValue(), true);
                }
            }
        }
        log << std::endl;
    }

    static std::set<std::string> ParseVarSet(const std::string& active_vars_str)
    {
        const auto list = SplitValueList(active_vars_str, false, ", \t", true);
        return std::set<std::string>(list.begin(), list.end());
    }

private:
    AnalyzerArguments args;
    std::set<std::string> activeVariables, limitVariables;
    std::shared_ptr<TFile> outputFile;
    PropertyConfigReader histConfig;
    std::shared_ptr<ModellingUncertaintyCollection> unc_collection;
    std::set<UncertaintySource> unc_sources_total;
};

} // namespace analysis

PROGRAM_MAIN(analysis::ProcessAnaTuple, analysis::AnalyzerArguments)
