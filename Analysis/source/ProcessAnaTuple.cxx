/*! Post-processing of the analysis results.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Core/include/ProgressReporter.h"
#include "hh-bbtautau/Analysis/include/AnaTuple.h"
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
    OPT_ARG(size_t, n_parallel, 10);
    run::Argument<std::vector<std::string>> input_friends{"input_friends", "", {}};
};

class ProcessAnaTuple : public EventAnalyzerCore {
public:
    using AnaData = ::analysis::EventAnalyzerData;
    using AnaDataCollection = ::analysis::EventAnalyzerDataCollection;
    using PlotsProducer = ::analysis::StackedPlotsProducer;

    ProcessAnaTuple(const AnalyzerArguments& _args) :
        EventAnalyzerCore(_args, _args.channel(), false), args(_args), activeVariables(ParseVarSet(args.vars())),
        tupleReader(args.input(), args.channel(), activeVariables, args.input_friends()),
        outputFile(root_ext::CreateRootFile(args.output() + "_full.root"))
    {
        histConfig.Parse(FullPath(ana_setup.hist_cfg));
        if(!ana_setup.unc_cfg.empty()) {
            ConfigReader config_reader;
            unc_collection = std::make_shared<ModellingUncertaintyCollection>();
            ModellingUncertaintyEntryReader unc_reader(*unc_collection);
            config_reader.AddEntryReader("UNC", unc_reader, true);
            config_reader.ReadConfig(FullPath(ana_setup.unc_cfg));
        }

        if(!args.input_friends().empty()) {
            std::cout << "Input TTree friends:\n";
            for(const std::string& file_name : args.input_friends())
                std::cout << '\t' << file_name << '\n';
            std::cout << std::endl;
        }

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

        if(args.shapes()) {
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
    }

    void Run()
    {
        const std::set<std::string> signal_names(ana_setup.signals.begin(), ana_setup.signals.end());
        const std::set<std::string> bkg_names(ana_setup.backgrounds.begin(), ana_setup.backgrounds.end());
        const auto samplesToDraw = PlotsProducer::CreateOrderedSampleCollection(
                    ana_setup.draw_sequence, sample_descriptors, cmb_sample_descriptors, ana_setup.signals,
                    ana_setup.data, args.channel());

        std::ofstream qcd_out(args.output() +"_QCD.txt");

        const std::vector<EventSubCategory> all_subCategories(sub_categories_to_process.begin(),
                                                              sub_categories_to_process.end());

        for(size_t n = 0; n * args.n_parallel() < all_subCategories.size(); ++n) {
            AnaDataCollection anaDataCollection(outputFile, channelId, activeVariables, histConfig.GetItems(),
                                                false, bkg_names, unc_collection);
            EventSubCategorySet subCategories;
            for(size_t k = 0; k < args.n_parallel() && n * args.n_parallel() + k < all_subCategories.size(); ++k) {
                const auto& subCategory = all_subCategories.at(n * args.n_parallel() + k);
                subCategories.insert(subCategory);
                std::cout << subCategory << " ";
            }
            std::cout << std::endl;

            std::cout << "\tCreating histograms..." << std::endl;
            ProduceHistograms(anaDataCollection, subCategories);

            std::cout << "\tProcessing combined samples and QCD... " << std::endl;
            for(const auto& subCategory : subCategories) {

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
                    for(const auto& subCategory : subCategories)
                        limitsInputProducer.Produce(args.output(), limit_setup.first, limit_setup.second, subCategory,
                                                    ana_setup.unc_sources, ana_setup.regions, mva_sel_aliases,
                                                    args.period());
                }
            }
            if(args.draw()) {
                std::cout << "\tCreating plots..." << std::endl;
                PlotsProducer plotsProducer(anaDataCollection, samplesToDraw, FullPath(ana_setup.plot_cfg),
                                            ana_setup.plot_page_opt);
                std::string pdf_prefix = args.output();
                if(n != 0)
                    pdf_prefix += "_part" + ToString(n + 1);
                plotsProducer.PrintStackedPlots(pdf_prefix, EventRegion::SignalRegion(), ana_setup.categories,
                                                subCategories, signal_names);
            }
        }

        std::cout << "Saving output file..." << std::endl;
    }

private:

    //struct AnaDataFiller : public TObject {
    struct AnaDataFiller : ROOT::Detail::RDF::RActionImpl<AnaDataFiller> {
        using Result_t = bool;
        using Hist = EventAnalyzerData::Entry::Hist;
        using Mutex = Hist::Mutex;
        using DataId = bbtautau::AnaTupleReader::DataId;
        using HistMap = std::map<DataId, Hist*>;
        using LorentzVectorM = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>;
        // template <typename T> using VecType = std::vector<T>;

        const bbtautau::AnaTupleReader* tupleReader;
        AnaDataCollection* anaDataCollection;
        const EventCategorySet* categories;
        const EventSubCategorySet* subCategories;
        std::map<SelectionCut,analysis::EllipseParameters> massWindowParams;
        bool use_kinFit{false}, use_svFit{false};

        const std::set<UncertaintySource>* unc_sources;
        std::string hist_name;
        bool is_mva_score, is_limit_var;
        std::shared_ptr<HistMap> histograms;
        std::shared_ptr<Mutex> mutex;
        Result_t result{true};

        AnaDataFiller(const bbtautau::AnaTupleReader& _tupleReader, AnaDataCollection& _anaDataCollection,
                      const EventCategorySet& _categories, const EventSubCategorySet& _subCategories,
                      std::map<SelectionCut,analysis::EllipseParameters> _massWindowParams,
                      bool _use_kinFit, bool _use_svFit,
                      const std::set<UncertaintySource>& _unc_sources,
                      const std::string& _hist_name, bool _is_limit_var) :
                tupleReader(&_tupleReader), anaDataCollection(&_anaDataCollection), categories(&_categories),
                subCategories(&_subCategories),
                massWindowParams(_massWindowParams), use_kinFit(_use_kinFit), use_svFit(_use_svFit),
                unc_sources(&_unc_sources), hist_name(_hist_name),
                is_mva_score(_hist_name == "mva_score"), is_limit_var(_is_limit_var),
                histograms(std::make_shared<HistMap>()), mutex(std::make_shared<Mutex>()) {}
        AnaDataFiller(AnaDataFiller&) = delete;
        AnaDataFiller(AnaDataFiller&&) = default;
        //AnaDataFiller& operator=(const AnaDataFiller&) = default;
        // virtual ~AnaDataFiller() {}

        void Initialize() {}
        void InitTask(TTreeReader *, unsigned int) {}


        std::shared_ptr<bool> GetResultPtr() const { return std::make_shared<bool>(true); }

        void Finalize() {}

        std::string GetActionName() {return "AnaDataFiller";}

        Result_t& PartialUpdate(unsigned int slot){ return result; }


        template<typename T>
        void Exec(unsigned int slot, std::vector<size_t> dataId_hash_vec, std::vector<double> weight_vec,
                  bbtautau::AnaTupleReader::category_storage category_storage, int vbf_tag_raw,
                  bool has_b_pair, LorentzVectorM SVfit_p4, LorentzVectorM MET_p4,
                  double m_bb, double m_tt_vis, int kinFit_convergence, T&& value) const
        {

            const std::map<DiscriminatorWP, size_t> bjet_counts = {{DiscriminatorWP::Loose,
                                                                                   category_storage.num_btag_loose},
                                                                                  {DiscriminatorWP::Medium,
                                                                                   category_storage.num_btag_medium},
                                                                                  {DiscriminatorWP::Tight,
                                                                                   category_storage.num_btag_tight}};

            boost::optional<DiscriminatorWP> vbf_tag;
            if(vbf_tag_raw > 0) vbf_tag = static_cast<DiscriminatorWP>(vbf_tag_raw);

            for (size_t i=0; i<dataId_hash_vec.size(); i++){

                auto dataId_hash = dataId_hash_vec.at(i);
                const auto& dataId = tupleReader->GetDataIdByHash(dataId_hash);
                static const EventCategory evtCategory_2j = EventCategory::Parse("2j");
                if((dataId.Get<EventCategory>() != evtCategory_2j)
                        || !(unc_sources->count(dataId.Get<UncertaintySource>()))) continue;

                auto weight = weight_vec.at(i);

                for(const auto& category : *categories) {
                    if(!category.Contains(static_cast<size_t>(category_storage.num_jets), bjet_counts, category_storage.is_vbf,
                                     category_storage.is_boosted,vbf_tag) ) continue;
                    EventSubCategory evtSubCategory;
                    if(has_b_pair) {
                        if(category.HasBoostConstraint() && category.IsBoosted()) {
                            if(use_svFit) {
                                const bool isInsideBoostedCut = cuts::hh_bbtautau_Run2::hh_tag::IsInsideBoostedMassWindow(SVfit_p4.M(), m_bb);
                                evtSubCategory.SetCutResult(SelectionCut::mh, isInsideBoostedCut);
                            }
                        } else {
                            if(!use_svFit && massWindowParams.count(SelectionCut::mh))
                                throw exception("Category mh inconsistent with the false requirement of SVfit.");
                            if(massWindowParams.count(SelectionCut::mh)) {
                                const bool cut_result = use_svFit
                                            //&& event.GetSVFitResults(ana_setup.allow_calc_svFit).has_valid_momentum
                                            && massWindowParams.at(SelectionCut::mh).IsInside(
                                            SVfit_p4.M(), m_bb);
                                evtSubCategory.SetCutResult(SelectionCut::mh, cut_result);
                            }
                            if(massWindowParams.count(SelectionCut::mhVis))
                                evtSubCategory.SetCutResult(SelectionCut::mhVis,massWindowParams.at(SelectionCut::mhVis)
                                        .IsInside(m_tt_vis, m_bb));

                            if(massWindowParams.count(SelectionCut::mhMET))
                                evtSubCategory.SetCutResult(SelectionCut::mhMET,massWindowParams.at(SelectionCut::mhMET)
                                        .IsInside((SVfit_p4+MET_p4).M(), m_bb));

                        }
                        if(use_kinFit)
                            evtSubCategory.SetCutResult(SelectionCut::KinematicFitConverged,
                                                        kinFit_convergence);
                    }




                    for(const auto& subCategory : *subCategories) {
                        if(!evtSubCategory.Implies(subCategory)) continue;

                        const auto& dataId_correct = dataId.Set(category).Set(subCategory);

                        Hist* hist = GetHistogram(dataId_correct);
                        if(hist) {
                            auto x = value;
                            /* if(is_mva_score) {
                                const auto& dataId = tupleReader->GetDataIdByHash(dataId_hash);
                                x = static_cast<T>(tupleReader->GetNormalizedMvaScore(dataId, static_cast<float>(x)));
                            }*/

                            std::lock_guard<Hist::Mutex> lock(hist->GetMutex());
                            hist->Fill(x, weight);
                        }
                    }
                }

            }
            //(int) slot;
        }

        void Merge(TList*) {}

    private:
        Hist* GetHistogram(const DataId& dataId) const
        {
            std::lock_guard<Mutex> lock(*mutex);
            auto iter = histograms->find(dataId);
            if(iter != histograms->end())
                return iter->second;
            //const auto& dataId_temp = tupleReader->GetDataIdByHash(dataId_hash);
            //const auto& dataId = dataId_temp.Set(evtCategory).Set(evtSubCategory);
            Hist* hist = nullptr;
            //if(categories->count(dataId.Get<EventCategory>())
            //        && subCategories->count(dataId.Get<EventSubCategory>())
            //        &&
            //        && (is_limit_var || dataId.Get<UncertaintyScale>() == UncertaintyScale::Central)) {
            hist = &anaDataCollection->Get(dataId).GetHistogram(hist_name)();
            // }
            (*histograms)[dataId] = hist;
            return hist;
        }
    };

    template <typename T> using VecType = ROOT::VecOps::RVec<T>;

    void ProduceHistograms(AnaDataCollection& anaDataCollection, const EventSubCategorySet& subCategories)
    {
        const auto has_column = [](bbtautau::AnaTupleReader::RDF& df, const std::string& column_name) {
            auto columns = df.GetColumnNames();
            auto iter = std::find(columns.begin(), columns.end(), column_name);
            return iter != columns.end();
        };

        const auto is_defined_column = [](bbtautau::AnaTupleReader::RDF& df, const std::string& column_name) {
            auto columns = df.GetDefinedColumnNames();
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
            std::cout << hist_name << " ";
            const std::string df_hist_name = hist_name == "mva_score" ? "all_mva_scores" : hist_name;
            const std::vector<std::string> branches = {"dataIds", "all_weights", "category_storage", "vbf_tag_raw",
                                                       "has_b_pair", "SVfit_p4", "MET_p4", "m_bb", "m_tt_vis",
                                                       "kinFit_convergence",df_hist_name};
            //const std::vector<std::string> branches = {"category_storage", df_hist_name};
            AnaDataFiller filter(tupleReader, anaDataCollection, ana_setup.categories, subCategories,
                                 ana_setup.massWindowParams, ana_setup.use_kinFit, ana_setup.use_svFit,
                                 ana_setup.unc_sources, hist_name, true) ; //limitVariables.count(hist_name));
            auto df = get_df(hist_name);
            ROOT::RDF::RResultPtr<bool> result;
            //if(filter.is_mva_score)
            //   result = df.Book< bbtautau::AnaTupleReader::category_storage,
                //result = df.Fill<VecType<size_t>, VecType<double>, bbtautau::AnaTupleReader::category_storage,
            //              VecType<float>>(std::move(filter), branches);
            //else if(bbtautau::AnaTupleReader::BoolBranches.count(df_hist_name))
            if(bbtautau::AnaTupleReader::BoolBranches.count(df_hist_name))
                //result = df.Book< bbtautau::AnaTupleReader::category_storage,
                result = df.Book<std::vector<size_t>, std::vector<double>, bbtautau::AnaTupleReader::category_storage,
                        int, bool, LorentzVectorM, LorentzVectorM, double, double, int,
                        bool>(std::move(filter), branches);
            else if(bbtautau::AnaTupleReader::IntBranches.count(df_hist_name))
                //result = df.Book< bbtautau::AnaTupleReader::category_storage,
                result = df.Book<std::vector<size_t>, std::vector<double>, bbtautau::AnaTupleReader::category_storage,
                        int, bool, LorentzVectorM, LorentzVectorM, double, double, int,
                        int>(std::move(filter), branches);
            else if(is_defined_column(df, hist_name))
                //result = df.Book< bbtautau::AnaTupleReader::category_storage,
                result = df.Book<std::vector<size_t>, std::vector<double>, bbtautau::AnaTupleReader::category_storage,
                        int, bool, LorentzVectorM, LorentzVectorM, double, double, int,
                        double>(std::move(filter), branches);
            else
                //result = df.Book<bbtautau::AnaTupleReader::category_storage,
                result = df.Book<std::vector<size_t>, std::vector<double>, bbtautau::AnaTupleReader::category_storage,
                        int, bool, LorentzVectorM, LorentzVectorM, double, double, int,
                        float>(std::move(filter), branches);
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
    bbtautau::AnaTupleReader tupleReader;
    std::shared_ptr<TFile> outputFile;
    PropertyConfigReader histConfig;
    std::shared_ptr<ModellingUncertaintyCollection> unc_collection;
};

} // namespace analysis

PROGRAM_MAIN(analysis::ProcessAnaTuple, analysis::AnalyzerArguments)
