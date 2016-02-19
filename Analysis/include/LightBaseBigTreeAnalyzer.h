/*! Definition of LightBaseBigTreeAnalyzer class, the base class for separate studies on big trees.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "Analysis/source/FlatTreeProducer_etau.cxx"
#include "Analysis/source/FlatTreeProducer_mutau.cxx"
#include "Analysis/source/FlatTreeProducer_tautau.cxx"

namespace analysis {

class LightBaseBigTreeAnalyzer : public FlatTreeProducer_etau, public FlatTreeProducer_mutau,
                                 public FlatTreeProducer_tautau {
public:
    LightBaseBigTreeAnalyzer(const std::string& channelName, const std::string& inputFileName,
                             const std::string& outputFileName, const std::string& configFileName,
                             const std::string& _prefix = "none", size_t _maxNumberOfEvents = 0)
        : BaseFlatTreeProducer(inputFileName, outputFileName, configFileName, _prefix, _maxNumberOfEvents),
          FlatTreeProducer_etau(inputFileName, outputFileName, configFileName, _prefix, _maxNumberOfEvents),
          FlatTreeProducer_mutau(inputFileName, outputFileName, configFileName, _prefix, _maxNumberOfEvents),
          FlatTreeProducer_tautau(inputFileName, outputFileName, configFileName, _prefix, _maxNumberOfEvents)
    {
        std::istringstream ss_channel(channelName);
        ss_channel >> channel;
        if(channel != Channel::ETau && channel != Channel::MuTau && channel != Channel::TauTau)
            throw exception("Unsupported channel ") << channel;
        writeFlatTree = false;
    }

    virtual analysis::BaseAnalyzerData& GetAnaData() override
    {
        if(channel == Channel::ETau)
            return HHbbetau_FlatTreeProducer::GetAnaData();
        if(channel == Channel::MuTau)
            return HHbbmutau_FlatTreeProducer::GetAnaData();
        return HHbbtautau_FlatTreeProducer::GetAnaData();
    }

    virtual void ProcessEvent(std::shared_ptr<const analysis::EventDescriptor> _event) override
    {
        BaseAnalyzer::ProcessEvent(_event);

        const SelectionResults& selection = ApplyBaselineSelection();

        AnalyzeSelection(selection);
        if(channel == Channel::ETau)
            AnalyzeETauSelection(HHbbetau_FlatTreeProducer::selection);
        else if(channel == Channel::MuTau)
            AnalyzeMuTauSelection(HHbbmutau_FlatTreeProducer::selection);
        else
            AnalyzeTauTauSelection(HHbbtautau_FlatTreeProducer::selection);
    }

protected:
    virtual void AnalyzeSelection(const SelectionResults& selection) {}
    virtual void AnalyzeETauSelection(const SelectionResults_etau& selection) {}
    virtual void AnalyzeMuTauSelection(const SelectionResults_mutau& selection) {}
    virtual void AnalyzeTauTauSelection(const SelectionResults_tautau& selection) {}

    virtual analysis::SelectionResults& ApplyBaselineSelection() override
    {
        if(channel == Channel::ETau)
            return HHbbetau_FlatTreeProducer::ApplyBaselineSelection();
        else if(channel == Channel::MuTau)
            return HHbbmutau_FlatTreeProducer::ApplyBaselineSelection();
        return HHbbtautau_FlatTreeProducer::ApplyBaselineSelection();
    }

    virtual void CalculateTriggerWeights(const analysis::Candidate& higgs) override
    {
        if(channel == Channel::ETau)
            HHbbetau_FlatTreeProducer::CalculateTriggerWeights(higgs);
        else if(channel == Channel::MuTau)
            HHbbmutau_FlatTreeProducer::CalculateTriggerWeights(higgs);
        else
            HHbbtautau_FlatTreeProducer::CalculateTriggerWeights(higgs);
    }

    virtual void CalculateIsoWeights(const analysis::Candidate& higgs) override
    {
        if(channel == Channel::ETau)
            HHbbetau_FlatTreeProducer::CalculateIsoWeights(higgs);
        else if(channel == Channel::MuTau)
            HHbbmutau_FlatTreeProducer::CalculateIsoWeights(higgs);
        else
            HHbbtautau_FlatTreeProducer::CalculateIsoWeights(higgs);
    }

    virtual void CalculateIdWeights(const analysis::Candidate& higgs) override
    {
        if(channel == Channel::ETau)
            HHbbetau_FlatTreeProducer::CalculateIdWeights(higgs);
        else if(channel == Channel::MuTau)
            HHbbmutau_FlatTreeProducer::CalculateIdWeights(higgs);
        else
            HHbbtautau_FlatTreeProducer::CalculateIdWeights(higgs);
    }

    virtual void CalculateDMWeights(const analysis::Candidate& higgs) override
    {
        if(channel == Channel::ETau)
            HHbbetau_FlatTreeProducer::CalculateDMWeights(higgs);
        else if(channel == Channel::MuTau)
            HHbbmutau_FlatTreeProducer::CalculateDMWeights(higgs);
        else
            HHbbtautau_FlatTreeProducer::CalculateDMWeights(higgs);
    }

    virtual analysis::Candidate SelectTau(size_t id, cuts::ObjectSelector* objectSelector,
                                          root_ext::AnalyzerData& _anaData,
                                          const std::string& selection_label) override
    {
        if(channel == Channel::ETau)
            return HHbbetau_FlatTreeProducer::SelectTau(id, objectSelector, _anaData, selection_label);
        else if(channel == Channel::MuTau)
            return HHbbmutau_FlatTreeProducer::SelectTau(id, objectSelector, _anaData, selection_label);
        return HHbbtautau_FlatTreeProducer::SelectTau(id, objectSelector, _anaData, selection_label);
    }

    virtual analysis::Candidate SelectSignalTau(size_t id, cuts::ObjectSelector* objectSelector,
                                                root_ext::AnalyzerData& _anaData,
                                                const std::string& selection_label) override
    {
        if(channel == Channel::ETau)
            return HHbbetau_FlatTreeProducer::SelectSignalTau(id, objectSelector, _anaData, selection_label);
        else if(channel == Channel::MuTau)
            return HHbbmutau_FlatTreeProducer::SelectSignalTau(id, objectSelector, _anaData, selection_label);
        return HHbbtautau_FlatTreeProducer::SelectSignalTau(id, objectSelector, _anaData, selection_label);
    }

    virtual void FillFlatTree(const analysis::SelectionResults& selection) override
    {
        if(channel == Channel::ETau)
            HHbbetau_FlatTreeProducer::FillFlatTree(selection);
        else if(channel == Channel::MuTau)
            HHbbmutau_FlatTreeProducer::FillFlatTree(selection);
        else
            HHbbtautau_FlatTreeProducer::FillFlatTree(selection);
    }

private:
    Channel channel;
};

} // namespace analysis
