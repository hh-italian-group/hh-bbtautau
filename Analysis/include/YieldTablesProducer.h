/*! Definition of BaseEventAnalyzer class, the base class for event analyzers.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

namespace analysis {

class YieldTablesProducer {
public:
    void PrintTables(const std::string& name_suffix, const std::wstring& sep)
    {
        std::wofstream of(args.outputFileName() + "_" + name_suffix + ".csv");

        static const std::set< std::pair<std::string, EventSubCategory> > interesting_histograms = {
            { EventAnalyzerData::m_sv_Name(), EventSubCategory::NoCuts },
            { EventAnalyzerData::m_sv_Name(), EventSubCategory::MassWindow },
            { EventAnalyzerData::m_ttbb_kinfit_Name(), EventSubCategory::KinematicFitConverged },
            { EventAnalyzerData::m_ttbb_kinfit_Name(), EventSubCategory::KinematicFitConvergedWithMassWindow }
        };

        for(const auto& hist_entry : interesting_histograms)
            PrintTables(of, sep, hist_entry.first, hist_entry.second, false, true);

        of.flush();
        of.close();
    }

    void PrintTables(std::wostream& of, const std::wstring& sep, const std::string& hist_name,
                     EventSubCategory subCategory, bool includeOverflow, bool includeError)
    {
        of << std::wstring(hist_name.begin(), hist_name.end());

        std::wstring table_name_suffix = L"";
        if(includeOverflow && includeError)
            table_name_suffix = L" with overflow and error";
        else if(includeOverflow && !includeError)
            table_name_suffix = L" with overflow";
        else if(!includeOverflow && includeError)
            table_name_suffix = L" with error";
        of << table_name_suffix << sep;

        for (EventCategory eventCategory : EventCategoriesToProcess())
            of << eventCategory << sep;
        of << std::endl;

        for (const DataCategory* dataCategory : dataCategoryCollection.GetAllCategories()) {
            of << std::wstring(dataCategory->title.begin(), dataCategory->title.end()) << sep;
            for (EventCategory eventCategory : EventCategoriesToProcess()) {
                const EventAnalyzerDataMetaId_noRegion_noName meta_id(eventCategory, subCategory,
                                                                     EventEnergyScale::Central);

                if( TH1D* histogram = GetSignalHistogram(meta_id, dataCategory->name, hist_name) ) {
                    const PhysicalValue integral = Integral(*histogram, includeOverflow);
                    of << integral.ToString<wchar_t>(includeError, false) << sep;
                }
                else
                    of << "not found" << sep;
            }
            of << std::endl;
        }
        of << std::endl << std::endl;
    }
};

} // namespace analysis
