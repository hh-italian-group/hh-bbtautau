import ROOT

ROOT.gInterpreter.Declare('''

using IdType = unsigned long long;
static std::set<IdType> eventIds;

void AddEventId(IdType evt)
{
    static std::mutex m;
    std::lock_guard<std::mutex> lock(m);
    eventIds.insert(evt);
}

void CollectEventIds(const std::string& file_name)
{
    ROOT::RDataFrame df("tauTau", file_name);
    df.Foreach(AddEventId, {"evt"});
}

bool PassSel(IdType evt)
{
    return eventIds.count(evt) > 0;
}

''')


def skim(df, anaTupleFile):
    ROOT.CollectEventIds(anaTupleFile)
    return df.Filter('PassSel(event)')
