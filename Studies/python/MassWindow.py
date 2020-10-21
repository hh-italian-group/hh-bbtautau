import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import ROOT

ROOT.gInterpreter.ProcessLine("""
using LorentzVectorXYZ = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>;
using LorentzVectorM = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>;
using LorentzVectorE = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double>>;
using LorentzVector = LorentzVectorE;
""")

def weighted_quantile(values, quantiles, sample_weight=None,
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, math.sqrt(variance))


year=2018
channel = "tauTau"
folder = "/eos/user/k/kandroso/cms-it-hh-bbtautau/anaTuples/2020-05-07/"
signal_name = 'GluGluSignal_NonRes_kl1'
#signal_name = 'VBFSignal_NonRes_VBFHH-CV_1_C2V_1_C3_0'
#signal_name = '(GluGluSignal_NonRes_kl1|VBFSignal_NonRes_VBFHH-CV_1_C2V_1_C3_0)'
#signal_name = 'TTTo.*'
#signal_name = 'DY_.*'
file_names = {
    'DeepFlavour' : "{}/{}_{}_tuple.root".format(folder, year, channel),
    'HHbtag' : "{}/{}_hhbtag_{}_tuple.root".format(folder, year, channel)
}


m_bb = {}
weights = {}
m_bb_hist = {}
hist_model = ROOT.RDF.TH1DModel("m_bb", "m_bb", 40, 0, 200)

for tagger, file_name in file_names.items():
    id_dfs = LoadIdFrames(file_name, id_collections = ['data', 'sample'])
    sample_match = getSampleIdMatchLine(id_dfs['sample'], signal_name)
    print(sample_match)
    df_all = ROOT.RDataFrame(channel, file_name)
    df_all = df_all.Filter("{} && is_central_es == true".format(sample_match))
    for b_index in range(1, 3):
        df_all = df_all.Define("b{}_p4".format(b_index),
                               "LorentzVectorM(b{0}_pt, b{0}_eta, b{0}_phi, b{0}_m)".format(b_index))
    df_all = df_all.Define("m_bb", "(b1_p4 + b2_p4).mass()")
    df_pure = df_all.Filter("b1_hadronFlavour == 5 && b2_hadronFlavour == 5")

    for idx,df in enumerate([df_all, df_pure]):
        pure = idx == 1
        m_bb_hist[(tagger, pure)] = df.Histo1D(hist_model, "m_bb", "weight")
        a = df.AsNumpy(["m_bb", "weight"])
        m_bb[(tagger, pure)] = np.array(a["m_bb"])
        weights[(tagger, pure)] = np.array(a["weight"])

for tagger in file_names:
    purity = (np.sum(weights[(tagger, True)]) / np.sum(weights[(tagger, False)])) * 100
    print("{}: signal purity = {:.1f}%".format(tagger, purity))
rel = (np.sum(weights[("HHbtag", True)])/np.sum(weights[("DeepFlavour", True)]) - 1)*100
print("relative increase of the pure signal yeild: {:.1f}%".format(rel))

old_conf = [79.66, 150.62]
new_conf = [85.01, 137.69]

conf_level=0.68
q = [(1 - conf_level)/2, (1 + conf_level)/2 ]
for tagger in file_names:
    purity = False
    x = m_bb[(tagger, purity)]
    w = weights[(tagger, purity)]
    conf_int = weighted_quantile(x, q, w)
    avg, std = weighted_avg_and_std(x, weights=w)
    print("{}: m_bb mean={:.2f} std={:.2f}  68% conf_interval=[{:.2f}, {:.2f}], conf_interval_width={:.2f}" \
          .format(tagger, avg, std, conf_int[0], conf_int[1], conf_int[1] - conf_int[0]))
    old_conf_cnt = np.sum(w[(x > old_conf[0]) & (x < old_conf[1])]) / np.sum(w) * 100
    new_conf_cnt = np.sum(w[(x > new_conf[0]) & (x < new_conf[1])]) / np.sum(w) * 100
    print("{}: in old conf = {:.2f}% in new conf = {:.2f}%".format(tagger, old_conf_cnt, new_conf_cnt))


canvas = ROOT.TCanvas("", "", 400, 400)
purity = False
h1 = m_bb_hist[("DeepFlavour", purity)]
h2 = m_bb_hist[("HHbtag", purity)]
h2.Draw()
h1.Draw("SAME")
h1.SetLineColor(ROOT.kBlue)
h2.SetLineColor(ROOT.kRed)
h2.SetStats(0)
h1.SetStats(0)
canvas.Draw()
