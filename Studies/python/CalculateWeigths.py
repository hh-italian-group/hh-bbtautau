import sys
import pandas as pd
import numpy as np
import numba
from scipy import interpolate

from InputsProducer import SampleType as st

ggHH_NonRes = st.ggHH_NonRes -7
VBFHH_NonRes = st.VBFHH_NonRes -7
ggHH_Res = st.ggHH_Res -7
VBFHH_Res = st.VBFHH_Res -7

@numba.njit
def calculeteNEvents(X, Z):
    mass_node = np.where(Z[:, 0, 3] == -1, Z[:, 0, 2], Z[:, 0, 3])
    mass_node_unique = np.unique(mass_node)
    counts = np.zeros((3,3,4,2,len(mass_node_unique)))

    for n in range(X.shape[0]):
        year = int(Z[n, 0, -2])#var_pos['sample_year']]
        channel = int(Z[n, 0, -1])#var_pos['channelId']]
        sample_type = int(Z[n,0,0])
        spin = int(Z[n,0, 1] > 0)
        mn = mass_node[n]
        mn_index = np.argmax(mass_node_unique == mn)

        counts[year - 2016, channel, sample_type - 6, spin,  mn_index] +=1
    weights = np.where(counts > 0, 1. / counts, 0)
    return counts, weights, mass_node_unique, mass_node

def CreateXSTable(X, Z, radion_xs, radion_br, vbf_radion_xs, graviton_xs, graviton_br, vbf_graviton_xs):

    counts, weights, mass_node_unique, mass_node = calculeteNEvents(X, Z)

    df_xs_radion = pd.read_table(radion_xs)
    df_xs_vbf_radion = pd.read_table(vbf_radion_xs)
    df_br_radion = pd.read_table(radion_br)

    df_xs_graviton = pd.read_table(graviton_xs)
    df_xs_vbf_graviton = pd.read_table(vbf_graviton_xs)
    df_br_graviton = pd.read_table(graviton_br)

    f_br_radion = interpolate.interp1d(df_br_radion['mR(GeV)'], df_br_radion['hh'], kind='cubic')
    f_xs_radion = interpolate.interp1d(df_xs_radion['mH(GeV)'], df_xs_radion['CrossSection(pb)'], kind='cubic')
    f_xs_vbf_radion = interpolate.interp1d(df_xs_vbf_radion['mH(GeV)'], df_xs_vbf_radion['CrossSection(pb)'], kind='cubic')

    f_br_graviton = interpolate.interp1d(df_br_graviton['mG(GeV)'], df_br_graviton['hh'], kind='cubic')
    f_xs_graviton = interpolate.interp1d(df_xs_graviton['mG(GeV)'], df_xs_graviton['XS(pb)'], kind='cubic')
    f_xs_vbf_graviton = interpolate.interp1d(df_xs_vbf_graviton['mG(GeV)'], df_xs_vbf_graviton['XS(pb)'], kind='cubic')

    xs_all = np.ones((4,2,len(mass_node_unique)))
    indices = mass_node_unique >= 250
    res_masses = np.copy(mass_node_unique[indices])
    if len(res_masses) > 0 and res_masses[0] == 250 :
        res_masses[0] = 251

    xs_radion = f_xs_radion(res_masses)
    xs_vbf_radion = f_xs_vbf_radion(res_masses)
    br_radion = f_br_radion(res_masses)

    xs_graviton = f_xs_graviton(res_masses)
    xs_vbf_graviton = f_xs_vbf_graviton(res_masses)
    br_graviton = f_br_graviton(res_masses)

    xs_all[ggHH_Res,0, indices] = xs_radion * br_radion
    xs_all[VBFHH_Res,0, indices] = xs_vbf_radion * br_radion

    xs_all[ggHH_Res,1, indices] = xs_graviton * br_graviton
    xs_all[VBFHH_Res,1, indices] = xs_vbf_graviton * br_graviton

    return xs_all, counts, weights, mass_node_unique, mass_node

def CreateWeights(X, Z):
    xs_all, counts, weights, mass_node_unique, mass_node = CreateXSTable(X, Z, '../config/xs_br/radion_ggF_xs.csv',
                                                                        '../config/xs_br/radion_br.csv',
                                                                        '../config/xs_br/radion_VBF_xs.csv',
                                                                        '../config/xs_br/graviton_ggF_xs.csv',
                                                                        '../config/xs_br/graviton_br.csv',
                                                                        '../config/xs_br/graviton_VBF_xs.csv')

    for sample_type in range(4):
        for spin in range(2):
            for mass in range(xs_all.shape[2]):
                weights[:,:,sample_type,spin,mass] *= xs_all[sample_type, spin, mass]

    # Radion = Graviton
    indices = mass_node_unique >= 250
    res_process = [ggHH_Res, VBFHH_Res]
    non_res_process = [ggHH_NonRes, VBFHH_NonRes]
    all_process = [res_process, non_res_process]
    for year in range(3):
        for channel in range(3):
            for prod_mode in res_process:
                radion = np.sum(weights[year, channel, prod_mode, 0, :] * counts[year, channel, prod_mode, 0, :])
                graviton = np.sum(weights[year, channel, prod_mode, 1, :] * counts[year, channel, prod_mode, 1, :])
                if radion == 0 or graviton == 0: continue
                C = radion / graviton
                weights[year, channel, prod_mode, 1, :] *= C

    # VBF = ggF
    for year in range(3):
        for channel in range(3):
            for x in all_process:
                VBF = np.sum(weights[year, channel, x[1], :, :] * counts[year, channel, x[1], :, :])
                ggF = np.sum(weights[year, channel,  x[0], :, :] * counts[year, channel,  x[0], :, :])
                if VBF == 0 or ggF == 0: continue
                E = VBF / ggF
                weights[year, channel, x[0], :, :] *= E

    # Res = Non Res
    for year in range(3):
        for channel in range(3):
            non_res = np.sum(weights[year, channel, non_res_process, :, :] * counts[year, channel, non_res_process, :, :])
            res = np.sum(weights[year, channel, res_process, :, :] * counts[year, channel, res_process, :, :])
            if non_res == 0 or res == 0: continue
            D = non_res / res
            weights[year, channel, res_process, :, :] *= D

    # eTau = muTau = tauTau
    for year in range(3):
        ch_cnt = np.zeros(3)
        for channel in range(3):
            ch_cnt[channel] = np.sum(weights[year, channel, :, :, :] * counts[year, channel, :, :, :])

        ref_ch = np.amax(ch_cnt)
        if ref_ch == 0: continue

        for channel in range(3):
            if ch_cnt[channel] != 0 :
                A = ref_ch / ch_cnt[channel]
                weights[year, channel, :, :, :] *= A

    # all years
    #https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiLUM
    luminosity = [35.92, 41.53, 59.74]
    year_cnt = np.zeros(3)
    for year in range(3):
         year_cnt[year] = np.sum(weights[year, :, :, :, :] * counts[year, :, :, :, :])

    for year in range(3):
        if year_cnt[year] != 0 :
            A = luminosity[year] / year_cnt[year]
            weights[year, :, :, :, :] *= A

    #sum all = n_evt
    weights_count = np.sum(weights * counts)
    C =  X.shape[0] / weights_count

    weights *= C

    return weights, mass_node_unique, mass_node

@numba.njit
def ConvertToVector(X,Z, weights, mass_node_unique, mass_node):
    weight_vec = np.zeros(X.shape[0])
    for n in range(X.shape[0]):
        year = int(Z[n, 0, -2])
        channel = int(Z[n, 0, -1])
        sample_type = int(Z[n,0,0])
        spin = int(Z[n,0, 1] > 0)
        mn = mass_node[n]
        mn_index = np.argmax(mass_node_unique == mn)

        weight_vec[n] = weights[year - 2016, channel, sample_type - 6, spin,  mn_index]
    return weight_vec

def CreateSampleWeigts(X,Z):
    weights, mass_node_unique, mass_node = CreateWeights(X, Z)
    return ConvertToVector(X,Z, weights, mass_node_unique, mass_node)

@numba.njit
def CrossCheckWeights(Z, X, w, g_r, res_non_res, check_channel, year, channel):
    weights, mass_node_unique, mass_node = CreateWeights(X, Z)
    w = ConvertToVector(X, Z, weights, mass_node_unique, mass_node)
    w_1 = []
    w_2 = []
    w_3 = []
    if res_non_res == True:
        for n in range(X.shape[0]):
            if (Z[n, 0, 0] == 6 or Z[n, 0, 0] == 7) and Z[n, 0, -2] == year :
                w_1.append(w[n])
            elif (Z[n, 0, 0] == 8 or Z[n, 0, 0] == 9) and Z[n, 0, -2] == year :
                w_2.append(w[n])

    elif g_r == True:
        for n in range(X.shape[0]):
            if Z[n, 0, 1] == 0 and Z[n, 0, -2] == year:
                w_1.append(w[n])
            elif Z[n, 0, 1] == 2 and Z[n, 0, -2] == year:
                w_2.append(w[n])

    elif check_channel == True:
        for n in range(X.shape[0]):
            if Z[n, 0, -1] == 0 and Z[n, 0, -2] == year:
                w_1.append(w[n])
            elif Z[n, 0, -1] == 1 and Z[n, 0, -2] == year:
                w_2.append(w[n])
            elif Z[n, 0, -1] == 2 and Z[n, 0, -2] == year:
                w_3.append(w[n])

    return np.sum(np.array(w_1)), np.sum(np.array(w_2)), np.sum(np.array(w_3))
