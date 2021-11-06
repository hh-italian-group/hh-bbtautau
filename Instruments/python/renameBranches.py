import ROOT

branches = [ 'b1_pt', 'b1_m', 'b1_resolution', 'b2_pt', 'b2_m', 'b2_resolution', 'MET_pt', 'MET_phi',
             'SVfit_valid', 'SVfit_pt', 'SVfit_eta', 'SVfit_phi', 'SVfit_m',
             'SVfit_pt_error', 'SVfit_eta_error', 'SVfit_phi_error', 'SVfit_m_error', 'SVfit_mt', 'SVfit_mt_error',
             'kinFit_convergence', 'kinFit_m', 'kinFit_chi2', 'MT2' ]

columns_to_exclude = branches

def rename(df):
    for b in branches:
        df = df.Define(b + '_orig', b)
    return df
