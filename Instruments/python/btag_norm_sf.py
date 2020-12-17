import ROOT
import argparse
import json
from numpy import *



def add_dic(r_factor, r_factor_name, r_factors):
    print ("r factor related to %s is %.4f " % (r_factor_name ,r_factor))
    if isnan(r_factor):
        r_factor = 1
    r_factors[r_factor_name]=r_factor

def evaluate_r(input_dir,channel, year, unc_source, r_factors):
    file = input_dir+year+"_"+channel+"_"+unc_source+".root"
    d = ROOT.RDataFrame(channel, file)
    not_data = '! is_data'
    d=d.Filter(not_data)
    if(unc_source == "Central"):
        print "\n ******** Central ********* \n"
        w_central_before = d.Define( "num", "weight").Sum('num')
        w_central_after = d.Define( "den", "weight*weight_btag_IterativeFit").Sum('den')
        r_factor=w_central_before.GetValue()/w_central_after.GetValue()
        add_dic(r_factor, "Central",r_factors)

        w_btag_lf_Up_before = d.Define( "num_btag_lf_Up", "weight*unc_btag_lf_Up").Sum('num_btag_lf_Up')
        w_btag_lf_Up_after = d.Define( "den_btag_lf_Up", "weight*unc_btag_lf_Up*weight_btag_IterativeFit*unc_btag_lf_Up").Sum('den_btag_lf_Up')
        r_factor_btag_lf_Up = w_btag_lf_Up_before.GetValue()/w_btag_lf_Up_after.GetValue()
        add_dic(r_factor_btag_lf_Up, "btag_lf_Up",r_factors)

        w_btag_lf_Down_before = d.Define( "num_btag_lf_Down", "weight*unc_btag_lf_Down").Sum('num_btag_lf_Down')
        w_btag_lf_Down_after = d.Define( "den_btag_lf_Down", "weight*unc_btag_lf_Down*weight_btag_IterativeFit*unc_btag_lf_Down").Sum('den_btag_lf_Down')
        r_factor_btag_lf_Down = w_btag_lf_Down_before.GetValue()/w_btag_lf_Down_after.GetValue()
        add_dic(r_factor_btag_lf_Down, "btag_lf_Down",r_factors)

        w_btag_hf_Up_before = d.Define( "num_btag_hf_Up", "weight*unc_btag_hf_Up").Sum('num_btag_hf_Up')
        w_btag_hf_Up_after = d.Define( "den_btag_hf_Up", "weight*unc_btag_hf_Up*weight_btag_IterativeFit*unc_btag_hf_Up").Sum('den_btag_hf_Up')
        r_factor_btag_hf_Up = w_btag_hf_Up_before.GetValue()/w_btag_hf_Up_after.GetValue()
        add_dic(r_factor_btag_hf_Up, "btag_hf_Up",r_factors)

        w_btag_hf_Down_before = d.Define( "num_btag_hf_Down", "weight*unc_btag_hf_Down").Sum('num_btag_hf_Down')
        w_btag_hf_Down_after = d.Define( "den_btag_hf_Down", "weight*unc_btag_hf_Down*weight_btag_IterativeFit*unc_btag_hf_Down").Sum('den_btag_hf_Down')
        r_factor_btag_hf_Down = w_btag_hf_Down_before.GetValue()/w_btag_hf_Down_after.GetValue()
        add_dic(r_factor_btag_hf_Down, "btag_hf_Down",r_factors)

        w_btag_hfstats1_Up_before = d.Define( "num_btag_hfstats1_Up", "weight*unc_btag_hfstats1_Up").Sum('num_btag_hfstats1_Up')
        w_btag_hfstats1_Up_after = d.Define( "den_btag_hfstats1_Up", "weight*unc_btag_hfstats1_Up*weight_btag_IterativeFit*unc_btag_hfstats1_Up").Sum('den_btag_hfstats1_Up')
        r_factor_btag_hfstats1_Up = w_btag_hfstats1_Up_before.GetValue()/w_btag_hfstats1_Up_after.GetValue()
        add_dic(r_factor_btag_hfstats1_Up, "btag_hfstats1_Up",r_factors)

        w_btag_hfstats1_Down_before = d.Define( "num_btag_hfstats1_Down", "weight*unc_btag_hfstats1_Down").Sum('num_btag_hfstats1_Down')
        w_btag_hfstats1_Down_after = d.Define( "den_btag_hfstats1_Down", "weight*unc_btag_hfstats1_Down*weight_btag_IterativeFit*unc_btag_hfstats1_Down").Sum('den_btag_hfstats1_Down')
        r_factor_btag_hfstats1_Down = w_btag_hfstats1_Down_before.GetValue()/w_btag_hfstats1_Down_after.GetValue()
        add_dic(r_factor_btag_hfstats1_Down, "btag_hfstats1_Down",r_factors)

        w_btag_hfstats2_Up_before = d.Define( "num_btag_hfstats2_Up", "weight*unc_btag_hfstats2_Up").Sum('num_btag_hfstats2_Up')
        w_btag_hfstats2_Up_after = d.Define( "den_btag_hfstats2_Up", "weight*unc_btag_hfstats2_Up*weight_btag_IterativeFit*unc_btag_hfstats2_Up").Sum('den_btag_hfstats2_Up')
        r_factor_btag_hfstats2_Up = w_btag_hfstats2_Up_before.GetValue()/w_btag_hfstats2_Up_after.GetValue()
        add_dic(r_factor_btag_hfstats2_Up, "btag_hfstats2_Up",r_factors)

        w_btag_hfstats2_Down_before = d.Define( "num_btag_hfstats2_Down", "weight*unc_btag_hfstats2_Down").Sum('num_btag_hfstats2_Down')
        w_btag_hfstats2_Down_after = d.Define( "den_btag_hfstats2_Down", "weight*unc_btag_hfstats2_Down*weight_btag_IterativeFit*unc_btag_hfstats2_Down").Sum('den_btag_hfstats2_Down')
        r_factor_btag_hfstats2_Down = w_btag_hfstats2_Down_before.GetValue()/w_btag_hfstats2_Down_after.GetValue()
        add_dic(r_factor_btag_hfstats2_Down, "btag_hfstats2_Down",r_factors)

        w_btag_lfstats1_Up_before = d.Define( "num_btag_lfstats1_Up", "weight*unc_btag_lfstats1_Up").Sum('num_btag_lfstats1_Up')
        w_btag_lfstats1_Up_after = d.Define( "den_btag_lfstats1_Up", "weight*unc_btag_lfstats1_Up*weight_btag_IterativeFit*unc_btag_lfstats1_Up").Sum('den_btag_lfstats1_Up')
        r_factor_btag_lfstats1_Up = w_btag_lfstats1_Up_before.GetValue()/w_btag_lfstats1_Up_after.GetValue()
        add_dic(r_factor_btag_lfstats1_Up, "btag_lfstats1_Up",r_factors)

        w_btag_lfstats1_Down_before = d.Define( "num_btag_lfstats1_Down", "weight*unc_btag_lfstats1_Down").Sum('num_btag_lfstats1_Down')
        w_btag_lfstats1_Down_after = d.Define( "den_btag_lfstats1_Down", "weight*unc_btag_lfstats1_Down*weight_btag_IterativeFit*unc_btag_lfstats1_Down").Sum('den_btag_lfstats1_Down')
        r_factor_btag_lfstats1_Down = w_btag_lfstats1_Down_before.GetValue()/w_btag_lfstats1_Down_after.GetValue()
        add_dic(r_factor_btag_lfstats1_Down, "btag_lfstats1_Down",r_factors)

        w_btag_lfstats2_Up_before = d.Define( "num_btag_lfstats2_Up", "weight*unc_btag_lfstats2_Up").Sum('num_btag_lfstats2_Up')
        w_btag_lfstats2_Up_after = d.Define( "den_btag_lfstats2_Up", "weight*unc_btag_lfstats2_Up*weight_btag_IterativeFit*unc_btag_lfstats2_Up").Sum('den_btag_lfstats2_Up')
        r_factor_btag_lfstats2_Up = w_btag_lfstats2_Up_before.GetValue()/w_btag_lfstats2_Up_after.GetValue()
        add_dic(r_factor_btag_lfstats2_Up, "btag_lfstats2_Up",r_factors)

        w_btag_lfstats2_Down_before = d.Define( "num_btag_lfstats2_Down", "weight*unc_btag_lfstats2_Down").Sum('num_btag_lfstats2_Down')
        w_btag_lfstats2_Down_after = d.Define( "den_btag_lfstats2_Down", "weight*unc_btag_lfstats2_Down*weight_btag_IterativeFit*unc_btag_lfstats2_Down").Sum('den_btag_lfstats2_Down')
        r_factor_btag_lfstats2_Down = w_btag_lfstats2_Down_before.GetValue()/w_btag_lfstats2_Down_after.GetValue()
        add_dic(r_factor_btag_lfstats2_Down, "btag_lfstats2_Down",r_factors)

        w_btag_cferr1_Up_before = d.Define( "num_btag_cferr1_Up", "weight*unc_btag_cferr1_Up").Sum('num_btag_cferr1_Up')
        w_btag_cferr1_Up_after = d.Define( "den_btag_cferr1_Up", "weight*unc_btag_cferr1_Up*weight_btag_IterativeFit*unc_btag_cferr1_Up").Sum('den_btag_cferr1_Up')
        r_factor_btag_cferr1_Up = w_btag_cferr1_Up_before.GetValue()/w_btag_cferr1_Up_after.GetValue()
        add_dic(r_factor_btag_cferr1_Up, "btag_cferr1_Up",r_factors)

        w_btag_cferr1_Down_before = d.Define( "num_btag_cferr1_Down", "weight*unc_btag_cferr1_Down").Sum('num_btag_cferr1_Down')
        w_btag_cferr1_Down_after = d.Define( "den_btag_cferr1_Down", "weight*unc_btag_cferr1_Down*weight_btag_IterativeFit*unc_btag_cferr1_Down").Sum('den_btag_cferr1_Down')
        r_factor_btag_cferr1_Down = w_btag_cferr1_Down_before.GetValue()/w_btag_cferr1_Down_after.GetValue()
        add_dic(r_factor_btag_cferr1_Down, "btag_cferr1_Down",r_factors)

        w_btag_cferr2_Up_before = d.Define( "num_btag_cferr2_Up", "weight*unc_btag_cferr2_Up").Sum('num_btag_cferr2_Up')
        w_btag_cferr2_Up_after = d.Define( "den_btag_cferr2_Up", "weight*unc_btag_cferr2_Up*weight_btag_IterativeFit*unc_btag_cferr2_Up").Sum('den_btag_cferr2_Up')
        r_factor_btag_cferr2_Up = w_btag_cferr2_Up_before.GetValue()/w_btag_cferr2_Up_after.GetValue()
        add_dic(r_factor_btag_cferr2_Up, "btag_cferr2_Up",r_factors)

        w_btag_cferr2_Down_before = d.Define( "num_btag_cferr2_Down", "weight*unc_btag_cferr2_Down").Sum('num_btag_cferr2_Down')
        w_btag_cferr2_Down_after = d.Define( "den_btag_cferr2_Down", "weight*unc_btag_cferr2_Down*weight_btag_IterativeFit*unc_btag_cferr2_Down").Sum('den_btag_cferr2_Down')
        r_factor_btag_cferr2_Down = w_btag_cferr2_Down_before.GetValue()/w_btag_cferr2_Down_after.GetValue()
        add_dic(r_factor_btag_cferr2_Down, "btag_cferr2_Down",r_factors)

    if(unc_source=="JES"):
        print "\n ******** JES ********* \n"
        w_before_JetReduced_Absolute_Up= d.Filter("unc_source==44").Filter("unc_scale==+1").Define("num_JetReduced_Absolute_Up", "weight").Sum('num_JetReduced_Absolute_Up')
        w_after_JetReduced_Absolute_Up= d.Filter("unc_source==44").Filter("unc_scale==+1").Define( "den_JetReduced_Absolute_Up", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_Absolute_Up')
        r_factor_JetReduced_Absolute_Up=w_before_JetReduced_Absolute_Up.GetValue()/w_after_JetReduced_Absolute_Up.GetValue()
        add_dic(r_factor_JetReduced_Absolute_Up, "JetReduced_Absolute_Up",r_factors)

        w_before_JetReduced_Absolute_Down = d.Filter("unc_source==44").Filter("unc_scale==-1").Define("num_JetReduced_Absolute_Down", "weight").Sum('num_JetReduced_Absolute_Down')
        w_after_JetReduced_Absolute_Down = d.Filter("unc_source==44").Filter("unc_scale==-1").Define( "den_JetReduced_Absolute_Down", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_Absolute_Down')
        r_factor_JetReduced_Absolute_Down=w_before_JetReduced_Absolute_Down.GetValue()/w_after_JetReduced_Absolute_Down.GetValue()
        add_dic(r_factor_JetReduced_Absolute_Down, "JetReduced_Absolute_Down",r_factors)


        w_before_JetReduced_Absolute_year_Up= d.Filter("unc_source==45").Filter("unc_scale==+1").Define("num_JetReduced_Absolute_year_Up", "weight").Sum('num_JetReduced_Absolute_year_Up')
        w_after_JetReduced_Absolute_year_Up= d.Filter("unc_source==45").Filter("unc_scale==+1").Define( "den_JetReduced_Absolute_year_Up", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_Absolute_year_Up')
        r_factor_JetReduced_Absolute_year_Up=w_before_JetReduced_Absolute_year_Up.GetValue()/w_after_JetReduced_Absolute_year_Up.GetValue()
        add_dic(r_factor_JetReduced_Absolute_year_Up, "JetReduced_Absolute_year_Up",r_factors)

        w_before_JetReduced_Absolute_year_Down = d.Filter("unc_source==45").Filter("unc_scale==-1").Define("num_JetReduced_Absolute_year_Down", "weight").Sum('num_JetReduced_Absolute_year_Down')
        w_after_JetReduced_Absolute_year_Down = d.Filter("unc_source==45").Filter("unc_scale==-1").Define( "den_JetReduced_Absolute_year_Down", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_Absolute_year_Down')
        r_factor_JetReduced_Absolute_year_Down=w_before_JetReduced_Absolute_year_Down.GetValue()/w_after_JetReduced_Absolute_year_Down.GetValue()
        add_dic(r_factor_JetReduced_Absolute_year_Down, "JetReduced_Absolute_year_Down",r_factors)


        w_before_JetReduced_BBEC1_Up= d.Filter("unc_source==46").Filter("unc_scale==+1").Define("num_JetReduced_BBEC1_Up", "weight").Sum('num_JetReduced_BBEC1_Up')
        w_after_JetReduced_BBEC1_Up= d.Filter("unc_source==46").Filter("unc_scale==+1").Define( "den_JetReduced_BBEC1_Up", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_BBEC1_Up')
        r_factor_JetReduced_BBEC1_Up=w_before_JetReduced_BBEC1_Up.GetValue()/w_after_JetReduced_BBEC1_Up.GetValue()
        add_dic(r_factor_JetReduced_BBEC1_Up, "JetReduced_BBEC1_Up",r_factors)

        w_before_JetReduced_BBEC1_Down = d.Filter("unc_source==46").Filter("unc_scale==-1").Define("num_JetReduced_BBEC1_Down", "weight").Sum('num_JetReduced_BBEC1_Down')
        w_after_JetReduced_BBEC1_Down = d.Filter("unc_source==46").Filter("unc_scale==-1").Define( "den_JetReduced_BBEC1_Down", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_BBEC1_Down')
        r_factor_JetReduced_BBEC1_Down=w_before_JetReduced_BBEC1_Down.GetValue()/w_after_JetReduced_BBEC1_Down.GetValue()
        add_dic(r_factor_JetReduced_BBEC1_Down, "JetReduced_BBEC1_Down",r_factors)


        w_before_JetReduced_BBEC1_year_Up= d.Filter("unc_source==47").Filter("unc_scale==+1").Define("num_JetReduced_BBEC1_year_Up", "weight").Sum('num_JetReduced_BBEC1_year_Up')
        w_after_JetReduced_BBEC1_year_Up= d.Filter("unc_source==47").Filter("unc_scale==+1").Define( "den_JetReduced_BBEC1_year_Up", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_BBEC1_year_Up')
        r_factor_JetReduced_BBEC1_year_Up=w_before_JetReduced_BBEC1_year_Up.GetValue()/w_after_JetReduced_BBEC1_year_Up.GetValue()
        add_dic(r_factor_JetReduced_BBEC1_year_Up, "JetReduced_BBEC1_year_Up",r_factors)

        w_before_JetReduced_BBEC1_year_Down = d.Filter("unc_source==47").Filter("unc_scale==-1").Define("num_JetReduced_BBEC1_year_Down", "weight").Sum('num_JetReduced_BBEC1_year_Down')
        w_after_JetReduced_BBEC1_year_Down = d.Filter("unc_source==47").Filter("unc_scale==-1").Define( "den_JetReduced_BBEC1_year_Down", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_BBEC1_year_Down')
        r_factor_JetReduced_BBEC1_year_Down=w_before_JetReduced_BBEC1_year_Down.GetValue()/w_after_JetReduced_BBEC1_year_Down.GetValue()
        add_dic(r_factor_JetReduced_BBEC1_year_Down, "JetReduced_BBEC1_year_Down",r_factors)



        w_before_JetReduced_EC2_Up= d.Filter("unc_source==48").Filter("unc_scale==+1").Define("num_JetReduced_EC2_Up", "weight").Sum('num_JetReduced_EC2_Up')
        w_after_JetReduced_EC2_Up= d.Filter("unc_source==48").Filter("unc_scale==+1").Define( "den_JetReduced_EC2_Up", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_EC2_Up')
        r_factor_JetReduced_EC2_Up=w_before_JetReduced_EC2_Up.GetValue()/w_after_JetReduced_EC2_Up.GetValue()
        add_dic(r_factor_JetReduced_EC2_Up, "JetReduced_EC2_Up",r_factors)

        w_before_JetReduced_EC2_Down = d.Filter("unc_source==48").Filter("unc_scale==-1").Define("num_JetReduced_EC2_Down", "weight").Sum('num_JetReduced_EC2_Down')
        w_after_JetReduced_EC2_Down = d.Filter("unc_source==48").Filter("unc_scale==-1").Define( "den_JetReduced_EC2_Down", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_EC2_Down')
        r_factor_JetReduced_EC2_Down=w_before_JetReduced_EC2_Down.GetValue()/w_after_JetReduced_EC2_Down.GetValue()
        add_dic(r_factor_JetReduced_EC2_Down, "JetReduced_EC2_Down",r_factors)


        w_before_JetReduced_EC2_year_Up= d.Filter("unc_source==49").Filter("unc_scale==+1").Define("num_JetReduced_EC2_year_Up", "weight").Sum('num_JetReduced_EC2_year_Up')
        w_after_JetReduced_EC2_year_Up= d.Filter("unc_source==49").Filter("unc_scale==+1").Define( "den_JetReduced_EC2_year_Up", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_EC2_year_Up')
        r_factor_JetReduced_EC2_year_Up=w_before_JetReduced_EC2_year_Up.GetValue()/w_after_JetReduced_EC2_year_Up.GetValue()
        add_dic(r_factor_JetReduced_EC2_year_Up, "JetReduced_EC2_year_Up",r_factors)

        w_before_JetReduced_EC2_year_Down = d.Filter("unc_source==49").Filter("unc_scale==-1").Define("num_JetReduced_EC2_year_Down", "weight").Sum('num_JetReduced_EC2_year_Down')
        w_after_JetReduced_EC2_year_Down = d.Filter("unc_source==49").Filter("unc_scale==-1").Define( "den_JetReduced_EC2_year_Down", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_EC2_year_Down')
        r_factor_JetReduced_EC2_year_Down=w_before_JetReduced_EC2_year_Down.GetValue()/w_after_JetReduced_EC2_year_Down.GetValue()
        add_dic(r_factor_JetReduced_EC2_year_Down, "JetReduced_EC2_year_Down",r_factors)


        w_before_JetReduced_FlavorQCD_Up= d.Filter("unc_source==50").Filter("unc_scale==+1").Define("num_JetReduced_FlavorQCD_Up", "weight").Sum('num_JetReduced_FlavorQCD_Up')
        w_after_JetReduced_FlavorQCD_Up= d.Filter("unc_source==50").Filter("unc_scale==+1").Define( "den_JetReduced_FlavorQCD_Up", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_FlavorQCD_Up')
        r_factor_JetReduced_FlavorQCD_Up=w_before_JetReduced_FlavorQCD_Up.GetValue()/w_after_JetReduced_FlavorQCD_Up.GetValue()
        add_dic(r_factor_JetReduced_FlavorQCD_Up, "JetReduced_FlavorQCD_Up",r_factors)

        w_before_JetReduced_FlavorQCD_Down = d.Filter("unc_source==50").Filter("unc_scale==-1").Define("num_JetReduced_FlavorQCD_Down", "weight").Sum('num_JetReduced_FlavorQCD_Down')
        w_after_JetReduced_FlavorQCD_Down = d.Filter("unc_source==50").Filter("unc_scale==-1").Define( "den_JetReduced_FlavorQCD_Down", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_FlavorQCD_Down')
        r_factor_JetReduced_FlavorQCD_Down=w_before_JetReduced_FlavorQCD_Down.GetValue()/w_after_JetReduced_FlavorQCD_Down.GetValue()
        add_dic(r_factor_JetReduced_FlavorQCD_Down, "JetReduced_FlavorQCD_Down",r_factors)


        w_before_JetReduced_HF_Up= d.Filter("unc_source==51").Filter("unc_scale==+1").Define("num_JetReduced_HF_Up", "weight").Sum('num_JetReduced_HF_Up')
        w_after_JetReduced_HF_Up= d.Filter("unc_source==51").Filter("unc_scale==+1").Define( "den_JetReduced_HF_Up", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_HF_Up')
        r_factor_JetReduced_HF_Up=w_before_JetReduced_HF_Up.GetValue()/w_after_JetReduced_HF_Up.GetValue()
        add_dic(r_factor_JetReduced_HF_Up, "JetReduced_HF_Up",r_factors)

        w_before_JetReduced_HF_Down = d.Filter("unc_source==51").Filter("unc_scale==-1").Define("num_JetReduced_HF_Down", "weight").Sum('num_JetReduced_HF_Down')
        w_after_JetReduced_HF_Down = d.Filter("unc_source==51").Filter("unc_scale==-1").Define( "den_JetReduced_HF_Down", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_HF_Down')
        r_factor_JetReduced_HF_Down=w_before_JetReduced_HF_Down.GetValue()/w_after_JetReduced_HF_Down.GetValue()
        add_dic(r_factor_JetReduced_HF_Down, "JetReduced_HF_Down",r_factors)


        w_before_JetReduced_HF_year_Up= d.Filter("unc_source==52").Filter("unc_scale==+1").Define("num_JetReduced_HF_year_Up", "weight").Sum('num_JetReduced_HF_year_Up')
        w_after_JetReduced_HF_year_Up= d.Filter("unc_source==52").Filter("unc_scale==+1").Define( "den_JetReduced_HF_year_Up", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_HF_year_Up')
        r_factor_JetReduced_HF_year_Up=w_before_JetReduced_HF_year_Up.GetValue()/w_after_JetReduced_HF_year_Up.GetValue()
        add_dic(r_factor_JetReduced_HF_year_Up, "JetReduced_HF_year_Up",r_factors)

        w_before_JetReduced_HF_year_Down = d.Filter("unc_source==52").Filter("unc_scale==-1").Define("num_JetReduced_HF_year_Down", "weight").Sum('num_JetReduced_HF_year_Down')
        w_after_JetReduced_HF_year_Down = d.Filter("unc_source==52").Filter("unc_scale==-1").Define( "den_JetReduced_HF_year_Down", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_HF_year_Down')
        r_factor_JetReduced_HF_year_Down=w_before_JetReduced_HF_year_Down.GetValue()/w_after_JetReduced_HF_year_Down.GetValue()
        add_dic(r_factor_JetReduced_HF_year_Down, "JetReduced_HF_year_Down",r_factors)


        w_before_JetReduced_RelativeBal_Up= d.Filter("unc_source==53").Filter("unc_scale==+1").Define("num_JetReduced_RelativeBal_Up", "weight").Sum('num_JetReduced_RelativeBal_Up')
        w_after_JetReduced_RelativeBal_Up= d.Filter("unc_source==53").Filter("unc_scale==+1").Define( "den_JetReduced_RelativeBal_Up", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_RelativeBal_Up')
        r_factor_JetReduced_RelativeBal_Up=w_before_JetReduced_RelativeBal_Up.GetValue()/w_after_JetReduced_RelativeBal_Up.GetValue()
        add_dic(r_factor_JetReduced_RelativeBal_Up, "JetReduced_RelativeBal_Up",r_factors)

        w_before_JetReduced_RelativeBal_Down = d.Filter("unc_source==53").Filter("unc_scale==-1").Define("num_JetReduced_RelativeBal_Down", "weight").Sum('num_JetReduced_RelativeBal_Down')
        w_after_JetReduced_RelativeBal_Down = d.Filter("unc_source==53").Filter("unc_scale==-1").Define( "den_JetReduced_RelativeBal_Down", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_RelativeBal_Down')
        r_factor_JetReduced_RelativeBal_Down=w_before_JetReduced_RelativeBal_Down.GetValue()/w_after_JetReduced_RelativeBal_Down.GetValue()
        add_dic(r_factor_JetReduced_RelativeBal_Down, "JetReduced_RelativeBal_Down",r_factors)


        w_before_JetReduced_RelativeSample_year_Up= d.Filter("unc_source==54").Filter("unc_scale==+1").Define("num_JetReduced_RelativeSample_year_Up", "weight").Sum('num_JetReduced_RelativeSample_year_Up')
        w_after_JetReduced_RelativeSample_year_Up= d.Filter("unc_source==54").Filter("unc_scale==+1").Define( "den_JetReduced_RelativeSample_year_Up", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_RelativeSample_year_Up')
        r_factor_JetReduced_RelativeSample_year_Up=w_before_JetReduced_RelativeSample_year_Up.GetValue()/w_after_JetReduced_RelativeSample_year_Up.GetValue()
        add_dic(r_factor_JetReduced_RelativeSample_year_Up, "JetReduced_RelativeSample_year_Up",r_factors)

        w_before_JetReduced_RelativeSample_year_Down = d.Filter("unc_source==54").Filter("unc_scale==-1").Define("num_JetReduced_RelativeSample_year_Down", "weight").Sum('num_JetReduced_RelativeSample_year_Down')
        w_after_JetReduced_RelativeSample_year_Down = d.Filter("unc_source==54").Filter("unc_scale==-1").Define( "den_JetReduced_RelativeSample_year_Down", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_RelativeSample_year_Down')
        r_factor_JetReduced_RelativeSample_year_Down=w_before_JetReduced_RelativeSample_year_Down.GetValue()/w_after_JetReduced_RelativeSample_year_Down.GetValue()
        add_dic(r_factor_JetReduced_RelativeSample_year_Down, "JetReduced_RelativeSample_year_Down",r_factors)

        w_before_JetReduced_Total_Up= d.Filter("unc_source==55").Filter("unc_scale==+1").Define("num_JetReduced_Total_Up", "weight").Sum('num_JetReduced_Total_Up')
        w_after_JetReduced_Total_Up= d.Filter("unc_source==55").Filter("unc_scale==+1").Define( "den_JetReduced_Total_Up", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_Total_Up')
        r_factor_JetReduced_Total_Up=w_before_JetReduced_Total_Up.GetValue()/w_after_JetReduced_Total_Up.GetValue()
        add_dic(r_factor_JetReduced_Total_Up, "JetReduced_Total_Up",r_factors)

        w_before_JetReduced_Total_Down = d.Filter("unc_source==55").Filter("unc_scale==-1").Define("num_JetReduced_Total_Down", "weight").Sum('num_JetReduced_Total_Down')
        w_after_JetReduced_Total_Down = d.Filter("unc_source==55").Filter("unc_scale==-1").Define( "den_JetReduced_Total_Down", "weight*weight_btag_IterativeFit").Sum('den_JetReduced_Total_Down')
        r_factor_JetReduced_Total_Down=w_before_JetReduced_Total_Down.GetValue()/w_after_JetReduced_Total_Down.GetValue()
        add_dic(r_factor_JetReduced_Total_Down, "JetReduced_Total_Down",r_factors)


        w_before_JetReduced_Total_withJES_Up= d.Filter("unc_source==55").Filter("unc_scale==+1").Define("num_JetReduced_Total_withJES_Up", "weight").Sum('num_JetReduced_Total_withJES_Up')
        w_after_JetReduced_Total_withJES_Up= d.Filter("unc_source==55").Filter("unc_scale==+1").Define( "den_JetReduced_Total_withJES_Up", "weight*weight_btag_IterativeFit_withJES").Sum('den_JetReduced_Total_withJES_Up')
        r_factor_JetReduced_Total_withJES_Up=w_before_JetReduced_Total_withJES_Up.GetValue()/w_after_JetReduced_Total_withJES_Up.GetValue()
        add_dic(r_factor_JetReduced_Total_withJES_Up, "JetReduced_Total_withJES_Up",r_factors)

        w_before_JetReduced_Total_withJES_Down = d.Filter("unc_source==55").Filter("unc_scale==-1").Define("num_JetReduced_Total_withJES_Down", "weight").Sum('num_JetReduced_Total_withJES_Down')
        w_after_JetReduced_Total_withJES_Down = d.Filter("unc_source==55").Filter("unc_scale==-1").Define( "den_JetReduced_Total_withJES_Down", "weight*weight_btag_IterativeFit_withJES").Sum('den_JetReduced_Total_withJES_Down')
        r_factor_JetReduced_Total_withJES_Down=w_before_JetReduced_Total_withJES_Down.GetValue()/w_after_JetReduced_Total_withJES_Down.GetValue()
        add_dic(r_factor_JetReduced_Total_withJES_Down, "JetReduced_Total_withJES_Down",r_factors)

    if(unc_source=="LES"):

        print "\n ******** LES ********* \n"

        w_before_TauES_Up= d.Filter("unc_source==1").Filter("unc_scale==+1").Define("num_TauES_Up", "weight").Sum('num_TauES_Up')
        w_after_TauES_Up= d.Filter("unc_source==1").Filter("unc_scale==+1").Define( "den_TauES_Up", "weight*weight_btag_IterativeFit").Sum('den_TauES_Up')
        r_factor_TauES_Up=w_before_TauES_Up.GetValue()/w_after_TauES_Up.GetValue()
        add_dic(r_factor_TauES_Up, "TauES_Up",r_factors)

        w_before_TauES_Down = d.Filter("unc_source==1").Filter("unc_scale==-1").Define("num_TauES_Down", "weight").Sum('num_TauES_Down')
        w_after_TauES_Down = d.Filter("unc_source==1").Filter("unc_scale==-1").Define( "den_TauES_Down", "weight*weight_btag_IterativeFit").Sum('den_TauES_Down')
        r_factor_TauES_Down=w_before_TauES_Down.GetValue()/w_after_TauES_Down.GetValue()
        add_dic(r_factor_TauES_Down, "TauES_Down",r_factors)


        w_before_TauES_DM0_Up= d.Filter("unc_source==96").Filter("unc_scale==+1").Define("num_TauES_DM0_Up", "weight").Sum('num_TauES_DM0_Up')
        w_after_TauES_DM0_Up= d.Filter("unc_source==96").Filter("unc_scale==+1").Define( "den_TauES_DM0_Up", "weight*weight_btag_IterativeFit").Sum('den_TauES_DM0_Up')
        r_factor_TauES_DM0_Up=w_before_TauES_DM0_Up.GetValue()/w_after_TauES_DM0_Up.GetValue()
        add_dic(r_factor_TauES_DM0_Up, "TauES_DM0_Up",r_factors)

        w_before_TauES_DM0_Down = d.Filter("unc_source==96").Filter("unc_scale==-1").Define("num_TauES_DM0_Down", "weight").Sum('num_TauES_DM0_Down')
        w_after_TauES_DM0_Down = d.Filter("unc_source==96").Filter("unc_scale==-1").Define( "den_TauES_DM0_Down", "weight*weight_btag_IterativeFit").Sum('den_TauES_DM0_Down')
        r_factor_TauES_DM0_Down=w_before_TauES_DM0_Down.GetValue()/w_after_TauES_DM0_Down.GetValue()
        add_dic(r_factor_TauES_DM0_Down, "TauES_DM0_Down",r_factors)


        w_before_TauES_DM1_Up= d.Filter("unc_source==97").Filter("unc_scale==+1").Define("num_TauES_DM1_Up", "weight").Sum('num_TauES_DM1_Up')
        w_after_TauES_DM1_Up= d.Filter("unc_source==97").Filter("unc_scale==+1").Define( "den_TauES_DM1_Up", "weight*weight_btag_IterativeFit").Sum('den_TauES_DM1_Up')
        r_factor_TauES_DM1_Up=w_before_TauES_DM1_Up.GetValue()/w_after_TauES_DM1_Up.GetValue()
        add_dic(r_factor_TauES_DM1_Up, "TauES_DM1_Up",r_factors)

        w_before_TauES_DM1_Down = d.Filter("unc_source==97").Filter("unc_scale==-1").Define("num_TauES_DM1_Down", "weight").Sum('num_TauES_DM1_Down')
        w_after_TauES_DM1_Down = d.Filter("unc_source==97").Filter("unc_scale==-1").Define( "den_TauES_DM1_Down", "weight*weight_btag_IterativeFit").Sum('den_TauES_DM1_Down')
        r_factor_TauES_DM1_Down=w_before_TauES_DM1_Down.GetValue()/w_after_TauES_DM1_Down.GetValue()
        add_dic(r_factor_TauES_DM1_Down, "TauES_DM1_Down",r_factors)

        w_before_TauES_DM10_Up= d.Filter("unc_source==98").Filter("unc_scale==+1").Define("num_TauES_DM10_Up", "weight").Sum('num_TauES_DM10_Up')
        w_after_TauES_DM10_Up= d.Filter("unc_source==98").Filter("unc_scale==+1").Define( "den_TauES_DM10_Up", "weight*weight_btag_IterativeFit").Sum('den_TauES_DM10_Up')
        r_factor_TauES_DM10_Up=w_before_TauES_DM10_Up.GetValue()/w_after_TauES_DM10_Up.GetValue()
        add_dic(r_factor_TauES_DM10_Up, "TauES_DM10_Up",r_factors)

        w_before_TauES_DM10_Down = d.Filter("unc_source==98").Filter("unc_scale==-1").Define("num_TauES_DM10_Down", "weight").Sum('num_TauES_DM10_Down')
        w_after_TauES_DM10_Down = d.Filter("unc_source==98").Filter("unc_scale==-1").Define( "den_TauES_DM10_Down", "weight*weight_btag_IterativeFit").Sum('den_TauES_DM10_Down')
        r_factor_TauES_DM10_Down=w_before_TauES_DM10_Down.GetValue()/w_after_TauES_DM10_Down.GetValue()
        add_dic(r_factor_TauES_DM10_Down, "TauES_DM10_Down",r_factors)

        w_before_TauES_DM11_Up= d.Filter("unc_source==99").Filter("unc_scale==+1").Define("num_TauES_DM11_Up", "weight").Sum('num_TauES_DM11_Up')
        w_after_TauES_DM11_Up= d.Filter("unc_source==99").Filter("unc_scale==+1").Define( "den_TauES_DM11_Up", "weight*weight_btag_IterativeFit").Sum('den_TauES_DM11_Up')
        r_factor_TauES_DM11_Up=w_before_TauES_DM11_Up.GetValue()/w_after_TauES_DM11_Up.GetValue()
        add_dic(r_factor_TauES_DM11_Up, "TauES_DM11_Up",r_factors)

        w_before_TauES_DM11_Down = d.Filter("unc_source==99").Filter("unc_scale==-1").Define("num_TauES_DM11_Down", "weight").Sum('num_TauES_DM11_Down')
        w_after_TauES_DM11_Down = d.Filter("unc_source==99").Filter("unc_scale==-1").Define( "den_TauES_DM11_Down", "weight*weight_btag_IterativeFit").Sum('den_TauES_DM11_Down')
        r_factor_TauES_DM11_Down=w_before_TauES_DM11_Down.GetValue()/w_after_TauES_DM11_Down.GetValue()
        add_dic(r_factor_TauES_DM11_Down, "TauES_DM11_Down",r_factors)

        w_before_EleFakingTauES_DM0_Up= d.Filter("unc_source==108").Filter("unc_scale==+1").Define("num_EleFakingTauES_DM0_Up", "weight").Sum('num_EleFakingTauES_DM0_Up')
        w_after_EleFakingTauES_DM0_Up= d.Filter("unc_source==108").Filter("unc_scale==+1").Define( "den_EleFakingTauES_DM0_Up", "weight*weight_btag_IterativeFit").Sum('den_EleFakingTauES_DM0_Up')
        r_factor_EleFakingTauES_DM0_Up=w_before_EleFakingTauES_DM0_Up.GetValue()/w_after_EleFakingTauES_DM0_Up.GetValue()
        add_dic(r_factor_EleFakingTauES_DM0_Up, "EleFakingTauES_DM0_Up",r_factors)

        w_before_EleFakingTauES_DM0_Down = d.Filter("unc_source==108").Filter("unc_scale==-1").Define("num_EleFakingTauES_DM0_Down", "weight").Sum('num_EleFakingTauES_DM0_Down')
        w_after_EleFakingTauES_DM0_Down = d.Filter("unc_source==108").Filter("unc_scale==-1").Define( "den_EleFakingTauES_DM0_Down", "weight*weight_btag_IterativeFit").Sum('den_EleFakingTauES_DM0_Down')
        r_factor_EleFakingTauES_DM0_Down=w_before_EleFakingTauES_DM0_Down.GetValue()/w_after_EleFakingTauES_DM0_Down.GetValue()
        add_dic(r_factor_EleFakingTauES_DM0_Down, "EleFakingTauES_DM0_Down",r_factors)

        w_before_EleFakingTauES_DM1_Up= d.Filter("unc_source==109").Filter("unc_scale==+1").Define("num_EleFakingTauES_DM1_Up", "weight").Sum('num_EleFakingTauES_DM1_Up')
        w_after_EleFakingTauES_DM1_Up= d.Filter("unc_source==109").Filter("unc_scale==+1").Define( "den_EleFakingTauES_DM1_Up", "weight*weight_btag_IterativeFit").Sum('den_EleFakingTauES_DM1_Up')
        r_factor_EleFakingTauES_DM1_Up=w_before_EleFakingTauES_DM1_Up.GetValue()/w_after_EleFakingTauES_DM1_Up.GetValue()
        add_dic(r_factor_EleFakingTauES_DM1_Up, "EleFakingTauES_DM1_Up",r_factors)

        w_before_EleFakingTauES_DM1_Down = d.Filter("unc_source==109").Filter("unc_scale==-1").Define("num_EleFakingTauES_DM1_Down", "weight").Sum('num_EleFakingTauES_DM1_Down')
        w_after_EleFakingTauES_DM1_Down = d.Filter("unc_source==109").Filter("unc_scale==-1").Define( "den_EleFakingTauES_DM1_Down", "weight*weight_btag_IterativeFit").Sum('den_EleFakingTauES_DM1_Down')
        r_factor_EleFakingTauES_DM1_Down=w_before_EleFakingTauES_DM1_Down.GetValue()/w_after_EleFakingTauES_DM1_Down.GetValue()
        add_dic(r_factor_EleFakingTauES_DM1_Down, "EleFakingTauES_DM1_Down",r_factors)

        w_before_MuFAkingTauES_Up= d.Filter("unc_source==108").Filter("unc_scale==+1").Define("num_MuFAkingTauES_Up", "weight").Sum('num_MuFAkingTauES_Up')
        w_after_MuFAkingTauES_Up= d.Filter("unc_source==108").Filter("unc_scale==+1").Define( "den_MuFAkingTauES_Up", "weight*weight_btag_IterativeFit").Sum('den_MuFAkingTauES_Up')
        r_factor_MuFAkingTauES_Up=w_before_MuFAkingTauES_Up.GetValue()/w_after_MuFAkingTauES_Up.GetValue()
        add_dic(r_factor_MuFAkingTauES_Up, "MuFAkingTauES_Up",r_factors)

        w_before_MuFAkingTauES_Down = d.Filter("unc_source==108").Filter("unc_scale==-1").Define("num_MuFAkingTauES_Down", "weight").Sum('num_MuFAkingTauES_Down')
        w_after_MuFAkingTauES_Down = d.Filter("unc_source==108").Filter("unc_scale==-1").Define( "den_MuFAkingTauES_Down", "weight*weight_btag_IterativeFit").Sum('den_MuFAkingTauES_Down')
        r_factor_MuFAkingTauES_Down=w_before_MuFAkingTauES_Down.GetValue()/w_after_MuFAkingTauES_Down.GetValue()
        add_dic(r_factor_MuFAkingTauES_Down, "MuFAkingTauES_Down",r_factors)


parser = argparse.ArgumentParser(description='Create bash command')
parser.add_argument('--ch', required=False, type=str, default= "all", help= "channel")
parser.add_argument('--year', required=False, type=str, default= "all", help= "year")
parser.add_argument('--unc_sources_group', required=False, type=str, default= "all", help="unc sources groups")
parser.add_argument('--input-dir', required=False, type=str, default="/mnt/data/Dottorato/anaTuples/2020-12-01/", help=" anatUples directory")
parser.add_argument('-n', required=False, type=bool, default=False, help=" don't write the file")

args = parser.parse_args()
channels=[]
if args.ch!="all":
    for ch in (args.ch).split(","):
      channels.append(ch)
else:
      channels.append("eTau")
      channels.append("muTau")
      channels.append("tauTau")
years=[]
if args.year!="all":
    for year in (args.year).split(","):
      years.append(year)
else:
      years.append("2016")
      years.append("2017")
      years.append("2018")

unc_sources=[]
if args.unc_sources_group!="all":
    for unc in (args.unc_sources_group).split(","):
      unc_sources.append(unc)
else:
      unc_sources.append("Central")
      unc_sources.append("LES")
      unc_sources.append("JES")

files = []
trees = []
for year in years:
    for channel in channels:
        r_factors={}
        for unc in unc_sources:
            evaluate_r(args.input_dir,channel,year,unc,r_factors)
        if(args.n==False):
            with open("btag_correction_factors_"+channel+"_"+year+".json", "w") as write_file:
                json.dump(r_factors, write_file)
