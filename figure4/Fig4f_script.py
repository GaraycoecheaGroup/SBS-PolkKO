# script to generate Fig 4f from fig_4f_start_data.tsv
# python version used: 3.9.18

from matplotlib import pyplot as plt

import numpy as np  # v1.26.4
import pandas as pd   # v1.5.3
import statsmodels.api as sm   # v0.14.4

"""
functions to retrieve the result statistics from the fit model
"""
def get_param_result(fit,param,CI_alpha=0.05):
    res = {}
    res['estimate'] = fit.params[param]
    res['pval'] = fit.pvalues[param]
    ci = tuple(fit.conf_int(alpha=CI_alpha).loc[param,:].values)
    res['lower'] = ci[0]
    res['upper'] = ci[1]
    return res

def get_model_results(fit_model,CI_alpha=0.05):
    results = {}
    terms = fit_model.params.keys()
    for term in terms:
        results[term] = get_param_result(fit_model,term,CI_alpha=CI_alpha)
    return results

"""
functions to arrange data for plotting as a stacked barplot
"""

def get_x_spread(n_groups,groupsize,within=1,extra_between=0.5,start=-0.5):
    res = []
    cur = start
    for g in range(n_groups):
        cur+=extra_between
        for n in range(groupsize):
            res.append(cur)
            cur+=within
        
    return res

def pad_missing(values, groupsize,position):
    """
    given y/height for 1 bar category, pad them with 0s
    depending on the number of categories and
    this categories position in each group
    """
    padded = []
    for value in values:
        cur_group = [0]*groupsize
        cur_group[position] = value
        padded += cur_group
    return padded

def get_yerr(est,lower,upper):
    conf_int = np.array([lower,upper]).T
    yerr = np.c_[est-conf_int[:,0],conf_int[:,1]-est ].T
    return yerr

def get_label_locs(x_values, n_groups,groupsize):
    """
    given x values, and of a number of grouped bars:
    get x-values at the middle of the groups
    """
    result = []
    cur_index=0
    
    for _ in range(n_groups):
        res = np.mean(x_values[cur_index:cur_index+groupsize])
        cur_index += groupsize
        result.append(res)
    return result

def get_cat_locs(x_values,groupsize,position):
    """
    given x_values of number of grouped bars. get x pos of one of the group categories
    e.g. the x value of the first number of each group
    """
    val_locs = np.arange(position,len(x_values)+position,groupsize)

    vals = np.array(x_values)[val_locs]
    return vals

if __name__ == "__main__":

    # load input table
    sbs_all_data_fn = './fig_4f_start_data.tsv'
    sbs_all_data = pd.read_csv(sbs_all_data_fn,sep='\t')


    """
    estimate the effect of each genotype with binomial regression
    """

    # outcome / dependent variable: mutation count in total coverage
    endog = sbs_all_data[['count','coverage']]
    
    # predictors / dependent variables for the full model including interaction terms
    exog_full = sbs_all_data[['background','PolK','Xpa','Xpc','Xpa_PolK','Xpc_PolK']]

    # fit binomial regression using the dependent and independent variables
    sbs_all_m2_res = sm.GLM(endog,exog_full,
    family=sm.families.Binomial(link=sm.families.links.Identity())).fit()
    m2_result = pd.DataFrame(get_model_results(sbs_all_m2_res)).T

    m2_result = m2_result.iloc[[0,1,3,2,5,4]]
    m2_result

    # convert the model outcome values from mutation rate per base to mutation rate per Mbase
    m2_perMb_result = m2_result.copy()
    m2_perMb_result[['estimate','lower','upper']] = m2_perMb_result[['estimate','lower','upper']]*1000000

    """
    plot the observed data and model output together in a stacked bar plot
    """

    # get the x axis locations and labels
    x_values = get_x_spread(6,2,within=0.75,extra_between=0.25)
    label_xs = get_label_locs(x_values,6,2)
    labels = ['WT','PolK -/-', 'Xpa -/-', 'Xpc -/-', 'PolK-Xpa -/-','PolK-Xpc -/-']
    labels = ['WT','Xpa -/-','Xpc -/-', 'Polk-/-', 'Polk-Xpa -/-', 'Polk-Xpc -/-']

    # get the model estimates, structure them for the barplot
    vwt,vwtl,vwtu = m2_perMb_result.loc['background',['estimate','lower','upper']]
    vpolk,vpolkl,vpolku = m2_perMb_result.loc['PolK',['estimate','lower','upper']]
    vxpa,vxpal,vxpau = m2_perMb_result.loc['Xpa',['estimate','lower','upper']]
    vxpc,vxpcl,vxpcu = m2_perMb_result.loc['Xpc',['estimate','lower','upper']]
    vxpap,vxpapl,vxpapu = m2_perMb_result.loc['Xpa_PolK',['estimate','lower','upper']]
    vxpcp,vxpcpl,vxpcpu = m2_perMb_result.loc['Xpc_PolK',['estimate','lower','upper']]

    wt = pad_missing([vwt]*6,2,1)
    wt_lower = pad_missing([vwtl]*6,2,1)
    wt_upper = pad_missing([vwtu]*6,2,1)
    wt_yerr = get_yerr(wt, wt_lower,wt_upper)

    polk = pad_missing([0,0,0,vpolk,vpolk,vpolk],2,1)
    polk_lower = pad_missing([0,0,0,vpolkl,vpolkl,vpolkl],2,1)
    polk_upper = pad_missing([0,0,0,vpolku,vpolku,vpolku],2,1)
    polk_yerr = get_yerr(polk, polk_lower,polk_upper)

    xpa = pad_missing([0,vxpa,0,0,vxpa,0],2,1)
    xpa_lower = pad_missing([0,vxpal,0,0,vxpal,0],2,1)
    xpa_upper = pad_missing([0,vxpau,0,0,vxpau,0],2,1)
    xpa_yerr = get_yerr(xpa, xpa_lower,xpa_upper)

    xpc = pad_missing([0,0,vxpc,0,0,vxpc],2,1)
    xpc_lower = pad_missing([0,0,vxpcl,0,0,vxpcl],2,1)
    xpc_upper = pad_missing([0,0,vxpcu,0,0,vxpcu],2,1)
    xpc_yerr = get_yerr(xpc, xpc_lower,xpc_upper)

    xpa_polk = pad_missing([0,0,0,0,vxpap,0],2,1)
    xpa_polk_lower = pad_missing([0,0,0,0,vxpapl,0],2,1)
    xpa_polk_upper = pad_missing([0,0,0,0,vxpapu,0],2,1)
    xpa_polk_yerr = get_yerr(xpa_polk, xpa_polk_lower,xpa_polk_upper)

    xpc_polk = pad_missing([0,0,0,0,0,vxpcp],2,1)
    xpc_polk_lower = pad_missing([0,0,0,0,0,vxpcpl],2,1)
    xpc_polk_upper = pad_missing([0,0,0,0,0,vxpcpu],2,1)
    xpc_polk_yerr = get_yerr(xpc_polk, xpc_polk_lower,xpc_polk_upper)

    wt_polk = np.add(wt,polk)
    wt_polk_xpa = np.add(wt_polk,xpa)
    wt_polk_xpc = np.add(wt_polk,xpc)

    # get the observed data and per-genotype means, structure them for barplot
    observed_means = sbs_all_data.groupby('genotype').mean(numeric_only=True)
    observed_means = observed_means.iloc[[1,2,4,0,3,5]]
    obs_means = pad_missing(observed_means['SBS per Mb'],2,0)

    obs_x_values = get_cat_locs(x_values,2,0)
    obs_data = []
    for i,gt in enumerate(['WT','Xpa','Xpc','Polk','Xpa-Polk','Xpc-Polk']):
        gt_counts = sbs_all_data[sbs_all_data['genotype']==gt]['SBS per Mb'].values
        x_vals = [obs_x_values[i]]*len(gt_counts)
        obs_data.append((x_vals,gt_counts))


    # plot the data in a barplot
    fig,ax = plt.subplots(figsize=(6,6))
    ax.bar(x_values,obs_means,label='observed',color='#3471eb',
        linewidth=0.9,edgecolor='black',width=0.6)
    ax.bar(x_values,wt,label='background',yerr=wt_yerr,color='grey',
        linewidth=0.9,edgecolor='black',width=0.6)
    ax.bar(x_values,polk,label='PolK',bottom=wt,yerr=polk_yerr,color='#D36135',
        linewidth=0.9,edgecolor='black',width=0.6)
    ax.bar(x_values,xpa,label='Xpa',bottom=wt_polk,yerr=xpa_yerr,color='#7FB069',
        linewidth=0.9,edgecolor='black',width=0.6)
    ax.bar(x_values,xpc,label='Xpc',bottom=wt_polk,yerr=xpc_yerr,color='#ECE4B7',
        linewidth=0.9,edgecolor='black',width=0.6)

    ax.bar(x_values,xpa_polk,label='Xpa-PolK',bottom=wt_polk_xpa,yerr=xpa_polk_yerr,color='#a434eb',
        linewidth=0.9,edgecolor='black',width=0.6)
    ax.bar(x_values,xpc_polk,label='Xpc-PolK',bottom=wt_polk_xpc,yerr=xpc_polk_yerr,color='#eb3446',
        linewidth=0.9,edgecolor='black',width=0.6)
    for x,y in obs_data:
        ax.scatter(x,y,color='black',s=8)

    ax.legend()
    ax.set_xlabel('genotype')
    ax.set_ylabel('SBS per Mb')

    plt.xticks(label_xs,rotation=90)
    ax.set_xticklabels(labels)
    plt.tight_layout()
    plt.show()