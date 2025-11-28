
# print(rsn_timeseries.shape)
from statsmodels.graphics.tsaplots import plot_acf
from statsmodels.tsa import stattools
from scipy.interpolate import CubicSpline
import numpy as np
import sys
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import nolds

def exp_decay(lag, A, tau):
    return A * np.exp(-lag / tau)

def intrinsic_timescale_expfit(signal, max_lag = 1000):
    # 1. 计算 ACF（不含 lag=0）
    autocorr = stattools.acf(signal, nlags = max_lag, fft=True)
    lags = np.arange(1, max_lag + 1)
    y = autocorr[1:]  # 排除 lag=0

    # 2. 拟合指数模型 A * exp(-t / tau)
    try:
        popt, _ = curve_fit(exp_decay, lags, y, p0=(1.0, 10.0), bounds=(0, np.inf))
        A_fit, tau_fit = popt
        return tau_fit
    except RuntimeError:
        # 拟合失败时返回 NaN
        return np.nan

def intrinic_timescale00(ica):
    autocorr = stattools.acf(ica,nlags = 1000000)
    # newx = np.arange(0,200,0.01)
    # spl = np.interp(newx,np.arange(0,200,1),autocorr[:200])
    return np.sum(autocorr[1:list(autocorr<0.0001).index(True)])

def intrinic_timescale01(ica):
    autocorr = stattools.acf(ica,nlags = 1000000)
    return np.sum(autocorr[1:list(autocorr<0.1).index(True)])

def intrinic_timescale05(ica):
    autocorr = stattools.acf(ica,nlags = 1000000)
    return np.sum(autocorr[1:list(autocorr<0.5).index(True)])

def acfstrength00(ica):
    autocorr = stattools.acf(ica,nlags = 1000000)
    newx = np.arange(0,len(autocorr),0.01)
    spl = np.interp(newx, np.arange(0,len(autocorr),1),autocorr[:len(autocorr)])
    return newx[np.where(spl < 0 )[0][0]]

def acfstrength01(ica):
    autocorr = stattools.acf(ica,nlags = 1000000)
    newx = np.arange(0,len(autocorr),0.01)
    spl = np.interp(newx, np.arange(0,len(autocorr),1),autocorr[:len(autocorr)])
    return newx[np.where(spl < 0.1 )[0][0]]

def acfstrength05(ica):
    autocorr = stattools.acf(ica,nlags = 1000000)
    newx = np.arange(0,len(autocorr),0.01)
    spl = np.interp(newx, np.arange(0,len(autocorr),1),autocorr[:len(autocorr)])
    return newx[np.where(spl < 0.5 )[0][0]]

def ac1(ica):
    autocorr = stattools.acf(ica,nlags = 1000000)
    return autocorr[1]

def hurst(ts, long = True):
    if long:
        return nolds.hurst_rs(ts)
    
    ts = list(ts)
    N = len(ts)
    if N < 20:
        raise ValueError("Time series is too short! input series ought to have at least 20 samples!")

    max_k = int(np.floor(N/2))
    R_S_dict = []
    for k in range(10,max_k+1):
        R,S = 0,0
        # split ts into subsets
        subset_list = [ts[i:i+k] for i in range(0,N,k)]
        if np.mod(N,k)>0:
            subset_list.pop()
            #tail = subset_list.pop()
            #subset_list[-1].extend(tail)
        # calc mean of every subset
        mean_list=[np.mean(x) for x in subset_list]
        for i in range(len(subset_list)):
            cumsum_list = pd.Series(subset_list[i]-mean_list[i]).cumsum()
            R += max(cumsum_list)-min(cumsum_list)
            S += np.std(subset_list[i])
        R_S_dict.append({"R":R/len(subset_list),"S":S/len(subset_list),"n":k})
    
    log_R_S = []
    log_n = []
    # print(R_S_dict)
    for i in range(len(R_S_dict)):
        R_S = (R_S_dict[i]["R"]+np.spacing(1)) / (R_S_dict[i]["S"]+np.spacing(1))
        log_R_S.append(np.log(R_S))
        log_n.append(np.log(R_S_dict[i]["n"]))

    Hurst_exponent = np.polyfit(log_n,log_R_S,1)[0]
    return Hurst_exponent

def CR_RAD(x, tau=1, do_abs=True):
    """ 
    Compute the rescaled auto-density (RAD), a metric for inferring the distance to criticality.

    Parameters
    -----------------
    x : array_like
        Input 1D time series.
    tau : int, optional
        Embedding and differencing delay in time steps (default is 1).
    do_abs : bool, optional
        Whether to center the time series at 0 and take the absolute value (default is True).
    Returns
    -------
    f : float
        The computed RAD feature value.
    """
    x = np.asarray(x,dtype=np.float64).flatten()  # Ensure 1D array
    if do_abs:
        x = x - np.median(x)
        x = np.abs(x)
    if tau >= len(x):
        raise ValueError("tau is too large for the length of the time series")
   # Delay embedding (m = 2)
    x1 = x[:-tau]
    y = x[tau:]
    # Median split
    median_val = np.median(x1)
    sub_medians = x1 < median_val
    super_median_sd = np.std(x1[~sub_medians], ddof=1)
    sub_median_sd = np.std(x1[sub_medians], ddof=1)
    # Properties of the auto-density
    sigma_dx = np.std(y - x1, ddof=1)
    density_difference = (1.0 / super_median_sd) - (1.0 / sub_median_sd)
    f = sigma_dx * density_difference
    return f


