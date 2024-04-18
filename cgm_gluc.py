import numpy as np

def mg_dL_to_mmol_L(x):
    return x*0.0555
def mmol_L_to_mg_dL(x):
    return x*18.0188

def eA1c(x, mmol_L = False):
    if mmol_L:
        return (np.nanmean(x['glucose']) + 2.59)/1.59
    else:
        return (np.nanmean(x['glucose']) + 46.7)/28.7

def TR_min(x, sd = 1, sr=5, range = []):
    up = np.nanmean(x['glucose']) + sd*np.nanstd(x['glucose'])
    dw = np.nanmean(x['glucose']) - sd*np.nanstd(x['glucose'])
    TIR = len(x[(x['glucose']<= up) & (x['glucose']>= dw)])/96
    TBR = len(x[x['glucose'] < dw])/96
    TAR = len(x[x['glucose'] > up])/96 
    return TIR, TBR, TAR

def SDRC(x, t=5):
    return np.nanstd(np.abs(x['glucose'].diff()/t))

def GMI(x):
    return 12.71 + (4.70587*np.nanmean(x['glucose']))

def Jindex(x):
    return 0.324*((np.nanmean(x['glucose'])+np.nanstd(x['glucose']))**2)
    