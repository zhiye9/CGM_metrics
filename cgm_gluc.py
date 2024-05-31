import numpy as np

def eA1c(x, mmol_L = False):
    """Estimated A1C values with avarage glucose values by using the formula (mean(glucose) + 46.7)/28.7 with mg/dL and (mean(glucose) + 2.59)/1.59 with mmol/L.
        Parameters:
            x: CGM values
            mmol_L: bool, default = False
                If False, the x should be in mg/dL. Otherwise, it should be in mmol/L
        Returns:
            mg/dL value
        ----------------
        References:
           Nathan, David M., et al. "Translating the A1C assay into estimated average glucose values." Diabetes care 31.8 (2008): 1473-1478. 
    """
    if mmol_L:
        return (np.nanmean(x) + 2.59)/1.59
    else:
        return (np.nanmean(x) + 46.7)/28.7

def TR_min(x, sd = 1., sr=5., TIR = False, range = [70, 180]):
    """Compute time of CGM values in and out of the target range in minutes.
        Parameters:
            x: CGM values
            sd: float, default = 1
                Standard deviation value for the range
            sr: float, default = 5
                Sampling rate
            TIR: bool, default = False
                If False, the range will be calculated by using the mean and standard deviation. Otherwise, the range will be given by the user.
            range: list, default = [70, 180]
                The standard target range of CGM 
        Returns:
            TIR: Time in Range (TIR)
            TBR: Time Below Range (TBR)
            TAR: Time Above Range (TAR)
            TOR: Time out of Range (TOR)
        ----------------
        References:
           Gabbay, Monica Andrade Lima, et al. "Time in range: a new parameter to evaluate blood glucose control in patients with diabetes." Diabetology & metabolic syndrome 12 (2020): 1-8.
    """
    if TIR:
        up = range[1]
        dw = range[0]
    else:
        up = np.nanmean(x) + sd*np.nanstd(x)
        dw = np.nanmean(x) - sd*np.nanstd(x)
    TIR = len(x[(x<= up) & (x>= dw)])
    TBR = len(x[x < dw])
    TAR = len(x[x > up])
    return TIR, TBR, TAR, TAR + TBR

def TR_percent(x, sd = 1, sr=5, TIR = False, range = [70, 180]):
    """Compute time of CGM values in and out of the target range in percentage.
    """
    if TIR:
        up = range[1]
        dw = range[0]
    else:
        up = np.nanmean(x) + sd*np.nanstd(x)
        dw = np.nanmean(x) - sd*np.nanstd(x)
    TIR = len(x[(x<= up) & (x>= dw)]) / len(x)
    TBR = len(x[x < dw]) / len(x)
    TAR = len(x[x > up]) / len(x)
    return TIR, TBR, TAR, TAR + TBR

def ADRR(df, mmol_L = False):
    """Compute the average daily risk range (ADRR) by normalizing the glucose values according to [1][2][3] firstly, then compute the ADDR according to [4].  
        Parameters:
            df: Data frame with columns 'glucose' 'Time' and 'Day'.
            mmol_L: bool, default = False
                If False, the x should be in mg/dL. Otherwise, it should be in mmol/L

        Returns:
            ADRR: Average Daily Risk Range
            LBGI: Low Blood Glucose (Risk) Index
            HBGI: High Blood Glucose (Risk) Index
        ----------------
        References:
           [1] Kovatchev, Boris P., et al. "Symmetrization of the blood glucose measurement scale and its applications." Diabetes care 20.11 (1997): 1655-1658.
           [2] Kovatchev, Boris P., et al. "Assessment of risk for severe hypoglycemia among adults with IDDM: validation of the low blood glucose index." Diabetes care 21.11 (1998): 1870-1875.
           [3] Kovatchev, Boris P., et al. "Risk analysis of blood glucose data: Az quantitative approach to optimizing the control of insulin dependent diabetes." Computational and Mathematical Methods in Medicine 3.1 (2000): 1-10.
           [4] Kovatchev, Boris P., et al. "Evaluation of a new measure of blood glucose variability in diabetes." Diabetes care 29.11 (2006): 2433-2438.
    """
    if mmol_L:
        df['glucose'] = df['glucose'].apply(lambda x: 1.794*(np.power(np.log(x), 1.026) - 1.861))
    else:
        df['glucose'] = df['glucose'].apply(lambda x: 1.509*(np.power(np.log(x), 1.084) - 5.381))

    # compute rh and rl
    df['rh'] = 10*np.square(df['glucose'].apply(lambda x: max(x, 0)))
    df['rl'] = 10*np.square(df['glucose'].apply(lambda x: min(x, 0)))

    HGRI = np.nanmean(df['rh'])
    LBGI = np.nanmean(df['rl'])

    # compute the max rh and rl for each day and ADRR
    ADRR = np.nanmean(df.groupby('Day')['rh'].max() + df.groupby('Day')['rl'].max())

    return ADRR, HGRI, LBGI

def COGI(x, range = [[70, 180], [3.9, 10]], overall_weights = [0.5, 0.35, 0.15], mmol_L = False):
    """Compute the continous glucose monitoring index (COGI) for the given CGM values.
    Parameters:
        x: CGM values
        range: list, default = [[70, 180], [3.9, 10]]
            The standard target range of CGM
        overall_weights: list, default = [0.5, 0.35, 0.15]
            Weights for TIR, TBR and SD
        mmol_L: bool, default = False
            If False, the x should be in mg/dL. Otherwise, it should be in mmol/L
    Returns:
        COGI: Continous Glucose Monitoring Index
    ----------------
    References:
       Leelarathna, Lalantha, et al. "Evaluating glucose control with a novel composite continuous glucose monitoring index." Journal of diabetes science and technology 14.2 (2020): 277-283.
    """
    SD = np.nanstd(x)
    if mmol_L:
        TIR, TBR, _, _ = TR_percent(x, TIR = True, range = range[1])
        return TIR*100*overall_weights[0] + (1 - np.clip(TBR / 0.15, 0, 1))*100*overall_weights[1] + (100 - np.clip((SD - 1) / (6 - 1) * 100, 0, 100))
    else:
        TIR, TBR, _, _ = TR_percent(x, TIR = True, range = range[0])
        return TIR*100*overall_weights[0] + (1 - np.clip(TBR / 0.15, 0, 1))*100*overall_weights[1] + (100 - np.clip((SD - 18) / (108 - 18) * 100, 0, 100))*overall_weights[2]

def CONGA(x, hours = 1):
    """Compute the continuous overall net glycemic action (CONGA) for the given CGM values.
    Parameters:
        x: CGM values
    Returns:
        CONGA: Continuous Overall Net Glycemic Action
    ----------------
    References:
       McDonnell, C. M., et al. "A novel approach to continuous glucose analysis utilizing glycemic variation." Diabetes technology & therapeutics 7.2 (2005): 253-263.
    """


    return np.nanstd(x.diff())

def MODD(x):
    """Compute the mean of daily differences (MODD) for the given CGM values.
    Parameters:
        x: CGM values
    Returns:
        MODD: Mean of Daily Differences
    ----------------
    References:
       Molnar, G. D., W. F. Taylor, and M. M. Ho. "Day-to-day variation of continuously monitored glycaemia: a further measure of diabetic instability." Diabetologia 8.5 (1972): 342-348.
    """
    return np.nanmean(np.abs(x.diff()))

def GMI(x):
    return 12.71 + (4.70587*np.nanmean(x))

def Jindex(x):
    return 0.324*((np.nanmean(x)+np.nanstd(x))**2)
    
def SDRC(x, t=5):
    return np.nanstd(np.abs(x.diff()/t))
