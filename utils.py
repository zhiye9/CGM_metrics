import numpy as np

def mg_dL_to_mmol_L(x):
    """Convert mg/dL to mmol/L.
        Parameters:
            x:mg/dL value
        Returns:
            mmol/L value
        ----------------
        References:
           Riemsma, Rob, et al. "Integrated sensor-augmented pump therapy systems [the MiniMed® Paradigm™ Veo system and the Vibe™ and G4® PLATINUM CGM (continuous glucose monitoring) system] for managing blood glucose levels in type 1 diabetes: a systematic review and economic evaluation." Health Technology Assessment (Winchester, England) 20.17 (2016): 1. 
    """
    return x*0.0555

def mmol_L_to_mg_dL(x):
    return x*18.0188

def separate_day_time(df):
    """Separate the time into two columns by space.
        Parameters:
            x: CGM data frame
        Returns:
            Date and time
    """
    df['Day'] = df['time'].dt.date
    df['Time'] = df['time'].dt.time
    return df




