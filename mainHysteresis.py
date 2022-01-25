from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import numpy as np

def plot_params():
    plt.rc("axes", axisbelow=True)
    plt.tick_params(direction="in", length=5, colors="k", width=0.75)
    plt.grid(True, color="silver", linestyle="solid",
            linewidth=0.75, alpha=0.75)

def EnvTwoHyst(Dnum,Vnum,Dexp,Vexp):
    # peaks_pos_num, _ = find_peaks(Vnum, prominence=1, width=10)
    # peaks_neg_num, _ = find_peaks(-Vnum, prominence=1, width=10)
    # peaks_pos_exp, _ = find_peaks(Vexp, prominence=1, width=10)
    # peaks_neg_exp, _ = find_peaks(-Vexp, prominence=1, width=10)
    
    peaks_pos_num, _ = find_peaks(Vnum, prominence=1, width=15)
    peaks_neg_num, _ = find_peaks(-Vnum, prominence=1, width=15)
    peaks_pos_exp, _ = find_peaks(Vexp, prominence=1, width=15)
    peaks_neg_exp, _ = find_peaks(-Vexp, prominence=1, width=15)
    
    DnumEnv = np.concatenate([np.flip(Dnum[peaks_neg_num]), Dnum[peaks_pos_num]])
    VnumEnv = np.concatenate([np.flip(Vnum[peaks_neg_num]), Vnum[peaks_pos_num]])
    DexpEnv = np.concatenate([np.flip(Dexp[peaks_neg_exp]), Dexp[peaks_pos_exp]])
    VexpEnv = np.concatenate([np.flip(Vexp[peaks_neg_exp]), Vexp[peaks_pos_exp]])
    
    return DnumEnv, VnumEnv, DexpEnv, VexpEnv

def CompEnvHyst(Dnum,Vnum,Dexp,Vexp,Lunits,Funits):
    peaks_pos_num, _ = find_peaks(Vnum, prominence=1, width=100)
    peaks_neg_num, _ = find_peaks(-Vnum, prominence=1, width=100)
    peaks_pos_exp, _ = find_peaks(Vexp, prominence=1, width=100)
    peaks_neg_exp, _ = find_peaks(-Vexp, prominence=1, width=100)
       
    plt.figure(), plot_params()
    plt.plot(Vnum)
    plt.plot(peaks_pos_num, Vnum[peaks_pos_num],'xr', peaks_neg_num, Vnum[peaks_neg_num])
    plt.show()
    
    DnumEnv = np.concatenate([np.flip(Dnum[peaks_neg_num]), Dnum[peaks_pos_num]])
    VnumEnv = np.concatenate([np.flip(Vnum[peaks_neg_num]), Vnum[peaks_pos_num]])
    DexpEnv = np.concatenate([np.flip(Dexp[peaks_neg_exp]), Dexp[peaks_pos_exp]])
    VexpEnv = np.concatenate([np.flip(Vexp[peaks_neg_exp]), Vexp[peaks_pos_exp]])
    
    plt.figure(), plot_params()
    plt.plot(DexpEnv,VexpEnv,'-k',linewidth=1.3,label="Experimental")
    plt.plot(DnumEnv,VnumEnv,'--r',linewidth=1.3,label="Numérica")
    plt.xlabel("Desplazamiento "+Lunits)
    plt.ylabel("Fuerza "+Funits)
    plt.legend()
    