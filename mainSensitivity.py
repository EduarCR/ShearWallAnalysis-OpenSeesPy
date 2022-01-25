# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 18:34:01 2022

@author: EduarCR
"""

from SFI_MVLEM2D_M3R10 import modelo as sfimvlem2D
from SALib.plotting.bar import plot as barplot
from SALib.sample import saltelli
from SALib.analyze import sobol
import matplotlib.pyplot as plt
import mainHysteresis as hysEdu
import hysteresis as hys
import seaborn as sns
import pandas as pd
import numpy as np
import itertools
import time
import os

Folder = 'Gráficas'
os.makedirs(Folder, exist_ok=True)

# Parámetros de figuras
def plot_params():
    plt.rc("axes", axisbelow=True)
    plt.rcParams.update({'font.size': 8})
    plt.rc('font', family='serif')
    plt.tick_params(direction="in", length=5, colors="k", width=0.75)
    plt.grid(True, color="silver", linestyle="solid",
            linewidth=0.75, alpha=0.75)


sns.set_context("paper", font_scale=0.9, rc={"lines.linewidth": 2.5})

# Función de desempeño
# --------------------
def funcDes(lw,R0,cR1,cR2,xcrn,xcrp,ru,rt,nu,alfa):
    
    numHys, expHys, Nfases, Dnum_Vnum, Dexp_Vexp = sfimvlem2D(lw,R0,cR1,cR2,xcrn,xcrp,ru,rt,nu,alfa)
    
    # Chequear sensibilidad
    
    lpSteps = [2]*int(Nfases) # 2 ciclos por Nfases
    Diff, Diffs = hys.compareHys(numHys, expHys)  # Función que calcula las diferencias entre curvas
    
    # Gráficas
    # --------
    # Envolventes numéricas y experimentales
    DnumEnv, VnumEnv, DexpEnv, VexpEnv = hysEdu.EnvTwoHyst(Dnum=Dnum_Vnum[:,0], Vnum=Dnum_Vnum[:,1],
                                                           Dexp=Dexp_Vexp[:,0],Vexp=Dexp_Vexp[:,1])
    
    plt.figure(figsize=(4, 3.5)), plot_params()
    plt.plot(DexpEnv,VexpEnv,'-k',linewidth=1.3,label="Experimental")
    plt.plot(DnumEnv,VnumEnv,'--r',linewidth=1.3,label="Numérica")
    plt.xlabel("Desplazamiento [mm]"), plt.ylabel("[Fuerza [kN]")
    plt.legend(loc="lower right")
    plt.savefig("Gráficas/EnvAvg "+str(round(Diff,2))+" "+str(lw)+"lw"+".png",dpi=1200, bbox_inches='tight')
    plt.show()

    # Hystéresis
    plt.figure(figsize=(4, 3.5)), plot_params()
    expHys.plot(color='black', linewidth = 0.9, linestyle='--', label = 'Experimental')
    numHys.plot(color="red", linewidth=0.9, linestyle="--", label = 'Numérica')
    plt.xlabel("Desplazamiento [mm]"), plt.ylabel("[Fuerza [kN]")
    plt.legend()
    plt.savefig("Gráficas/CompHys "+str(round(Diff,2))+" "+str(lw)+"lw"+".png",dpi=1200, bbox_inches='tight')
    plt.show()

    # Diferencia entre hystéresis
    plt.figure(figsize=(4, 3.5)), plot_params()
    plt.stem(np.arange(len(Diffs)), Diffs, 'black')
    plt.xlabel('Ciclo (#)')
    plt.ylabel('Diferencia promedio')
    plt.savefig("Gráficas/compHysDiff "+str(round(Diff,2))+" "+str(lw)+"lw"+".png",dpi=1200, bbox_inches='tight')
    plt.show()
    
    return Diff, np.array([DnumEnv,VnumEnv]), len(DnumEnv)

start_time = time.time()

# Análisis de sensibilidad
# ------------------------

# Definimos la cantidad de iteraciones y el problema
Ns = 3
problem = {
    "num_vars": 9, 
    "names": ["R0","cR1","cR2","xcrn","xcrp","ru","rt","nu","alfa"], 
    "bounds": [[10, 20],[0.925,0.925*0.02],[0.002, 0.15],
               [1.035,1.035*0.07],[10000,10000*0.07],[5,50], [1.2,1.2*0.07],
               [0.1, 1.5],[0.01, 0.05]],
    'dists': ["unif","norm",'unif',"norm",'norm',"unif","norm","unif","unif"]
}

param_values = saltelli.sample(problem, 2**Ns)

# Evaluar el modelo + Salida gráfica
# ----------------------------------
hystExp = np.loadtxt('Hyst. M3R10.txt') #Cargar datos experimentales
Vexp = hystExp[:,1]
Dexp = hystExp[:,0]

lws = [1.5]    # Longitud del muro a evaluar

Y = np.zeros([param_values.shape[0],len(lws)])
envNum = np.zeros((50,2,param_values.shape[0]))
pos = np.zeros(param_values.shape[0])

for j,lw in enumerate(lws):
    
    for i, X in enumerate(param_values):
        Y[i,j], EnvNum, lennum= funcDes(lw,X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7],X[8])
        envNum[0:lennum,:,i] = np.transpose(EnvNum)
        pos1 = np.where(envNum[:,0,i]==0)
        pos[i] = np.array(pos1)[0,0]
    
    
    plt.figure(figsize=(3.5, 2.7)), plot_params()
    DnumEnv, VnumEnv, DexpEnv, VexpEnv = hysEdu.EnvTwoHyst(Dnum=Dexp, Vnum=Vexp,
                                                           Dexp=Dexp,Vexp=Vexp)
      
    for k in range(param_values.shape[0]):
        if k!=159 and k!=155 and k!=126:
            plt.plot(envNum[0:int(pos[k]),0,k],envNum[0:int(pos[k]),1,k],'.-k',linewidth=0.7)
        # plt.scatter(envNum[0:int(pos[k]),0,k],envNum[0:int(pos[k]),1,k])
    
    plt.plot(DexpEnv,VexpEnv,'-r',linewidth=1.3,label="Experimental")
    plt.xlabel("Desplazamiento [mm]"), plt.ylabel("[Fuerza [kN]")
    plt.savefig("Envolventes_1.png",dpi=1200, bbox_inches='tight')
    plt.show()    
    
    sensitivity = sobol.analyze(problem, Y[:,j], print_to_console=True)
    total, first, second = sensitivity.to_df()
    
    dictt=dict(itertools.islice(sensitivity.items(), 4))
    result_sens =pd.DataFrame(dictt,index=["R0","cR1","cR2","xcrn","xcrp","ru","rt","nu","alfa"])
    
    rows, columns = 1, 2
    grid = plt.GridSpec(rows, columns, wspace=0.3, hspace=0.25)
    plt.figure(figsize=(6, 2.7))
    plt.subplot(grid[0]), plot_params()
    plt.bar(x=["R0","cR1","cR2","xcrn","xcrp","ru","rt","nu","alfa"], height=sensitivity["ST"], color="black", label="ST")
    plt.xticks(rotation=45)
    plt.legend()
    plt.subplot(grid[1]), plot_params()
    plt.bar(x=["R0","cR1","cR2","xcrn","xcrp","ru","rt","nu","alfa"], height=sensitivity["S1"], color="black", label="S1")
    plt.xticks(rotation=45)
    plt.legend()
    plt.savefig("SobolMod "+" "+str(lw)+"lw"+".png",dpi=1200, bbox_inches='tight')
    plt.show()   
    
    fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(14,5))
    ax1 = barplot(total, ax=ax1)
    ax2 = barplot(first, ax=ax2)
    ax3 = barplot(second, ax=ax3)
    
    plt.savefig("Sobol "+" "+str(lw)+"lw"+".png",dpi=1200, bbox_inches='tight')
    np.savetxt("Sensitivity "+ str(lw)+"lw"+".txt",np.transpose([sensitivity["ST"], sensitivity["S1"]]))
    plt.show()
       
    plt.figure(figsize=(3.5,2.7)), plot_params()
    plt.step(np.arange(1,len(Y)+1,1),Y,'black')
    plt.xlabel("Iteraciones"), plt.ylabel("Diferencia total promedio")
    plt.savefig("Gráficas/DifTotal "+" "+str(lw)+"lw"+".png",dpi=1200, bbox_inches='tight')
    np.savetxt("Diff_total "+ str(lw)+"lw"+".txt",Y[:,j])
    plt.show()

print(
        f"\nTotal elapsed time is {(time.time() - start_time):.3f} seconds.\n")
