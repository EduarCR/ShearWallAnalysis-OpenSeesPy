# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 18:14:55 2022

@author: EduarCR
"""

# -------------------------------------------------------------------
# Modelo 1  :   Simulación de comportamiento cíclico usando SFI_MVLEM
# Espécimen :   M3R10 (Ortega et al., 2020)
# Fecha     :   05/01/2022
# Update    :   24/01/2022

# +-----------------------------+
# |         Librerias           |
# +-----------------------------+
import openseespy.postprocessing.Get_Rendering as opsplt
from scipy.signal import find_peaks
import openseespy.opensees as ops
import mainAnalysis as man
import hysteresis as hys
import numpy as np
import os

# +-----------------------------+
# |          Unidades           |
# +-----------------------------+
# Unidades base
# -------------
m, mN, sec = 1.0, 1.0, 1.0

# Otras unidades
# --------------
mm  = m*1E3
N   = mN/1E6
kN  = mN/1E3
Pa, kPa, MPa, GPa = N/m**2, 1E3*N/m**2, 1E6*N/m**2, 1E9*N/m**2
g   = 9.80665*m/(sec**2)

# Conversión de unidades
# ----------------------
inch  = 0.0254*m
ft    = 12*inch
kip   = N/0.000225
ksi   = kip/(inch**2)
psi   = ksi/1E3

def modelo(lw,R0,cR1,cR2,xcrn,xcrp,ru,rt,nu,alfa):
    # +-------------------------------------------------------------+
    # |   Inicio de generación de modelo (Units: mN, m, sec, KPa)   |
    # +-------------------------------------------------------------+
    NombreModelo = 'M3R10-SFI_MVLEM'
    os.makedirs(NombreModelo, exist_ok=True)

    # Crear modelo 2D
    ops.wipe()
    ops.model('basic','-ndm', 2, '-ndf',3)

    # +-----------------------------------------+
    # |     Geometría, nodos, restricciones     |
    # +-----------------------------------------+
    # Geometría del Muro
    # ------------------
    hw = 2.4*m    # altura
    tw = 0.1*m    # espesor
    lw = lw*m     # longitud
    lb = 0.2*m    # longitud de elmento de borde

    # Nodos
    # -----
    # node nodeId xCrd yCrd zCrd
    ops.node(1,0,0)
    ops.node(2,0,0.48)
    ops.node(3,0,0.96)
    ops.node(4,0,1.44)
    ops.node(5,0,1.92)
    ops.node(6,0,hw)
    
    # Nodos adicionales (protocolo de desplazamiento aplicado 
    # a 1.5 m de la corona del muro)
    ops.node(7,0,hw+1.5/2)  
    ops.node(8,0,hw+1.5)

    # restricciones
    ops.fix(1,1,1,1) # Nodo 1 restingido

    # Nodo y gdl de control
    IDctrNode  = 6
    IDctrDOF   = 1

    # +-----------------------------------------+
    # |               MATERIALES                |
    # +-----------------------------------------+
    # ACERO
    # -----
    # acero en X (transversal)
    fyX = 590*MPa       # fy
    bx  = 0.0029        # strain hardening

    # acero en Y alma (longitudinal)
    fyYw = 590*MPa      # fy
    byw  = 0.0029       # strain hardening

    # acero en Y EB (longitudinal)
    fyYb = 460*MPa      # fy
    byb  = 0.0023       # strain hardening

    # otros parámetros
    Esy = 200000*MPa    # Young's modulus
    Esx = Esy           # Young's modulus
    R0 = R0             # Valor inicial de parámetro de curvatura
    cR1 = cR1           # Parámetro de degradación de curvatura
    cR2 = cR2           # Parámetro de degradación de curvatura 

    # Definición de materiales de ACERO
    # uniaxialMaterial('SteelMPF', matTag, fyp, fyn, E0, bp, bn, *params, a1=0.0, a2=1.0, a3=0.0, a4=1.0)
    ops.uniaxialMaterial('SteelMPF',1,fyX,fyX,Esx,bx,bx,R0,cR1,cR2)         # acero en X
    ops.uniaxialMaterial('SteelMPF',2,fyYw,fyYw,Esy,byw,byw,R0,cR1,cR2)     # acero en Y web
    ops.uniaxialMaterial('SteelMPF',3,fyYb,fyYb,Esy,byb,byb,R0,cR1,cR2)     # acero en Y EB


    # CONCRETO
    # --------
    # Inconfinado
    fpc = 24.5*MPa      # peak compressive stress
    ec0 = -0.002        # strain at peak compressive stress
    ft = 1.5*MPa        # peak tensile stress
    et = 0.0004         # strain at peak tensile stress
    Ec = 19304*MPa      # Young's modulus   
    xcrnu = xcrn        # cracking strain - compression
    xcrp = xcrp         # cracking strain - tension
    ru = ru             # shape parameter - compression
    rt = rt             # shape parameter - tension

    # Definiendo materiales de CONCRETO
    # uniaxialMaterial('ConcreteCM', matTag, fpcc, epcc, Ec, rc, xcrn, ft, et, rt, xcrp, '-GapClose', GapClose=0)
    ops.uniaxialMaterial('ConcreteCM',4,-fpc,ec0,Ec,ru,xcrnu,ft,et,rt,xcrp,'-GapClose',0)       # concreto inconfinado

    # +-----------------------------------------+
    # |        Parámetros material FSAM         |
    # +-----------------------------------------+
    # Cuantías de refuerzo
    rouXw = 0.0027989   # X web
    rouXb = rouXw       # X EB
    rouYw = rouXw       # Y web
    rouYb = 0.026127    # Y EB

    # Parámetros de mecanismo de resstencia al cortante
    nu = nu            # Concrete friction coefficient (0.0<ν<1.5)
    alfadow = alfa     # Stiffness coefficient of reinforcement dowel action (0.0<alfadow<0.05)

    # nDMaterial('FSAM', matTag, rho, sXTag, sYTag, concTag, rouX, rouY, nu, alfadow)
    ops.nDMaterial('FSAM', 6,    0.0, 1,     2,     4,       rouXw,rouYw,nu, alfadow) # Alma (concreto inconfinado)
    ops.nDMaterial('FSAM', 7,    0.0, 1,     3,     4,       rouXb,rouYb,nu, alfadow) # EB (concreto inconfinado)

    # +-----------------------------------------+
    # |     Definiendo elementos SFI_MVLEM      |
    # +-----------------------------------------+
    if lw==0.6: m1=3
    elif lw==0.8: m1=4
    elif lw==1.0: m1=5
    elif lw==1.2: m1=6
    elif lw==1.4: m1=7
    elif lw==1.5: m1=6
    elif lw==1.6: m1=8
    elif lw==1.8: m1=9
    elif lw==2.0: m1=10
    elif lw==2.2: m1=11
    elif lw==2.4: m1=12
    elif lw==2.6: m1=13
    elif lw==2.8: m1=14
    elif lw==3.0: m1=15
    
    thick = np.ones(m1)*tw
    width = np.concatenate([[lb],np.ones(m1-2)*(lw-2*lb)/(m1-2),[lb]])
    mat = np.concatenate([[7],np.ones(m1-2)*6,[7]])


    # element('SFI_MVLEM', eleTag, *eleNodes, m, c, '-thick', *thick, '-width', *widths, '-mat', *mat_tags)
    ops.element('SFI_MVLEM', 1,  1, 2, m1, 0.4, '-thick', *thick, '-width', *width, '-mat', *mat)
    ops.element('SFI_MVLEM', 2,  2, 3, m1, 0.4, '-thick', *thick, '-width', *width, '-mat', *mat)
    ops.element('SFI_MVLEM', 3,  3, 4, m1, 0.4, '-thick', *thick, '-width', *width, '-mat', *mat)
    ops.element('SFI_MVLEM', 4,  4, 5, m1, 0.4, '-thick', *thick, '-width', *width, '-mat', *mat)
    ops.element('SFI_MVLEM', 5,  5, 6, m1, 0.4, '-thick', *thick, '-width', *width, '-mat', *mat)

    # Agregar elemento superior rígido
    # --------------------------------
    ops.geomTransf('Linear', 1)
    A=0.2*0.2
    E=200000*MPa*1e10
    Iz=0.2*0.2**3/12

    ops.element('ModElasticBeam2d' , 11 , *[6,7], A, E, Iz , 1, 1e3, 1, 1)
    ops.element('ModElasticBeam2d' , 12 , *[7,8], A, E, Iz , 1, 1e3, 1, 1)

    # opsplt.plot_model('nodes','elements')

    # +-----------------------------------------+
    # |     Análisis por cargas de gravedad     |
    # +-----------------------------------------+

    P = 0.184*mN  # Carga Axial
    ops.timeSeries("Linear", 1)
    ops.pattern('Plain',1,1)
    ops.load(IDctrNode,0.0,-P,0.0)	# vertical load node 11

    man.run_gravity_analysis(steps=10)

    # Calcular reaccion vertical
    # --------------------------
    ops.reactions()
    baseNodes = [1]
    print("Reaccion vertical",ops.nodeReaction(baseNodes[0],2)/kN,' kN')


    # +-----------------------------------------+
    # |      Recorders - salida del modelo      |
    # +-----------------------------------------+
    # Nodal recorders
    ops.recorder('Node', '-file', NombreModelo+"/Dtop.out", '-time', '-node', IDctrNode, '-dof',IDctrDOF, 'disp')
    ops.recorder('Node', '-file', NombreModelo+"/React.out", '-time', '-node', 1, '-dof',1, 'reaction')


    # ------------------------------------------+
    # |           Análisis cíclico              |
    # ------------------------------------------+
    ops.wipeAnalysis()
    ops.loadConst("-time", 0.0)

    # Aplicar carga lateral con serie de tiempo lineal
    # ------------------------------------------------
    ops.timeSeries('Linear',2)
    ops.pattern('Plain',2,2)

    Plateral = 1.0           # Carga lateral de referencia
    ops.load(8,Plateral,0.0,0.0)

    # Cargar protocolo de desplazamiento
    # ----------------------------------
    stop = 2686
    protDsteps = np.loadtxt('PROT_M3R10_2.txt')
    protDsteps = protDsteps[0:stop]
    baseNodes = [1]
    Fact=1

    cyclic_output = man.run_cyclic_analysis_txt(dispSteps=protDsteps, ctrlNode=IDctrNode,
                            baseNodes=baseNodes, dof=IDctrDOF, Fact=Fact)


    # +----------------------------------------+
    # |              Resultados                |
    # -----------------------------------------+

    # Resultados analíticos
    # ---------------------
    Dnum=np.array(cyclic_output["rel_disp"])*1e3
    Vnum=np.array(cyclic_output["force"])*1e3
    Dnum_Vnum = np.column_stack([Dnum,Vnum])
    numHys = hys.Hysteresis(Dnum_Vnum)
    peaks_pos_num, _ = find_peaks(Dnum)
    
    # Condiciones para poder realizar comparaciones entre curvas si el análisis falla
    Npeaks =    len(peaks_pos_num)
    if Npeaks%2==0 and len(Dnum)>2682: # 12 fases
        Nfases=Npeaks/2
    elif Npeaks%2==0 and len(Dnum)>2304: # 11 fases
        Nfases=Npeaks/2
    elif Npeaks%2==0 and len(Dnum)>1920: # 10 fases
        Nfases=Npeaks/2
    elif Npeaks%2==0 and len(Dnum)>1535: # 9 fases
        Nfases=Npeaks/2
    elif Npeaks%2==0 and len(Dnum)>1152: # 8 fases
        Nfases=Npeaks/2
    elif Npeaks%2==0 and len(Dnum)>832: # 7 fases
        Nfases=Npeaks/2
    elif Npeaks%2==0 and len(Dnum)>575: # 6 fases
          Nfases=Npeaks/2
    elif Npeaks%2==0 and len(Dnum)>382: # 5 fases
          Nfases=Npeaks/2
    elif Npeaks%2==0 and len(Dnum)>252: # 4 fases
          Nfases=Npeaks/2
    else:
        Nfases=(Npeaks-1)/2
    
    # Cargando curva experimental
    hystExp = np.loadtxt('Hyst. M3R10.txt')
    Vexp = hystExp[:,1]
    Dexp = hystExp[:,0]
    Dexp_Vexp = np.column_stack([Dexp,Vexp])
    
    tam1=np.size(Dnum_Vnum,0)
    tam2=np.size(Dexp_Vexp,0)
    tam=min(tam1,tam2)
    Dexp_Vexp=Dexp_Vexp[0:tam,:]
    expHys = hys.Hysteresis(Dexp_Vexp)

    
    return numHys,expHys,Nfases,Dnum_Vnum,Dexp_Vexp
