import openseespy.opensees as ops
import numpy as np
import time

# +-----------------------------------------------------+
# |          Función para análisis de gravedad          |
# +-----------------------------------------------------+
def run_gravity_analysis(steps=10):
    print("==> ... Running Gravity Analysis ...\n")
    ops.wipeAnalysis()
    ops.system("BandGeneral")
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0E-12, 10, 3)
    ops.algorithm("Newton")  # KrylovNewton
    ops.integrator("LoadControl", 1/steps)
    ops.analysis("Static")
    ops.analyze(steps)
    print("\n==> Gravity Analysis Completed!\n")
    # Set the gravity loads to be constant & reset the time in the domain
    ops.loadConst("-time", 0.0)

# +---------------------------------------------------------------------------+
# |       Función para análisis Pushover controlado por desplazamiento        |
# +---------------------------------------------------------------------------+
def run_pushover_analysis(ctrlNode, baseNodes, dof, Dincr, max_disp, IOflag=True):

    start_time = time.time()
   
    print("\n==> ... Running Pushover Analysis ...\n")

    testType = "NormDispIncr"  # EnergyIncr
    tolInit = 1.0e-6
    iterInit = 1000  # Set the initial Max Number of Iterations
    algorithmType = "KrylovNewton"  # Set the algorithm type

    ops.system("BandGeneral")
    ops.constraints('Penalty',1e20,1e20)   #Transformation
    ops.numberer("RCM")
    ops.test(testType, tolInit, iterInit)
    ops.algorithm(algorithmType)
    # Change the integration scheme to be displacement control
    #                                     node      dof  init Jd min max
    ops.integrator("DisplacementControl", ctrlNode, dof, Dincr)
    ops.analysis("Static")

    if IOflag:
        print(f"Single Pushover: Push node {ctrlNode} to {max_disp} metros.\n")

    # Set some parameters
    tCurrent = ops.getTime()
    currentStep = 1

    outputs = {
        "time": [],
        "rel_disp": [],
        "force": []
    }

    nodeList = []
    for node in baseNodes:
        nodeList.append(f"- ops.nodeReaction({node}, dof) ")

    nodeList = "".join(nodeList)
    currentDisp = ops.nodeDisp(ctrlNode, dof)
    ok = 0

    test = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    # algorithm = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}
    algorithm = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'NewtonLineSearch'}
    
    while ok == 0 and currentDisp < max_disp:
        ops.reactions()
        ok = ops.analyze(1)
        tCurrent = ops.getTime()
        currentDisp = ops.nodeDisp(ctrlNode, dof)

        if IOflag:
            print(
                f"Current displacement ==> {ops.nodeDisp(ctrlNode, dof):.3f} metros")

        if ok != 0:
            for i in test:
                for j in algorithm:
                    if j < 4:
                        ops.algorithm(algorithm[j], '-initial')
                            
                    else:
                        ops.algorithm(algorithm[j])
                    
                    print("\n==> Probando con:", test[i], algorithm[j]) 
                    ops.test(test[i], tolInit/0.1, iterInit)            
                    ok = ops.analyze(1) 
                    
                    if ok == 0:
                        print("==> That worked...\n")
                        break
                if ok==0:
                    break

        print('Step:',currentStep)
        currentStep += 1
        tCurrent = ops.getTime()

        outputs["time"].append(tCurrent)
        outputs["rel_disp"].append(ops.nodeDisp(ctrlNode, dof))
        outputs["force"].append(eval(nodeList))

    # Print a message to indicate if analysis completed succesfully or not
    if ok == 0:
        print("==> Pushover Analysis Completed SUCCESSFULLY!\n")
    else:
        print("==> Pushover Analysis Completed FAILED!\n")

    print(
        f"Pushover elapsed time is {(time.time() - start_time):.3f} seconds.\n")

    return outputs

# +--------------------------------------------------------------+
# |      Función para generar protocolo de desplazamiento        |
# +--------------------------------------------------------------+
def GeneratePeaks(Dmax,Dincr,CycleType,Fact):
    Disp = 0
    
    Dmax = Dmax*Fact
    dSteps = []
    
    if Dmax<0:
        dx = -Dincr
    else:
        dx = Dincr
    
    NstepsPeak = int(np.abs(Dmax)/Dincr)
    dSteps.append(Disp)
    
    for i in range(0, NstepsPeak):
        Disp = Disp+dx
        dSteps.append(Disp)
        
    if CycleType != "Push":
        for i in range(0, NstepsPeak):
            Disp = Disp-dx
            dSteps.append(Disp)
            
        if CycleType != "Half":
            for i in range(0, NstepsPeak):
                Disp = Disp-dx
                dSteps.append(Disp)
                
            for i in range(0, NstepsPeak):
                Disp = Disp+dx
                dSteps.append(Disp)
    
    np.savetxt('Dsteps.txt',dSteps)
    return dSteps

# +---------------------------------------------------------------------------+
# |      Función para análisis cíclico cuaciestático por desplazamiento       |
# +---------------------------------------------------------------------------+
def run_cyclic_analysis(ctrlNode, baseNodes, dof, iDmax, DincrStatic, CycleType, Fact, Ncycles):

    start_time = time.time()
    
    print("\n==> ... Running Cyclic Analysis ...\n")

    testType = "NormDispIncr"  # EnergyIncr
    tolInit = 1.0e-6
    iterInit = 1000  # Set the initial Max Number of Iterations
    algorithmType = "KrylovNewton"  # Set the algorithm type

    # ops.system('BandGeneral')
    ops.system('FullGeneral')
    ops.constraints('Penalty',1e20,1e20)
    # ops.constraints('Transformation')
    ops.numberer("RCM")
    ops.test(testType, tolInit, iterInit)
    ops.algorithm(algorithmType)

    # Set some parameters
    tCurrent = ops.getTime()
    currentStep = 1

    outputs = {
        "time": [],
        "rel_disp": [],
        "force": []
    }

    nodeList = []
    for node in baseNodes:
        nodeList.append(f"- ops.nodeReaction({node}, dof) ")

    nodeList = "".join(nodeList)
    ok = 0

    test = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    # algorithm = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}
    algorithm = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'NewtonLineSearch'}
    
    for DincrID,Dmax in  enumerate(iDmax):
        iDstep = GeneratePeaks(Dmax,DincrStatic[DincrID],CycleType,Fact)
        
        for i in range(Ncycles):
            D0 = 0.0
            for Dstep in iDstep:
                D1 = Dstep
                Dincr = D1-D0
                ops.reactions()
                ops.integrator("DisplacementControl", ctrlNode, dof, Dincr)
                ops.analysis("Static")
                ok = ops.analyze(1)
                tCurrent = ops.getTime()

                if ok != 0:
                    for i in test:
                        for j in algorithm:
                            if j < 4:
                                ops.algorithm(algorithm[j], '-initial')
                                    
                            else:
                                ops.algorithm(algorithm[j])
                            
                            print("\n==> Probando con:", test[i], algorithm[j])  
                            ops.test(test[i], tolInit/0.1, iterInit)                    
                            ok = ops.analyze(1) 
                            
                            if ok == 0:
                                print("==> That worked...\n")
                                break
                        if ok==0:
                            break

                print('Step:',currentStep)
                currentStep += 1
                tCurrent = ops.getTime()

                outputs["time"].append(tCurrent)
                outputs["rel_disp"].append(ops.nodeDisp(ctrlNode, dof))
                outputs["force"].append(eval(nodeList))
                
                D0=D1
                
                if ok != 0:
                    break
                
            if ok != 0:
                break
            
        if ok != 0:
            break   

    # Print a message to indicate if analysis completed succesfully or not
    if ok == 0:
        print("\n==> Cyclic Analysis Completed SUCCESSFULLY!\n")
    else:
        print("\n==> Cyclic Analysis Completed FAILED!\n")

    # print(
    #     f"Cyclic analysis elapsed time is {(time.time() - start_time):.3f} seconds.\n")

    return outputs

# +---------------------------------------------------------------------------+
# |      Función para análisis cíclico cuaciestático por desplazamiento       |
# +---------------------------------------------------------------------------+
def run_cyclic_analysis_txt(dispSteps, ctrlNode, baseNodes, dof, Fact):

    start_time = time.time()
    
    print("\n==> ... Running Cyclic Analysis ...\n")

    testType = "NormDispIncr"  # EnergyIncr
    tolInit = 1.0e-6
    iterInit = 1000  # Set the initial Max Number of Iterations
    algorithmType = "KrylovNewton"  # Set the algorithm type

    ops.system('BandGeneral')
    ops.constraints('Penalty',1e20,1e20)
    # ops.constraints('Transformation')
    ops.numberer("RCM")
    ops.test(testType, tolInit, iterInit)
    ops.algorithm(algorithmType)

    # Set some parameters
    tCurrent = ops.getTime()
    currentStep = 1

    outputs = {
        "time": [],
        "rel_disp": [],
        "force": []
    }

    nodeList = []
    for node in baseNodes:
        nodeList.append(f"- ops.nodeReaction({node}, dof) ")

    nodeList = "".join(nodeList)
    ok = 0
    D0=0.0
    
    test = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algorithm = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}
    # algorithm = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'NewtonLineSearch'}
    
    for Dstep in dispSteps:
        
        D1 = Dstep*Fact
        Dincr = D1-D0
        ops.integrator("DisplacementControl",  ctrlNode, dof, Dincr)
        ops.analysis("Static")
        ok=ops.analyze(1)
                        
        if ok!=0:
            for i in test:
                for j in algorithm:
                    if j < 4:
                        ops.algorithm(algorithm[j],'-initial')
                    else:
                        ops.algorithm(algorithm[j])
                        
                    print('\n\n==> Probando con:',test[i],algorithm[j])
                    ops.test(test[i],tolInit/0.1, iterInit,0)
                    ok = ops.analyze(1)
                    
                    if ok==0: 
                        print('==> That worked...\n')
                        break
                if ok==0: 
                    break
                                                
        # print('Step:',currentStep)
        currentStep += 1
        tCurrent = ops.getTime()

        outputs["time"].append(tCurrent)
        outputs["rel_disp"].append(ops.nodeDisp(ctrlNode, dof))
        outputs["force"].append(eval(nodeList))
                
        if ok!=0: 
            break
            
        D0=D1 # move to next step

    if ok == 0:
        print("==> Cyclic Analysis Completed SUCCESSFULLY!\n")
    else:
        print("==> Cyclic Analysis Completed FAILED!\n")
    
    print(f"Cyclic analysis elapsed time is {(time.time() - start_time):.3f} seconds.\n")

    return outputs