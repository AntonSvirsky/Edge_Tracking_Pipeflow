import netCDF4 as nc
import numpy as np
import os
import sys
import subprocess
import configparser
import matplotlib.pyplot as plt
import opf_analysis   as oa
import bisec_algo_aux as aux
import time
import logging

global main_dir, utils_dir, state_file, config_file
global IC_dir, Bisec_dir, FwdInt_dir
global sbatch
######################################################################################
def state_IC():
    logger.info(f'IC {IC_dir}')
    s_IC_1() # Initiate bisection folder   
    return
#------------------------------------------------------------------------------------    
def s_IC_1():
    ## Initiate bisection folder
    #convert IC states to mat_vec.cdf
    Lz, nx, ny ,nz = aux.load_p2m_settings(config_file) #load importing parameters
    for idx in [1,2]: #create input file for extraction script
        input_str  = ["1","1","0","0","1","0",Lz,nx,ny,nz, "1", "1",
                      f"./s_{idx}p.cdf.dat"]
        with open(os.path.join(IC_dir,f"input{idx}.inp"), "w") as f:
            f.writelines(s + '\n' for s in input_str)
    subprocess.run([f"cp {utils_dir}/mv_ic_ext.sh {IC_dir}/mv_ic_ext.sh"], shell=True) #copy the script to dir
    subprocess.run([os.path.join(IC_dir,"mv_ic_ext.sh")], shell=True)                  #RUN the script
    subprocess.run(f'rm {os.path.join(IC_dir,"mv_ic_ext.sh")}', shell=True)            #remove the script

    d = oa.distanceL2_mv_files('mv_1p.cdf','mv_2p.cdf',IC_dir,sub_HPf = False) #compute distance
    
    #create a bisection dir with the inital status fule
    subprocess.run([f"mkdir -p {Bisec_dir}"], shell=True)
    subprocess.run([f"cp {IC_dir}/s_1p.cdf.dat {Bisec_dir}/s_1p_0.cdf.dat"], shell=True)
    subprocess.run([f"cp {IC_dir}/s_2p.cdf.dat {Bisec_dir}/s_2p_0.cdf.dat"], shell=True)

    #create bisection status files
    with open(os.path.join(Bisec_dir,f"status"), "w") as f:
        f.writelines("s_1p_0.cdf.dat"+ "\n") #acrive 1p state
        f.writelines("s_2p_0.cdf.dat"+ "\n") #active 2p state
        f.writelines(f"{d}"+ "\n")           #distance between states
        f.writelines("0"+ "\n")              #number of bisection
        f.writelines("none"+ "\n")           #b_state
        #f.writelines("1p_old"+ "\n")         #indicates wheter 1p was replaced
        #f.writelines("2p_old"+ "\n")         #indicates wheter 2p was replaced

    #update main state
    aux.updateMainState('S','Bis',state_file) 
    return
######################################################################################
def state_Bis():
    logger.info(f'Bis {Bisec_dir}')              #OUT
    b_state = aux.get_Bi_state(Bisec_dir) 
    logger.info(f'b_state {b_state}')            #OUT
    
    match b_state:
        case 'none':     s_Bis_1() # Check bisection status
        case 'pre-run':  s_Bis_2() # Start bisection Run
        case 'running':  s_Bis_3() # Check if simulation running
        case 'stoped':   s_Bis_4() # Classification of finished run
        case 'finished': s_Bis_5() # Initiate fwdInt dir
        case 'unknown':  s_Bis_6() # Continue integration of unknown state

    return
#------------------------------------------------------------------------------------
def s_Bis_1():
    ## Check bisection status. (Is another bisection required or can continue to FwdInt?).
    #load current state
    bs_1p,bs_2p,dist,b_M = aux.get_bisection_status(Bisec_dir)
    #load targer distance bor bisection L_bs
    m_config_File = configparser.ConfigParser()
    m_config_File.read(config_file)
    L_bs = np.double(m_config_File['State distance']['L_bs']) #bisection distance treshold
    force_change = m_config_File['State distance']['force_change'].lower() == "true" #if required to change representative state 
    
    #check distance and that 
    if dist > L_bs:
        logger.info(f'd([{bs_1p}],[{bs_2p}]) = {dist:.3e} > {L_bs:.3e} --> continue bisection')
        aux.update_Bi_state("pre-run", Bisec_dir)
        #write distance in log
        with open(os.path.join(Bisec_dir,f"bisection.log"), "a") as f:
            f.writelines(f"{b_M:04} : d([{bs_1p}] ,[{bs_2p}]) = {dist}>{L_bs}"+ "\n")
    elif force_change and (bs_1p == 's_1p_0.cdf.dat' or bs_2p == 's_2p_0.cdf.dat'):
        logger.info(f'did not find replacemnt for both states {bs_1p} {bs_2p} --> continue bisection')
        aux.update_Bi_state("pre-run", Bisec_dir)
        #write distance in log
        with open(os.path.join(Bisec_dir,f"bisection.log"), "a") as f:
            f.writelines(f"{b_M:04} : d([{bs_1p}] ,[{bs_2p}]) = {dist}<{L_bs}"+ "\n")
            f.writelines(f"{b_M:04} : {bs_1p} or {bs_2p} did not change"+ "\n")
    else:
        logger.info(f'd([{bs_1p}],[{bs_2p}]) = {dist:.3e} <= {L_bs:.3e} --> bisection finished')
        #write distance in log
        with open(os.path.join(Bisec_dir,f"bisection.log"), "a") as f:
            f.writelines(f"{b_M:04} : d([{bs_1p}] ,[{bs_2p}]) = {dist}<{L_bs} ---> bisection finished"+ "\n")
        aux.update_Bi_state("finished", Bisec_dir)
    return 
#------------------------------------------------------------------------------------
def s_Bis_2():
    ## Start Bisection Run
    #load current state
    bs_1p,bs_2p,dist,b_M = aux.get_bisection_status(Bisec_dir)
    # Create avg state between bs_1p and bs_2p using the 'avg_state.sh' script
    with open(os.path.join(Bisec_dir,f"input0.inp"), "w") as f:
        f.writelines(bs_1p + '\n') 
        f.writelines(bs_2p + '\n')
    subprocess.run(f"cp {utils_dir}/avg_state.sh {Bisec_dir}/avg_state.sh", shell=True) #copy the script to dir
    subprocess.run(f"{Bisec_dir}/avg_state.sh", shell=True) #RUN the script
    subprocess.run(f'rm {Bisec_dir}/avg_state.sh', shell=True) #remove the script    
    
    # Create a simulation for the avg_run
    subprocess.run(f"cp -r {utils_dir}/bisec_sim {Bisec_dir}/bisec_sim_{b_M:04}", shell=True) #copy the script to dir
    subprocess.run(f"cp {Bisec_dir}/s_avg.cdf.dat {Bisec_dir}/bisec_sim_{b_M:04}/state.cdf.in", shell=True) #copy the script to dir
    
    #Run bisection sim
    subprocess.run(f"cp {utils_dir}/sim_run.sh {Bisec_dir}/bisec_sim_{b_M:04}/sim_run.sh", shell=True) #copy the script to run sim
    subprocess.run(f"{Bisec_dir}/bisec_sim_{b_M:04}/sim_run.sh", shell=True) #RUN sim
    
    aux.update_Bi_state("running", Bisec_dir)
    return 
#------------------------------------------------------------------------------------
def s_Bis_3():
    ## Check if simulation running
    bs_1p,bs_2p,dist,b_M = aux.get_bisection_status(Bisec_dir)
    isRunning, save_rate,max_step,s_files_count = aux.get_sim_status(f"{Bisec_dir}/bisec_sim_{b_M:04}") 
    max_s_files = max_step/save_rate
    if isRunning:
        logger.info(f'Sim running with: save_rate = {save_rate}, max_step = {max_step}. Currenty: {s_files_count} / {max_s_files} state files')    
        #sleep for c_timer_min
        c_timer = aux.get_cTimer(config_file)
        logger.info(f'Sleeping for {c_timer} min')
        time.sleep(60*c_timer)
        logger.info(f'Waking up') 
    else:
        logger.info(f'Sim sroped with: save_rate = {save_rate}, max_step = {max_step}. Currenty: {s_files_count} / {max_s_files} state files')    
        aux.update_Bi_state("stoped", Bisec_dir)
        return 
#------------------------------------------------------------------------------------
def s_Bis_4():
    ## Classification of finished run
    bs_1p,bs_2p,dist,b_M = aux.get_bisection_status(Bisec_dir)
    #Classify the simulation run and update bisection state accordingly
    s_files_count = oa.count_cdf_files(f"{Bisec_dir}/bisec_sim_{b_M:04}/", end_with = 'cdf.dat')
    # Create input file for the extraction script
    Lz, nx, ny ,nz = aux.load_p2m_settings(config_file) #load importing parameters
    input_str  = ["1","0","0","0","1","0",Lz,nx,ny,nz,"1","2","./",
                  "0",str(s_files_count-1),"1"]
    #extract mat_vec from stateXXXX.cdf.dat files
    oa.convertToMatVec(f"bisec_sim_{b_M:04}" ,input_str, Bisec_dir, utils_dir)
    mv_files_count = oa.count_cdf_files(f"{Bisec_dir}/bisec_sim_{b_M:04}/mat_vec_out/", end_with = 'cdf')
    logger.info(f'Extraction complete: {mv_files_count} mat_vec.cdf files') #OUT
    
    #Classify the simulation run and update bisection state accordingly
    state_def = aux.get_state_def(config_file)
    cls,t = oa.classificationPlot(f"bisec_sim_{b_M:04}",state_def,Bisec_dir, Bisec_dir,Bisec_dir)
    
    #acknolage classification and update state
    
    #b_cls = cls[-1] if np.var(cls[-3:]) == 0 else 0  #clasify based on the last 3 snaps. if not equal classify as 0. 
    #b_cls = cls[-1] #classify based on last snapshot
    b_cls = aux.get_classification(cls,t,config_file)
    
    if b_cls == 0: #case that no classification was obtained
        logger.info(f'Unable to find classification for run') #OUT 
        with open(os.path.join(Bisec_dir,f"bisection.log"), "a") as f:
            f.writelines(f"{b_M:04} : [{bs_1p}] ,[{bs_2p}] ---> s_avg.cdf.dat ---> unknown ----> continue"+ "\n")
        aux.update_Bi_state("unknown", Bisec_dir)
        return
    
    bs_new_name = f's_{b_cls}p_{b_M+1:04}.cdf.dat'
    logger.info(f'Clasified as {bs_new_name}') #OUT
    subprocess.run(f"mv {Bisec_dir}/s_avg.cdf.dat {Bisec_dir}/{bs_new_name}", shell=True) #rename the classified state

    #update logger
    with open(os.path.join(Bisec_dir,f"bisection.log"), "a") as f:
        f.writelines(f"{b_M:04} : [{bs_1p}] ,[{bs_2p}] ---> s_avg.cdf.dat ---> [{bs_new_name}]"+ "\n")
    
    #replace the approporiate active 1p or 2p  for the bisected run
    if b_cls == 1:   
        bs_1p = bs_new_name
    elif b_cls == 2: 
        bs_2p = bs_new_name
    else:            logger.info('Error! no valid classification obtained!') #OUT ERROR
    #extract mat_vec representaion of new pair and compute distance
    selected_states = [bs_1p,bs_2p]
    for idx in [0,1]:
        input_str  = ["1","1","0","0","1", "0",Lz,nx,ny,nz,"1", "1", 
                      selected_states[idx]]
        with open(os.path.join(Bisec_dir,f"input{idx+1}.inp"), "w") as f:
            f.writelines(s + '\n' for s in input_str)
    
    subprocess.run([f"cp {utils_dir}/mv_ic_ext.sh {Bisec_dir}/mv_ic_ext.sh"], shell=True) #copy the script to dir
    subprocess.run([os.path.join(Bisec_dir,"mv_ic_ext.sh")], shell=True)                  #RUN the script
    subprocess.run(f'rm {os.path.join(Bisec_dir,"mv_ic_ext.sh")}', shell=True)            #remove the script
    
    d = oa.distanceL2_mv_files('mv_1p.cdf','mv_2p.cdf',Bisec_dir,sub_HPf = False) #compute distance
    logger.info(f'Dist = {d}') #OUT
    
    #write distance in log
    with open(os.path.join(Bisec_dir,f"bisection.log"), "a") as f:
        f.writelines(f"{b_M:04} : d([{bs_1p}] ,[{bs_2p}]) = {d}"+ "\n")
    
    #update status file "none"
    str_write = [bs_1p,bs_2p,f"{d}",f"{b_M+1:04}","none"]
    with open(os.path.join(Bisec_dir,f"status"), "w") as f:
        f.writelines(s + "\n" for s in str_write)
    return 
#------------------------------------------------------------------------------------
def s_Bis_5():
    ## Initiate fwdInt folder
    #create a fwdInt dir with two simulations, with IC as the active bisection states indicated by the status file of bisection_N
    subprocess.run([f"mkdir -p {FwdInt_dir}"], shell=True)
    fwdInt_1p_ic,fwdInt_2p_ic,dist,_ = aux.get_bisection_status(Bisec_dir)

    #Create simulation folders
    subprocess.run(f"cp -r {utils_dir}/fwdInt_sim {FwdInt_dir}/fwdInt_1p", shell=True)
    subprocess.run(f"cp -r {utils_dir}/fwdInt_sim {FwdInt_dir}/fwdInt_2p", shell=True) 
    # copy IC from bisection folder
    subprocess.run([f"cp {Bisec_dir}/{fwdInt_1p_ic} {FwdInt_dir}/fwdInt_1p/state.cdf.in"], shell=True)
    subprocess.run([f"cp {Bisec_dir}/{fwdInt_2p_ic} {FwdInt_dir}/fwdInt_2p/state.cdf.in"], shell=True)
    
    # Create fwdInt status 
    with open(os.path.join(FwdInt_dir,f"status"), "w") as f:
        f.writelines(f"{Bisec_dir}"+ "\n")    #bisection dir
        f.writelines(f"{fwdInt_1p_ic}"+ "\n") #acrive 1p state
        f.writelines(f"{fwdInt_2p_ic}"+ "\n") #active 2p state
        f.writelines(f"{dist}"+ "\n")         #inital distance between states
        f.writelines("none"+ "\n")            #fwdInt_1_state
        f.writelines("none"+ "\n")            #fwdInt_2_state

    aux.updateMainState('S','Fwd',state_file) 
    return 
#-------------------------------------------------------------------------------------
def s_Bis_6():
    ## Continue integration for bisection with of unknown classification
    bs_1p,bs_2p,dist,b_M = aux.get_bisection_status(Bisec_dir)
    s_cnt = oa.count_cdf_files(f"{Bisec_dir}/bisec_sim_{b_M:04}", end_with = 'cdf.dat')
    
    # change inital condition
    subprocess.run(f"cp -r {utils_dir}/bisec_sim {Bisec_dir}/bisec_sim_{b_M+1:04}", shell=True) #create new sim dir
    subprocess.run(f"cp {Bisec_dir}/bisec_sim_{b_M:04}/state{s_cnt-1:04}.cdf.dat {Bisec_dir}/bisec_sim_{b_M+1:04}/state.cdf.in", shell=True) #set final state as new IC
    
    logger.info(f'{Bisec_dir}/bisec_sim_{b_M:04}/state{s_cnt-1:04}.cdf.dat --> {Bisec_dir}/bisec_sim_{b_M+1:04}/state.cdf.in"') 
    
    #Run bisection sim
    subprocess.run(f"cp {utils_dir}/sim_run.sh {Bisec_dir}/bisec_sim_{b_M+1:04}/sim_run.sh", shell=True) #copy the script to run sim
    subprocess.run(f"{Bisec_dir}/bisec_sim_{b_M+1:04}/sim_run.sh", shell=True) #RUN sim
    
    #update status file "running"
    str_write = [bs_1p,bs_2p,f"{dist}",f"{b_M+1:04}","running"]
    with open(os.path.join(Bisec_dir,f"status"), "w") as f:
        f.writelines(s + "\n" for s in str_write)
    return 
######################################################################################
def state_Fwd():
    logger.info(f'Fwd {FwdInt_dir}')
    flag_stop   = 0;
    flag_finish = 0;
    for fwd_run in [1,2]:
        fwd_state = aux.get_fwdInt_state(fwd_run,FwdInt_dir)
        match fwd_state:
            case "none"     : s_Fwd_1(fwd_run) # Run fwdInt sim
            case "running"  : s_Fwd_2(fwd_run) #Check if simulation running
            case "stoped"   : flag_stop   += 1
            case "finished" : flag_finish += 1   

    if   flag_stop == 2:   s_Fwd_3() # Get new pair of states 
    elif flag_finish == 2: s_Fwd_4() # Make IC_{N+1} folder
    return
#------------------------------------------------------------------------------------
def s_Fwd_1(idx):
    ## Run fwdInt sims
    logger.info(f'Starting fwdInt_{idx} run')
    subprocess.run(f"cp {utils_dir}/sim_run.sh {FwdInt_dir}/fwdInt_{idx}p/sim_run.sh", shell=True) #copy the script to run sim
    subprocess.run(f"{FwdInt_dir}/fwdInt_{idx}p/sim_run.sh", shell=True) #RUN fwd Sim
    aux.update_fwdInt_state(idx,"running", FwdInt_dir)
    return 
#------------------------------------------------------------------------------------
def s_Fwd_2(idx):
    ## Check if simulation running
    fwdIn_sim_dir = f"{FwdInt_dir}/fwdInt_{idx}p" 
    isRunning, save_rate,max_step,s_files_count = aux.get_sim_status(fwdIn_sim_dir) 
    max_s_files = max_step/save_rate
    if isRunning:
        logger.info(f'fwdInt sim {idx} running with: save_rate = {save_rate}, max_step = {max_step}. Currenty: {s_files_count} / {max_s_files} state files')    
        #sleep for c_timer_min
        c_timer = aux.get_cTimer(config_file)
        logger.info(f'Sleeping for {c_timer} min')
        time.sleep(60*c_timer)
        logger.info(f'Waking up') 
    else:
        logger.info(f'fwdInt sim {idx} sroped with: save_rate = {save_rate}, max_step = {max_step}. Currenty: {s_files_count} / {max_s_files} state files')
        aux.update_fwdInt_state(idx,"stoped", FwdInt_dir)
    return 
#------------------------------------------------------------------------------------
def s_Fwd_3():
    ## Get new pair of states 
    #extract mat_vec from stateXXXX.cdf.dat files for each
    for i in [1,2]:
        sf_count = oa.count_cdf_files(f"{FwdInt_dir}/fwdInt_{i}p/", end_with = 'cdf.dat')
        
        # Create input file for the extraction script
        Lz, nx, ny ,nz = aux.load_p2m_settings(config_file) #load importing parameters
        input_str  = ["1","0","0","0","1", "0",Lz, nx, ny ,nz, "1", "2","./",
                      "0",str(sf_count-1),"1" ]
        #copy and run extraction script
        oa.convertToMatVec(f"fwdInt_{i}p" ,input_str, FwdInt_dir, utils_dir)
        mv_files_count = oa.count_cdf_files(f"{FwdInt_dir}/fwdInt_{i}p/mat_vec_out/", end_with = 'cdf')
        logger.info(f'FwdInt {i} extraction complete: {mv_files_count} mat_vec.cdf files')
    #get classification plots
    state_def = aux.get_state_def(config_file)
    cls1,t2 = oa.classificationPlot("fwdInt_1p",state_def,FwdInt_dir, FwdInt_dir,FwdInt_dir)
    cls1,t2 = oa.classificationPlot("fwdInt_2p",state_def,FwdInt_dir, FwdInt_dir,FwdInt_dir)
    #find distance between states
    sf_1_count = oa.count_cdf_files(f"{FwdInt_dir}/fwdInt_1p/mat_vec_out", end_with = 'cdf')
    sf_2_count = oa.count_cdf_files(f"{FwdInt_dir}/fwdInt_2p/mat_vec_out", end_with = 'cdf')
    if sf_1_count!= sf_2_count: logger.info('Warning: diffrent simulation lengths')
        
    sf_N = np.min([sf_1_count,sf_2_count])
    dist_fwdInt = np.zeros(sf_N)
    for j in range(sf_N):
        f1 = f"fwdInt_1p/mat_vec_out/mat_vec{j:04}.cdf"
        f2 = f"fwdInt_2p/mat_vec_out/mat_vec{j:04}.cdf"
        dist_fwdInt[j] = oa.distanceL2_mv_files(f1,f2,FwdInt_dir,sub_HPf = True)
    #load targer distance bor bisection L_bs and L_fwd
    m_config_File = configparser.ConfigParser()
    m_config_File.read(config_file)
    L_bs  = np.double(m_config_File['State distance']['L_bs'])
    L_fwd = np.double(m_config_File['State distance']['L_fwd'])
    
    #find new states as the first to cross distance
    #larger_dist = np.where(dist_fwdInt > L_fwd)[0]
    larger_dist = np.where(dist_fwdInt < L_fwd)[0] #chose the last point bellow the treshold
    new_idx = larger_dist[-1] if np.any(larger_dist) else int(np.ceil(sf_count*0.2))
    logger.info(f"fwdInt states diverge at idx = {new_idx}, with d = {dist_fwdInt[new_idx]:.4e}") #OUT
    
    #Generate plot for distance and inital and final states
    oa.GenDistPlots(FwdInt_dir,dist_fwdInt,new_idx,state_def,L_bs,L_fwd,FwdInt_dir)
    oa.fwdInt_initalFinalStatesPlot(FwdInt_dir,dist_fwdInt,new_idx)
    
    #add the selected states to the status file
    with open(os.path.join(FwdInt_dir,f"status"),'r') as f:
        lines = f.readlines()
    with open(os.path.join(FwdInt_dir,f"status"), "w") as f:
        f.writelines(lines)
        f.writelines(f"fwdInt_1p/state{new_idx:04}.cdf.dat"+ "\n") 
        f.writelines(f"fwdInt_2p/state{new_idx:04}.cdf.dat"+ "\n")
        f.writelines(str(dist_fwdInt[new_idx])+ "\n")

    #update status file to indicate job finished
    aux.update_fwdInt_state(1,"finished", FwdInt_dir)
    aux.update_fwdInt_state(2,"finished", FwdInt_dir)
    return 
#------------------------------------------------------------------------------------
def s_Fwd_4():
    ## Make IC_{N+1} folder
    #make new IC folder
    New_IC_dir = os.path.join(main_dir,f'IC_{m_N+1 :04}')
    subprocess.run([f"mkdir -p {New_IC_dir}"], shell=True)
    
    with open(os.path.join(FwdInt_dir,f"status")) as f:
        lines = f.read().splitlines()
        New_1p_ic = lines[6]  #new 1 puff state
        New_2p_ic = lines[7]  #new 2 puff state

    #cope state files to IC folder
    subprocess.run(f"cp {FwdInt_dir}/{New_1p_ic} {New_IC_dir}/s_1p.cdf.dat", shell=True)
    subprocess.run(f"cp {FwdInt_dir}/{New_2p_ic} {New_IC_dir}/s_2p.cdf.dat", shell=True) 
    aux.updateMainState('S','Adv',state_file)
    return
######################################################################################
def state_Adv():
    s_Adv_1()
    return
#------------------------------------------------------------------------------------
def s_Adv_1():
    logger.info(f'Advancing to next iteration {m_N} --> {m_N+1}')
    aux.updateMainState('N',m_N+1,state_file)
    globals()['IC_dir']     = os.path.join(main_dir,f'IC_{m_N+1 :04}')
    globals()['Bisec_dir']  = os.path.join(main_dir,f'Bisec_{m_N+1 :04}')
    globals()['FwdInt_dir'] = os.path.join(main_dir,f'FwdInt_{m_N+1 :04}')
    aux.updateMainState('S','IC',state_file)
    return
######################################################################################
def getMainState():
    #initialize the script
    m_state_File = configparser.ConfigParser()
    m_state_File.read(state_file)
    m_state  = str(m_state_File['State']['S'])
    m_N      = int(m_state_File['State']['N'])
    return m_state,m_N 
######################################################################################

## Main: 

##parameters
#main_dir    = '../bisec_R2050_L180_01'
#main_dir    = '/storage/ph_anna/antonsv/openpipe/runs/bisec_R2050_L150_01'
#main_dir    = '/storage/ph_anna/antonsv/openpipe/runs/bisec_R2050_L150_02'
main_dir = sys.argv[1]
m_N_max     = 100

#init logger
logger = logging.getLogger(__name__)
date_strftime_format = '%Y-%m-%y %H:%M:%S'
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format="%(asctime)s %(message)s", datefmt=date_strftime_format)

#set global dirs
utils_dir   = os.path.join(main_dir,'utils')
state_file  = os.path.join(main_dir,'main.state')
config_file = os.path.join(main_dir,'bisec.config')

m_state,m_N = getMainState()

globals()['IC_dir']     = os.path.join(main_dir,f'IC_{m_N:04}')
globals()['Bisec_dir']  = os.path.join(main_dir,f'Bisec_{m_N:04}')
globals()['FwdInt_dir'] = os.path.join(main_dir,f'FwdInt_{m_N:04}')

while m_N<=m_N_max:
    logger.info(f'MAIN: N={m_N}, S={m_state}') 
    match m_state:
        case 'IC':  state_IC()
        case 'Bis': state_Bis()
        case 'Fwd': state_Fwd()
        case 'Adv': state_Adv()
    m_state,m_N = getMainState()
    time.sleep(1)
                
