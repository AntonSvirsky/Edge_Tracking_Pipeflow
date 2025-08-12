import os
import configparser
import numpy as np
import opf_analysis   as oa
import time

def updateMainState(field,new_state,state_file):
    #update main state file
    m_state_File = configparser.ConfigParser()
    m_state_File.read(state_file)
    m_state_File['State'][field] = str(new_state)
    with open(state_file, 'w') as f:
      m_state_File.write(f)
    return
def load_p2m_settings(config_file):
    m_config_File = configparser.ConfigParser()
    m_config_File.read(config_file)
    cset = m_config_File['Conversion settings']
    Lz, nx, ny ,nz = cset['Lz'],cset['nx'],cset['ny'],cset['nz']
    return Lz, nx, ny ,nz

def get_sim_status(sim_dir):
    isRunning = os.path.isfile(f"{sim_dir}/RUNNING") or (not os.path.isfile(f"{sim_dir}/state0000.cdf.dat"))
    with open(f"{sim_dir}/main.info") as f:
            lines = f.read().splitlines()
    save_rate = int(lines[18].split("=")[1])
    max_step  = int(lines[20].split("=")[1])
    s_files_count = oa.count_cdf_files(sim_dir, end_with = 'cdf.dat')
    
    return isRunning, save_rate,max_step,s_files_count

def get_cTimer(config_file):
    m_config_File = configparser.ConfigParser(); 
    m_config_File.read(config_file)
    c_timer_min =  int(m_config_File['Misc']['c_timer_min'])
    return c_timer_min

def get_state_def(config_file):
    m_config_File = configparser.ConfigParser(); 
    m_config_File.read(config_file)
    s_def = m_config_File['State deffention']
    state_def = {
                "F_th"      : np.double(s_def["F_th"]),
                "min_gap"   : np.double(s_def["min_gap"]),
                "small_gap" : np.double(s_def["small_gap"]),
                "p0_Ft_lim" : np.fromstring(s_def["p0_Ft_lim"], dtype=float, count=2, sep = ','),
                "p1_Ft_lim" : np.fromstring(s_def["p1_Ft_lim"], dtype=float, count=2, sep = ','),
                "p2_Ft_lim" : np.fromstring(s_def["p2_Ft_lim"], dtype=float, count=2, sep = ','),
                "p3_Ft_lim" : np.fromstring(s_def["p3_Ft_lim"], dtype=float, count=2, sep = ',')
                }
    return state_def
#----------------------------------------------------------------------
def get_Bi_state(Bisec_dir):
    with open(os.path.join(Bisec_dir,f"status"),'r') as f:
        lines = f.read().splitlines()
    return lines[4]

def update_Bi_state(state, Bisec_dir):
    with open(os.path.join(Bisec_dir,f"status"),'r') as f:
        lines = f.readlines()
    lines[4] = state + "\n"
    with open(os.path.join(Bisec_dir,f"status"),'w') as f:
        f.writelines(lines)
    return
def get_bisection_status(Bisec_dir):
    with open(os.path.join(Bisec_dir,f"status")) as f:
        lines = f.read().splitlines()
    bs_1p = lines[0]             #Current 1 puff state
    bs_2p = lines[1]             #Current 2 puff state
    dist  = np.double(lines[2])  #Distance between bs_1p and bs_2p
    b_M   = int(lines[3])        #index of the current bisection
    return bs_1p,bs_2p,dist,b_M

def get_classification(cls,t,config_file):
    #get minimum time period to concider classification
    m_config_File = configparser.ConfigParser(); 
    m_config_File.read(config_file)
    min_class_time =  float(m_config_File['State deffention']['min_class_time'])
    currennt_class = cls[0]; #current considered classification
    time_in_class  = 0.0;    #time spent in currennt_class
    for i in range(1,len(cls)):
        if cls[i] == currennt_class: 
            time_in_class += t[i] - t[i-1]
        else: 
            currennt_class,time_in_class = cls[i], 0.0;

        #classification is accepted if time spent in classification is suficently long
        if (time_in_class >= min_class_time) and currennt_class != 0:
            print(f'Selected class {currennt_class} after time {time_in_class}')
            return currennt_class
    #if no valid classification was obtained, return 0 (no class)
    currennt_class = 0
    return currennt_class
    

#----------------------------------------------------------------------
def update_fwdInt_state(fwdInt_num,state, FwdInt_dir):
    with open(os.path.join(FwdInt_dir,f"status"),'r') as f:
        lines = f.readlines()
    if fwdInt_num == 1: lines[4] = state + "\n";
    if fwdInt_num == 2: lines[5] = state + "\n";   
    with open(os.path.join(FwdInt_dir,f"status"),'w') as f:
        f.writelines(lines)
    return

def get_fwdInt_state(fwdInt_num,FwdInt_dir):
    with open(os.path.join(FwdInt_dir,f"status"),'r') as f:
        lines = f.read().splitlines()
    if fwdInt_num == 1: return lines[4]
    if fwdInt_num == 2: return lines[5]
    return 0