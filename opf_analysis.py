#--V1.0 17/10/2023
#-- Anton

"""
    Openpipe analysis functions
    Contains:
        1) Q-U extraction and loading 
        2) Q-U Classification
        3) State distance
        4) Display
"""

import netCDF4 as nc
import numpy as np
import os
import matplotlib.pyplot as plt
import subprocess

###########################################################################
### 1) Q-U extraction and loading:
# - Code for extraction of Q-U signals from mat_vecXXXX.cdf.dat files 
# - outputed by 'openpipeflow', saving them as .npz and loading them

# -----------------------------------------------------------------------

def convertToMatVec(name,input, run_dir, utils_dir, pressure = False):
    """
    Extracts mat_vec.cdf files from raw stateXXXX.cdf.dat files for openpipe 
    simulation '(run_dir)/(name)' using the appropriate prim2matlab.out 
    program - '(utils_dir)/p2m' with the general extraction script
    '(utils_dir)/mv_exp.sh'. 

    Inputs:
        name(str):       run name
        run_type(str):   type of run to chose appropriate utils file (='R2050_L180')
        input(str arr):  input for prim2matlab.out program
        run_dir(str):    directory where runs are located
        utils_dir(str):  utility directory
        presure(bool):   If False, extracts velocity data, else pressure
    Output:
        1
    Result:
        Cration of 'mat_vec_out' folder, containing the mat_vec.cdf files in
        the simulation directory '(run_dir)/(name)/mat_vec_out'
    """
    rpath =  os.path.join(run_dir,name) #run directory
    
    #create imput file for script from input_str
    with open(os.path.join(rpath,"input.inp"), "w") as f:
        f.writelines(s + '\n' for s in input)
    
    #Copy the extraction script togeter with the apropriate prim2matlab.out file to the run dir 
    if not pressure: #use diffrent script for presure and velocity
        subprocess.run([f"cp {utils_dir}/mv_series_ext.sh {rpath}/mv_series_ext.sh"], shell=True)
    else:
        subprocess.run([f"cp {utils_dir}/mv_pressure_ext.sh {rpath}/mv_series_ext.sh"], shell=True)
    
    subprocess.run([f"cp {utils_dir}/p2m.out {rpath}/prim2matlab.out"], shell=True)
    #Run the extraction script
    subprocess.run([os.path.join(rpath,"mv_series_ext.sh")], shell=True)
    #subprocess.run(f"rm {rpath}/prim2matlab.out", shell=True)
    #subprocess.run(f"rm {rpath}/mv_series_ext.sh", shell=True)

    return 

# -----------------------------------------------------------------------
def load_QU(path):
#load existing Q-U file
    with np.load(path) as data:
        Q,U,z,t = data['Q'],data['U'],data['z'], data['t']
    return Q,U,z,t
#------------------------------------------------------------------------    

def count_cdf_files(folder, end_with = ".cdf"):
    file_list = []
    for f in os.listdir(folder): 
        if f.endswith(end_with):
            file_list.append(f)
    return len(file_list)

#------------------------------------------------------------------------    
def compute_QU(mat_vec_dir):
#extracting and compute Q-U files from .cdf.dat files
    #Count number of snaphots
    N = count_cdf_files(mat_vec_dir)
    print(f'Total of {N} cdf.dat files')
    #load meta data
    ds = nc.Dataset(os.path.join(mat_vec_dir,f'mat_vec{0:04}.cdf'))
    metadata = ds.__dict__
    n_xd,n_yd,n_zd = ds.dimensions['xd'].size,ds.dimensions['yd'].size,ds.dimensions['zd'].size
    x,y,z = ds['x'][:],ds['y'][:],ds['z'][:]
    [Q,U] = np.zeros(shape=(2,N,n_zd))
    #load Q fields for all runs
    t  = np.zeros(N)
    for idx in range(N):
        fn =  os.path.join(mat_vec_dir,f'mat_vec{idx:04}.cdf')
        ds = nc.Dataset(fn) 
        [u,v,w]  = ds['A'][:,:,:,:]
        Q[idx,:] = np.trapz(np.trapz(u**2 + v**2,y,axis = 2),x,axis = 1); 
        U[idx,:] = w[:,n_xd//2,n_yd//2]
        t[idx]   = ds.__dict__['t']
        print(f'Loading {fn} t = {t[idx]}')
    return Q,U,z,t
#------------------------------------------------------------------------
def get_QU(name,run_dir = '..',post_dir = './post', echo = False):
    #get Q-U file of simulation by loading or extracting from .cdf.dat files
    mat_vec_dir      = f'{run_dir}/{name}/mat_vec_out';  #directory where the mat_vec_files are located
    Q_save_file_name = f'{post_dir}/{name}_QU.npz'       #path to saved QU file
    if os.path.isfile(Q_save_file_name): #check if saved file exist
        if echo: print(f'Loading Q file {Q_save_file_name}')
        Q,U,z,t = load_QU(Q_save_file_name)
    else:
        if echo: print(f'Computing Q-U file in {mat_vec_dir}')
        Q,U,z,t = compute_QU(mat_vec_dir)
        if echo: print(f'Saving Q-U file as {Q_save_file_name}')
        np.savez(Q_save_file_name,Q = Q,U = U, z = z, t = t) #save Q-U file
    return Q,U,z,t

#############################################################################

### 2) Q-U Classification:
# - Code for classification of flow state based on Q-U signals

#------------------------------------------------------------------------
def find_mu(q,x):
    # Finds CoM using the phase of the first Fourier mode on a q signal 
    L = x[-1] + x[1]
    N = len(q)
    c = (L/N)*np.sum(q*np.exp(2j*np.pi*(x/L)))
    mu = x[int(np.angle(c)*N/(2*np.pi))]       
    return mu
#------------------------------------------------------------------------
def find_mu_over_t(Q,x):
    # Finds CoM using the phase of the first Fourier mode on a q signal 
    # For each time point q = Q[i,:]
    T = Q.shape[0]
    mu = np.zeros(T)
    for i in range(T):
        mu[i] = find_mu(Q[i,:],x)
    return mu
#-----------------
def find_mu_vel(mu,z,t):
    #finds the velocity of CoM
    L = np.max(z)
    v_mu = np.gradient(mu,t)
    valid_i   = np.where(np.abs(v_mu) < L/2) #exclude transition points
    v0   = np.mean(v_mu[valid_i])
    return v0
#-------------------------------------------------------
def CenterOfMassAlign(F,z):
    # Align to CoM of f
    mu = find_mu_over_t(F,z)
    F_com = np.zeros_like(F)
    for i in range(F.shape[0]):
        F_com[i] = np.roll(F[i], F.shape[1]//2 - int(round(mu[i])))
    return F_com
#------------------------------------------------------------------------
def CenterOfMassAlign_QU(Q,U,z):
    # Align to CoM of f
    mu = find_mu_over_t(Q,z)
    Q_com = np.zeros_like(Q)
    U_com = np.zeros_like(U)
    for i in range(Q.shape[0]):
        Q_com[i] = np.roll(Q[i], Q.shape[1]//2 - int(round(mu[i])))
        U_com[i] = np.roll(U[i], U.shape[1]//2 - int(round(mu[i])))
    return Q_com, U_com
#------------------------------------------------------------------------
def GalileanAlign_QU(Q,U,z,t):
    # Align to to the velocity of CoM
    mu   = find_mu_over_t(Q,z)
    v    = find_mu_vel(mu,z,t)
    #inital shift to CoM at t = 0 
    Q_s = np.roll(Q, Q.shape[1]//2 - int(round(mu[0])),axis = 1)
    U_s = np.roll(U, U.shape[1]//2 - int(round(mu[0])),axis = 1)
    
    Q_com = GalShift(Q_s,z,t,v) 
    U_com = GalShift(U_s,z,t,v) 
    return Q_com, U_com
#------------------------------------------------------------------------
def GetMoment(F,z,m):
    #Get m-th moment of F for each time point f=F[i,:]
    T = F.shape[0]
    F_com = CenterOfMassAlign(F,z)
    z_com = z - z[z.shape[0]//2]
    moment = np.zeros(T)
    for i in range(T):
        moment[i] = np.trapz((F_com[i]*z_com**m),z_com)/np.trapz(F_com[i],z_com)
    return moment
#------------------------------------------------------------------------
def GetPuffAndGaps(f_com,z):
    ## Find size of patches and gaps from f_com field (f centered to CoM)
    dz = np.abs(z[1] - z[0])
    dg = np.diff(f_com,prepend = 0)
    rise_idx = np.where(dg ==  1)[0]
    fall_idx = np.where(dg == -1)[0]
    puffs = np.array([])
    gaps = np.array([])
    for r_i in rise_idx:
        for f_i in fall_idx:
            if f_i> r_i:
                puffs = np.append(puffs, f_i-r_i);
                break;
    for f_i in fall_idx:
        for r_i in rise_idx:
            if r_i> f_i:
                gaps = np.append(gaps,r_i-f_i);
                break;
    return puffs*dz, gaps*dz
#------------------------------------------------------------------------
def puffCount(F_com,z,min_gap = 1.0, small_gap = 0.5):
    #find number of turbulent patches for each snapshot f_com = F_com[i,:] in CoM alined F (F_com)
    #small gap indicates a gap to small to be considered two-puffs but large enogt to not be considered one-puff
    T = F_com.shape[0]
    count = np.zeros(T)
    for i in range(T):
        _, gaps =  GetPuffAndGaps(F_com[i],z)
        count[i] = np.sum(gaps>min_gap) + 1
        count[i] += 0.1*np.sum((gaps>small_gap)*(gaps<min_gap))
        #count 
    return count  
#------------------------------------------------------------------------    
def classifyQ(Q,z,state_def):
    #classify the state at each time step
    F_th      = state_def['F_th']      #Treshold of Q: F:=Q>F_th
    min_gap   = state_def['min_gap']   #Minimum space [in units of z] to be considerd as gap
    small_gap = state_def['small_gap'] #Samaller value of gap, for which no longer a single puff (not yet two)
    #Tresholds on Ft:= trapz(F,z), such that Ft_lim[0] > Ft > Ft_lim[1]
    p0_Ft_lim = state_def['p0_Ft_lim'] #Decay 
    p1_Ft_lim = state_def['p1_Ft_lim'] #Single puff
    p2_Ft_lim = state_def['p2_Ft_lim'] #Two puffs
    p3_Ft_lim = state_def['p3_Ft_lim'] #More then two puffs
    
    F     = np.array(Q>F_th,dtype = np.float32)   #Turbulent patches
    F_com = CenterOfMassAlign(F,z)                #Center of mass F
    Pc    = puffCount(F_com,z, min_gap = min_gap, small_gap = small_gap) #Number of puffs
    Ft    = np.trapz(F,z,axis = 1)                #Turbulent fraction
    
    cls_0p = (Ft>p0_Ft_lim[0])*(Ft<=p0_Ft_lim[1])         #Decay      cls = -1
    cls_1p = (Ft>p1_Ft_lim[0])*(Ft<=p1_Ft_lim[1])*(Pc==1.0) #One puff   cls =  1
    cls_2p = (Ft>p2_Ft_lim[0])*(Ft<=p2_Ft_lim[1])*(Pc==2.0) #Two puffs  cls =  2
    cls_3p = (Ft>p3_Ft_lim[0])*(Ft<=p3_Ft_lim[1])*(Pc>2.0)  #More ouffs cls =  3

    cls = (-1)*cls_0p + (1)*cls_1p + (2)*cls_2p + (3)*cls_3p #classification
    return cls
#------------------------------------------------------------------------
def getClassDetailed(Q,z,state_def):
    #extract classification series
    F     = np.array(Q>state_def["F_th"],dtype = np.float32) 
    F_com = CenterOfMassAlign(F,z)
    Q_com = CenterOfMassAlign(Q,z)
    Ft    = np.trapz(F,z,axis = 1)
    Pc    = puffCount(F_com,z, min_gap = state_def["min_gap"],small_gap = state_def["small_gap"])
    cls   = classifyQ(Q,z,state_def)
    return Ft,Pc,cls,F_com,Q_com
###########################################################################

### 3) State distance:
# -Computing the distance between flow states stored in mat_vecXXXX.cdf.dat files 

#------------------------------------------------------------------------
def Norm_S(s):
    #compute amplitude of velocity field [u,v,w]
    return np.sqrt(s[0]**2 + s[1]**2 + s[2]**2)
#------------------------------------------------------------------------
def HP_flow(x,y,z):
    #base HP flow profile
    U  = np.zeros(shape=(3,np.size(z),np.size(y),np.size(x)))
    zm,ym,xm = np.meshgrid(z,x,y,indexing = 'ij')
    U[2,:,:,:] = 1 - (xm**2 + ym**2)
    U[U<0] = 0
    return U
#------------------------------------------------------------------------
def computeQU_fromS(s0,x,y,z):
    #compute Q from velocity [u,v,w]
    n_xd,n_yd = np.size(x),np.size(y)
    [u,v,w] = s0
    Q = np.trapz(np.trapz(u**2 + v**2,y,axis = 2),x,axis = 1)
    U = w[:,n_xd//2,n_yd//2]
    return Q,U
#------------------------------------------------------------------------
def centState(s,x,y,z):
    # center state to the CoM of Q
    Q,_ = computeQU_fromS(s,x,y,z)
    mu = int(np.round(find_mu(Q,z)))
    return np.roll(s,len(z)//2 - mu,axis = 1)
#------------------------------------------------------------------------
def int_s_XYZ(f,x,y,z):
    #integrate over the entire domain
    return np.trapz(np.trapz(np.trapz(f,x=z,axis = 0),x=y,axis = 0),x=x,axis = 0)
#------------------------------------------------------------------------
# def distanceL2(s1,s2,x,y,z, align = True):
#     #Compute L2 distance between velocity fields s1,s2
#     # old verison
#     if align:
#         a1,a2 = centState(s1,x,y,z),centState(s2,x,y,z) 
#     else:
#         a1,a2 = s1,s2
#     Norm = np.sqrt(int_s_XYZ(Norm_S(a1),x,y,z))*np.sqrt(int_s_XYZ(Norm_S(a2),x,y,z))
#     return int_s_XYZ(Norm_S(a1 - a2),x,y,z)/Norm

def distanceL2(s1,s2,x,y,z, align = True):
    #Compute L2 distance between velocity fields s1,s2
    # Update version with z shifting
    if align:
        a1,a2 = centState(s1,x,y,z),centState(s2,x,y,z)
        Norm = np.sqrt(int_s_XYZ(Norm_S(a1),x,y,z))*np.sqrt(int_s_XYZ(Norm_S(a2),x,y,z))
        D,_ =  meas_dist_z_inv(s1, s2, x, y, z)
    else:
        a1,a2 = s1,s2
        Norm = np.sqrt(int_s_XYZ(Norm_S(a1),x,y,z))*np.sqrt(int_s_XYZ(Norm_S(a2),x,y,z))
        D    = int_s_XYZ(Norm_S(a1 - a2),x,y,z)
    return D/Norm
#------------------------------------------------------------------------
def meas_dist_z_inv(u1_in, u2_in, x, y, z , z_shift_max = 5,weights = np.array([1.0,1.0,1.0])):
    #compute a 'z shift invariant' distance between two velocity fields fields u1(z,x,y) and u2(x,y,z)
    #the fields u1 and u2 are assumed to be in the shape (3, nz,nx,ny)
    #z_shift_max defines the maximum shift (in z units) considered
    #weights = the 'weigth' given to each component of u for distance computation
    
    u1 = np.copy(u1_in)
    u2 = np.copy(u2_in)
    #allpy weight to each component
    for k in range(3):
        u1[k] = weights[k]*u1[k]
        u2[k] = weights[k]*u2[k]
    
    #create a vector of index shifts 
    z_i_shift_max = int(np.ceil(z_shift_max/(z[1] - z[0])))
    shift_vector = np.arange(-z_i_shift_max, z_i_shift_max+1)
    N_shifts = len(shift_vector)
    
    #creates a matrix with a row for each shifted u2
    u2_shifts = np.expand_dims(u2, 0)
    u2_shifts = np.repeat(u2_shifts, N_shifts, axis=0)
    for i in range(N_shifts):
        u2_shifts[i] = np.roll(u2_shifts[i], shift_vector[i], axis = 1)

    #compute distance for each shift
    u_diff = np.sum((u2_shifts - u1)**2, axis = 1)
    dist_shifts = np.zeros(N_shifts)
    dist_shifts = np.trapz(np.trapz(np.trapz(u_diff, x = z, axis = 1), x = x, axis = 1), x = y, axis = 1)
    
    #minimal distance is chosen
    dist  = np.min(dist_shifts)
    shift = shift_vector[np.argmin(dist_shifts)]

    return dist,shift
#------------------------------------------------------------------------
def getState(fileName,dir):
    #load specific state file "fileName" from directory "dir"
    fn =  os.path.join(dir,fileName)
    ds = nc.Dataset(fn)
    return ds['A'][:]
#------------------------------------------------------------------------
def getGrid(fileName,dir):
    #load grid 
    fn =  os.path.join(dir,fileName)
    ds = nc.Dataset(fn)
    x,y,z = ds['x'][:],ds['y'][:],ds['z'][:]
    return x,y,z
#------------------------------------------------------------------------
def distanceL2_mv_files(file_1,file_2,directory,sub_HPf = False, align = True):
    #Compute L2 distance between velocity fields saved in files file_1,file_2 in dir
    #sub_HPf=True if HP flow need's to be subtracted. 
    x,y,z = getGrid(file_1,directory)
    s1 = getState(file_1,directory) 
    s2 = getState(file_2,directory) 
    if sub_HPf:
        U = HP_flow(x,y,z) 
        s1 = s1 - U
        s2 = s2 - U
    return distanceL2(s1,s2,x,y,z,align = align)

#############################################################################
### 4) Display:
#------------------------------------------------------------------------
def GalShift(F,z,t,v):
    # Galilean shift of the signal
    Fs = np.zeros_like(F)
    dt,dz = np.diff(t)[0],np.diff(z)[0]
    for i in range(F.shape[0]):
        Fs[i,:] = np.roll(F[i,:], -int(v*i*dt/dz))  
    return Fs
#------------------------------------------------------------------------
def classificationPlot(name,state_def,run_dir, post_dir,fig_dir):
    """
    Creates classification figure for run '(run_dir)/(name)' based on state files
    for which mat_vec.cdf files have been extracted. The classification is done with 
    tresholds defined by 'state_def'.

    Inputs:
        name(str):       run name
        state_def(dict): deffention of tresholds for the diffrent states
        run_dir(str):    directory where runs are located
        post_dir(str):   directory where the post-proccesing data (Q-U fields) are saved
        fig_dir(str):    directory for the output figures
    Output:
        1
    Result:
        Cration classification figure (fig_dir)/(name)_class.png
    """    
    Q,_,z,t = get_QU(name,run_dir = run_dir, post_dir = post_dir)
    Ft,Pc,cls,F_com,Q_com = getClassDetailed(Q,z,state_def)
    zm, tm = np.meshgrid(z, t)
    fig, ax = plt.subplots(figsize = (6,6),nrows = 4, sharex = True)
    #Turbulance
    #ax[0].pcolormesh(tm,zm,F_com,shading='auto', cmap = 'hot', alpha = 1)
    ax[0].pcolormesh(tm,zm,np.log(Q_com),shading='auto', cmap = 'cool'  , alpha = 1)
    ax[0].pcolormesh(tm,zm,F_com        ,shading='auto', cmap = 'binary', alpha = 0.2)
    ax[0].contour(tm,zm,F_com, colors = 'black',linewidths = 0.5, alpha = 1)
    
    ax[0].set_ylabel('$z [R]$')
    
    #Turubulent fraction
    ax[1].plot(t,Ft)
    
    max_Ff = np.max([np.max(Ft),state_def['p3_Ft_lim'][0]+10])
    ax[1].fill_between(t, state_def['p0_Ft_lim'][0], y2=state_def['p0_Ft_lim'][1] , color = 'red'    , alpha = 0.1, label = 'decay')
    ax[1].fill_between(t, state_def['p1_Ft_lim'][0], y2=state_def['p1_Ft_lim'][1] , color = 'blue'   , alpha = 0.1, label = 'x1 puff')
    ax[1].fill_between(t, state_def['p2_Ft_lim'][0], y2=state_def['p2_Ft_lim'][1] , color = 'green'  , alpha = 0.1, label = 'x2 puff')
    ax[1].fill_between(t, state_def['p3_Ft_lim'][0], y2=max_Ff                    , color = 'magenta', alpha = 0.1, label = '>2 puff')
    
    ax[1].set_ylabel(f"Total Turb. \n (th = {state_def['F_th']:.0e})")
    ax[1].grid()
    ax[1].legend(ncols = 4)
    
    #Puff count
    ax[2].plot(t,Pc)
    
    ax[2].fill_between(t, 0.75, y2=1.25 , color = 'blue', alpha = 0.1  , label = 'x1 puff')
    ax[2].fill_between(t, 1.75, y2=2.25 , color = 'green', alpha = 0.1 , label = 'x2 puff')
    
    ax[2].set_ylabel(f"Puff Count \n (gap = {state_def['min_gap']})")
    ax[2].set_ylim([0,3])
    ax[2].grid()
    #ax[2].legend(ncols = 2)
    
    #Classification
    ax[3].plot(t,cls)
    
    ax[3].fill_between(t, -1.25, y2=-0.75 , color = 'red'    , alpha = 0.1, label = 'decay')
    ax[3].fill_between(t,  0.75, y2= 1.25 , color = 'blue'   , alpha = 0.1, label = 'x1 puff')
    ax[3].fill_between(t,  1.75, y2= 2.25 , color = 'green'  , alpha = 0.1, label = 'x2 puff')
    ax[3].fill_between(t,  2.75, y2= 3.25 , color = 'magenta', alpha = 0.1, label = '>2 puff')
    
    ax[3].set_ylabel('Classification')
    ax[3].set_ylim([-1,3])
    ax[3].grid()
    #ax[3].legend(ncols = 4)
    
    ax[3].set_xlabel('$t [R/U_{cl}]$')
    
    plt.tight_layout()
    fig.suptitle(f'Classification {name}', y = 1.01)
    plt.savefig(os.path.join(fig_dir,f"{name}_class.png"),bbox_inches='tight')
    return cls,t
#------------------------------------------------------------------------
def GenDistPlots(FwdInt_dir,dist_fwdInt,new_idx,state_def,L_bs,L_fwd,fig_dir):
    #distance plots for fwd integration
    #get QU fields
    Q1,U1,z1,t1 = get_QU("fwdInt_1p",FwdInt_dir, FwdInt_dir)
    Q2,U2,z2,t2 = get_QU("fwdInt_2p",FwdInt_dir, FwdInt_dir)
    
    #shift time
    t1 = t1 - t1[0]
    t2 = t2 - t2[0]
    
    Ft1,Pc1,cls1,F1_com,Q1_com = getClassDetailed(Q1,z1,state_def)
    Ft2,Pc2,cls2,F2_com,Q2_com = getClassDetailed(Q2,z2,state_def)
    
    
    fig, ax = plt.subplots(figsize = (6,6),nrows = 3, sharex = True)
    
    #distance
    ax[0].plot(t1,dist_fwdInt,'.', label = '$d(1p,2p)$')
    ax[0].plot(t1[new_idx],dist_fwdInt[new_idx],'*', label = 'selected')
    
    ax[0].axhline(L_bs,  color = 'green', label = '$L_{bs}$' )
    ax[0].axhline(L_fwd, color = 'red'  , label = '$L_{fws}$')
    
    ax[0].set_ylabel(f"$L_2$ dist.")
    ax[0].grid()
    ax[0].legend(ncols = 4)
    
    # Turubulent fraction
    ax[1].plot(t1,Ft1, label = 'fwdInt_1p')
    ax[1].plot(t2,Ft2, label = 'fwdInt_2p')
    
    max_Ff = np.max([np.max(Ft2),state_def['p3_Ft_lim'][0]+10])
    ax[1].fill_between(t2, state_def['p0_Ft_lim'][0], y2=state_def['p0_Ft_lim'][1] , color = 'red'    , alpha = 0.1)
    ax[1].fill_between(t2, state_def['p1_Ft_lim'][0], y2=state_def['p1_Ft_lim'][1] , color = 'blue'   , alpha = 0.1)
    ax[1].fill_between(t2, state_def['p2_Ft_lim'][0], y2=state_def['p2_Ft_lim'][1] , color = 'green'  , alpha = 0.1)
    ax[1].fill_between(t2, state_def['p3_Ft_lim'][0], y2=max_Ff                    , color = 'magenta', alpha = 0.1)
    
    ax[1].set_ylabel(f"Total Turb. \n (th = {state_def['F_th']:.0e})")
    ax[1].grid()
    ax[1].legend(ncols = 4)
    
    #Classification
    ax[2].plot(t1,cls1)
    ax[2].plot(t2,cls2)
    
    ax[2].fill_between(t2, -1.25, y2=-0.75 , color = 'red'    , alpha = 0.1)
    ax[2].fill_between(t2,  0.75, y2= 1.25 , color = 'blue'   , alpha = 0.1)
    ax[2].fill_between(t2,  1.75, y2= 2.25 , color = 'green'  , alpha = 0.1)
    ax[2].fill_between(t2,  2.75, y2= 3.25 , color = 'magenta', alpha = 0.1)
    
    ax[2].set_ylabel('Classification')
    ax[2].set_ylim([-1,3])
    ax[2].grid()
    ax[2].set_xlabel('$t [R/U_{cl}]$')
    
    plt.tight_layout()
    fig.suptitle(f'Forward integration', y = 1.01)
    plt.savefig(os.path.join(fig_dir,f"fwdInt_dist.png"),bbox_inches='tight')
    plt.close(fig)
    return
#------------------------------------------------------------------------
def fwdInt_initalFinalStatesPlot(FwdInt_dir,dist_fwdInt,new_idx):
    Q1,U1,z1,t1 = get_QU("fwdInt_1p",FwdInt_dir, FwdInt_dir)
    Q2,U2,z2,t2 = get_QU("fwdInt_2p",FwdInt_dir, FwdInt_dir)
    Q1_com,U1_com = CenterOfMassAlign_QU(Q1,U1,z1)
    Q2_com,U2_com = CenterOfMassAlign_QU(Q2,U2,z2)
    fig, ax = plt.subplots(figsize = (9,3),ncols = 3,nrows = 2, sharex = True)
    plot_names = ['inital', 'selected','final']
    for i,idx in enumerate([0,new_idx,-1]):
        ax[0,i].plot(z1,Q1[idx] , label = '$Q_1$')
        ax[0,i].plot(z2,Q2[idx] , label = '$Q_2$')
        ax[1,i].plot(z1,U1[idx] , label = '$U_1$')
        ax[1,i].plot(z2,U2[idx] , label = '$U_2$')   
        
        ax[0,i].set_title(f'{plot_names[i]}: t = {t1[idx]:.1f} \n d = {dist_fwdInt[idx]:.2e} ')
        ax[0,i].set_ylim([-0.001,0.025])
        ax[0,i].legend()
        
        #ax[1,i].set_ylim([-0.5,  0.05])
        ax[1,i].legend()
    plt.tight_layout()
    plt.savefig(os.path.join(FwdInt_dir,f"fwdInt_QU_states.png"),bbox_inches='tight')
    plt.close(fig)
    return

#############################################################################
### 5) Fowrard Trajectory analysis:
#------------------------------------------------------------------------
def load_fwdTrajectory(main_dir,state=1):
    #loads the forawrds drajectories saved in main_dir, combined at the divergence time as set in 
    #the status faile in each fwdInt_XXXXX directory. 
    # state indicaates wheter to take 1p (state = 1) or 2p (state = 2) data
    print(f'Loading {state} puff data from {main_dir}')
    #get ordered FwdInt dirs
    fwdDirs = sorted([name for name in os.listdir(main_dir) if name.startswith('FwdInt_')])
    
    #load meta data and dimensions
    mat_vec_dir_meta =  os.path.join(main_dir,fwdDirs[0],'fwdInt_1p/mat_vec_out') #regular variant
    ds = nc.Dataset(os.path.join(mat_vec_dir_meta,f'mat_vec{0:04}.cdf'))
    metadata = ds.__dict__
    
    # load metadata
    n_xd,n_yd,n_zd = ds.dimensions['xd'].size,ds.dimensions['yd'].size,ds.dimensions['zd'].size
    x,y,z = ds['x'][:],ds['y'][:],ds['z'][:]
    Re = ds.Re
    
    for idx in range(0,len(fwdDirs)):
        #loading
        mat_vec_dir =  os.path.join(main_dir,fwdDirs[idx],f'fwdInt_{state}p/mat_vec_out') #regular variant
        mat_presure_dir =  os.path.join(main_dir,fwdDirs[idx],f'fwdInt_{state}p/mat_vec_pressure') #regular variant
        #mat_vec_dir =  os.path.join(main_dir,fwdDirs[idx],'fwdInt_1p/mat_vec_hiRes') #high resolution
        print(f'{idx} - Selected {mat_vec_dir} ')
    
        #find the last fwdInt step before bisection
        with open(os.path.join(main_dir,fwdDirs[idx],f"status")) as f:
            lines = f.read().splitlines()
            New_1p_ic = lines[6]  #new 1 puff state
            last_idx_1p = int(New_1p_ic.split('/')[1][5:9]) #the index of the state
    
        #set the step count 
        N_T_max = last_idx_1p #set the last posible time point
        N_T_load = np.min([count_cdf_files(mat_vec_dir), N_T_max])
        print(f'{idx} - Total of {count_cdf_files(mat_vec_dir)} cdf.dat files, taking {N_T_load} for trajectory')
        
        #load velocity snapshots
        U_field_load = np.zeros(shape = (N_T_load,3,n_zd,n_xd,n_yd))
        t_load = np.zeros(N_T_load)
        for n_mat in range(N_T_load):
            fn =  os.path.join(mat_vec_dir,f'mat_vec{n_mat:04}.cdf')
            ds = nc.Dataset(fn) 
            [u,v,w]  = ds['A'][:,:,:,:]
            U_field_load[n_mat] = np.array([w,u,v])
            t_load[n_mat] = ds.__dict__['t']
        #print(f'{fwdDirs[idx]}: Imported velocity ({U_field_load.shape}); T = {t_load[-1] - t_load[0]}')
        
        P_field_load = np.zeros(shape = (N_T_load,1,n_zd,n_xd,n_yd))
        if os.path.isdir(mat_presure_dir):
            for n_mat in range(N_T_load):
                fn =  os.path.join(mat_presure_dir,f'mat_vec{n_mat:04}.cdf')
                ds = nc.Dataset(fn) 
                pres  = ds['A'][:,:,:,:]
                P_field_load[n_mat] = np.array(pres) 
        else:
            print(f'{idx} - Pressure data missing in {fwdDirs[idx]}');
        #merge fields
        if idx == 0:
            U_field = np.copy(U_field_load)
            P_field = np.copy(P_field_load)
            t = np.copy(t_load)
        else:
            U_field = np.append(U_field, U_field_load[1:], axis = 0)
            P_field = np.append(P_field, P_field_load[1:], axis = 0)
            t = np.append(t,t_load[1:] + t[-1]) #make times continuse
            
    return U_field, P_field, t , x, y, z , metadata
#------------------------------------------------------------------------
def align_to_front_UP(U_field,P_field,x,y,z,q_th = 0.002):
    #align U and P fields to front defined by q>q_th
    q_time = np.trapz(np.trapz((U_field[:,1]**2 + U_field[:,2]**2), axis = 2, x = x), axis = 2, x = y)
    N_T = q_time.shape[0]
    #upstream from pos
    front_z   = np.zeros(N_T, dtype = np.float64);
    front_idx = np.zeros(N_T, dtype = int);
    for nn in range(N_T):
        #first shift tom maximum and then find fron (to correct for periodic bc)
        shift_to_center = np.argmax(q_time[nn])
        q_cent = np.roll(q_time[nn], len(z)//2 - shift_to_center)
        front_idx[nn] = (int(np.where(q_cent > q_th)[0][0]) - (len(z)//2 - shift_to_center))%len(z)
        front_z[nn] = z.data[front_idx[nn]] 
    #galilean shift
    U_align = np.zeros_like(U_field)
    P_align = np.zeros_like(P_field)
    #shift to front
    for nn in range(0,N_T):
        shift = len(z)//2 - front_idx[nn]
        U_align[nn] = np.roll(U_field[nn], shift, axis = 1)
        P_align[nn] = np.roll(P_field[nn], shift, axis = 1)
    return U_align, P_align
#------------------------------------------------------------------------
def align_to_front(U_field,x,y,z,q_th = 0.002):
    #align U and P fields to front defined by q>q_th
    q_time = np.trapz(np.trapz((U_field[:,1]**2 + U_field[:,2]**2), axis = 2, x = x), axis = 2, x = y)
    N_T = q_time.shape[0]
    #upstream from pos
    front_z   = np.zeros(N_T, dtype = np.float64);
    front_idx = np.zeros(N_T, dtype = int);
    for nn in range(N_T):
        #first shift tom maximum and then find fron (to correct for periodic bc)
        shift_to_center = np.argmax(q_time[nn])
        q_cent = np.roll(q_time[nn], len(z)//2 - shift_to_center)
        front_idx[nn] = (int(np.where(q_cent > q_th)[0][0]) - (len(z)//2 - shift_to_center))%len(z)
        front_z[nn] = z.data[front_idx[nn]] 
    #galilean shift
    U_align = np.zeros_like(U_field)
    #shift to front
    for nn in range(0,N_T):
        shift = len(z)//2 - front_idx[nn]
        U_align[nn] = np.roll(U_field[nn], shift, axis = 1)
    return U_align