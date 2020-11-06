'''
Script to test functionality from namelist creation to run and postprocessing.

It is recommended to run this in IPython.
'''

import numpy as np
import matplotlib.pyplot as plt
plt.ion()
import omfit_eqdsk, omfit_gapy
import sys

# Make sure that package home is added to sys.path
import sys
sys.path.append('../')
import aurora

# read in default Aurora namelist
namelist = aurora.default_nml.load_default_namelist()

# Use gfile and statefile in local directory:
geqdsk = omfit_eqdsk.OMFITgeqdsk('example.gfile')
inputgacode = omfit_gapy.OMFITgacode('example.input.gacode')

# save kinetic profiles on a rhop (sqrt of norm. pol. flux) grid
kp = namelist['kin_profs']
kp['Te']['rhop'] = kp['ne']['rhop'] = np.sqrt(inputgacode['polflux']/inputgacode['polflux'][-1])
kp['ne']['vals'] = inputgacode['ne']*1e13 # 1e19 m^-3 --> cm^-3
kp['Te']['vals'] = inputgacode['Te']*1e3  # keV --> eV

# set impurity species and sources rate
imp = namelist['imp'] = 'Ar'
namelist['source_type'] = 'const'
namelist['Phi0'] = 1e24

# Now get aurora setup
asim = aurora.core.aurora_sim(namelist, geqdsk=geqdsk)

# set time-independent transport coefficients (flat D=1 m^2/s, V=-2 cm/s)
D_z = 1e4 * np.ones(len(asim.rvol_grid))  # cm^2/s
V_z = -2e2 * np.ones(len(asim.rvol_grid)) # cm/s

# Here the issue is that I have to copy code from the run_aurora section into here
nz_init = np.zeros((len(self.rvol_grid),int(self.Z_imp+1)))

if D_z.ndim==2:
    # set all charge states to have the same transport
    D_z = np.tile(np.atleast_3d(D_z),(1,1,self.Z_imp+1))  # include elements for neutrals
    V_z = np.tile(np.atleast_3d(V_z),(1,1,self.Z_imp+1))
    
    # unless specified, set transport coefficients for neutrals to 0
    D_z[:,:,0] = 0.0
    V_z[:,:,0] = 0.0

if D_z.ndim==1:
    # D_z was given as time-independent
    D_z = np.tile(np.atleast_3d(D_z[:,None]),(1,1,self.Z_imp+1))  # include elements for neutrals
    V_z = np.tile(np.atleast_3d(V_z[:,None]),(1,1,self.Z_imp+1))
    times_DV = [1.] # dummy, no time dependence


return (len(asim.time_out),  # number of times at which simulation outputs results
            times_DV,
            D_z, V_z, # cm^2/s & cm/s    #(ir,nt_trans,nion)
            asim.par_loss_rate,  # time dependent
            asim.source_rad_prof,# source profile in radius
            asim.S_rates, # ioniz_rate,
            asim.R_rates, # recomb_rate,
            asim.rvol_grid, asim.pro_grid, asim.qpr_grid,
            asim.mixing_radius, asim.decay_length_boundary,
            asim.time_grid, asim.saw_on,
            asim.source_time_history, # source profile in time
            asim.save_time, asim.sawtooth_erfc_width, # dsaw width  [cm, circ geometry]
            asim.wall_recycling,
            asim.source_div_fraction, # divbls [fraction of source into divertor]
            asim.tau_div_SOL_ms * 1e-3, asim.tau_pump_ms *1e-3, asim.tau_rcl_ret_ms *1e-3,  #[s] 
            asim.rvol_lcfs, asim.bound_sep, asim.lim_sep, asim.prox_param,
            nz_init, alg_opt, evolneut)