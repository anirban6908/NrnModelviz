from bmtk.builder.networks import NetworkBuilder
import ateam.sim.setup as sim
import os
from ateamopt.utils import utility
import shutil
import h5py
import numpy as np
import matplotlib.pyplot as plt
from ateamopt.morph_handler import Morph_handler

# function to calculate spike triggered average and covariance for one channels
def cal_STA_STC_onech(data,spikes,win,dt): # data: one_dimensional data
    lwin = win[0]
    hwin = win[1]
    t_STA = np.arange(lwin,hwin,dt)
    
    for i in range(0,len(spikes)):
        indx = range(int(spikes[i,]/dt)+int(lwin/dt),int(spikes[i,]/dt)+int(hwin/dt))   #-1ms, +3ms
        if i==0:
            arr=np.array(data[indx])
        else:
            if len(data)>int(spikes[i,]/dt)+int(hwin/dt):
                arr = np.vstack((arr, np.array(data[indx])))
    return t_STA,np.mean(arr,axis=0),np.std(arr,axis=0)

# function to calculate spike triggered average and covariance for multiple channels
def cal_STA_STC(data,spikes,win,dt): # ecp_data: channels*times; output: channels*ecp_win_times
    
    num_channels = data.shape[0]  # number of channels
    times=int((win[1]-win[0])/dt)
    
    #Calculate STA and STC for multiple channels
    STA = np.array([[0.0]*times]*num_channels)     #spike triggered average
    STC = np.array([[0.0]*times]*num_channels)     #spike triggered standard deviation 

    for i in range(num_channels):
        t_STA,STA[i,:],STC[i,:]=cal_STA_STC_onech(data[i,:],spikes,win,
                 dt)
    
    return (t_STA,STA,STC)



bmtk_config_template = os.path.join(os.getcwd(),"bmtk_config/default_config.json")
morph_path = os.path.join(os.getcwd(),'sim_components/reconstruction.swc')


sim_folder = 'bmtk_simulation'
sm = sim.SimManager.from_template(config_template=bmtk_config_template, overwrite=True,
                                          sim_folder=sim_folder)

network_name = 'snm'
net = NetworkBuilder(network_name)
        
cell_name = 'aa_model'
bmtk_morph_path = os.path.join(os.getcwd(),"bmtk_config/morphologies/{}.swc".\
                               format(cell_name))
utility.create_filepath(bmtk_morph_path)
shutil.copyfile(morph_path,bmtk_morph_path)
morph_handler = Morph_handler(morph_path)
morph_data,morph_apical,morph_axon = morph_handler.\
                            get_morph_coords()
theta_x,theta_y,theta_z = morph_handler.calc_euler_angle\
                                (morph_data,morph_apical)
node_props = {'cell_name' : cell_name,
              'model_type' : 'biophysical',
              'model_template' : 'ctdb:Biophys1.hoc',
              'model_processing' : 'aibs_allactive_ani_directed',
              'dynamics_params' : 'param_aibs_%s.json'%cell_name,
              'morphology': '%s.swc'%cell_name,
              'rotation_angle_xaxis':[theta_x],
              'rotation_angle_yaxis':[theta_y],
              'rotation_angle_zaxis':[theta_z]
              }

net.add_nodes(**node_props)

protocol_name = 'DC_stim'
stim_amp = 200*1e-3
delay=270
duration=1000
total_duration=2270 
dt = 0.1
protocol_dict = {
                protocol_name:{"stimuli":
                    [
                        {
                            "amp": stim_amp,
                            "delay": delay,
                            "duration": duration,
                            "totduration": total_duration,
                        }
                    ]
                }
            }   

sm.add_current_clamp_input(protocol_name, protocol_dict[protocol_name]['stimuli'][0], 
                           loop_delay = 0)

sm.add_network(net)
net.build()
sm.save_network_files(use_abs_paths=True)
sm.add_ecp_report(locs=[[2,0,2]])
sm.sim_time = total_duration
sm.sim_timestep = dt
sm.run_bionet()

ecp_report_file = os.path.join('{}/output'.format(sim_folder),
                   sm.config['reports']['ecp_report']['file_name'])
spike_report_file = os.path.join('{}/output'.format(sim_folder),
                   sm.config['output']['spikes_file'])

with h5py.File(ecp_report_file,'r') as ecp_data:
    ecp_rec0=ecp_data['data']
    ecp_rec0=1e3 * np.array(ecp_rec0)  # Change the unit from mV to uV
Ve = np.transpose(ecp_rec0)   #ecp_data: channels*times
t = np.arange(Ve.shape[1])*dt

fig,ax = plt.subplots()
ax.plot(t,np.squeeze(Ve),color='k')
plt.show()

with h5py.File(spike_report_file,'r') as spike_data:
    spike_times = np.array(spike_data['spikes/timestamps'],
                            dtype=np.float)
win=[-1,3]
[t_STA,oneCol_VeSTA,oneCol_VeSTC] = cal_STA_STC(Ve,spike_times,win,dt)  

oneCol_VeSTA = np.squeeze(oneCol_VeSTA)
oneCol_VeSTC = np.squeeze(oneCol_VeSTC)
oneCol_VeSTA_plus = oneCol_VeSTA + oneCol_VeSTC
oneCol_VeSTA_minus = oneCol_VeSTA - oneCol_VeSTC


fig,ax = plt.subplots()
ax.plot(t_STA,oneCol_VeSTA,color='k')
ax.fill_between(t_STA,oneCol_VeSTA_minus,
                oneCol_VeSTA_plus,alpha=0.2,color='k')
plt.show()