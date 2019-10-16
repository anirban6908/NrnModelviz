from ateamopt.utils import utility
from collections import defaultdict
import re
import numpy as np
import matplotlib as mpl
   
def param_dict_to_AIBS_param(param_dict,ena=53,ek=-107,\
                             temp=34,v_init= -80,
                             expand_params=True,
                             no_apical = False):

    pass_params = ['e_pas','g_pas','cm','Ra']
    bpopt_section_map = utility.bpopt_section_map
    bpopt_section_map_inv = utility.bpopt_section_map_inv

    aibs_format_params = defaultdict(list)
   
   
    
    aibs_format_params['passive'].append({'ra' : param_dict['Ra.all']})
    aibs_format_params['fitting'].append({
            'junction_potential' : -14.0,
             'sweeps' : []
             })

    ena_sectionlist = []
    ek_sectionlist = []

    sect_reject_list = ['all']
    if no_apical:sect_reject_list.append('apical')
    
    if expand_params:
        # If section == 'all' distribute to all section names

        param_all_entries = {}
        param_del_entries = []
        sect_reject_list_bpopt = [bpopt_section_map_inv[sect_] for sect_ in sect_reject_list]
        temp_sect_map = [sect for sect in bpopt_section_map_inv \
                         if sect not in sect_reject_list_bpopt]

        for param,val in param_dict.items():
            param_name,sect = param.split('.')
            if sect == 'all':
                for sec_ in temp_sect_map:
                    param_all_entries.update({'%s.%s'%(param_name,sec_):\
                                              val})
                param_del_entries.append(param)

        # delete the all entries and repopulate with the four section names
        param_dict = utility.remove_entries_dict(param_dict,param_del_entries)
        param_dict.update(param_all_entries)


    for param,val in param_dict.items():
        param_name,sect = param.split('.')
        if param_name in pass_params:
           mech = ''

        else:
            mech = param_name.split('_',1)[1]
            
        if bool(re.search('Na',mech,re.IGNORECASE)):
            ena_sectionlist.append(sect)
        elif  bool(re.search('K',mech,re.IGNORECASE)):
            ek_sectionlist.append(sect)
        aibs_format_params['genome'].append(
                {
                  'section' : bpopt_section_map_inv[sect],
                  'name'    : param_name,
                  'value'   : str(val),
                  'mechanism': mech
                })
    
        
    ena_sectionlist = list(set(ena_sectionlist))
    ek_sectionlist = list(set(ek_sectionlist))
    
    
    temp_sect_map = [sect for sect in bpopt_section_map_inv.keys() \
                         if sect not in sect_reject_list]

    erev_list = []
    for sect_ in temp_sect_map:
        temp_dict = {}
        if sect_ in ena_sectionlist:
            temp_dict['ena'] = utility.rev_potential['ena']

        if sect_ in ek_sectionlist:
            temp_dict['ek'] = utility.rev_potential['ek']

        if bool(temp_dict):
            temp_dict['section'] = bpopt_section_map_inv[sect_]
            erev_list.append(temp_dict.copy())


    aibs_format_params['conditions'].append({'celsius' : temp,
                      "erev": erev_list,
                      "v_init": v_init
                      })
  
    return aibs_format_params

def AIBS_param_to_bpopt_param(AIBS_param,ena=53,ek=-107,\
                             temp=34,v_init= -80,):

    bpopt_section_map = utility.bpopt_section_map
    model_params_bpopt = list()
    mechs_bpopt = defaultdict(list)
    
    ena_sectionlist = []
    ek_sectionlist = []
    
    for key, values in AIBS_param.items():
        if key == 'genome':
            for j in range(len(values)):
                 iter_dict_release = {'param_name':AIBS_param[key][j]['name']}
                 iter_dict_release['sectionlist'] = bpopt_section_map\
                             [AIBS_param[key][j]['section']]
                 iter_dict_release['type'] = 'section'
                 iter_dict_release['value'] = float(AIBS_param[key][j]['value'])
                 iter_dict_release['dist_type'] = 'uniform'
                 if AIBS_param[key][j]['mechanism'] != '':
                     iter_dict_release['mech'] = AIBS_param[key][j]['mechanism']
                     iter_dict_release['type'] = 'range'
                     mechs_bpopt[bpopt_section_map\
                         [AIBS_param[key][j]['section']]].append(\
                                 AIBS_param[key][j]['mechanism'])
                 else:
                     mechs_bpopt[bpopt_section_map\
                         [AIBS_param[key][j]['section']]].append('pas')
                
                 model_params_bpopt.append(iter_dict_release)
                 if bool(re.search('Na',AIBS_param[key][j]['mechanism'],
                                   re.IGNORECASE)):
                     ena_sectionlist.append(bpopt_section_map\
                                    [AIBS_param[key][j]['section']])
                 elif  bool(re.search('K',AIBS_param[key][j]['mechanism'],
                                      re.IGNORECASE)):
                    ek_sectionlist.append(bpopt_section_map\
                                    [AIBS_param[key][j]['section']])
                
                
    ena_sectionlist = list(set(ena_sectionlist))
    ek_sectionlist = list(set(ek_sectionlist))
    rev_potential = {'ena' : ena_sectionlist,'ek' : ek_sectionlist}
    
    
    for rev,sect_list in rev_potential.items():
        for sect in sect_list:
            iter_dict_release =  {'param_name':rev, 'sectionlist':sect, 
                                  'dist_type': 'uniform', 'type':'section'}
            if rev == 'ena':
                iter_dict_release['value'] = ena
            elif rev == 'ek':
                iter_dict_release['value'] = ek
            model_params_bpopt.append(iter_dict_release)
    
    model_params_bpopt.append({"param_name": "celsius",
                 "type": "global","value": temp})
    model_params_bpopt.append({"param_name": "v_init",
                 "type": "global","value": v_init})

    for mech_sect,mech_list in mechs_bpopt.items():
        mechs_bpopt[mech_sect] = list(set(mech_list))
    return model_params_bpopt,mechs_bpopt

def get_morph_points(morph_obj,theta,axis_of_rot,reject_axon=True):
        
        all_x, all_y, all_z = [], [], []
        
        for n in morph_obj.morph.compartment_list:
            if reject_axon and n['type'] == 2:
                continue
            nx,ny,nz = morph_obj.shift_origin(np.array([n['x'],n['y'],n['z']]))
            [nx_rot,ny_rot,nz_rot] = morph_obj.rotate3D_point([nx,ny,nz],
                            theta,axis_of_rot)
            for c in morph_obj.morph.children_of(n):        
                cx,cy,cz = morph_obj.shift_origin(np.array([c['x'],c['y'],c['z']]))
                [cx_rot,cy_rot,cz_rot] = morph_obj.rotate3D_point([cx,cy,cz],
                            theta,axis_of_rot)
                
                all_x.extend((nx_rot, cx_rot))
                all_y.extend((ny_rot, cy_rot))
                all_z.extend((nz_rot, cz_rot))
        
        return all_x,all_y,all_z
    
    
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


def convert_mpl_color_to_rgb_tuple(named_mpl_color):
    color = mpl.colors.to_rgb(named_mpl_color)
    color=tuple(color_*255 for color_ in color)
    return color