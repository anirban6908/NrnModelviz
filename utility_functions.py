from ateamopt.utils import utility
from collections import defaultdict
import re
        
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
            
        if bool(re.search('Na',param_name,re.IGNORECASE)):
            ena_sectionlist.append(sect)
        elif  bool(re.search('K',param_name,re.IGNORECASE)):
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
