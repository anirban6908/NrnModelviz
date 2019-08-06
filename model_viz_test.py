import sys
from PyQt5 import QtCore, QtGui, uic
from PyQt5.QtWidgets import QFileDialog
import pyqtgraph as pg
import numpy as np
import os
import math
import os
import bluepyopt.ephys as ephys
from ateamopt.bpopt_evaluator import Bpopt_Evaluator
from ateamopt.utils import utility
import utility_functions as uf
import ateam.sim.setup as sim
from ateamopt.morph_handler import Morph_handler
import shutil
import h5py
from bmtk.builder.networks import NetworkBuilder
import efel
from ateamopt.analysis.analysis_module import get_spike_shape
import pyqtgraph.opengl as gl


#pg.setConfigOption('background', 'w')


qtCreatorFile = "Model_viz.ui" # Enter file here.

Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

class NrnModelViz(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
#        self.sim_components = os.path.join(os.getcwd(),'sim_components')
        self.bpopt_config = os.path.join(os.getcwd(),'bluepyopt_config')
        self.bpopt_sim = os.path.join(os.getcwd(),'bluepyopt_simulation')
        self.bmtk_config = os.path.join(os.getcwd(),'bmtk_config')
        self.bmtk_sim = os.path.join(os.getcwd(),'bmtk_simulation')
        self.model = None
        self.RunSimulation.clicked.connect(self.run_bpopt_sim)
        self.closeWindow.clicked.connect(self.closeApp)
        self.plotMultipatchRec.clicked.connect(self.plot_comp_recording)
        self.plotEAP.clicked.connect(self.run_bmtk_sim)
        self.eap_STA.clicked.connect(self.plot_eap_STA)
        self.plotAPshape.clicked.connect(self.plot_ap_shape)
        self.vizMorph.clicked.connect(self.viz_morph)

    def import_model(self):
        options = QFileDialog.Options()
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*)", options=options)
        param_unknown = utility.load_json(fileName)
        bpopt_mech_filename = os.path.join(self.bpopt_config,
                           'mech_bpopt_aa_model.json')
        bpopt_param_filename = os.path.join(self.bpopt_config,
                           'param_bpopt_aa_model.json')
        if isinstance(param_unknown,dict):
            if 'genome' in list(param_unknown.keys()): # Parameter in AIBS format
                bpopt_param,bpopt_mech = uf.AIBS_param_to_bpopt_param(\
                                                      param_unknown)
                utility.create_filepath(bpopt_param_filename)
                utility.save_json(bpopt_param_filename,bpopt_param)
                utility.save_json(bpopt_mech_filename,bpopt_mech)
            else:
                aibs_param_dict = uf.param_dict_to_AIBS_param(param_unknown,
                                                  no_apical=self.no_apical)
                aibs_param_filename = os.path.join(self.bmtk_config,
                           'biophysical_neuron_templates/param_aibs_aa_model.json')

                utility.create_filepath(aibs_param_filename)
                utility.save_json(aibs_param_filename,aibs_param_dict)

                bpopt_param,bpopt_mech = uf.AIBS_param_to_bpopt_param(\
                                                      aibs_param_dict)
                utility.create_filepath(bpopt_param_filename)
                utility.save_json(bpopt_param_filename,bpopt_param)
                utility.save_json(bpopt_mech_filename,bpopt_mech)

    def import_morphology(self):
        options = QFileDialog.Options()
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*)", options=options)
        self.morph_path = fileName
        self.no_apical = utility.check_swc_for_apical(self.morph_path)

    def prepare_stim_protocol(self,stim_amp,delay=270,duration=1000,
                              total_duration=2270):

        protocol_name = 'DC_stim'
        axon_sec_len = 30.0
        soma_loc = ephys.locations.NrnSeclistCompLocation(
            name='soma',
            seclist_name='somatic',
            sec_index=0,
            comp_x=0.5)

        dend_elec_pos = float(self.intra_electrode_dendrite.value())
        axon_elec_pos = float(self.intra_electrode_axon.value())
        somav_recording = ephys.recordings.CompRecording(
            name='%s.soma.v' %
            protocol_name,
            location=soma_loc,
            variable='v')

        recordings = [somav_recording]
        dend_type = 'basal' if self.no_apical else 'apical'
        extra_recording_def = {dend_type:dend_elec_pos,
                               'axonal':axon_elec_pos}

        for rec_key,rec_val in extra_recording_def.items():
            if rec_val != 0:

                if rec_key == 'axonal':
                    sec_index = math.floor(rec_val/axon_sec_len)
                    location = ephys.locations.NrnSeclistCompLocation(
                        name='loc_%s'%rec_key,
                        sec_index=sec_index,
                        comp_x=0.5,
                        seclist_name=rec_key)
                else:
                    location = ephys.locations.NrnSomaDistanceCompLocation(
                    name='loc_%s'%rec_key,
                    soma_distance=rec_val,
                    seclist_name=rec_key)
                var = 'v'
                recording = ephys.recordings.CompRecording(
                name='%s.%s.%s' % (protocol_name, rec_key, var),
                location=location,
                variable=var)
                recordings.append(recording)

        stimuli = []
        stimuli.append(ephys.stimuli.NrnSquarePulse(
            step_amplitude=stim_amp,
            step_delay=float(delay),
            step_duration=float(duration),
            location=soma_loc,
            total_duration=float(total_duration)))

        stim_protocol =  ephys.protocols.SweepProtocol(
                    protocol_name,
                    stimuli,
                    recordings)
        return stim_protocol,protocol_name

    def run_bpopt_sim(self):
        stim_amp = float(self.stimAmp.value())*1e-3 # in nA
        delay,duration,total_duration=270.0,1000.0,2270.0
        stim_protocol,protocol_name = self.prepare_stim_protocol(stim_amp)

        param_path = os.path.join(self.bpopt_config,'param_bpopt_aa_model.json')
        mech_path= os.path.join(self.bpopt_config,'mech_bpopt_aa_model.json')

        release_params = {}
        axon_type = str(self.axonType.currentText())
        if axon_type == 'Stub axon':
            morphology = ephys.morphologies.\
                NrnFileMorphology(self.morph_path, stub_axon=True)
        else:
            morphology = ephys.morphologies.\
                NrnFileMorphology(self.morph_path, do_replace_axon=True)
        sim = ephys.simulators.NrnSimulator()
        eval_handler = Bpopt_Evaluator(protocol_path=None,
                           feature_path=None,
                           morph_path=None, param_path=param_path,
                           mech_path=mech_path)
        mechanisms = eval_handler.define_mechanisms()
        parameters = eval_handler.define_parameters()

        model_aa = ephys.models.CellModel('aa_model',
                       morph=morphology, mechs=mechanisms,
                       params=parameters)

        responses = stim_protocol.run(model_aa,release_params,sim)
        time_response = responses['{}.soma.v'.format(protocol_name)]['time']
        voltage_response =responses['{}.soma.v'.format(protocol_name)]['voltage']

        resp_filename = os.path.join(self.bpopt_sim,'model_resp.pkl')
        utility.create_filepath(resp_filename)
        utility.save_pickle(resp_filename,responses)

        self.simview.clear()

        self.simview.plot(time_response,voltage_response,
                          pen=pg.mkPen('b', width=3))
        self.simview.setRange(xRange=[delay-50,delay+duration+50])


    def plot_comp_recording(self):
        self.PatchView.clear()
        delay,duration,total_duration=270.0,1000.0,2270.0
        resp_filename = os.path.join(self.bpopt_sim,'model_resp.pkl')
        responses = utility.load_pickle(resp_filename)
        color_dict = {'soma':'b','apical': 'r','axonal':'g'}
        self.PatchView.addLegend()
        for resp_key,resp_val in responses.items():
            comp_name = resp_key.split('.')[1]
            if comp_name != 'soma':
                t_ = resp_val['time']
                v_ = resp_val['voltage']
                self.PatchView.plot(t_,v_,
                         pen=pg.mkPen(color_dict[comp_name], width=3),
                         name = comp_name)
        midpoint = delay+duration/2
        self.PatchView.setRange(xRange=[midpoint-50,midpoint+50])
#        self.PatchView.setRange(xRange=[delay-50,delay+duration+50])

    def run_bmtk_sim(self):
        stim_amp = float(self.stimAmp.value())*1e-3 # in nA
        delay,duration,total_duration=270.0,1000.0,2270.0
        elec_pos_x = float(self.electrode_x.value())
        elec_pos_y = float(self.electrode_y.value())
        elec_pos_z = float(self.electrode_z.value())
        elec_pos = [elec_pos_x,elec_pos_y,elec_pos_z]
        bmtk_config_template = os.path.join(self.bmtk_config,
                                            "default_config.json")

        sim_folder =os.path.basename(self.bmtk_sim)
        sm = sim.SimManager.from_template(config_template=bmtk_config_template, overwrite=True,
                                          sim_folder=sim_folder)
        network_name = 'snm'
        cell_name = 'aa_model'

        net = NetworkBuilder(network_name)

        bmtk_morph_path = os.path.join(self.bmtk_config,'morphologies/{}.swc'.\
                               format(cell_name))
        utility.create_filepath(bmtk_morph_path)
        shutil.copyfile(self.morph_path,bmtk_morph_path)
        morph_handler = Morph_handler(self.morph_path)
        morph_data,morph_apical,morph_axon = morph_handler.\
                                    get_morph_coords()
        theta_x,theta_y,theta_z = morph_handler.calc_euler_angle\
                                        (morph_data,morph_apical)

        axon_type = str(self.axonType.currentText())
        if axon_type == 'Stub axon':
            axon_processing = 'aibs_allactive_ani_directed'
        else:
            axon_processing = 'aibs_allactive_bpopt_axon_directed'

        node_props = {'cell_name' : cell_name,
                      'model_type' : 'biophysical',
                      'model_template' : 'ctdb:Biophys1.hoc',
                      'model_processing' : axon_processing,
                      'dynamics_params' : 'param_aibs_%s.json'%cell_name,
                      'morphology': '%s.swc'%cell_name,
                      'rotation_angle_xaxis':[theta_x],
                      'rotation_angle_yaxis':[theta_y],
                      'rotation_angle_zaxis':[-theta_z]
                      }

        net.add_nodes(**node_props)

        protocol_name = 'DC_stim'
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
        dt = 0.01
        sm.add_network(net)
        net.build()
        sm.save_network_files(use_abs_paths=True)
        sm.add_ecp_report(locs=[elec_pos])
        sm.sim_time = total_duration
        sm.sim_timestep = dt
        sm.run_bionet()

        ecp_report_file = os.path.join('{}/output'.format(sim_folder),
                   sm.config['reports']['ecp_report']['file_name'])


        with h5py.File(ecp_report_file,'r') as ecp_data:
            ecp_rec=ecp_data['data']
            ecp_rec=1e3 * np.array(ecp_rec)  # Change the unit from mV to uV
        Ve = np.transpose(ecp_rec)   #ecp_data: channels*times
        t = np.arange(Ve.shape[1])*dt
        self.EAPview.clear()
        self.EAPview.plot(t,np.squeeze(Ve),
                          pen=pg.mkPen('y', width=2))
        self.EAPview.setRange(xRange=[delay-50,delay+duration+50])

    def plot_eap_STA(self):
        sim_folder =os.path.basename(self.bmtk_sim)
        config_path = os.path.join(self.bmtk_sim,'config.json')
        sm = sim.SimManager(config_path= config_path,sim_folder=sim_folder)
        spike_report_file = os.path.join('{}/output'.format(sim_folder),
                           sm.config['output']['spikes_file'])
        ecp_report_file = os.path.join('{}/output'.format(sim_folder),
                   sm.config['reports']['ecp_report']['file_name'])
        win=[-1,3]

        with h5py.File(spike_report_file,'r') as spike_data:
            spike_times = np.array(spike_data['spikes/timestamps'],
                                    dtype=np.float)

        with h5py.File(ecp_report_file,'r') as ecp_data:
            ecp_rec=ecp_data['data']
            ecp_rec=1e3 * np.array(ecp_rec)  # Change the unit from mV to uV
        Ve = np.transpose(ecp_rec)   #ecp_data: channels*times

        [t_STA,oneCol_VeSTA,oneCol_VeSTC] = uf.cal_STA_STC(Ve,spike_times,
                    win,sm.sim_timestep)

        oneCol_VeSTA = np.squeeze(oneCol_VeSTA)
        oneCol_VeSTC = np.squeeze(oneCol_VeSTC)
        oneCol_VeSTA_plus = oneCol_VeSTA + oneCol_VeSTC
        oneCol_VeSTA_minus = oneCol_VeSTA - oneCol_VeSTC

        self.EAPview.clear()
        self.EAPview.plot(t_STA,oneCol_VeSTA,
                          pen=pg.mkPen('y', width=2))
        fill_brush = pg.mkBrush(200, 0, 0,150)
        sta_low = pg.PlotCurveItem(t_STA, oneCol_VeSTA_minus,
                                   brush=fill_brush)
        sta_high = pg.PlotCurveItem(t_STA, oneCol_VeSTA_plus,
                                   brush=fill_brush)
        fill = pg.FillBetweenItem(sta_low,
                sta_high,brush=fill_brush)
        self.EAPview.addItem(fill)
        self.EAPview.setRange(xRange=[t_STA[0]-.1,t_STA[-1]+.1])

    def plot_ap_shape(self):
        delay,duration,total_duration=270.0,1000.0,2270.0
        protocol_name = 'DC_stim'
        resp_filename = os.path.join(self.bpopt_sim,'model_resp.pkl')
        responses = utility.load_pickle(resp_filename)
        time = responses['{}.soma.v'.format(protocol_name)]['time']
        voltage =responses['{}.soma.v'.format(protocol_name)]['voltage']

        # Prepare sweep for eFEL
        sweep = {}
        sweep['T'] = time
        sweep['V'] = voltage
        sweep['stim_start'] = [delay]
        sweep['stim_end'] = [delay+duration]
        sweeps = [sweep]

        # Extract experimental spike times
        feature_results = efel.getFeatureValues(sweeps, ['peak_time'])
        spike_times = feature_results[0]['peak_time']

        prefix_pad,posfix_pad,res = 2,5,0.05
        AP_shape_time = np.arange(-prefix_pad,posfix_pad, res)
        AP_shape_voltage = np.zeros(AP_shape_time.size)
        AP_shape_voltage = get_spike_shape(time,voltage,
                        spike_times,AP_shape_time,
                        AP_shape_voltage)
        AP_shape_voltage /= len(spike_times)

        self.APshapeview.clear()

        self.APshapeview.plot(AP_shape_time,AP_shape_voltage,
                          pen=pg.mkPen('b', width=3))
        self.APshapeview.setRange(xRange=[AP_shape_time[0]-.1,
                                      AP_shape_time[-1]+.1])

    def viz_morph(self):
        morph_handler = Morph_handler(self.morph_path)
        morph_data,morph_apical,morph_axon = morph_handler.\
                            get_morph_coords()
        theta,axis_of_rot = morph_handler.calc_rotation_angle\
                                        (morph_data,morph_apical)
        morph_x,morph_y,morph_z = uf.get_morph_points(morph_handler,theta,axis_of_rot)
        morph_pos = np.array([morph_x,morph_y,morph_z]).T
        
        morph_z_center = (np.max(morph_z) + np.min(morph_z))/2 
        morph_pos[:,2] -= morph_z_center
#        self.morphview.clear()
        for item_ in self.morphwidget.items:
            self.morphwidget.items.remove(item_)
        
        size_vec = 1.5*np.ones(morph_pos.shape[0])
#        color_arr = np.array([255,105,180,255])/255.0
        color_arr = (1, 0,0,.5)
        color_mat = np.tile(color_arr,(morph_pos.shape[0],1))
        sc = gl.GLScatterPlotItem(pos=morph_pos,
                    color = color_arr,
                    size=size_vec,
                    pxMode=False)
        soma_pos = np.reshape(np.array([0,0,-morph_z_center]),(1,3))
        sc_soma = gl.GLScatterPlotItem(pos=soma_pos,
                    color = (1, 0,0,1),
                    size=15,
                    pxMode=False)
#        sc= pg.ScatterPlotItem(x=morph_x,y=morph_z,pen=pg.mkPen(None),
#                        brush=pg.mkBrush(255,105,180,200),size=1.5)
#        sc.addPoints(x=[0],y=[0],pen=pg.mkPen(None),
#                        brush=pg.mkBrush(255,105,180,255),size=15)
#    
#        self.morphview.addItem(sc)
        
        self.morphwidget.addItem(sc)
        self.morphwidget.addItem(sc_soma)
        self.morphwidget.setCameraPosition(distance=250,
#                                           azimuth = 45,
                                           elevation=30)

    def closeApp(self):
        sys.exit()


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = NrnModelViz()
    window.show()
    sys.exit(app.exec_())
