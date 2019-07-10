import sys
from PyQt5 import QtCore, QtGui, uic
import pyqtgraph as pg
import numpy as np
import glob
import os
import bluepyopt as bpopt
import bluepyopt.ephys as ephys
from ateamopt.bpopt_evaluator import Bpopt_Evaluator
from ateamopt.utils import utility


pg.setConfigOption('background', 'w')

 
qtCreatorFile = "Model_viz.ui" # Enter file here.
 
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)
 
class NrnModelViz(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        self.sim_components = os.path.join(os.getcwd(),'sim_components')
        self.morph_path = glob.glob(self.sim_components+'/*.swc')
        if len(self.morph_path) > 1:
            raise Exception
        else:
            self.morph_path = self.morph_path[0]
        self.bpopt_config = os.path.join(os.getcwd(),'bluepyopt_config')
        self.bmtk_config = os.path.join(os.getcwd(),'bmtk_config')
        self.RunSimulation.clicked.connect(self.run_bpopt_sim)
    
    @staticmethod
    def prepare_stim_protocol(stim_amp,delay=270,duration=1000,
                              total_duration=2270):
        
        protocol_name = 'DC_stim'
        
        
        
        soma_loc = ephys.locations.NrnSeclistCompLocation(
            name='soma',
            seclist_name='somatic',
            sec_index=0,
            comp_x=0.5)
        somav_recording = ephys.recordings.CompRecording(
            name='%s.soma.v' %
            protocol_name,
            location=soma_loc,
            variable='v')

        recordings = [somav_recording]
        
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
        stim_protocol,protocol_name = self.prepare_stim_protocol(stim_amp)
        param_path = os.path.join(self.bpopt_config,'parameters.json')
        mech_path= os.path.join(self.bpopt_config,'mechanism.json')
        release_param_path = os.path.join(self.bpopt_config,'release_param.json')
        release_params = utility.load_json(release_param_path)
        
        morphology = ephys.morphologies.\
            NrnFileMorphology(self.morph_path, stub_axon=True)
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
        
        self.simview.clear()
#        a=np.random.randn(10)
#        b=np.random.randn(10)
        self.simview.plot(time_response,voltage_response,
                          pen=pg.mkPen('k', width=2))
        
 
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = NrnModelViz()
    window.show()
    sys.exit(app.exec_())