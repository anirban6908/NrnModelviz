import sys
from PyQt5 import QtCore, QtGui, uic
import pyqtgraph as pg
import numpy as np
import glob
import os
import bluepyopt as bpopt
import bluepyopt.ephys as ephys
from ateamopt.bpopt_evaluator import Bpopt_Evaluator


pg.setConfigOption('background', 'w')

 
qtCreatorFile = "Model_viz.ui" # Enter file here.
 
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)
 
class NrnModelViz(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        self.sim_components = os.path.join(os.getcwd(),'sim_components')
        self.morph_path = glob.glob(self.sim_components+'*.swc')
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
            step_delay=delay,
            step_duration=duration,
            location=soma_loc,
            total_duration=total_duration))
        
        
        stim_protocol = {protocol_name : ephys.protocols.SweepProtocol(
                    protocol_name,
                    stimuli,
                    recordings)}
        return stim_protocol
    
    def run_bpopt_sim(self):
        stim_amp = float(self.stimAmp.toPlainText()) # in pA
        stim_protocol = self.prepare_stim_protocol(stim_amp)
        
        morphology = ephys.morphologies.\
            NrnFileMorphology(self.morph_path, do_replace_axon=True)
        sim = ephys.simulators.NrnSimulator()
        eval_handler = Bpopt_Evaluator(
                   morph_path, param_path,
                   mech_path,timed_evaluation = False)
        
        self.simview.clear()
        a=np.random.randn(10)
        b=np.random.randn(10)
        self.simview.plot(a,b,pen=pg.mkPen('k', width=2))
        
 
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = NrnModelViz()
    window.show()
    sys.exit(app.exec_())