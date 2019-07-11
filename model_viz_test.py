import sys
from PyQt5 import QtCore, QtGui, uic
from PyQt5.QtWidgets import QFileDialog
import pyqtgraph as pg
import numpy as np
import glob
import os
import bluepyopt as bpopt
import bluepyopt.ephys as ephys
from ateamopt.bpopt_evaluator import Bpopt_Evaluator
from ateamopt.utils import utility
import utility_functions as uf


pg.setConfigOption('background', 'w')


qtCreatorFile = "Model_viz.ui" # Enter file here.

Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

class NrnModelViz(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        self.sim_components = os.path.join(os.getcwd(),'sim_components')
#        self.morph_path = glob.glob(self.sim_components+'/*.swc')
#        if len(self.morph_path) > 1:
#            raise Exception
#        else:
#            self.morph_path = self.morph_path[0]
        self.bpopt_config = os.path.join(os.getcwd(),'bluepyopt_config')
        self.bmtk_config = os.path.join(os.getcwd(),'bmtk_config')
        self.model = None
        self.RunSimulation.clicked.connect(self.run_bpopt_sim)
        self.closeWindow.clicked.connect(self.closeApp)

    def import_model(self):
        options = QFileDialog.Options()
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*)", options=options)
        param_dict = utility.load_json(fileName)
        aibs_param_dict = uf.param_dict_to_AIBS_param(param_dict)
        aibs_param_filename = os.path.join(self.bmtk_config,
                   'biophysical_neuron_templates/param_aibs_aa_model.json')
        utility.create_filepath(aibs_param_filename)
        utility.save_json(aibs_param_filename,aibs_param_dict)

    def import_morphology(self):
        options = QFileDialog.Options()
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*)", options=options)
        self.morph_path = fileName

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
        delay,duration,total_duration=270,1000,2270
        stim_protocol,protocol_name = self.prepare_stim_protocol(stim_amp)
        
        param_path = os.path.join(self.bpopt_config,'parameters.json')
        mech_path= os.path.join(self.bpopt_config,'mechanism.json')
        release_param_path = os.path.join(self.bpopt_config,'release_param.json')
        release_params = utility.load_json(release_param_path)
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

        self.simview.clear()

        self.simview.plot(time_response,voltage_response,
                          pen=pg.mkPen('k', width=4))
        self.simview.setRange(xRange=[delay-50,delay+duration+50])

    def closeApp(self):
        sys.exit()


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = NrnModelViz()
    window.show()
    sys.exit(app.exec_())
