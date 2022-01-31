from PyQt5 import QtWidgets, QtCore
import sys
import FortranDriver
sys.path.append("GUI-resources")
from ui_mainwindow import Ui_MainWindow


class MainWindow(QtWidgets.QMainWindow):

    # Custom signal to communicate with real-time plot widgets
    # This needs to go here and not in the constructor. See:
    # https://stackoverflow.com/questions/2970312/pyqt4-qtcore-pyqtsignal-object-has-no-attribute-connect
    # for an explanation
    timestep_update = QtCore.pyqtSignal(int)

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.md = FortranDriver.MDInterface()
        # Timestep
        self.tau = 0

        # Keep track of the open plotting windows
        self.window_list = []

        # Register the signal handlers
        self.register_handlers()

        # Start with the pause button disabled until we start the simulation
        self.ui.pause_button.setEnabled(False)

        # Print the input values to the QTextBrowser widget
        self.param_str = self.md.format_params()
        self.ui.input_textbrowser.setPlainText(self.param_str)
        self.initialise_simulation()

    def register_handlers(self):
        # Set up all signal handlers not already setup in Qt Designer
        self.ui.start_button.clicked.connect(self.start_sim)
        self.ui.pause_button.clicked.connect(self.pause_sim)
        self.ui.restart_button.clicked.connect(self.restart_sim)

        # Add signals so spin-boxes change the underlying simulation parameters
        self.ui.lj_eps_spinbox.editingFinished.connect(self.update_parameters)
        self.ui.field_spinbox.editingFinished.connect(self.update_parameters)
        self.ui.temp_spinbox.editingFinished.connect(self.update_parameters)
        self.ui.density_spinbox.editingFinished.connect(self.update_parameters)

        # Checkbox to toggle NEMD field
        self.ui.nemd_checkbox.stateChanged.connect(self.toggle_ne_field)

        # Now initialise (but don't start) a timer to update the plot
        self.sim_timer = QtCore.QTimer()
        self.sim_timer.setInterval(8)

        # Now do a longer timer for updating the UI elements (but not the plot)
        self.gui_timer = QtCore.QTimer()
        self.gui_timer.setInterval(60)
        # Connect the timer's "timeout" (finished) event to our update function
        self.sim_timer.timeout.connect(self.update_plot_data)
        self.gui_timer.timeout.connect(self.update_GUI_elements)

        # Now connect the file editing dialog to our open and close menu buttons
        self.ui.open_input_file.triggered.connect(self.open_input_file)
        self.ui.save_input_file.triggered.connect(self.save_to_file)
        self.ui.actionQuit.triggered.connect(self.clean_exit)
        finish = QtWidgets.QAction("Quit", self)
        finish.triggered.connect(self.closeEvent)

    def initialise_simulation(self):

        # Fix the X and Y ranges so they don't constantly shift throughout the simulation
        self.ui.plot_window.clear()
        bounds = self.md.box_bounds
        self.ui.plot_window.setXRange(0, bounds[0])
        self.ui.plot_window.setYRange(0, bounds[1])
        self.ui.plot_window.setBackground('w')
        self.pos_data = self.ui.plot_window.plot(self.md.x, self.md.y, pen=None, symbol='o')

        # Now do the g(2) radial-distribution function
        compute_rdf = self.md.compute_rdf()
        r = compute_rdf['r']
        rdf = compute_rdf['rdf']
        self.ui.rdf_window.setXRange(0, max(r)+1)
        self.ui.rdf_window.setYRange(0, max(rdf)+1)
        self.ui.rdf_window.setBackground('w')
        self.rdf_data = self.ui.rdf_window.plot(r[1:], rdf[1:], color='k')

        # Initialise the N, V and T labels
        self.ui.npart_label.setText(f"N particles = {self.md.npart}")
        self.ui.volume_label.setText(f"Vol = {self.md.vol:.2f}")
        self.ui.temperature_label.setText(f"Temperature = {self.md.temp:.2f}")

        # Finally, set the simulation controls to the correct value
        self.ui.lj_eps_spinbox.setValue(self.md.eps)
        self.ui.field_spinbox.setValue(self.md.fieldstrength)
        self.ui.temp_spinbox.setValue(self.md.temp)
        self.ui.density_spinbox.setValue(self.md.reduced_density)

    def update_parameters(self):
        # Get the widget which sent this signal, as well as its new value
        sender = self.sender()
        value = sender.value()

        # Now change the appropriate simulation parameter
        if sender == self.ui.lj_eps_spinbox:
            self.md.kf = value

        elif sender == self.ui.field_spinbox:
            self.md.fieldstrength = value

        # These spinboxes control initial parameters, and require the simulation to be restarted after
        # changing
        elif sender == self.ui.temp_spinbox:
            self.md.temp = value

        elif sender == self.ui.density_spinbox:
            self.md.reduced_density = value
            # The density is controlled by changing the box size and keeping NPART constant, so
            # we need to adjust the plot window's range
            bounds = self.md.box_bounds
            self.ui.plot_window.setXRange(0, bounds[0])
            self.ui.plot_window.setYRange(0, bounds[1])
            self.ui.volume_label.setText(f"Vol = {self.md.vol:.2f}")

        else:
            print("Unknown sender")
            pass
        # Finally, update the input values in the QTextBrowser widget
        self.param_str = self.md.format_params()
        self.ui.input_textbrowser.setPlainText(self.param_str)
        self.ui.input_textbrowser.repaint()

    def toggle_ne_field(self, state):
        # Toggles the nonequilibrium field (on or off) based on the status of nemd_checkbox
        if state == QtCore.Qt.Checked:
            self.md.toggle_nemd()
        else:
            self.md.toggle_nemd()

        # Finally, update the input values in the QTextBrowser widget
        self.param_str = self.md.format_params()
        self.ui.input_textbrowser.setPlainText(self.param_str)
        self.ui.input_textbrowser.repaint()

    def update_plot_data(self):
        # First, run an MD timestep
        self.md.run(1)

        # Only plot the fluid particles for now
        self.pos_data.setData(self.md.x, self.md.y)  # Update the data.

        # Now do the g(2) radial-distribution function
        compute_rdf = self.md.compute_rdf()
        r = compute_rdf['r']
        rdf = compute_rdf['rdf']
        self.ui.rdf_window.setXRange(0, max(r)+1)
        self.ui.rdf_window.setYRange(0, max(rdf)+1)
        self.rdf_data.setData(r, rdf)

    def update_GUI_elements(self):
        # Update the temperature and volume labels
        temp = self.md.temp
        vol = self.md.vol
        self.ui.temperature_label.setText(f"Temperature = {temp:.2f}")
        self.ui.volume_label.setText(f"Volume = {vol:.2f}")
        self.ui.npart_label.setText(f"N particles = {self.md.npart}")

    def start_sim(self):
        """ Start the simulation.

            The timer has already been initialised and linked to the update function, so we only need to
            start the timer here."""
        self.sim_timer.start()
        self.gui_timer.start()

        self.ui.start_button.setEnabled(False)
        self.ui.pause_button.setEnabled(True)

        # Also want to disable the density spinbox, since it makes no sense to change it
        # while the simulation is running
        self.ui.density_spinbox.setEnabled(False)

        # Finally, run an MD timestep
        self.md.run(1)

    # Plotting routines
    def pause_sim(self):
        """ Pause the simulation.

            This is simplest to achieve by simply stopping the timer temporarily, so the plot stops
            updating. It will start back up again when the timer is restarted."""
        if self.sim_timer.isActive():
            self.sim_timer.stop()
            self.gui_timer.stop()

        self.ui.start_button.setEnabled(True)
        self.ui.pause_button.setEnabled(False)

    def restart_sim(self):
        """ Restart the simulation by stopping the timer and reinitialising parameters."""
        if self.sim_timer.isActive():
            self.sim_timer.stop()
            self.gui_timer.stop()
        self.tau = 0

        # Re-enable buttons which can't be changed while the simulation is running
        self.ui.density_spinbox.setEnabled(True)
        self.ui.start_button.setEnabled(True)

        self.ui.plot_window.clear()
        self.ui.rdf_window.clear()

        self.md.setup()
        self.initialise_simulation()

    # I/O Control routines
    def open_input_file(self):
        self.sim_timer.stop()
        input_file, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open File", "", "Input Files (*.in)")
        if input_file:
            self.md.read_from_file(input_file)
            # Update the input values in the QTextBrowser widget
            self.param_str = self.md.format_params()
            self.ui.input_textbrowser.setPlainText(self.param_str)
            self.ui.input_textbrowser.repaint()
            self.restart_sim()

    def save_to_file(self):
        self.sim_timer.stop()
        self.gui_timer.stop()
        output_file, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save File", "", "Input Files (*.in)")
        if output_file:
            out_string = self.md.format_params()
            with open(output_file, 'w') as ofp:
                ofp.write(out_string)

    # Clean exit
    def clean_exit(self):
        # Close all open windows when the main window is closed
        for window in self.window_list:
            window.close()

        self.close()

    def closeEvent(self, event):
        # Close all open windows when the main window is closed
        for window in self.window_list:
            window.close()

        event.accept()
