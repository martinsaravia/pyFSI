class fsiModel:
    def __init__(self, control, beam, flow):
        # Store the beam and flow objects
        self._beam = beam
        self._flow = flow
        # Store the control dict
        self._control = control

    # Getters
    def control(self):
        return self._control

    def solid(self):
        return self._beam

    def flow(self):
        return self._flow
