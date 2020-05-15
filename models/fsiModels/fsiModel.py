class fsiModel:
    def __init__(self, execution, control, beam, flow):
        # Store objects
        self._beam = beam
        self._flow = flow
        self._execution = execution
        self._control = control

    # Getters
    def control(self):
        return self._control

    def execution(self):
        return self._execution

    def solid(self):
        return self._beam

    def flow(self):
        return self._flow
