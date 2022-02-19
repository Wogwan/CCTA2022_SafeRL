import cbf_sosp

class LEARNER():
    def __init__(self,env,eng):
        self.env = env
        self.eng = eng
        # Build barrier function

    def current_cbf(self):
        coeffs = cbf_sosp.build_barrier(self, self.eng)
        return coeffs
