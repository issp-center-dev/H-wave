from requests.structures import CaseInsensitiveDict

class solver_base():
    def __init__(self, param_ham, info_log, info_mode, param_mod=None):
        self.param_mod = param_mod
        self.param_ham = param_ham
        self.info_log = info_log
        info_mode_param = info_mode.get("param", CaseInsensitiveDict({}))

        if param_mod is not None:
            # initial values
            para_init = CaseInsensitiveDict({
                "Nsite": 0,
                "Ne": 0,
                "Ncond": 0,
                "2Sz": None,
                "Mix": 0.5,
                "EPS": 6,
                "Print": 0,
                "IterationMax": 20000,
                "RndSeed": 1234,
                "CDataFileHead": "zvo",
                "CParaFileHead": "zqp",
            })

            for k, v in para_init.items():
                self.param_mod.setdefault(k, v)

            # set initial values
            for key in ["CDataFileHead", "CParaFileHead"]:
                if self.param_mod[key] is not None and type(self.param_mod[key]) == type([]):
                    self.param_mod[key] = str(self.param_mod[key][0])

            for key in ["nsite", "ne", "2Sz", "ncond", "eps", "IterationMax", "Print", "RndSeed"]:
                if self.param_mod[key] is not None and type(self.param_mod[key]) == type([]):
                    self.param_mod[key] = int(self.param_mod[key][0])

            #The type of mix is float.
            for key in ["mix"]:
                if self.param_mod[key] is not None and type(self.param_mod[key]) == type([]):
                    self.param_mod[key] = float(self.param_mod[key][0])

            # overwrite by info_mode_param
            for key, value in info_mode_param.items():
                self.param_mod[key] = value
        else:
            self.param_mod = CaseInsensitiveDict(info_mode_param)

        # canonicalize
        self.param_mod["EPS"] = pow(10, -self.param_mod["EPS"])

    def solve(self, path_to_output):
        pass

    def get_results(self):
        return (self.physics, self.Green)

    def save_results(self, info_outputfile, green_info):
        pass