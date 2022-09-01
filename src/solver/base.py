import logging
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
            })

            for k, v in para_init.items():
                self.param_mod.setdefault(k, v)

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

        if self._check_info_mod(info_mode) != 0:
            exit(1)
        if self._check_param_mod(self.param_mod) != 0:
            exit(1)

        # canonicalize
        self.param_mod["EPS"] = pow(10, -self.param_mod["EPS"])

    def _check_info_mod(self, info_mode):
        exit_code = 0
        fix_list = {"mode": ["UHF", "UHFk"], "flag_fock": [True, False]}
        logger = logging.getLogger("qlms").getChild(self.name)
        for key, _fix_list in fix_list.items():
            if info_mode[key] not in _fix_list:
                logger.warning("mode.{} in mode section is incorrect: {}.".format(key, _fix_list))
                exit_code += 1
        return exit_code

    def _check_param_mod(self, param_mod):
        error_code = 0
        min_list = {"T": 0, "2Sz": -param_mod["Nsite"],
                      "Nsite": 1, "Ncond": 1, "IterationMax":0,
                    "Mix":0, "print_step": 1,
                    "EPS": 0}
        max_list = {"2Sz": param_mod["Nsite"],
                    # "Ncond": param_mod["Nsite"], # UHF and UFHk have different maximums...
                    "Mix":1}

        keyset = set(min_list.keys())
        keyset = keyset.union(set(max_list.keys()))

        logger = logging.getLogger("qlms").getChild(self.name)

        for key in keyset:
            if key not in keyset:
                logger.warning("mode.param.{} is not defined.".format(key))
                error_code += 1

        if param_mod["2Sz"] is not None:
            if ((param_mod["Ncond"] % 2 == 0 and param_mod["2Sz"] % 2 == 1) or
                (param_mod["Ncond"] % 2 == 1 and param_mod["2Sz"] % 2 == 0)):
                logger.error("mode.param.2Sz must be even(odd) when Ncond is even(odd).")
                error_code += 1

        for key, value in min_list.items():
            if key in param_mod and param_mod[key] < value:
                logger.error("value of mode.param.{} must be greater than {}.".format(key, value))
                error_code += 1

        for key, value in max_list.items():
            if key in param_mod and param_mod[key] > value:
                logger.error("value of mode.param.{} must be smaller than {}.".format(key, value))
                error_code += 1

        return error_code

    def solve(self, path_to_output):
        pass

    def get_results(self):
        return (self.physics, self.Green)

    def save_results(self, info_outputfile, green_info):
        pass
