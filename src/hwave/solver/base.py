import logging
from requests.structures import CaseInsensitiveDict

# from pprint import pprint

class solver_base():
    def __init__(self, param_ham, info_log, info_mode, param_mod=None):
        logger = logging.getLogger("qlms").getChild(self.name)

        self.param_mod = param_mod
        self.param_ham = param_ham
        self.info_log = info_log
        info_mode_param = info_mode.get("param", CaseInsensitiveDict({}))

        # constants
        self.threshold = 1.0e-12  # cutoff on elements of green function

        # default values
        para_init = CaseInsensitiveDict({
            "Nsite":  0,
            "Ncond":  0,
            "2Sz":    None,
            "Mix":    0.5,
            "EPS":    6,
            "IterationMax": 20000,
            "RndSeed": 1234,
            "T":      0.0,
            "ene_cutoff": 1.0e+2,
        })

        range_list = {
            "T":     [ 0, None ],
            # "2Sz":   [ -param_mod["Nsite"], param_mod["Nsite"] ],
            # "Nsite": [ 1, None ],
            "Ncond": [ 1, None ],
            "IterationMax": [ 0, None ],
            "Mix":   [ 0.0, 1.0 ],
            # "print_step": [ 1, None ],
            "EPS":   [ 0, None ],
        }

        if param_mod is not None:
            for k, v in para_init.items():
                self.param_mod.setdefault(k, v)

            # integer parameters
            for key in ["Nsite", "Ncond", "2Sz", "EPS", "IterationMax", "RndSeed"]:
                if self.param_mod[key] is not None and type(self.param_mod[key]) == type([]):
                    self.param_mod[key] = int(self.param_mod[key][0])
            # float parameters
            for key in ["Mix", "T", "ene_cutoff"]:
                if self.param_mod[key] is not None and type(self.param_mod[key]) == type([]):
                    self.param_mod[key] = float(self.param_mod[key][0])

            # overwrite by info_mode_param
            for key, value in info_mode_param.items():
                self.param_mod[key] = value
        else:
            self.param_mod = CaseInsensitiveDict(para_init)
            for key, value in info_mode_param.items():
                self.param_mod[key] = value

        if "2Sz" in self.param_mod and self.param_mod["2Sz"] == "free":
            self.param_mod["2Sz"] = None

        if self._check_info_mode(info_mode) != 0:
            logger.error("Parameter check failed for info_mode.")
            exit(1)
        if self._check_param_range(self.param_mod, range_list) != 0:
            logger.error("Parameter range check failed for param_mod.")
            exit(1)
        if self._check_param_mod(self.param_mod) != 0:
            logger.error("Parameter check failed for param_mod.")
            exit(1)

        # canonicalize
        self.param_mod["EPS"] = pow(10, -self.param_mod["EPS"])

        # debug
        # pprint(self.param_mod)

    def _check_info_mode(self, info_mode):
        logger = logging.getLogger("qlms").getChild(self.name)

        fix_list = {
            "mode": ["UHFr", "UHFk"],
            # "flag_fock": [True, False]
        }

        exit_code = 0
        for key, _fix_list in fix_list.items():
            if key not in info_mode:
                logger.warning("mode.{} is not defined.".format(key))
                exit_code += 1
            elif info_mode[key] not in _fix_list:
                logger.warning("mode.{} in mode section is incorrect: {}.".format(key, _fix_list))
                exit_code += 1
        return exit_code

    def _check_param_mod(self, param_mod):
        logger = logging.getLogger("qlms").getChild(self.name)

        error_code = 0
        if "2Sz" in param_mod and param_mod["2Sz"] is not None:
            ncond_eo = param_mod["Ncond"] % 2
            twosz_eo = param_mod["2Sz"] % 2
            if ((ncond_eo == 0 and twosz_eo == 1) or
                (ncond_eo == 1 and twosz_eo == 0)):
                logger.error("mode.param.2Sz must be even(odd) when Ncond is even(odd).")
                error_code += 1
        return error_code

    def _check_param_range(self, param_mod, range_list):
        logger = logging.getLogger("qlms").getChild(self.name)

        error_code = 0
        for key, [ vmin, vmax ] in range_list.items():
            if key not in param_mod:
                logger.warning("mode.param.{} is not defined.".format(key))
                error_code += 1
            elif param_mod[key] is not None:
                if vmin is not None and param_mod[key] < vmin:
                    logger.warning("mode.param.{} must be greater than {}.".format(key, vmin))
                    error_code += 1
                if vmax is not None and param_mod[key] > vmax:
                    logger.warning("mode.param.{} must be smaller than {}.".format(key, vmax))
                    error_code += 1
        return error_code

    def solve(self, path_to_output):
        pass

    def get_results(self):
        return (self.physics, self.Green)

    def save_results(self, info_outputfile, green_info):
        pass
