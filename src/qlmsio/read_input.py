import logging
import numpy as np
import re
from requests.structures import CaseInsensitiveDict

logger = logging.getLogger("qlms").getChild("read_input")


class QMLSInput():
    def __init__(self, file_name_list, solver_type="UHF"):
        self.param = CaseInsensitiveDict()
        self.file_names = self._read_file_names(file_name_list)
        # Get ModPara file
        self.param["modpara"] = self._read_para("modpara")
        self.param["calcmod"] = self._read_para("calcmod", start_line=0)
        for k, v in self.param["calcmod"].items():
            self.param["calcmod"][k] = int(v[0])
        self.param["locspin"] = self._read_para("locspin")
        for k, v in self.param["locspin"].items():
            self.param["locspin"][k] = int(v[0])
        self.ham_param = CaseInsensitiveDict()
        self.ham_param["Transfer"] = self._read_ham("trans", value_type="complex")
        self.ham_param["CoulombInter"] = self._read_ham("coulombinter")
        self.ham_param["CoulombIntra"] = self._read_ham("coulombintra")
        self.ham_param["PairHop"] = self._read_ham("parihop")
        self.ham_param["Hund"] = self._read_ham("hund")
        self.ham_param["Exchange"] = self._read_ham("exchange")
        self.ham_param["Interall"] = self._read_ham("interall", value_type="complex")
        self.ham_param["Initial"] = self._read_ham("initial", value_type="complex")
        # TODO Check Pair(Hermite or not)
        # TODO Check site is larger than defined lattics size
        # TODO Add validation function (ex.:Check site is smaller than defined lattics size)
        self.output = CaseInsensitiveDict()
        self.output["OneBodyG"] = self._read_green("onebodyg")
        self.output["TwoBodyG"] = self._read_green("twobodyg")

    def get_param(self, key):
        if key == "param":
            return self.param
        elif key == "ham":
            return self.ham_param
        elif key == "output":
            return self.output
        else:
            # Add error message
            logger.error("Get_param: key must be param or ham or output.")
            return None

    def _read_file_names(self, file_name_list):
        file_names = CaseInsensitiveDict()
        with open(file_name_list, "r") as f:
            lines = f.readlines()
            for line in lines:
                words = line.split()
                file_names[words[0]] = words[1]
        # TODO Check essential files
        return file_names

    def _read_para(self, file_key, start_line=5):
        file_name = self.file_names[file_key]
        value = CaseInsensitiveDict()
        with open(file_name, "r") as f:
            lines = f.readlines()
            match_pattern = r'^\s*(#.*|)$'
            for line in lines[start_line:]:
                match = re.match(match_pattern, line)
                if not match:
                    words = line.split()
                    value[words[0]] = words[1:]
        return value

    def _read_ham(self, file_key, value_type="real"):
        ham_info = {}
        if file_key in self.file_names:
            file_name = self.file_names[file_key]
            if value_type == "real":
                with open(file_name, "r") as f:
                    lines = f.readlines()
                    for line in lines[5:]:
                        words = line.split()
                        list = tuple([int(i) for i in words[:-1]])
                        value = float(words[-1])
                        if list in ham_info:
                            ham_info[list] += value
                        else:
                            ham_info[list] = value
            if value_type == "complex":
                with open(file_name, "r") as f:
                    lines = f.readlines()
                    for line in lines[5:]:
                        words = line.split()
                        list = tuple([int(i) for i in words[:-2]])
                        value = float(words[-2]) + 1j * float(words[-1])
                        if list in ham_info:
                            ham_info[list] += value
                        else:
                            ham_info[list] = value
        else:
            return None
        return ham_info

    def _read_green(self, file_key):
        list = []
        if file_key in self.file_names:
            file_name = self.file_names[file_key]
            with open(file_name, "r") as f:
                lines = f.readlines()
                for line in lines[5:]:
                    list.append(line.split())
        else:
            return None
        return list


if __name__ == '__main__':
    qml_input = QMLSInput("namelist.def")
