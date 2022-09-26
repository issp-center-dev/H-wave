from functools import wraps
import time

class PerfDB:
    def __init__(self):
        self._db_count = {}
        self._db_value = {}
    def __del__(self):
        print("--------------------------------------------------------------------------------")
        print("Statistics")
        print("  function                 :  total elapsed  : average elapsed : ncalls")
        print("--------------------------------------------------------------------------------")
        for item in self._db_count.keys():
            print("  {:32s} : {:10.3f} msec : {:10.3f} msec : {:6d}".format(
                item,
                self._db_value[item] * 1000,
                self._db_value[item] * 1000 / self._db_count[item],
                self._db_count[item]
            ))
        print("--------------------------------------------------------------------------------")
    def put(self, name, value):
        self._db_count[name] = self._db_count.get(name, 0) + 1
        self._db_value[name] = self._db_value.get(name, 0) + value

_perf_db_data = PerfDB()

def do_profile(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        # start time
        _start = time.perf_counter()

        # exec function body
        result = func(*args, **kwargs)

        # end time
        _end = time.perf_counter()

        # elapsed time
        _elapsed = _end - _start

        #print("{}: elapsed time {:.3f} ms".format(func.__name__, _elapsed * 1000))
        _perf_db_data.put(func.__module__ + '.' + func.__name__, _elapsed)

        return result
    return wrapper
