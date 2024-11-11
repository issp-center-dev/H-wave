"""Performance profiling utilities.

This module provides utilities for profiling function execution times and
collecting performance statistics.
"""

from functools import wraps
import time


class PerfDB:
    """Database for storing performance profiling data.
    
    Stores execution counts and total elapsed time for profiled functions.
    Prints summary statistics on deletion.
    
    Attributes
    ----------
    _db_count : dict
        Number of calls per function
    _db_value : dict
        Total elapsed time per function
    """

    def __init__(self):
        """Initialize empty performance database."""
        self._db_count = {}
        self._db_value = {}

    def __del__(self):
        """Print summary statistics when object is deleted."""
        if len(self._db_count) == 0:
            return
        print("--------------------------------------------------------------------------------")
        print("Statistics")
        print("  function                         :  total elapsed  : average elapsed : ncalls")
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
        """Add a timing measurement.

        Parameters
        ----------
        name : str
            Function name
        value : float
            Elapsed time in seconds
        """
        self._db_count[name] = self._db_count.get(name, 0) + 1
        self._db_value[name] = self._db_value.get(name, 0) + value


_perf_db_data = PerfDB()


def do_profile(func):
    """Decorator for profiling function execution time.

    Parameters
    ----------
    func : callable
        Function to profile

    Returns
    -------
    callable
        Wrapped function that measures and records execution time
    """
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
