#-----------------------------------------------------------------------
# This file is part of Nazca.
#
# Nazca is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# Nazca is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with Nazca.  If not, see <http://www.gnu.org/licenses/>.
#
# 2019 (c) Ronald Broeke
#-----------------------------------------------------------------------
"""
Set up logging to stdout or a file of the output from Nazca.

Logging is based on the standard Python logging concepts.
"""

import os
import sys
import inspect
import datetime
import logging as log  # rename to not confuse this Nazca logging module with the python one in later nd.logging calls.
import matplotlib.pyplot as plt
from . import cfg
from nazca import __version__

nazca_dir = os.path.dirname(__file__)
logging_module_dir = os.path.abspath(log.__file__)

def clear_handlers(logger):
    logger.handlers.clear()
    logger.filters.clear()
    #if logger.handlers:
        #logger.removeHandler(logger.handlers[0])


class ContextFilter(log.Filter):
    """Filter to inject contextual information into the log."""
    cnt = 0
    def filter(self, record):
        ContextFilter.cnt += 1

        line = ""
        callerfilename = ""
        record.callerfilenameline = ""

        if record.levelno >= cfg.logging_level_linenumber:  # only add line numbers for logging_level_linenumber and up.
            n = 2
            stack = inspect.stack()
            nmax = len(stack)
            flag = False
            while n < 10:
                callerframerecord = stack[n]  # n=0 represents this, n=1 represents line at caller, etc
                frame = callerframerecord[0]
                info = inspect.getframeinfo(frame)
                filename = os.path.basename(info.filename)
                n += 1
                if ContextFilter.cnt > 3 and logging_module_dir not in info.filename:
                    flag = True
                if flag and nazca_dir not in info.filename:
                    line = f"{info.lineno}"
                    callerfilename = filename
                    record.callerfilenameline = f" file/#:{callerfilename}/{line}:"
                    break
                if n == nmax:  # Iterator out of range
                    break
        else:
            record.cnt = f"{ContextFilter.cnt:d}"
        return True


loggers = [] # keep track of loggers to not initialize more than once.

logger = log.getLogger('nazca')
#logger.filters = []  # clear filters
clear_handlers(logger=logger)

#logname = 'nazca.log'
#filehandler = log.FileHandler(logname, 'w')
#filehandler.setFormatter(formatter0)
#logger.addHandler(filehandler)


formatter0 = log.Formatter('%(levelname)-7s:%(callerfilenameline)s %(message)s')
streamhandler = log.StreamHandler(sys.stdout)
streamhandler.setFormatter(formatter0)

f = ContextFilter()
streamhandler.addFilter(f)

logger.addHandler(streamhandler)
logger.setLevel(log.INFO)

#f = ContextFilter()
#logger.addFilter(f)



def nolog2stdout():
    """Switch off the logger's stdout stream.

    By default the stdout is on. It is recommended to always use a log file and
    use the logfile() method instead with stdout=False setting.

    Returns:
        None
    """
    logger.removeHandler(log.streamhandler)


def logfile(name='nazca.log', level='INFO', formatter=None, stdout=True, create=False):
    """Create a file handler for logging to file.

    To start logging as soon as possible the logger should be placed directly
    after import nazca and before loading e.g. a PDK.

    Args:
        name (str): filename of the logfile. If the name does not
            and with '.log' it will be added automatically.
        level (str): minimum log level (case insensitive) reported:
            DEBUG, INFO (default), WARNING, ERROR, CRITICAL
        formatter (str): format string, if None a default formatter will be used
        stdout (bool): set False to switch off a copy of the log to stdout
            (default=True)
        create (bool):

    Example:
         Create a filehandler with the filename of your python file and switch
         off streamhandler to standard output (stdout)::

             import nazca as nd
             nd.logfile(name=__file__, stdout=False)

             # rest of code

    Returns:
        logger: logger
    """
    global globname, datetimestamp

    name = os.path.basename(name)
    if not name.endswith('.log'):
        name = f"{name}.log"
    globname = name

    if formatter is None:
        formatter = formatter0

    if not create:
        logger0 = logger
        if not stdout:
            logger0.removeHandler(streamhandler)
    else:
        logger0 = log.getLogger('new logger')
        clear_handlers(logger0)

    if logger0 in loggers:
        return logger0
    else:
        loggers.append(logger0)

    filehandler = log.FileHandler(name, 'w')
    filehandler.setFormatter(formatter0)

    f = ContextFilter()
    filehandler.addFilter(f)
    logger0.addHandler(filehandler)

    levels = {
        'DEBUG': log.DEBUG,
        'INFO': log.INFO,
        'WARNING': log.WARNING,
        'ERROR': log.ERROR,
        'CRITICAL': log.CRITICAL,
    }
    level = level.upper()
    if level not in levels.keys():
        raise Exception(f"logging level '{level}' not recognized, falling back to 'INFO'.")
    else:
        logger0.setLevel(levels[level])

    # Add header info the the logger:
    print(f"...creating log file '{name}'")
    if not create:
        now = datetime.datetime.now()
        logger0.info('datetime: {}'.format(now.strftime("%Y-%m-%d %H:%M")))
        logger0.info('logging level: {}'.format(level))
        logger0.info('Nazca Design version: {}'.format(__version__))
        datetimestamp = now
    return logger0


def newlogfile(name='nazca.log', level='INFO', formatter=None, stdout=True):
    return logfile(name=name, level=level, formatter=formatter, stdout=stdout, create=True)


# Log "bare" exceptions, but only in development version
dir_path = os.path.dirname(os.path.abspath(__file__))
if os.path.exists(dir_path[:-5] + 'tests'):
    print_except = logger.exception
else:
    def print_except(e):
        pass


def summary(filename="", plot=True):
    """Summarize errors and warnings in the logfile.

    Adds log info entry for number of warning and errors in the logfile.

    Args:
        filename (str): logfile, default the main logger.
        plot (bool): show warnings + errors in a bar chart, default=True

    Retrun:
        None
    """
    from matplotlib.ticker import MaxNLocator
    global globname, datetimestamp

    if filename == "":
        filename = globname
    warnings = 0
    errors = 0
    with open(filename, 'r') as F:
        for line in F:
            if "ERROR" in line:
                errors += 1
            elif "WARNING" in line:
                warnings += 1
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_title(f"Logfile '{globname}' - {datetimestamp.strftime('%Y-%m-%d %H:%M')}")
    level = ['WARNINGS', 'ERRORS']
    occurance = [warnings, errors]
    barlist = ax.bar(level, occurance)
    barlist[0].set_color('y')
    barlist[1].set_color('r')
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.show()
    logger.info(f"total number of WARNINGS: {warnings}")
    logger.info(f"total number of ERROR: {errors}")
