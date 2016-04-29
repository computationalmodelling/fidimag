from . import llg
from . import llg_stt
from . import llg_stt_cpp
from . import baryakhtar

KNOWN_DRIVERS = {'llg': llg.LLG,
                 'llg_stt': llg_stt.LLG_STT,
                 'llg_stt_cpp': llg_stt_cpp.LLG_STT_CPP,
                 'llbar': baryakhtar.LLBar,
                 'llbar_full': baryakhtar.LLBarFull}


def Sim(*args, **kwargs):
    """
    def Sim(*args, **kwargs)
    
    This function returns a simulation object.

    By default, it will return a simulation object for the LLG driver;
    other available drivers are:
        llg_stt - LLG w. spin transfer torque
        llg_stt_cpp - LLG w. spin transfer torque 
        llbar - Landau-Lifshitz-Baryakhtar equation
        llbar_full
    
    To use the STT driver, for example, pass as an argument driver="llg_stt".

    """
    driver = 'llg'

    if 'driver' in kwargs:
        driver = kwargs.pop('driver')

    if driver not in KNOWN_DRIVERS:
        raise NotImplementedError("""Driver '{}' is not implemented.
                                  Valid choices: one of '{}'.""".format(driver, KNOWN_DRIVERS.keys()))

    return KNOWN_DRIVERS[driver](*args, **kwargs)
