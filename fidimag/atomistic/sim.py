from .llg import LLG
from .sllg import SLLG
from .llg_stt import LLG_STT
from .llg_stt_cpp import LLG_STT_CPP


KNOWN_DRIVERS = {'llg': LLG,
                 'sllg': SLLG,
                 'llg_stt': SLLG_STT,
                 'llg_stt_cpp': LLG_STT_CPP
		}


def Sim(*args, **kwargs):

    driver = 'llg'

    if 'driver' in kwargs:
        driver = kwargs.pop('driver')

    if driver not in KNOWN_DRIVERS:
        raise NotImplementedError("""Driver '{}' is not implemented.
                                  Valid choices: one of '{}'.""".format(
                                      driver, KNOWN_DRIVERS.keys()))

    return KNOWN_DRIVERS[driver](*args, **kwargs)
