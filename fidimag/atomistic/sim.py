from .llg import LLG
from .sllg import SLLG
from .llg_stt_slonczewski_type import LLG_STT_Slonczewski


KNOWN_DRIVERS = {'llg': LLG,
                 'sllg': SLLG,
                 'slonczewski': LLG_STT_Slonczewski
		}


def Sim(*args, **kwargs):

    driver = 'llg'

    if kwargs.has_key('driver'):
        driver = kwargs['driver']
        kwargs.pop('driver')

    if driver not in KNOWN_DRIVERS:
        raise NotImplementedError("""Driver '{}' is not implemented.
                                  Valid choices: one of '{}'.""".format(
                                      driver, KNOWN_DRIVERS.keys()))

    return KNOWN_DRIVERS[driver](*args, **kwargs)
