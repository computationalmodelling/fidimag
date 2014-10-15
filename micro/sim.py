from llg import LLG
from llg_stt import LLG_STT
from baryakhtar import LLBar
from baryakhtar import LLBarFull

KNOWN_DRIVERS = {'llg': LLG,
                 'llg_stt': LLG_STT, 
                 'llbar': LLBar,
                 'llbar_full': LLBarFull}

def Sim(*args, **kwargs):
    
    driver = 'llg'
    
    if kwargs.has_key('driver'):
        driver = kwargs['driver']
        kwargs.pop('driver')
    
    if driver not in KNOWN_DRIVERS:
        raise NotImplementedError("""Driver '{}' is not implemented.
                                  Valid choices: one of '{}'.""".format(driver, KNOWN_DRIVERS.keys()))
    
    return KNOWN_DRIVERS[driver](*args, **kwargs)


if __name__=='__main__':
    pass