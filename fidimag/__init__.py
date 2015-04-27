try:
    from . import extensions
except ImportError as e:
    import os
    cwd = os.getcwd()
    FIDIMAG_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(FIDIMAG_DIR)
    print "Could not find fidimag extensions. Trying to build them now..."
    os.system("make build")
    print "Building extensions done."
    os.chdir(cwd)
from . import util
from . import pc
from . import micro
