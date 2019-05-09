import glob
import os

def glob_files(path, excludes, extension="*.cpp"):
    files = []
    for srcfile in glob.glob(os.path.join(path, extension)):
        filename = os.path.basename(srcfile)
        if not filename in tuple(excludes):
            files.append(srcfile)
    return files


class BuildError(Exception):
    pass



def get_user_module_sources(folder):
    print('Found User Module: {}'.format(folder))
    user_sources = glob.glob(folder + '/*.pyx')
    print('\tFound Cython sources: {}'.format(user_sources))

    if len(user_sources) != 1:
        raise BuildError("User Modules are only allowed one Cython .pyx file")

    filename_string = user_sources[0].split('/')[-1][:-4]
    if filename_string != module_name:
        print(filename_string, module_name)
        raise BuildError("The Cython source file in {} must match the folder name - i.e. it must be {}.pyx".format(module_name, module_name))
    cfilename = filename_string + '.cpp'
    print(cfilename)
    user_sources += glob_files(folder, excludes=[cfilename])
    
    return user_sources
