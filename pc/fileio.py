import os
import numpy as np
from types import TupleType, StringType

class DataSaver(object):
    
    comment_symbol = '# '

    def __init__(self, sim, filename):
        
        self.sim = sim
        self.filename = filename
        
        precision=12
        charwidth = 18
        self.float_format = "%" + str(charwidth)+'.'+str(precision)+ "g "
        self.string_format = "%" + str(charwidth) + "s "

        self.entities = {
            'time': {'unit': '<s>',
                        'get': lambda sim: sim.t,
                        'header': 'time'},
            'm': {'unit': '<>',
                  'get': lambda sim: sim.compute_average(),
                  'header': ('m_x', 'm_y', 'm_z')}
            }

        self.save_head=False
        self.entity_order = self.default_entity_order()
        
    def default_entity_order(self):
        keys = self.entities.keys()
        # time needs to go first
        keys.remove('time')
        return ['time'] + sorted(keys)
    
    def update_entity_order(self):
        self.entity_order = self.default_entity_order()

    def headers(self):
        """return line one and two of ndt data file as string"""
        line1 = [self.comment_symbol]
        line2 = [self.comment_symbol]
        for entityname in self.entity_order:
            colheaders = self.entities[entityname]['header']
            # colheaders can be a 3-tuple ('mx','my','mz'), say
            # or a string ('time'). Avoid iterating over string:
            if isinstance(colheaders, str):
                colheaders = [colheaders]
            for colhead in colheaders:
                line1.append(self.string_format % colhead)
                line2.append(self.string_format % \
                    self.entities[entityname]['unit'])
        return "".join(line1) + "\n" + "".join(line2) + "\n"

    def save(self):
        """Append data (spatial averages of fields) for current configuration"""
        
        if not self.save_head:
            f = open(self.filename, 'w')
            # Write header
            f.write(self.headers())
            f.close()
            self.save_head=True
        
        # open file
        with open(self.filename, 'a') as f:
            f.write(' ' * len(self.comment_symbol))  # account for comment symbol width
            for entityname in self.entity_order:
                value = self.entities[entityname]['get'](self.sim)
                if isinstance(value, np.ndarray):
                    if len(value) == 3:  # 3d vector
                        for i in range(3):
                            f.write(self.float_format % value[i])
                    else:
                        msg = "Can only deal with 3d-numpy arrays so far, but shape is %s" % value.shape
                        raise NotImplementedError(msg)
                elif isinstance(value, float) or isinstance(value, int):
                    f.write(self.float_format % value)
                else:
                    msg = "Can only deal with numpy arrays, float and int so far, but type is %s" % type(value)
                    raise NotImplementedError(msg)

            f.write('\n')


class DataReader(object):

    # open ndt file
    def __init__(self, filename):
        self.filename = filename
        # if file exists, cowardly stop
        if not os.path.exists(filename):
            raise RuntimeError("Cannot see file '%s'" % self.filename)
        # immediatey read file
        self.reload()

    def reload(self):
        """Read Table data file"""

        try:
            self.f = open(self.filename, 'r')
        except IOError:
            raise RuntimeError("Cannot see file '%s'" % self.filename)

        line1 = self.f.readline()
        line2 = self.f.readline()
        headers = line1.split()
        units = line2.split()

        assert len(headers) == len(units)

        # use numpy to read remaining data
        try:
            self.data = np.loadtxt(self.f)
        except ValueError:
            raise RuntimeError("Cannot load data from file '{}'. Maybe the file was incompletely written?".format(self.f))
        self.f.close()

        # some consistency checks: must have as many columns as
        # headers (disregarding the comment symbol)
        # for the case that only one line data 
        if len(self.data) == self.data.size:
            assert self.data.size == len(headers) - 1
            self.data.shape=(1, self.data.size)
        else:
            assert self.data.shape[1] == len(headers) - 1

        datadic = {}
        # now wrap up data conveniently
        for i, entity in enumerate(headers[1:]):
            datadic[entity] = self.data[:, i]

        self.datadic = datadic

    def entities(self):
        """Returns list of available entities"""
        return self.datadic.keys()


    def __getitem__(self, entity):
        """
        Given the entity name, return the data as a 1D numpy array.
        If multiple entity names (separated by commas) are given
        then a 2D numpy array is returned where the columns represent
        the data for the entities.
        """
        if isinstance(entity, StringType):
            res = self.datadic[entity]
        elif isinstance(entity, TupleType):
            res = [self.datadic[e] for e in entity]
        else:
            raise TypeError("'entity' must be a string or a tuple. "
                            "Got: {} ({})".format(entity, type(entity)))
        return res


