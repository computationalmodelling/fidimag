#! /usr/bin/python

import math
import os
import sys
import numpy
import struct


class OMF2:

    def __init__(self, file_name, filter=False):
        self.file_name = file_name
        self.read()
        if filter:
            pass

    def read(self):
        f = open(self.file_name, 'rb')
        line = f.readline()
        # print line
        if line != "# OOMMF OVF 2.0\n":
            print self.file_name + ": NOT OOMMF OVF 2.0!"
            sys.exit(0)

        while not line.startswith("# Begin: Data Binary"):
            line = f.readline()

            if line.startswith("# xnodes:"):
                self.xnodes = int(line.split(":")[1])
            elif line.startswith("# ynodes:"):
                self.ynodes = int(line.split(":")[1])
            elif line.startswith("# znodes:"):
                self.znodes = int(line.split(":")[1])
            elif line.startswith("# ystepsize:"):
                self.ystepsize = float(line.split(":")[1])

        msb = f.read(8)

        if struct.unpack('d', msb)[0] != 123456789012345.0:
            print 'check value error!'
            return

        count = 3 * self.xnodes * self.ynodes * self.znodes
        data = f.read(8 * count)
        self.data = numpy.frombuffer(data)
        self.data = numpy.reshape(self.data, (-1, 3))
        f.close()
        return

    def get_mag(self, id_x, id_y, id_z, comp='x'):
        """
        return x, y or z component of magnetisation at index (i,j,k)
        """

        index = self.xnodes * self.ynodes * id_z + self.xnodes * id_y + id_x

        id_comp = ord(comp) - ord('x')

        return self.data[index][id_comp]

    def get_all_mag(self, comp='x'):
        """
        return x, y or z component of magnetisation of all nodes 
        """
        index = ord(comp) - ord('x')
        return self.data[:,index]
	
    def get_all_mags(self, order='xyz'):
        if order=='xyz':
            d = self.data.copy()
            d.shape=(-1)
            return d
        elif order=='xxx':
            return np.array([self.data[:,0], self.data[:,1], self.data[:,2]])
