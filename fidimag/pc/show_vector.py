from __future__ import division

import numpy as np
import visual as vis


class VisualSpin(object):

    def __init__(self, sim):
        self.sim = sim
        mesh = sim.mesh
        self.mesh = mesh
        self.spin = sim.spin
        self.nxyz = mesh.nxyz

    def init(self):
        vis.scene.autoscale = False
        self.vspins = []
        tmp = np.reshape(self.spin, (self.nxyz, 3), order='F')
        for i in range(self.nxyz):
            a = vis.arrow(
                pos=self.mesh.pos[i], axis=tmp[i], color=vis.color.yellow)
            self.vspins.append(a)

    def update(self):
        tmp = np.reshape(self.spin, (self.nxyz, 3), order='F')
        for i in range(len(tmp)):
            self.vspins[i].axis = tmp[i]
