from __future__ import division
import os
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg  as la


class EigenProblem(object):

	def __init__(self, mesh, m0, mu_s=1.0, J=1, D=0, K=0, H=None, gamma=1.0):
		self.mesh = mesh
		self.m0 = m0
		self.mu_s = mu_s
		self.nx = mesh.nx
		self.ny = mesh.ny
		self.nxyz = mesh.nxyz
		self.ngx = mesh.neighbours_x
		self.ngy = mesh.neighbours_y
		self.ngxy = mesh.neighbours_xy

		self.H = H
		self.K = K/self.mu_s
		self.J = J/self.mu_s
		self.D = D/self.mu_s

		assert len(m0)/3 == mesh.nxyz

		#self.M = sp.csr_matrix((2*self.nxyz, 2*self.nxyz), dtype=np.float)
		self.M = sp.lil_matrix((2*self.nxyz, 2*self.nxyz), dtype=np.float)

		self.compute_theta_phi()
		self.build_matrix()

		self.solve_eigen()


	def compute_theta_phi(self):
		xyz = self.m0
		
		xyz.shape = (3, -1)
		r_xy = np.sqrt(xyz[0, :] ** 2 + xyz[1, :] ** 2)
		theta = np.arctan2(r_xy, xyz[2, :])
		phi = np.arctan2(xyz[1, :], xyz[0, :])
		xyz.shape = (-1,)

		self.theta = theta
		self.phi = phi
		self.st = np.sin(theta)
		self.ct = np.cos(theta)
		self.sp = np.sin(phi)
		self.cp = np.cos(phi)


	def add_H0_H(self):
		st = self.st
		ct = self.ct
		sp = self.sp
		cp = self.cp
		for i in range(self.nxyz):
			j = i+self.nxyz
			tmp = self.H[0]*cp[i]*st[i] + self.H[1]*sp[i]*st[i] + self.H[2]*ct[i]
			self.M[i,j] += tmp
			self.M[j,i] -= tmp

	def add_H0_K(self):
		st = self.st
		ct = self.ct
		sp = self.sp
		cp = self.cp
		for i in range(self.nxyz):
			j = i+self.nxyz
			tmp = 2*self.K*cp[i]*cp[i]*st[i]*st[i]
			self.M[i,j] += tmp
			self.M[j,i] -= tmp


	def add_H0_J(self):
		st = self.st
		ct = self.ct
		sp = self.sp
		cp = self.cp
		t = self.theta
		p = self.phi

		for i in range(self.nxyz):
			
			total = 0
			for k in self.ngxy[i]:
				total += self.J*(ct[i]*ct[k]+st[i]*st[k]*np.cos(p[i]-p[k]))

			j = i+self.nxyz

			self.M[i,j] += total
			self.M[j,i] -= total


	def add_H0_D(self):
		st = self.st
		ct = self.ct
		sp = self.sp
		cp = self.cp
		p = self.phi

		for i in range(self.nxyz):
			total = 0

			for j in self.ngx[i]:
				total += self.D*np.sign(i-j)*(st[i]*sp[i]*ct[j]-st[j]*sp[j]*ct[i])

			for j in self.ngy[i]:
				total += self.D*np.sign(i-j)*(-st[i]*cp[i]*ct[j]+st[j]*cp[j]*ct[i])

			k = i+self.nxyz

			self.M[i,k] += total
			self.M[k,i] -= total

	def add_h_J(self):
		st = self.st
		ct = self.ct
		sp = self.sp
		cp = self.cp
		p = self.phi

		for i in range(self.nxyz):

			k = i+self.nxyz

			for j in self.ngxy[i]:
				#hy
				self.M[i,j] -= self.J*ct[i]*np.sin(p[i]-p[j])
				self.M[i,k] -= self.J*(np.cos(p[i]-p[j])*ct[i]*ct[j]+st[i]*st[j])

				#hx
				self.M[k, j] += np.cos(p[i]-p[j])
				self.M[k, j+self.nxyz] += ct[j]*np.sin(p[i]-p[j])

	def add_h_D(self):
		st = self.st
		ct = self.ct
		sp = self.sp
		cp = self.cp
		p = self.phi

		for i in range(self.nxyz):

			k = i+self.nxyz

			for j in self.ngx[i]:
				#hy
				self.M[i,j] += self.D*np.sign(i-j)*st[i]*cp[j]
				self.M[i,k] -= self.D*(st[i]*sp[j]*ct[j]-st[j]*sp[i]*ct[i])

				#hx
				self.M[k, j+self.nxyz] += self.D*np.sign(i-j)*st[j]*cp[i]

			for j in self.ngy[i]:
				#hy
				self.M[i,j] += self.D*np.sign(i-j)*st[i]*sp[j]
				self.M[i,k] -= self.D*(-st[i]*cp[j]*ct[j]+ct[j]*cp[i]*ct[i])

				#hx
				self.M[k, j+self.nxyz] += self.D*np.sign(i-j)*st[j]*sp[i]


	def build_matrix(self):
		if self.H is not None:
			self.add_H0_H()
		if self.K !=0:
			self.add_H0_K()
		if self.J !=0:
			self.add_H0_J()
			self.add_h_J()
		if self.D !=0:
			self.add_H0_D()
			self.add_h_D()
		print 'build_matrix: Done!'

	def solve_eigen(self):

		vals, vecs = la.eigs(self.M, k=10, ncv=self.nxyz,  which='SI')
		print vals
		#print vecs[:,0]
		#print self.M.toarray()

		#w,v = np.linalg.eig(self.M.toarray())
		#print w


if __name__ == '__main__':
	pass
