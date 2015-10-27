from __future__ import division
import os
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg  as la


class EigenProblem(object):

	def __init__(self, mesh, m0, mu_s=1.0, J=1, D=0, Kz=0, Kx=0, H=None, gamma=1.0):
		self.mesh = mesh
		self.m0 = m0
		self.mu_s = mu_s
		self.nx = mesh.nx
		self.ny = mesh.ny
		self.n = mesh.n
		
		self.ngx = mesh.neighbours_x
		self.ngy = mesh.neighbours_y
		self.ngz = mesh.neighbours_z
		self.ngs = mesh.neighbours

		self.H = H
		self.Kx = Kx/self.mu_s
		self.Kz = Kz/self.mu_s
		self.J = J/self.mu_s
		self.D = D/self.mu_s

		assert len(m0)/3 == mesh.n

		#self.M = sp.csr_matrix((2*self.n, 2*self.n), dtype=np.float)
		self.M = sp.lil_matrix((2*self.n, 2*self.n), dtype=np.float)

		self.compute_theta_phi()
		self.build_matrix()

		#self.solve_eigen()


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
		for i in range(self.n):
			j = i+self.n
			total = self.H[0]*cp[i]*st[i] + self.H[1]*sp[i]*st[i] + self.H[2]*ct[i]
			self.M[i,j] -= total
			self.M[j,i] += total

	def add_H0_K(self):
		st = self.st
		ct = self.ct
		sp = self.sp
		cp = self.cp
		for i in range(self.n):
			j = i+self.n
			total = 2*self.Kx*(cp[i]*st[i])**2+2*self.Kz*(ct[i])**2
			self.M[i,j] -= total
			self.M[j,i] += total


	def add_H0_J(self):
		st = self.st
		ct = self.ct
		sp = self.sp
		cp = self.cp
		t = self.theta
		p = self.phi

		for i in range(self.n):
			
			total = 0
			for j in self.ngs[i]:
				total += self.J*(ct[i]*ct[j]+st[i]*st[j]*np.cos(p[i]-p[j]))

			k = i+self.n

			self.M[i,k] -= total
			self.M[k,i] += total


	def add_H0_D(self):
		st = self.st
		ct = self.ct
		sp = self.sp
		cp = self.cp
		p = self.phi

		for i in range(self.n):
			total = 0

			for j in self.ngx[i]:
				total += self.D*np.sign(j-i)*(st[j]*sp[j]*ct[i]-st[i]*sp[i]*ct[j])

			for j in self.ngy[i]:
				total += self.D*np.sign(j-i)*(st[i]*cp[i]*ct[j]-st[j]*cp[j]*ct[i])

			for j in self.ngz[i]:
				total += self.D*np.sign(j-i)*(st[i]*st[j]*np.sin(p[i]-p[j]))

			k = i + self.n

			self.M[i,k] -= total
			self.M[k,i] += total

	def add_h_K(self):
		st = self.st
		ct = self.ct
		sp = self.sp
		cp = self.cp
		for i in range(self.n):
			j = i+self.n

			#hu
			self.M[j,i] -= 2*self.Kx*(cp[i]*ct[i])**2 + 2*self.Kz*st[i]**2
			self.M[j,j] += 2*self.Kx*sp[i]*cp[i]*ct[i]
			
			#hv
			self.M[i,i] += -2*self.Kx*sp[i]*cp[i]*ct[i]
			self.M[i,j] += 2*self.Kx*sp[i]*sp[i]
			

	def add_h_J(self):
		st = self.st
		ct = self.ct
		sp = self.sp
		cp = self.cp
		p = self.phi

		for i in range(self.n):

			k = i+self.n

			for j in self.ngs[i]:

				#hu
				self.M[k, j] -= self.J*(np.cos(p[i]-p[j])*ct[i]*ct[j]+st[i]*st[j])
				self.M[k, j+self.n] -= self.J*ct[i]*np.sin(p[i]-p[j])

				#hv
				self.M[i, j] -= self.J*ct[j]*np.sin(p[i]-p[j])
				self.M[i, j+self.n] += self.J*np.cos(p[i]-p[j])
				


	def add_h_D(self):
		st = self.st
		ct = self.ct
		sp = self.sp
		cp = self.cp
		p = self.phi

		for i in range(self.n):

			k = i+self.n

			for j in self.ngx[i]:
				#hu
				self.M[k,j] -= self.D*np.sign(j-i)*(st[j]*sp[i]*ct[i]-st[i]*sp[j]*ct[j])
				self.M[k,j+self.n] += self.D*np.sign(j-i)*st[i]*cp[j]

				#hv
				self.M[i, j] += self.D*np.sign(j-i)*st[j]*cp[i]

			for j in self.ngy[i]:
				#hu
				self.M[k, j] -= self.D*np.sign(j-i)*(st[i]*cp[j]*ct[j]-st[j]*cp[i]*ct[i])
				self.M[k,j+self.n] += self.D*np.sign(j-i)*st[i]*sp[j]
				
				#hv
				self.M[i, j+self.n] += self.D*np.sign(j-i)*st[j]*sp[i]

			for j in self.ngz[i]:
				#hu
				self.M[k,j] -= self.D*np.sign(j-i)*ct[i]*ct[j]*np.sin(p[i]-p[j])
				self.M[k,j+self.n] += self.D*np.sign(j-i)*ct[i]*np.cos(p[i]-p[j])

				#hv
				self.M[i, j] += self.D*np.sign(j-i)*ct[j]*np.cos(p[i]-p[j])
				self.M[i, j+self.n] += self.D*np.sign(j-i)*np.sin(p[i]-p[j])


	def build_matrix(self):
		if self.H is not None:
			self.add_H0_H()
			print 'build_matrix for applied field: Done!'
		if self.J !=0:
			self.add_H0_J()
			self.add_h_J()
			print 'build_matrix for exchange: Done!'
			
		if self.Kx!=0 or self.Kz !=0:
			self.add_H0_K()
			self.add_h_K()
			print 'build_matrix for anisotropy: Done!'

		if self.D !=0:
			self.add_H0_D()
			self.add_h_D()
			print 'build_matrix DMI: Done!'

	def solve_sparse(self, w0=1e-8, n=20):
		
		w, v = la.eigs(self.M, k=2*n,  sigma = 0+w0*1j,  maxiter=10000, OPpart='r')

		freqs_all = np.imag(w)

		freqs = []
		vs = []
		for i, f in enumerate(freqs_all):
			if f > 0:
				freqs.append(f)
				vs.append(v[:,i])

		return np.array(freqs), np.array(vs)

	def solve_dense(self):

		M = self.M.toarray()

		print 'max:::::',np.max(np.abs(M.transpose()+M))

		w,v = np.linalg.eig(M)

		freqs_all = np.imag(w)

		idx = freqs_all.argsort()
		ids = [i for i in idx if freqs_all[i]>0]

		freqs = freqs_all[ids]
		vectors = np.transpose(v[:,ids])

		return freqs, vectors
		


if __name__ == '__main__':
	pass
