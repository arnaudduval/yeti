from lib.__init__ import *

class solver():
	def __init__(self):
		self._nbIter       = 100
		self._thresholdPCG = 1e-10
		return
	
	def CG(self, Afun, b):
		x = np.zeros(np.shape(b))
		r = b; normb = np.max(np.abs(r))
		resPCG = [1.0]
		if normb <= self._thresholdPCG: return
		rsold = np.dot(r, r); p = r

		for i in range(self._nbIter):
			Ap = Afun(p)
			alpha = rsold/np.dot(p, Ap)
			x += alpha*p
			r -= alpha*Ap

			resPCG.append(np.max(np.abs(r))/normb)
			if (resPCG[-1]<=self._thresholdPCG): break

			rsnew = np.dot(r, r)
			p = r + rsnew/rsold*p
			rsold = rsnew

		return x, resPCG
	
	def PCG(self, Afun, Pfun, b):
		x = np.zeros(np.shape(b))
		r = b; normb = np.max(np.abs(r))
		resPCG = [1.0]
		if normb <= self._thresholdPCG: return
		z = Pfun(r)
		rsold = np.dot(r, z); p = z

		for i in range(self._nbIter):
			Ap = Afun(p)
			alpha = rsold/np.dot(p, Ap)
			x += alpha*p
			r -= alpha*Ap

			resPCG.append(np.max(np.abs(r))/normb)
			if (resPCG[-1]<=self._thresholdPCG): break

			z = Pfun(r)
			rsnew = np.dot(r, z)
			p = z + rsnew/rsold*p
			rsold = rsnew

		return x, resPCG
	
	def BiCGSTAB(self, Afun, b):
		x = np.zeros(np.shape(b))
		r = b; normb = np.max(np.abs(r))
		resPCG = [1.0]
		if normb <= self._thresholdPCG: return
		rhat = r; p = r
		rsold = np.dot(r, rhat)

		for i in range(self._nbIter):
			Ap = Afun(p)
			alpha = rsold/np.dot(Ap, rhat)
			s = r - alpha*Ap
			As = Afun(s)
			omega = np.dot(As, s)/np.dot(As, As)
			x += alpha*p + omega*s
			r = s - omega*As

			resPCG.append(np.max(np.abs(r))/normb)
			if (resPCG[-1]<=self._thresholdPCG): break

			rsnew = np.dot(r, rhat)
			beta = (alpha/omega)*(rsnew/rsold)
			p = r + beta(p - omega*Ap)
			rsold = rsnew

		return x, resPCG
	
	def PBiCGSTAB(self, Afun, Pfun, b):
		x = np.zeros(np.shape(b))
		r = b; normb = np.max(np.abs(r))
		resPCG = [1.0]
		if normb <= self._thresholdPCG: return
		rhat = r; p = r
		rsold = np.dot(r, rhat)

		for i in range(self._nbIter):
			ptilde = Pfun(p)
			Aptilde = Afun(ptilde)
			alpha = rsold/np.dot(Aptilde, rhat)
			s = r - alpha*Aptilde
			stilde = Pfun(s)
			Astilde = Afun(stilde)
			omega = np.dot(Astilde, s)/np.dot(Astilde, Astilde)
			x += alpha*ptilde + omega*stilde
			r = s - omega*Astilde

			resPCG.append(np.max(np.abs(r))/normb)
			if (resPCG[-1]<=self._thresholdPCG): break

			rsnew = np.dot(r, rhat)
			beta = (alpha/omega)*(rsnew/rsold)
			p = r + beta*(p - omega*Aptilde)
			rsold = rsnew

		return x, resPCG