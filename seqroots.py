import itertools as it
from Bio import SeqIO
import numpy as np

def parse_seq(subseq,MODE,output=None,idx=None):
	for seq in SeqIO.parse(fasta,'fasta'):
		if subseq in seq.seq:
			if MODE == 'count':
				return 1
			else:
				output[idx] = subseq
				return None
	return 0

class seqroots:

	global aa
	aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

	def __init__(self,fasta,**kwargs):
		self.fasta = fasta
		self.set_kwargs(kwargs)
		self.scipy = False

	def set_kwargs(self,kwargs):
		allowed_keys = ['MODE','nthreads','min_subseq_len','max_subseq_len']
		default_keys = ['count', 1, 2, 3]
		self.__dict__.update((key, value) for key,value in zip(allowed_keys,default_keys))
		self.__dict__.update((key, value) for key, value in kwargs.items() if key in allowed_keys)

	def import_libs(self):
		if self.nthreads > 1:
			try:
				global Parallel, delayed, has_shareable_memory
				from joblib import Parallel, delayed
			except ImportError:
				print "can't import joblib. Will use 1 thread"
				self.nthreads = 1
		if self.MODE == 'count':
			global timer, plt, font, sy, curve_fit
			from timeit import default_timer as timer
			import matplotlib as mpl ; mpl.use('Agg')
			font = {'family' : 'normal','weight' : 'normal','size'   : 15}
			mpl.rc('font', **font) ; import matplotlib.pyplot as plt
			try:
				import scipy as sy
				from scipy.optimize import curve_fit
				self.scipy = True
			except ImportError:
				print "can't import scipy curve_fit. Will skip this step"

		elif self.MODE == 'save':
			global os, shutil
			import os, shutil
	
	@staticmethod
	def fit(y):
		nbins = np.arange(len(y))+2
		x = np.asarray(nbins).ravel()
		def func(x,a,b): return a*20**(b*x)
		p0 = sy.array([1,1])
		coeffs, matcov = curve_fit(func, x, y, p0) 
		yaj = func(x, coeffs[0], coeffs[1])
		plt.plot(x,y,'--',x,yaj, 'r-')
		plt.xlabel('Subsequence Length')
		if len(y) > 2:
			plt.text(3,0.3, r'y = %.2f*$20^{%.2fx}$' % (coeffs[0],coeffs[1]), fontdict=font)
		return None
	
	def compute_scaling(self,y,pltname):
		plt.clf()
		try:
			self.fit(y)
		except TypeError:
			print "need more points to run fitting algorithm"
			return None
		plt.savefig(pltname)
		return None

	def save_memmap_readable(self,N,output,outname):
		print "writing...", outname
		output = np.array(output)
		output = output[np.where(output != '')]
		np.savetxt(outname,output, delimiter=" ", fmt="%s")
		return None

	def run(self):
		self.import_libs()

		if self.MODE == "count":
			t = [] ; n = []
			for N in range(self.min_subseq_len,self.max_subseq_len+1):
				print "working on", N
				start = timer()
				if self.nthreads == 1:
					n.append(float(sum([parse_seq(''.join(subseq),self.MODE) for subseq in it.combinations(aa,N)]))/(len(aa)**N))
				elif self.nthreads > 1:
					n.append(float(sum(Parallel(n_jobs=self.nthreads)(delayed(parse_seq)(''.join(subseq),self.MODE) for subseq in it.combinations(aa,N))))/(len(aa)**N))
				end = timer()
				t.append(end-start)
				print N, n[-1], end-start
			if self.scipy:
				#self.compute_scaling(t,'time_scaling.png')
				self.compute_scaling(n,'number_scaling.png')
		
		elif self.MODE == "save":
			N = self.max_subseq_len
			folder = './joblib_memmap'
			try:
				os.mkdir(folder)
			except OSError:
                                pass
			output_filename_memmap = os.path.join(folder, 'output_memmap')
			output = np.memmap(output_filename_memmap, mode='w+',shape=(int(20**N),),dtype="S"+str(N))
			outname = 'N'+str(N)+'.npy'
			if os.path.exists(outname):
				overwrite = raw_input("'"+outname+"' exists. Overwrite? [y/n] ")
				if overwrite == 'n': 
					print "exiting..."
					exit()
			print "working..."
			if self.nthreads > 1:
				Parallel(n_jobs=self.nthreads)(delayed(parse_seq)(''.join(subseq),self.MODE,output,idx) for idx,subseq in enumerate(it.combinations(aa,N)))
			else:
				[parse_seq(''.join(subseq),self.MODE,output,idx) for idx,subseq in enumerate(it.combinations(aa,N))]
			try:
				shutil.rmtree(folder)
			except:
				print 'Could not clean-up automatically'
			self.save_memmap_readable(N,output,outname)

if __name__ == "__main__":
	fasta = "/home/dillion/Dropbox/james_project/uniprot_sprot.fasta"
	s = seqroots(fasta, nthreads=2, MODE='save', min_subseq_len = 2, max_subseq_len = 3)
	s.run()
