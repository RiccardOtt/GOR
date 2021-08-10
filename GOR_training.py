import sys
import numpy as np


def get_pssm(pssm_file):
	seq_prof = []
	for line in pssm_file:
		line = line.split()
		try:
			line[0]
		except:
			pass
		else:
			if line[0].isdigit():
				seq_prof.append(line[22:-2])

	x = np.array(seq_prof, dtype=np.float64)
	x /= 100
#	print(seq_prof)
	return x


def matrices(SeqProf, dssp, H, E, C, R, ProbSS):
	Padding = np.zeros(shape=(8,20))
	SS = dssp.read().splitlines()[1]
#	SS = SS.replace('-','C')
	PadSS = '-'*8 + SS + '-'*8
#	print(SS)
#	print(SeqProf)
	PSeqProf = np.vstack((Padding,SeqProf,Padding))
#	print(PSeqProf)
	for i in range(8,len(PSeqProf)-8):
		j = i-8
		k = i+8
		if PadSS[i] == 'H':
			H += PSeqProf[j:k+1]
			ProbSS[0] += 1
		if PadSS[i] == 'E':
			E += PSeqProf[j:k+1]
			ProbSS[1] += 1
		if PadSS[i] == '-':
			C += PSeqProf[j:k+1]
			ProbSS[2] += 1
		R += PSeqProf[j:k+1]

	return H,E,C,R,ProbSS




if __name__ == '__main__':
	fileid = sys.argv[1]
	H = np.zeros(shape=(17,20))
	E = np.zeros(shape=(17,20))
	C = np.zeros(shape=(17,20))
	R = np.zeros(shape=(17,20))
	ProbSS = np.zeros(3)
	with open(fileid) as filein:
		for id in filein:
			id = id.rstrip()
#			print(id)
			profile_file = '/home/riccardo/Documents/Documents/LB2/Castrense/project/jpred4.pssm/' + id + '.pssm'
			dssp_file = '/home/riccardo/Documents/Documents/LB2/Castrense/project/jpred4.dssp/' + id + '.dssp'
			try:
				prof = open(profile_file)
				dssp = open(dssp_file)
			except: continue

			else:
				profile = get_pssm(prof)
				H,E,C,R,ProbSS = matrices(profile,dssp,H,E,C,R,ProbSS)


	tot_res = ProbSS[0] + ProbSS[1] + ProbSS[2]
	H /= tot_res
	E /= tot_res
	C /= tot_res
	R /= tot_res
	ProbSS /= tot_res
	print(H) #,E,C,R,ProbSS)

	np.save('GOR_tr_H',H)
	np.save('GOR_tr_E',E)
	np.save('GOR_tr_C',C)
	np.save('GOR_tr_R',R)
	np.save('GOR_tr_ProbSS',ProbSS)
