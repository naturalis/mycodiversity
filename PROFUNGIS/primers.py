#!/usr/bin/python3
"""
Class Primer contains all methods needed to get the primer the user needs.
If it's not in the database, 
"""
class Primer:
	def __init__(self):
		self.primerfile = "./deps/primer.data"
		self.primerdict = self.do_primerdict(self.primerfile)

	@staticmethod
	def do_primerdict(primerfile):
		primerlist=[]
		primerdict={}
		with open(primerfile,'r') as primers:
			for line in primers:
				if not line.isspace():
					primerlist.append(line.strip("\n").split(':'))		
		for primer in primerlist:
			primerdict[primer[0]] = [primer[1],primer[2]]

		return primerdict

	def check_input(self, name, primer, fragment):
		goodchars = ['A','C','G','T','R','Y','W','D','S','B','V','K','H']
		goodfrags = ['its1','its2']
		for char in primer:
			if char not in goodchars:
				return False, "Sequence invalid"
		if name in self.primerdict:
			return False, "Name already present"
		present_primers = [x[0] for x in self.primerdict.values()]
		if primer in present_primers:
			present = list(self.primerdict)[present_primers.index(primer)]
			return False, "Primer already present ({primer})".format(primer=present)
		if len(primer) < 6:
			return False, "Primer too short, please add a sequence of >6 bases"

		if fragment.lower() not in goodfrags:
			return False, "Fragment not recognized. Please use ITS1 or ITS2"

		return True,"Primer added"

	def add_primer(self, name, primer, fragment):
		proper, message = self.check_input(name, primer, fragment)
		if proper:
			with open(self.primerfile,'a') as writefile:
				writefile.write("%s:%s:%s\n"%(name, primer,fragment.upper()))
			self.primerdict[name]=primer
		return proper, message
	def get_fragment(self, name):
		try:
			return self.primerdict[name][1]
		except KeyError:
			return False

	def get_primer(self, name):
		try:
			return self.primerdict[name][0]
		except KeyError:
			return False
	def get_revcomp(self, name):
		complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'R':'Y','Y':'R', 'S':'S',
					  'W':'W', 'K':'M','B':'T','V':'B','D':'H','H':'D'}
		reverse_complement = "".join(complement.get(base, base) for base in reversed(self.primerdict[name][0]))
		return reverse_complement


#primertje = Primer()
#print(primertje.get_primer("ITS4"))
#print(primertje.get_revcomp("ITS4"))
#print(primertje.get_fragment("ITS3"))
