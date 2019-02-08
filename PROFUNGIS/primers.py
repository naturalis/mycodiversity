#!/usr/bin/python3
"""
Class Primer contains all methods needed to get the primer the user needs.
"""
class Primer:
	def __init__(self):
		"""\
		Init method gathers the primer datafile and uses it in do_primerdict
		"""
		self.primerfile = "./deps/primer.data"
		self.primerdict = self.do_primerdict(self.primerfile)

	@staticmethod
	def do_primerdict(primerfile):
		"""/
		Method loads the data of the primer datafileself.
		This file is structured as name:sequence:amplicon
		The method reads the datafile and condenses the data to a dictionary
		with the primer name as key and sequence and amplicon as values
		"""
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
		"""/
		Method checks if primer input by user is in the right formatself.
		It checks all characters in the primer for invalid characters
		It checks if the name of the primer is already added before
		It checks if the sequence has already been added before under a
		different name
		It checks the length of the primer (>6)
		It checks if the provided amplicon is valid (ITS1 or ITS2)
		If none are violated, the method returns True. It returns False
		for the first violation.
		"""
		goodchars = ['A','C','G','T','R','Y','W','D','S','B','V','K','H','M']
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
		"""/
		Method adds primers to the primer datafile. It gets a verdict
		on validity from the check_input method. If check_input returned True
		the primer is written to the primer datafile and to the dictionary
		containing the datafile on runtime.
		"""
		proper, message = self.check_input(name, primer, fragment)
		if proper:
			with open(self.primerfile,'a') as writefile:
				writefile.write("%s:%s:%s\n"%(name, primer,fragment.upper()))
			self.primerdict[name]=primer
		return proper, message
	def get_fragment(self, name):
		"""/
		Method returns the amplicon of a primer based on the primer name.
		If the primer name does not exist, it returns False
		"""
		try:
			return self.primerdict[name][1]
		except KeyError:
			return False

	def get_primer(self, name):
		"""/
		Method returns sequence of a primer based on primer name. If the name
		does not exist, the method returns False.
		"""
		try:
			return self.primerdict[name][0]
		except KeyError:
			return False

	def get_revcomp(self, name):
		"""/
		Method returns reverse complement of a sequence based on the complement
		dictionary containing all valid nucleotides and wildcards
		"""
		complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'R':'Y','Y':'R',
					  'S':'S', 'W':'W', 'K':'M','B':'T','V':'B','D':'H','H':'D'}
		reverse_complement = "".join(complement.get(base, base) for base in reversed(self.primerdict[name][0]))
		return reverse_complement
