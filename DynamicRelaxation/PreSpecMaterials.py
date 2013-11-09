class Materials():

    def __init__(self):    
	self.E = 0
	self.I = 0
	self.A = 0
	self.C = 0

    def steel(self, E = 2500, I = 1000, A = 1000, C = 0): 
	self.E = E
	self.I = I
	self.A = A
	self.C = C
	
	return ['STEEL', E, A, I, C]

    def gfrp(self, E = 2100, I = 1000, A = 1000, C = 0): 
	self.E = E
	self.I = I
	self.A = A
	self.C = C

	return ['GFRP', E, A, I, C]

    def membrane(self, E = 2100, I = 1000, A = 1000, C = 0.9): 
	self.E = E
	self.I = I
	self.A = A
	self.C = C

	return ['MEMBRANE', E, A, I, C]