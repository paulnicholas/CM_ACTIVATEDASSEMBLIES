class Materials():

    def __init__(self):    
        self.E = -1
        self.A = -1
        self.I = -1
        self.C = -1
        self.D = -1
        self.name = -1

    def steel(self, NAME = 'STEEL', E = 2500, I = 1000, A = 1000, C = 1, D = 1): 
        self.E = E
        self.A = A
        self.I = I
        self.C = C
        self.D = D
        self.name = NAME
        return [self.name, self.E, self.A, self.I, self.C, self.D]

    def gfrp(self, NAME = 'GFRP', E = 2100, I = 1000, A = 1000, C = 1, D = 1): 
        self.E = E
        self.A = A
        self.I = I
        self.C = C
        self.D = D
        self.name = NAME
        return [self.name, self.E, self.A, self.I, self.C, self.D]

    def rattan(self, NAME = 'RATTAN', E = 2100, I = 1000, A = 1000, C = 1, D = 1): 
        self.E = E
        self.A = A
        self.I = I
        self.C = C
        self.D = D
        self.name = NAME
        return [self.name, self.E, self.A, self.I, self.C, self.D]
    
    def pvc(self, NAME = 'PVC', E = 2100, I = 1000, A = 1000, C = 1, D = 1):
        self.E = E
        self.A = A
        self.I = I
        self.C = C
        self.D = D
        self.name = NAME
        return [self.name, self.E, self.A, self.I, self.C, self.D]
        
    def membrane(self, NAME = 'MEMBRANE', E = 2100, I = 1000, A = 1000, C = 1.9, D = 1): 
        self.E = E
        self.A = A
        self.I = I
        self.C = C
        self.D = D
        self.name = NAME
        return [self.name, self.E, self.A, self.I, self.C, self.D]

    def wood(self, NAME = 'WOOD', E = 2100, I = 1000, A = 1000, C = 1.9, D = 1): 
        self.E = E
        self.A = A
        self.I = I
        self.C = C
        self.D = D
        self.name = NAME
        return [self.name, self.E, self.A, self.I, self.C, self.D]