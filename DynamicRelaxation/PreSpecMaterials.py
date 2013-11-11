class Materials():
    
    # Materials properties are in N and mm
    # E in N/mm2
    # I in mm^4
    # A in mm^2
    # D in mm
    
    def __init__(self):    
        self.E = -1
        self.A = -1
        self.I = -1
        self.C = -1
        self.D = -1
        self.name = -1

    def steel(self, NAME = 'STEEL', E = 500, I = 0.78539, A = 3.1415, C = 1, D = 1): 
        self.E = E
        self.A = A
        self.I = I
        self.C = C
        self.D = D
        self.name = NAME
        return [self.name, self.E, self.A, self.I, self.C, self.D]

    def gfrp(self, NAME = 'GFRP', E = 40000, I = 0.78539, A = 3.1415, C = 1, D = 1): 
        self.E = E
        self.A = A
        self.I = I
        self.C = C
        self.D = D
        self.name = NAME
        return [self.name, self.E, self.A, self.I, self.C, self.D]

    def rattan(self, NAME = 'RATTAN', E = 100, I = 0.78539, A = 3.1415, C = 1, D = 1): 
        self.E = E
        self.A = A
        self.I = I
        self.C = C
        self.D = D
        self.name = NAME
        return [self.name, self.E, self.A, self.I, self.C, self.D]
    
    def pvc(self, NAME = 'PVC', E = 3000, I = 0.78539, A = 3.1415, C = 1, D = 1):
        self.E = E
        self.A = A
        self.I = I
        self.C = C
        self.D = D
        self.name = NAME
        return [self.name, self.E, self.A, self.I, self.C, self.D]
        
    def membrane(self, NAME = 'MEMBRANE', E = 10, I = 0.78539, A = 3.1415, C = 0.5, D = 1): 
        self.E = E
        self.A = A
        self.I = I
        self.C = C
        self.D = D
        self.name = NAME
        return [self.name, self.E, self.A, self.I, self.C, self.D]
    
    def wood(self, NAME = 'WOOD', E = 7000, I = 0.78539, A = 3.1415, C = 1, D = 1): 
        self.E = E
        self.A = A
        self.I = I
        self.C = C
        self.D = D
        self.name = NAME
        return [self.name, self.E, self.A, self.I, self.C, self.D]