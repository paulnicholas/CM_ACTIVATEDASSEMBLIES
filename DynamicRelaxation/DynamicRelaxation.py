# FileName 'DynamicRelaxation.py'

import math

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class Vector3D(object):
    
    def __init__(self,*args):
        if len(args) == 0:
            self.init(0,0,0)
        elif len(args) == 1:
            v = args[0]
            self.init(v.x(),v.y(),v.z())
        elif len(args) == 3:
            self.init(args[0],args[1],args[2])
        
    def init(self,x,y,z):
        self._x = x
        self._y = y
        self._z = z
        
    def __str__(self): return '['+str(self._x)+', '+str(self._y)+', '+str(self._z)+']'
    def x(self): return self._x
    def y(self): return self._y
    def z(self): return self._z
        
    def set(self,*args):
        if len(args)==1:
            v = args[0]
            set(v.x(),v.y(),v.z()); return self 
        elif len(args)==3: 
            self._x = args[0]; self._y = args[1]; self._z = args[2]; return self
        
    def setX(self,x):self._x = x; return self 
    def setY(self,y):self._y = y; return self 
    def setZ(self,z):self._z = z; return self 
    
    
    def __add__(self, p ): return Vector3D(self.x() + p.x(),self.y() + p.y(),self.z() + p.z())
    def addBy(self, a,b,c ): return Vector3D(self._x + a, self._y + b, self._z + c)
    
    def __sub__(self, p ): return Vector3D(self.x() - p.x(),self.y() - p.y(),self.z() - p.z())
    def subtractBy(self, a, b, c): return Vector3D(self._x - a, self._y - b, self._z - c)
    
    def __mul__(self,  v ): return Vector3D( self._x * v.x(),self._y * v.y(), self._z * v.z())
    def multiplyBy(self,  f ): return Vector3D(self._x * f, self._y * f, self._z * f)
    
    def __div__(self,  v ): return Vector3D( self._x / v.x(),self._y / v.y(), self._z / v.z())
    def divideBy(self,  f ): return Vector3D(self._x / f, self._y / f, self._z/ f)
    
    def distanceTo(self, p ): return math.sqrt( self.distanceSquaredTo( p ) )
    def distanceSquaredTo(self, p ): d = self - p; return d.dot(d)
    
    
    def dot(self, p ): return self._x * p.x() + self._y * p.y() + self._z * p.z()
    def length(self): return math.sqrt( self.lengthSquared() )
    def lengthSquared(self): return self.dot(self) 
    
    def normalize(self): return Vector3D(self).divideBy(self.length())
      
    def clear(self):self._x = 0; self._y = 0; self._z = 0; return self
    
    
    def __str__(self): return "(" + str(self._x) + ", " + str(self._y) + ", " + str(self._z) + ")" 
    
    def __neg__(self): return Vector3D(-self.x(),-self.y(),-self.z())
    def cross(self, p ): return Vector3D(     self.y() * p.z() - self.z() * p.y(), self.z() * p.x() - self.x() * p.z(), self.x() * p.y() - self.y() * p.x() )
    #def cross(self, p ): return Vector3D(     self.y() * p.z() - self.z() * p.y(), self.x() * p.z() - self.z() * p.x(), self.x() * p.y() - self.y() * p.x() )
    
    def isZero(self): return self._x == 0 and self._y == 0 and self._z == 0;
    
    def distanceSegmentSq(self, p, q) :
        v = q - p
        w = self - p
        
        c1 = w.dot(v)
        if ( c1 <= 0 ):
            return self.distanceSquaredTo(p);
        
        c2 = v.dot(v)
        if ( c2 <= c1 ):
            return self.distanceSquaredTo(q);
            
        b = c1 / c2
        Pb =  p + v.multiplyBy(b);
        return self.distanceSquaredTo(Pb);
    
    
    
    def distanceLineSq(self, p, q) :
        v = q - p
        w = self - p
        c1 = w.dot(v);
        c2 = v.dot(v);
        b = c1 / c2;
        Pb =  p + v.multiplyBy(b);
        return self.distanceSquaredTo(Pb);
    
    
    #     * project on plane     * assume normalized direction     * @param dir
    def project(self, dir) : t = -self.dot(dir); return Vector3D(self + dir.multiplyBy(t))
    
    
    #     * project on vector     * assume normalized direction     * @param dir
    def projectDir(self, dir) : t = self.dot(dir); return Vector3D(dir.multiplyBy(t))
    
class Particle(object):
    
    def __init__(self,  m ):
        self.position = Vector3D();
        self.velocity = Vector3D();
        self.force = Vector3D();
        self.mass = Vector3D(m,m,m);
        self.fixed = False;
        self.constraint = Vector3D();
        self.constrained = False;
    
    def __str__(self):
        return 'p: '+str(self.position)+', v: '+str(self.velocity)+ ', fixed: '+str(self.fixed)+', constrained: '+str(self.constrained)
    
    def distanceTo(self, p ):
      return self.position.distanceTo( p.position )
    
    def makeFixed(self):
        self.fixed = True;
        self.velocity.clear()
    
    def isFixed(self):
        return self.fixed
        
    def isFree(self):
        return not self.fixed
        
    def isConstrained(self):
        return self.constrained
        
    def makeFree(self):
        self.fixed = False
    
    def massAverage(self):
        return (self.mass.x() + self.mass.y() + self.mass.z())/3
        
    def setMass(self, m ):
        self.mass = Vector3D(m,m,m)
        
    def constrained(self, nx, ny, nz) :
        self.constraint = Vector3D(nx,ny,nz).normalize()
        self.constrained = True
    
    def getConstraint(self) :
        return self.constraint
        
    def constraint(self):
        self.force.project(constraint)
        
    def update(self):
        pass
        
    def reset(self):
        self.position.clear()
        self.velocity.clear()
        self.force.clear()
        #self.mass = Vector3D(1,1,1)

class ParticleSystem(object):

    RUNGE_KUTTA = 0
    MODIFIED_EULER = 1
    VERLET = 2
    EULER = 3
    
    #DEFAULT_GRAVITY = Vector3D()
    #DEFAULT_DRAG = 0.001  
    
    def __init__(self, gravity = Vector3D(), drag = 0.0001):
        self.init()
        self.gravity = gravity
        self.drag = drag
        
    def __str__(self):
        return "integrator   : "+str(self.integrator)+"\n"+"constraints  : "+str(self.hasConstrainedParticles)+"\n"+"gravity  : "+str(self.gravity)+"\n"+"drag     : "+str(self.drag)
        
    def init(self):
        self.integrator = VerletIntegrator( self )
        
        self.hasConstrainedParticles = False;
        self.prevKinetic = 0
        
        self.particles = []
        self.springs = []
        self.attractions = []
        self.bendings = []
        self.elastics = []
        self.cables = []
        self.loads = []
        self.customForces = []
        
    def setIntegrator(self, integrator ) :
        if integrator == ParticleSystem.RUNGE_KUTTA:
            self.integrator =  RungeKuttaIntegrator( self )
        elif integrator == ParticleSystem.MODIFIED_EULER:
            self.integrator =  ModifiedEulerIntegrator( self )
        elif integrator == ParticleSystem.VERLET:
            self.integrator =  VerletIntegrator( self )
        elif integrator == ParticleSystem.EULER:
            self.integrator =  EulerIntegrator( self )
        
        
    def setGravity(self, x, y, z ) :
        self.gravity.set( x, y, z )
        
    def setGravity(self, g ) :
        self.gravity.set( 0, g, 0 )
        
    def setDrag( self, d ) :
        drag = d
    
    def step( self,  t=1.0 ):  
        self.integrator.step( t )
        
    def stepRelax( self, t=1.0 ) :
        self.integrator.step( t )
        
        tmp = self.computeKinetic()
        
        if (tmp < self.prevKinetic - 10e-04) :
            
            for p in self.particles:
                p.velocity.clear()
            self.prevKinetic = 0
            
        else: self.prevKinetic = tmp
    
    def computeKinetic(self) :
        kinetic = 0;
        
        for p in self.particles:
            kinetic += .5 * p.mass.dot(p.velocity * p.velocity)
        return kinetic
    
    
    
    def makeParticle(self, *args) :
        if len(args)==0:        
            return self.makeParticle( 1.0, 0, 0, 0 );
        #def makeParticle(self, v) :
        elif len(args)==1:
            v = args[0]
            p = Particle( 1.0 )
            p.position.set( v )
            self.particles.append( p )
            return p
        #def makeParticle( self, x, y, z ) :
        elif len(args)==3:
            p = Particle( 1.0 )
            p.position.set( args[0], args[1], args[2] )
            self.particles.append( p )
            return p
        #def makeParticle( self, mass, x, y, z ) :
        elif len(args)==4:
            mass = args[0]; 
            p = Particle( mass )
            p.position.set( args[1], args[2], args[3] )
            self.particles.append( p )
            return p
        
    def makeSpring(self, a, b, k, l0, damp = 0 ) :
        s = Spring( a, b, k, l0, damp )
        self.springs.append( s )
        return s
    
    def makeAttraction(self, a, b, k, minDistance ) :
        m = Attraction( a, b, k, minDistance )
        self.attractions.append( m )
        return m
    
    def makeBending(self, a, b, c, E, I, S ) :
        m = Bending( a, b, c, E , I , S)
        self.bendings.append( m )
        return m
    
    def makeElastic( self, a,  b,  E,  A,  l0,  t0 ) :
        m = Elastic( a, b, E , A , l0, t0)
        self.elastic.append( m )
        return m
    
    def makeCable( self, a,  b,  E,  A,  l0,  t0 ) :
        m = Cable( a, b, E , A , l0, t0)
        self.cables.append( m )
        return m
    
    def makeLoad( self, a, v, f) :
        m = Load( a, v, f)
        self.loads.append( m )
        return m
    
    def addCustomForce(self,  f ):  
        self.customForces.append( f )
    
    def clear() :
        particles = []
        springs = []
        attractions = []
        bendings = []
        elastics = []
        cables = []
        loads = []
    
    def applyForces(self) :
        
        if not self.gravity.isZero() : 
            for p in self.particles:  p.force = p.force + self.gravity
            
        # apply drag
        for p in self.particles:  
            p.force = p.force + p.velocity.multiplyBy(-self.drag) 
        
        for  f in self.springs: f.apply()
        for  f in self.attractions: f.apply()
        for  f in self.bendings: f.apply()
        for  f in self.elastics: f.apply()
        for  f in self.cables: f.apply()
        for  f in self.loads: f.apply()
        for  f in self.customForces: f.apply()
      
        if (self.hasConstrainedParticles) :
            for p in self.particles: 
                p.update()
                p.constraint()
    
    def clearForces(self) :
        for p in self.particles: p.force.clear()
    
    def numberOfParticles(self):    return len(self.particles)  
    def numberOfSprings(self):        return len(self.springs)
    def numberOfAttractions(self):  return len(self.attractions)
    def numberOfBendings(self):        return len(self.bendings)
    def numberOfElastic(self):      return len(self.elastics)
    def numberOfCables(self):        return len(self.cables)
    def numberOfLoads(self):        return len(self.loads)
    def numberOfCustomForces(self):    return len(self.customForces)
    
    def getParticle( self, i ):        return self.particles[i]  
    def getSpring( self, i ):        return self.springs[i]  
    def getAttraction( self, i ):   return self.attractions[i]
    def getBending( self, i ):        return self.bendings[i]
    def getElastic( self, i ):        return self.elastics[i]
    def getCable( self, i ):        return self.cables[i]
    def getLoad( self, i ):            return self.loads[i]
    def getCustomForce( self, i ):    return self.customForces[i]  
    
    def removeParticle( self, i ):      return self.particles.remove( i )  
    def removeSpring( self, i ):        return self.springs.remove( i )
    def removeAttraction( self, i  ):    return self.attractions.remove( i )
    def removeBending( self, i  ):        return self.bendings.remove( i )
    def removeElastic( self, i  ):        return self.elastics.remove( i )
    def removeCable( self, i  ):        return self.cables.remove( i )
    def removeLoad( self, i  ):            return self.loads.remove( i )
    def removeCustomForce( self, i ):    return self.customForces.remove( i )  
    
    def removeParticle( self, p ):        self.particles.pop( p ); return self.particles
    def removeSpring( self, a ) :        self.springs.pop( a ) ; return self.springs
    def removeAttraction( self, s ):    self.attractions.pop( s ); return self.attractions
    def removeBending( self, s ):        self.bendings.pop( s ); return self.bendings
    def removeElastic( self, s ):        self.elastics.pop( s ); return self.elastics
    def removeCable( self, s ):            self.cables.pop( s ); return self.cables
    def removeLoad( self, s ):            self.loads.pop( s ); return self.loads
    def removeCustomForce( self, f ):   self.customForces.pop( f ); return self.customForces

#@@@@@@@@@@@@@@@@

class Integrator:
    
    def step(self, t):
        pass

class VerletIntegrator(Integrator):
    
    # number of iterations for constraints
    NUM_ITERATIONS = 10
    EPSILON = 1e-06
    
    def __init__(self,particleSystem):
        self._ps = particleSystem
        
    def step(self, t ):
        
        self._ps.clearForces()
        self._ps.applyForces()
        
        halftt = 0.5/(t*t)
        halft = 0.5/t
        
        for p in self._ps.particles:
            if p.isFree():
                
                ax = p.force.x() / p.mass.x()
                ay = p.force.y() / p.mass.y()
                az = p.force.z() / p.mass.z()
                
                # update position
                p.position = p.position + p.velocity.divideBy(t)
                p.position = p.position + Vector3D( ax*halftt, ay*halftt, az*halftt )
                
                # half step with previous acceleration
                p.velocity = p.velocity + Vector3D( ax*halft, ay*halft, az*halft )
                
                # compute forces with new position
                self._ps.clearForces()
                self._ps.applyForces()
                
                # compute acceleration in new position
                ax = p.force.x() / p.mass.x()
                ay = p.force.y() / p.mass.y()
                az = p.force.z() / p.mass.z()
                
                # half step with new acceleration
                p.velocity = p.velocity + Vector3D( ax*halft, ay*halft, az*halft )
                
        self.satisfyConstraints()
        
    def satisfyConstraints(self) :
        
        for j in range(VerletIntegrator.NUM_ITERATIONS):
            
            for f in self._ps.springs:
                
                if (f.k == 0) :
                    # apply constraints
                    
                    # measure distance
                    a = f.getOneEnd().position
                    b = f.getTheOtherEnd().position
                    
                    d1 = a - b
                    d2 = d1.length() # potential to optimize with taylor expansion
                    diff = (d2-f.l0())/d2;
                    d1.multiplyBy(.5 * diff);
                    f.getOneEnd().position.set(a - d1);
                    f.getTheOtherEnd().position.set(a + d1);
                    
                    if diff<VerletIntegrator.EPSILON:
                        continue

class EulerIntegrator(Integrator):
    
    def __init__(self,particleSystem):
        self._ps = particleSystem
    
    def step(self, t ):
        
        self._ps.clearForces();
        self._ps.applyForces();
        
        for p in self._ps.particles:
            if p.isFree():
                p.velocity = p.velocity + Vector3D( p.force.x()/(p.mass.x()*t), p.force.y()/(p.mass.y()*t), p.force.z()/(p.mass.z()*t) )
                p.position = p.position + p.velocity.divideBy(t)

#@@@@@@@@@@@@@@@@

class Force:
    
    def __init__(self):
        self.on = True;
      
    def apply(self):
        pass
          
    def turnOff(self):
        self.on = False
    
    def turnOff(self):
          self.on = True
          
    def isOn(self):
        return self.on
    def isOff(self):
        return not self.on

class Spring(Force):
    
    def __init__( self, A,  B,  k, l0, damp ):
        self.k = k
        self.damping = damp
        self.l0 = l0
        self.a = A
        self.b = B
        self.on = True
        
        
    def getOneEnd(self) : return self.a
    def getTheOtherEnd(self)  : return self.b
    
    def currentLength(self) : return self.a.position.distanceTo( self.b.position )
    
    def setStrength( self, k )  : self.k = k
    
    def setDamping( self, d )  : self.damping = d
    
    def setRestLength( self, l )  : self.l0 = l
    
    def setA( self, p )  : self.a = p
    def setB( self, p )  : self.b = p
    
    def apply(self):
        
        if ( self.on and ( self.a.isFree() or self.b.isFree() ) ):
            
            a2b = self.a.position - self.b.position
            
            a2bDistance = math.sqrt( a2b.dot(a2b) )
        
            if ( a2bDistance == 0 ):
                a2b = Vector3D()
            else :
                a2b = a2b.multiplyBy(a2bDistance)
                
            # spring force is proportional to how much it stretched
            
            springForce = -( a2bDistance - self.l0 ) * self.k
            
            
            # want velocity along line b/w a & b, damping force is proportional to this
            
            Va2b = self.a.velocity - self.b.velocity
                    
            dampingForce = - self.damping * a2b.dot(Va2b)
            
            
            # forceB is same as forceA in opposite direction
            
            r = springForce + dampingForce
            
            a2b = a2b.multiplyBy(r)
            
            if ( self.a.isFree() ):
                self.a.force = self.a.force + a2b
            if (self. b.isFree() ):
                self.b.force = self.b.force - a2b

class Elastic(Force) :
    # E module of elastic material   
    # A section area of elastic material
    # l0 initial length 
    # t0 prestress 
    
    
    # Particle a, b;
    
    def __init__ (self,  _a,  _b,  _E,  _A,  _l0,  _t0 = 0 ) :
        self.E = _E
        self.A = _A
        self.l0 = _l0
        self.t0 = _t0
        self.a = _a
        self.b = _b
        self.on = true
        self.fixTension = False
        
    def getOneEnd(self)  : return self.a
    def getTheOtherEnd(self)  :      return self.b
    
    def setPartA(self, p )  :      self.a = p
    def setPartB( self, p )  :      self.b = p
    
    def currentLength(self) :      return self.a.position.distanceTo( self.b.position )
    
    def getRestLength(self)  :  return self.l0
    def setRestLength( self, _l ) :  self.l0 = _l 
    
    def getPrestress(self)  :  return self.t0
    def setPrestress( self, _t ) : self.t0 = _t 
    
    
    def getE(self)  :  return self.E
    def setE(self, _E)  :  self.E = _E
    
    def getA(self)  :  return self.A
    def setA(self, _A)  : self.A = _A
    
    
    def apply(self):
        
        if ( self.on and ( self.a.isFree() or self.b.isFree() ) ) :
        
            a2b = self.a.position - self.b.position            
            a2bDistance = math.sqrt( a2b.dot(a2b) )
            
            if ( a2bDistance == 0 ):
                a2b = Vector3D()
            else :
                a2bX .divideBy( a2bDistance)
            
            # elastic force is proportional to how much it stretched 
            # t = EA/l0 * (l-l0) + t0  
            elasticForce = 0
            
            # if tension not fixed
            if (not fixTension):
                elasticForce = -( a2bDistance - self.l0 ) * self.E*self.A/self.l0 + self.t0
            # if tension fixed prestress only
            else: elasticForce = t0
            
            
            a2b = a2b.multiplyBy(elasticForce)
            
            if ( self.a.isFree() ):
                self.a.force = self.a.force + a2b
            
            # forceB is same as forceA in opposite direction
            if ( self.b.isFree() ):
                self.b.force = self.b.force - a2b

class Cable(Elastic) :
    
    def __init__( self, _a,  _b,  _E,  _A,  _l0,  _t0 = 0 ):
        super(_a,_b,_E,_A,_l0,_t0)
        self.slack = False
        
    def apply():
        if ( self.on and ( self.a.isFree() or self.b.isFree() ) ):
            
            a2b = self.a.position - self.b.position
            
            
            a2bDistance = math.sqrt( a2b.dot(a2b) )
            
            # check if cable is slack
            if ( a2bDistance >= l0 ) :
                slack = False;
                
                a2b.divideBy(a2bDistance)
                
                # elastic force is proportional to how much it stretched 
                # t = EA/l0 * (l-l0) + t0  
                
                elasticForce = 0;
                
                # if tension not fixed
                if (not fixTension):
                    elasticForce = -( a2bDistance - self.l0 ) * self.E*self.A/self.l0 + self.t0
                # if tension fixed prestress only
                else: elasticForce = t0
                
                
                a2b = a2b.multiplyBy( elasticForce)
                
                if ( self.a.isFree() ):
                    self.a.force = self.a.force + a2b
                
                # forceB is same as forceA in opposite direction
                if ( self.b.isFree() ):
                    self.b.force = self.b.force - a2b
            else :
                slack = True;

class Bending(Force):
    
    def __init__(  self, a,  b,  c,  E,  I,  S )    :
        self.init(a,b,c,a.distanceTo(b),b.distanceTo(c),E,I,S);
    
    def init( self, _a,  _b,  _c,  _L0ab,  _L0bc,  _E,  _I,  _S ) :
    #def init( self, _a,  _b,  _c,  _L0ab,  _L0bc,  _E,  _S,  _I ) :
        self.a = _a
        self.b = _b
        self.c = _c
        self.on = True
        self.EI = _E*_I
        self.ES = _E*_S
        # compute rest distances
        self.L0ab = _L0ab
        self.L0bc = _L0bc 
    
    def setA(  self, p ) : self.a = p
    def setB(  self, p ) : self.b = p
    def setC(  self, p ) : self.c = p
    
    def getOneEnd( self): return self.a
    def getTheMiddle( self): return self.b
    def getTheOtherEnd( self): return self.c
    
    def setEI( self, E, I )    : self.EI = E*I
    def setES( self, E, S )    : self.ES = E*S
    
    def apply( self):
        
        if ( self.on and ( self.a.isFree() or self.b.isFree() or self.c.isFree()) ):
            
            # compute current distances
            Lab = self.a.distanceTo(self.b);
            Lbc = self.b.distanceTo(self.c);
            Lac = self.a.distanceTo(self.c);
            
            # compute direction vectors
            ab = self.b.position - self.a.position
            bc = self.c.position - self.b.position
            
            # compute tension
            Tab = ab.multiplyBy(self.ES*(1/self.L0ab - 1/Lab));
            Tbc = bc.multiplyBy(self.ES*(1/self.L0bc - 1/Lbc));
            
            # compute moment
            Fab = Vector3D()
            Fbc = Vector3D()
            
            if (self.a.position.distanceLineSq(self.b.position,self.c.position) > 10e-8) :
                Mb = ab.cross(bc).multiplyBy(2*self.EI/(Lac*Lab*Lbc));
                Fab = ab.cross(Mb).multiplyBy(1.0/(Lab*Lab));
                Fbc = bc.cross(Mb).multiplyBy(1.0/(Lbc*Lbc));
            
            # else zero
            
            # apply
            
            if ( self.a.isFree() ) :
                self.a.force = self.a.force + Tab
                self.a.force = self.a.force + Fab
            
            if ( self.b.isFree() ) :
                self.b.force = self.b.force - Tab + Tbc 
                self.b.force = self.b.force - Fab - Fbc 
            
            if ( self.c.isFree()) :
                self.c.force = self.c.force - Tbc
                self.c.force = self.c.force + Fbc

class Attraction(Force):
    
    
    def __init__( self,  _a,  _b,  _k,  _distanceMin ):
        self.a = _a
        self.b = _b
        self.k = _k
        self.on = True
        self.distanceMin = _distanceMin
        self.distanceMinSquared = _distanceMin*_distanceMin

    def setA( self, p ):
        self.a = p

    def setB( self, p ):
        self.b = p

    def getMinimumDistance(self):
        return self.distanceMin

    def setMinimumDistance( self, d ):
        self.distanceMin = d
        self.distanceMinSquared = d*d


    def setStrength( self, _k ):
        self.k = _k
        
    def  getStrength(self):
        return self.k

    def getOneEnd(self):
        return self.a

    def getTheOtherEnd(self):
        return self.b

    def apply(self):
        if ( self.on and ( self.a.isFree() or self.b.isFree() ) ):
            
            a2b = self.a.position - self.b.position
            
            a2bDistanceSquared = a2b.dot(a2b)
            
            if ( a2bDistanceSquared < self.distanceMinSquared ) :
                a2bDistanceSquared = self.distanceMinSquared;

            force = self.k * self.a.massAverage() * self.b.massAverage() / a2bDistanceSquared;

            length = math.sqrt( a2bDistanceSquared );
            
            # make unit vector
            
            a2b = a2b.dividedBy(length) 
            
            # multiply by force 
            
            a2b = a2b.multiplyBy(force) 
            
            # apply
            
            if ( self.a.isFree() ):
                self.a.force = self.a.force - a2b
            if ( self.b.isFree() ):
                self.b.force = self.b.force + a2b

class Load(Force):
    
    def __init__( self, _a, _d, _force):
        self.a = _a
        self.force = _force
        # normalize direction
        self.dir = _d.multiplyBy(1/_d.length())
        self.on = True
    
    def setParticle( self, p ):
        self.a = p
    
    
    def setForce( self, f ):
        self.force = f
    
    def getForce(self):
        return self.force
    
    def setDirection( self, v ):
        # normalize direction
        self.dir = v.multiplyBy(1/v.length())
    
    def getParticle(self):
        return self.a
    
    
    def apply(self):
        
        load = self.dir.multiplyBy(self.force) 
        
        if ( self.on and self.a.isFree() ):
                self.a.force = self.a.force + load

