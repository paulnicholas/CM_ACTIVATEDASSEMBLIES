
# Particle-based simulations

## Different stages of implementations and their advantages/consequences

The base is a simple spring system (3DOF) using verlet integrator for stability. Particles have mass. it includes Dynamic Relaxation (DR) (killing kinetic energy). There is a possibility to constrain particles.

 
### Stage 1: Engineering upgrade

Added cables with prestress, prestress in members, Young modulus, etc. There is a possibility to have a cable which goes through more than two points, to model membrane edge cables.

*advantages* The first part is absolutely needed, the second one is not too difficult and might be interesting.

### Stage 2: Including bending rod model 

Include the Adriaenssens bending Rod model (3DOF).

*advantages* the model is reliable and very fast.

*inconvenients* The model is still 3DOF and therefore there is no possibility to fix the orientation of the bending rod since the particles can rotate. No reaction resulting from rod torsion is modeled.

### Stage 3: New model of bending rods for torsion

According to S.Adriaenssens, there is a new model that takes into account torsion. 

*advantages* More realistic, but not known yet, she has to send the paper...

*inconvenients* It might be a simplification that still does not allow for controlling rotations at nodes

### Stage 4: Upgrading to 6DOF

To gain more control on the local moments, an upgrade of the model to 6DOF might be interesting. Another possibility would be to work with quaternions. This should remain a topic of discussion for next steps.

*advantages* larger modeling capabilities.

*inconvenients* most probably slower than the 3DOF, unless we do clever tricks with the quaternions.

