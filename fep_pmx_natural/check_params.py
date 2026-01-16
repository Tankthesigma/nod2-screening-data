#!/usr/bin/env python
"""Check global parameters in alchemical system."""
from openmm import XmlSerializer

with open('C:/Users/vasud/nod2-screening-data/fep_pmx_natural/wt_complex/alchemical_system.xml', 'r') as f:
    system = XmlSerializer.deserialize(f.read())

print('Forces and their global parameters:')
for i in range(system.getNumForces()):
    force = system.getForce(i)
    fname = force.__class__.__name__
    if hasattr(force, 'getNumGlobalParameters'):
        n = force.getNumGlobalParameters()
        if n > 0:
            params = [force.getGlobalParameterName(j) for j in range(n)]
            print(f'  {i}: {fname} - {params}')
