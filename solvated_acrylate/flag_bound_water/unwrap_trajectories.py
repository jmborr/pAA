from __future__ import print_function

from tqdm import tqdm
import os
from tempfile import mkstemp

script_template = """package require pbctools
mol new _PDB_
animate read dcd _WRAPPED_ beg 0 end -1 waitfor all
pbc unwrap -all
animate write dcd _UNWRAPPED_ beg 1 end -1 waitfor all
"""

datadir = '/data/nfs/camm/users/ggz'
pdb_file = 'acrylate.pdb'

script_handle, script_file = mkstemp()
print('script_file =', script_file)

wrappeds = 'PAA.T270K.2NS.dcd  PAA.T280K.2NS.dcd  PAA.T290K.2NS.dcd  PAA.T300K.2NS.dcd  PAA.T310K.2NS.dcd'.split()
for wrapped in wrappeds:
    print(wrapped)
    script = script_template.replace('_PDB_', pdb_file)
    script = script.replace('_WRAPPED_', os.path.join(datadir,wrapped))
    unwrapped, extension = os.path.splitext(wrapped)
    unwrapped += '_unwrapped.dcd'
    script = script.replace('_UNWRAPPED_', unwrapped)
    open(script_file, 'w').write(script)
    os.system('sleep 9s; vmd -dispdev text -e {} &'.format(script_file))

os.system('/bin/rm {}'.format(script_file))
