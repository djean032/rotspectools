import subprocess
import sys

from pathlib import Path
from subprocess import PIPE, Popen


base_path = Path(sys.argv[0]).resolve().parent
spfit = str(base_path) + '\spfit.exe'
piform = str(base_path) + '\piform.exe'
base_file = 'cyanomethcycloprop_gs'
res_file = 'cyanomethcycloprop_gs.res'
subprocess.run([spfit], shell=True)
print(base_file)