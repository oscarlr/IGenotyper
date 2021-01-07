#!/bin/env python
import os
#from lsf.lsf import Lsf
#from IGenotyper.helper import non_emptyfile

from IGenotyper.command_lines.clt import CommandLine

class PlotTools(CommandLine):
    def __init__(self,files,cpu,sample):
        CommandLine.__init__(self,files,cpu,sample)

    def run_pygenometracks(self,config,plotfn):
        args = [config,plotfn]
        command = ("CONDA_BASE=$(conda info --base) \n"
                   "source ${CONDA_BASE}/etc/profile.d/conda.sh \n"
                   "conda activate pygenometracks \n"
                   "pyGenomeTracks --tracks %s --region igh:1-1193129 -o %s" % tuple(args))
        print command
        self.run_command(command,plotfn)

