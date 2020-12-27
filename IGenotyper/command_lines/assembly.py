#!/bin/env python
import os
from lsf.lsf import Lsf

from IGenotyper.command_lines.clt import CommandLine

class Assembly(CommandLine):
    def __init__(self,files,cpu,sample):
        CommandLine.__init__(self,files,cpu,sample)
        
    def run_assembly_scripts(self,assembly_scripts):
        if not self.cpu.cluster:
            for script in assembly_scripts:
                command = "sh %s" % script
                os.system(command)
        else:
            hpc = Lsf()
            for job in assembly_scripts:
                hpc.config(cpu=self.cpu.threads,walltime=self.cpu.walltime,
                           memory=self.cpu.mem,queue=self.cpu.queue)
                hpc.submit("%s" % job)
            hpc.wait()

