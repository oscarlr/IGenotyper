#!/bin/env python
import os
#from lsf.lsf import Lsf
#from IGenotyper.helper import non_emptyfile

from IGenotyper.command_lines.clt import CommandLine

class VariantTools(CommandLine):
    def __init__(self,files,cpu,sample):
        CommandLine.__init__(self,files,cpu,sample)

    def run_kalign(self,fastafn,clufn):
        command = "kalign -i %s -f clu -o %s" % (fastafn,clufn)
        self.run_command(command,clufn)

