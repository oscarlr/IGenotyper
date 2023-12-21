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
        self.run_command(command,plotfn)

    def rplot_gene_cov(self):
        args = [self.files.scripts,
                self.files.gene_cov,
                self.files.plot_gene_cov,
                self.files.plot_sv_gene_cov]
        command = "Rscript %s/rplot_gene_cov.R %s %s %s" % tuple(args)
        self.run_command(command,self.files.plot_sv_gene_cov)
