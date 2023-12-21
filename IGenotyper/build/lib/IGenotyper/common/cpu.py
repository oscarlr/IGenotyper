#!/bin/env python

class CpuManager:
    def __init__(self,threads=None,mem=None,cluster=None,queue=None,walltime=None):
        self.threads = threads
        self.mem = mem
        self.cluster = cluster
        self.queue = queue
        self.walltime = walltime
