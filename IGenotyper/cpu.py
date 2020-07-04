#!/bin/env python

class CpuManager:
    def __init__(self,threads,mem,cluster,queue,walltime):
        self.threads = threads
        self.mem = mem
        self.cluster = cluster
        self.queue = queue
        self.walltime = walltime
