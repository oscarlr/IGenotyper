#!/bin/env python
import distutils.spawn

def check_tools():
    tools = [
        "ccs",
        "bam2fastq",
        "blasr",
        "samtools"
    ]

    tools = []
    
    missing_tools = []

    for tool in tools:
        tool_path = distutils.spawn.find_executable(tool)
        if tool_path == None:
            missing_tools.append(tool)

    return missing_tools
