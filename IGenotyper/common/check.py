#!/usr/bin/env python3
import shutil

def check_tools(command_name):
    tools_by_command = {
        "phase": [
            "bamCoverage",
            "minimap2",
            "Rscript",
            "samtools",
            "whatshap",
        ],
        "assembly": ["canu", "minimap2", "samtools"],
        "detect": ["kalign"],
        "alleles": [],
    }
    tools = tools_by_command.get(command_name, [])
    return [tool for tool in tools if shutil.which(tool) is None]
