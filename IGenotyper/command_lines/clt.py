#!/usr/bin/env python3
import os
import subprocess


def non_emptyfile(path):
    return os.path.isfile(path) and os.path.getsize(path) > 0


class CommandLine:
    """Run an external pipeline command once and fail clearly on errors."""

    def __init__(self, files, cpu, sample):
        self.files = files
        self.cpu = cpu
        self.sample = sample

    def run_command(self, command, output_file):
        print("-----------------")
        print("Checking %s" % output_file)
        if not non_emptyfile(output_file):
            print("\tRunning command... \n%s" % command)
            subprocess.check_call(command, shell=True, executable="/bin/bash")
            if not non_emptyfile(output_file):
                raise RuntimeError(
                    "Command completed without creating expected output: %s"
                    % output_file
                )
        print("-----------------")
