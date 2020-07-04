#!/bin/env python
import os

def create_folders(folders):
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder)

def non_emptyfile(checkfile):
    return os.path.isfile(checkfile) and os.path.getsize(checkfile) > 0