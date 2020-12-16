#!/bin/env python
import sys
import argparse
import importlib
from IGenotyper.check import check_tools

def main():
    commands = [
        "phase",
        "assembly",
        "detect"
    ]

    if len(sys.argv) < 2:
        sys.exit("Please run one of the following commands: \n%s" % "\n".join(commands))

    command_name = sys.argv[1]

    if command_name not in commands:
        sys.exit("Please run one of the following commands: \n%s" % "\n".join(commands))

    missing_tools = check_tools()
    if len(missing_tools) != 0:
        sys.exit("Install tools: %s" % ",".join(missing_tools))

    parser = argparse.ArgumentParser(description='Process IGH capture data')
    subparsers = parser.add_subparsers()

    command = importlib.import_module('.%s' % command_name, "IGenotyper.commands")
    subparser = subparsers.add_parser(command_name)
    command.add_arguments(subparser)

    args = parser.parse_args(sys.argv[1:])
    command.main(args)

if __name__ == "__main__":
    main()
