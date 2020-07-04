#!/bin/env python
import sys
import argparse
import importlib

def main():
    commands = [
        "phase"
    ]

    if len(sys.argv) < 2:
        sys.exit("Please run one of the following commands: \n%s" % "\n".join(commands))

    command_name = sys.argv[1]

    if command_name not in commands:
        sys.exit("Please run one of the following commands: \n%s" % "\n".join(commands))

    parser = argparse.ArgumentParser(description='Process IGH capture data')
    subparsers = parser.add_subparsers()

    command = importlib.import_module('.%s' % command_name, "IGenotyper.commands")
    subparser = subparsers.add_parser(command_name)
    command.add_arguments(subparser)

    args = parser.parse_args(sys.argv[1:])
    command.main(args)

if __name__ == "__main__":
    main()