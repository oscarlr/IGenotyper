#!/usr/bin/env python3
import json
import os
import re
import subprocess
import tempfile
from pathlib import Path
from shlex import quote


def non_emptyfile(path):
    return os.path.isfile(path) and os.path.getsize(path) > 0


def signatures(paths):
    return {str(p): [os.stat(p).st_size, os.stat(p).st_mtime_ns, os.stat(p).st_ctime_ns] for p in paths}


def stage_command(command, replacements):
    """Replace complete shell path words, including quoted paths, in our commands.

    Never replace substrings (e.g. x.bam inside x.bam.bai). All output paths must
    be declared; derived sidecars share the staged parent's basename/directory.
    """
    words = {}
    for original, staged in replacements.items():
        for spelling in (quote(original), '"' + original + '"', "'" + original + "'"):
            words[spelling] = quote(staged)
    pattern = r'(?<![\w./-])(?:' + '|'.join(re.escape(w) for w in sorted(words, key=len, reverse=True)) + r')(?![\w./-])'
    return re.sub(pattern, lambda match: words[match.group()], command)


class CommandLine:
    """Stage outputs, validate successful commands, publish files then a receipt.

    Legacy outputs have no proof of completion and are rebuilt. Failed commands
    leave prior files untouched. Multi-output commands get one receipt, published
    only after every file has been replaced; an interrupted publication is untrusted.
    Concurrent writers to the same output directory are not supported.
    """

    def __init__(self, files, cpu, sample):
        self.files = files
        self.cpu = cpu
        self.sample = sample

    def run_command(self, command, output_file, validator=None, inputs=(), action=None):
        outputs = [str(output_file)] if isinstance(output_file, (str, os.PathLike)) else list(map(str, output_file))
        receipt = Path(outputs[0] + '.success.json')
        expected = {'schema': 2, 'command': command, 'inputs': signatures(inputs)}

        def validate(paths):
            for path in paths:
                if not non_emptyfile(path):
                    raise RuntimeError('Missing or empty command output: %s' % path)
            for path in paths:
                if path.endswith('.bam'):
                    import pysam
                    pysam.quickcheck(path)
                    if path + '.bai' in paths:
                        with pysam.AlignmentFile(path, 'rb', check_sq=False) as bam:
                            bam.check_index()
            if validator:
                validator(paths)

        try:
            state = json.loads(receipt.read_text())
            if state == dict(expected, outputs=signatures(outputs)):
                validate(outputs)
                return
        except (OSError, ValueError, RuntimeError):
            pass

        receipt.unlink(missing_ok=True)
        with tempfile.TemporaryDirectory(prefix='.igenotyper-', dir=receipt.parent) as work:
            parents = {}
            staged = []
            for path in outputs:
                parent = str(Path(path).parent.resolve())
                directory = Path(work).resolve() / str(parents.setdefault(parent, len(parents)))
                directory.mkdir(exist_ok=True)
                staged.append(str(directory / Path(path).name))
            staged_command = stage_command(command, dict(zip(outputs, staged)))
            print('Running command:\n%s' % staged_command)
            if action is None:
                subprocess.check_call('set -euo pipefail; ' + staged_command, shell=True, executable='/bin/bash')
            else:
                action(staged)
            validate(staged)
            # Inputs must not have changed while the command was running.
            if signatures(inputs) != expected['inputs']:
                raise RuntimeError('Command inputs changed during execution')
            for temporary, final in zip(staged, outputs):
                os.replace(temporary, final)
            temporary_receipt = Path(work) / 'success.json'
            temporary_receipt.write_text(json.dumps(dict(expected, outputs=signatures(outputs))))
            os.replace(temporary_receipt, receipt)
