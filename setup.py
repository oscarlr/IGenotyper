from setuptools import setup

setup(
    name='IGenotyper',
    version='1.1',
    packages=['IGenotyper', 'IGenotyper.commands','IGenotyper.commands.msa'],
    url='oscarlr.github.io',
    license='',
    author='Oscar Rodriguez',
    author_email='oscar.rodriguez@icahn.mssm.edu',
    description='',
    package_data={'IGenotyper': ['data/*',
                                 'scripts/*']},
    entry_points = {
        'console_scripts': ['IG = IGenotyper.main:main']
    }
)
