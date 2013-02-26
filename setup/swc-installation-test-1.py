#!/usr/bin/env python

"""Test script to check required Python version.

Execute this code at the command line by typing:

  python swc-installation-test-1.py

How to get a command line:

- On OSX run this with the Terminal application.

- On Windows, go to the Start menu, select 'Run' and type 'cmd'
(without the quotes) to run the 'cmd.exe' Windows Command Prompt.

- On Linux, either use your login shell directly, or run one of a
  number of graphical terminals (e.g. 'xterm', 'gnome-terminal', ...).

For some screen shots, see:

  http://software-carpentry.org/setup/terminal.html

Run the script and follow the instructions it prints at the end.  If
you see an error saying that the 'python' command was not found, than
you may not have any version of Python installed.  See:

  http://www.python.org/download/releases/2.7.3/#download

for installation instructions.

This test is separate to avoid Python syntax errors parsing the more
elaborate `swc-installation-test-2.py`.
"""

import sys as _sys


__version__ = '0.1'


def check():
    if _sys.version_info < (2, 6):
        print('check for Python version (python):')
        print('outdated version of Python: ' + _sys.version)
        return False
    return True


if __name__ == '__main__':
    if check():
        print('Passed')
    else:
        print('Failed')
        print('Install a current version of Python!')
        print('http://www.python.org/download/releases/2.7.3/#download')
        _sys.exit(1)
