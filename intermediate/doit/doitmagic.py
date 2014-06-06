"""doitmagic provices a simple magic for running
tasks through doit in the ipython notebook

Usage:

%%doit doit_args

def task_example():
    return { 'actions' : ['echo "Hello world!"'] }

"""

# This file is copyright 2014 by Rob Beagrie: see
# https://github.com/gvwilson/sqlitemagic/blob/master/LICENSE
# for the license.
# Inspired by https://github.com/tkf/ipython-sqlitemagic
# and Greg Wilson's sqlitemagic elsewhere in the SWC repo

from tempfile import NamedTemporaryFile
from subprocess import Popen, PIPE
import signal
import time
import sys
import os
from IPython.core.magic import Magics, magics_class, cell_magic
from IPython.utils import py3compat

@magics_class
class DoitMagic(Magics):
    '''Provide the 'doit' calling point.'''

    def __init__(self, shell):
        """
        Parameters
        ----------
        shell : IPython shell

        """
        super(DoitMagic, self).__init__(shell)
        self._temp_file = NamedTemporaryFile()

    @cell_magic
    def doit(self, doit_args, cell):
        with NamedTemporaryFile(delete=False, suffix='.py') as tmp_file:
            tmp_name = tmp_file.name
            tmp_file.write(cell)

        cur_dir = os.getcwd()
        doit_args = doit_args.split()
        if doit_args:
            doit_command = [doit_args.pop(0)]
        else:
            doit_command = []

        cmd = ['doit']
        cmd += doit_command
        cmd += [ '-d', cur_dir, '-f', tmp_name]
        cmd += doit_args

        p = Popen(cmd, stdout=PIPE, stderr=PIPE)

        try:
            out, err = p.communicate(cell)
        except KeyboardInterrupt:
            try:
                p.send_signal(signal.SIGINT)
                time.sleep(0.1)
                if p.poll() is not None:
                    print("Process is interrupted.")
                    return
                p.terminate()
                time.sleep(0.1)
                if p.poll() is not None:
                    print("Process is terminated.")
                    return
                p.kill()
                print("Process is killed.")
            except OSError:
                pass
            except Exception as e:
                print("Error while terminating subprocess (pid=%i): %s" \
                    % (p.pid, e))
            return

        out = py3compat.bytes_to_str(out)
        err = py3compat.bytes_to_str(err)

        sys.stdout.write(out)
        sys.stdout.flush()
        sys.stderr.write(err)
        sys.stderr.flush()

        os.remove(tmp_name)

def load_ipython_extension(ipython):
    ipython.register_magics(DoitMagic)
