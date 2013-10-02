"""
ears.py : a simple unit testing library for teaching in the IPython Notebook.

ears.run() looks for all the functions defined in the calling stack frame
(hopefully the top level interpreter session) whose names begin with the
characters 'test_' and calls them in an undetermined order, collecting and
reporting results.

Usage:

    import ears

    def test_pass(): pass
    def test_fail(): assert False, 'Error message'
    def test_error(): 1/0 # zero division error

    ears.run()
"""

import sys
import inspect
import traceback

def run(prefix='test_', verbose=False):
    """
    Look for test functions defined by caller, execute, and report.
    """
    # Collect functions defined in calling context.
    caller_defs = inspect.stack()[1][0].f_globals
    test_functions = dict([(n, caller_defs[n]) for n in caller_defs
                           if n.startswith(prefix) and callable(caller_defs[n])])
    setup = caller_defs.get('setup', None)
    teardown = caller_defs.get('teardown', None)

    # Execute and record.
    passes = []
    fails = []
    errors = []
    for (name, test) in test_functions.iteritems():
        if verbose:
            print name
        if setup is not None:
            setup()
        try:
            test()
            passes.append((name, None))
            sys.stdout.write('.')
        except AssertionError as e:
            fails.append((name, traceback.format_exc()))
            sys.stdout.write('f')
        except Exception as e:
            errors.append((name, traceback.format_exc()))
            sys.stdout.write('E')
        if teardown is not None:
            teardown()

    # Report.
    print
    print '{0} pass, {1} fail, {2} error'.format(len(passes),
                                                 len(fails),
                                                 len(errors))
    for (title, group) in (('fail', fails),
                           ('error', errors)):
        for (name, exc) in group:
            print '{0}\n{1}: {2}'.format('-'*40, title, name)
            print exc

if __name__ == '__main__':

    def test_pass():
        pass

    def test_fail():
        assert False, 'Error message'

    def test_error():
        1/0

    run()
