# Run this as
#
#  python swc_installation_test.py
#
# If if says nothing, everything is fine!

import os

python_modules = ['nose']

tools = ['bash',
         ('easy_install', 'Python setuptools'),
         ('hg', 'Mercurial'),
         'make',
         ('nosetests', 'Python nose'),
         'sqlite3']

# Check Python modules/packages

def check_python_modules(module_names):
    for module_name in module_names:
        try:
            __import__(module_name)
        except ImportError:
            print "Python module '%s' is missing" % module_name


# Check command line tools

def check_command_line_tools(tools):
    shell_path = os.environ['PATH'].split(':')
    for tool in tools:
        if isinstance(tool, basestring):
            command = tool
            package = None
        else:
            command, package = tool
        found = False
        for directory in shell_path:
            filename = os.path.join(directory, command)
            if os.access(filename, os.X_OK) and not os.path.isdir(filename):
                found = True
                break
        if not found:
            if package is None:
                print "Command line tool '%s' is missing" % command
            else:
                print "Command line tool '%s' " \
                    "from package '%s' is missing" % (command, package)


# Run all the checks

def main():
    check_python_modules(python_modules)
    check_command_line_tools(tools)


if __name__ == '__main__':
    main()
