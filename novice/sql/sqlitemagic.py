"""sqlitemagic provices a simple magic for interacting with SQLite
databases stored on disk.

Usage:

%%sqlite filename.db
select personal, family from person;

produces:

Alan|Turing
Grace|Hopper
"""

# This file is copyright 2013 by Greg Wilson: see
# https://github.com/gvwilson/sqlitemagic/blob/master/LICENSE
# for the license.
# Inspired by https://github.com/tkf/ipython-sqlitemagic.

import sqlite3
from IPython.core.magic import Magics, magics_class, cell_magic
from IPython.display import display, HTML

@magics_class
class SqliteMagic(Magics):
    '''Provide the 'sqlite' calling point.'''

    @cell_magic
    def sqlite(self, filename, query):
        connection = sqlite3.connect(filename)
        cursor = connection.cursor()
        try:
            cursor.execute(query)
            results = cursor.fetchall()
            display(HTML(self.tablify(results)))
        except Exception, e:
            import sys
            print >> sys.stderr, "exception", e
        cursor.close()
        connection.close()

    def tablify(self, rows):
        return '<table>\n' + '\n'.join(self.rowify(r) for r in rows) + '\n</table>'

    def rowify(self, row):
        return '<tr>' + ''.join('<td>' + str(r) + '</td>' for r in row) + '</tr>'

def load_ipython_extension(ipython):
    ipython.register_magics(SqliteMagic)
