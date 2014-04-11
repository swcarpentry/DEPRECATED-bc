# EPUB Version of Software Carpentry Lessons

This document explains how to build the `epub` version of Software Carpentry
lessons and others ebook formats (e.g. `mobi` and `azw`).

## Dependencies

- [lxml](http://lxml.de/)
- [Jekyll](http://jekyllrb.com/)
- [Pandoc](http://johnmacfarlane.net/pandoc/)

  We recommend use pandoc >= 1.13 due a [issue with
  tables](https://github.com/jgm/pandoc/issues/1341). The build process won't
  fail if you use an older version but will hard to read tables.
- [Calibre](http://calibre-ebook.com/)

  Only need to build the `mobi` and `azw` formats.

## Build Steps

To build `epub`:

~~~
$ make epub
~~~

After build `epub` you can build `mobi` and/or `azw`:

~~~
$ make -f epub.mk mobi
$ make -f epub.mk azw
~~~
