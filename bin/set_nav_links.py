#!/usr/bin/python
"""
This set the navigation links for lessons.
"""

import os
import os.path
import re
import warnings

def write_nav(afile_, alines, aprev, anext):
    """Write updated navigation into file"""
    end_of_yaml = False
    for line in alines:
        # Locate YAML
        if '---\n' == line:
            if end_of_yaml:
                if aprev:
                    afile_.write("prev: {}\n".format(aprev))
                if anext:
                    afile_.write("next: {}\n".format(anext))
                end_of_yaml = False
            else:
                end_of_yaml = True

        # Ignoring previous navigation
        if "prev:" in line or "next:" in line:
            pass
        else:
            afile_.write(line)

def get_nav(afile_name):
    """Discovery the previous and next lesson."""
    prev_lesson = None
    next_lesson = None

    file_name_head, file_name_tail = os.path.split(afile_name)
    if file_name_head != '':
        list_files = os.listdir(file_name_head)

        pattern = re.compile('^(?P<number>\d\d)-.*\.md')
        this_lesson = pattern.match(file_name_tail)
        if this_lesson:
            this_lesson_number = int(this_lesson.group("number"))

            for file_ in list_files:
                other_lesson = pattern.match(file_)
                if other_lesson:
                    other_lesson_number = int(other_lesson.group("number"))
                    if other_lesson_number == this_lesson_number - 1:
                        prev_lesson = file_
                    if other_lesson_number == this_lesson_number + 1:
                        next_lesson = file_

    return prev_lesson, next_lesson

def set_nav(afile_name):
    """Set the navigation links into file."""
    if os.path.isfile(afile_name):
        prev_lesson, next_lesson = get_nav(afile_name)

        with open(afile_name, 'r') as file_:
            lines = file_.readlines()

        with open(afile_name, 'w') as file_:
            write_nav(file_, lines, prev_lesson, next_lesson)
    else:
        warnings.warn("{} isn't a regular file.".format(afile_name), 
                UserWarning)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
            description="Script to set navigation links")

    parser.add_argument("file_name", nargs="+",
            help="The files to set the navigation links")

    parsed_args = parser.parse_args()

    for file_name in parsed_args.file_name:
        try:
            set_nav(file_name)
        except FileNotFoundError as err:
            print(err)
