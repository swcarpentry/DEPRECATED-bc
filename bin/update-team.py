"""
Update team.md
"""

import sys
import os.path
import subprocess
if sys.version_info.major == 3:
    import urllib.request as url
elif sys.version_info.major == 2:
    import urllib2 as url
import warnings
import yaml


TEMPLATE = """---
layout: lesson
root: .
title: Our Team
---
Many thanks to all of the people who have contributed to Software Carpentry over the years.

"""

LINK2BOOTCAMPS = "workshops_saved.yml"
URLOFLINK2BOOTCAMPS = "https://raw.githubusercontent.com/swcarpentry/site/master/config/workshops_saved.yml"

def contributor_name_is_valid(contributor_name, origin):
    """Check if contributor name is valid"""
    names = contributor_name.split()
    if len(names) < 2:
        warnings.warn("{} from {} is invalid".format(contributor_name, origin))

def get_workshops():
    """Get metadata from workshops"""
    if os.path.isfile(LINK2BOOTCAMPS):
        with open(LINK2BOOTCAMPS, "r") as file_:
            lines = ''.join(file_.readlines())
    else:
        lines = url.urlopen(URLOFLINK2BOOTCAMPS).read()
        with open(LINK2BOOTCAMPS, "w") as file_:
            file_.write(lines.decode("utf-8"))
    return yaml.load(lines)

def get_contributors_from_workshops():
    """Create a set of contributors from workshops"""
    contributors = set()
    for workshop in get_workshops():
        if "instructor" in workshop:
            for instructor in workshop["instructor"]:
                contributors.add(instructor)
        if "helper" in workshop:
            for helper in workshop["helper"]:
                contributors.add(helper)

    # Check if contributor name is valid
    for contributor in contributors:
        contributor_name_is_valid(contributor, "workshop list")

    return contributors

def get_contributors_from_git():
    """Create a set of contributors from bc repo"""
    contributors = set()
    with open("AUTHORS", "r") as file_:
        file_.readline()
        for line in file_:
            contributors.add(line.split('<')[0].strip())

    # Check if contributor name is valid
    for contributor in contributors:
        contributor_name_is_valid(contributor, "git")

    return contributors

def write_team(contributors):
    """Write team.md"""
    i = 0
    with open("team.md", "w") as file_:
        file_.write(TEMPLATE)
        file_.write("<table width=\"100%\" class=\"team\">")
        sorted_contributors = sorted(contributors, key=lambda item: item.split(" ")[-1])
        # Print contributor name
        for contributor in sorted_contributors:
            uni_contributor = contributor
            if i % 3 == 0:
                cell = "<tr class=\"team\">\n<td class=\"team\">{0}</td>".format(uni_contributor)
            elif i % 3 == 2:
                cell = "<td class=\"team\">{0}</td>\n</tr>".format(uni_contributor)
            else:
                cell = "<td class=\"team\">{0}</td>".format(uni_contributor)
            file_.write(cell)
            i = i + 1
        # Print end of table
        if i % 3 == 0:
            file_.write("</table>")
        else:
            file_.write("</tr></table>")

if __name__ == "__main__":
    contributors = get_contributors_from_workshops()
    contributors = contributors.union(get_contributors_from_git())
    write_team(contributors)
