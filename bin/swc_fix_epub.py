#!/usr/bin/python

"""
Fix XHTML from EPUB generate with pandoc
"""

import sys
import os.path

# Need lib to parse XML and we can't use xml since it didn't  fine-grained
# control over serializing namespaced elements
import lxml.etree as etree

def find_glossary(file_name):
    tree = etree.parse(file_name)
    links = tree.findall('.//{http://www.w3.org/1999/xhtml}a')
    for link in links:
        if link.text == "Glossary":
            return link.attrib["href"]

def generate_name(string):
    """Generate the name attribute for the entries at the glossary."""
    return string.lower().replace(' ', '-').replace("'", '')

def fix_glossary(file_name):
    tree = etree.parse(file_name)
    links = tree.findall('.//{http://www.w3.org/1999/xhtml}strong')
    for link in links:
        fixed_link = etree.Element('{http://www.w3.org/1999/xhtml}a',
            name=generate_name(link.text))
        fixed_link.text = link.text
        link.append(fixed_link)
        link.text = ""

    tree.write(file_name, xml_declaration=True, encoding="UTF-8")

def fix_link2glossary(file_name, glossary_file_name):
    tree = etree.parse(file_name)
    links = tree.findall('.//{http://www.w3.org/1999/xhtml}a')
    for link in links:
        if 'href' in link.attrib and link.attrib['href'].startswith('#g:'):
            link.attrib['href'] = link.attrib['href'].replace('#g:',
                    '{0}#'.format(glossary_file_name))

    tree.write(file_name, xml_declaration=True, encoding="UTF-8")

def main(path_to_files):
    gloss = find_glossary(os.path.join(path_to_files, "ch001.xhtml"))
    fix_glossary(os.path.join(path_to_files, gloss))
    for file_name in os.listdir(path_to_files):
        if file_name.endswith(".xhtml"):
            fix_link2glossary(os.path.join(path_to_files, file_name), gloss)

if __name__ == "__main__":
    if len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        print("Usage: ./fix_epub.py path/to/xhtmls")

