#!/usr/bin/env python

'''Test script to check whether index.html is valid.

Prints out warnings and errors when the header in index.html is malformed,
or when the information for 

Checks for:
    1. There should be the right number of categories
    2. Categories are allowed to appear only once
    3. Contact email should be valid (letters + @ + letters + . + letters)
    4. Latitute/longitude should be 2 floating point numbers separated by comma
    5. startdate/enddate should be valid dates
    6. country should be a string with no spaces
    7. instructor list should be a valid Python/Ruby list
    8. Template header should not exist
    9. humandate should have three-letter month and four-letter year
    10. layout should be 'bootcamp'
    11. root must be '.'
    12. humantime should have 'am' or 'pm' or both
    13. address, venue should be non-empty
    14. registration should be 'open' or 'restricted'
'''


import sys
import re

import yaml
from collections import Counter

__version__ = '0.4'

# See docstring for definitions and restrictions for each category
CATEGORIES = set(['layout', 'root', 'venue', 'address', 'country',
        'humandate', 'humantime', 'startdate', 'enddate', 'latlng',
        'registration', 'instructor', 'contact'])
REGISTRATIONS = set(['open', 'restricted', 'closed'])

EMAIL_PATTERN = '[^@]+@[^@]+\.[^@]+'
HUMANTIME_PATTERN = '((0?[0-9]|1[0-1]):[0-5][0-9](am|pm)(-|to)(0?[0-9]|1[0-1]):[0-5][0-9](am|pm))|((0?[0-9]|1[0-9]|2[0-3]):[0-5][0-9](-|to)(0?[0-9]|1[0-9]|2[0-3]):[0-5][0-9])'

def check_layout(layout):
    '''Checks whether layout equals "bootcamp".'''
    return layout == 'bootcamp'

def check_root(root):
    '''Checks root - can only be "."'''
    return root == '.'

def check_country(country):
    '''A valid country has no spaces, is one string, isn't empty'''
    return (country is not None) and (' ' not in country)


def check_humandate(date):
    '''A valid human date starts with a three-letter month and ends with four-letter year,
    Example: "Feb 18-20, 2525"
    other example: "Feb 18 and 20, 2014"
    '''
    if "," not in date:
        return False

    month_dates, year = date.split(",")

    # The first three characters of month_dates are not empty
    month = month_dates[:3]
    if any(char == " " for char in month):
        return False

    # But the fourth character is empty ("February" is illegal)
    if month_dates[3] != " ":
        return False

    # year contains *only* numbers
    try:
        int(year)
    except:
        return False

    return True

def check_humantime(time):
    '''A valid humantime contains at least one number'''
    return bool(re.match(HUMANTIME_PATTERN, time.replace(" ","")))

def check_date(this_date):
    '''A valid date is YEAR-MONTH-DAY, example: 2014-06-30'''
    from datetime import date
    # yaml automatically loads valid dates as datetime.date
    return isinstance(this_date, date)

def check_latitude_longitude(latlng):
    '''A valid latitude/longitude listing is two floats, separated by comma'''
    try:
        # just one of them has to break
        lat, lng = latlng.split(',')
        float(lat)
        float(lng)
    except ValueError:
        return False
    return True

def check_registration(registration):
    '''Legal registrations are defined in REGISTRATIONS'''
    return registration in REGISTRATIONS

def check_instructor(instructor):
    '''Checks whether instructor list is of format ['First instructor', 'Second instructor', ...']'''
    # yaml automatically loads list-like strings as lists
    return isinstance(instructor, list) and len(instructor) > 0

def check_email(email):
    '''A valid email has letters, then an @, followed by letters, followed by a dot, followed by letters.'''
    return bool(re.match(EMAIL_PATTERN, email))

def check_validity(data, function, error):
    '''Wrapper-function around the various check-functions.'''
    valid = function(data)
    if valid:
        return True
    else:
        sys.stderr.write(error)
        sys.stderr.write('\tOffending entry is: "{0}"\n'.format((data)))
        return False

def check_categories(first_categories, second_categories, message):
    result = first_categories - second_categories
    if result:
        sys.stderr.write(message)
        sys.stderr.write('\tOffending entries: {0}\n'.format(result))
        return False
    return True

def check_double_categories(seen_categories, message):
    category_counts = Counter(seen_categories)
    double_categories = [category for category in category_counts if category_counts[category] > 1]
    if double_categories:
        sys.stderr.write(message)
        sys.stderr.write('\t"{0}" appears more than once.\n'.format(double_categories))
        return False
    return True

def get_header(index_fh):
    '''Parses index.html file, returns just the header'''
    # We stop the header once we see the second '---'
    header_counter = 0
    header = []
    this_categories = []
    for line in index_fh:
        line = line.rstrip('\n\r')
        if line == '---':
            header_counter += 1
            continue
        if header_counter != 2:
            header.append(line)
            this_categories.append(line.split(":")[0].strip())

        if "This page is a template for bootcamp home pages." in line:
            sys.stderr.write('WARN:\tYou seem to still have the template header in your index.html. Please remove that.\n')
            sys.stderr.write('\tLook for: "<!-- Remove the block below. -->" in the index.html.\n')
            break # we can stop here - for now, just check header and template header

    return yaml.load("\n".join(header)), this_categories

def check_file(index_fh):
    '''Gets header from index.html, calls all other functions and checks file for validity.
    Returns True when 'index.html' has no problems and False when there are problems.
    '''
    header_data, seen_categories = get_header(index_fh)

    if not header_data:
        sys.stderr.write('ERROR:\tCan\'t find header in given file "{0}". Please check path, is this the bc index.html?\n'.format((filename)))
        sys.exit(1)

    file_is_valid = True

    handlers_and_messages = { 'layout' : (check_layout, 'ERROR:\tLayout isn\'t "bootcamp".\n'),
        'root' : (check_root, 'ERROR:\troot can only be ".".\n'), 
        'country' : (check_country, 'ERROR:\tCountry seems to be invalid. Please check whether there are spaces inside the country-name.\n'),
        'humandate' : (check_humandate, 'ERROR:\tCategory "humandate" seems to be invalid. Please use a three-letter month like "Jan" and four-letter year like "2025".\n'),
        'humantime' : (check_humantime, 'ERROR:\t"humantime" doesn\'t include numbers.\n'),
        'startdate' : (check_date, 'ERROR:\t"startdate" seems to be invalid. Must be of format year-month-day, i.e., 2014-01-31.\n'),
        'enddate' : (check_date, 'ERROR:\t"enddate" seems to be invalid. Must be of format year-month-day, i.e., 2014-01-31.\n'),
        'latlng' : (check_latitude_longitude, 'ERROR:\tLatitude/Longitude seems to be invalid. Check whether it\'s two floating point numbers, separated by a comma.\n'),
        'registration' : (check_registration, 'ERROR:\tregistration can only be {0}.\n'.format(REGISTRATIONS)), 
        'instructor' : (check_instructor, 'ERROR:\tInstructor string isn\'t a valid list of format ["First instructor", "Second instructor",..].\n'),
        'contact' : (check_email,'ERROR:\tEmail seems to be invalid.\n') }

    # look through all header entries
    for category in handlers_and_messages:
        handler_function, error_message = handlers_and_messages[category]
        file_is_valid &= check_validity(header_data[category], handler_function, error_message)

    # Do we have double categories?
    file_is_valid &= check_double_categories(seen_categories, 'ERROR:\tThere are categories appearing twice or more.\n')

    # Check whether we have missing or too many categories
    seen_categories = set(seen_categories)
    file_is_valid &= check_categories(seen_categories, CATEGORIES, 'ERROR:\tThere are superfluous categories.\n')
    file_is_valid &= check_categories(CATEGORIES, seen_categories, 'ERROR:\tThere are missing categories.\n')

    return file_is_valid

if __name__ == '__main__':
    args = sys.argv
    if len(args) > 2:
        sys.stderr.write('Usage: python swc-index-test.py\nAlternative: python swc-index-test.py path/to/index.html\n')
        sys.exit(0)
    elif len(args) == 1:
        filename = '../index.html'
    else:
        filename = args[1]

    sys.stderr.write('Testing file "{0}".\n'.format(filename))

    with open(filename) as index_fh:
        file_is_valid = check_file(index_fh)

    if file_is_valid:
        sys.stderr.write('Everything seems to be in order.\n')
        sys.exit(0)
    else:
        sys.stderr.write('There were problems, please see above.\n')
        sys.exit(1)
