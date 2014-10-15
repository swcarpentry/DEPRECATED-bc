"""
Test suite for ``validate_index.py``

Workshop metadata is stored in yaml header and PyYaml
strip all the strings.
"""

from io import StringIO
from datetime import date
import swc_index_validator as validator

def make_file(text):
    try: # this happens in Python3
        f = StringIO(text)
    except TypeError: # this happens in Python2
        f = StringIO(unicode(text))
    return f


def test_check_layout():
    assert validator.check_layout("workshop")

def test_check_layout_fail():
    assert not validator.check_layout("lesson")

def test_check_root():
    assert validator.check_root(".")

def test_check_root_fail():
    assert not validator.check_root("setup")

def test_check_contry():
    assert validator.check_country("Country")

def test_check_contry_none():
    assert not validator.check_country(None)

def test_check_contry_two_words():
    assert not validator.check_country("Some Country")

def test_check_humandate():
    assert validator.check_humandate("Feb 18-20, 2525")

def test_check_humandate_fail():
    assert not validator.check_humandate("February 18-20, 2525")

def test_check_humandate_chars():
    assert not validator.check_humandate("XXX SomeDay, Year")

def test_check_humantime():
    assert not validator.check_humantime("09:00am")

def test_check_euro_humantime():
    assert validator.check_humantime("09:00-17:00")

def test_check_humantime_fail():
    assert not validator.check_humantime("09:00")

def test_check_humantime_only_am():
    assert not validator.check_humantime("am")

def test_check_humantime_without_spaces():
    assert validator.check_humantime("9:00am-5:00pm")

def test_check_humantime_with_spaces():
    assert validator.check_humantime("9:00am - 5:00pm")

def test_check_humantime_with_extra_spaces():
    assert validator.check_humantime("9:00 am - 5:00 pm")

def test_check_humantime_with_to():
    assert validator.check_humantime("9:00am to 5:00pm")

def test_check_humantime_with_to_and_spaces():
    assert validator.check_humantime("9:00 am to 5:00 pm")

def test_check_humantime_without_am_pm():
    assert validator.check_humantime("9:00-17:00")

def test_check_humantime_without_am_pm_with_to():
    assert validator.check_humantime("9:00 to 17:00")

def test_check_date():
    assert validator.check_date(date(2525, 2, 20))

def test_check_date_fail():
    assert not validator.check_date("Feb 18-20, 2525")

def test_check_latitude_longitude():
    assert validator.check_latitude_longitude("0.0,0.0")

def test_check_latitude_longitude_chars():
    assert not validator.check_latitude_longitude("foo,bar")

def test_check_registration_open():
    assert validator.check_registration("open")

def test_check_registration_restricted():
    assert validator.check_registration("restricted")

def test_check_registration_closed():
    assert validator.check_registration("closed")

def test_check_registration_fail():
    assert not validator.check_registration("close")

def test_check_instructors():
    assert validator.check_instructors(["John Doe", "Jane Doe"])

def test_check_instructor_only_one():
    assert validator.check_instructors(["John Doe"])

def test_check_instructor_empty():
    assert not validator.check_instructors([])

def test_check_instructor_string():
    assert not validator.check_instructors("John Doe")

def test_check_helpers():
    assert validator.check_helpers(["John Doe", "Jane Doe"])

def test_check_instructor_only_one():
    assert validator.check_helpers(["John Doe"])

def test_check_instructor_empty():
    assert validator.check_helpers([])

def test_check_instructor_string():
    assert not validator.check_helpers("John Doe")

def test_check_email():
    assert validator.check_email("user@box.com")

def test_check_email_obfuscate():
    assert not validator.check_email("user AT box DOT com")

def test_check_email_not_default():
    assert not validator.check_email('admin@software-carpentry.org')

def test_check_eventbrite_9_digits():
    assert validator.check_eventbrite('1' * 9)

def test_check_eventbrite_10_digits():
    assert validator.check_eventbrite('1' * 10)

def test_check_not_eventbrite_8_digits():
    assert not validator.check_eventbrite('1' * 8)

def test_check_not_eventbrite_empty():
    assert not validator.check_eventbrite('')

def test_check_not_eventbrite_non_digits():
    assert not validator.check_eventbrite('1' * 8 + 'a')

def test_check_with_enddate():
    header_sample = """---
layout: workshop
root: .
venue: Euphoric State University
address: 123 College Street, Euphoria
country: United-States
humandate: Feb 17-18, 2020
humantime: 9:00 am - 4:30 pm
startdate: 2020-06-17
enddate: 2020-06-18
latlng: 41.7901128,-87.6007318
registration: restricted
instructor: ["Grace Hopper", "Alan Turing"]
helper: [ ]
contact: alan@turing.com
---"""

    assert validator.check_file(make_file(header_sample))

def test_check_without_enddate():
    header_sample = """---
layout: workshop
root: .
venue: Euphoric State University
address: 123 College Street, Euphoria
country: United-States
humandate: Feb 17-18, 2020
humantime: 9:00 am - 4:30 pm
startdate: 2020-06-17
latlng: 41.7901128,-87.6007318
registration: restricted
instructor: ["Grace Hopper", "Alan Turing"]
contact: alan@turing.com
helper: [ "John von Neumann" ]
---"""

    assert validator.check_file(make_file(header_sample))
