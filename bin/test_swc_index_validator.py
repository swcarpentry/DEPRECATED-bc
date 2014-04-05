"""
Test suite for ``validate_index.py``

Bootcamp metadata is stored in yaml header and PyYaml
strip all the strings.
"""

from io import StringIO
from datetime import date
import swc_index_validator

def test_check_layout():
    assert swc_index_validator.check_layout("bootcamp")

def test_check_layout_fail():
    assert not swc_index_validator.check_layout("lesson")

def test_check_root():
    assert swc_index_validator.check_root(".")

def test_check_root_fail():
    assert not swc_index_validator.check_root("setup")

def test_check_contry():
    assert swc_index_validator.check_country("Country")

def test_check_contry_none():
    assert not swc_index_validator.check_country(None)

def test_check_contry_two_words():
    assert not swc_index_validator.check_country("Some Country")

def test_check_humandate():
    assert swc_index_validator.check_humandate("Feb 18-20, 2525")

def test_check_humandate_fail():
    assert not swc_index_validator.check_humandate("February 18-20, 2525")

def test_check_humandate_chars():
    assert not swc_index_validator.check_humandate("XXX SomeDay, Year")

def test_check_humantime():
    assert not swc_index_validator.check_humantime("09:00am")

def test_check_euro_humantime():
    assert swc_index_validator.check_humantime("09:00-17:00")

def test_check_humantime_fail():
    assert not swc_index_validator.check_humantime("09:00")

def test_check_humantime_only_am():
    assert not swc_index_validator.check_humantime("am")

def test_check_humantime_without_spaces():
    assert swc_index_validator.check_humantime("9:00am-5:00pm")

def test_check_humantime_with_spaces():
    assert swc_index_validator.check_humantime("9:00am - 5:00pm")

def test_check_humantime_with_extra_spaces():
    assert swc_index_validator.check_humantime("9:00 am - 5:00 pm")

def test_check_humantime_with_to():
    assert swc_index_validator.check_humantime("9:00am to 5:00pm")

def test_check_humantime_with_to_and_spaces():
    assert swc_index_validator.check_humantime("9:00 am to 5:00 pm")

def test_check_humantime_without_am_pm():
    assert swc_index_validator.check_humantime("9:00-17:00")

def test_check_humantime_without_am_pm_with_to():
    assert swc_index_validator.check_humantime("9:00 to 17:00")

def test_check_date():
    assert swc_index_validator.check_date(date(2525, 2, 20))

def test_check_date_fail():
    assert not swc_index_validator.check_date("Feb 18-20, 2525")

def test_check_latitude_longitude():
    assert swc_index_validator.check_latitude_longitude("0.0,0.0")

def test_check_latitude_longitude_chars():
    assert not swc_index_validator.check_latitude_longitude("foo,bar")

def test_check_registration_open():
    assert swc_index_validator.check_registration("open")

def test_check_registration_restricted():
    assert swc_index_validator.check_registration("restricted")

def test_check_registration_closed():
    assert swc_index_validator.check_registration("closed")

def test_check_registration_fail():
    assert not swc_index_validator.check_registration("close")

def test_check_instructor():
    assert swc_index_validator.check_instructor(["John Doe", "Jane Doe"])

def test_check_instructor_only_one():
    assert swc_index_validator.check_instructor(["John Doe"])

def test_check_instructor_empty():
    assert not swc_index_validator.check_instructor([])

def test_check_instructor_string():
    assert not swc_index_validator.check_instructor("John Doe")

def test_check_email():
    assert swc_index_validator.check_email("user@box.com")

def test_check_email_obfuscate():
    assert not swc_index_validator.check_email("user AT box DOT com")

def test_check_file():
    header_sample = """---
layout: bootcamp
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
contact: admin@software-carpentry.org
---"""

    try: # this happens in Python3
        file_ = StringIO(header_sample)
    except TypeError: # this happens in Python2
        file_ = StringIO(unicode(header_sample))

    assert swc_index_validator.check_file(file_)
