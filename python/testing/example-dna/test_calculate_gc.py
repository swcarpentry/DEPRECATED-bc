from nose.tools import assert_equal, assert_almost_equal, assert_true, \
    assert_false, assert_raises, assert_is_instance

from calculate_gc import calculate_gc

def test_only_G_and_C():
    '''
    Sequence of only G's and C's has fraction 1.0
    '''
    fixture = 'GGCGCCGGC'
    result = calculate_gc(fixture)
    assert_equal(result, 1.0)

def test_half():
    '''
    Sequence with half G and C has fraction 0.5
    '''
    fixture = 'ATGC'
    result = calculate_gc(fixture)
    assert_equal(result, 0.5)

def test_lower_case():
    '''
    Sequence with lower case letters
    '''
    fixture = 'atgc'
    result = calculate_gc(fixture)
    assert_equal(result, 0.5)

def test_not_DNA():
    '''
    Raise TypeError if not DNA
    '''
    fixture = 'qwerty'
    assert_raises(TypeError, calculate_gc, fixture)
