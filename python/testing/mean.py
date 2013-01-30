def mean(numlist):
    try :
        total = sum(numlist)
        length = len(numlist)
    except TypeError :
        raise TypeError("The list was not numbers.")
    except ZeroDivisionError :
        raise ZeroDivisionError("Does your list have elements in it?")
    except :
        print "Something unknown happened with the list."
    return float(total)/float(length)
