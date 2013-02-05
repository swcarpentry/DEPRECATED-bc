def mean(numlist):
    try :
        total = sum(numlist)
        length = len(numlist)
    except TypeError :
        raise TypeError("The list was not numbers.")
    except :
        print "Something unknown happened with the list."
    return total/length
