def mean(numlist):
    try:
        total = sum(numlist)
        length = len(numlist)
    except ValueError:
        print "The number list was not a list of numbers."
    except:
        print "There was a problem evaluating the number list."
    return total/length

