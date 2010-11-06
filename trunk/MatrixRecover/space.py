from math import *

def nchoosek(n, k):
    return factorial(n)/(factorial(n-k)*factorial(k))


def spaceRatio(n, l):
    matrixSpace = 2**(n*l)/factorial(n)
    constraintSpace = (n+1)**(nchoosek(l,2)+1)
    print '|D|= ', constraintSpace
    ratio = matrixSpace/constraintSpace
    print ratio
    return ratio


