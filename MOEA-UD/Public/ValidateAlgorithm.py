"""
Validate algorithm.
"""

def validateAlgorithm(algorithm):
    """Returns True if a valid algorithm is given, returns False otherwise"""
    valid_algorithms = ['MOEA-UD']
    if algorithm in valid_algorithms:
        return True
    else:
        return False
