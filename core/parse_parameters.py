def parseParameters(parFile):

    par = {}
    with open(parFile) as f:
        for line in f:
            (key, val) = line.split()
            par[key] = float(val)
        
    return par