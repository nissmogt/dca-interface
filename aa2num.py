def aa2num(aa):
    """convert aa into num"""
    alphabet = "ARNDCQEGHILKMFPSTWYV-"
    states = len(alphabet)
    a2n = {}
    for a, n in zip(alphabet, range(states)):
        a2n[a] = n
    if aa in a2n:
        return a2n[aa]
    else:
        return a2n['-']
