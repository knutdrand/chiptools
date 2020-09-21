import numpy as np

def alignscore(lines, read_length):
    a = np.zeros(read_length+1, dtype="int")
    scores = (int(line.split()[-2].split(":")[-1])
              for line in lines if not line.startswith("@"))
    for score in scores:
        a[score]+=1
    return a
        
              
