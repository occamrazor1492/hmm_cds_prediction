from pomegranate import *
import numpy as np

with open("temp.txt", "r") as myfile:
    data = myfile.read().replace('\n','')

hmm = HiddenMarkovModel()
hmm.add_states()

trans, ems = hmm.forward_backward(data)
print(trans)
