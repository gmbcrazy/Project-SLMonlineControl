import numpy as np

filepath = "I:\\Analyzed_x\\20230522_mouse9-1_p8-results\suite2p\plane0\ops.npy"
ops = np.load(filepath, allow_pickle=True)
#print(ops.item())


list = ops.item()
dictionary = {}

for key,value in list.items():
    dictionary[key] = str(value)

#for key in dictionary:
    #print(key)
print(dictionary['timing'])    