import os, sys
import random, h5py
import numpy as np
import mxnet as mx
from collections import namedtuple


#Initialize the Program
alphabet = "NACGT"
vocab_size = 5
batch_size = 10000              #Original 10000
embedding_size = 4              #Original 8
time_steps = 101
category = 2
max_data_size = 15200000        #Original 1.5 Mil
char_to_int = dict((c, i) for i, c in enumerate(alphabet))
int_to_char = dict((i, c) for i, c in enumerate(alphabet))
CODING = True

#Read the Input File

def read_input_file(file_path, label=-1):
    x_data = []
    y_data = []

    file_read = open(file_path, "r")
    for line in file_read:
        data = [int(i) for i in line.strip()]
        x_data.append(data)
        y_data.append(label)
        #print(x_data[-1], y_data[-1])
        if len(x_data) == max_data_size:
            break
    file_read.close()
    print("Sequences Read: ", len(x_data))
    return np.array(x_data), np.array(y_data)

x_data_pos, y_data_pos = read_input_file("./gene_range_"+sys.argv[1]+"_codon.txt", 1)
x_data_nogene, y_data_nogene = read_input_file("./intragenic_"+sys.argv[1]+"_codon.txt", 0)
x_data_coding, y_data_coding = read_input_file("./coding_"+sys.argv[1]+"_codon.txt", 0)

prefix = sys.argv[2]
if sys.argv[1] == "stop":
    check_point = 200
else:
    check_point = 1

Batch = namedtuple("Batch", ["data"])
sym, arg_params, aux_params = mx.model.load_checkpoint(prefix, check_point)

all_layers = sym.get_internals()
print(all_layers.list_outputs()[:])

#layer_sym = all_layers["logistic_output"]
#print(layer_sym)
#devices = [mx.gpu(i) for i in range(8)]
#model_loaded = mx.mod.Module(symbol=layer_sym, context=devices, data_names=['data'], label_names=['label'])
#exit()

def return_model(layer_name):
    softmax_sym = all_layers[layer_name]
    print(softmax_sym)
    
    devices = [mx.gpu(i) for i in range(8)]
    model_loaded = mx.mod.Module(symbol=sym, context=devices, data_names=["data"], label_names=["label"])
    model_loaded.bind(for_training=False, label_shapes=None, data_shapes=[("data", (batch_size, time_steps))])
    model_loaded.set_params(arg_params, aux_params)
    
    return model_loaded

trained_model = return_model("logistic_output")
#exit()

def compute_accuracy(x_data, y_data):

    prediction = []
    for index in range(0, len(x_data)-batch_size+1, batch_size):
        trained_model.forward(Batch([mx.nd.array(x_data[index:index+batch_size])]))
        results = trained_model.get_outputs()[0].asnumpy()
        prediction.extend(results)
        
    prediction = np.array(prediction).ravel()
    index = len(prediction)
    print(y_data.shape, prediction.shape)
    #return prediction

    result = [int(prediction[i] >= float(sys.argv[3])) for i in range(index)]
    accuracy = (1.0 * np.sum(np.equal(result, y_data[:index])) / index)
    print("Accuracy: ", np.count_nonzero(y_data[:index]), np.count_nonzero(result))
    return accuracy


accuracy = compute_accuracy(x_data_pos, y_data_pos)
print("ROC Adjusted True Gene Accuracy: ", accuracy)

accuracy = compute_accuracy(x_data_nogene, y_data_nogene)
print("ROC Adjusted Outside Gene Accuracy: ", accuracy)

accuracy = compute_accuracy(x_data_coding, y_data_coding)
print("ROC Adjusted Inside Gene Accuracy: ", accuracy)



