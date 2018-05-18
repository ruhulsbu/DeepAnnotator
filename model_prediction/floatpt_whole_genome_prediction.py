import os, sys, gzip
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
max_data_size = 11000# 20000000        #Original 1.5 Mil
char_to_int = dict((c, i) for i, c in enumerate(alphabet))
int_to_char = dict((i, c) for i, c in enumerate(alphabet))
CODING = True
MAX_GENOME = 100


def read_input_file(file_path, label=-1):
    file_count = 0
    file_name = []
    line_number = []

    x_data = []
    count_line = 0
    file_read = gzip.open(file_path, "rt")

    for line in file_read:
        if line[0] == '>':
            file_count += 1
            if file_count > MAX_GENOME:
                break

            file_name.append(line.strip())
            line_number.append(count_line)
            print(file_name[-1], line_number[-1])
            continue

        data = [int(i) for i in line.strip()]
        x_data.append(data)
        count_line += 1
        #print(x_data[-1], y_data[-1])
        #if len(x_data) == max_data_size:
        #    break

    file_read.close()
    print("Sequences Read: ", len(x_data))
    return file_name, line_number, np.array(x_data)


file_name, line_number, x_data = read_input_file(sys.argv[1], 1)
mod = batch_size - (len(x_data) % batch_size)
big_array = np.zeros((mod,time_steps))
x_data = np.append(x_data, big_array, axis=0)
print("Reformatted Input: ", x_data.shape)
line_number.append(-1)


def return_model(prefix, check_point, layer_name):

    Batch = namedtuple("Batch", ["data"])
    sym, arg_params, aux_params = mx.model.load_checkpoint(prefix, check_point)

    all_layers = sym.get_internals()
    #print(all_layers.list_outputs()[:])

    softmax_sym = all_layers[layer_name]
    print(softmax_sym)
    
    devices = [mx.gpu(i) for i in range(8)]
    model_loaded = mx.mod.Module(symbol=sym, context=devices, data_names=["data"], label_names=["label"])
    model_loaded.bind(for_training=False, label_shapes=None, data_shapes=[("data", (batch_size, time_steps))])
    model_loaded.set_params(arg_params, aux_params)
    
    return model_loaded


def compute_accuracy(x_data, trained_model):
    Batch = namedtuple("Batch", ["data"])
    prediction = []

    for index in range(0, len(x_data)-batch_size+1, batch_size):
        trained_model.forward(Batch([mx.nd.array(x_data[index:index+batch_size])]))
        results = trained_model.get_outputs()[0].asnumpy()
        prediction.extend(results)
        
    prediction = np.array(prediction).ravel()
    result = np.array([int(prediction[i] >= 0.50) for i in range(len(prediction))])
    print("Prediction: ", result.shape)
    return prediction


model_name = ["gene_prediction_lstm_model", "general_lstm_model_start", "general_lstm_model_stop"]
check_point = [1, 1, 200]
predict = []
#file_write = gzip.open("predicted_whole_sequences.gz", "wt")


for i in range(len(model_name)):
    trained_model = return_model(model_name[i], check_point[i], "logistic_output")
    y_data = compute_accuracy(x_data, trained_model)
    predict.append(y_data)


count = 0
organism = 0
print_string = [[], [], []]

for i in range(0, len(y_data)):
    if i == line_number[organism]:
        if count > 0:
            for k in range(len(print_string)):
                file_write.write(",".join(x for x in print_string[k])+"\n")
                print_string[k] = []
            file_write.close()
        
        file_write = gzip.open(sys.argv[2] + "floatpt_whole_sequences_"+file_name[organism][1:], "wt")
        file_write.write(file_name[organism]+"\n")
        organism += 1
        count = 0

    count = count + 1
    for k in range(len(print_string)):
        print_string[k].append(str(predict[k][i])) 
    
if count > 0:
    for k in range(len(print_string)):
        file_write.write(",".join(x for x in print_string[k])+"\n")
file_write.close()

print("Writing Done! ", i)
