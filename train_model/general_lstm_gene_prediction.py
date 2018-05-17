import os, sys, gzip
import random, h5py
import numpy as np
import mxnet as mx

#Initialize the Program
alphabet = "NACGT"
vocab_size = 5
batch_size = 10000		#Original 10000
embedding_size = 4		#Original 8
time_steps = 101
category = 2
max_data_size = 10000000	#Original 1.5 Mil
char_to_int = dict((c, i) for i, c in enumerate(alphabet))
int_to_char = dict((i, c) for i, c in enumerate(alphabet))
CODING = True

#Read the Input File

def read_input_file(file_path, label=-1):
    x_data = []
    y_data = []

    file_read = gzip.open(file_path, "rt")
    for line in file_read:
        line = str(line.strip())
        data = [int(line[i]) for i in range(len(line))]        

        x_data.append(data)
        y_data.append(label)
        #print(x_data[-1], y_data[-1])
        if len(x_data) == max_data_size:
            break
    file_read.close()
    print("Sequences Read: ", len(x_data))
    return np.array(x_data), np.array(y_data)

x_data_pos, y_data_pos = read_input_file("./coding_region_sequences.gz", 1)
x_data_neg, y_data_neg = read_input_file("./intrag_region_sequences.gz", 0)

print(x_data_pos.shape, x_data_neg.shape)

np.random.shuffle(x_data_pos)
np.random.shuffle(x_data_neg)
#exit()
#Process Data Index

train_index = int((len(x_data_pos) / batch_size) * 0.60 * batch_size)
eval_index = train_index + int((len(x_data_pos) / batch_size) * 0.20 * batch_size)
test_index = eval_index + int((len(x_data_pos) / batch_size) * 0.20 * batch_size)
print("train, eval, test = ", (train_index, eval_index, test_index))

#Process Negative Data

x_train = x_data_neg[0:train_index]
y_train = y_data_neg[0:train_index]

x_eval = x_data_neg[train_index:eval_index]
y_eval = y_data_neg[train_index:eval_index]

x_test = x_data_neg[eval_index:test_index]
y_test = y_data_neg[eval_index:test_index]

#Process Positive Data

x_train = np.append(x_train, x_data_pos[0:train_index], axis=0)
y_train = np.append(y_train, y_data_pos[0:train_index], axis=0)

x_eval = np.append(x_eval, x_data_pos[train_index:eval_index], axis=0)
y_eval = np.append(y_eval, y_data_pos[train_index:eval_index], axis=0)

x_test = np.append(x_test, x_data_pos[eval_index:test_index], axis=0)
y_test = np.append(y_test, y_data_pos[eval_index:test_index], axis=0)

print("Sanity Check: ", np.sum(y_train), np.sum(y_eval), np.sum(y_test))

#Create Data Iterator

train_iter = mx.io.NDArrayIter(data=x_train, label=y_train, \
                               data_name='data', label_name='label', \
                               batch_size=batch_size, shuffle=True)
print(train_iter.data)

eval_iter = mx.io.NDArrayIter(data=x_eval, label=y_eval, \
                               data_name='data', label_name='label', \
                               batch_size=batch_size, shuffle=True)
print(eval_iter.data)

test_iter = mx.io.NDArrayIter(data=x_test, label=y_test, \
                               data_name='data', label_name='label', \
                               batch_size=batch_size, shuffle=False)
print(test_iter.data)

###############################################################################

unroll_dim = time_steps #Length of sequence
sequence_dim = 2*unroll_dim #Sum of two sequence length

num_hidden = 128 #Output unit of lstm layer
num_embed = embedding_size #Dimension of embeddings
num_label = category #dense network output

max_features = vocab_size #Number of events
#multiply = 100 #Control the size of data

session = mx.sym.Variable(name='data')
session_src = mx.sym.slice_axis(data=session, axis=1, begin=0, end=unroll_dim, name='source')
#session_tar = mx.sym.slice_axis(data=session, axis=1, begin=unroll_dim, end=sequence_dim, name='target')
embed_weight = mx.sym.Variable(name='weight')
label = mx.sym.Variable(name='label')

#Shared Embedding Layer
#Original Dropout 0.40
embed_seq_src = mx.sym.Embedding(data=session_src, \
                                input_dim=max_features, \
                                output_dim=num_embed, \
                                weight=embed_weight, \
                                name='vocab_embed_src')
"""
embed_seq_tar = mx.sym.Embedding(data=session_tar, \
                                input_dim=max_features, \
                                output_dim=num_embed, \
                                weight=embed_weight, \
                                name='vocab_embed_tar')
"""
#Shared LSTM Layer

fused_lstm_src = mx.rnn.FusedRNNCell(num_hidden=num_hidden, \
                                mode='lstm', \
                                prefix='lstm_src_', \
                                forget_bias=True, \
				bidirectional=True, \
                                dropout=0.20)	#Original 0.50
fused_lstm_src.reset()

lstm_outputs_src, _ = fused_lstm_src.unroll(length=unroll_dim, \
                                inputs=embed_seq_src, \
                                merge_outputs=True, \
                                layout="NTC")

fused_lstm_tar = mx.rnn.FusedRNNCell(num_hidden=num_hidden, \
                                mode='lstm', \
                                prefix='lstm_tar_', \
                                #params = fused_lstm_src.params, \
                                forget_bias=True, \
				bidirectional=True, \
                                dropout=0.20)	#Original 0.50
fused_lstm_tar.reset()

lstm_outputs_tar, _ = fused_lstm_tar.unroll(length=unroll_dim, \
                                inputs=lstm_outputs_src, \
                                merge_outputs=False, \
                                layout="NTC")

#Logistic Regression Layer

out_lstm1 = mx.sym.Reshape(data=lstm_outputs_tar[-1], shape=(int(batch_size/8), -1), name="reshape_lstm_src")
#out_lstm2 = mx.sym.Reshape(data=lstm_outputs_tar[-1], shape=(int(batch_size/8), -1), name="reshape_lstm_tar")

norm_lstm1 = mx.sym.L2Normalization(data=out_lstm1, mode='instance', name="norm_lstm_src") 
#norm_lstm2 = mx.sym.L2Normalization(data=out_lstm2, mode='instance', name="norm_lstm_tar")
#norm_product = norm_lstm1 * norm_lstm2

#scalar_sum = mx.sym.sum_axis(data=norm_product, axis=1, name="sum_dot_product")
#scalar_output = mx.sym.Reshape(data=scalar_sum, shape=(int(batch_size/8), 1), name="reshape_scalar")

#shared_lstm_model = mx.sym.LogisticRegressionOutput(data=scalar_output, label=label, name="logistic")

dense_output = mx.sym.FullyConnected(data=out_lstm1, num_hidden=num_label, name="dense_net")
shared_lstm_model = mx.sym.SoftmaxOutput(data=dense_output, label=label, name="softmax_label")

dense_output = mx.sym.FullyConnected(data=out_lstm1, num_hidden=1, name="dense_net")
logistic_model = mx.sym.LogisticRegressionOutput(data=dense_output, label=label, name="logistic")

devices = [mx.gpu(i) for i in range(8)]
shared_model = mx.mod.Module(symbol=logistic_model, \
                                 data_names=['data'], \
                                 label_names=['label'], \
                                 context=devices)

#Train the Model

import logging
logging.getLogger().setLevel(logging.logMultiprocessing)  # logging to stdout

root_logger = logging.getLogger()
str_handler = logging.FileHandler("gene_logging_logregnet.txt", "w")
root_logger.addHandler(str_handler)
root_logger.setLevel(logging.INFO)

learn_rate = 0.001
opt_params = {
    'learning_rate': learn_rate,
    'wd': 0.00001
}

test_train = mx.io.NDArrayIter(data=x_train, label=y_train, \
                               data_name='data', label_name='label', \
                               batch_size=batch_size, shuffle=False)

test_eval = mx.io.NDArrayIter(data=x_eval, label=y_eval, \
                               data_name='data', label_name='label', \
                               batch_size=batch_size, shuffle=False)

max_accuracy = 0.0
saved_model = None

def check_point_analysis():
    def _callback(iter_no, sym, arg, aux):
        global max_accuracy
        global saved_model

        print("Check Point: ", iter_no)
        output = shared_model.predict(test_eval).asnumpy()
        y_prediction = np.round(output).ravel()
        print("Evaluating Accuracy at Thr=0.50: ", output.shape, (1.0 * np.sum(np.equal(y_prediction, y_eval)) / len(y_eval)))

        output = shared_model.predict(test_iter).asnumpy()
        y_prediction = np.round(output).ravel()
        accuracy = 1.0 * np.sum(np.equal(y_prediction, y_test)) / len(y_test)
        print("Testing Accuracy at Thr=0.50: ", output.shape, (1.0 * np.sum(np.equal(y_prediction, y_test)) / len(y_test)))

        output = shared_model.predict(test_train).asnumpy()
        y_prediction = np.round(output).ravel()
        print("Training Accuracy at Thr=0.50: ", output.shape, (1.0 * np.sum(np.equal(y_prediction, y_train)) / len(y_train)))
    
        if max_accuracy < accuracy:
            max_accuracy = accuracy
            saved_model = shared_model
    return _callback

check_point = 1 #200
for i in range(check_point):
    shared_model.fit(train_iter,  					# train data
                      eval_data=eval_iter,  				# validation data
                      optimizer='adam', 
                      optimizer_params=opt_params,  
                      #initializer=mx.init.Xavier(factor_type="in", magnitude=2.34),
                      eval_metric=mx.metric.CompositeEvalMetric(metrics=["acc"]), 
                      #batch_end_callback=mx.callback.Speedometer(max_data_size, 1), 
                      epoch_end_callback=check_point_analysis(),
                      num_epoch=50)  
    """
    print("Check Point: ", i)
    output = shared_model.predict(test_eval).asnumpy()
    y_prediction = np.round(output).ravel()
    print("Evaluating Accuracy at Thr=0.50: ", output.shape, (1.0 * np.sum(np.equal(y_prediction, y_eval)) / len(y_eval)))

    output = shared_model.predict(test_iter).asnumpy()
    y_prediction = np.round(output).ravel()
    accuracy = 1.0 * np.sum(np.equal(y_prediction, y_test)) / len(y_test)
    print("Testing Accuracy at Thr=0.50: ", output.shape, (1.0 * np.sum(np.equal(y_prediction, y_test)) / len(y_test)))

    output = shared_model.predict(test_train).asnumpy()
    y_prediction = np.round(output).ravel()
    print("Training Accuracy at Thr=0.50: ", output.shape, (1.0 * np.sum(np.equal(y_prediction, y_train)) / len(y_train)))
    
    if sys.argv[1] == "stop":
        if accuracy > 0.90:
            learn_rate = 0.0002
    if sys.argv[1] == "start":
        if accuracy > 0.80:
            learn_rate = 0.0001
    """
print("Done!")

#Save the Model

prefix = "gene_prediction_lstm_model"
saved_model.save_checkpoint(prefix=prefix, epoch=check_point)
netgraph = mx.viz.plot_network(logistic_model, title="gene_prediction_lstm_architecture",save_format='pdf')
netgraph.render()

#Evaluate Model

output = saved_model.predict(test_train).asnumpy()
#y_prediction = np.round(output).ravel()
#print("\nTraining Accuracy at Thr=0.50: ", output.shape, (1.0 * np.sum(np.equal(y_prediction, y_train)) / len(y_train)))

from sklearn.metrics import roc_curve, auc, f1_score
fpr, tpr, thr = roc_curve(y_true=y_train, y_score=output)

score = []
for i in range(len(thr)):
    tp = tpr[i]
    fn = 1 - tpr[i]
    fp = fpr[i]
    tn = 1 - fpr[i]
    
    maximize = tp + tn
    
    #print(tpr, fnr, fpr, tnr, maximize)
    #print(maximize)
    score.append(maximize)

index = np.argmax(score)
print(index, fpr[index], tpr[index], thr[index], score[index])

result = [int(output[i] >= thr[index]) for i in range(len(output))]
print("ROC adjusted Accuracy: ", (1.0 * np.sum(np.equal(result, y_train)) / len(y_train)))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr)

plt.plot(fpr[index], tpr[index], 'ro')
plt.xlabel('False +VE Rate')
plt.ylabel('True +VE Rate')
plt.savefig("gene_prediction_general_lstm_roc.png")

# Test Accuracy

output = saved_model.predict(test_iter).asnumpy()
result = [int(output[i] >= thr[index]) for i in range(len(output))]
print("Final Testing Accuracy: ", (1.0 * np.sum(np.equal(result, y_test)) / len(y_test)))

from sklearn import metrics
#print "Confusion Matrix: "
#print(metrics.confusion_matrix(y_true=y_train, y_pred=y_prediction))
print("Classification Report: ")
print(metrics.classification_report(y_true=y_test, y_pred=result))

