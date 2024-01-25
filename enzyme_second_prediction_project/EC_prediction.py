import os
import argparse
import subprocess
import sys
import torch
import torch.nn as nn
import utility as u
import models as md
import numpy as np
import time
import pandas as pd


parser = argparse.ArgumentParser(description='ESM-ECP second (Enzyme Commission Predicting)')
parser.add_argument('fasta_data', type=str, help='fasta file path')
parser.add_argument('-o','--output', type=str, help='output file path', default=f'{os.getcwd()}/result')
parser.add_argument('-g', '--cuda_num', help='Give cuda number, such as "-g 0,1", if no gpu, cpu will be used', default=None)
parser.add_argument('-b', '--batch_size', type=int, default=2)
args = parser.parse_args()


input = args.fasta_data
output = args.output
batch_size = args.batch_size


if args.cuda_num is not None:
    os.environ['CUDA_VISIBLE_DEVICES'] = args.cuda_num
    device=torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
else:
    device = 'cpu'

if os.path.exists(output) is False:
    subprocess.call([f"mkdir {output}"], shell=True)
else:
    raise RuntimeError(f"Output folder '{output}' already exists. Aborting to avoid overwriting.")


print('Process 1: Embedding')
start =  time.time()
id, data = u.data_embedding_vector(input, batch_size, device)
torch.save(data, f'{output}/data.pt')
end =  time.time()
time_str = u.consumption_time(start, end)
print(f'Process 1: Embedding finished. Consumption time: {time_str}.')

print('Process 2: IsEnzyme prediction')
start =  time.time()
net = md.MLP_esm3B_1()
net.load_state_dict(torch.load("./models_param/esm2_3B_MLP_1.pt"))
probs, labels = u.prediction(net, data, batch_size, device)
probs = [np.array(n).max() for i in probs for n in i]
labels = np.concatenate(labels).tolist()
IsEnzyme = u.transfer_name(labels, u.enzyme_dic)
result_IsEnzyme = pd.DataFrame({"id": id, "IsEnzyme":IsEnzyme, "prob":probs})
result_IsEnzyme.to_csv(f'{output}/result_IsEnzyme.csv',index=False)
end =  time.time()
time_str = u.consumption_time(start, end)
print(f'Process 2: IsEnzyme prediction finished. Consumption time: {time_str}.')

print('Process 3: First EC prediction')
start =  time.time()
bool_list = result_IsEnzyme['IsEnzyme'] == 'enzyme'
bool_tensor = torch.tensor(bool_list)
enzyme_data = data[bool_tensor]
net = md.MLP_esm3B_2()
net.load_state_dict(torch.load("./models_param/esm2_3B_MLP_2.pt"))
probs, labels = u.prediction(net, enzyme_data, batch_size, device)
probs = [np.array(n).max() for i in probs for n in i]
labels = np.concatenate(labels).tolist()
First_EC = u.transfer_name(labels, u.first_enzyme_dic)
result_IsEnzyme.loc[bool_list,'IsEnzyme'] = First_EC
result_IsEnzyme.loc[bool_list,'prob'] = probs
result_First_EC = result_IsEnzyme
result_First_EC.rename(columns={'IsEnzyme':'First_EC'}, inplace=True)
result_First_EC.to_csv(f'{output}/result_First_EC.csv',index=False, columns=['id','First_EC','prob'])
end =  time.time()
time_str = u.consumption_time(start, end)
print(f'Process 3: First EC prediction finished. Consumption time: {time_str}.')

print('Process 4: Second EC prediction')
start =  time.time()
First_EC_value = result_First_EC['First_EC'].tolist()
Second_EC_value = []
Second_EC_prob = []
for i, n in enumerate(First_EC_value):
    if n == 'non_enzyme':
        Second_EC_value.append(n)
        Second_EC_prob.append(0)
    elif n == '1':
        net = md.MLP_esm3B_3()
        net.load_state_dict(torch.load("./models_param/esm2_3B_MLP_3.pt"))
        Second_EC_data = data[[i]]
        probs, labels = u.prediction(net, Second_EC_data, batch_size, device)
        probs = [np.array(na).max() for ia in probs for na in ia]
        labels = labels[0]
        Second_EC_value.append(u.transfer_name(labels, u.second_enzyme_first_dic)[0])
        Second_EC_prob.append(probs[0])
    elif n == '2':
        net = md.MLP_esm3B_4()
        net.load_state_dict(torch.load("./models_param/esm2_3B_MLP_4.pt"))
        Second_EC_data = data[[i]]
        probs, labels = u.prediction(net, Second_EC_data, batch_size, device)
        probs = [np.array(na).max() for ia in probs for na in ia]
        labels = labels[0]
        Second_EC_value.append(u.transfer_name(labels, u.second_enzyme_second_dic)[0])
        Second_EC_prob.append(probs[0])
    elif n == '3':
        net = md.MLP_esm3B_5()
        net.load_state_dict(torch.load("./models_param/esm2_3B_MLP_5.pt"))
        Second_EC_data = data[[i]]
        probs, labels = u.prediction(net, Second_EC_data, batch_size, device)
        probs = [np.array(na).max() for ia in probs for na in ia]
        labels = labels[0]
        Second_EC_value.append(u.transfer_name(labels, u.second_enzyme_third_dic)[0])
        Second_EC_prob.append(probs[0])
    elif n == '4':
        net = md.MLP_esm3B_6()
        net.load_state_dict(torch.load("./models_param/esm2_3B_MLP_6.pt"))
        Second_EC_data = data[[i]]
        probs, labels = u.prediction(net, Second_EC_data, batch_size, device)
        probs = [np.array(na).max() for ia in probs for na in ia]
        labels = labels[0]
        Second_EC_value.append(u.transfer_name(labels, u.second_enzyme_fourth_dic)[0])
        Second_EC_prob.append(probs[0])
    elif n == '5':
        net = md.MLP_esm3B_7()
        net.load_state_dict(torch.load("./models_param/esm2_3B_MLP_7.pt"))
        Second_EC_data = data[[i]]
        probs, labels = u.prediction(net, Second_EC_data, batch_size, device)
        probs = [np.array(na).max() for ia in probs for na in ia]
        labels = labels[0]
        Second_EC_value.append(u.transfer_name(labels, u.second_enzyme_fiveth_dic)[0])
        Second_EC_prob.append(probs[0])
    elif n == '6':
        net = md.MLP_esm3B_8()
        net.load_state_dict(torch.load("./models_param/esm2_3B_MLP_8.pt"))
        Second_EC_data = data[[i]]
        probs, labels = u.prediction(net, Second_EC_data, batch_size, device)
        probs = [np.array(na).max() for ia in probs for na in ia]
        labels = labels[0]
        Second_EC_value.append(u.transfer_name(labels, u.second_enzyme_sixth_dic)[0])
        Second_EC_prob.append(probs[0])
    elif n == '7':
        net = md.MLP_esm3B_9()
        net.load_state_dict(torch.load("./models_param/esm2_3B_MLP_9.pt"))
        Second_EC_data = data[[i]]
        probs, labels = u.prediction(net, Second_EC_data, batch_size, device)
        probs = [np.array(na).max() for ia in probs for na in ia]
        labels = labels[0]
        Second_EC_value.append(u.transfer_name(labels, u.second_enzyme_seventh_dic)[0])
        Second_EC_prob.append(probs[0])
result_Second_EC = pd.DataFrame({"id": id, "Second_EC":Second_EC_value, "prob":Second_EC_prob})
result_Second_EC.to_csv(f'{output}/result_Second_EC.csv',index=False)
end =  time.time()
time_str = u.consumption_time(start, end)
print(f'Process 3: Second EC prediction finished. Consumption time: {time_str}.')

print('All processes finished!')