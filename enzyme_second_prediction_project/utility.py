import esm
from esm import model
import torch.nn as nn
import torch 
import time
from Bio import SeqIO
import os
from torch.utils import data
import numpy as np

def data_embedding_vector(file, batch_size, device):
    data = []
    record_id = []
    records = SeqIO.parse(file, "fasta")
    for record in records:
        record_id.append(record.id)
        data.append((record.id, record.seq))

    model, alphabet = esm.pretrained.esm2_t36_3B_UR50D()    
    batch_converter = alphabet.get_batch_converter()
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_tokens = batch_tokens[:,:1024]

    model = model.to(device)
    if device != 'cpu' and torch.cuda.device_count() > 1:
        print("Let's use", torch.cuda.device_count(), "GPUs!")
        model = nn.DataParallel(model)
    batch_tokens = batch_tokens.split(batch_size)
    sequence_representations = []
    model.eval()
    with torch.no_grad():
        for batch_token in batch_tokens:
            batch_token = batch_token.to(device)            
            results = model(batch_token, repr_layers=[36], return_contacts=False)
            token_representations = results["representations"][36].cpu()
            sequence_representation = []
            for i in range(token_representations.shape[0]):
                seq_len = 1024 - ((batch_token[i] == 1) | (batch_token[i] == 2)).sum()
                sequence_representation.append(token_representations[i, 1 : min(seq_len,1024)].mean(0).reshape(1,-1))
            sequence_representations.append(torch.cat(sequence_representation, dim=0))
    final_result = torch.cat(sequence_representations, dim=0) 
    return record_id, final_result


def prediction(net, data, batch_size, device):
    iters = torch.split(data, split_size_or_sections=batch_size, dim=0)
    model = net.to(device)
    if device != 'cpu' and torch.cuda.device_count() > 1:
        print("Let's use", torch.cuda.device_count(), "GPUs!")
        model = nn.DataParallel(model)
    model.eval()
    with torch.no_grad():
        labels = []
        probs = []
        for X in iters:
            X = X.to(device)
            preds = net(X).detach()
            preds = preds.to('cpu')
            prob = nn.Softmax(dim=1)
            preds_prob = prob(preds)
            preds_label = np.argmax(preds_prob.numpy(),axis=1)
            probs.append(preds_prob)
            labels.append(preds_label)
    return  probs, labels

def consumption_time(start, end):
    run_time = round(end-start)
    hour = run_time // 3600
    minute = (run_time - 3600 * hour) // 60
    second = run_time - 3600 * hour - 60 * minute
    time = [hour, 'h', minute, 'm', second, 's']
    time_str = ''.join(map(str, time))
    return time_str

def transfer_name(data, dic):
    replaced_list = [dic[item] for item in data]
    return replaced_list

enzyme_dic = {0: 'non_enzyme', 1: 'enzyme'}
first_enzyme_dic = {0: "1", 
                    1: "2",
                    2: "3",
                    3: "4",
                    4: "5",
                    5: "6",
                    6: "7"}
second_enzyme_first_dic = { 0: '1.1',
                            1: '1.2',
                            2: '1.3',
                            3: '1.4',
                            4: '1.5',
                            5: '1.6',
                            6: '1.7',
                            7: '1.8',
                            8: '1.9',
                            9: '1.10',
                            10: '1.11',
                            11: '1.12',
                            12: '1.13',
                            13: '1.14',
                            14: '1.15',
                            15: '1.16',
                            16: '1.17',
                            17: '1.18',
                            18: '1.21',
                            19: '1.97',
                            20: '1.others'} # 将 '1.20': 63, '1.23': 12这两类归为others
second_enzyme_second_dic = {0: '2.1',
                            1: '2.2',
                            2: '2.3',
                            3: '2.4',
                            4: '2.5',
                            5: '2.6',
                            6: '2.7',
                            7: '2.8',
                            8: '2.others'} # 将 '2.9': 241, '2.10': 25这两类归为others

second_enzyme_third_dic = { 0: '3.1',
                            1: '3.2',
                            2: '3.4',
                            3: '3.5',
                            4: '3.6',
                            5: '3.others'} # 将 '3.3': 54, '3.7': 194, '3.8': 56, '3.9': 29, '3.10': 1, '3.11': 71, '3.12': 2, '3.13': 263, 这八类归为others
second_enzyme_fourth_dic = {0: '4.1',
                            1: '4.2',
                            2: '4.3',
                            3: '4.4',
                            4: '4.6',
                            5: '4.98',
                            6: '4.others'} # 将 '4.5': 6, '4.7': 2, '4.8': 7, '4.99': 292这三类归为others
second_enzyme_fiveth_dic = {0: '5.1',
                            1: '5.2',
                            2: '5.3',
                            3: '5.4',
                            4: '5.6',
                            5: '5.others'} # 将 '5.5': 180, '5.99': 59这两类归为others
second_enzyme_sixth_dic = { 0: '6.1',
                            1: '6.2',
                            2: '6.3',
                            3: '6.5',
                            4: '6.7',
                            5: '6.others'} # 将 '6.4': 78, '6.6': 50这两类归为others
second_enzyme_seventh_dic ={0: '7.1',
                            1: '7.2',
                            2: '7.3',
                            3: '7.4',
                            4: '7.5',
                            5: '7.6'}




    


