import argparse
from enformer_pytorch import Enformer
import torch
import time
import numpy as np
from pyfaidx import Fasta
import glob
import matplotlib.pyplot as plt
import pandas as pd
import warnings
import os
warnings.filterwarnings("ignore")

SEQ_LEN = 196608
BIN_SIZE= 128
SEQ_LEN_DIV_BY2 = SEQ_LEN//2

one_hot_embed = torch.zeros(256, 4)
one_hot_embed[ord('a')] = torch.Tensor([1., 0., 0., 0.])
one_hot_embed[ord('c')] = torch.Tensor([0., 1., 0., 0.])
one_hot_embed[ord('g')] = torch.Tensor([0., 0., 1., 0.])
one_hot_embed[ord('t')] = torch.Tensor([0., 0., 0., 1.])
one_hot_embed[ord('n')] = torch.Tensor([0., 0., 0., 0.])
one_hot_embed[ord('A')] = torch.Tensor([1., 0., 0., 0.])
one_hot_embed[ord('C')] = torch.Tensor([0., 1., 0., 0.])
one_hot_embed[ord('G')] = torch.Tensor([0., 0., 1., 0.])
one_hot_embed[ord('T')] = torch.Tensor([0., 0., 0., 1.])
one_hot_embed[ord('N')] = torch.Tensor([0., 0., 0., 0.])
one_hot_embed[ord('.')] = torch.Tensor([0.25, 0.25, 0.25, 0.25])

def cast_list(t):
    return t if isinstance(t, list) else [t]
def torch_fromstring(seq_strs):
    batched = not isinstance(seq_strs, str)
    seq_strs = cast_list(seq_strs)
    np_seq_chrs = list(map(lambda t: np.fromstring(t, dtype = np.uint8), seq_strs))
    seq_chrs = list(map(torch.from_numpy, np_seq_chrs))
    return torch.stack(seq_chrs) if batched else seq_chrs[0]
def str_to_one_hot(seq_strs):
    seq_chrs = torch_fromstring(seq_strs)
    return one_hot_embed[seq_chrs.long()]
def get_res(one_hot, target_length, head = "mouse"): 
    return enformer(one_hot, target_length =target_length,head=head)

def get_enformer_prediction_for_one_gene(gene_id, gtf, seqs, SEQ_LEN=SEQ_LEN, offset=0): 
    SEQ_LEN = SEQ_LEN +offset 
    SEQ_LEN_DIV_BY2 = SEQ_LEN//2
    select = gtf[gtf.gene_id == gene_id]
    #print(select)
    chromo = select.chr.values[0]
    n_seq = len(seqs)
    start = select.start.values[0] - 1 # python index starts from 0
    end = select.end.values[0] - 1 # python index starts from 0
    if start > SEQ_LEN_DIV_BY2:
        cut_start = start - SEQ_LEN_DIV_BY2
        pred_loc = (SEQ_LEN_DIV_BY2 + 1)//BIN_SIZE
    else:
        cut_start = 0
        pred_loc  = (start + 1)//BIN_SIZE
    cut_end = min(n_seq, start + SEQ_LEN_DIV_BY2) 
    one_hot = str_to_one_hot(np.array(seqs)[cut_start:cut_end])
  #  print(one_hot.shape)
    target_length= (cut_end-cut_start+1)//BIN_SIZE
    pred = get_res(one_hot, target_length=target_length, head="mouse").detach().numpy()[0]
    print("target_length: ", target_length, " pred_loc", pred_loc, " pred shape:", pred.shape, 
          " pred(kidney 1486): ",  pred[pred_loc, 1486])
    print("pred (kidney 1486) +/-3: ", pred[pred_loc-3:pred_loc+3, 1486])
    return pred[pred_loc,:]

def get_enformer_prediction_per_chrome(chrome):
    gene_list = gtf[gtf.chr == chrome].gene_id.values[0:2]
    #print(chrome, gene_list)
    res = {}
    i=0
    for gene_id in gene_list:
        print(i, "------", gene_id)
        pred = get_enformer_prediction_for_one_gene(gene_id, gtf, seqs)
        res[gene_id] = pred
        i=i+1
    return res


parser = argparse.ArgumentParser()
parser.add_argument("chrome")
parser.add_argument("ref_path")


args = parser.parse_args()
chrome = args.chrome
ref_path =  args.ref_path

#chrome ="1"
#ref_path="/fastscratch/chenm/A_J.39.fa"

print("chr:",  chrome)
print("ref_path:", ref_path)
gtf_output_filepath = "/projects/compsci/vmp/USERS/chenm/mahoney/enformer/gtf_data/"
gtf_filepath = gtf_output_filepath + os.path.basename(ref_path).replace("fa", "gtf")
print("gtf_path:", gtf_filepath)
gtf= pd.read_csv(gtf_filepath)
seqs = Fasta(ref_path)[chrome]

enformer = Enformer.from_pretrained('EleutherAI/enformer-official-rough')
enformer.eval()
res = get_enformer_prediction_per_chrome(chrome)

output_filepath = "/projects/compsci/vmp/USERS/chenm/mahoney/enformer/results/"
output_file = output_filepath + os.path.basename(ref_path).replace(".fa", "_chr") + chrome + ".csv"
pd.DataFrame(res).transpose().to_csv(output_file)
