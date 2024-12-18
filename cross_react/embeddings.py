
from Bio import SeqIO
from transformers import AutoTokenizer, AutoModel, pipeline
from esm.models.esm3 import ESM3
from esm.sdk.api import ESMProtein, SamplingConfig
from esm.utils.structure.protein_chain import ProteinChain
import pandas as pd
import torch



class Embeddings:
    def __init__(self, model, tokenizer, device="cuda", **kwargs):
        self.model = AutoModel.from_pretrained(model, **kwargs)
        self.tokenizer = AutoTokenizer.from_pretrained(tokenizer, **kwargs)

        if device is None:
            self.device = "cuda" if torch.cuda.is_available() else "cpu"
        else:
            self.device = device

        self.esm2 = pipeline("feature-extraction", framework="pt", model=self.model,
                                          tokenizer=self.tokenizer, device=self.device)

        #there is only one option here
        self.esm3=ESM3.from_pretrained("esm3-open").to(device=self.device)

    #TODO this is not correct
    def get_esm2_embeddings(self, fasta):
        embeddings = {}
        sequences={}

        for record in SeqIO.parse(fasta, "fasta"):
            sequences[record.id] = str(record.seq)

        record_names = list(sequences.keys())
        for i in range(len(record_names)):
            embeddings[record_names[i]]= self.esm2(sequences[record_names[i]], return_tensors = "pt")
        return embeddings

    def get_esm3_embeddings(self, pdbs, names):
        if len(names)!=len(pdbs):
            raise ValueError("Length of names and pdbs do not match")

        embeddings = {}

        for name, struct in zip(pdbs, names):
            protein_chain = ProteinChain.from_pdb(struct)
            protein = ESMProtein.from_protein_chain(protein_chain)
            protein_tensor = self.esm3.encode(protein)
            output = self.esm3.forward_and_sample(protein_tensor, SamplingConfig(return_per_residue_embeddings=True))
            embeddings[name] = output.per_residue_embedding

        return embeddings

    def split_features(self, embeddings, window=8, stride=1, dim=0):
        split_features = embeddings.unfold(dim, window, stride).transpose(2, 1)
        split_features = torch.mean(split_features, 1, False)
        return split_features

    def compare_split_features(self, feat1, feat2, order=2):
        if order == 2:
            distances = torch.cdist(feat1, feat2, p=order, compute_mode='donot_use_mm_for_euclid_dist')
        else:
            distances = torch.cdist(feat1, feat2, p=order)
        return distances

    def save(self, features, dest):
        pd.to_pickle(features, dest)


