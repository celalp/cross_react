import os

from Bio import Seq
from biotite.sequence import sasa


class Component:
    def __init__(self, name, sequence, structure=None, esm2_embed=None,
                 esm3_embed=None):
        self.name = name
        self.sequence = sequence
        self.structure = structure
        self.esm2_embed = esm2_embed
        self.esm3_embed = esm3_embed
        self.homologs={}

    def get_full_seq(self, allergen):
        blast_results=allergen.blast_db.Blast(self.sequence)
        blast_results=blast_results.sort_values("bitscore", ascending=False)
        best_match=blast_results["saccver"][blast_results["bitscore"]==blast_results["bitscore"].max()].to_list()

        self.proteome_sequence=[]
        for item in best_match:
            self.proteome_sequence.append(Seq.Seq(allergen.proteome.fetch(item)))

        return self

    def get_homologs(self, allergen, return_seq=True):
        blast_results = allergen.blast_db.Blast(self.sequence)
        blast_results = blast_results.sort_values("bitscore", ascending=False)
        best_match = blast_results["saccver"][blast_results["bitscore"] == blast_results["bitscore"].max()].to_list()

        if return_seq:
            homologs={}
            for item in best_match:
                homologs[item]=(Seq.Seq(allergen.proteome.fetch(item)))
        else:
            homologs=best_match

        self.homologs[allergen.name]=homologs
        return self

    def get_sasa(self, **kwargs):
        structure=self.structure.get_structure()
        str_sasa=sasa(structure, **kwargs)
