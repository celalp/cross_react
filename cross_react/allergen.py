import os

import pysam
import pandas as pd
from Bio import SeqIO
from biotite.structure.io.pdb import PDBFile

from .blast import Blast
from .component import Component

parser=PDBParser()

class Allergen:
    def __init__(self, name,  path, blast_db=None, components=None, proteome=None):
        """
        This class represent the whole allergen (e.g. peanut), it contains organism level information as well as all the
        componenets (see component class).
        :param name: name of the allergen
        :param latin_name: latin name of the allergen
        :param short_name: short name of the allergen
        :param path: path for the folder that contains all the files for that allergen
        :param blast_db: blast db created from the proteome
        :param components: list of Component objects
        :param proteome: fasta file for the proteome
        """
        self.name = name
        self.path = path
        self.proteome=proteome
        self.components=components
        self.blast_db=blast_db
        self.reduced_components=None


    def get_components(self, fasta, esm2, esm3, structures):
        components = []
        with open(fasta, "rb") as f:
            for record in SeqIO.parse(f, "fasta"):
                name = record.id
                structure = os.path.join(structures, name + ".pdb")
                esm2 = esm2[name]
                esm3 = esm3[name]
                components[name] = Component(name, structure, esm2, esm3)
        return components

    def build_from_dir(self):
        contents = os.listdir(self.path)
        self.proteome=pysam.FastaFile(os.path.join(self.path, "protein.faa"))

        if "blast" in contents:
            self.blast_db=Blast(db=os.path.join(self.path, "blast/"), db_name="blast", dbtype="p")
        else:
            if self.proteome is not None:
                db_dir=os.makedirs(db=os.path.join(self.path, "blast/"), db_name="blast", exist_ok=True)
                blast = Blast(db_dir, db=db_dir, dbtype="p")
                blast.create_db(fasta=self.proteome, output_path=os.path.join(db_dir, 'blast'),
                            db_name="blast")
            self.blast_db = blast

        component_fasta=os.path.join(self.path, "allergens.fa")
        if "allergens_uniq.fa" in contents:
            reduced_components_fasta=os.path.join(self.path, "allergens_uniq.fa")
        else:
            reduced_components_fasta=None

        structures_folder=os.path.join(self.path, "allergen_structures")
        esm2_embeddings=pd.read_pickle(os.path.join(self.path, "esm2_embeddings.pkl"))
        esm3_embeddings=pd.read_pickle(os.path.join(self.path, "esm3_embeddings.pkl"))

        self.components={}
        with open(component_fasta, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                name=record.id
                sequence=record.seq
                structure = PDBFile.read(os.path.join(structures_folder, name+".pdb"))
                esm2=esm2_embeddings[name]
                esm3=esm3_embeddings[name]
                self.components[name]=Component(name=name, sequence=sequence,
                                                structure=structure, esm2_embed=esm2,
                                                esm3_embed=esm3)
        #TODO maybe a mask
        self.reduced_components = {}
        if reduced_components_fasta is not None:
            with open(reduced_components_fasta, "r") as f:
                for record in SeqIO.parse(f, "fasta"):
                    name=record.id
                    structure=parser.get_structure(name, os.path.join(structures_folder, name+".pdb"))
                    sequence = record.seq
                    esm2=esm2_embeddings[name]
                    esm3=esm3_embeddings[name]
                    self.reduced_components[name] = Component(name=name, sequence=sequence,
                                                      structure=structure, esm2_embed=esm2,
                                                      esm3_embed=esm3)

    def save(self, path):
        pd.to_pickle(self, path)

    def __str__(self):
        "{} with {} components".format(self.name, len(self))

    def __len__(self):
        print(len(self.components.keys()))





