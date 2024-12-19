# this will contain the code to create msa alignment between component claseses it will take a list of components
from functools import partial

import biotite.sequence.align as align
from biotite.sequence import ProteinSequence
from biotite.structure import superimpose, rmsd, rmspd
from biotite.structure.io.pdb import PDBFile

matrix=align.SubstitutionMatrix.std_protein_matrix()

class SequenceAlignment:
    def __init__(self, matrix=matrix, **kwargs):
        self.matrix = matrix
        self.aligner=partial(align.align_multiple(matrix=matrix, **kwargs))

    def __call__(self, components):
        sequences=[]
        for component in components:
            sequences.append(str(component.sequence))
        self.alignments=self.aligne(sequences=sequences)
        self.distances=self.alignments[3]
        return self

class StructureAlignment:
    # this calculates rmsd between structures this can be a whole structure or a part
    def __init__(self, structure1, structure2):
        self.structure1=structure1.get_structure()
        self.structure2=structure2.get_structure()


    def align(self, **kwargs):
        aligned, transformation = superimpose(self.structure1, self.structure2, **kwargs)
        self.aligned=aligned
        self.transformation=transformation
        return self

    def rmsd(self):
        struct_rmsd=rmsd(self.structure1, self.structure2)
        self.rmsd=struct_rmsd
        return self

