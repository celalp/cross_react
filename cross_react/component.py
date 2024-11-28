

class Component:
    def __init__(self, name, sequence, structure=None, esm2_embed=None,
                 esm3_embed=None, canonical=False, homologs=None, full_seq=None):
        self.name = name
        self.sequence = sequence
        self.structure = structure
        self.esm2_embed = esm2_embed
        self.esm3_embed = esm3_embed
        self.canonical = canonical
        self.homologs = homologs
        self.full_seq = full_seq

    def get_full_seq(self, blast_db):
        pass

    def get_homolog(self, blast_db):
        pass

    def get_esm2_embed(self, pipeline):
        pass

    def get_esm3_embed(self, pipeline):
        pass

    def get_foldseeek_seq(self, foldseek):
        pass

    def sasa(self):
        pass
