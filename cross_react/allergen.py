import os
from blast import Blast



class Allergen:
    def __init__(self, name, latin_name, short_name, path, blast_db=None, components=None, proteome=None):
        """

        :param name:
        :param latin_name:
        :param short_name:
        :param path:
        :param blast_db:
        :param components:
        :param proteome:
        """
        self.name = name
        self.short_name = short_name
        self.latin_name=latin_name
        self.path = path
        self.blast_db=blast_db
        self.proteome=proteome
        self.components=components

    def create_db(self):
        pass

    def reduce(self):
        pass

    def __str__(self):
        print("{} with latin name {}".format(self.name, self.latin_name))

    def __len__(self):
        print(len(self.components))





