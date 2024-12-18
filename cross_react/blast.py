import json
import os
import subprocess as sub
import warnings
from shutil import which

import pandas as pd


class Blast:
    def __init__(self, path=None, db=None, db_name="blast", dbtype="n"):
        """
        initiate a Blast class instance
        :param path: path of the executable if none will check $PATH
        :type path: str
        :param db: path and name of the blast database if exists if not it can be created using create_db
        :type db: str
        :param dbtype: type of the database n for nucleotide p for protein
        :type dbtype: str
        """
        execs = ["blastn", "blastp", "blastx", "tblastn", "tblastx", "makeblastdb"]
        if path is not None:
            for ex in execs:
                full_path = os.path.join(path, ex)
                if not os.path.exists(full_path):
                    raise FileNotFoundError("There was a problem finding executable {} please check your blast "
                                            "installation".format(exec))
        else:
            for ex in execs:
                if which(ex) is None:
                    raise EnvironmentError("{} does not seem to be installed, have you added blast to your $PATH?")

        if db is not None:
            if not os.path.exists(db):
                raise FileNotFoundError("Your blast database path does not exist")
            else:
                self.db = db+"/"+db_name
                if dbtype == "n":
                    self.db_type = "n"
                elif dbtype == "p":
                    self.db_type = "p"
                else:
                    raise ValueError("You can only have a nucleotide 'n' or a protein 'p' database")
        else:
            self.db = None

    def create_db(self, fasta, output_path, dbname, overwrite=True, arg_dict=None):
        """
        create a blast databse and stor in self.db
        :param dbtype: database type n for nucleotide and p for protein
        :type dbtype: str
        :param fasta: path of the fasta file only fasta is implemented
        :type fasta: str
        :param output_path: output path for the database this is different from the databse name
        :type output_path: str
        :param dbname: database name so self.db will be output_path/dbname
        :type dbname: str
         :param overwrite: if there is already a self.db you can override this just edits the class instance value
        dooes not touch the databse
        :type overwrite: bool
        :param arg_dict: a dictionary of arguments, if left empty will use default values see blast documentation
        :type arg_dict: dict
        :return: nothing just puts the new database path in self.db after database creation
        :rtype: None
        """

        if not os.path.isfile(fasta):
            raise FileNotFoundError("{} does not exists".format(fasta))

        if self.db is not None and not overwrite:
            raise FileExistsError("There is already a database for this class instance you "
                                  "can create another instance")

        arguments = self.__parse_args__(arg_dict)
        dbname = os.path.join(output_path, dbname)

        command = "makeblastdb -dbtype {} -input_type fasta -in {} -out {} {}".format(self.dbtype, fasta, dbname, arguments)
        makedb = sub.Popen(command.split(" "), stdout=sub.PIPE, stderr=sub.PIPE)
        makedb_out, makedb_err = makedb.communicate()

        if makedb.returncode != 0:
            print(makedb_err)
        else:
            print(makedb_out)
            self.db = dbname
            self.db_type = type

    def BLAST(self, seq, output_type="tabular", exec="blastn", arg_dict=None, cols=None):
        """
        blast a sequence
        :param exec: which execuctable to use, blastn, blastp, blastx, tblastn or tblastx see blast documentation for details
        :type exec:
        :param output_type: tabular or json if tabular is selected you can include a list of column names to be added to the table
        see blast documentation for details
        :type output_type: str
        :param seq: Bio.Seq object if multiple queries are needed you can call this function multiple times
        :type seq: Bio.Seq
        :param arg_dict: keyword arguments for the respective blast program if None default values will be used see blast documentations for details
        :type arg_dict: dict
        :return: a pandas dataframe if tabular dict if json
        :rtype: pd.DataFrame or dict
        """

        if exec in ["blastn", "tblastn", "tblastx"] and self.db_type == "p":
            raise ValueError("You are trying to use a protein database for a query that needs nucleotide info")

        if exec in ["blastp", "blastx"] and self.db_type == "n":
            raise ValueError("You are trying to use a nucleotide database for a query that needs protein info")

        command = [exec, "-db", self.db]
        arguments = self.__parse_args__(arg_dict)
        command = command + arguments
        if output_type == "tabular":
            if cols is None:
                cols=["qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart",
                      "qend", "sstart", "send", "evalue", "bitscore"]

                command_cols = "-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore'"
            else:
                command_cols=" ".join(cols)
                command_cols="-outfmt"+"'"+"6 "+command_cols+"'"
            command.append(command_cols)

        if output_type == "json":
            command.append("-outfmt 15")



        fasta_seq=">query\n"+str(seq)
        with open("tmp.fasta", "w") as tmp_fasta:
            tmp_fasta.write(fasta_seq)
        command.append("-query tmp.fasta")

        command = " ".join(command)
        blast = sub.Popen(command, shell=True, stdin=sub.PIPE, stdout=sub.PIPE, stderr=sub.STDOUT)
        results, err = blast.communicate()

        if blast.returncode != 0:
            warnings.warn("Blast run resulted in an error please see the error in the output")
            return err
        else:
            if output_type=="tabular":
                parsed_results = self.__parse_output__(out_type=output_type, results=results,
                                                       cols=cols)

        os.remove("tmp.fasta")
        return parsed_results


    def __parse_args__(self, arg_dict):
        """
        take a dict of arguments to be appended to the blast subprocess see blast documentation for available features
        :param arg_dict: a dictionary of argument key is the flat and value is the value if no value is needed for the flag
        it can be a 0 length string or None type
        :type arg_dict: dict
        :return: a list of strings to be passed to subprocess.run
        :rtype: list
        """
        arguments = []
        if arg_dict is not None:
            for arg in arg_dict.keys():
                arguments = arguments + ["-" + arg, arg_dict[arg]]
        else:
            arguments=[""]
        return arguments

    def __parse_output__(self, out_type, results, cols=None):
        """
        parse blast output, this depends on the out_type which there are several
        :param out_type: the kind of blast output
        :type out_type str
        :param results: subprocess PIPE output
        :type results: bytes
        :return: depends on the inputs
        :rtype:
        """
        if out_type == "tabular":
            parsed=pd.DataFrame([item.split("\t") for item in results.decode().split("\n")][:-1])
            parsed.columns=cols
        elif out_type == "json":
            parsed = json.loads(results)

        return parsed
