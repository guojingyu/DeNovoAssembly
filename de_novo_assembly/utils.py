#!/usr/bin/env python

"""
This script is to provide utility functions, such as functions for text
based file IO

Author : Jingyu Guo
"""

from Bio import SeqIO
import os.path
import sys
import logging
import datetime

# Configure logging
timestamp = datetime.datetime.now().strftime("%Y-%m-%d")
logging.basicConfig(filename="../" + timestamp + "_denovo_assembly.log",
                    stream=sys.stdout,
                    level=logging.DEBUG)
stderrLogger=logging.StreamHandler()
stderrLogger.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
logger = logging.getLogger(__name__)
logger.addHandler(stderrLogger)

def get_sequences_from_fasta_file(file, include_reverse_complement= True):
    """
    read fasta file and return all relevant sequences
    :param file: fasta file
    :param include_reverse_complement: if to include the reverse complement
    sequences
    :return: a dictionary of Seq objects (biopython)
    """
    if os.path.isfile(file):
        try:
            record_dict = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
            logger.info(datetime.datetime.now().strftime("%H:%M:%S") +
                " Read the input file: " + file)
            if len(record_dict) == 0:
                logger.error(clock_now() +
                             "Read : no sequence input in file: " + file)
                raise RuntimeError("Read : no sequence input in file: " + file)
            logger.info(datetime.datetime.now().strftime("%H:%M:%S") +
                        " A number of " + str(len(record_dict.keys())) +
                        " records read in.")
            return record_dict
        except IOError as err:
            print(err.args)
            print("Reading fasta file failed. Please check file :" + file)
            logger.error(clock_now() + " FileRead: " + err.message)
            raise err

    else:
        raise RuntimeError("The provided file not exist :" + file)

def include_reverse_complement(record_dict):
    """
    the methods to add reverse complement sequences into to-be-assembled
    sequences
    :param record_dict: dictionary contains fasta reads
    :return: a dictionary containing also the reverse complement sequences
    """
    reversed_record_dict = {}
    for seq_id, seq in record_dict.iteritems():
        try:
            reversed_record_dict[seq_id + "_rev_com"] = seq.reverse_complement()
        except Exception as exp:
            print "Reverse Complement : failed to convert the sequence: " + seq_id
            logger.error(clock_now() + " Reverse Complement : failed to convert the sequence: " + seq_id)
            logger.error(clock_now() + " Reverse Complement :" + exp.message)
            # skip and not raised

    logger.info(clock_now() +
                " A number of " + str(len(reversed_record_dict.keys())) +
                " reverse complement records added.")
    record_dict.update(reversed_record_dict)
    return record_dict

def output_assembly(sequences, output_file_name="./output.txt",
                    longest_assembly_only=False):
    """
    Output function
    :param sequences_in_list: a list of strings as sequences
    :param output_file_name_base: an output file path and name
    :param longest_assembly_only: true or false
    :return:
    """
    try:
        with open(output_file_name, "w") as output_handle:
            if longest_assembly_only:
                sequence = max(sequences, key=len)
                output_handle.write(sequence + '\n')
            else:
                for sequence in sorted(sequences,key=len,reverse=True):
                    output_handle.write(sequence + '\n')

        logger.info(clock_now() +
                    " Writing to output file " + output_file_name +
                    " successfully.")

    except IOError as err:
        print("Writing to output file failed. Please check path :" +
              output_file_name)
        logger.error(clock_now() + " Write output : failed to write to "
                                   "output." + output_file_name + "\n" + err.message)
        raise err
    except RuntimeError as err:
        print("Writing to output file failed"+
                     err.message)
        logger.error(clock_now() + " Write output : failed to write to output." +
                     err.message)
        raise err


def clock_now():
    return datetime.datetime.now().strftime("%H:%M:%S")


if __name__ == "__main__":
    # run the file reading
    fa_test_file = "../data/coding_challenge_data_set.fasta"
    get_sequences_from_fasta_file(fa_test_file)
