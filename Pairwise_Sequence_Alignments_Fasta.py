## import necessary python modules
from Bio import SeqIO, pairwise2
#Biopython modules to use for sequences and alignments
import os, sys, re, getopt, itertools, glob
#import system operator modules to enable the use of certain functions
from glob import glob
#import glob function from the glob module

def pairwise_to_dict(file):
        '''Takes a fasta file containing two sequences and makes a single "best" alignment, which is saved into a
           dictionary. The dictionary key is a combination of both fasta headers.'''
        aln_dict = {}
        #create an empty dictionary
        sequences = list(SeqIO.parse(file, "fasta"))
        #create a variable and parse into a list fasta sequences. this function takes in two parameters. the first parameter is the sequence
        #the second parameter indicates what type of file format which is "fasta" in this case
        for aln in pairwise2.align.globalxx(sequences[0].seq, sequences[1].seq, one_alignment_only = True):
            #this function uses a for loop to align sequences in a pairwise format using the pairwise function from Biopython
            #this global match function doesnt penalises gaps. for identical matches, it gives a score of 1 and 0 for non-identical (xx)
            aln_dict[sequences[0].id + '_' + sequences[1].id] = (pairwise2.format_alignment(*aln))
            #and only gives the only one best possible alignment for the two sequences
        return aln_dict
    #return the dictionary which should contain both headers and an alignemnt for two sequences and also a score of every key of the dictionary
        


## start of the main options block
def main():
"""Documentations... defining friendly user iinterface for use"""
    workingdirectory = 0
    #set variable to 0
    pathname = 0 
    #set variable to 0
        
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hw:p:', ['help=', 'workingdirectory=','pathname=',])
        #try user options 
    except getopt.GetoptError:
    #in case of an error
        usage()
        sys.exit(2)
    
## define all optionality for the script to allow for user input
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit(2)
        elif opt in ('-w', '--workingdirectory'):
            workingdirectory = arg
        elif opt in ('-p', '--pathname'):
            pathname = arg
        else:
            usage()
            sys.exit(2)
            
            
            
    ## define protocols if information required is not given
    ## for example, if 'primer seq' is not given (i.e. set to 0), a default sequence is taken.
    ## A message is written to stdout to tell the user.
    
    if workingdirectory == 0:
    #if variable is set to 0
        workingdirectory = os.getcwd()
        #assign the variable to the current working directory
        sys.stderr.write("Using default current directory %s \n" % workingdirectory)  
        #write stdout message to inform user
        

    if pathname == 0:
    #if variable is set to 0
        pathname = os.getcwd() + r"/Documents/Problem_2"
        #assign variable the current working with the name of folder
        sys.stderr.write("Using default current directory and foldername %s \n" % pathname)
        #write stdout message to inform user if incorrect

   
    os.chdir(pathname)
    
    
    
    for a, b in itertools.combinations(glob("*.fasta"), 2):
        if not "*combined*" in a and not "*combined*" in b:
            con_cat_test = os.path.splitext(a)[0] + "_" + os.path.splitext(b)[0] + "_combined.fasta"
            if not os.path.exists(con_cat_test):
                os.system("cat " + a + " " + b + " > " + os.path.splitext(a)[0] + "_" + os.path.splitext(b)[0] + "_combined.fasta")
    
   
    
    
    with open('Sequence_To_PairwiseAlignment.csv', 'w') as f_out:
        f_out.write('headers,score,dif_in_bp\n')
        for files in glob("*combined.fasta"):
            #using the for loop
            #search for combined.fasta files using glob indicating the path directory where the combined.fasta is
            for k, v in pairwise_to_dict(files).items():
                #parse the keys and values of the glob files into the function pairwise_to_dict
                print(k, '\n', v, '\n')
                head_a, head_b = k.split('_')[0], k.split('_')[1]
                seq_a, match, seq_b = v.split('\n')[0],v.split('\n')[1],v.split('\n')[2]
                match_count = match.count('|')
                f_out.write("{0},{1},{2}\n".format(k, match_count, (len(seq_a)- match_count)))
                #print("Sequence {0} varies with sequence {1} by {2} bp.".format(head_a, head_b, (len(seq_a)- match_count)))
            
            
            
        
    os.chdir("..")
    os.chdir("..")
## start of help menu
def usage():  
    print(""" 
______________________________________________________________________________
                        Pairwise_Sequence_Alignments_Fasta
                               by Kevin A Agyeman
______________________________________________________________________________
G.
==============================================================================
-h	: help
-w	: working directory (current)
-p	: pathname which is the working directory and the foldername      

______________________________________________________________________________""")



    
if __name__ == "__main__":
    main()





