# This class uses the Altschul-Erickson method for dinucleotide shuffling
# P. Clote, Oct 2003
## Code is derived from From github.com/wassermanlab/BiasAway
#Original publication: http://mbe.oxfordjournals.org/content/2/6/526.abstract
import sys,string,random,re
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import os


def dinuc_shuffle_fasta(input_file,output_file=''):
    """
    Dinucleotide shuffles an input fasta file
    If an output_file is specified, the shuffled file
    will be written to the output file.
    Else, output will go to stdout
    """
    
    seq_parser = SeqIO.parse(input_file,"fasta")
    if output_file != '':
        #abs_output_file = os.path.join(os.path.dirname(input_file),output_file)
        #print abs_output_file
        output_file_handle =open(output_file,'w')
        
        
    for record in seq_parser:
        #Writes dinucleotide shuffled version of each input sequence
        if output_file == '':
            generate_sequences(seq_parser,1,None)
        else:
            generate_sequences(seq_parser,1,output_file_handle)
            output_file_handle.close()
    seq_parser.close()
    
    return output_file



def computeCountAndLists(s):
  #WARNING: Use of function count(s,'UU') returns 1 on word UUU
  #since it apparently counts only nonoverlapping words UU
  #For this reason, we work with the indices.

  # NOTE: One cannot use function "count(s,word)" to count the number
  # of occurrences of dinucleotide word in string s, since the built-in
  # function counts only nonoverlapping words, presumably in a left to
  # right fashion.

  
  #Initialize lists and mono- and dinucleotide dictionaries
  List = {} #List is a dictionary of lists
  List['A'] = []; List['C'] = [];
  List['G'] = []; List['T'] = [];
  nuclList   = ["A","C","G","T"]
  s       = s.upper()
  s       = s.replace("T","T")
  nuclCnt    = {}  #empty dictionary
  dinuclCnt  = {}  #empty dictionary
  for x in nuclList:
    nuclCnt[x]=0
    dinuclCnt[x]={}
    for y in nuclList:
      dinuclCnt[x][y]=0

  #Compute count and lists
  nuclCnt[s[0]] = 1
  nuclTotal     = 1
  dinuclTotal   = 0
  for i in range(len(s)-1):
    x = s[i]; y = s[i+1]
    List[x].append( y )
    nuclCnt[y] += 1; nuclTotal  += 1
    dinuclCnt[x][y] += 1; dinuclTotal += 1
  assert (nuclTotal==len(s))
  assert (dinuclTotal==len(s)-1)
  return nuclCnt,dinuclCnt,List
 
 
def chooseEdge(x,dinuclCnt):
  numInList = 0
  for y in ['A','C','G','T']:
    numInList += dinuclCnt[x][y]
  z = random.random()
  denom=dinuclCnt[x]['A']+dinuclCnt[x]['C']+dinuclCnt[x]['G']+dinuclCnt[x]['T']
  numerator = dinuclCnt[x]['A']
  if z < float(numerator)/float(denom):
    dinuclCnt[x]['A'] -= 1
    return 'A'
  numerator += dinuclCnt[x]['C']
  if z < float(numerator)/float(denom):
    dinuclCnt[x]['C'] -= 1
    return 'C'
  numerator += dinuclCnt[x]['G']
  if z < float(numerator)/float(denom):
    dinuclCnt[x]['G'] -= 1
    return 'G'
  dinuclCnt[x]['T'] -= 1
  return 'T'



def connectedToLast(edgeList,nuclList,lastCh):
  D = {}
  for x in nuclList: D[x]=0
  for edge in edgeList:
    a = edge[0]; b = edge[1]
    if b==lastCh: D[a]=1
  for i in range(2):
    for edge in edgeList:
      a = edge[0]; b = edge[1]
      if D[b]==1: D[a]=1
  ok = 0
  for x in nuclList:
    if x!=lastCh and D[x]==0: return 0
  return 1
 

def eulerian(s):
  nuclCnt,dinuclCnt,List = computeCountAndLists(s)
  #compute nucleotides appearing in s
  nuclList = []
  for x in ["A","C","G","T"]:
    if x in s: nuclList.append(x)
  #compute numInList[x] = number of dinucleotides beginning with x
  numInList = {}
  for x in nuclList:
    numInList[x]=0
    for y in nuclList:
      numInList[x] += dinuclCnt[x][y]
  #create dinucleotide shuffle L 
  firstCh = s[0]  #start with first letter of s
  lastCh  = s[-1]
  edgeList = []
  for x in nuclList:
    if x!= lastCh: edgeList.append( [x,chooseEdge(x,dinuclCnt)] )
  ok = connectedToLast(edgeList,nuclList,lastCh)
  return ok,edgeList,nuclList,lastCh



def shuffleEdgeList(L):
  n = len(L); barrier = n
  for i in range(n-1):
    z = int(random.random() * barrier)
    tmp = L[z]
    L[z]= L[barrier-1]
    L[barrier-1] = tmp
    barrier -= 1
  return L



def dinuclShuffle(s):
  ok = 0
  while not ok:
    ok,edgeList,nuclList,lastCh = eulerian(s)
  nuclCnt,dinuclCnt,List = computeCountAndLists(s)

  #remove last edges from each vertex list, shuffle, then add back
  #the removed edges at end of vertex lists.
  for [x,y] in edgeList: List[x].remove(y)
  for x in nuclList: shuffleEdgeList(List[x])
  for [x,y] in edgeList: List[x].append(y)

  #construct the eulerian path
  L = [s[0]]; prevCh = s[0]
  for i in range(len(s)-2):
    ch = List[prevCh][0] 
    L.append( ch )
    del List[prevCh][0]
    prevCh = ch
  L.append(s[-1])
  t = string.join(L,"")
  return t



def generate_sequences(seqs, nfold, output_file_handle=None):
    # Module allowing the generation of sequences by using a di-nucleotide
    # shuffling of the given sequences
    cpt = 1
    #bg_gc_list = []
    bg_lengths = []
    for record in seqs:
        seq = record.seq.__str__()
        descr = ''
        #descr = "Background sequence for {0:s}".format(record.name)
        for n in range(0, nfold):
            new_sequence = ""
            for sequence in re.split('(N+)', seq):
                if re.match('N', sequence):
                    new_sequence += sequence
                elif sequence:
                    new_sequence += dinuclShuffle(sequence)
            new_seq = SeqRecord(Seq(new_sequence, generic_dna),
                            id="dinucleotide_shuffled_seq_{0:s}".format(record.name),
                                        description=descr)
            if (output_file_handle == None):
                print new_seq.format("fasta"),
                bg_lengths.append(len(new_sequence))
            else:
                output_file_handle.write(new_seq.format("fasta"))
                bg_lengths.append(len(new_sequence))
            cpt += 1
        return bg_lengths
