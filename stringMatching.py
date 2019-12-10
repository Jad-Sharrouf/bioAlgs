import urllib.request
import math
import numpy as np
import matplotlib.pyplot as plt
import os.path
import time
from operator import attrgetter
import matplotlib.pyplot as plt

#FileName = "cullpdb_pc90_res3.0_R1.0_d191010_chains39545.gz"
FileName = input("Enter the path to the text file containing the list of protein IDs & their chains: ")

#each line is a pattern
patternFile = input("Enter the path to the text file containing the patterns: ")
patterns = open(patternFile,"r").read().split()

myFile = open(FileName,'r')
lines = myFile.readlines()

#I've included the file cullpdb_pc90_res3.0_R1.0_d191010_chains39545.gz that contains 39k proteins
#you can change the number of proteins down below (each line is a protein)
lines = lines[:5]

proteins = []
for line in lines:
    if len(line.strip()) != 0:
        proteins.append(line.split()[0])

proteins = list(filter(lambda a: len(a)<=5, proteins))

myFile.close()

#preprocessing time for boyer moore
preBM = 0
class last_occurrence(object):
    def __init__(self, pattern, alphabet):
        self.occurrences = dict()
        for letter in alphabet:
            self.occurrences[letter] = pattern.rfind(letter)

    def __call__(self, letter):
        return self.occurrences[letter]

    def __str__(self):
        return str(self.occurrences)

bmTime = 0
def boyer_moore(text, pattern):   
    global preBM
    global bmTime
    
    s = time.time()
    
    alphabet = set(text)
    last = last_occurrence(pattern, alphabet)
    e = time.time()
    preBM += (e-s)
    #print(last)
    
    s1 = time.time()
    m = len(pattern) 
    n = len(text)
    occurrences = []
    
    i = 0
    while(i <= n-m): 
        j = m-1
  
        while j>=0 and pattern[j] == text[i+j]: 
            j -= 1
  
        if j<0: 
            occurrences.append(i)
            i += (m-last(text[i+m]) if i+m<n else 1) 
        else: 
            i += max(1, j-last(text[i+j])) 
    
    e1 = time.time()
    bmTime += (e1-s1)
    return occurrences 

preKMP = 0
def prefix(P):
    m = len(P)
    pi = [None]*m
    pi[0] = 0
    
    i = 0
    j = 1
    while j<m:
        if P[i] == P[j]:
            pi[j] = i+1
            i += 1
            j += 1
        else:
            if i == 0:
                pi[j] = 0
                j += 1
            else:
                i = pi[i-1]
    return pi


kmpTime = 0
def KMP(T,P):
    global preKMP
    global kmpTime
    
    s = time.time()
    pi = prefix(P)
    e = time.time()
    preKMP += (e-s)
    
    s2 = time.time()
    occurrences = []

    n = len(T)
    m = len(P)
    
    i = 0
    j = 0
    while i < n:
        if P[j] == T[i]:
            i+=1
            j+=1
            
        if j == m:
            occurrences.append(i-j)
            j = pi[j-1]
        elif i < n and P[j] != T[i]:
            if j != 0:
                j = pi[j-1]
            else:
                i+=1
                
    e2 = time.time()
    kmpTime += (e2-s2)
    return occurrences

def Naive(T,P):
    occurrences = []
    
    n = len(T)
    m = len(P)
    
    for i in range(n-m+1): 
        j = 0
          
        while(j < m): 
            if (T[i + j] != P[j]): 
                break
            j += 1
  
        if (j == m):  
            occurrences.append(i)

    return occurrences
    
def suffixArrays(array, query, seq, lo=0, hi=None):
    if lo < 0:
        print(lo,'must be non-negative!')
    if hi is None: 
        hi = len(array)
    while lo < hi: 
        mid = (lo+hi)//2
        if seq[array[mid]:] < query:
            lo = mid+1
        else:
            hi = mid

    def match_at(i):
        return seq[i: i + len(query)] == query

    first = lo
    while first > 0 and match_at(array[first - 1]):
        first -= 1

    last = lo
    while match_at(array[last]):
        last += 1

    return array[first:last]


leafEnd = -1
class Node:

    def __init__(self, leaf):
        self.children = {}
        self.leaf = leaf
        self.suffixIndex = None
        self.start = None
        self.end = None
        self.suffixLink = None

    def __eq__(self, node):
        atg = attrgetter('start', 'end', 'suffixIndex')
        return atg(self) == atg(node)

    def __ne__(self, node):
        atg = attrgetter('start', 'end', 'suffixIndex')
        return atg(self) != atg(node)

    def __getattribute__(self, name):
        if name == 'end':
            if self.leaf:
                return leafEnd
        return super(Node, self).__getattribute__(name)


class SuffixTree:
    def __init__(self, data):
        self._string = data
        self.lastNewNode = None
        self.activeNode = None

        self.activeEdge = -1
        self.activeLength = 0

        self.remainingSuffixCount = 0
        self.rootEnd = None
        self.splitEnd = None
        self.size = -1 
        self.root = None

    def edge_length(self, node):
        return node.end - node.start + 1

    def walk_down(self, current_node):
        length = self.edge_length(current_node)
        if (self.activeLength >= length):
            self.activeEdge += length
            self.activeLength -= length
            self.activeNode = current_node
            return True
        return False

    def new_node(self, start, end=None, leaf=False):
        node = Node(leaf)
        node.suffixLink = self.root
        node.start = start
        node.end = end
        node.suffixIndex = -1
        return node

    def extend_suffix_tree(self, pos):
        global leafEnd

        leafEnd = pos
        self.remainingSuffixCount += 1

        self.lastNewNode = None
        while(self.remainingSuffixCount > 0):
            if (self.activeLength == 0):
                self.activeEdge = pos 
        
            if (self.activeNode.children.get(self._string[self.activeEdge]) is None):
                self.activeNode.children[self._string[self.activeEdge]] = self.new_node(pos, leaf=True)
               
                if (self.lastNewNode is not None):
                    self.lastNewNode.suffixLink = self.activeNode
                    self.lastNewNode = None

            else:
                _next = self.activeNode.children.get(self._string[self.activeEdge])
                if self.walk_down(_next): 
                    continue

                if (self._string[_next.start + self.activeLength] == self._string[pos]):
                    if((self.lastNewNode is not None) and (self.activeNode != self.root)):
                        self.lastNewNode.suffixLink = self.activeNode
                        self.lastNewNode = None
                    self.activeLength += 1
                    break
                    
                self.splitEnd = _next.start + self.activeLength - 1

                split = self.new_node(_next.start, self.splitEnd)
                self.activeNode.children[self._string[self.activeEdge]] = split

                split.children[self._string[pos]] = self.new_node(pos, leaf=True)
                _next.start += self.activeLength
                split.children[self._string[_next.start]] = _next
              
                if (self.lastNewNode is not None):
                    self.lastNewNode.suffixLink = split
                self.lastNewNode = split
            self.remainingSuffixCount -= 1
            if ((self.activeNode == self.root) and (self.activeLength > 0)):
                self.activeLength -= 1
                self.activeEdge = pos - self.remainingSuffixCount + 1
            elif (self.activeNode != self.root):
                self.activeNode = self.activeNode.suffixLink

    def walk_dfs(self, current):
        start, end = current.start, current.end
        yield self._string[start: end + 1]

        for node in current.children.values():
            if node:
                yield from self.walk_dfs(node)

    def build_suffix_tree(self):
        self.size = len(self._string)

        self.rootEnd = -1
        self.root = self.new_node(-1, self.rootEnd)
        self.activeNode = self.root
        for i in range(self.size):
            self.extend_suffix_tree(i)

    def print_dfs(self):
        for sub in self.walk_dfs(self.root):
            print(sub)
            
class Base:
    def __init__(self, tree):
        self.tree = tree
        self.main_string = tree._string
        self.root = tree.root


class FindOccurrence(Base):
    def __init__(self, tree, sub_string, findall=False):
        super(FindOccurrence, self).__init__(tree)
        self.sub_string = sub_string
        self.latest_index = 0
        self.findall = findall
        self.continue_flag = False
        self.sub_length = len(sub_string)

    def traverse(self, node, sub_string):
        if sub_string:
            item = next(((char, child) for char, child in node.children.items() if sub_string.startswith(char)), None)

            if item:
                char, child = item
                start, end = child.start, child.end
                if self.main_string[start: end + 1].startswith(sub_string):
                    if self.findall:
                        return self.find_all_match(child, len(sub_string))
                    return start - (self.sub_length - len(sub_string))
                return self.traverse(child, sub_string[end - start + 1:])
            else:
                return -1
        if self.findall:
            return self.find_all_match(node, len(sub_string))
        return node.start - (self.sub_length - len(sub_string))

    def check(self):
        if self.root is None:
            return -1
        if not isinstance(self.sub_string, str):
            return -1
        if not self.sub_string:
            return 0

        return self.traverse(self.root, self.sub_string)

    def find_all_match(self, node, sub_length):

        def inner(node, traversed_edges):
            for char, child in node.children.items():
                if child.leaf:
                    yield child.start - traversed_edges
                else:
                    start, end = child.start, child.end
                    sub_length = end - start + 1
                    yield from inner(child, traversed_edges + sub_length)

        if node.leaf:
            first = node.start - (self.sub_length - sub_length)
            return [first, *inner(node, self.sub_length)]
        else:
            return list(inner(node, self.sub_length))



naiveTime = 0
saTime = 0
preSATime = 0
stTime = 0
preSTTime = 0
for p in proteins:
    protein = p[:-1].lower()
    chainID = p[-1]
    
    url = 'https://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=fastachain&compression=NO&structureId='+ protein +'&chainId='+chainID
    response = ""
    try:
        response = urllib.request.urlopen(url)
    except Exception:
        continue
    
    data = response.read()    
    seq = data.decode('utf-8')
    
    seq = seq[seq.find('\n')+1:seq.rfind('\n')]
    seq = "".join(seq.split())
    
    ssa = time.time()
    suffix = [seq[i:] for i in range(len(seq))]
    suffix_sorted = sorted([seq[i:] for i in range(len(seq))])
    suffix_array = [suffix.index(ss) for ss in suffix_sorted]
    esa = time.time()
    preSATime += (esa-ssa)
    
    sst = time.time()
    T = seq + "$"
    tree = SuffixTree(T)
    tree.build_suffix_tree()
    est = time.time()
    preSTTime += (est-sst)
    
    for p in patterns:
        s0 = time.time()
        naive = Naive(seq, p) 
        e0 = time.time()
        naiveTime += (e0-s0)
        
        bm = boyer_moore(seq, p) 
        
        kmp = KMP(seq, p) 
        
        s3 = time.time()
        SA = []
        try:
            SA = suffixArrays(suffix_array, p, seq)
        except:
            pass
            
        e3 = time.time()
        saTime += (e3-s3)
        
        s4 = time.time()
        st = FindOccurrence(tree, p, findall=True).check()
        e4 = time.time()
        stTime += (e4-s4)
        
        #Printing the occurrence count and their indices
        print(p,"-->",len(naive),"occurrences.")
        for x in naive:
            print(x,end=" ")
        print("\n")

        print(p,"-->",len(bm),"occurrences.")
        for x in bm:
            print(x,end=" ")
        print("\n")

        print(p,"-->",len(kmp),"occurrences.")
        for x in kmp:
            print(x,end=" ")
        print("\n")

        print(p,"-->",len(SA),"occurrences.")
        for x in SA:
            print(x,end=" ")
        print("\n")

        if st != -1:
            print(p,"-->",len(st),"occurrences.")
            for x in st:
                print(x,end=" ")
        else:
            print(p,"--> 0 occurrences.")
        print("\n")