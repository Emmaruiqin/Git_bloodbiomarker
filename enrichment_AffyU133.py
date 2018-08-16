from mybio import analysis
from glob import glob
import os
from collections import defaultdict

def verify(gene):
    gene = gene.strip()
    try:
        int(gene)
        return gene
    except ValueError:
        return False

for file in glob('filter_probe.txt'):
    genelist = []
    name = os.path.splitext(file)[0]
    print(name)
    head, *data = open(file)
    head = head.strip().split('\t')
    num1 = head.index('Entrez Gene') if 'Entrez Gene' in head else head.index('GeneID')
    FC_name = [i for i in head if 'FC' in i and 'logFC' not in i]
    # FC_name = [i for i in head if 'FC' in i and 'logFC' not in i and 'rma-gc-scale' in i]
    num2 = head.index(FC_name[0])
    alltar = {'all': set(), 'up': set(), 'down': set()}
    generegulation = defaultdict(set)
    for line in data:
        ll = line.strip().split('\t')
        genes = ll[num1].split('///')

        if float(ll[num2]) > 1:
            regu = 'up'
        elif float(ll[num2]) < 1:
            regu = 'down'

        for gene in genes:
            gene = gene.strip()
            if verify(gene):
                alltar['all'].add(gene)
                alltar[regu].add(gene)
                generegulation[gene].add(regu)

    generegulation1 = {}
    for key, value in generegulation.items():
        if len(value) > 1:
            generegulation1[key] = "other"
        else:
            generegulation1[key] = value.pop()
    for regulation, genes in alltar.items():
        if 'all' in regulation:
            analysis(genes, name + '-' + regulation, generegulation1, org='hsa')
