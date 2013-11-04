
from Bio import EUtils
from Bio.EUtils import HistoryClient


client = HistoryClient.HistoryClient()
records = client.search('Q8R5B6', db='protein')
print records.dbids

#result = client.post(EUtils.DBIds("protein", "4579714"))
#related = result.neighbor_links("protein")
#related_dbids = related.linksetdbs["protein_protein"].dbids
#proteins = client.post(related_dbids)
#len(proteins)
#infile = proteins.efetch(retmode = "text", rettype = "fasta")
#fasta = infile.read()
#print fasta[:788]
