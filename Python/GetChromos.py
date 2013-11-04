#
# Copyright John Reid 2006
#
import urllib
import os

# Details of parts of the ensembl url
urlprefix="ftp://ftp.ensembl.org/pub/"
data_dir="/data/fasta/dna/"
file_part=".may.dna.chromosome."
ext=".fa.gz"

# Where the files will be stored
local_dir="c:\\data\\ensembl\\chromosones\\zipped\\"

# Get one file
def get_chromo(version, filename):
    filepath=local_dir+filename
    if not os.path.isfile(filepath):
        print version, " ", filename
        urllib.urlretrieve(
            urlprefix+version+data_dir+filename,
            filepath)


# Get all files for one version of one species
def get_chromos(version, name, number, extras):
    for e in extras:
        get_chromo(version, name+file_part+e+ext)
    for n in range(1,number+1):
        get_chromo(version, name+file_part+str(n)+ext)



#current_chicken Gallus_gallus.WASHUC1
get_chromos(
    "current_chicken",
    "Gallus_gallus.WASHUC1",
    24,
    ["Un", "W", "Z", "26", "27", "28", "32"])

#current_human Homo_sapiens.NCBI35
get_chromos(
    "current_human",
    "Homo_sapiens.NCBI35",
    22,
    ["X", "Y"])

#current_mouse Mus_musculus.NCBIM33
get_chromos(
    "current_mouse",
    "Mus_musculus.NCBIM33",
    19,
    ["X", "Y"])

#current_chimp Pan_troglodytes.CHIMP1
get_chromos(
    "current_chimp",
    "Pan_troglodytes.CHIMP1",
    23,
    ["X", "Y"])

#current_rat Rattus_norvegicus.RGSC3.4
get_chromos(
    "current_rat",
    "Rattus_norvegicus.RGSC3.4",
    20,
    ["X"])

#current_tetraodon Tetraodon_nigroviridis.TETRAODON7
get_chromos(
    "current_tetraodon",
    "Tetraodon_nigroviridis.TETRAODON7",
    21,
    [])

#current_dog Canis_familiaris.BROADD1
get_chromos(
    "current_dog",
    "Canis_familiaris.BROADD1",
    38,
    ["X", "Un"])

#current_zebrafish Danio_rerio.ZFISH4
get_chromos(
    "current_zebrafish",
    "Danio_rerio.ZFISH4",
    25,
    [])

#current_fugu Fugu_rubripes.FUGU2
urllib.urlretrieve(
    "ftp://ftp.ensembl.org/pub/current_fugu/data/fasta/dna/Fugu_rubripes.FUGU2.may.dna.scaffold.fa.gz",
    local_dir+"Fugu_rubripes.FUGU2.may.dna.scaffold.fa.gz")
