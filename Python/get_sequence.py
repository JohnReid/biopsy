
import biopsy

try: remome
except NameError:
    remome = biopsy.Remome.load( 'c:/data/ReMos/100/100.filtered' )


# get the remos position and sequence
as = remome.get_aligned_sequences()[0]
cs = as.centre_sequence
print 'Looking at: %s' % str( cs )
region = as.get_sequence_info(cs).region
print region
remos = remome.get_remos_for( as )
remo = remos[0]
seqs = remo.get_sequences( cs )
assert( len( seqs ) == 1 )
location = seqs[0].location
print str( location )
bases = seqs[0].unmasked_sequence
print bases


# where does the gene start in the chromosome?
gene = str( cs.gene_id )
genome = str( cs.db_id )
gene_locations = biopsy.gene_locations( genome )
gene_location = gene_locations[gene]


# get the right part of the chromosome
chr = biopsy.get_chromosome_file( genome, gene_location.chromosome )
genomic_seq = biopsy.get_chromosome_sequence(
        chr,
        gene_location.start - location.end,
        location.end - location.start + 1
)
print genomic_seq
