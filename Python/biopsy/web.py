
_transfac_major_version = 10
_transfac_minor_version = 2

def transfac_url( table_link ):
    return (
            'http://www.biobase-international.com/cgi-bin/biobase/transfac/%d.%d/bin/getTFProf.cgi?%s'
            %
            (
                    _transfac_major_version,
                    _transfac_minor_version,
                    table_link,
            )
    )

def ensembl_url( gene_id, species = 'Mus_musculus' ):
    return (
            'http://www.ensembl.org/%s/geneview?gene=%s'
            %
            (
                    species,
                    gene_id,
            )
    )
