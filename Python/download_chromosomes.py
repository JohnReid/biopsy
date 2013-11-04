
import os.path, sys, ftplib, re, chromosome, gzip

_ensembl_ftp = 'ftp.ensembl.org'
_top_dir = '/pub'
_chromosome_re = re.compile( "\.dna\.chromosome\..*.fa.gz$" )

def gettext(ftp, filename, outfile=None):
    # fetch a text file
    if outfile is None:
        outfile = sys.stdout
    # use a lambda to add newlines to the lines read from the server
    ftp.retrlines("RETR " + filename, lambda s, w=outfile.write: w(s+"\n"))

def getbinary(ftp, filename, outfile=None):
    # fetch a binary file
    if outfile is None:
        outfile = sys.stdout
    ftp.retrbinary("RETR " + filename, outfile.write)


def download_chromosomes(
        software_version,
        ncbi_build,
        ncbi_tag,
        organism_name = 'mouse',
        organism = 'mus_musculus',
):
    print 'Software version: %d' % software_version
    print 'NCBI build: %d' % ncbi_build
    print 'NCBI tag: %s' % ncbi_tag
    print 'Organism name: %s' % organism_name
    print 'Organism: %s' % organism
    print

    data_dir = '%s/%s_core_%d_%d%s' % (
            chromosome._base_data_dir,
            organism,
            software_version,
            ncbi_build,
            ncbi_tag
    )
    print 'Downloading data to directory: %s' % data_dir
    if not os.path.exists( data_dir ):
        os.makedirs( data_dir )

    ftp = ftplib.FTP( _ensembl_ftp )
    ftp.login()
    d = '%s/release-%d/%s-%d.%d%s/data/fasta/dna' % (
            _top_dir,
            software_version,
            organism_name,
            software_version,
            ncbi_build,
            ncbi_tag
    )
    #print d
    ftp.cwd( d )

    dir_listing = []
    ftp.dir( dir_listing.append )
    for l in dir_listing:
        fields = l.strip().split()
        if len( fields ) > 8 and _chromosome_re.search( fields[8] ):
            filename = fields[8]
            local_filename = '%s/%s' % (data_dir, filename)
            if not os.path.exists( local_filename ):
                print 'Downloading %s to %s' % ( filename, local_filename )
                getbinary(
                        ftp,
                        filename,
                        open(
                                local_filename,
                                'wb'
                        )
                )
            unzipped_filename = local_filename.strip( '.fa.gz' )
            if not os.path.exists( unzipped_filename ):
                print 'Unzipping: unzipped_filename'
                zipped_file = gzip.GzipFile( local_filename )
                zipped_file.readline() # strip header
                unzipped_file = open( unzipped_filename, 'w' )
                for l in zipped_file:
                    unzipped_file.write( l.strip() )

    return ftp


for sv, build, tag in [
        ( 28, 33, 'd' ),
        ( 31, 33, 'g' ),
        ( 32, 34, '' ),
]:
    ftp = download_chromosomes(
            sv,
            build,
            tag,
            'mouse',
            'mus_musculus'
    )
