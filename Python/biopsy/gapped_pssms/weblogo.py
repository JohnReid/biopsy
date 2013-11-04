#
# Copyright John Reid 2006
#

def dna_seqs_from_dist( dist, num_seqs = 100 ):
    """Generate a sequence of sequences that have a similar base distribution to
    dist

    dist: distribution over pssms
    num_seqs: number of sequences, the more the closer the match to the dist.
    """
    K = len( dist )
    seqs = [ ]
    for i in xrange( num_seqs ):
        seqs.append( '' )
        for j in xrange( K ):
            if i < num_seqs * dist[j,0]:
                seqs[i] += 'a'
            elif i < num_seqs * (dist[j,0] + dist[j,1]):
                seqs[i] += 'c'
            elif i < num_seqs * (dist[j,0] + dist[j,1] + dist[j,2]):
                seqs[i] += 'g'
            else:
                seqs[i] += 't'
    return seqs

def weblogo_data_from_dist( dist ):
    """Data for weblogo from pssm distribution"""
    import weblogolib, corebio.seq
    return weblogolib.LogoData.from_counts(
            corebio.seq.unambiguous_dna_alphabet,
            dist * 100
    )

def weblogo_from_dist( dist, filename = 'logo.eps' ):
    """Generate a weblogo from a pssm distribution"""
    import weblogolib
    data = weblogo_data_from_dist( dist )
    options = weblogolib.LogoOptions()
    options.size = weblogolib.LogoSize(
            stack_width = 5.4*12,
            stack_height = 5.4*12*5
    )
    options.color_scheme = weblogolib.std_color_schemes[ "classic" ]
    format = weblogolib.LogoFormat( data, options )
    weblogolib.eps_formatter(
            data,
            format,
            open( filename, 'w' )
    )

def format_weblogo_from_dist( dist, basename, ext, convert_args = '' ):
    """Generate a weblogo from a pssm distribution in format defined by ext"""
    import os.path
    d = os.path.dirname( basename )
    if '' != d and not os.path.exists( d ): os.makedirs( d )
    eps_file = basename + '.eps'
    converted_file = basename + '.' + ext
    if eps_file == converted_file:
        raise RuntimeError( 'Extension should not be same as eps' )
    weblogo_from_dist( dist, eps_file )
    convert_format( eps_file, converted_file, convert_args )
    os.remove( eps_file )


def convert_format( source, dest, convert_args = '' ):
    """Converts an image from one format to another

    Uses imagemagick convert program which must be installed

    Args:
            source: input image file
            dest: output image file
    """
    import os
    command = 'convert.exe %s "%s" "%s"' % ( convert_args, source, dest )
    status = os.system( command )
    if status:
        raise RuntimeError(
                'Could not convert "%s" to "%s".\nCommand: %s\nStatus: %d'
                % (
                        source,
                        dest,
                        command,
                        status
                )
        )

if '__main__' == __name__:
    import numpy, weblogolib
    dist = numpy.array(
            [
                    [ 0.25, 0.25, 0.25, 0.25 ],
                    [ 0.25, 0.25, 0.25, 0.25 ],
                    [ 0.25, 0.25, 0.00, 0.50 ],
                    [ 1.00, 0.00, 0.00, 0.00 ],
            ],
            numpy.float64
    )
    print '\n'.join( dna_seqs_from_dist( dist, num_seqs = 4 ) )
    weblogo_from_dist( dist, filename = 'test_logo.eps' )
    convert_format( 'test_logo.eps', 'test_logo.png' )
