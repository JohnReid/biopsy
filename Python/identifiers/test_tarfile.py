import urllib2, gzip, csv, StringIO
url = 'ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz'
compressed_data = urllib2.urlopen(url).read()
compressed_stream = StringIO.StringIO(compressed_data)
gzipfile = gzip.GzipFile(fileobj=compressed_stream)
for l in csv.reader(gzipfile):
    print l
    break
