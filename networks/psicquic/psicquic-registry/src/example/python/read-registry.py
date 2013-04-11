import urllib2

urlStr = 'http://www.ebi.ac.uk/Tools/webservices/psicquic/registry/registry?action=ACTIVE&format=txt'
try:
    fileHandle = urllib2.urlopen(urlStr)
    content = fileHandle.read()
    fileHandle.close()
except IOError:
    print 'Cannot open URL' % urlStr
    content = ''

lines = content.splitlines()

for line in lines:
    fields = line.split('=')
    print fields[0] + ' ---> ' + fields[1]

