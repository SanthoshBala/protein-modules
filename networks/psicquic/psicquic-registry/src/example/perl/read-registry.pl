#!/usr/bin/perl
use strict;
use LWP::Simple;

# Registry URL that lists all ACTIVE PSICQUIC services
# Documentation: http://code.google.com/p/psicquic/wiki/Registry
my $url = 'http://www.ebi.ac.uk/Tools/webservices/psicquic/registry/registry?action=ACTIVE&format=txt';

my $content = get $url;
die "Couldn't get $url" unless defined $content;

# Now list all services
my @lines = split(/\n/, $content);
for my $line (@lines){
  my @flds = split(/=/, $line);
  my ($serviceName, $serviceUrl) = @flds;
  print "$serviceName  --->  $serviceUrl\n";
}