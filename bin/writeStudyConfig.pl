#!/usr/bin/perl                                                                                                                                                                                            
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my ($help, $name, $file, $outputFile, $sourceIdType, @inputProtocolAppNodes, $protocol, @protocolParams, $profileSetName);

&GetOptions('help|h!' => \$help,
            'name=s' => \$name,
            'file=s' => \$file,
            'outputFile=s' =>\$outputFile,
            'sourceIdType=s' => \$sourceIdType,
            'inputProtocolAppNodes=s' => \@inputProtocolAppNodes,
            'protocol=s' => \$protocol,
            'profileSetName=s' => \$profileSetName,
            'protocolParams=s' => \@protocolParams # string of form param|paramValue                                                                                                                       
            );

&usage() if($help);

unless($file && $outputFile && $name && $protocol) {
    &usage();
}

sub usage {
    print STDERR "writeStudyConfig --name=s --file=s --outputFile <FILE> --sourceIdType=s --protocol=s [--inputProtocolAppNodes=list] [--protocolParams=list]\n";
    exit;
}

my @columns = ('Name', 'File Name', 'Source Id Type', 'Input ProtocolAppNodes', 'Protocol', 'ProtocolParams', 'ProfileSet');

open(OUT, "> $outputFile") or die "Cannot open $outputFile for writing: $!";

print OUT join("\t", @columns) . "\n";

print OUT "$name\t$file\t$sourceIdType\t".join(';', @inputProtocolAppNodes)."\t$protocol\t".join(';', @protocolParams)."\t" . $profileSetName . "\n";

close OUT;
