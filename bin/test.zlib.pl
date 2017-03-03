#!/usr/bin/perl -w 
use IO::Zlib;
    $fh = IO::Zlib->new("file.gz", "wb9");
    if (defined $fh) {
        print $fh "bar\n";
        $fh->close;
    }
