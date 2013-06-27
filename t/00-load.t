#!perl -T
use 5.006;
use strict;
use warnings FATAL => 'all';
use Test::More;

plan tests => 1;

BEGIN {
    use_ok( 'CRISPR_db' ) || print "Bail out!\n";
}

diag( "Testing CRISPR_db $CRISPR_db::VERSION, Perl $], $^X" );
