#!/acrm/usr/local/bin/perl

$tm_file = $ARGV[0];

open(TMFILE, "$tm_file") || die "ERROR: Can't open $tm_file FILE: $!\n";

@tm_file = <TMFILE>;

for($i = 0; $i <@tm_file; $i++)
{
    $structure = substr($tm_file[$i], 0, 4);
    $structurevoidfile = "$structure.voids";
    ($j,$j,$j,$j,$j,$j,$j,$size,@j) = stat($structurevoidfile);
    if($size == 0)
    {
        print STDERR "$structurevoidfile is empty\n";
    }
    else
    {
#        print "perl void_tm_printout.pl $structurevoidfile  $tm_file $structure\n";
        system("perl void_tm_printout.pl $structurevoidfile $tm_file $structure");
    }
}
