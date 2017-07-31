#!/usr/local/bin/perl 

$pdbprep = "/acrm/data/pdb/pdb";
$pdbext = ".ent";
$domprep = "/acrm/data/dompdb/";
$domext = "";
$file = $ARGV[0];
$output_dir = $ARGV[1];

open(FILE, "$file") || die "Oops ... can't open file: $!\n";

while(<FILE>)
{
    @fields = split;
    $var = $fields[0];
    $pdbcode = substr($var, 0,4);
    $chain = substr($var, 4,1);
    #$domain = substr($var, 5,1);
                          
    $fileloc = "$pdbprep$pdbcode$pdbext";
    $fileloc = "\L$fileloc";
    print "$fileloc\n";
    CalcVoid($fileloc, $pdbcode);

}

#####################################################

sub CalcVoid
{
    my($location, $nme) = @_;
    #print "/acrm/home/andrew/avp/avp -p 0.5 -R $location > $output_dir$nme.voids\n"; 
    system("/acrm/home/andrew/avp/avp -p 0.5 -R $location > $output_dir$nme.voids");
}

#####################################################

