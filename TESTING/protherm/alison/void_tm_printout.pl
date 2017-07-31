#!/acrm/usr/local/bin/perl

$pdbprep = "/acrm/data/pdb/pdb";
$pdbext  = ".ent";
$countpdb = "/home/bsm/martin/bin/countpdb";
$getresol = "/home/bsm/martin/cprogs/jiffies/getresol";

$void_file = $ARGV[0];
$tm_file = $ARGV[1];
$pdbstructure = $ARGV[2];

open(TMFILE, "$tm_file") || die "ERROR: Can't open $tm_file FILE: $!\n";
$printed = 0;

$pdbcode = "\L$pdbstructure";
$pdbfile = "$pdbprep$pdbcode$pdbext";

while(<TMFILE>)
{
    
    @fields = split;

    $pdbname = $fields[0];
    $pdbtemp = $fields[1];

    if($pdbname eq $pdbstructure)
    {
	$file_tm = "$pdbname\t$pdbtemp\t";
        $printed++;
        last;
    }
}

open(VOIDFILE, "$void_file") || die "ERROR: Can't open $void_file FILE: $!\n";


while(<VOIDFILE>)
{
    if(/Total\s+void\s+volume:\s+([\d\.]+)/)
    {
	$total_void_volume = $1;
        $printed++;
    }
    
    if(/Largest\s+void\s+volume:\s+([\d\.]+)/)
    {
	$largest_void_volume = $1;
        $printed++;
    }
}

if($printed == 3)
{
    $rescount = GetRescount($pdbfile);
    $resol = GetResol($pdbfile);

    print $file_tm;
    print "$rescount\t";
    print "$resol\t";
    print "$total_void_volume\t";
    print "$largest_void_volume\n";
}


sub GetRescount
{
    my($pdbfile) = @_;
    my($result, @fields);

    $result = `$countpdb $pdbfile`;
    chomp $result;
    @fields = split(/\s+/, $result);
    return($fields[3]);
}

sub GetResol
{
    my($pdbfile) = @_;
    my($result, @fields);
    $result = `$getresol $pdbfile`;
    chomp $result;
    @fields = split(/[\s\,\/]+/, $result);
    $fields[1] =~ s/A//;
    return($fields[1]);
}
