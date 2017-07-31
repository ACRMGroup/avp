#!/usr/bin/perl
$pdbprep = "/acrm/data/pdb/pdb";
$pdbext = ".ent";

while(<>)
{
    chomp;
    @fields = split;
    $fields[0] =~ tr/A-Z/a-z/;
    $pdb = $pdbprep . $fields[0] . $pdbext;

    print "$_\n" if(HasWater($pdb));
}

sub HasWater
{
    my($pdb) = @_;
    my($haswater) = 0;

    if(open(PDB, $pdb))
    {
        while($line = <PDB>)
        {
            if($line =~ /^HETATM/)
            {
                if(($line = ~ /HOH/) ||
                   ($line = ~ /OH2/) ||
                   ($line = ~ /OHH/) ||
                   ($line = ~ /DOD/) ||
                   ($line = ~ /OD2/) ||
                   ($line = ~ /ODD/) ||
                   ($line = ~ /WAT/))
                {
                    $haswater = 1;
                    last;
                }
            }
        }
        close PDB;
    }
    return($haswater);
}
