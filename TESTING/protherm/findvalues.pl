#!/usr/bin/perl
while(<>)
{
    chomp;
    @fields = split;

    $origpdb = $fields[@fields-3];
    $mutpdb = $fields[@fields-2];
    $tm = $fields[@fields-1];

    if($mutpdb eq "NULL")
    {
        $pdb = $origpdb;
    }
    else
    {
        $pdb = $mutpdb;
    }
    if($pdb ne "NULL")
    {
        print "$pdb $tm\n";
    }
}


