#!/usr/bin/perl

$_ = <>;
chomp;
($pdb_prev, $tm) = split;
$count = 1;
$tmtot = $tm;                                

while(<>)
{
    chomp;
    ($pdb, $tm) = split;
    if($pdb eq $pdb_prev)
    {
        $count++;
        $tmtot += $tm;
    }
    else
    {
        printf "$pdb_prev %f\n", $tmtot/$count;
        $count = 1;
        $pdb_prev = $pdb;
        $tmtot = $tm;
    }
}
printf "$pdb %.2f\n", $tmtot/$count;
