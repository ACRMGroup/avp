#!/usr/local/bin/perl 


$step   = 15;
$rotate = '/home/andrew/bin/rotate';
$avp    = '/home/andrew/voids/avp';

$pdb = $ARGV[0];

close STDIN;
$|=1;

$Sx      = 0;
$SxSq    = 0;
$NValues = 0;
$mean    = 0;
$sd      = 0;
$minval  = 999999999.0;
$maxval  = 0.0;

for($x=0; $x<360; $x+=$step)
{
    for($y=0; $y<360; $y+=$step)
    {
        for($z=0; $z<360; $z+=$step)
        {
            $results = `$rotate -x $x -y $y -z $z $pdb | $avp -q -r`;
            @lines = split('\n', $results);
            ($junk,$vol) = split(': ', $lines[$#lines]);
            print "$vol\n";
            CalcExtSD($vol, 0, \$Sx, \$SxSq, \$NValues, \$mean, \$sd);
            if($vol > $maxval)
            {
                $maxval = $vol;
            }
            if($vol < $minval)
            {
                $minval = $vol;
            }
        }
    }
}

CalcExtSD(0, 1, \$Sx, \$SxSq, \$NValues, \$mean, \$sd);

print "Minimum: $minval Mean: $mean Maximum: $maxval StdDev: $sd\n";

#*************************************************************************
#   void CalcExtSD(REAL val, int action, REAL *Sx, REAL *SxSq, 
#                  int *NValues, REAL *mean, REAL *SD)
#   ----------------------------------------------------------
#   Calculate the mean and standard deviation from a set of numbers. 
#   The routine is called with each value to be sampled and the action 
#   required is specified:
#
#   Input:   val     int       The value to be sampled
#            action  short     0: Sample the value
#   1: Calculate & return mean and SD
#   2: Clear the sample lists
#   Output:  mean    *REAL     The returned mean
#            SD      *REAL     The returned standard deviation
#   I/O:     Sx      *REAL     Sum of values
#            SxSq    *REAL     Sum of values squared
#            NValues *int      Number of values
#
#   The output values are only set when action==1
#
#   This is the same as CalcSD except that the Sx, SxSq and NValues
#   variables are kept outside the function instead of being static
#   within the function
#
#   13.10.93 Original based on CalcSD   By: ACRM
#   22.06.94 Fixed for only one value supplied
#   22.11.99 Translated from C to Perl
#
sub CalcExtSD
{
    my($val, $action, $Sx, $SxSq, $NValues, $mean, $SD) = @_;

    if($action==0)
    {
        ($$NValues)++;
        $$SxSq += ($val * $val);
        $$Sx   += $val;
    }
    elsif($action==1)
    {
        $$mean = $$SD = 0.0;
        $$mean = ($$Sx) / ($$NValues)  if($$NValues > 0);
        $$SD   = sqrt(($$SxSq - (($$Sx) * ($$Sx)) / ($$NValues)) / ($$NValues - 1))  if($$NValues > 1);
    }
    else
    {
        $$SxSq    = 0.0;
        $$Sx      = 0.0;
        $$NValues = 0;
    }
}

sub PrintResults
{
    CalcExtSD($LoopLength, 0, \$Sx, \$SxSq, 
              \$NValues, \$mean, \$SD);

    CalcExtSD($sa, 1, \$Sx, \$SxSq, \$NValues, \$mean, \$SD);
    printf "   Loops:   mean=%.3f SD=%.3f\n", $mean,  $SD;
}

