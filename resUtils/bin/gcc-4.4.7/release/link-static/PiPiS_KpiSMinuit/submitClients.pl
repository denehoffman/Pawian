#!/usr/bin/perl -w
my $queue="bsub -q pwa";
#my $queue="bsub -l nodes=4:ppn=1 -q pwa";
#my $queue="bsub -q pwaLarge";
print $queue."\n";
for ($i=1;$i<=8;$i++)
{
    my $client="Client$i";
    print $client."\n";
    my $submit="$queue -o batch$client.out \"../resReactionApp -c resReactionDefault.cfg --mode client >&! log.$client\"";
    print $submit."\n";
    system($submit); 
}
