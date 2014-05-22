#!/usr/bin/perl
use strict;
open FILE, "<", "./debug/posra_test.txt" or die $!;
my @image_paths = <FILE>;
close(FILE);
my $counter = 0;
my $total = 0;
foreach (@image_paths) {
      if (substr($_, 0, 1) ne "#") {
            print $_;
            $counter += test($_);
            print "\n";
            ++$total
      }
}
print "$counter / $total tests passed.";

sub test {
      my $res = 0;
      my @io = split(/\|/, $_[0]);
      $io[0] =~s/ /\\ /g;
      my $input = "./src/osra -f can -r 150 " . $io[0] . " |";
      open POSRA, $input or die $!;
      my @output = <POSRA>;
      close(POSRA);
      foreach (@output) {
            print "\t$_";
      }
      my $i = ((scalar @output) < 2) ? 0 : 1;
      my $fail = ((scalar @output) > 0) ? 0 : 1;
      chomp($output[0]);
      chomp($io[(scalar @io) - 1]);
      if (($output[0] eq $io[(scalar @io) - 1]) and !$fail) {
            print "\t[  \e[32mOK\e[0m  ]\n";
            ++$res;
      } else {
            print "\t[ \e[31mFAIL\e[0m ]\n";
      }
      return $res;
}
