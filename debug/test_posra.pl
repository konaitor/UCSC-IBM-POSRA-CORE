#!/usr/bin/perl
use strict;
open FILE, "<", "./debug/posra_tests.txt" or die $!;
my @image_paths = <FILE>;
close(FILE);
my $counter = 0;
foreach (@image_paths) {
      if (substr($_, 0, 1) ne "#") {
            print $_;
            $counter += test($_);
            print "\n";
      }
}
my $total = scalar(@image_paths);
print "$counter / $total tests passed.";

sub test {
      my $res = 0;
      my @io = split(" ", $_[0]);
      $io[0] =~s/ /\\ /g;
      my $input = "./src/osra " . "\"$io[0]\"" . " |";
      open POSRA, $input or die $!;
      my @output = <POSRA>;
      close(POSRA);
      foreach (@output) {
            print "\t$_";
      }
      my $i = ((scalar @output) < 2) ? 0 : 1;
      chomp($output[$i]);
      chomp($io[1]);
      if ($output[$i] eq $io[1]) {
            print "\t[  \e[32mOK\e[0m  ]\n";
            ++$res;
      } else {
            print "\t[ \e[31mFAIL\e[0m ]\n";
      }
      return $res;
}
