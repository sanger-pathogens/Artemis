#!/usr/local/bin/perl5.6.1

use strict;

BEGIN 
{

  my $configFile = "/nfs/pathdb/dev/go-cgi/config.pl";

  if (-f $configFile) 
  {
    require $configFile;
  }

}

use Getopt::Long;

my $associations;
my $definition;

&GetOptions("assoc" => \$associations,
            "def"   => \$definition);

use GO::AppHandle;

if(@ARGV == 0) 
{
 print "USAGE: go -def -assoc acc\n"; exit 0;
}

my $acc_num = $ARGV[0];

my $dbname = $ENV{GO_DBNAME};
my $dbport = $ENV{GO_DBPORT};
my $dbhost = $ENV{GO_DBHOST};
my $dbuser = $ENV{GO_DBUSER};

if (not $dbname =~ /go\d/) 
{
  print STDERR "GO database name uncorrect, must be go1 or go2 !!";
  exit 1;
}

my $apph = GO::AppHandle->connect (
                                   -dbname => $dbname,
                                   -dbport => $dbport,
                                   -dbhost => $dbhost,
                                   -dbuser => $dbuser,
                                  ) 
or die "can't connect to GO database, $dbname!!!\n";

#my @accs = qw(O00221);
my @accs;
push(@accs, $acc_num);

my @pqlist = map { {acc=>$_} } @accs;
my $term_l = $apph->get_terms({products=>[@pqlist]});

foreach my $term (@$term_l) 
{
  printf "%s; %s; ", $term->acc, $term->name;
  if($associations)
  {
    foreach my $assoc (@{$term->selected_association_list || []}) 
    {
      my $gp = $assoc->gene_product;
      my $ev_l = $assoc->evidence_list || [];
      printf "%s; %s; %s; evidence=%s; %s;", 
                        $gp->full_name, $gp->acc, $gp->symbol,
                        join('; ', map {$_->code} @$ev_l), $gp->as_str;

      foreach my $syn (@{$gp->synonym_list || []})
      {
        print " Synonym: $syn;";
      }
    }
  }

  if($definition)
  {
    printf " %s;", $term->definition;
  }

# printf "Synonyms:%s\n", join(", ", @{$term->synonym_list});

  print "\n";

}

