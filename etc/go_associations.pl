#!/usr/local/bin/perl5.6.1

use strict;

BEGIN 
{

#  my $configFile = "/nfs/pathdb/dev/go-cgi/config.pl";
  my $configFile = "/nfs/pathdb/amigo/conf/config.pl";

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
my $dbauth = $ENV{GO_DBAUTH};

if (not $dbname =~ /^go/) {
  print STDERR "GO database name uncorrect, must start 'go...' !!";
  exit 1;
}

my $apph = GO::AppHandle->connect (
                                   -dbname => $dbname,
                                   -dbport => $dbport,
                                   -dbhost => $dbhost,
                                   -dbuser => $dbuser,
                                   -dbauth => $dbauth,
                                  ) 
    or die "can't connect to GO database, $dbname!!!\n";


#my @accs = qw(O00221);
my @accs;
push(@accs, $acc_num);

my @pqlist = map { {acc=>$_} } @accs;
my $term_l = $apph->get_terms({products=>[@pqlist]});

foreach my $term (@$term_l) 
{
  my $type = $term->term_type;
  if($type =~ m/_component/i)
  {
    printf "/GO_component=\"";
  }
  elsif($type =~ m/_function/i)
  {
    printf "/GO_function=\"";
  }
  elsif($type =~ m/_process/i)
  {
    printf "/GO_process=\"";
  }

  printf "%s (%s);", $term->acc, $term->name;
  if($associations)
  {
    foreach my $assoc (@{$term->selected_association_list || []}) 
    {
      my $gp = $assoc->gene_product;
      my $ev_l = $assoc->evidence_list || [];
      printf " %s; %s:%s", $assoc->evidence_as_str, $gp->speciesdb, $gp->acc;

      foreach my $syn (@{$gp->synonym_list || []})
      {
        print " ($syn)";
      }
      printf ";";

      printf " db_xref=";
      foreach my $ref (@{$term->dbxref_list || []})
      {
        printf "%s", $ref->as_str;
      }
      printf ";";
    }
  }

  if($definition)
  {
    printf " %s;", $term->definition;
  }

# printf "Synonyms:%s\n", join(", ", @{$term->synonym_list});

  print "\n";

}

