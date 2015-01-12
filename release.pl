#!/usr/bin/perl



#dynamically generate a C program (with a baked hash or tag)
#and compile it for release
#see here http://stackoverflow.com/questions/3442874/in-git-how-can-i-write-the-current-commit-hash-to-a-file-in-the-same-commit
#include <stdio.h>
#int main(int numargs,char** args)
#	{
#	if(numargs==2)
#		{
#		if(strcmp(args[1],"desc")==0)
#			{
#			printf("desc\n");
#			}
#		else if(strcmp(args[1],"tag")==0)
#			{
#			printf("tag\n");
#			}
#		else
#			{
#			printf("tag\n");
#			}	
#		}
#	else
#		{
#		printf("tag\n");
#		}
#	return 0;
#	}
my $whichGit=`which git`;
$whichGit=trim($whichGit);
my $whichGcc=`which gcc`;
$whichGcc=trim($whichGcc);
if(length($whichGcc)==0)
	{
	die "Error, fail to find gcc!";
	}
if(length($whichGit)==0)
	{
	die "Error, fail to find git!";
	}

my $commitHash=`$whichGit rev-parse HEAD`;
$commitHash=trim($commitHash);
my $commit=$commitHash;
my $commitDesc=`$whichGit describe`;
$commitDesc=trim($commitDesc);
print "The commit is $commit\n";
my $cCode="";
$cCode.="#include <stdio.h>\n";
$cCode.="int main(int numargs,char** args)\n";
$cCode.="{\n";
$cCode.="if(numargs==2) {\n";
$cCode.="if(strcmp(\"desc\",args[1])==0) printf(\"".$commitDesc."\\n\");\n";
$cCode.="else if(strcmp(\"hash\",args[1])==0) printf(\"".$commitHash."\\n\");\n";
$cCode.="else printf(\"".$commitDesc."\\n\");\n";
$cCode.="}\n";
$cCode.="else { printf(\"".$commitDesc."\\n\"); }\n";
$cCode.="return 0;\n";
$cCode.="}\n";
my $tmpPath=$commit.".c";
open(WRITER,">$tmpPath");
print WRITER $cCode;
close(WRITER);
print "Wrote to $tmpPath\n";
my $gccCmd="$whichGcc $tmpPath -o release_info";
`$gccCmd`;
unlink($tmpPath);


#http://perlmaven.com/trim
sub  trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };
