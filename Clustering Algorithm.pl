#!usr/bin/perl


#AIM: To Cluster Images using correlation coefficients between the images.
$location="/home/yogesh/"
opendir(DIR,"$location.corrnew") or die "the directory could not be found $!";   #directory containing the files containing the correlation coefficients between images
open(FIL,">Clusters-version0.txt") or die "cant open the file $!";   #Output file to store the clustered image names
my @files=readdir(DIR);
my $count=0;
my $total=1;
my $sum=0;
my %Head_hash;  #contains the member to head associations. used to avoid cases of creating a new cluster on a member to add a new element.
my $maxcor=0;

my %checkHash; #way to avoid duplicates for clustering elements.
my %testHash;   #Contains all the unique elements
my %Cluster; #contain all the clustered image names and its members.
my %hash;   #the hash which stores the img_name;img_name=>correlation_coefficient

foreach $file (@files)  #Read contents from all correlation files.
{

    if($file =~/corrstage4-0/)    #we perform clustering for stage 4 for the Drosophila images.
    {
        open(FILE,"/home/yogesh/Desktop/corrnew/$file") or die "cant open the file $!";
        print "$file"."\npress...";
        #$m=<STDIN>;
        my @file_contents=<FILE>;
        foreach $line (@file_contents)  #The format present in the file is: img_name1 gene1 stage Img_name2 gene2 corr_coef
        {
            my @words=split(/\t/,$line);
            #print @words;
            #print $words[0].";".$words[1].";".$words[2].";".$words[3].";".$words[4].";".$words[5];
            if($words[3] =~/insitu/)
            {
                $checkHash{$words[3]}=1;
                $checkHash{$words[0]}=1;
                $testHash{$words[3]}=1;
                $testHash{$words[0]}=1;
                $hash{$words[0]}[0]{$words[3]}[0]=$words[5];
                $hash{$words[0]}[0]{$words[3]}[1]=$words[4];
                $hash{$words[0]}[1]=$words[1];
         
                $total++;
                $sum=$sum+$words[5];
                
            }    
        }    
    }

    #$m=<STDIN>;
   # print FIL "\nHASHED HEADERS";
    foreach $img (keys %hash)
    {
 #       print FIL $img."\n";
       # print $img."\n";       
    }
    
    #Determining the threshold for Clustering. The correlation values above the threshold implies that the images belong to same group.
    my $avg=$sum/$total;
    
  #  print  "\n\nTHE count is $count out of $total\nMAX value $maxcor \n max images are $maximg1 and $maximg2 \n average is $avg\n";
}   
   
   
   #####################
   #CLUSTERING ALGORITHM.
   #####################
   
   # UNION_FIND structure $Head_hash is a member to Head Association structure. Every HEAD has its head value to be -1. Every member has its head value to be its Parent.
   # Depending upon the strength of the correlation coefficients, We fill in this structure with the image names and cluster them by tracing the Head for every member.
   # structure:  $Head_hash{$member}=$head;
   
   # Initialize Union_Find Structure.
   foreach $key (keys %testHash)
   {
    $Head_hash{$key}[0]="-0.1";
    $Head_hash{$key}[1]="-0.2";
   }
   
   $i=0;
   foreach $key1 (keys %testHash)
   {
        foreach $key2 (keys %testHash)
        {
	    print $i."\n";
	    $i++;
            if($key1!~$key2)
            {
                if($hash{$key1}[0]{$key2}[0]>0.5 || $hash{$key2}[0]{$key1}[0]>0.5) 
                {
                    my $value;
                    
                    if($hash{$key1}[0]{$key2}[0]>0.5)
                    {
                        $value=$hash{$key1}[0]{$key2}[0];
                    }
                    if($hash{$key2}[0]{$key1}[0]>0.5)
                    {
                        $value=$hash{$key2}[0]{$key1}[0];
                    }
                    
                    if($Head_hash{$key1}[0] eq "-0.1") #If key1 is new
                    {
			print "both $key1 and $key2 are new\n";
			
                        if($Head_hash{$key2}[0] eq "-0.1")  #both key1 and key2 are new =>make one the Head and other follows it
                        {
                            $Head_hash{$key1}[0]=$key2;
                            $Head_hash{$key1}[1]=$value;
                            $Head_hash{$key2}[0]="-1";                        
                        }
                        else    #key1 is new and key2 may be a Head or a Member => make key1 a member of key2
                        {
                           $Head_hash{$key1}[0]=$key2;
                           $Head_hash{$key1}[1]=$value;
                        }
                    }
                    else
                    {
                        if($Head_hash{$key2}[0] eq "-0.1")   #key2 is new and key1 may be a Head or a Member => make key2 a member of key1
                        {
                            $Head_hash{$key2}[0]=$key1;
                            $Head_hash{$key2}[1]=$value;
                            
                        }
                        else   #BOTH ARE NOT NEW.. THEY COULD BE MEMBERS OR HEADS
                        {
			    
                            if($Head_hash{$key1}[0]!~/-1/)  #key1 is a member
                            {
                                if($Head_hash{$key2}[0]!~/-1/)  #key2 is a member BOTH ARE MEMBERS
                                {
				    print "Both $key1 and $key2 are members \n";
                                    if($Head_hash{$key1}[0]>$value) #if key1's correlation coefficient to its Head is greater than that between key1 and key2 
                                    {
                                        if($Head_hash{$key2}[0]<$value) #if key2's correlation coefficient to its Head is lesser than that between key1 and key2 => key2 follows key1
                                        {
                                            $Head_hash{$key2}[0]=$key1;
                                            $Head_hash{$key2}[1]=$value;
                                        }
                                           #if key1's and key2's correlation coefficient to their Head are greater than that between key1 and key2 =>dont bother about the new association
                                    }
                                    else   #if key1's correlation coefficient to its Head is lesser than that between key1 and key2
                                    {
                                        if($Head_hash{$key2}[0]>$value)  #if key2's correlation coefficient to its Head is greater than that between key1 and key2 =>key1 follows key2
                                        {
                                            $Head_hash{$key1}[0]=$key2;
                                            $Head_hash{$key1}[1]=$value;
                                        }
                                        else  #if key1's and key2's correlation coefficient to their Head are lesser than that between key1 and key2 =>choose the one with larger coeff to its head and the other follows this one.
                                        {
                                            if($Head_hash{$key1}[1]> $Head_hash{$key2}[1]) #key2 follows key1
                                            {
                                                $Head_hash{$key2}[0]=$key1;
                                                $Head_hash{$key2}[1]=$value;
                                            }
                                            else
                                            {
                                                $Head_hash{$key1}[0]=$key2;
                                                $Head_hash{$key1}[1]=$value;
                                            }
                                        }
                                    }
                                }
                                else  #key1 is a member and key2 is a head . 
                                {
                                    if($Head_hash{$key1}[1]>$value)  #if member(key1) to its head correlation is greater than that of the key1-key2 correlation => key2 follows key1
                                    {
                                        $Head_hash{$key2}[0]=$key1;
                                        $Head_hash{$key2}[1]=$value;
                                    }
                                    else               #else member(key1) to its head correlation is lesser than that of the key1-key2 correlation => key1 follows key2. 
                                    {
                                        $Head_hash{$key1}[0]=$key2;
                                        $Head_hash{$key1}[1]=$value;
                                    }
                                }
                            }
                            else #key1 is a head
                            {
                                if($Head_hash{$key2}[0]!~/-1/) #key1 is a head and key2 is a member
                                {
                                    if($Head_hash{$key2}[1]>$value) #if member(key2) has greater correlation coefficient to its head than the key1-key2 correlation coefficient then head(key1) follows key2
                                    {
                                        $Head_hash{$key1}[0]=$key2;
                                        $Head_hash{$key1}[1]=$value;
                                    }
                                    else   #if member correlation to its head is lesser than the correlation between the key1-key2 then the member(key2) now follows the new key1
                                    {
                                        $Head_hash{$key2}[0]=$key1;
                                        $Head_hash{$key2}[1]=$value;
                                    }
                                }
                                else  #BOTH ARE HEADS Choose a head to follow the other.
                                {
				    print "Both $key1 and $key2 are heads \n";
                                    $Head_hash{$key1}[0]=$key2;
                                    $Head_hash{$key1}[1]=$value;
                                }
                            }
                        }
                    }
                }
				
            }
        }
    }
   
 #Ithe correlation is below the threshold.
 # if they are new make them as heads or dont bother
   foreach $key (keys %testHash)
    {
	if($Head_hash{$key}[0] eq "-0.1")
	{
	$Head_hash{$key}[0]="-1";
	$Head_hash{$key}[1]="-2";
	}
    }
   
   foreach $key (keys %Head_hash)
   {
    print "key: $key and value : $Head_hash{$key}[0]\n";
   }
   print "press...\n";
   $my=<STDIN>;
    #Enumerating the Clusters
    # In order to enunciate the clusters, we need to find heads of all the nodes. This is done by FIND operations.
    # We traverse the head for every member till the head becomes -1. The head becomes the member on the next iteration.
    # The traversal is recursive and we need to prevent the recursion{Takes Too long!}
    # we can create a reverse hash{maps Values to the Key}.
    # Here the Reverse_Hash adds the members to their Heads{which are obtained from the Head_hash}
    # It also adds the members of the members{if any(obtained from previous additions}} to the Head.
   my %Reverse_hash;
   foreach $key (keys %testHash)
   {
    if($Head_hash{$key}[0] ne "-1")
    {
	push @{$Reverse_hash{$Head_hash{$key}[0]}},$key;
    }
    if(exists $Reverse_hash{$key})
    {
	push @{$Reverse_hash{$Head_hash{$key}[0]}},@{$Reverse_hash{$key}};
	delete $Reverse_hash{$key};
    }
   }
   
   foreach $key (keys %Reverse_hash)
   {
     if($key ne "-1")
     {
	foreach $member (@{$Reverse_hash{$key}})
	{
	    if(exists $Reverse_hash{$member})
	    {
	       push @{$Reverse_hash{$Head_hash{$key}[0]}},@{$Reverse_hash{$key}};
	       delete $Reverse_hash{$key};    
	    }
	}
     }
   }
   
   foreach $key (keys %Reverse_hash)
   {
     if($key ne "-1")
     {
	foreach $member (@{$Reverse_hash{$key}})
	{
	    if(exists $Reverse_hash{$member})
	    {
	       push @{$Reverse_hash{$Head_hash{$key}[0]}},@{$Reverse_hash{$key}};
	       delete $Reverse_hash{$key};    
	    }
	}
     }
   }
   foreach $key (keys %Reverse_hash)
   {
     if($key ne "-1")
     {
	foreach $member (@{$Reverse_hash{$key}})
	{
	    if(exists $Reverse_hash{$member})
	    {
	       push @{$Reverse_hash{$Head_hash{$key}[0]}},@{$Reverse_hash{$key}};
	       delete $Reverse_hash{$key};    
	    }
	}
     }
   }
   foreach $key (keys %Reverse_hash)
   {
     print "}\n\n$key:\n{";
     foreach $member (@{$Reverse_hash{$key}})
     {
	print "$member,";
     }
   }
