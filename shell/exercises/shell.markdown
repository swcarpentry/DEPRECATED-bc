# Shell

Let's try out your new shell skills on some real data.

The file `1000gp.vcf` is a small sample (1%) of a very large text file
containing human genetics data. Specifically, it describes genetic variation in
three African individuals sequenced as part of the [1000 Genomes
Project](http://www.1000genomes.org). 

## Exercise Part 1 (setup)

* If you had forgotten where you downloaded the file, how would you locate the
path of all files with that name on the computer (using the shell)?  

  **Hint:**
  > `$ man find`

  **Answer:**
  > `$ find / -name "1000gp.vcf"`

* It's usually a good idea to use an empty directory as a workspace so that
  other files don't get in the way (or accidentally get overwritten or deleted).
  Create a new subdirectory directory named "sandbox", move our data file there,
  and make the directory your current working directory (sandbox should be the
  last part of the path given when you type `pwd`. 


  **Answer:**
  > ```bash
  > $ mkdir sandbox
  > $ mv /home/orion/Downloads/1000gp.vcf sandbox
  > $ cd sandbox
  > ```


## Exercise Part 2 (analysis)

* The data file you downloaded is a line-based text file. The "vcf" extension
  lets us know that it's in a specific text format, namely "Variant Call
  Format". The file starts with a bunch of comment lines (they start with "#" or
  "##"), and then a large number of data lines. The human genome can be thought
  of as an encyclopedia, where each chromosome is a volume. Each volume is just
  a long string of characters, but rather than the english alphabet, the genome
  uses just the characters "A", "C", "G", and "T". This VCF file lists the
  differences between the three African individuals and a standard "individual"
  called the reference (actually based upon a few different people). Each line
  in the file corresponds to a difference. The line tells us the position of the
  difference (chromosome and position), the genetic sequence in the reference,
  and the corresponding sequence in each of the three Africans. Research is
  ongoing to understand the full effects of these genetic differences; some
  cause diseases such as Tay-Sachs and Hemophilia, while others determine your
  blood type and eye color.

  Before we start processing the file, let's get a high-level view of the file
  that we're about to work with.

  What's the file size (in kilo-bytes), and how many lines are in the file?


  **Hint:**
  > There's an option to `ls` that will print the file sizes in a more
  > human-friendly format.


  **A hint about the number of lines:**
  > `$ man wc`


  **Answer:**
  > We should get a file size around 3.6 MB with:
  > `$ ls -lh 1000gp.vcf`
  > Alternatively, the command `du` can be used to achieve a similar result:
  > `$ du -h 1000gp.vcf`
  > 
  > We find there are 45034 lines with:
  > `$ wc -l 1000gp.vcf`


* Because this file is so large, you're going to almost always want to pipe
("|") the result of any command to `less` (a simple text viewer, type 'q' to
exit) or `head` (to print the first 10 lines) so that you don't accidentally
print 45,000 lines to the screen.

  Let's start by printing the first 5 lines to see what it looks like.  

  **Answer:**
  > `$ head -5 1000gp.vcf`

* That isn't very interesting; it's just a bunch of the comments at the
beginning of the file (they all start with "#")! Print the first 20 lines to see
more of the file.

  **Answer:**
  > `$ head -20 1000gp.vcf`


* Okay, so now we can see the basic structure of the file. A few comment lines
  that start with "#" or "##" and then a bunch of lines of data that contain all
  the data and are pretty hard to understand. Each line of data contains the
  same number of fields, and all fields are separated with TABs. These fields
  are:

  1. the chromosome (which volume the difference is in)
  2. the position (which character in the volume the difference starts at)
  3. the ID of the difference
  4. the sequence in the reference human(s)

  The rest of the columns tell us, in a rather complex way, a bunch of
  additional information about that position, including: the predicted sequence
  for each of the three Africans and how confident the scientists are that these
  sequences are correct.

  To start analyzing the actual data, we have to remove the header. How can we
  print the first 10 non-header lines (those that _don't_ start with a "#")?

  **Hint:**
    $ man grep

  **Hint:**
  > You can use a pipe ("|") to connect the output of `grep` to the input of
  > `head`.

  **Hint:**
  In `grep` regular expressions, the carat '^' character matches the start of a
  line and the dollar sign '$' matches the end of a line. Thus, the following
  will print all non-blank lines from `file`:
    $ grep -v "^$" file

  **Our answer:**
  >     $ grep -v "^#" 1000gp.vcf | head
  > 
  > Why are neither of these correct?
  >     $ grep -v "#" 1000gp.vcf | head
  >     $ grep -v "^##" 1000gp.vcf | head

* How many lines of data are in the file (rather than counting the number of
  header lines and subtracting, try just counting the number of data lines)?

  **Hint:**
  > Instead of piping to `head`, try piping to `wc`.

  **Our Answer:**
  >     $ grep -v "^#" 1000gp.vcf | wc -l
  >
  > should print `45024`

* Where these differences are located can be important. If all the differences
  between two encyclopedias were in just the first volume, that would be
  interesting. The first field of each data line is the name of the chromosome
  that the difference occurs on (which volume we're on). Print the first 10
  chromosomes, one per line.

  **Hint:**
  > You can extract a column from a tab-delimited text file using the `cut`
  > command.

  **Hint:**
  > Use `grep` to print only non-comment lines, and `cut` to extract the
  > chromosome column.

  **Our Answer:**
  >     $ grep -v "^#" 1000gp.vcf | cut -f 1 | head

* As you should have observed, the first 10 lines are on numbered chromosomes.
  Every normal cell in your body has 23 pairs of chromosomes, 22 pairs of
  "autosomal" chromosomes (these are numbered 1-22) and a pair of sex
  chromosomes (two Xs if you're female, an X and a Y if you're male). If you've
  heard of the genetics company [23andMe](https://www.23andme.com), the 23
  refers to these 23 pairs of chromosomes. 

  Let's look at which chromosomes these variations are on. Print a list of the
  chromosomes that are in the file (each chromosome name should only be printed
  once, so you should only print 23 lines).

  **Hint:**
  > You need to remove all the duplicate lines from your previous answer.

  **Hint:**
  > `sort` has an option that should make this easier.

  **Our Answer:**
  >     $ grep -v "^#" 1000gp.vcf | cut -f 1 | sort -u


* Rather than using `sort` to print unique results, a common pipeline is to
  first sort and then pipe to another UNIX command, `uniq`. The `uniq` command
  takes _sorted_ input and prints only unique lines, but it provides more
  flexibility than just using `sort` by itself. Keep in mind, if the input isn't
  sorted, `uniq` won't work properly.

  Using `sort` and `uniq`, print the number of times each chromosome occurs in
  the file.

  **Hint:**
  >     $ man uniq

  **Hint:**
  > Instead of using `sort` to remove duplicates, just use it to sort and pipe
  > the result to `uniq`.

  **Our Answer:**
  >     $ grep -v "^#" 1000gp.vcf | cut -f 1 | sort | uniq -c


* Add to your previous solution to list the chromosomes from most frequently
  observed to least frequently observed.

  **Hint:**
  >     $ man sort

  **Hint:**
  > Make sure you're sorting in descending order. By default, `sort` sorts in
  > ascending order.

  **Our Answer:**
  >     $ grep -v "^#" 1000gp.vcf | cut -f 1 | sort | uniq -c | sort -n -r
  >
  > should output the following:
  >
  >     3721 2
  >     3387 1
  >     3224 4
  >     3219 3
  >     2894 5
  >     2860 6
  >     2527 8
  >     2525 7
  >     2203 10
  >     2166 11
  >     2032 12
  >     1865 9
  >     1656 13
  >     1409 14
  >     1362 16
  >     1304 X
  >     1275 18
  >     1265 15
  >     1097 17
  >      993 20
  >      814 19
  >      661 21
  >      565 22
    
* The autosomal chromosomes (1-22) are named according to their size. The
  largest of them is chromosome 1, while the smallest is chromosome 22. Does it
  look like differences occur relatively randomly across the genome, or are some
  chromosomes more different than you'd expect at random (very roughly taking
  their sizes into account)?

  It's worth noting that the chromosomes were numbered by the sizes of the
  actual molecules, not how much of them had been sequenced.

  Wikipedia has a nice table of chromosome sizes and how much of each has been
  sequenced (and you can sort it):
  http://en.wikipedia.org/wiki/Human_chromosome#Human_chromosomes

  Notice anything?

  **Our Hypothesis:**
  > Since variation can only be found in the known sequence, the order you
  > printed corresponds closely to ordering by the number of bases sequenced
  > (rather than the total number of bases). 
  > 
  > Given this, it seems like differences occur relatively randomly across the
  > genome. We see more differences on longer chromosomes, fewer on shorter,
  > without any striking outliers.

* This is great, but biologists might also like to see the chromosomes ordered
  by their number (not dictionary order), since different chromosomes have
  different attributes and this ordering allows them to find a specific
  chromosome more easily.

  **Hint:**
  > A lot of the power of `sort` comes from the fact that you can specify which
  > fields to sort on, and the order in which to sort them. In this case you
  > only need to sort on one field.

  **Answer:**
  >     $ grep -v "^#" 1000gp.vcf | cut -f 1 | sort | uniq -c | sort -k 2n


## Exercise Part 3 (scripts and svn)
* Wonderful! Now we have a (long) command for printing chromosome statistics
  from our `1000gp.vcf` file. Using `nano`, create a new file, "chrom_stats.sh",
  with just your answer to the previous question in it.

  **Answer:**
  >  Type the following to open a new file:
  >     $ nano chrom_stats.sh
  >  Type in the command. Type ^o to save and ^x (where ^ means the control key).

* Just to be illustrate the flexibility of the shell, try creating the same file
  directly from the shell (without a text editor). Once you do, you can use
  `cat` to make sure the contents of the file are exactly what you expect.

  **Hint:**
  > You can use `echo` to print something and `>` to redirect to a file.

  **Hint:**
  > Since our long command has double-quotes in it, you either need to use
  > single-quotes or escape these with back-slashes.

  **Answer:**
    $ echo 'grep -v "^#" 1000gp.vcf | cut -f 1 | sort | uniq -c | sort -k 2n' > chrom_stats.sh
    $ cat chrom_stats.sh

* Now, execute your new script to print the chromosome statistics.

  **Hint:**
  > You may have to change the permissions to allow you to execute it.

  **Hint:**
  > It's good form to only make permissions as permissive as necessary. So,
  > rather than allow everyone to execute the file, it is better to just allow
  > you to execute it.

  **Answer:**
    $ chmod u+x chrom_stats.sh
    $ ./chrom_stats.sh

  > Note that it is `u+x` instead of just `+x` or `a+x`. This only adds the
  > ability for the owner to execute it, whereas the other two options would
  > allow anyone to execute it.

* We'd like to be able to use this script in the future with arbitrary VCF
  files, instead of just our `1000gp.vcf` file. Edit the script so that it takes
  VCF-formatted text input on stdin and prints out chromosome statistics on
  stdout. This is simpler than you might think.

  **Hint:**
  > If `grep` isn't given an input file, it will read from stdin.

  **Answer:**
  > Change
  > `grep -v "^#" 1000gp.vcf | ...`
  > to
  > `grep -v "^#" | ...`
  > 
  > Since this is in a file instead of the shell prompt, we aren't showing the
  > "$" at the beginning of the line.

* Now that we have a script that reads from stdin and prints to stdout, how do
  we run it on the `1000gp.vcf` file to get the same output as before?

  **Hint:**
  > The `cat` command is used to print files to stdout.

  **Hint:**
  > You can pipe the output of `cat` directly into our script.

  **Hint:**
  > Just like before, in order to tell the shell that the `chrom_stats.sh` file
  > we want to execute is the one in our current directory, we need to use
  > `./chrom_stats.sh`.

  **Answer:**
  `$ cat 1000gp.vcf | ./chrom_stats.sh`

* Finally, add a copy of this file to your folder in the class SVN repository.
  1. `cp chrom_stats.sh /path/to/repo/participants/user/`
  2. Add the file to subversion version control
  3. Commit your changes

**Fin.** 
Comments, questions, and suggestions are encouraged and appreciated.
Thanks to Tommy Guy, Jon Pipitone, Greg Wilson, and Elango Cheran for their help
with these exercises.

