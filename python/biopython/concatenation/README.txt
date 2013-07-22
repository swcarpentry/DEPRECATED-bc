This is a (admittedly fake) sample sheet that describes an experiment where six samples (controls vs. biopsies) for three cartoon characters.

Since we had six samples and eight flowcell lanes, we wanted one and a third lanes per sample.
Our sequencing center director told us that we needed at least two samples per lane or the machine would become confused.

So we have a situation where two of our samples are split over two lanes and the rest of our samples are split over three lanes, and we desire to concatenate all of the samples form Donald Duck's biopsy into one (possibly large) data file.    After the fact, we need to check that the concatenated data are valid--that we got all the reads, that there are the number of reads that we exepect, and that the number of reads from read 1 match the number of reads from read 2.

This is the classic case for automation--you can write a script to concatenate the data by hand, or you can solve the more general problem of concatenating all the data form the same sample and changing the names of the files to things that are easier for you and your collaborators to understand.


