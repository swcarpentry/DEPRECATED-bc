Notes on teaching git remote with toy (occasion-built) repositories.

Aron Ahmadia created a toy repository for 2013-06-tufts and
walked the audience through forking, committing, and pushing 
to the twenty forks of the toy repository.

I'm not a fan of Fork and Spoon; if we could use the boot-camps
repository for examples more this would be great.

I tried the same approach at 2013-06-chicago and 2013-07-notredame
and asked for pull requests to update the upstream repository
with content.  I suggested that there might be a universe in 
which students submit their homework by pull request.

By the end of the class at Notre Dame, most of the class had 
sent pull requests, but it looks like only a few updated their
forks with all the changes.

Requirements:
Instructor sets up a toy repository with placeholder content.

Teaching assistant has access to the toy repository and can approve
pull requests (after the instructor approves the first one or two)
so that the instructor's toy repository gets commits from most
of the class.

Hint:  It is helpful to have a TA give the instructor push access to 
the TA's fork of the repository.  This allows the instructor to do 
exactly as the students -- updating a repository (the TA's)  that is 
NOT upstream.

origin    https://github.com/TA/testrepo-YYYY-MM-PLACE
upstream  https://github.com/INSTRUCTOR/testrepo-YYYY-MM-PLACE

the placeholders USERNAME and TESTREPO should be replaced by 
the instructor's username and test repository name.

We used 
https://github.com/ahmadia/bio-pipeline
https://github.com/wltrimbl/testrepo-2013-06-chicago
https://github.com/wltrimbl/testrepo-2013-07-notredame

W. Trimble 2013-07-31

