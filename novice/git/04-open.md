---
layout: lesson
root: ../..
title: Open Science
---
<div class="objectives" markdown="1">

#### Objectives
* Learn how to distribute open source software:
  - Choosing an appropriate open source license
  - Choosing an appropriate hosting repository
* Learn how to distribute open data:
  - Understand licensing concerns for data and metadata
  - Choosing an appropriate data repository

</div>


Knowing how to effectively publish and distribute open source software 
and open data is becoming as important to scientific research as publishing
papers -- indeed, it is already required by many of the most prestigious
journals. In this lesson we focus on the two key components to publishing
data or source code: licensing and repositories. 


## Open Source Software ## 

### Licensing Software ###

Open source licenses assist the creator of a creative work in waiving
some of the rights and privileges which they are automatically granted
under [_copyright_ law](http://en.wikipedia.org/wiki/Copyright).


Broadly speaking, there are two kinds of open license for
software: [copyleft](http://en.wikipedia.org/wiki/Copyleft)
licenses such as the [GNU General Public
Licenses](http://opensource.org/licenses/GPL-3.0) (GPL),  and
[permissive](http://en.wikipedia.org/wiki/Permissive_free_software_licence)
licenses such as the [MIT](http://opensource.org/licenses/MIT) and
[BSD](http://opensource.org/licenses/BSD-2-Clause) licenses.  All of these
licenses allow unrestricted sharing and modification of programs, but
copyleft licenses are [infective](../../gloss.html#infective-license):
anyone who distributes a modified version of the code (or anything
that includes GPL'd code) must make *their* code freely available as
well. Lacking this clause, code under permissive licenses can be more
immediately used in commercial software.

#### How to apply a license ####

Before releasing open source software you should confirm with your 
employer that you are the current copyright holder (in academic settings,
faculty tend to control their own copyrights while the copyrights of work
done by staff often belong to the university).

Software licenses are typically applied by including a plain-text file
with name such as `LICENSE` or `COPYING` in the project directory.
Some projects will place the full text of the license in comments at
the top of every source file, while others may only declare the choice
of license by an abbreviation and/or a link to the license terms. 

The legal text for most open source licenses can be found from the [Open
Source Initiative](http://opensource.org/), which maintains a list of
open source licenses which have gone through their approval process.
[tl;drLegal](http://www.tldrlegal.com/) explains many of them in plain
English.  

When selecting a license, be sure that your choice is consistent with
the terms of any software you may be reusing or modifying (usually by 
adopting the license already in use). Note that many licenses have 
multiple versions which are not necessarily compatible, so be sure to
be explicit. 


------------------

[Software Carpentry](http://software-carpentry.org/license.html) uses
CC-BY for its lessons and the MIT License for its code in order to
encourage the widest possible re-use.  Again, the most important thing
is for the `LICENSE` file in the root directory of your project to state
clearly what your license is.  You may also want to include a file called
`CITATION` or `CITATION.txt` that describes how to reference your project;
the one for Software Carpentry states:

<div class="file" markdown="1">
~~~
To reference Software Carpentry in publications, please cite both of the following:

Greg Wilson: "Software Carpentry: Lessons Learned". arXiv:1307.5448, July 2013.

@online{wilson-software-carpentry-2013,
  author      = {Greg Wilson},
  title       = {Software Carpentry: Lessons Learned},
  version     = {1},
  date        = {2013-07-20},
  eprinttype  = {arxiv},
  eprint      = {1307.5448}
}
~~~
</div>

### Hosting & Distributing Software ###

Open Source research software is best distributed through the use of a
dedicated code repository or academic data archive.  Most (but not all)
code repositories are built around the use of a version control system
such as `git` or `subversion`, which creates some barrier to entry
(fortunately you've just completed the `git` SWC lessons!)

Public hosting services such as [GitHub](http://github.com),
[BitBucket](http://bitbucket.org), [Google Code](http://code.google.com),
or [SourceForge](http://sourceforge.net) are feature rich,
user friendly and widely adopted options. All provide free
hosting for open-source projects (and usually a limited number
of free private projects as well).  See other recommendations
for code repositories from the [Journal for Open Research
Software](http://openresearchsoftware.metajnl.com/about/editorialPolicies#custom-0).

Researchers may also choose to distribute software through dedicated
language repositories such as CRAN (R). These language-specific
repositories host only code that is ready for use and will usually make
it easier for other users to install your software. These repositories
also archive versions as they are released, but typically do not
require using version management software. These repositories often
have stricter criteria than the public hosting services described
above, so be sure to consult the appropriate policies (e.g. [CRAN
policies](cran.r-project.org/web/packages/policies.html)) before
proceeding.  Many projects host their daily development on public hosting
services while also distributing releases through a system such as this.


It has been common practice for researchers to host software they develop
on computer servers managed by their lab, department, or institution.
Experience has shown that software and other resources hosted in this
fashion has a much higher rate of link rot, where changes to websites,
changing jobs, or other factors make it unlikely that these resources are
still available years later.  These options also typically lack many of
the features dedicated software repositories provide. Online supplement
sections of journals are also not ideal mechanisms to distribute
software, for many of the same reasons. It is best to simply link from
your publication or personal website to the permanent software repository.

make it easy to run a Github-like environment on a private server but
may not be as well suited for long-term hosting as the larger dedicated
hosting services.  -->

Some scientific data archives will also host software.  Because
these archives are backed by long-term redundant archiving
(e.g. [CLOCKSS](http://clockss.org/)) and permanent identifiers
(e.g. [DOIs](http://en.wikipedia.org/wiki/Digital_object_identifier)),
they offer a more long-term archival storage solution (see Archiving
Data, below).  The data repositories [zenodo](https://zenodo.org) and
[figshare](http://figshare.com) currently have automated [integration
with Github](http://collaborate.mozillascience.org/projects/codemeta)
to facilitate this.


## Open Data ##

Learners should know how to publish open data effectively, whether
or not they choose to do so in any particular circumstance. -->

### Licensing Data ###

Unlike software or other creative works, data are considered facts and
generally not subject to copyright.  Many academic data repositories
underscore this by requiring a public-domain declaration such as
[Creative Commons Zero](https://creativecommons.org/about/cc0) (or CC0,
not technically a license) for all data that they host (see [Panton
Principles](http://pantonprinciples.org) of open data.) Even when placing
data or other work in the public domain it is preferable to use a standard
declaration such as CC0, since writing an internationally valid legal
document is a task best left to the relevant experts. 

Data formats, descriptions, or databases are considered creative works
and are frequently accompanied by a copyright statement. Creative Commons
provides a suite of licenses to waive various aspects of copyright in
order to facilitate open reuse.  The most permissive of these is the
_Attribution_ or CC-BY license. Alternatives may restrict commercial use
(NC: non-commercial), restrict derivative products (ND: no-derivatives),
or include the copyleft clause (SA: share-alike) similar to the GPL. Different
licenses offer any combination of the latter three clauses on top of the 
default BY clause.

It is worth noting that only the CC-BY license
is considered compatible with the widely recognized Budapest Open Access
Initiative [definition](en.wikipedia.org/wiki/Budapest_Open_Access_Initiative).
Several studies have shown that researchers choose the more restrictive
variations by default and are unaware of the limitations this may
place on uses they condone, such as education.


### Archiving & Distributing Data ###

Many journals now require authors to deposit all data supporting published
results into a scientific data repository. As with software repositories,
data repositories are better suited for sharing data than hosting on one's own
website or in a journal's supplemental online materials.

Scientific data repositories may be divided into two types: those
accepting only published data accompanying a scientific article
(e.g. [Dryad](http://datadryad.org), and those that also accept data
that is not (or not yet) associated with any particular publication
([Zenodo](http://zenodo.org), [figshare](http://figshare.com)). Some
repositories focus on narrow subject areas or data types, while others
are more general purpose.  Consult the policies of your journals,
discipline-specific literature on data archiving, and the policies of the
data archives themselves in finding a good match. The [recommendations
from _Nature_](http://www.nature.com/sdata/data-policies/repositories)
are one good place to start.


Data repositories provide many advantages, including: 


- **Permanent identifiers:** Though widely touted as making your data
'citable', permanent identifiers are designed to avoid link rot that
results from changing URLs (hence the name). The [Digital Object
Identifier](http://en.wikipedia.org/wiki/Digital_object_identifier),
or DOI is the best known because of its association with scientific
publications.[^1]   An object with a DOI number can be found by entering
the number into a central registry, [http://doi.org](http://doi.org),
regardless of the URL address currently hosting it. Repositories must pay
a small fee for each DOI. If a repository fails to update the records
allowing the DOI to resolve to the correct resource, the DOI provider
may refuse to sell them additional DOIs.

- **Metadata & data discovery**: Data repositories collect basic metadata
such as author and subject information. This facilitates search and
discovery of relevant datasets. DOI-based repositories submit much
of this information in a standardized format to the central registry
at DataCite, which allows tools and researchers to search across all
DataCite repositories at once.

- **Data management** Data repositories are well equipped to
provide redundant and reliable access to data over the long
term. Data can be updated or corrected while maintaining links
to the original versions. Looks great on [Data management
plans](http://www.nsf.gov/eng/general/dmp.jsp).

[^1]: Technically Data DOIs are different than scientific publication
DOIs, in that the former are administered by DataCite and the latter by
CrossRef, and as such include slightly different metadata and protocols.

#### Special cases ####

Data security concerns are not a good reason to be lazy about data archiving.
Sensitive data (e.g. human experimental subjects) should always be
dealt with as such, following appropriate anonymization and/or security
protocols defined before the data is collected. Many repositories have
explicit mechanisms in place to to handle sensitive data appropriately.
Storing sensitive data on personal machines without clear security policies in 
place may be inappropriate.

Rapidly updated, streaming, or very large datasets (usually >2-10
GB) still pose challenges for most general purpose scientific data
repositories.

<div class="keypoints" markdown="1">

## Key Points ##

* Open source licenses include both permissive (BSD, MIT) and copyleft
(GPL) style licenses. Anyone distributing software with code taken from
or modified from code under a GPL style-license must make their derivative
source code available under the same terms.

* Open data should be placed in the public domain using the CC0
declaration (copyright not being applicable to facts).

* Dedicated software repositories such as [GitHub](http://github.com),
[BitBucket](http://bitbucket.org), [Google Code](http://code.google.com),
or [SourceForge](http://sourceforge.net) are preferable to self hosting software.

- Other creative works, including data descriptions and publications, can
use Creative Commons licenses to facilitate reuse. The most permissive
license, CC-BY, corresponds with community definitions of Open Access,
while others are more restrictive.

- Dedicated scientific data repositories, such as those integrated with
DataCite (e.g. any that provide DOIs), are the preferable mechanism for
data archiving.

</div>

<div class="challenge" markdown="1">
Find out whether you are allowed to apply an open license to your software.
Can you do this unilaterally,
or do you need permission from someone in your institution?
If so, who?
</div>

<div class="challenge" markdown="1">
Find out whether you are allowed to host your work openly on a public forge.
Can you do this unilaterally,
or do you need permission from someone in your institution?
If so, who?
</div>
