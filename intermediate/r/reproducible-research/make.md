# Make - automate document writing

`Make` is a simple but powerful tool for managing a build process in a language-independent manner. 

`Make` is a domain-specific language for encoding a dependence graph.

By writing a make file you are instructing the program make to synthesize files from dependencies.

`Make` is widely available in many platforms and can easily be installed in others where it's not available as default (e.g. Windows).

---

# Basic idea of make

file-to-generate : dependency (or dependencies)
    instruction on how to generate the file

Note: The instructions should be indented by a single `<tab>`.


# Process multiple files at once
`$<` are the input files. I.e the .md file in a pandoc command
`$@` is the output file.


# Example make files

* [Makefile from Ethan White's paper](https://github.com/weecology/data-sharing-paper/blob/master/makefile)

## References

* Some notes from [here](http://matt.might.net/articles/intro-to-make/)
