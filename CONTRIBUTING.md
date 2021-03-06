How to contribute
=================

Contributions to Astrochem are welcome. Please follow these guidelines
if you would like to contribute.

Getting started
---------------

* Make sure you have a [GitHub account](https://github.com/signup/free).
* [Fork](http://guides.github.com/activities/forking/) the astrochem repository.

Making changes
--------------

* [Create a topic branch](http://guides.github.com/introduction/flow/) for
  your commits.
* Make commits of logical units.
* Check for unnecessary whitespace with `git diff --check` before
  committing.
* Make sure your commit messages are in the proper format:

```
Short (50 chars or less) summary of changes

More detailed explanatory text, if necessary.  Wrap it to
about 72 characters or so.  In some contexts, the first
line is treated as the subject of an email and the rest of
the text as the body.  The blank line separating the
summary from the body is critical (unless you omit the body
entirely); tools like rebase can get confused if you run
the two together. If the commit fixes a bug, the description
should contain the bug number, for example: "fixes #3".

```	

* Make sure you have added the necessary tests for your changes.
* Run the tests with `make check` to assure nothing else was
  accidentally broken.

Submitting changes
------------------

* [Rebase](https://help.github.com/articles/about-git-rebase/) your
  work on Astrochem master branch to ease the merging and help us
  keeping a clean git history. If you made many small, trivial commits
  during the development, please
  [squash](http://gitready.com/advanced/2009/02/10/squashing-commits-with-rebase.html)
  them into logical units with an interactive rebasing.
* [Push](https://help.github.com/articles/pushing-to-a-remote/) your
  changes to the topic branch in your fork of the repository.
* Submit a
  [pull request](https://help.github.com/articles/creating-a-pull-request/)
  to the astrochem repository.
* Wait for your changes to be reviewed.

Coding style
------------

Format your source code in the same style than the rest of astrochem
code; a mixture of formating styles tends to make the code difficult
to read. For C code, we use the
[GNU style](http://www.gnu.org/prep/standards/standards.html#Formatting).
For Python code, we follow the recommendations in the
[Style Guide for Python code](http://legacy.python.org/dev/peps/pep-0008/).
