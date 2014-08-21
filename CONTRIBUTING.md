How to contribute
=================

Contributions to Astrochem are welcome. Please follow these guidelines
if you would like to contribute.

Getting started
---------------

* Make sure you have a [GitHub account](https://github.com/signup/free).
* Fork the astrochem repository on GitHub.

Making changes
--------------

* Create a topic branch from where you want to base your work.
* Make commits of logical units.
* Check for unnecessary whitespace with `git diff --check` before
  committing.
* Make sure your commit messages are in the proper format:

```
Summary of the commit (less than 70 chars)

A detailed description of the commit. It can be omitted for
trivial changes. The description can span over several lines. If
the commit fixes a bug, the description should contain the bug
number, for example: "fixes #3".
```	

* Make sure you have added the necessary tests for your changes.
* Run the tests with `make check` to assure nothing else was
  accidentally broken.

Submitting changes
------------------

* Rebase your work on Astrochem master branch to ease the merging and
  help us keeping a clean git history. If you made many small, trivial
  commits during the development, please squash them into logical
  units with an interactive rebasing.
* Push your changes to a topic branch in your fork of the repository.
* Submit a pull request to the Astrochem repository.

Coding style
------------

Format your source code in the same style than the rest of astrochem
code; a mixture of formating styles tends to make the code difficult
to read. For C code, we use the
[GNU style](http://www.gnu.org/prep/standards/standards.html#Formatting).
For Python code, we follow the recommendations in the
[Style Guide for Python code](http://legacy.python.org/dev/peps/pep-0008/).
