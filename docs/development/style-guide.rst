.. _astropy-style-guide:

******************************************************************
Astropy Narrative Style Guide: A Writing Resource for Contributors
******************************************************************

The purpose of this style guide is to provide the Astropy community with a set
of style and formatting guidelines that can be referenced when writing Astropy
documentation. Following the guidelines offered in this style guide will bring
greater consistency and clarity to Astropy's documentation, supporting its
mission to develop a common core package for Astronomy in Python and foster an
ecosystem of interoperable astronomy packages.

This style guide is organized alphabetically by writing topic, with usage
examples in each section, and tone and formatting guidelines at the end.

Abbreviations
=============

Place abbreviations such as i.e. and e.g. within parentheses, where they are
followed by a comma. Alternatively, consider using "that is" and “for example”
instead, preceded by an em dash or semicolon and followed by a comma, or
contained within em dashes.

Examples
--------
* The only way to modify the data in a frame is by using the ``data`` attribute
  directly and not the aliases for components on the frame (i.e., the following
  will not work).
* There are no plans to support more complex evolution (e.g., non-inertial
  frames or more complex evolution), as that is out of scope for the ``astropy``
  core.
* Once you have a coordinate object you can access the components of that
  coordinate — for example, RA or Dec — to get string representations of the
  full coordinate.

For general use and scientific terms, use abbreviations only when the
abbreviated term is well-known and widely used within the astronomy community.
For less common scientific terms, or terms specific to a given field, write out
the term or link to a resource of explanation. A good rule of thumb to follow
when deciding whether or not something should be abbreviated is: when in doubt,
write it out.

Examples
--------
* 1D, 2D, etc. is preferred over one-dimensional, two-dimensional, etc.
* Units such as SI and CGS can be abbreviated as is more commonly seen in the
  scientific community.
* White dwarf should be written out fully instead of abbreviated as WD.
* Names of organizations or other proper nouns that employ acronyms should be
  written as their known acronym, but with a hyperlink to a website or resource
  for reference, for instance, `CODATA <https://codata.org/>`_.

Capitalization
==============

Capitalize all proper nouns (names) in plain text, except when referring to
package/code names, in which case use lowercase and double backticks. Astropy
capitalized refers to The Astropy Project, while ``astropy`` lowercase and in
backticks refers to the core package.

Examples
--------
* Follow Astropy guidelines for contributing code.
* Affiliated packages are astronomy-related software packages that are not part
  of the ``astropy`` core package.
* Provide a code example along with details of the operating system and the
  Python, ``numpy``, and ``astropy`` versions you are using.

In Documentation materials, title case capitalization is preferred in headings,
meaning capitalize first, last, and all major words in the heading, but
lowercase articles (the, a, an), prepositions (at, to, up, down, with, in,
etc.), and common coordinating conjunctions (and, but, for, or). Sentence case
capitalization is acceptable for longer example headings.

Examples
--------
* Building and Installing
* Frames without Data
* Checklist for Contributing Code
* Astropy Guidelines
* Importing ``astropy`` and Subpackages
* Example: Use velocity to compute sky position at different epochs

In Tutorials and other learning materials, title case capitalization is
preferred in headings of structured introductory/template sections, but within
the tutorial, sentence case (i.e., capitalize first word and proper nouns only)
is acceptable for longer headings designating different learning/code sections.

Contractions
============

Do not use contractions in formal documentation material.

Examples
--------
* If you are making changes that impact ``astropy`` performance, consider adding
  a performance benchmark.
* You do not need to include a changelog entry.

In all other materials, avoid use of contractions only when the tense can be
confused, such as in the case of “she is gone” versus “she has gone,” etc.

.. _Hyphenation:

Hyphenation
===========

Phrasal adjectives/compound modifiers placed before a noun should be hyphenated
to avoid confusion.

Examples
--------
* Astronomy-related software packages.
* Astropy provides sustainable, high-level education to the astronomy community.

Hyphenated compound words should contain hyphens in plain text, but no hyphens
in code.

Example
-------
* Do not forget to double-check your formatting.

Numbers
=======

For numbers followed by a unit or as part of a name, use the numeral.

Examples
--------
* 1 arcminute
* 32 degrees
* Gaia data release 2 catalog
* 1D, 2D, etc. is preferred over one-dimensional, two-dimensional, etc.

For all other whole numbers, follow Associated Press (AP) style: spell out
numbers one through nine, and use numerals for 10 and higher, with numeral-word
combinations for millions, billions, and trillions.

Examples
--------
* There are two ways to build Astropy documentation.
* Follow these 11 steps.
* Measuring astrometry for about 2 billion stars.

For casual expressions, spell out the number.

Example
-------
* A picture is worth a thousand words.

Punctuation
===========

For consistency across Astropy materials, non-U.S. punctuation will be edited
to reflect American punctuation preferences.

**Parentheses**: punctuation belonging to parenthetical material will be placed
inside of closing parentheses, with the exception of commas to denote a small
pause coming after parenthetical material, and periods when parenthetical
material is included within another sentence.

Examples
--------
* (For full contributor guidelines, see our documentation.)
* Once you open a pull request (which should be opened against the ``master``
  branch), please make sure to include the following.
* In some cases, most of the required functionality is contained in a single
  class (or a few classes).

**Quotation marks**: periods and commas will be placed inside of closing
quotation marks, whether double or single.

Examples
--------
* Chief among these terms is the concept of a “coordinate system.”
* Because of the likelihood of confusion between these meanings of “coordinate
  system,” `~astropy.coordinates` avoids this term wherever possible.

**Hyphens vs. En Dashes vs. Em Dashes**

Hyphens (-) should be used for phrasal adjectives and compound words (see
`Hyphenation`_ above).

En dashes (– longer) should be used for number ranges (dates, times, pages) or
to replace the words “to” or “through,” without spaces around the dash.

Examples
--------
* See chapters 14–18.
* We have blocked off March 2019–May 2019 to develop a new version.

Em dashes (— longest) can be used in place of commas, parentheses, or colons to
set off amplifying or explanatory elements. In Astropy materials, follow
Associated Press (AP) style, which calls for spaces on either side of each em
dash.

Examples
--------
* Several types of input angles — array, scalar, tuple, string — can be used in
  the creation of an Angle object.
* The creation of an Angle object supports a variety of input angle types —
  array, scalar, tuple, string, etc.

Spelling
========

For consistency across Astropy materials, non-U.S. spelling will be edited to
reflect American spelling preferences.

Example
-------
* Cross-matching catalog coordinates (versus catalogue)

Time and Date
=============

Use numerals when exact times are expressed. Use the 24-hour system to express
exact times. For consistency across Astropy materials, all instances of exact
times will be edited to reflect 24-hour time system preferences.

Example
-------
* The presentation starts at 15:00.

Express specific dates as numerals in ISO 8601 format, year-month-day.

Example
-------
* Data from the Gaia mission was released on 2018-04-25.

A Note About Voice and Tone
===========================

Across all Astropy materials in narrative sections, please follow these voice
and tone guidelines.

Write in the present tense.

Example
-------
* In the following section, we are going to make a plot...
* To test if your version of ``astropy`` is running correctly...

Use the first-person inclusive plural.

Example
-------
* We did this the long way, but next we can try it the short way...

Use the generic pronoun “you” instead of “one.”

Example
-------
* You can access any of the attributes on a frame by...

Always avoid extraneous or belittling words such as “obviously,” “easily,”
“simply,” “just,” or “straightforward.” Avoid extraneous phrases like, “we just
have to do one more thing.”

Avoid words or phrases that create worry in the mind of the reader. Instead,
use positive language that establishes confidence in the skills being learned.

Examples
--------
* As a best practice...
* One recommended way to...
* An important note to remember is...

Along these lines, use "warning" directives only to note limitations in the
code, not implied limitations in the skills or knowledge of the reader.

Documentation vs. Tutorials vs. Guides
--------------------------------------

Documentation
^^^^^^^^^^^^^
Tone: academic and slightly more formal.

* Use title case capitalization in section headings.
* Do not use contractions.

Tutorials
^^^^^^^^^
Tone: academic but less formal and more friendly.

* Use title case capitalization in introductory/template headings, switch to
  sentence case capitalization for learning/example section headings.
* Section headings should use the imperative mood to form a command or request
  (e.g., “Download the data”).
* Contractions can be used as long as the tense is clear.

Guides
^^^^^^
Tone: academic but less formal and more friendly.

* Use title case capitalization in introductory/template headings, switch to
  sentence case capitalization for learning/example section headings.
* Contractions can be used as long as the tense is clear.

Formatting Guidelines
=====================

Astropy documentation is written in reStructuredText using the Sphinx
documentation generator. When formatting the different sections of your
documentation files, please follow these guidelines to maintain consistency in
section heading hierarchy across Astropy's RST files.

Section headings in reStructuredText files are created by underlining (and
optionally overlining) the section title with a punctuation character the same
length as the text.

Examples
--------

::

  *************************
  This is a Chapter Heading
  *************************

::

  This is a Section Heading
  =========================

Although there are no formally assigned characters to create heading level
hierarchy, as the hierarchy rendering is determined from the succession of
headings, here is a suggested convention to follow when formatting Astropy
documentation files:

# with overline, for parts
* with overline, for chapters
=, for sections
-, for subsections
^, for subsubsections
", for paragraphs

These guidelines follow Sphinx's recommendation in the `Sections
<https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html#sections>`_
chapter of its reStructuredText Primer and Python's convention in the `7.3.6.
Sections <https://devguide.python.org/documenting/#sections>`_ part of its style
guide.

Other Writing Resources
=======================

Some other resources that may be useful when writing Astropy documentation are:

* Python's `Style Guide
  <https://devguide.python.org/documenting/#style-guide>`_
* Sphinx's `reStructuredText Primer
  <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_
* `Quick reStructuredText
  <https://docutils.sourceforge.io/docs/user/rst/quickref.html>`_
