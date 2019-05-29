.. _astropy-style-guide:

******************************************************************
Astropy Narrative Style Guide: A Writing Resource for Contributors
******************************************************************

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

For general use and scientific terms, use the full term instead of the
abbreviation, for example, write out white dwarf instead of WD, etc.

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

In Documentation materials, up style capitalization is preferred in headings,
meaning capitalize first, last, and all major words in the heading, but
lowercase articles (the, a, an), prepositions (at, to, up, down, with, in,
etc.), and common coordinating conjunctions (and, but, for, or). Down style is
acceptable for longer “example” headings.

Examples
--------
* Building and Installing
* Frames without Data
* Checklist for Contributing Code
* Astropy Guidelines
* Importing ``astropy`` and Sub-packages
* Example: Use velocity to compute sky position at different epochs

In Tutorials and other learning materials, up style capitalization is preferred
in headings of structured introductory/template sections, but within the
tutorial, down style (i.e., capitalize first word and proper nouns only) is
acceptable for longer headings designating different learning/code sections.

Contractions
============

Do not use contractions in formal documentation material.

Examples
--------
* If you are making changes that impact ``astropy`` performance, consider adding
  a performance benchmark.
* You do not need to include a changelog entry.

In all other materials, avoid use of contractions only when tense can be
confused, such as in the case of “she is gone” versus “she has gone,” etc.

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
* Since most of the functionality of ``astropy`` resides in sub-packages, it is
  best to import the desired sub-package.

Numbers
=======

For numbers followed by a unit or as part of a name, use the numeral.

Examples
--------
* 1 arcminute
* 32 degrees
* Gaia data release 2 catalog

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

.. note::

    Please note that en dashes and em dashes do not render on GitHub and appear
    as hyphens. This section is included for reference only.

Hyphens (-) should be used for phrasal adjectives and compound words (see
Hyphenation above).

En dashes (– longer) should be used for number ranges (dates, times, pages) or
to replace the words “to” or “through,” without spaces around the dash.

Examples
--------
* See chapters 14–18
* We’ve blocked off March 2019–May 2019 to develop a new version.

Em dashes (— longest) can be used in place of commas, parentheses, or colons to
set off amplifying or explanatory elements. In Astropy materials, follow AP
style, which calls for spaces on either side of each em dash.

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

Use numerals when exact times are expressed, followed by *ante meridiem* or
*post meridiem* abbreviated in lowercase with periods, or use the twenty-four-
hour system.

Examples
--------
* The meeting is set for 9:30 a.m.
* The presentation starts at 15:00.

Express specific dates as numerals in ISO 8601 format, year-month-day.

Example
-------
* Data from the Gaia mission was released on 2018-04-25.

A Note About Style and Tone
===========================

Across all Astropy materials in narrative sections, please write in:

* The present tense. For example:
  * In the following section, we are going to make a plot…
  * To test if your version of ``astropy`` is running correctly…

* The first-person inclusive plural. For example:
  * We did this the long way, but next let’s try the short way…

* Use the generic pronoun “you” instead of “one.” For example:
  * You can access any of the attributes on a frame by...

Always avoid extraneous or belittling words such as “obviously,” “easily,”
“simply,” “just,” or “straightforward.” Avoid extraneous phrases like, “we just
have to do one more thing.”

Documentation vs. Tutorials vs. Guides
--------------------------------------

Documentation
~~~~~~~~~~~~~
Tone: academic and slightly more formal.

* Use up style capitalization in section headings.
* Do not use contractions.

Tutorials
~~~~~~~~~
Tone: academic but less formal and more friendly.

* Use up style capitalization in introductory/template headings, switch to down
  style capitalization for learning/example section headings.
* Section headings should use the imperative mood to form a command or request
  (e.g., “Download the data”).
* Contractions can be used as long as the tense is clear.

Guides
~~~~~~
Tone: academic but less formal and more friendly.

* Use up style capitalization in introductory/template headings, switch to down
  style capitalization for learning/example section headings.
* Contractions can be used as long as the tense is clear.
