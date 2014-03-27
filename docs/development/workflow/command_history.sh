#
# Setup
#
# Start by cloning the repo from your GitHub account if you need to...
git clone https://github.com/mwcraig/astropy.git apy-bugs
# ...then navigate to that clone in a terminal, and...
cd apy-bugs
# ...rename the origin remote make its meaning clearer and...
git remote rename origin mwcraig
# ...remember to set up the official remote.
git remote add astropy git://github.com/astropy/astropy.git
# get the latest from the official repo
git fetch astropy
#
# Make a new branch for this fix
#
git branch i-1761 astropy/master # branch just for this work
git checkout i-1761               # work on this branch
git push --set-upstream mwcraig i-1761 # Hey GitHub, get ready to receive my work
#
# Make a python environment for this work -- conda version shown
#
conda create -n apy-1761 --clone root   # copy my default environment
source activate apy-1761                # switch to this environment
# next step DOES NOT WORK in python 3
python setup.py develop                 # install astropy
#
# Check setup by running tests
#
cd astropy/coordinates/tests    # get ready to run coordinate tests
py.test                         # test
ls                              # what is here?
#
# Write a test in test_arrays.py to expose this bug
# 
# After edit, re-test
#
py.test test_arrays.py          # Hopefully this succeeds
#
# Commit changes to local git
#
git status                      # what did we do? Remind me, git
git add test_arrays.py          # get ready to commit
# always include a message with your commits
git commit -m'Add test for array coordinate length (issue #1761)'
git status                      # Best to be informed!
#
# Edit test_arrays.py
#
cd ..                           # up a level to coordinates
#
# Code fix made in editor
#
#
# Test change
#
py.test tests/test_arrays.py   # Did our fix actually fix the problem?
py.test                        # do all of the coordinate tests pass?
cd ../..                       # up to the top level to check all of the tests
python setup.py test           # grab a coffee now; this takes a while
#
# Success! Commit our change
#
git status                     # check what we have changed, if you forgot
git add astropy/coordinates/coordsystems.py  # Copy and paste path from git status
                                             # unless you like typing
# check git staus again if you want...or not.
# commit, with message!
git commit -m"
    > Add len() to coordinates
    >
    > Closes #1761"
git push                       # or not--can wait until you are completely done
#
# Oops! Forgot to include some tests!
#
cd astropy/coordinates/tests   # back to add more tests
#
# edit test_arrays.py to add tests to test_array_len
#
# Now re-test
#
py.test test_arrays.py  # do these tests pass with the new tests?
#   yes! now move up a couple of levels to check all tests
#
#   Do I really have to re-test everything? It certainly never hurts and
#   potentially saves time for the maintainers if you have introduced a bug.
#
cd ../../..          # back to the top level
python setup.py test # grab another coffee....
#
# Success! Time to commit
#
git status           # optional, but useful for figuring out what to add
git add astropy/coordinates/tests/test_arrays.py # use path git status supplies
# Commit, with a message...
git commit -m"Add tests of len() for scalar coordinate and length 1 coordinate"
#
# Edit the changelog in CHANGES.rst
#
# Then commit changes
#
git add CHANGES.rst
git commit -m"Add changelog entry for 1761"
#
# Push changes to my GitHub account
#
git push   # make sure all of our changes are in our github account
#
# Go to GitHub to make pull request. Remember to switch to this branch before 
# pushing the pull request button.
#