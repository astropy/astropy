#
# Setup
#
# Start by cloning the repo from your GitHub account if you need to...
$ git clone https://github.com/mwcraig/astropy.git apy-bugs
Cloning into 'apy-bugs'...                                                      
remote: Counting objects: 44980, done.                                          
remote: Compressing objects: 100% (11772/11772), done.                          
remote: Total 44980 (delta 32930), reused 44980 (delta 32930)                   
Receiving objects: 100% (44980/44980), 18.13 MiB | 3.57 MiB/s, done.            
Resolving deltas: 100% (32930/32930), done.                                     
Checking connectivity... done.  

# ...then navigate to that clone in a terminal, and...
$ cd apy-bugs
# ...rename the origin remote make its meaning clearer and...
$ git remote rename origin mwcraig
# ...remember to set up the official remote.
$ git remote add astropy git://github.com/astropy/astropy.git

# get the latest from the official repo
$ git fetch astropy
remote: Counting objects: 37, done.                                             
remote: Compressing objects: 100% (30/30), done.                                
remote: Total 37 (delta 10), reused 23 (delta 7)                                
Unpacking objects: 100% (37/37), done.                                          
From git://github.com/astropy/astropy                                           
 * [new branch]      master     -> astropy/master                               
 * [new branch]      stable     -> astropy/stable                               
 * [new branch]      v0.1.x     -> astropy/v0.1.x                               
 * [new branch]      v0.2.x     -> astropy/v0.2.x                               
 * [new branch]      v0.3.x     -> astropy/v0.3.x  
#
# Make a new branch for this fix
#
$ git branch i-1761 astropy/master # branch just for this work
Branch i-1761 set up to track remote branch master from astropy.

$ git checkout i-1761               # work on this branch
Switched to branch 'i-1761'                                                     
Your branch is up-to-date with 'astropy/master'. 

# let my github account know about this branch
$ git push --set-upstream mwcraig i-1761 # Hey GitHub, get ready to receive my work
Branch i-1761 set up to track remote branch i-1761 from mwcraig.                
Everything up-to-date  

#
# Make a python environment for this work -- conda version shown
#
$ conda create -n apy-1761 --clone root   # copy my default environment
Output not shown

$ source activate apy-1761                # switch to this environment
prepending /Users/mcraig/anaconda/envs/apy-1761/bin to PATH

# next step DOES NOT WORK in python 3
$ python setup.py develop   # install astropy, python 2
# using python 3? do this instead: python3 setup.py install
Configured with: --prefix=/Applications/Xcode.app/Contents/Developer/usr --with-
gxx-include-dir=/usr/include/c++/4.2.1                                          
running develop                                                                 
running egg_info                                                                
writing requirements to astropy.egg-info/requires.txt                           
writing astropy.egg-info/PKG-INFO
...
Using /Users/mcraig/anaconda/envs/apy-1761/lib/python2.7/site-packages          
Finished processing dependencies for astropy==0.4.dev6873  

#
# Check setup by running tests
#
$ cd astropy/coordinates/tests    # get ready to run coordinate tests
$ py.test                         # test
============================= test session starts ==============================
platform darwin -- Python 2.7.6 -- pytest-2.4.2

Running tests with Astropy version 0.4.dev6877.
Running tests in /Users/mcraig/Development/astronomy/apy-bugs/astropy/coordinates/tests.

Platform: Darwin-13.0.0-x86_64-i386-64bit

Executable: /Users/mcraig/anaconda/envs/apy-1761/bin/python

Full Python Version:
2.7.6 |Anaconda 1.8.0 (x86_64)| (default, Nov 11 2013, 10:49:09)
[GCC 4.0.1 (Apple Inc. build 5493)]

encodings: sys: ascii, locale: US-ASCII, filesystem: utf-8, unicode bits: 15
byteorder: little
float info: dig: 15, mant_dig: 15

Numpy: 1.7.1
Scipy: 0.13.0
Matplotlib: 1.3.1
h5py: 2.2.0

plugins: capturelog
collected 103 items / 2 skipped

test_angles.py .......................
test_angular_separation.py ..
test_api.py ............
test_arrays.py .............
test_distance.py ........
test_formatting.py ....................
test_matching.py ....
test_name_resolve.py ss
test_transformations.py ...............
accuracy/test_fk4_no_e_fk4.py .
accuracy/test_fk4_no_e_fk5.py .
accuracy/test_galactic_fk4.py .
accuracy/test_icrs_fk5.py .

=================== 101 passed, 4 skipped in 152.23 seconds ===================

$ ls   # what is here?
__init__.py                test_api.pyc                                         
__init__.pyc               test_arrays.py                                       
__pycache__                test_distance.py                                     
accuracy                   test_formatting.py                                   
test_angles.py             test_matching.py                                     
test_angular_separation.py test_name_resolve.py                                 
test_api.py                test_transformations.py  
#
# Write a test in test_arrays.py to expose this bug
# 
# After edit, re-test
#
$ py.test test_arrays.py # Hopefully this FAILS--we are trying to expose bug
                                                                                
Running tests with Astropy version 0.4.dev6873.                                 
Running tests in test_arrays.py.                                                
                                                                                
Platform: Darwin-13.0.0-x86_64-i386-64bit                                       
                                                                                
Executable: /Users/mcraig/anaconda/envs/apy-1761/bin/python                     
                                                                                
Full Python Version:                                                            
2.7.6 |Anaconda 1.8.0 (x86_64)| (default, Nov 11 2013, 10:49:09)                
[GCC 4.0.1 (Apple Inc. build 5493)]                                             
                                                                                
encodings: sys: ascii, locale: US-ASCII, filesystem: utf-8, unicode bits: 15    
byteorder: little                                                               
float info: dig: 15, mant_dig: 15                                               
                                                                                
Numpy: 1.7.1                                                                    
Scipy: 0.13.0                                                                   
Matplotlib: 1.3.1                                                               
h5py: 2.2.0                                                                     
                                                                                
plugins: capturelog                                                             
collected 13 items                                                              
                                                                                
test_arrays.py ............F                                                    
                                                                                
=================================== FAILURES ===================================
________________________________ test_array_len ________________________________
                                                                                
    def test_array_len():                                                       
        from .. import ICRS                                                     
                                                                                
        input_length = 5                                                        
        ra = np.linspace(0, 360, input_length)                                  
        dec = np.linspace(0, 90, input_length)                                  
                                                                                
        c = ICRS(ra, dec, unit=(u.degree, u.degree))                            
                                                                                
>       assert len(c) == input_length                                           
E       TypeError: object of type 'ICRS' has no len()                           
                                                                                
test_arrays.py:291: TypeError                                                   
===================== 1 failed, 12 passed in 2.49 seconds ======================

# Good, failure, as expected.
#
# what did we do? Remind me, git
#
$ git status                     
On branch i-1761                                                                
Your branch is up-to-date with 'mwcraig/i-1761'.                                
                                                                                
Changes not staged for commit:                                                  
  (use "git add <file>..." to update what will be committed)                    
  (use "git checkout -- <file>..." to discard changes in working directory)     
                                                                                
        modified:   test_arrays.py                                              
                                                                                
no changes added to commit (use "git add" and/or "git commit -a")       

#
# How about more detail, git?
$ git diff
diff --git a/astropy/coordinates/tests/test_arrays.py b/astropy/coordinates/test
index 2785b59..7eecfbb 100644                                                   
--- a/astropy/coordinates/tests/test_arrays.py                                  
+++ b/astropy/coordinates/tests/test_arrays.py                                  
@@ -278,3 +278,14 @@ def test_array_indexing():                                 
     assert c2.equinox == c1.equinox                                            
     assert c3.equinox == c1.equinox                                            
     assert c4.equinox == c1.equinox                                            
+                                                                               
+def test_array_len():                                                          
+    from .. import ICRS                                                        
+                                                                               
+    input_length = 5                                                           
+    ra = np.linspace(0, 360, input_length)                                     
+    dec = np.linspace(0, 90, input_length)                                     
+                                                                               
+    c = ICRS(ra, dec, unit=(u.degree, u.degree))                               
+                                                                               
+    assert len(c) == input_length 

#
# Stage, then commit
$ git add test_arrays.py          # get ready to commit
# always include a message with your commits
$ git commit -m'Add test for array coordinate length (issue #1761)'
[i-1761 23ba4ce] Add test for array coordinate length (issue #1761)             
 1 file changed, 11 insertions(+)   

#
# Now fix the bug...edit ../coordsystems.py
#
# Code fix made in editor
#
# Test change
#
$ py.test test_arrays.py   # Did our fix actually fix the problem?
============================= test session starts ==============================
platform darwin -- Python 2.7.6 -- pytest-2.4.2                                 
                                                                                
Running tests with Astropy version 0.4.dev6874.                                 
Running tests in test_arrays.py.                                                
                                                                                
Platform: Darwin-13.0.0-x86_64-i386-64bit                                       
                                                                                
Executable: /Users/mcraig/anaconda/envs/apy-1761/bin/python                     
                                                                                
Full Python Version:                                                            
2.7.6 |Anaconda 1.8.0 (x86_64)| (default, Nov 11 2013, 10:49:09)                
[GCC 4.0.1 (Apple Inc. build 5493)]                                             
                                                                                
encodings: sys: ascii, locale: US-ASCII, filesystem: utf-8, unicode bits: 15    
byteorder: little                                                               
float info: dig: 15, mant_dig: 15                                               
                                                                                
Numpy: 1.7.1                                                                    
Scipy: 0.13.0                                                                   
Matplotlib: 1.3.1                                                               
h5py: 2.2.0                                                                     
                                                                                
plugins: capturelog                                                             
collected 13 items                                                              
                                                                                
test_arrays.py .............                                                    
                                                                                
========================== 13 passed in 2.60 seconds ===========================

# Great! We fixed it!
#
# Run all of the coordinate tests to make sure we broke nothing...
$ py.test                        # do all of the coordinate tests pass?
============================= test session starts ==============================
platform darwin -- Python 2.7.6 -- pytest-2.4.2

Running tests with Astropy version 0.4.dev6877.
Running tests in /Users/mcraig/Development/astronomy/apy-bugs/astropy/coordinates/tests.

Platform: Darwin-13.0.0-x86_64-i386-64bit

Executable: /Users/mcraig/anaconda/envs/apy-1761/bin/python

Full Python Version:
2.7.6 |Anaconda 1.8.0 (x86_64)| (default, Nov 11 2013, 10:49:09)
[GCC 4.0.1 (Apple Inc. build 5493)]

encodings: sys: ascii, locale: US-ASCII, filesystem: utf-8, unicode bits: 15
byteorder: little
float info: dig: 15, mant_dig: 15

Numpy: 1.7.1
Scipy: 0.13.0
Matplotlib: 1.3.1
h5py: 2.2.0

plugins: capturelog
collected 103 items / 2 skipped

test_angles.py .......................
test_angular_separation.py ..
test_api.py ............
test_arrays.py .............
test_distance.py ........
test_formatting.py ....................
test_matching.py ....
test_name_resolve.py ss
test_transformations.py ...............
accuracy/test_fk4_no_e_fk4.py .
accuracy/test_fk4_no_e_fk5.py .
accuracy/test_galactic_fk4.py .
accuracy/test_icrs_fk5.py .

=================== 101 passed, 4 skipped in 152.23 seconds ===================

# So far, so good, now check ALL tests
$ cd ../../..                       # up to the top level to check all of the tests
$ python setup.py test           # grab a coffee now; this takes a while
onfigured with: --prefix=/Applications/Xcode.app/Contents/Developer/usr --with-gxx-include-dir=/usr/include/c++/4.2.1
Freezing version number to astropy/version.py
running test
running build
running build_py
copying astropy/version.py -> build/lib.macosx-10.5-x86_64-2.7/astropy
copying astropy/coordinates/tests/test_arrays.py -> build/lib.macosx-10.5-x86_64-2.7/astropy/coordinates/tests
...
OUTPUT TRUNCATED

#
# Success! Commit our change
#
$ git status                     # check what we have changed, if you forgot
On branch i-1761                                                                
Your branch is ahead of 'mwcraig/i-1761' by 1 commit.                           
  (use "git push" to publish your local commits)                                
                                                                                
Changes not staged for commit:                                                  
  (use "git add <file>..." to update what will be committed)                    
  (use "git checkout -- <file>..." to discard changes in working directory)     
                                                                                
        modified:   astropy/coordinates/coordsystems.py                         
                                                                                
no changes added to commit (use "git add" and/or "git commit -a") 

# stage change
$ git add astropy/coordinates/coordsystems.py # Copy and paste path from git status
                                             # unless you like typing
# commit, with message!
$ git commit -m"
    > Add len() to coordinates
    >
    > Closes #1761"
[i-1761 9688f7e] Add len() to coordinates                                       
 1 file changed, 2 insertions(+)  

# copy changes to the i-1761 branch in my github account...
$ git push           # or not--can wait until you are completely done
Counting objects: 52, done.                                                     
Delta compression using up to 2 threads.                                        
Compressing objects: 100% (11/11), done.                                        
Writing objects: 100% (11/11), 1.01 KiB | 0 bytes/s, done.                      
Total 11 (delta 9), reused 0 (delta 0)                                          
To https://github.com/mwcraig/astropy.git                                       
   7fa981e..9688f7e  i-1761 -> i-1761  

#
# Oops! Forgot to include some tests!
#
$ cd astropy/coordinates/tests   # back to add more tests
#
# edit test_arrays.py to add tests to test_array_len
#
# Now re-test
#
$ py.test test_arrays.py  # do these tests pass with the new tests?
============================= test session starts ==============================
platform darwin -- Python 2.7.6 -- pytest-2.4.2                                 
                                                                                
Running tests with Astropy version 0.4.dev6874.                                 
Running tests in test_arrays.py.                                                
                                                                                
Platform: Darwin-13.0.0-x86_64-i386-64bit                                       
                                                                                
Executable: /Users/mcraig/anaconda/envs/apy-1761/bin/python                     
                                                                                
Full Python Version:                                                            
2.7.6 |Anaconda 1.8.0 (x86_64)| (default, Nov 11 2013, 10:49:09)                
[GCC 4.0.1 (Apple Inc. build 5493)]                                             
                                                                                
encodings: sys: ascii, locale: US-ASCII, filesystem: utf-8, unicode bits: 15    
byteorder: little                                                               
float info: dig: 15, mant_dig: 15                                               
                                                                                
Numpy: 1.7.1                                                                    
Scipy: 0.13.0                                                                   
Matplotlib: 1.3.1                                                               
h5py: 2.2.0                                                                     
                                                                                
plugins: capturelog                                                             
collected 13 items                                                              
                                                                                
test_arrays.py .............                                                    
                                                                                
========================== 13 passed in 2.60 seconds ===========================

#
#   yes! now move up a couple of levels to check all tests
#
#   Do I really have to re-test everything? It certainly never hurts and
#   potentially saves time for the maintainers if you have introduced a bug.
#
$ cd ../../..          # back to the top level
$ python setup.py test # grab another coffee....
OUTPUT OMITTED

#
# Success! Time to commit
#
$ git status           # optional, but useful for figuring out what to add
On branch i-1761                                                                
Your branch is up-to-date with 'mwcraig/i-1761'.                                
                                                                                
Changes not staged for commit:                                                  
  (use "git add <file>..." to update what will be committed)                    
  (use "git checkout -- <file>..." to discard changes in working directory)     
                                                                                
        modified:   astropy/coordinates/tests/test_arrays.py                                              
                                                                                
no changes added to commit (use "git add" and/or "git commit -a")

# 
# Stage...
$ git add astropy/coordinates/tests/test_arrays.py # use path git status supplies

# Commit, with a message...
$ git commit -m"Add tests of len() for scalar coordinate and length 1 coordinate"
[i-1761 ed92c2d] Add tests of len() for scalar coordinate and length 1 coordinat
e                                                                               
 1 file changed, 10 insertions(+), 5 deletions(-) 

# Almost done...
#
# Edit the changelog in CHANGES.rst
#
# Then commit changes
#
$ git add CHANGES.rst
$ git commit -m"Add changelog entry for 1761"
[i-1761 bf6fcfc] Add changelog entry for 1761                                   
 1 file changed, 2 insertions(+)                                                

#
# Push changes to my GitHub account
#
$ git push   # make sure all of our changes are in our github account
Counting objects: 47, done.                                                     
Delta compression using up to 2 threads.                                        
Compressing objects: 100% (8/8), done.                                          
Writing objects: 100% (8/8), 836 bytes | 0 bytes/s, done.                       
Total 8 (delta 6), reused 0 (delta 0)                                           
To https://github.com/mwcraig/astropy.git                                       
   9688f7e..bf6fcfc  i-1761 -> i-1761 

#
# Go to GitHub to make pull request. Remember to switch to this branch before 
# pushing the pull request button.
#