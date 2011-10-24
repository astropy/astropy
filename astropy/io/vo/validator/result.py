# Copyright (C) 2008-2010 Association of Universities for Research in Astronomy (AURA)

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

#     1. Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.

#     2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.

#     3. The name of AURA and its representatives may not be used to
#       endorse or promote products derived from this software without
#       specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.

from __future__ import absolute_import

# STDLIB
import hashlib
import httplib
import os
import cPickle as pickle
import shutil
import subprocess
import sys
import urllib2
import warnings

# VO
from .. import table
from .. import voexceptions
from .. import xmlutil
from ..util import IS_PY3K

class Result:
    def __init__(self, url, root='results'):
        self.url = url
        m = hashlib.md5()
        m.update(url)
        self._hash = m.hexdigest()
        self._root = root
        self._path = os.path.join(
            self._hash[0:2], self._hash[2:4], self._hash[4:])
        if not os.path.exists(self.get_dirpath()):
            os.makedirs(self.get_dirpath())
        self.load_attributes()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.save_attributes()

    def get_dirpath(self):
        return os.path.join(self._root, self._path)

    def get_htmlpath(self):
        return self._path

    def get_attribute_path(self):
        return os.path.join(self.get_dirpath(), "values.dat")

    def get_vo_xml_path(self):
        return os.path.join(self.get_dirpath(), "vo.xml")

    # ATTRIBUTES

    def load_attributes(self):
        path = self.get_attribute_path()
        if os.path.exists(path):
            try:
                with open(path, 'rb') as fd:
                    self._attributes = pickle.load(fd)
            except:
                shutil.rmtree(self.get_dirpath())
                os.makedirs(self.get_dirpath())
                self._attributes = {}
        else:
            self._attributes = {}

    def save_attributes(self):
        path = self.get_attribute_path()
        with open(path, 'wb') as fd:
            pickle.dump(self._attributes, fd)

    def __getitem__(self, key):
        return self._attributes[key]

    def __setitem__(self, key, val):
        self._attributes[key] = val

    def __contains__(self, key):
        return key in self._attributes

    # VO XML

    def download_xml_content(self):
        path = self.get_vo_xml_path()

        if 'network_error' not in self._attributes:
            self['network_error'] = None

        if os.path.exists(path):
            return

        def fail(reason):
            reason = str(reason)
            with open(path, 'wb') as fd:
                fd.write("FAILED: %s\n" % reason)
            self['network_error'] = reason

        r = None
        try:
            if IS_PY3K:
                r = urllib2.urlopen(self.url.decode('ascii'))
            else:
                r = urllib2.urlopen(self.url, timeout=10)
        except urllib2.URLError as e:
            if hasattr(e, 'reason'):
                reason = e.reason
            else:
                reason = e.code
            fail(reason)
            return
        except httplib.HTTPException as e:
            fail("HTTPException: %s" % str(e))
            return

        if r is None:
            fail("Invalid URL")
            return

        content = r.read()
        r.close()

        with open(path, 'wb') as fd:
            fd.write(content)

    def get_xml_content(self):
        path = self.get_vo_xml_path()
        if not os.path.exists(path):
            self.download_xml_content()
        with open(path, 'rb') as fd:
            content = fd.read()
        return content

    def validate_vo(self):
        orig_stdout = sys.stdout
        path = self.get_vo_xml_path()
        if not os.path.exists(path):
            self.download_xml_content()
        self['version'] = ''
        if self['network_error'] is not None:
            self['nwarnings'] = 0
            self['nexceptions'] = 0
            self['warnings'] = []
            self['xmllint'] = None
            self['warning_types'] = set()
            return

        nexceptions = 0
        nwarnings = 0
        t = None
        lines = []
        with open(path, 'rb') as input:
            with warnings.catch_warnings(record=True) as warning_lines:
                try:
                    t = table.parse(input, pedantic=False, filename=path)
                except ValueError as e:
                    lines.append(str(e))
                    nexceptions += 1
        lines = [str(x.message) for x in warning_lines] + lines

        if t is not None:
            self['version'] = version = t.version
        else:
            self['version'] = version = "1.0"

        if 'xmllint' not in self:
            # Now check the VO schema based on the version in
            # the file.
            success, stdout, stderr = xmlutil.validate_schema(path, version)
            self['xmllint'] = (success == 0)
            self['xmllint_content'] = stderr

        warning_types = set()
        for line in lines:
            w = voexceptions.parse_vowarning(line)
            if w['is_warning']:
                nwarnings += 1
            if w['is_exception']:
                nexceptions += 1
            warning_types.add(w['warning'])

        self['nwarnings'] = nwarnings
        self['nexceptions'] = nexceptions
        self['warnings'] = lines
        self['warning_types'] = warning_types

    def has_warning(self, warning_code):
        return warning_code in self['warning_types']

    def match_expectations(self):
        if self['expected'] == 'good':
            return (not self['network_error'] and
                    self['nwarnings'] == 0 and
                    self['nexceptions'] == 0)
        elif self['expected'] == 'incorrect':
            return (not self['network_error'] and
                    (self['nwarnings'] > 0 or
                     self['nexceptions'] > 0))
        elif self['expected'] == 'broken':
            return self['network_error'] is not None

    def validate_with_votlint(self, path_to_stilts_jar):
        filename = self.get_vo_xml_path()
        p = subprocess.Popen(
            "java -jar %s votlint validate=false %s" %
            (path_to_stilts_jar, filename),
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        if len(stdout) or p.returncode:
            self['votlint'] = False
        else:
            self['votlint'] = True
        self['votlint_content'] = stdout

def get_result_subsets(results, root):
    all_results      = []
    not_expected     = []
    fail_schema      = []
    schema_mismatch  = []
    fail_votlint     = []
    votlint_mismatch = []
    network_failures = []
    version_10       = []
    version_11       = []
    version_12       = []
    version_unknown  = []
    has_warnings     = []
    warning_set      = {}
    has_exceptions   = []
    exception_set    = {}

    for url in results:
        x = Result(url, root=root)
        all_results.append(x)
        if not x.match_expectations():
            not_expected.append(x)
        if x['xmllint'] is False:
            fail_schema.append(x)
        if x['xmllint'] is False and x['nwarnings'] == 0 and x['nexceptions'] == 0:
            schema_mismatch.append(x)
        if 'votlint' in x and x['votlint'] is False:
            fail_votlint.append(x)
            if x['nwarnings'] == 0 and x['nexceptions'] == 0 and x['network_error'] is None:
                votlint_mismatch.append(x)
        if x['network_error'] is not None:
            network_failures.append(x)
        version = x['version']
        if version == '1.0':
            version_10.append(x)
        elif version == '1.1':
            version_11.append(x)
        elif version == '1.2':
            version_12.append(x)
        else:
            version_unknown.append(x)
        if x['nwarnings'] > 0:
            has_warnings.append(x)
            for warning in x['warning_types']:
                if warning is not None and len(warning) == 3 and warning.startswith('W'):
                    warning_set.setdefault(warning, [])
                    warning_set[warning].append(x)
        if x['nexceptions'] > 0:
            has_exceptions.append(x)
            for exc in x['warning_types']:
                if exc is not None and len(exc) == 3 and exc.startswith('E'):
                    exception_set.setdefault(exc, [])
                    exception_set[exc].append(x)

    warning_set = list(warning_set.items())
    warning_set.sort()
    exception_set = list(exception_set.items())
    exception_set.sort()

    tables = [
        ('all', u'All tests', all_results),
        ('unexpected', u'Unexpected', not_expected),
        ('schema', u'Invalid against schema', fail_schema),
        ('schema_mismatch', u'Invalid against schema/Passed vo.table', schema_mismatch, ['ul']),
        ('fail_votlint', u'Failed votlint', fail_votlint),
        ('votlint_mismatch', u'Failed votlint/Passed vo.table', votlint_mismatch, ['ul']),
        ('network_failures', u'Network failures', network_failures),
        ('version1.0', 'Version 1.0', version_10),
        ('version1.1', 'Version 1.1', version_11),
        ('version1.2', 'Version 1.2', version_12),
        ('version_unknown', 'Version unknown', version_unknown),
        ('warnings', 'Warnings', has_warnings)]
    for warning_code, warnings in warning_set:
        warning_class = getattr(voexceptions, warning_code, None)
        if warning_class:
            warning_descr = warning_class.get_short_name()
            tables.append(
                (warning_code, '%s: %s' % (warning_code, warning_descr), warnings, ['ul', 'li']))
    tables.append(
        ('exceptions', 'Exceptions', has_exceptions))
    for exception_code, exceptions in exception_set:
        exception_class = getattr(voexceptions, exception_code, None)
        if exception_class:
            exception_descr = exception_class.get_short_name()
            tables.append(
                (exception_code, '%s: %s' % (exception_code, exception_descr), exceptions, ['ul', 'li']))

    return tables
