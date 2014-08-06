#!/usr/bin/env python

# I wasn't happy with any of the GitHub libraries for Python that I tried so I
# just used the GitHub API directly.  If someone would like to rewrite this
# using a library please be my guest

from __future__ import unicode_literals

import argparse
import base64
import getpass
import io
import json
import logging
import os
import re
import stat
import sys

try:
    from urllib.request import Request, urlopen
    from urllib.error import HTTPError
    from urllib.parse import urlencode
except ImportError:
    from urllib2 import Request, urlopen, HTTPError
    from urllib import urlencode

try:
    input = raw_input
except NameError:
    pass

# Because pkg_resources provides better version parsing than distutils
import pkg_resources


BASE_URL = 'https://api.github.com/repos/'

# This regex ensures that only the 'Conflicts:' section at the end of the
# commit message is matched (in case there are multiple 'Conflicts' sections
# which can happen).
CONFLICTS_RE = re.compile(r'((?:.|\n)+)\nConflicts:(\n.*)+', flags=re.M)

log = logging.getLogger()



class _MaxLevelFilter(logging.Filter):
    def __init__(self, maxlevel):
        self.maxlevel = maxlevel

    def filter(self, record):
        return record.levelno <= self.maxlevel


class GithubRequestError(Exception):
    pass


class GithubSuggestBackports(object):
    # Cache all the commits found for the given branch so we don't have to
    # re-request them for each pull request
    _cached_commits = []

    def __init__(self, owner, repo, branch, username=None, password=None):
        self.owner = owner
        self.repo = repo
        self.branch = branch
        if username is not None and password is not None:
            # We can't rely on urllib2 to handle basic authentication in the
            # normal way since GitHub requests don't always have
            # www-authenticate in the headers
            auth = ':'.join((username, password)).encode('ascii')
            self._auth = base64.b64encode(auth).decode('ascii')
        else:
            self._auth = None

    def _github_repo_request(self, *resource, **parameters):
        resource = tuple(str(r) for r in resource)
        url = BASE_URL + '/'.join((self.owner, self.repo) + resource)
        if parameters:
            url += '?' + urlencode(parameters)
        log.debug('Requesting ' + url)
        req = Request(url)
        if self._auth:
            req.add_header('Authorization', 'Basic ' + self._auth)
        try:
            f = urlopen(req)
            enc = f.headers.get_content_charset()
            content = f.read().decode(enc)
            response = json.loads(content)
        except HTTPError as e:
            response = json.loads(e.fp.read().decode('utf8'))
            if 'message' in response:
                raise GithubRequestError(response['message'])
            raise e
        return response

    def get_tags(self):
        return self._github_repo_request('tags')

    def get_milestones(self, state=None):
        parameters = {}
        if state is not None:
            parameters['state'] = state

        return self._github_repo_request('milestones', **parameters)

    def iter_issues(self, milestone=None, state=None):
        parameters = {}
        if milestone is not None:
            parameters['milestone'] = milestone
        if state is not None:
            parameters['state'] = state

        parameters['page'] = 1
        issues = []
        while True:
            if not issues:
                response = self._github_repo_request('issues', **parameters)
                if response:
                    issues.extend(response)
                    parameters['page'] += 1
                else:
                    raise StopIteration
            yield issues.pop(0)

    def iter_issue_events(self, issue, filter_=None, count=None):
        events = []
        page = 1
        while (count is None or count):
            # Events can be paginated
            if not events:
                next = self._github_repo_request('issues', issue, 'events',
                                                 page=page)
                if not next:
                    raise StopIteration
                if filter_ is not None:
                    next = filter(lambda e: e['event'] == filter_, next)
                events.extend(next)
                page += 1
                # Either continue to the next page of events or start popping
                # any found events off the first page
                continue
            yield events.pop(0)
            count -= 1

    def get_pull_request_merge_commit(self, pr):
        """Returns the full commit object of the merge commit for a pull
        request or `None` if the given PR has not been merged.

        This is different from the commit named by merge_commit_sha listed in a
        pull request in that it's the commit that actually goes into mainline
        branch. The commit listed in merge_commit_sha only seems to be an
        artifact of how GitHub implements pull requests.
        """

        events = list(self.iter_issue_events(pr, filter_='merged', count=1))
        if events:
            return self.get_commit(events[0]['commit_id'])

    def get_commits(self, sha):
        """Get the first page of commits in the tree starting at sha.

        Commits are returned 30 at a time and paginated according to sha. So in
        order to get the second page of commits it's necessary to use a
        subsequent call to get_commits using the sha of the last commit from
        the previous call (which will be the first commit listed in the second
        call).
        """

        return self._github_repo_request('commits', sha=sha)

    def get_commit(self, sha):
        """Return a single commit."""

        return self._github_repo_request('commits', sha)

    def iter_pull_requests(self, state=None):
        parameters = {}
        if state is not None:
            parameters['state'] = state

        parameters['page'] = 1
        prs = []
        while True:
            if not prs:
                response = self._github_repo_request('pulls', **parameters)
                if response:
                    prs.extend(response)
                    parameters['page'] += 1
                else:
                    raise StopIteration
            yield prs.pop(0)

    def get_pull_request(self, number):
        try:
            pr = self._github_repo_request('pulls', str(number))
        except GithubRequestError as e:
            if e.message == 'Not Found':
                return None
            raise
        return pr

    def find_unmerged_commit(self, commit, since=None):
        def expand_cache():
            if not self._cached_commits:
                # Initialize with the first page of commits from the bug fix
                # branch
                next_commits = self.get_commits(self.branch)
            else:
                last_commit = self._cached_commits[-1]
                if last_commit['commit']['committer']['date'] <= since:
                    return False
                next_commits = self.get_commits(last_commit['sha'])[1:]
            if next_commits:
                self._cached_commits.extend(next_commits)
                return True
            else:
                return False

        idx = 0

        while True:
            try:
                merged_commit = self._cached_commits[idx]
            except IndexError:
                # Try growing the list of commits; but if there are no more to be
                # found return None
                if expand_cache():
                    continue
                return None

            # For cherry-picks we can't rely on comparing the sha, but the
            # author and commit message should be close enough
            a = commit['commit']
            b = merged_commit['commit']

            # Remove conflicts from the cherry-picked commit's commit message;
            # conflicts can cause the message to be different where it
            # otherwise wouldn't have been, and we don't care if there were
            # conflicts so long as it was merged successfully
            b['message'] = CONFLICTS_RE.sub(r'\1', b['message'])

            if a['author'] == b['author'] and a['message'] == b['message']:
                return merged_commit

            idx += 1

    def get_next_milestone(self):
        """Get the next open milestone that has the same version prefix as the
        branch.  For example if the repo has milestones v0.2.1 and v0.2.2 and the
        branch is v0.2.x, this will return v0.2.1.
        """

        prefix = self.branch[:-1]
        milestones = [m for m in self.get_milestones(state='open')
                      if m['title'].startswith(prefix)]
        sort_key = lambda m: int(m['title'].rsplit('.', 1)[1])
        return sorted(milestones, key=sort_key)[0]

    _last_tag = None
    def get_last_tag(self):
        if self._last_tag is not None:
            return self._last_tag
        branch_ver = pkg_resources.parse_version(self.branch.lstrip('v'))
        tags = sorted(self.get_tags(),
                      key=lambda t: pkg_resources.parse_version(t['name']),
                      reverse=True)
        # Get the last tag that should be in this branch
        for tag in tags:
            tag_ver = pkg_resources.parse_version(tag['name'].lstrip('v'))
            if tag_ver[:2] == branch_ver[:2]:
                self._last_tag = tag
                return tag

        self._last_tag = False

    _last_tag_commit = None
    def get_last_tag_commit(self):
        if self._last_tag_commit is not None:
            return self._last_tag_commit
        last_tag = self.get_last_tag()
        if last_tag:
            last_tag_commit = self.get_commit(last_tag['commit']['sha'])
        else:
            last_tag_commit = False

        self._last_tag_commit = last_tag_commit
        return last_tag_commit


    def iter_suggested_prs(self):
        next_milestone = self.get_next_milestone()
        next_ms_num = next_milestone['number']
        log.info("Finding PRs in milestone {0} that haven't been merged into "
                 "{1}".format(next_milestone['title'], self.branch))
        log.info('Merge these into {0} by doing "git checkout {0}; git pull; '
                 'git cherry-pick -m 1 <SHA>"'.format(self.branch))
        last_tag_commit = self.get_last_tag_commit()
        last_tag_date = last_tag_commit['commit']['committer']['date']

        # Get the issue #s of all closed issues in the relevant milestone
        milestone_issues = set(issue['number'] for issue in
                               self.iter_issues(milestone=next_ms_num,
                                                state='closed'))

        # Now get all PRs and filter by whether or not they belong to the
        # milestone; requesting them all at once is still faster than
        # requesting one at a time. This would also be easier if the API
        # supported sorting on PR lists
        for pr in self.iter_pull_requests(state='closed'):
            if (pr['number'] not in milestone_issues or not pr['merged_at']):
                continue

            merge_commit = self.get_pull_request_merge_commit(pr['number'])

            if not self.find_unmerged_commit(merge_commit,
                                             since=last_tag_date):
                yield pr, merge_commit['sha']


def main(argv):
    parser = argparse.ArgumentParser(
        description='Find pull requests that need be backported to a bug fix '
                    'branch')
    parser.add_argument('owner', metavar='OWNER',
                        help='owner of the repository')
    parser.add_argument('repo', metavar='REPO', help='the repository name')
    parser.add_argument('branch', metavar='BRANCH',
                        help='the name of the bug fix branch (eg. v0.2.x)')
    parser.add_argument('-f', '--file', metavar='FILE',
                        help='save the cherry-pick script to a file; '
                             'otherwise it is written to stdout')
    parser.add_argument('--debug', action='store_true')

    args = parser.parse_args(argv)

    # Configure log
    log.setLevel(logging.DEBUG)
    stdout_handler = logging.StreamHandler(sys.stdout)
    if args.debug:
        stdout_handler.setLevel(logging.DEBUG)
    else:
        stdout_handler.setLevel(logging.INFO)
    stdout_handler.addFilter(_MaxLevelFilter(logging.INFO))
    log.addHandler(stdout_handler)
    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(logging.WARNING)
    log.addHandler(stderr_handler)


    log.info("Enter your GitHub username and password so that API requests "
             "aren't as severely rate-limited...")
    username = input('Username: ')
    password = getpass.getpass('Password: ')
    suggester = GithubSuggestBackports(args.owner, args.repo, args.branch,
                                       username, password)

    pr_format = '[#{0}][{1}]: {2}'
    suggestions = []
    for pr, sha in suggester.iter_suggested_prs():
        # If sys.stdout's default encoding has a limited codepage this blows up if
        # the PR title contains unencodable characters =_=
        title = pr['title'].encode('ascii',
                errors='replace').decode('ascii')
        log.info(pr_format.format(pr['number'], sha, title))
        suggestions.append((pr, sha))

    suggestions.sort(key=lambda p: p[0]['merged_at'])

    script_lines = [
        '#!/bin/bash',
        '# git commands:',
        'git checkout {0} || exit 1'.format(args.branch),
        'git pull upstream {0} || exit 1'.format(args.branch)
    ]
    for pr, sha in suggestions:
        script_lines.append('# ' + pr_format.format(pr['number'], sha,
                                                    pr['title']))
        script_lines.append('git cherry-pick -m 1 {0} || exit 1'.format(sha))

    if args.file:
        with io.open(args.file, 'w', encoding='utf8') as f:
            f.writelines(line + '\n' for line in script_lines)
        os.chmod(args.file, stat.S_IRWXU)
    else:
        for line in script_lines:
            log.info(line)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
