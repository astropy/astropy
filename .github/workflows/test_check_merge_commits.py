#!/usr/bin/env python3
"""
Unit tests for check_merge_commits.py script.

These tests verify that the merge commit detection logic works correctly.
"""
import unittest
from unittest.mock import patch, MagicMock
import sys
import os

# Add the script directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import functions from the script
from check_merge_commits import (
    get_commit_details,
    is_merge_from_main,
    get_merge_commits
)


class TestMergeCommitDetection(unittest.TestCase):
    """Test cases for merge commit detection functionality."""

    def test_is_merge_from_main_with_main_branch(self):
        """Test detection of merge from main branch."""
        commit_details = {
            'hash': 'abc123',
            'subject': "Merge branch 'main' into feature",
            'parents': ['parent1', 'parent2']
        }
        self.assertTrue(is_merge_from_main(commit_details))

    def test_is_merge_from_main_with_master_branch(self):
        """Test detection of merge from master branch."""
        commit_details = {
            'hash': 'abc123',
            'subject': "Merge branch 'master' into feature",
            'parents': ['parent1', 'parent2']
        }
        self.assertTrue(is_merge_from_main(commit_details))

    def test_is_merge_from_main_with_origin_main(self):
        """Test detection of merge from origin/main."""
        commit_details = {
            'hash': 'abc123',
            'subject': "Merge remote-tracking branch 'origin/main' into feature",
            'parents': ['parent1', 'parent2']
        }
        self.assertTrue(is_merge_from_main(commit_details))

    def test_is_merge_from_main_with_upstream_main(self):
        """Test detection of merge from upstream/main."""
        commit_details = {
            'hash': 'abc123',
            'subject': "Merge upstream/main into feature",
            'parents': ['parent1', 'parent2']
        }
        self.assertTrue(is_merge_from_main(commit_details))

    def test_is_merge_from_main_with_feature_branch(self):
        """Test that feature branch merges are not flagged."""
        commit_details = {
            'hash': 'abc123',
            'subject': "Merge branch 'feature-xyz' into another-feature",
            'parents': ['parent1', 'parent2']
        }
        self.assertFalse(is_merge_from_main(commit_details))

    def test_is_merge_from_main_with_case_insensitivity(self):
        """Test that detection is case-insensitive."""
        commit_details = {
            'hash': 'abc123',
            'subject': "Merge Branch 'Main' Into Feature",
            'parents': ['parent1', 'parent2']
        }
        self.assertTrue(is_merge_from_main(commit_details))

    def test_is_merge_from_main_with_regular_commit(self):
        """Test that regular commits are not flagged."""
        commit_details = {
            'hash': 'abc123',
            'subject': "Add new feature for data processing",
            'parents': ['parent1']
        }
        self.assertFalse(is_merge_from_main(commit_details))

    @patch('check_merge_commits.run_git_command')
    def test_get_merge_commits_no_merges(self, mock_git):
        """Test get_merge_commits when there are no merge commits."""
        mock_git.return_value = ""
        result = get_merge_commits()
        self.assertEqual(result, [])

    @patch('check_merge_commits.run_git_command')
    def test_get_merge_commits_with_merges(self, mock_git):
        """Test get_merge_commits when merge commits are present."""
        mock_git.return_value = "abc123 Merge branch 'main'\ndef456 Another merge commit"
        result = get_merge_commits()
        self.assertEqual(result, ['abc123', 'def456'])

    @patch('check_merge_commits.run_git_command')
    def test_get_commit_details(self, mock_git):
        """Test get_commit_details retrieves correct information."""
        # First call returns subject, second call returns parents
        mock_git.side_effect = [
            "Merge branch 'main' into feature",
            "parent1hash parent2hash"
        ]

        result = get_commit_details('abc123')

        self.assertEqual(result['hash'], 'abc123')
        self.assertEqual(result['subject'], "Merge branch 'main' into feature")
        self.assertEqual(result['parents'], ['parent1hash', 'parent2hash'])


class TestMergeCommitPatterns(unittest.TestCase):
    """Test various merge commit message patterns."""

    def test_pattern_with_quotes(self):
        """Test merge with quoted branch name."""
        patterns = [
            "Merge branch 'main' into feature",
            'Merge branch "main" into feature',
            "Merge branch main into feature"
        ]
        for pattern in patterns:
            commit = {'hash': 'x', 'subject': pattern, 'parents': ['a', 'b']}
            self.assertTrue(is_merge_from_main(commit),
                          f"Failed to detect: {pattern}")

    def test_pattern_with_remotes(self):
        """Test merge with remote references."""
        patterns = [
            "Merge remote-tracking branch 'origin/main'",
            "Merge branch 'main' of github.com:astropy/astropy",
            "Merge upstream/main",
            "Merge origin/master into HEAD"
        ]
        for pattern in patterns:
            commit = {'hash': 'x', 'subject': pattern, 'parents': ['a', 'b']}
            self.assertTrue(is_merge_from_main(commit),
                          f"Failed to detect: {pattern}")

    def test_pattern_non_main_branches(self):
        """Test that non-main branch merges are not flagged."""
        patterns = [
            "Merge branch 'feature-x' into develop",
            "Merge pull request #123 from user/branch",
            "Merge branch 'bugfix-456'",
            "Merge maintenance branch"
        ]
        for pattern in patterns:
            commit = {'hash': 'x', 'subject': pattern, 'parents': ['a', 'b']}
            self.assertFalse(is_merge_from_main(commit),
                           f"False positive for: {pattern}")


if __name__ == '__main__':
    unittest.main()
