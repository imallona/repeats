"""
Shared pytest configuration and path setup.

Adds workflow/scripts to sys.path so unit tests can import scripts directly.
"""
import sys
import os

SCRIPTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', 'workflow', 'scripts')
)
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)
