"""Shared pytest setup.

The validation tests import `validation.checks`, which lives at the repository root rather than
under python/, so the root has to be importable. Putting it here rather than in the test module
keeps the import statement in the test file honest about what it depends on.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
