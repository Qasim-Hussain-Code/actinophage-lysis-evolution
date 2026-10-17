#!/usr/bin/env python3
"""Write an isolated, outcome-blind component and matching feasibility audit."""
from pathlib import Path
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from actinolysis.amendment_feasibility import main

if __name__ == "__main__":
    main()
