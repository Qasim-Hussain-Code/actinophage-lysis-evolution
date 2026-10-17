"""Save the exact installed WSL Python dependency set."""
from pathlib import Path
import subprocess
import sys
root=Path(__file__).resolve().parents[1]
result=subprocess.run([sys.executable,'-m','pip','freeze'],check=True,capture_output=True,text=True)
(root/'requirements-lock.txt').write_text(result.stdout.replace('\r\n','\n'),encoding='utf-8')
