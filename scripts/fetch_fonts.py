"""Download the fonts required by the TU Darmstadt Typst template.

The template (``@preview/athena-tu-darmstadt-thesis``) uses Roboto (headings,
captions) and XCharter (body text).  Both are free fonts distributed on CTAN;
this script downloads them and extracts only the needed OpenType files into
``thesis/assets/fonts`` (git-ignored).  Cross-platform replacement for the
template's ``download_fonts.sh``.

    python scripts/fetch_fonts.py
"""

from __future__ import annotations

import io
import sys
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FONT_DIR = ROOT / "thesis" / "assets" / "fonts"

PACKAGES = {
    # target folder: (url, predicate selecting the files to extract)
    "roboto": (
        "https://mirrors.ctan.org/fonts/roboto.zip",
        lambda n: "/opentype/" in n and n.endswith(".otf")
        and (Path(n).name.startswith("Roboto-") or Path(n).name.startswith("RobotoSlab-")),
    ),
    "xcharter": (
        "https://mirrors.ctan.org/fonts/xcharter.zip",
        lambda n: "/opentype/" in n and n.endswith(".otf") and Path(n).name.startswith("XCharter-"),
    ),
}


def fetch(name: str, url: str, select) -> int:
    target = FONT_DIR / name
    if target.exists() and any(target.glob("*.otf")):
        print(f"{name}: already present in {target}")
        return 0
    print(f"{name}: downloading {url} ...", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "nbody-thesis-font-fetcher"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        data = resp.read()
    target.mkdir(parents=True, exist_ok=True)
    count = 0
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        for member in zf.namelist():
            if select(member):
                (target / Path(member).name).write_bytes(zf.read(member))
                count += 1
    print(f"{name}: extracted {count} font files to {target}")
    return count


def main() -> int:
    for name, (url, select) in PACKAGES.items():
        fetch(name, url, select)
    return 0


if __name__ == "__main__":
    sys.exit(main())
