"""Download one dated BUSCO lineage release and unpack it as <download_path>/lineages/<name>.

`busco --download` always fetches the newest release, so the dated tarball is
fetched directly and checked against its md5 (from busco-data's file_versions.tsv).

usage: get_busco_lineage.py NAME DATE MD5 DOWNLOAD_PATH
"""
import hashlib
import shutil
import sys
import tarfile
import urllib.request
from pathlib import Path

name, date, md5, download_path = sys.argv[1:]
url = f"https://busco-data.ezlab.org/v5/data/lineages/{name}.{date}.tar.gz"
lineages = Path(download_path) / "lineages"
lineages.mkdir(parents=True, exist_ok=True)
tarball = lineages / f"{name}.{date}.tar.gz"

with urllib.request.urlopen(url, timeout=120) as response, open(tarball, "wb") as out:
    shutil.copyfileobj(response, out)

digest = hashlib.md5(tarball.read_bytes()).hexdigest()
if digest != md5:
    sys.exit(f"{url}: md5 {digest}, expected {md5}")

with tarfile.open(tarball) as tar:
    tar.extractall(lineages, filter="data")
tarball.unlink()
