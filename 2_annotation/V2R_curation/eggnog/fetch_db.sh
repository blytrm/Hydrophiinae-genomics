#!/usr/bin/env bash
# download_eggnog_data.py wrote zero-byte files here, so fetch what emapper
# needs directly. curl -C - resumes a partial file; each one is checked against
# the size the server reports before anything is unpacked.
set -uo pipefail
DATA=/scratchdata1/users/a1864358/sanders_lab/prog/eggnog_db
URL=http://eggnog5.embl.de/download/emapperdb-5.0.2
mkdir -p "$DATA"
cd "$DATA"

fetch() {
    f=$1
    want=$(curl -sfI "$URL/$f" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2+0}')
    attempt=1
    while [ "$attempt" -le 8 ]; do
        have=0
        [ -f "$f" ] && have=$(stat -c%s "$f")
        if [ "$have" -eq "$want" ]; then
            echo "OK $f ($want bytes)"
            return 0
        fi
        echo "$f: $have / $want bytes -- attempt $attempt"
        curl -sS -C - -o "$f" "$URL/$f"
        attempt=$((attempt + 1))
    done
    echo "FAILED $f" >&2
    return 1
}

fetch eggnog.db.gz            || exit 1
fetch eggnog_proteins.dmnd.gz || exit 1
fetch eggnog.taxa.tar.gz      || exit 1

[ -s eggnog.db ]            || { echo "unpacking eggnog.db";  gunzip -c eggnog.db.gz > eggnog.db; }
[ -s eggnog_proteins.dmnd ] || { echo "unpacking dmnd";       gunzip -c eggnog_proteins.dmnd.gz > eggnog_proteins.dmnd; }
[ -s eggnog.taxa.db ]       || { echo "unpacking taxa";       tar xzf eggnog.taxa.tar.gz; }
echo "EGGNOG DB READY"
ls -la "$DATA"
