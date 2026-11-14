#!/usr/bin/env bash
set -euo pipefail

root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)
layer="$root/test/linecablemodels/prose"
lcm=/home/amartins/Documents/KUL/LineCableModels/worktrees/LineCableModels-datamodel/.github/prose
generated="$root/test/linecablemodels/.prose"
mkdir -p "$generated/styles/config/vocabularies/HarmonicPowerModels"
for style in "$lcm"/styles/*; do
    [[ $(basename -- "$style") == config ]] && continue
    ln -sfn -- "$style" "$generated/styles/$(basename -- "$style")"
done
for config in "$lcm"/styles/config/*; do
    [[ $(basename -- "$config") == vocabularies ]] && continue
    ln -sfn -- "$config" "$generated/styles/config/$(basename -- "$config")"
done
ln -sfn -- "$layer/accept.txt" "$generated/styles/config/vocabularies/HarmonicPowerModels/accept.txt"
# Preserve LCM's configuration. Change its styles path and add the vocabulary.
python3 - "$lcm/vale.ini" "$generated/vale.ini" "$generated/styles" <<'PY'
from pathlib import Path
import sys
source, target, styles = sys.argv[1:]
lines = Path(source).read_text().splitlines(keepends=True)
assert sum(line.startswith('StylesPath =') for line in lines) == 1
assert not any(line.startswith('Vocab =') for line in lines)
Path(target).write_text(''.join(
    f'StylesPath = {styles}\nVocab = HarmonicPowerModels\n'
    if line.startswith('StylesPath =') else line for line in lines))
PY
cd -- "$root"
# Run reports contain generated results. CSV files contain tables rather than prose.
# accept.txt configures Vale's vocabulary and contains no prose.
mapfile -d '' files < <(git ls-files --cached --others --exclude-standard -z -- \
    ext/HarmonicPowerModelsLineCableModelsExt test/linecablemodels examples/rhhc-variants \
    ':!:examples/rhhc-variants/output/**' ':!:*.csv' \
    ':!:test/linecablemodels/prose/fixture.md' \
    ':!:test/linecablemodels/prose/accept.txt')
set +e
vale --no-global --config="$generated/vale.ini" "${files[@]}"
territory_vale=$?
cspell --no-progress --config "$lcm/cspell.json" "${files[@]}"
territory_cspell=$?
set -e

# The complete fixture retains the citation, admitted compounds, and negative controls.
set +e
vale --no-global --config="$generated/vale.ini" --output=JSON "$layer/fixture.md" > "$generated/fixture.json"
fixture_status=$?
set -e
set +e
python3 - "$layer/fixture.md" "$generated/fixture.json" <<'PY'
import json
from pathlib import Path
import sys
fixture, report = sys.argv[1:]
lines = Path(fixture).read_text().splitlines()
start = lines.index('<!-- defined -->') + 1
stray = lines.index('<!-- stray -->') + 1
alerts = [alert for values in json.loads(Path(report).read_text()).values() for alert in values]
for alert in alerts:
    print(f"fixture.md:{alert['Line']}: {alert['Severity']} {alert['Check']}: {alert['Message']}")
for severity in ('error', 'warning', 'suggestion'):
    print(f"Fixture {severity}s: {sum(a['Severity'] == severity for a in alerts)}")
assert not [a for a in alerts if start < a['Line'] < stray], alerts
for row, line in enumerate(lines, 1):
    if row > stray and 'robust' in line:
        assert any(a['Line'] == row and 'robust' in a['Match'].lower() for a in alerts), (row, alerts)
print('Vale fixture: defined compounds passed; all four stray uses refused.')
PY
fixture_check=$?
cspell --no-progress --config "$lcm/cspell.json" "$layer/fixture.md"
fixture_cspell=$?
set -e
printf 'Territory: Vale exit %s, cspell exit %s\n' "$territory_vale" "$territory_cspell"
printf 'Fixture: Vale exit %s, assertions exit %s, cspell exit %s\n' \
    "$fixture_status" "$fixture_check" "$fixture_cspell"
(( territory_vale == 0 && territory_cspell == 0 && fixture_status <= 1 && \
   fixture_check == 0 && fixture_cspell == 0 ))
