import sys
from pathlib import Path
import pandas as pd
from openpyxl import Workbook, load_workbook

# usage: python make_excel.py <book_name> <out_dir> <tsv_dir>
book_name, out_dir, tsv_dir = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
out = out_dir / f"{book_name}.xlsx"

files = sorted(f for f in tsv_dir.glob("*.tsv") if not f.name.startswith("._")) # only the top level of tsv_dir
if not files:
    sys.exit(f"ERROR: no .tsv files in {tsv_dir}")

out_dir.mkdir(parents=True, exist_ok=True)

if out.exists():
    wb = load_workbook(out)                     # keep existing tabs
    default = None
else:
    wb = Workbook()                             # new workbook
    default = wb.active                         # its empty "Sheet", removed below

for f in files:
    sheet = f.stem[:31]                         # Excel sheet names are capped at 31 characters
    if sheet in wb.sheetnames:
        del wb[sheet]                           # replace this tab only
    df = pd.read_csv(f, sep="\t", keep_default_na=False, encoding="utf-8", encoding_errors="replace")
    ws = wb.create_sheet(title=sheet)
    ws.append(list(df.columns))
    for row in df.itertuples(index=False):
        ws.append(list(row))

if default is not None:
    wb.remove(default)

wb.save(out)
print(f"Updated {out} ({len(wb.sheetnames)} sheets)")