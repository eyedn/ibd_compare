###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_ilash_map.py
###############################################################################


from bisect import bisect_right
from pathlib import Path


def _format_cm(cm_value):
    formatted = f"{cm_value:.12f}".rstrip("0").rstrip(".")
    return formatted if formatted else "0"


def build_ilash_map(reference_map_path, ped_map_path, output_map_path):
    reference_rows = []
    with open(reference_map_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.strip().split()
            if len(fields) < 4:
                continue
            reference_rows.append(
                (fields[0], fields[1], float(fields[2]), int(fields[3]))
            )

    if not reference_rows:
        raise ValueError(f"No rows found in reference map: {reference_map_path}")

    reference_bps = [row[3] for row in reference_rows]
    reference_cms = [row[2] for row in reference_rows]

    output_map_path.parent.mkdir(parents=True, exist_ok=True)

    with open(ped_map_path, "r", encoding="utf-8") as ped_handle, open(
        output_map_path, "w", encoding="utf-8"
    ) as out_handle:
        for marker_idx, line in enumerate(ped_handle, start=1):
            if not line.strip():
                continue
            fields = line.strip().split()
            if len(fields) < 4:
                continue

            chrom = fields[0]
            snp_id = fields[1]
            bp = int(fields[3])

            if snp_id in {".", "0", ""}:
                snp_id = f"{chrom}_{bp}_{marker_idx}"

            pos = bisect_right(reference_bps, bp)
            if pos == 0:
                cm_value = reference_cms[0]
            elif pos >= len(reference_bps):
                cm_value = reference_cms[-1]
            elif reference_bps[pos - 1] == bp:
                cm_value = reference_cms[pos - 1]
            else:
                bp_lo = reference_bps[pos - 1]
                bp_hi = reference_bps[pos]
                cm_lo = reference_cms[pos - 1]
                cm_hi = reference_cms[pos]
                fraction = (bp - bp_lo) / (bp_hi - bp_lo)
                cm_value = cm_lo + fraction * (cm_hi - cm_lo)

            out_handle.write(
                f"{chrom}\t{snp_id}\t{_format_cm(cm_value)}\t{bp}\n"
            )
