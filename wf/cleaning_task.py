import logging
import subprocess

from latch import large_task
from latch.types import LatchFile


logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s",
    level=logging.INFO
)

@large_task
def cleaning_task(
    run_id: str,
    output_dir: str,
    singlecell_file: LatchFile,
    positions_file: LatchFile,
    fragments_file: LatchFile,
    deviations: int,
    ) -> LatchFile:

    _r_cmd = [
        "python",
        "/root/wf/clean.py",
        run_id,
        singlecell_file.local_path,
        positions_file.local_path,
        fragments_file.local_path,
        str(deviations),
    ]

    logging.info("cleaning...")
    subprocess.run(_r_cmd)
    out_table = f"{run_id}_fragments.tsv"

    _sort_cmd = [
        "sort",
        "-k1,1V",
        "-k2,2n",
        out_table
    ]
    
    logging.info("sorting...")
    subprocess.run(_sort_cmd, stdout=open(f"cleaned_{out_table}", "w"))

    _zip_cmd = [
        "bgzip",
        f"cleaned_{out_table}" 
    ]

    logging.info("zipping...")
    subprocess.run(_zip_cmd)
    out_file = f"cleaned_{out_table}.gz"

    local_location = f"/root/{out_file}"
    remote_location = f"latch:///cleaned/{output_dir}/{out_file}"

    return LatchFile(str(local_location), remote_location)

if __name__ == "__main__":
    cleaning_task(
        run_id="ds_D01033_NG01681",
        output_dir="ds_D01033_NG01681",
        singlecell_file = LatchFile("latch://13502.account/atac_outs/ds_D01033_NG01681/outs/ds_D01033_NG01681_singlecell.csv"),
        positions_file = LatchFile("latch://13502.account/spatials/demo/spatial/tissue_positions_list.csv"),
        fragments_file = LatchFile("latch://13502.account/atac_outs/ds_D01033_NG01681/outs/ds_D01033_NG01681_fragments.tsv.gz"),
        deviations=1,      
        )
