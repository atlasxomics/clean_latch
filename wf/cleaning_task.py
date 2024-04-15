import logging
import subprocess

from latch import large_task
from latch.types import LatchDir, LatchFile

from dataclasses import dataclass

logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s", level=logging.INFO
)


@dataclass
class Sample:
    run_id: str
    singlecell_file: LatchFile
    positions_file: LatchFile
    fragments_file: LatchFile
    output_dir: str
    deviations: int


@dataclass
class CleaningOutput:
    run_id: str
    cleaned_fragment_dir: LatchDir
    positions_file: LatchFile


@large_task
def cleaning_task(sample: Sample) -> CleaningOutput:
    _r_cmd = [
        "python",
        "/root/wf/clean.py",
        sample.run_id,
        sample.singlecell_file.local_path,
        sample.positions_file.local_path,
        sample.fragments_file.local_path,
        str(sample.deviations),
    ]

    logging.info("cleaning...")
    subprocess.run(_r_cmd)

    out_table = f"{sample.run_id}_fragments.tsv"
    out_metrics = f"{sample.run_id}_cleaning_metrics.csv"

    _sort_cmd = ["sort", "-k1,1V", "-k2,2n", out_table]

    logging.info("sorting...")
    subprocess.run(_sort_cmd, stdout=open(f"cleaned_{out_table}", "w"))

    _zip_cmd = ["bgzip", f"cleaned_{out_table}"]

    logging.info("zipping...")
    subprocess.run(_zip_cmd)
    out_file = f"cleaned_{out_table}.gz"

    subprocess.run(["mkdir", sample.output_dir])
    subprocess.run(["mv", out_file, out_metrics, sample.output_dir])

    local_dir = f"/root/{sample.output_dir}"
    remote_dir = f"latch:///cleaned/{sample.output_dir}"

    return CleaningOutput(
        run_id=sample.run_id,
        cleaned_fragment_dir=LatchDir(local_dir, remote_dir),
        positions_file=sample.positions_file
    )


if __name__ == "__main__":
    cleaning_task(
        sample=Sample(
            run_id="ds_D01033_NG01681",
            output_dir="hannah_ds_D01033_NG01681",
            singlecell_file = LatchFile("latch://13502.account/atac_outs/ds_D01033_NG01681/outs/ds_D01033_NG01681_singlecell.csv"),
            positions_file = LatchFile("latch://13502.account/spatials/demo/spatial/tissue_positions_list.csv"),
            fragments_file = LatchFile("latch://13502.account/atac_outs/ds_D01033_NG01681/outs/ds_D01033_NG01681_fragments.tsv.gz"),
            deviations=1, 
        )
    )
