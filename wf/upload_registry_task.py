import logging

from latch import small_task
from latch.types import LatchDir, LatchFile
from latch.registry.table import Table
from wf.cleaning_task import CleaningOutput

from typing import List

logging.basicConfig(format="%(levelname)s - %(asctime)s - %(message)s")


@small_task(retries=0)
def upload_registry_task(
    cleaned_outputs: List[CleaningOutput], table_id: str = "761"
):
    table = Table(table_id)
    try:
        with table.update() as updater:
            for c in cleaned_outputs:
                updater.upsert_record(
                    c.run_id,
                    cleaned_fragment_file=LatchFile(
                        f"{c.cleaned_fragment_dir.remote_path}/cleaned_{c.run_id}_fragments.tsv.gz"
                    ),
                    positions_file=c.positions_file
                )
    except TypeError:
        print("error")
        logging.warning(f"No table with id {table_id} found.")
        return
    finally:
        return


if __name__ == "__main__":
    upload_registry_task(
        cleaned_outputs=[
            CleaningOutput(
                run_id="D1213",
                cleaned_fragment_dir=LatchDir(
                    "latch://13502.account/cleaned/Natrajan_cleaned_12042022"
                ),
            )
        ],
        table_id="761",
    )
