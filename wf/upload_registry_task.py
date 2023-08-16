from latch import small_task
from latch.types import LatchDir
from latch.registry.table import Table

@small_task(retries=0)
def upload_registry_task(
    cleaned_fragment_dir: LatchDir,
    run_id: str,
    table_id: str = "390"
):
    table = Table(table_id)

    with table.update() as updater:
        updater.upsert_record(
            run_id,
            cleaned_fragment_file=cleaned_fragment_dir
        )

if __name__ == "__main__":
    upload_registry_task(
        run_id="D01214_NG02300",
        cleaned_fragment_dir=LatchDir("latch:///cleaned/Natrajan_cleaned_12042022")
    )
