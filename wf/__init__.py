"""Latch workflow for normalizing hot rows and columns in spatial ATAC-seq data
"""

from wf.cleaning_task import cleaning_task
from wf.upload_registry_task import upload_registry_task

from latch import workflow
from latch.resources.launch_plan import LaunchPlan
from latch.types import (
    LatchAuthor,
    LatchFile,
    LatchMetadata,
    LatchParameter,
    LatchRule
)

metadata = LatchMetadata(
    display_name="clean",
    author=LatchAuthor(
        name="James McGann",
        email="jpaulmcgann@gmail.com",
        github="https://github.com/jpmcga",
    ),
    repository="https://github.com/jpmcga/clean_latch",
    license="MIT",
    parameters={
        "run_id": LatchParameter(
            display_name="run id",
            description="ATX Run ID with optional prefix, default to \
                        Dxxxxx_NGxxxxx format.",
            batch_table_column=True,
            placeholder="Dxxxxx_NGxxxxx",
            rules=[
                LatchRule(
                    regex="^[^/].*",
                    message="run id cannot start with a '/'"
                )
            ]
        ),
        "output_dir": LatchParameter(
            display_name="output directory",
            description="Name of Latch subdirectory for downloaded file; files \
                will be saved to /cleaned/{output directory}.",
            batch_table_column=True,
            rules=[
                LatchRule(
                    regex="^[^/].*",
                    message="output directory name cannot start with a '/'"
                )
            ]
        ),
        "singlecell_file": LatchParameter(
            display_name="singlecell file",
            description="singlecell.csv file from cellranger output.",
            batch_table_column=True,
        ),
        "positions_file": LatchParameter(
            display_name="positions file",
            description="tissue_positions_list.csv file from spatial folder, \
                denoting on/off-tissue tixels.",
            batch_table_column=True, 
        ),
        "fragments_file": LatchParameter(
            display_name="fragments file",
            description="fragments.tsv.gz file from cellranger output",
            batch_table_column=True, 
        ),
        "deviations": LatchParameter(
            display_name="standard deviations",
            description="Number from standard deviations above the mean above \
            which row/column medians will be considered outliers.",
            batch_table_column=True, 
        ),
        "table_id": LatchParameter(
            display_name="Registry Table ID",
            description="Provide the ID of the Registry table. The cleaned fragment file will be populated in the table once the workflow finishes."
        )
    },
    tags=[],

)
@workflow(metadata)
def clean_workflow(
    run_id: str,
    output_dir: str,
    singlecell_file: LatchFile,
    positions_file: LatchFile,
    fragments_file: LatchFile,
    deviations: int=1,
    table_id: str="390"
) -> LatchFile:
    """Workflow for normalizing hot rows and columns in spatial ATAC-seq data

    ## Clean high fragment couns in spatial ATAC-seq data

    Latch workflow for remediating lane artifacts in spatial ATAC-seq 
    experiments, specifically for data generated via DBiT-seq
    (see [Deng, 2022](https://www.nature.com/articles/s41586-022-05094-1)).

   Workflow takes the following outputs from Cellragnger ATAC,
    - fragments.tsv.gz
    - singlecell.csv

    and returns a 'cleaned' fragments.tsv.gz where hot row/columns are 
    downsampled to be within x standard deviations of the mean of row/column
    medians.

    Output table is sorted for continuous chromosome blocks, and compressed
    with gbzip for ArchR.
    
    Output can be used for analysis with ArchR, Seurat, and other scATAC-seq
    packages.
    """
    
    cleaned_fragment_file = cleaning_task(
        run_id=run_id,
        output_dir=output_dir,
        singlecell_file=singlecell_file,
        positions_file=positions_file,
        fragments_file=fragments_file,
        deviations=deviations
    )

    upload_registry_task(
        run_id=run_id, 
        cleaned_fragment_file=cleaned_fragment_file, 
        table_id=table_id
    )

    return cleaned_fragment_file

LaunchPlan(
    clean_workflow,
    "test data",
    {
        "run_id" : "ds_D01033_NG01681",
        "output_dir" : "ds_D01033_NG01681",
        "singlecell_file" : LatchFile("latch:///cr_outs/ds_D01033_NG01681/outs/ds_D01033_NG01681_singlecell.csv"),
        "positions_file" : LatchFile("latch:///position_files/D01033/test.csv"),
        "fragments_file" : LatchFile("latch:///cr_outs/ds_D01033_NG01681/outs/ds_D01033_NG01681_fragments.tsv.gz"),
        "deviations" : 1,
    },
)
